################################################################################
#
# solid_dmft - A versatile python wrapper to perform DFT+DMFT calculations
#              utilizing the TRIQS software library
#
# Copyright (C) 2018-2020, ETH Zurich
# Copyright (C) 2021, The Simons Foundation
#      authors: A. Hampel, M. Merkel, and S. Beck
#
# solid_dmft is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# solid_dmft is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# solid_dmft (in the file COPYING.txt in this directory). If not, see
# <http://www.gnu.org/licenses/>.
#
################################################################################
"""
main DMFT cycle, DMFT step, and helper functions
"""


# system
from copy import deepcopy
from timeit import default_timer as timer
import numpy as np

# triqs
from triqs.operators.util.observables import S_op, N_op
from triqs.version import git_hash as triqs_hash
from triqs.version import version as triqs_version
from h5 import HDFArchive
import triqs.utility.mpi as mpi
from triqs.operators import c_dag, c, Operator
from triqs.gf import make_hermitian, fit_hermitian_tail, MeshReFreq, MeshImFreq, make_gf_from_fourier, iOmega_n
from triqs.gf.tools import inverse, make_zero_tail
from triqs_dft_tools.sumk_dft import SumkDFT

# own modules
from solid_dmft.version import solid_dmft_hash
from solid_dmft.version import version as solid_dmft_version
from solid_dmft.dmft_tools.observables import (calc_dft_kin_en, add_dmft_observables, calc_bandcorr_man, write_obs,
                                         add_dft_values_as_zeroth_iteration, write_header_to_file, prep_observables)
from solid_dmft.dmft_tools.solver import SolverStructure
from solid_dmft.dmft_tools import convergence
from solid_dmft.dmft_tools import formatter
from solid_dmft.dmft_tools import interaction_hamiltonian
from solid_dmft.dmft_tools import results_to_archive
from solid_dmft.dmft_tools import afm_mapping
from solid_dmft.dmft_tools import manipulate_chemical_potential as manipulate_mu
from solid_dmft.dmft_tools import initial_self_energies as initial_sigma
from solid_dmft.dmft_tools import greens_functions_mixer as gf_mixer
from solid_dmft.io_tools import verify_input_params, dict_to_h5

def _extract_quantity_per_inequiv(param_name, n_inequiv_shells, general_params):
    """
    For quantities that can be different for each inequivalent shell, this
    function checks if the quantity is a single value or a list of the correct
    length. If just a single value is given, this value is applied to each shell.

    Parameters
    ----------
    param_name : str
        name of the parameter to be checked
    n_inequiv_shells : int
        number of inequivalent shells
    general_params : dict
        general parameters as a dict

    Returns
    -------
    general_params : dict
        updated general parameters as a dict

    """


    def formatter(value):
        if isinstance(value, float):
            return '{:.2f}'.format(value)
        return str(value)

    param = general_params[param_name]
    if not isinstance(param, list):
        mpi.report('Assuming {} = '.format(param_name)
                + '{} for all correlated shells'.format(formatter(param)))
        general_params[param_name] = [param] * n_inequiv_shells
    elif len(param) == n_inequiv_shells:
        mpi.report(f'{param_name} list for correlated shells: '
                   + ', '.join([formatter(p) for p in param]))
    else:
        raise IndexError(f'{param_name} must be list of length '
                         f'n_inequiv_shells={n_inequiv_shells} '
                         'or single value but is ' + str(general_params[param_name]))

    return general_params

def _determine_block_structure(sum_k, general_params, advanced_params, solver_type_per_imp, dens_mat):
    """
    Determines block structrure and degenerate deg_shells
    computes first DFT density matrix to determine block structure and changes
    the density matrix according to needs i.e. magnetic calculations, or keep
    off-diag elements

    Parameters
    ----------
    sum_k : SumK Object instances

    Returns
    -------
    sum_k : SumK Object instances
        updated sum_k Object
    """
    mpi.report('\n *** Determination of block structure ***')

    # Adds spin degeneracies if not spin-orbit coupled
    if sum_k.SO == 0:
        sum_k.block_structure.deg_shells = [[['up_0', 'down_0']] for _ in range(sum_k.n_inequiv_shells)]

    # Evaluates degeneracies for ftps
    imp_ftps = [i for i, s in enumerate(solver_type_per_imp) if s == 'ftps']
    if imp_ftps:
        mock_sumk = deepcopy(sum_k)
        mock_sumk.analyse_block_structure(dm=dens_mat, threshold=general_params['block_threshold'], include_shells=imp_ftps)
        deg_orbs_ftps = [mock_sumk.deg_shells[icrsh] if icrsh in imp_ftps else None
                         for icrsh in range(sum_k.n_inequiv_shells)]
        mpi.report('Block structure written to "deg_orbs_ftps":')
        mpi.report(deg_orbs_ftps)
    else:
        deg_orbs_ftps = None

    # Only removes the off-diagonal terms for the selected impurities
    imp_to_analyze = [i for i, offdiag in enumerate(general_params['enforce_off_diag']) if not offdiag]
    mpi.report('using 1-particle density matrix and Hloc (atomic levels) to determine the block structure')
    sum_k.analyse_block_structure(dm=dens_mat, threshold=general_params['block_threshold'], include_shells=imp_to_analyze)

    # Applies manual selection of the solver struct
    if any(s is not None for s in advanced_params['pick_solver_struct']):
        mpi.report('selecting subset of orbital space for gf_struct_solver from input:')
        mpi.report(advanced_params['pick_solver_struct'][icrsh])
        sum_k.block_structure.pick_gf_struct_solver(advanced_params['pick_solver_struct'])

    # Applies the manual mapping to each inequivalent shell
    if any(s is not None for s in advanced_params['map_solver_struct']):
        sum_k.block_structure.map_gf_struct_solver(advanced_params['map_solver_struct'])
        if any(s is not None for s in advanced_params['mapped_solver_struct_degeneracies']):
            sum_k.block_structure.deg_shells = advanced_params['mapped_solver_struct_degeneracies']

    # if we want to do a magnetic calculation we need to lift up/down degeneracy
    if general_params['magnetic'] and sum_k.SO == 0:
        mpi.report('Magnetic calculation: removing the spin degeneracy from the block structure')
        for icrsh, deg_shells_site in enumerate(sum_k.block_structure.deg_shells):
            deg_shell_mag = []
            # find degenerate orbitals that do not simply connect different spin channels
            for deg_orbs in deg_shells_site:
                for spin in ['up', 'down']:
                    # create a list of all up / down orbitals
                    deg = [orb for orb in deg_orbs if spin in orb]
                    # if longer than one than we have two deg orbitals also in a magnetic calculation
                    if len(deg) > 1:
                        deg_shell_mag.append(deg)
            sum_k.block_structure.deg_shells[icrsh] = deg_shell_mag
    # for SOC we remove all degeneracies
    elif sum_k.SO == 1:
        sum_k.block_structure.deg_shells = [[] for _ in range(sum_k.n_inequiv_shells)]

    return sum_k, deg_orbs_ftps


def _calculate_rotation_matrix(general_params, sum_k):
    """
    Applies rotation matrix to make the DMFT calculations easier for the solver.
    Possible are rotations diagonalizing either the local Hamiltonian or the
    density. Diagonalizing the density has not proven really helpful but
    diagonalizing the local Hamiltonian has.
    Note that the interaction Hamiltonian has to be rotated if it is not fully
    orbital-gauge invariant (only the Kanamori fulfills that).

    Parameters
    ----------
    general_params : dict
        general parameters as a dict
    sum_k : SumkDFT object
        Sumk Object

    Returns
    -------
    sum_k : SumkDFT Object
        updated Sumk object with the new rotation matrices
    """

    # Extracts new rotation matrices from density_mat or local Hamiltonian
    if general_params['set_rot'] == 'hloc':
        q_diag = sum_k.eff_atomic_levels()
    elif general_params['set_rot'] == 'den':
        q_diag = sum_k.density_matrix(method='using_gf')
    else:
        raise ValueError('Parameter set_rot set to wrong value.')

    chnl = sum_k.spin_block_names[sum_k.SO][0]

    rot_mat = []
    for icrsh in range(sum_k.n_corr_shells):
        ish = sum_k.corr_to_inequiv[icrsh]
        eigvec = np.array(np.linalg.eigh(np.real(q_diag[ish][chnl]))[1], dtype=complex)
        if sum_k.use_rotations:
            rot_mat.append( np.dot(sum_k.rot_mat[icrsh], eigvec) )
        else:
            rot_mat.append( eigvec )

    sum_k.rot_mat = rot_mat
    # in case sum_k.use_rotations == False before:
    sum_k.use_rotations = True
    # sum_k.eff_atomic_levels() needs to be recomputed if rot_mat were changed
    if hasattr(sum_k, "Hsumk"): delattr(sum_k, "Hsumk")
    mpi.report('Updating rotation matrices using dft {} eigenbasis to maximise sign'.format(general_params['set_rot']))

    # Prints matrices
    mpi.report('\nNew rotation matrices')
    formatter.print_rotation_matrix(sum_k)

    return sum_k


def _chi_setup(sum_k, solver_params, map_imp_solver):
    """

    Parameters
    ----------
    sum_k : SumkDFT object
        Sumk object with the information about the correct block structure
    solver_params: solver params dict

    Returns
    -------
    ops_chi_measurement : list of one-particle operators to measure per impurity
    """

    ops_chi_measure = [None] * sum_k.n_inequiv_shells

    for icrsh in range(sum_k.n_inequiv_shells):
        isolv = map_imp_solver[icrsh]
        if 'measure_chi' not in solver_params[isolv] or solver_params[isolv]['measure_chi'] is None:
            continue

        n_orb = sum_k.corr_shells[icrsh]['dim']
        if solver_params[isolv]['measure_chi'] == 'SzSz':
            mpi.report(f'\nImp {icrsh} with solver #{isolv}: Setting up Chi(S_z(tau),S_z(0)) measurement')

            ops_chi_measure[icrsh] = S_op('z', spin_names=sum_k.spin_block_names[sum_k.SO],
                                          n_orb=n_orb, map_operator_structure=sum_k.sumk_to_solver[icrsh])
        elif solver_params[isolv]['measure_chi'] == 'NN':
            mpi.report(f'\nImp {icrsh} with solver #{isolv}: Setting up Chi(n(tau),n(0)) measurement')

            ops_chi_measure[icrsh] = N_op(spin_names=sum_k.spin_block_names[sum_k.SO],
                                          n_orb=n_orb, map_operator_structure=sum_k.sumk_to_solver[icrsh])

    return ops_chi_measure


def dmft_cycle(general_params, solver_params, advanced_params, dft_params,
               gw_params, n_iter, dft_irred_kpt_indices=None, dft_energy=None):
    """
    main dmft cycle that works for one shot and CSC equally

    Parameters
    ----------
    general_params : dict
        general parameters as a dict
    solver_params : dict
        solver parameters as a dict
    advanced_params : dict
        advanced parameters as a dict
    dft_params : dict
        dft parameters as a dict
    gw_params : dict
        gw parameters as a dict
    n_iter : int
        number of iterations to be executed
    dft_irred_kpt_indices: iterable of int
        If given, writes density correction for csc calculations only for
        irreducible kpoints

    Returns
    ---------
    observables : dict
        updated observable array for calculation
    """

    # Creates real- or imaginary-frequency mesh to store Green functions on
    if general_params['beta'] is not None:
        mpi.report('Running solid_dmft on imag-freq grid because "general.beta" specified')
        sumk_mesh = MeshImFreq(beta=general_params['beta'],
                               S='Fermion',
                               n_iw=general_params['n_iw'])
        broadening = None
    else:
        mpi.report('Running solid_dmft on real-freq grid because "general.beta" not specified')
        sumk_mesh = MeshReFreq(window=general_params['w_range'],
                               n_w=general_params['n_w'])
        broadening = general_params['eta']

    # Creates SumkDFT object
    sum_k = SumkDFT(hdf_file=general_params['jobname']+'/'+general_params['seedname']+'.h5',
                    mesh=sumk_mesh, use_dft_blocks=False, h_field=general_params['h_field'])

    iteration_offset = 0

    # if magnetic field is set we make sure that this is a magnetic calculation
    if general_params['h_field'] != 0.0:
        if not general_params['magnetic']:
            mpi.report('WARNING: Magnetic field set but general.magnetic is False. Setting it to True now.')
            general_params['magnetic'] = True

    # determine chemical potential for bare DFT sum_k object
    if mpi.is_master_node():
        archive = HDFArchive(general_params['jobname']+'/'+general_params['seedname']+'.h5', 'a')
        if 'DMFT_results' not in archive:
            archive.create_group('DMFT_results')
        if 'last_iter' not in archive['DMFT_results']:
            archive['DMFT_results'].create_group('last_iter')
        if 'DMFT_input' not in archive:
            archive.create_group('DMFT_input')
            archive['DMFT_input']['program'] = 'solid_dmft'
            archive['DMFT_input'].create_group('solver')
            archive['DMFT_input'].create_group('version')
            archive['DMFT_input']['version']['triqs_hash'] = triqs_hash
            archive['DMFT_input']['version']['triqs_version'] = triqs_version
            archive['DMFT_input']['version']['solid_dmft_hash'] = solid_dmft_hash
            archive['DMFT_input']['version']['solid_dmft_version'] = solid_dmft_version

        if 'iteration_count' in archive['DMFT_results']:
            iteration_offset = archive['DMFT_results/iteration_count']
            sum_k.chemical_potential = archive['DMFT_results/last_iter/chemical_potential_post']
            print(f'RESTARTING DMFT RUN at iteration {iteration_offset+1} using last self-energy')
        else:
            print('INITIAL DMFT RUN')
        print('#'*80, '\n')
    else:
        archive = None

    iteration_offset = mpi.bcast(iteration_offset)
    sum_k.chemical_potential = mpi.bcast(sum_k.chemical_potential)

    # Checks the general parameters that are impurity-dependent and ensures they are a list
    mpi.report('Checking parameters that are impurity-dependent')
    for key in ['U', 'J', 'U_prime', 'ratio_F4_F2', 'h_int_type', 'enforce_off_diag', 'dc_type']:
        general_params = _extract_quantity_per_inequiv(key, sum_k.n_inequiv_shells, general_params)
    # Same for advanced params
    for key in ['dc_U', 'dc_J', 'dc_fixed_occ', 'map_solver_struct', 'pick_solver_struct', 'mapped_solver_struct_degeneracies']:
        advanced_params = _extract_quantity_per_inequiv(key, sum_k.n_inequiv_shells, advanced_params)

    # Checks that all impurities have an associated solver
    # and creates a map between impurity index and solver index
    if len(solver_params) == 1 and solver_params[0]['idx_impurities'] is None:
        map_imp_solver = [0] * sum_k.n_inequiv_shells
    else:
        all_idx_imp = [i for entry in solver_params for i in entry['idx_impurities']]
        if sorted(all_idx_imp) != list(range(sum_k.n_inequiv_shells)):
            raise ValueError('All impurities must be listed exactly once in solver.idx_impurities'
                             f'but instead got {all_idx_imp}')

        map_imp_solver = []
        for iineq in range(sum_k.n_inequiv_shells):
            for isolver, entry in enumerate(solver_params):
                if iineq in entry['idx_impurities']:
                    map_imp_solver.append(isolver)
                    break
    solver_type_per_imp = [solver_params[map_imp_solver[iineq]]['type'] for iineq in range(sum_k.n_inequiv_shells)]
    mpi.report(f'\nSolver type per impurity: {solver_type_per_imp}')

    # Checks all parameters that need to be checked against info in SumkDFT object
    verify_input_params.verify_h5_dependent(sum_k, solver_type_per_imp, general_params)

    # Initializes chemical potential with mu_initial_guess if this is the first iteration
    if general_params['mu_initial_guess'] is not None and iteration_offset == 0:
        sum_k.chemical_potential = general_params['mu_initial_guess']
        mpi.report('\nInitial chemical potential set to {:.3f} eV\n'.format(sum_k.chemical_potential))

    dft_mu = sum_k.calc_mu(precision=general_params['prec_mu'], method=general_params['calc_mu_method'],
                           broadening=broadening)

    # calculate E_kin_dft for one shot calculations
    if not general_params['csc'] and general_params['calc_energies']:
        E_kin_dft = calc_dft_kin_en(general_params, sum_k, dft_mu)
    else:
        E_kin_dft = None

    # check for previous broyden data oterhwise initialize it:
    if mpi.is_master_node() and  general_params['g0_mix_type'] == 'broyden':
        if not 'broyler' in archive['DMFT_results']:
            archive['DMFT_results']['broyler'] = [{'mu' : [],'V': [], 'dV': [], 'F': [], 'dF': []}
                                                  for _ in range(sum_k.n_inequiv_shells)]

    # Generates a rotation matrix to change the basis
    if general_params['set_rot'] is not None:
        # calculate new rotation matrices
        sum_k = _calculate_rotation_matrix(general_params, sum_k)
    # Saves rotation matrix to h5 archive:
    if mpi.is_master_node() and iteration_offset == 0:
        archive['DMFT_input']['rot_mat'] = sum_k.rot_mat
    mpi.barrier()

    # Calculates density in default block structure
    local_gf_dft_corr = sum_k.extract_G_loc(broadening=broadening, transform_to_solver_blocks=False, with_Sigma=False)
    dens_mat_dft = [gf.density() for gf in local_gf_dft_corr]

    for iineq, gf in enumerate(local_gf_dft_corr):
        mpi.report(f'Total density for imp {iineq} from DFT: '
                   '{:10.6f}'.format(np.real(gf.total_density())))

    # determine block structure for solver
    det_blocks = None
    # load previous block_structure if possible
    if mpi.is_master_node():
        det_blocks = 'block_structure' not in archive['DMFT_input']
    det_blocks = mpi.bcast(det_blocks)

    # Previous rot_mat only not None if the rot_mat changed from load_sigma or previous run
    previous_rot_mat = None
    deg_orbs_ftps = None
    # determine block structure for GF and Hyb function
    if det_blocks and not general_params['load_sigma']:
        sum_k, deg_orbs_ftps = _determine_block_structure(sum_k, general_params, advanced_params, solver_type_per_imp, dens_mat_dft)
    # if load sigma we need to load everything from this h5 archive
    elif general_params['load_sigma']:
        #loading block_struc and rot_mat and deg_shells
        if mpi.is_master_node():
            with HDFArchive(general_params['path_to_sigma'], 'r') as old_calc:
                sum_k.block_structure = old_calc['DMFT_input/block_structure']
                sum_k.deg_shells = old_calc['DMFT_input/deg_shells']
                previous_rot_mat = old_calc['DMFT_input/rot_mat']
                if deg_orbs_ftps is not None:
                    deg_orbs_ftps = old_calc['DMFT_input/solver_struct_ftps']

            if not all(np.allclose(x, y) for x, y in zip(sum_k.rot_mat, previous_rot_mat)):
                print('WARNING: rot_mat in current run is different from loaded_sigma run.')
            else:
                previous_rot_mat = None

        sum_k.block_structure = mpi.bcast(sum_k.block_structure)
        sum_k.deg_shells = mpi.bcast(sum_k.deg_shells)
        previous_rot_mat = mpi.bcast(previous_rot_mat)
        deg_orbs_ftps = mpi.bcast(deg_orbs_ftps)

        # In a magnetic calculation, no shells are degenerate
        if general_params['magnetic'] and sum_k.SO == 0:
            sum_k.deg_shells = [[] for _ in range(sum_k.n_inequiv_shells)]
    else:
        # Master node checks if rot_mat stayed the same
        if mpi.is_master_node():
            sum_k.block_structure = archive['DMFT_input']['block_structure']
            sum_k.deg_shells = archive['DMFT_input/deg_shells']
            previous_rot_mat = archive['DMFT_input']['rot_mat']
            if not all(np.allclose(x, y) for x, y in zip(sum_k.rot_mat, previous_rot_mat)):
                print('WARNING: rot_mat in current step is different from previous step.')
                archive['DMFT_input']['rot_mat'] = sum_k.rot_mat
            else:
                previous_rot_mat = None
            if 'solver_struct_ftps' in archive['DMFT_input']:
                deg_orbs_ftps = archive['DMFT_input/solver_struct_ftps']

        sum_k.block_structure = mpi.bcast(sum_k.block_structure)
        sum_k.deg_shells = mpi.bcast(sum_k.deg_shells)
        previous_rot_mat = mpi.bcast(previous_rot_mat)
        deg_orbs_ftps = mpi.bcast(deg_orbs_ftps)

    # Compatibility with h5 archives from the triqs2 version
    # Sumk doesn't hold corr_to_inequiv anymore, which is in block_structure now
    if sum_k.block_structure.corr_to_inequiv is None:
        if mpi.is_master_node():
            sum_k.block_structure.corr_to_inequiv = archive['dft_input/corr_to_inequiv']
        sum_k.block_structure = mpi.bcast(sum_k.block_structure)

    # Determination of shell_multiplicity
    shell_multiplicity = [sum_k.corr_to_inequiv.count(icrsh) for icrsh in range(sum_k.n_inequiv_shells)]

    # print block structure and DFT input quantitites!
    formatter.print_block_sym(sum_k, dens_mat_dft, general_params)

    if general_params['magnetic']:
        sum_k.SP = 1
        if general_params['afm_order']:
            general_params = afm_mapping.determine(general_params, archive, sum_k.n_inequiv_shells)

    # Constructs interaction Hamiltonian and writes it to the h5 archive
    h_int, gw_params = interaction_hamiltonian.construct(sum_k, general_params, solver_type_per_imp, gw_params)
    if mpi.is_master_node():
        archive['DMFT_input']['h_int'] = h_int

    # If new calculation, writes input parameters and sum_k <-> solver mapping to archive
    if iteration_offset == 0:
        if mpi.is_master_node():
            archive['DMFT_input']['general_params'] = dict_to_h5.prep_params_for_h5(general_params)
            archive['DMFT_input']['solver_params'] = dict_to_h5.prep_params_for_h5(solver_params)
            archive['DMFT_input']['dft_params'] = dict_to_h5.prep_params_for_h5(dft_params)
            archive['DMFT_input']['advanced_params'] = dict_to_h5.prep_params_for_h5(advanced_params)

            archive['DMFT_input']['block_structure'] = sum_k.block_structure
            archive['DMFT_input']['deg_shells'] = sum_k.deg_shells
            archive['DMFT_input']['shell_multiplicity'] = shell_multiplicity
            if deg_orbs_ftps is not None:
                archive['DMFT_input']['solver_struct_ftps'] = deg_orbs_ftps
    mpi.barrier()

    solvers = [None] * sum_k.n_inequiv_shells
    for icrsh in range(sum_k.n_inequiv_shells):
        # Construct the Solver instances
        solvers[icrsh] = SolverStructure(general_params, solver_params[map_imp_solver[icrsh]],
                                         gw_params, advanced_params, sum_k, icrsh, h_int[icrsh],
                                         iteration_offset, deg_orbs_ftps)

    # store solver hash to archive
    if mpi.is_master_node():
        if 'version' not in archive['DMFT_input']:
            archive['DMFT_input'].create_group('version')
        for iineq, solver in enumerate(solvers):
            if f'solver {iineq}' not in archive['DMFT_input']['version']:
                archive['DMFT_input']['version'].create_group(f'solver {iineq}')
            archive['DMFT_input']['version'][f'solver {iineq}']['name'] = solver.solver_params['type']
            archive['DMFT_input']['version'][f'solver {iineq}']['hash'] = solver.git_hash
            archive['DMFT_input']['version'][f'solver {iineq}']['version'] = solver.version

    # Extracts local GF per *inequivalent* shell
    local_gf_dft = sum_k.extract_G_loc(broadening=broadening, with_Sigma=False, mu=dft_mu)

    # Determines initial Sigma and DC
    sum_k, solvers = initial_sigma.determine_dc_and_initial_sigma(general_params, gw_params, advanced_params, sum_k,
                                                                  archive, iteration_offset, local_gf_dft, solvers,
                                                                  solver_type_per_imp)

    sum_k = manipulate_mu.set_initial_mu(general_params, sum_k, iteration_offset, archive, broadening)


    # setup of measurement of chi(SzSz(tau) if requested
    ops_chi_measure = _chi_setup(sum_k, solver_params, map_imp_solver)

    mpi.report('\n {} DMFT cycles requested. Starting with iteration  {}.\n'.format(n_iter, iteration_offset+1))

    # Prepares observable and conv dicts
    observables = None
    conv_obs = None
    if mpi.is_master_node():
        observables = prep_observables(archive, sum_k)
        conv_obs = convergence.prep_conv_obs(archive)
    observables = mpi.bcast(observables)
    conv_obs = mpi.bcast(conv_obs)

    if mpi.is_master_node() and iteration_offset == 0:
        write_header_to_file(general_params, sum_k)
        observables = add_dft_values_as_zeroth_iteration(observables, general_params, solver_type_per_imp, dft_mu, dft_energy, sum_k,
                                                         local_gf_dft, shell_multiplicity)
        write_obs(observables, sum_k, general_params)
        # write convergence file
        convergence.prep_conv_file(general_params, sum_k)

    # The infamous DMFT self consistency cycle
    is_converged = False
    for it in range(iteration_offset + 1, iteration_offset + n_iter + 1):

        # remove h_field when number of iterations is reached
        if sum_k.h_field != 0.0 and general_params['h_field_it'] != -1 and it > general_params['h_field_it']:
            mpi.report('\nRemoving magnetic field now.\n')
            sum_k.h_field = 0.0
            # enforce recomputation of eff_atomic_levels
            delattr(sum_k, 'Hsumk')

        mpi.report('#'*80)
        mpi.report('Running iteration: {} / {}'.format(it, iteration_offset + n_iter))
        (sum_k, solvers,
         observables, is_converged) = _dmft_step(sum_k, solvers, it, general_params, solver_params, gw_params,
                                                 advanced_params, dft_params, map_imp_solver, solver_type_per_imp,
                                                 h_int, archive, shell_multiplicity, E_kin_dft,
                                                 observables, conv_obs, ops_chi_measure, dft_irred_kpt_indices, dft_energy, broadening,
                                                 is_converged, is_sampling=False)

        if is_converged:
            break

    if is_converged:
        mpi.report('*** Required convergence reached ***')
    else:
        mpi.report('** All requested iterations finished ***')
    mpi.report('#'*80)

    # Starts the sampling dmft iterations if requested
    if is_converged and general_params['sampling_iterations'] > 0:
        mpi.report('*** Sampling now for {} iterations ***'.format(general_params['sampling_iterations']))
        iteration_offset = it

        for it in range(iteration_offset + 1,
                        iteration_offset + 1 + general_params['sampling_iterations']):
            mpi.report('#'*80)
            mpi.report('Running iteration: {} / {}'.format(it, iteration_offset+general_params['sampling_iterations']))
            sum_k, solvers, observables, _ = _dmft_step(sum_k, solvers, it, general_params, solver_params, gw_params,
                                                        advanced_params, dft_params, map_imp_solver, solver_type_per_imp,
                                                        h_int, archive, shell_multiplicity, E_kin_dft,
                                                        observables, conv_obs, ops_chi_measure, dft_irred_kpt_indices, dft_energy, broadening,
                                                        is_converged=True, is_sampling=True)

        mpi.report('** Sampling finished ***')
        mpi.report('#'*80)

    mpi.barrier()

    # close the h5 archive
    if mpi.is_master_node():
        del archive

    return is_converged, sum_k


def _dmft_step(sum_k, solvers, it, general_params, solver_params, gw_params,
               advanced_params, dft_params, map_imp_solver, solver_type_per_imp, h_int,
               archive, shell_multiplicity, E_kin_dft, observables, conv_obs,
               ops_chi_measure, dft_irred_kpt_indices, dft_energy, broadening,
               is_converged, is_sampling):
    """
    Contains the actual dmft steps when all the preparation is done
    """

    # init local density matrices for observables
    density_tot = 0.0
    density_shell = np.zeros(sum_k.n_inequiv_shells)
    density_mat = [None] * sum_k.n_inequiv_shells
    density_mat_unsym = [None] * sum_k.n_inequiv_shells
    density_shell_pre = np.zeros(sum_k.n_inequiv_shells)
    density_mat_pre = [None] * sum_k.n_inequiv_shells

    mpi.barrier()

    if sum_k.SO:
        printed = ((np.real, 'real'), (np.imag, 'imaginary'))
    else:
        printed = ((np.real, 'real'), )

    # Extracts G local
    G_loc_all = sum_k.extract_G_loc(broadening=broadening)

    # Copies Sigma and G0 before Solver run for mixing later
    Sigma_freq_previous = [solvers[iineq].Sigma_freq.copy() for iineq in range(sum_k.n_inequiv_shells)]
    G0_freq_previous = [solvers[iineq].G0_freq.copy() for iineq in range(sum_k.n_inequiv_shells)]

    # looping over inequiv shells and solving for each site seperately
    for icrsh in range(sum_k.n_inequiv_shells):
        # copy the block of G_loc into the corresponding instance of the impurity solver
        # TODO: why do we set solvers.G_freq? Isn't that simply an output of the solver?
        solvers[icrsh].G_freq << G_loc_all[icrsh]

        density_shell_pre[icrsh] = np.real(solvers[icrsh].G_freq.total_density())
        mpi.report('\n *** Correlated Shell type #{:3d} : '.format(icrsh)
                   + 'Estimated total charge of impurity problem = {:.6f}'.format(density_shell_pre[icrsh]))
        density_mat_pre[icrsh] = solvers[icrsh].G_freq.density()
        mpi.report('Estimated density matrix:')
        for key, value in sorted(density_mat_pre[icrsh].items()):
            for func, name in printed:
                mpi.report('{}, {} part'.format(key, name))
                mpi.report(func(value))

        # dyson equation to extract G0_freq, using Hermitian symmetry
        solvers[icrsh].G0_freq << inverse(solvers[icrsh].Sigma_freq + inverse(solvers[icrsh].G_freq))

        # mixing of G0 if wanted from the second iteration on
        if it > 1:
            solvers[icrsh] = gf_mixer.mix_g0(solvers[icrsh], general_params, icrsh, archive,
                                             G0_freq_previous[icrsh], it, sum_k.deg_shells[icrsh])

        if isinstance(sum_k.mesh, MeshImFreq):
            solvers[icrsh].G0_freq << make_hermitian(solvers[icrsh].G0_freq)
        sum_k.symm_deg_gf(solvers[icrsh].G0_freq, ish=icrsh)

        if ((solver_type_per_imp[icrsh] == 'cthyb' and solvers[icrsh].solver_params['delta_interface'])
                or solver_type_per_imp[icrsh] == 'ctseg'):
            mpi.report('\n Using the delta interface for passing Delta(tau) and Hloc0 directly to the solver.')
             # prepare solver input
            sumk_eal = sum_k.eff_atomic_levels()[icrsh]
            solver_eal = sum_k.block_structure.convert_matrix(sumk_eal, space_from='sumk', ish_from=sum_k.inequiv_to_corr[icrsh])
            # fill Delta_time from Delta_freq sum_k to solver
            # for name, g0 in self.G0_freq:
            for name, g0 in solvers[icrsh].G0_freq:
                solvers[icrsh].Delta_freq[name] << iOmega_n - inverse(g0) - solver_eal[name]
                known_moments = make_zero_tail(solvers[icrsh].Delta_freq[name], 1)
                tail, err = fit_hermitian_tail(solvers[icrsh].Delta_freq[name], known_moments)
                # without SOC delta_tau needs to be real
                if not sum_k.SO == 1:
                    solvers[icrsh].Delta_time[name] << make_gf_from_fourier(solvers[icrsh].Delta_freq[name],
                                                                            solvers[icrsh].Delta_time.mesh, tail).real
                else:
                    solvers[icrsh].Delta_time[name] << make_gf_from_fourier(solvers[icrsh].Delta_freq[name],
                                                                           solvers[icrsh].Delta_time.mesh, tail)

                if solvers[icrsh].solver_params['diag_delta']:
                    for o1 in range(g0.target_shape[0]):
                        for o2 in range(g0.target_shape[0]):
                            if o1 != o2:
                                solvers[icrsh].Delta_time[name].data[:, o1, o2] = 0.0 + 0.0j

                # Make non-interacting operator for Hloc0
                Hloc_0 = Operator()
                for spin, spin_block in solver_eal.items():
                    for o1 in range(spin_block.shape[0]):
                        for o2 in range(spin_block.shape[1]):
                            # check if off-diag element is larger than threshold
                            if o1 != o2 and abs(spin_block[o1,o2]) < solvers[icrsh].solver_params['off_diag_threshold']:
                                continue
                            else:
                                # TODO: adapt for SOC calculations, which should keep the imag part
                                Hloc_0 += spin_block[o1,o2].real/2 * (c_dag(spin,o1) * c(spin,o2) + c_dag(spin,o2) * c(spin,o1))
                solvers[icrsh].Hloc_0 = Hloc_0

         # store solver to h5 archive
        if general_params['store_solver'] and mpi.is_master_node():
            if 'solver' not in archive['DMFT_input']:
                archive['DMFT_input'].create_group('solver')
            archive['DMFT_input/solver'].create_group('it_'+str(it))
            archive['DMFT_input/solver/it_'+str(it)]['S_'+str(icrsh)] = solvers[icrsh].triqs_solver

        # store DMFT input directly in last_iter
        if mpi.is_master_node():
            archive['DMFT_results/last_iter']['G0_freq_{}'.format(icrsh)] = solvers[icrsh].G0_freq
            if solver_type_per_imp[icrsh] == 'cthyb' and solvers[icrsh].solver_params['delta_interface']:
                archive['DMFT_results/last_iter']['Delta_time_{}'.format(icrsh)] = solvers[icrsh].Delta_time

        # setup of measurement of chi(SzSz(tau) if requested
        # TODO: move this into solver class?
        if ops_chi_measure[icrsh] is not None:
            solvers[icrsh].solver_params['measure_O_tau'] = (ops_chi_measure[icrsh], ops_chi_measure[icrsh])

        if (general_params['magnetic'] and general_params['afm_order'] and general_params['afm_mapping'][icrsh][0]):
            # If we do a AFM calculation we can use the init magnetic moments to
            # copy the self energy instead of solving it explicitly
            solvers = afm_mapping.apply(general_params, icrsh, sum_k.gf_struct_solver[icrsh], solvers)
        else:
            # Solve the impurity problem for this shell
            mpi.report('\nSolving the impurity problem for shell {} ...'.format(icrsh))
            mpi.barrier()
            start_time = timer()
            solvers[icrsh].solve(it=it)
            mpi.barrier()
            mpi.report('Actual time for solver: {:.2f} s'.format(timer() - start_time))

        # some printout of the obtained density matrices and some basic checks from the unsymmetrized solver output
        density_shell[icrsh] = np.real(solvers[icrsh].G_freq_unsym.total_density())
        density_tot += density_shell[icrsh]*shell_multiplicity[icrsh]
        density_mat_unsym[icrsh] = solvers[icrsh].G_freq_unsym.density()
        density_mat[icrsh] = solvers[icrsh].G_freq.density()
        formatter.print_local_density(density_shell[icrsh], density_shell_pre[icrsh],
                                      density_mat_unsym[icrsh], sum_k.SO)

        # update solver in h5 archive
        if general_params['store_solver'] and mpi.is_master_node():
            archive['DMFT_input/solver/it_'+str(it)]['S_'+str(icrsh)] = solvers[icrsh].triqs_solver

    # Done with loop over impurities

    if mpi.is_master_node():
        # Done. Now do post-processing:
        print('\n *** Post-processing the solver output ***')
        print('Total charge of all correlated shells : {:.6f}\n'.format(density_tot))

    solvers = gf_mixer.mix_sigma(general_params, sum_k.n_inequiv_shells, solvers, Sigma_freq_previous)

    # calculate new DC
    # for the hartree solver the DC potential will be formally set to zero as it is already present in the Sigma
    if general_params['dc'] and general_params['dc_dmft']:
        sum_k = initial_sigma.calculate_double_counting(sum_k, density_mat, general_params,
                                                        gw_params, advanced_params,
                                                        solver_type_per_imp)

    #The hartree solver computes the DC energy internally, set it in sum_k
    for icrsh in range(sum_k.n_corr_shells):
        iineq = sum_k.corr_to_inequiv[icrsh]
        if solver_type_per_imp[iineq] == 'hartree':
            sum_k.dc_energ[icrsh] = solvers[iineq].DC_energy

    # symmetrize Sigma over degenerate blocks
    for icrsh in range(sum_k.n_inequiv_shells):
        sum_k.symm_deg_gf(solvers[icrsh].Sigma_freq, ish=icrsh)
    # doing the dmft loop and set new sigma into sumk
    sum_k.put_Sigma([solvers[icrsh].Sigma_freq for icrsh in range(sum_k.n_inequiv_shells)])

    # saving previous mu for writing to observables file
    previous_mu = sum_k.chemical_potential
    sum_k = manipulate_mu.update_mu(general_params, sum_k, it, archive, broadening)

    # if we do a CSC calculation we need always an updated GAMMA file
    E_bandcorr = 0.0
    deltaN = None
    dens = None
    if general_params['csc']:
        # handling the density correction for fcsc calculations
        assert dft_irred_kpt_indices is None or dft_params['dft_code'] == 'vasp'
        deltaN, dens, E_bandcorr = sum_k.calc_density_correction(dm_type=dft_params['dft_code'],
                                                                 kpts_to_write=dft_irred_kpt_indices)
    elif general_params['calc_energies']:
        # for a one shot calculation we are using our own method
        E_bandcorr = calc_bandcorr_man(general_params, sum_k, E_kin_dft)

    # Writes results to h5 archive
    results_to_archive.write(archive, sum_k, general_params, solver_params, solvers, map_imp_solver, solver_type_per_imp, it,
                             is_sampling, previous_mu, density_mat_pre, density_mat, deltaN, dens)

    mpi.barrier()

    # calculate observables and write them to file
    if mpi.is_master_node():
        print('\n *** calculation of observables ***')
        observables = add_dmft_observables(observables,
                                           general_params,
                                           solver_params,
                                           map_imp_solver,
                                           solver_type_per_imp,
                                           dft_energy,
                                           it,
                                           solvers,
                                           h_int,
                                           previous_mu,
                                           sum_k,
                                           density_mat,
                                           shell_multiplicity,
                                           E_bandcorr)

        write_obs(observables, sum_k, general_params)

        # write the new observable array to h5 archive
        archive['DMFT_results']['observables'] = observables

    # Computes convergence quantities and writes them to file
    if mpi.is_master_node():
        conv_obs = convergence.calc_convergence_quantities(sum_k, general_params, conv_obs, observables,
                                                           solvers, G0_freq_previous, G_loc_all, Sigma_freq_previous)
        convergence.write_conv(conv_obs, sum_k, general_params)
        archive['DMFT_results']['convergence_obs'] = conv_obs
    conv_obs = mpi.bcast(conv_obs)

    mpi.report('*** iteration finished ***')

    # Checks for convergence
    is_now_converged = convergence.check_convergence(sum_k.n_inequiv_shells, general_params, conv_obs)
    if is_now_converged is None:
        is_converged = False
    else:
        # if convergency criteria was already reached don't overwrite it!
        is_converged = is_converged or is_now_converged

    # Final prints
    formatter.print_summary_observables(observables, sum_k.n_inequiv_shells,
                                        sum_k.spin_block_names[sum_k.SO])
    if general_params['calc_energies']:
        formatter.print_summary_energetics(observables)
    if general_params['magnetic'] and sum_k.SO == 0:
        # if a magnetic calculation is done print out a summary of up/down occ
        formatter.print_summary_magnetic_occ(observables, sum_k.n_inequiv_shells)
    formatter.print_summary_convergence(conv_obs, general_params, sum_k.n_inequiv_shells)

    return sum_k, solvers, observables, is_converged
