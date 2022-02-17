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
import os
from copy import deepcopy
from timeit import default_timer as timer
import numpy as np

# triqs
from triqs.operators.util.observables import S_op, N_op
from triqs.version import git_hash as triqs_hash
from h5 import HDFArchive
import triqs.utility.mpi as mpi
from triqs.gf import BlockGf, GfReFreq, make_hermitian
from triqs.gf.tools import inverse
from triqs_dft_tools.sumk_dft import SumkDFT

# own modules
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


def _determine_block_structure(sum_k, general_params, advanced_params):
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
    mpi.report('\n *** determination of block structure ***')

    # this returns a list of dicts (one entry for each corr shell)
    # the dict contains one entry for up and one for down
    # each entry is a square complex numpy matrix with dim=corr_shell['dim']
    # GfReFreq for ftps, GfImFreq else
    if general_params['solver_type'] in ['ftps']:
        zero_Sigma_w = [sum_k.block_structure.create_gf(ish=iineq, gf_function=GfReFreq,
                                                        window = general_params['w_range'],
                                                        n_points = general_params['n_w'])
                        for iineq in range(sum_k.n_inequiv_shells)]
        sum_k.put_Sigma(zero_Sigma_w)
        G_loc_all = sum_k.extract_G_loc(iw_or_w='w', broadening = general_params['eta'], transform_to_solver_blocks=False)
        dens_mat = [G_loc_all[iineq].density() for iineq in range(sum_k.n_inequiv_shells)]
    else:
        zero_Sigma_iw = [sum_k.block_structure.create_gf(ish=iineq,
                                                         beta=general_params['beta'],
                                                         n_points = general_params['n_iw'])
                         for iineq in range(sum_k.n_inequiv_shells)]
        sum_k.put_Sigma(zero_Sigma_iw)
        dens_mat = sum_k.density_matrix(method='using_gf', beta=general_params['beta'])

    original_dens_mat = deepcopy(dens_mat)

    # if we want to do a magnetic calculation we need to lift up/down degeneracy
    if not general_params['csc'] and general_params['magnetic'] and sum_k.SO == 0:
        mpi.report('magnetic calculation: removing the spin degeneracy from the block structure')
        for i, elem in enumerate(dens_mat):
            for key, value in elem.items():
                if key == 'up':
                    for a in range(len(value[:, 0])):
                        for b in range(len(value[0, :])):
                            if a == b:
                                dens_mat[i][key][a, b] = value[a, b]*1.1
                elif key == 'down':
                    for a in range(len(value[:, 0])):
                        for b in range(len(value[0, :])):
                            if a == b:
                                dens_mat[i][key][a, b] = value[a, b]*0.9
                else:
                    mpi.report('warning spin channels not found! Doing a PM calculation')

    # for certain systems it is needed to keep off diag elements
    # this enforces to use the full corr subspace matrix
    solver_struct_ftps = None
    if general_params['enforce_off_diag'] or general_params['solver_type'] in ['ftps']:
        if general_params['solver_type'] in ['ftps']:
            # first round to determine real blockstructure
            mock_sumk = deepcopy(sum_k)
            mock_sumk.analyse_block_structure(dm=dens_mat, threshold=general_params['block_threshold'])
            solver_struct_ftps = [None] * sum_k.n_inequiv_shells
            for icrsh in range(sum_k.n_inequiv_shells):
                solver_struct_ftps[icrsh] = mock_sumk.deg_shells[icrsh]
            mpi.report('Block structure written to "solver_struct_ftps":')
            mpi.report(solver_struct_ftps)

        mpi.report('enforcing off-diagonal elements in block structure finder')
        for dens_mat_per_imp in dens_mat:
            for dens_mat_per_block in dens_mat_per_imp.values():
                dens_mat_per_block += 2 * general_params['block_threshold']

    if not general_params['enforce_off_diag'] and general_params['block_suppress_orbital_symm']:
        mpi.report('removing orbital symmetries in block structure finder')
        for dens_mat_per_imp in dens_mat:
            for dens_mat_per_block in dens_mat_per_imp.values():
                dens_mat_per_block += 2*np.diag(np.arange(dens_mat_per_block.shape[0]))

    mpi.report('using 1-particle density matrix and Hloc (atomic levels) to '
               'determine the block structure')
    sum_k.analyse_block_structure(dm=dens_mat, threshold=general_params['block_threshold'])

    # Applies the manual mapping to each inequivalent shell
    if advanced_params['map_solver_struct'] != 'none':
        # TODO: prints for debug only, remove
        mpi.report(sum_k.gf_struct_solver[0])
        mpi.report(sum_k.gf_struct_sumk[0])
        sum_k.block_structure.map_gf_struct_solver([advanced_params['map_solver_struct']] * sum_k.n_inequiv_shells)
        if advanced_params['mapped_solver_struct_degeneracies'] != 'none':
            sum_k.block_structure.deg_shells = [advanced_params['mapped_solver_struct_degeneracies']] * sum_k.n_inequiv_shells

    return sum_k, original_dens_mat, solver_struct_ftps


def _calculate_rotation_matrix(general_params, sum_k):
    """
    Applies rotation matrix to make the DMFT calculations easier for the solver.
    Possible are rotations diagonalizing either the local Hamiltonian or the
    density. Diagonalizing the density has not proven really helpful but
    diagonalizing the local Hamiltonian has.
    Note that the interaction Hamiltonian has to be rotated if it is not fully
    orbital-gauge invariant (only the Kanamori fulfills that).
    """

    # Extracts new rotation matrices from density_mat or local Hamiltonian
    if general_params['set_rot'] == 'hloc':
        q_diag = sum_k.eff_atomic_levels()
    elif general_params['set_rot'] == 'den':
        q_diag = sum_k.density_matrix(method='using_gf', beta=general_params['beta'])
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


def _chi_setup(sum_k, general_params, solver_params):
    """

    Parameters
    ----------
    sum_k : SumkDFT object
        Sumk object with the information about the correct block structure
    general_paramters: general params dict
    solver_params: solver params dict

    Returns
    -------
    solver_params :  dict
        solver_paramters for the QMC solver
    Op_list : list of one-particle operators to measure per impurity
    """

    if general_params['measure_chi'] == 'SzSz':
        mpi.report('\nSetting up Chi(S_z(tau),S_z(0)) measurement')
    elif general_params['measure_chi'] == 'NN':
        mpi.report('\nSetting up Chi(n(tau),n(0)) measurement')

    Op_list = [None] * sum_k.n_inequiv_shells

    for icrsh in range(sum_k.n_inequiv_shells):
        n_orb = sum_k.corr_shells[icrsh]['dim']
        orb_names = list(range(n_orb))

        if general_params['measure_chi'] == 'SzSz':
            Op_list[icrsh] = S_op('z',
                                  spin_names=sum_k.spin_block_names[sum_k.SO],
                                  orb_names=orb_names,
                                  map_operator_structure=sum_k.sumk_to_solver[icrsh])
        elif general_params['measure_chi'] == 'NN':
            Op_list[icrsh] = N_op(spin_names=sum_k.spin_block_names[sum_k.SO],
                                  orb_names=orb_names,
                                  map_operator_structure=sum_k.sumk_to_solver[icrsh])

    solver_params['measure_O_tau_min_ins'] = general_params['measure_chi_insertions']

    return solver_params, Op_list


def dmft_cycle(general_params, solver_params, advanced_params, dft_params,
               n_iter, dft_irred_kpt_indices=None, dft_energy=None):
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
    observables : dict
        current observable array for calculation
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

    # create Sumk object
    # TODO: use_dft_blocks=True yields inconsistent number of blocks!
    if general_params['csc']:
        sum_k = SumkDFT(hdf_file=general_params['seedname']+'.h5', use_dft_blocks=False,
                        h_field=general_params['h_field'])
    else:
        sum_k = SumkDFT(hdf_file=general_params['jobname']+'/'+general_params['seedname']+'.h5',
                        use_dft_blocks=False, h_field=general_params['h_field'])
        # This is a quick-and-dirty feature for a magnetic field with spin-orbit coupling
        # TODO: replace by more elegant implementation, e.g. new field in h5 giving spin
        if general_params['energy_shift_orbitals'] != 'none':
            assert np.allclose(sum_k.n_orbitals, sum_k.n_orbitals[0]), 'Energy shift in projector formalism not implemented'
            assert len(general_params['energy_shift_orbitals']) == sum_k.n_orbitals[0]
            sum_k.hopping += np.diag(general_params['energy_shift_orbitals'])

    iteration_offset = 0

    # determine chemical potential for bare DFT sum_k object
    if mpi.is_master_node():
        archive = HDFArchive(general_params['jobname']+'/'+general_params['seedname']+'.h5', 'a')
        if 'DMFT_results' not in archive:
            archive.create_group('DMFT_results')
        if 'last_iter' not in archive['DMFT_results']:
            archive['DMFT_results'].create_group('last_iter')
        if 'DMFT_input' not in archive:
            archive.create_group('DMFT_input')
            archive['DMFT_input'].create_group('solver')
            archive['DMFT_input'].create_group('version')
            archive['DMFT_input']['version']['triqs_hash'] = triqs_hash
        if 'iteration_count' in archive['DMFT_results']:
            iteration_offset = archive['DMFT_results/iteration_count']
            print('previous iteration count of {} '.format(iteration_offset)
                  + 'will be added to total number of iterations')
            sum_k.chemical_potential = archive['DMFT_results/last_iter/chemical_potential_post']
    else:
        archive = None

    iteration_offset = mpi.bcast(iteration_offset)
    sum_k.chemical_potential = mpi.bcast(sum_k.chemical_potential)

    # Incompatabilities for SO coupling
    if sum_k.SO == 1:
        if not general_params['csc'] and general_params['magnetic'] and general_params['afm_order']:
            raise ValueError('AFM order not supported with SO coupling')

    # need to set sigma immediately here, otherwise mesh in unclear for sumK
    # Initializes empty Sigma for calculation of DFT density even if block structure changes later
    if general_params['solver_type'] in ['ftps']:
        zero_Sigma_w = [sum_k.block_structure.create_gf(ish=iineq, gf_function=GfReFreq,
                                                        window = general_params['w_range'],
                                                        n_points = general_params['n_w'])
                        for iineq in range(sum_k.n_inequiv_shells)]
        sum_k.put_Sigma(zero_Sigma_w)
    else:
        zero_Sigma_iw = [sum_k.block_structure.create_gf(ish=iineq,
                                                         beta=general_params['beta'],
                                                         n_points = general_params['n_iw'])
                         for iineq in range(sum_k.n_inequiv_shells)]
        sum_k.put_Sigma(zero_Sigma_iw)

    # Sets the chemical potential of the DFT calculation
    # Either directly from general parameters, if given, ...
    if general_params['dft_mu'] != 'none':
        dft_mu = general_params['dft_mu']
        # Initializes chemical potential with dft_mu if this is the first iteration
        if iteration_offset == 0:
            sum_k.chemical_potential = dft_mu
            mpi.report('\n chemical potential set to {:.3f} eV\n'.format(sum_k.chemical_potential))
    # ... or with sum_k.calc_mu
    else:
        if general_params['solver_type'] in ['ftps']:
            dft_mu = sum_k.calc_mu(precision=general_params['prec_mu'],
                                   iw_or_w='w', broadening=general_params['eta'])
        else:
            dft_mu = sum_k.calc_mu(precision=general_params['prec_mu'], iw_or_w='iw')

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
    if general_params['set_rot'] != 'none':
        # calculate new rotation matrices
        sum_k = _calculate_rotation_matrix(general_params, sum_k)
    # Saves rotation matrix to h5 archive:
    if mpi.is_master_node() and iteration_offset == 0:
        archive['DMFT_input']['rot_mat'] = sum_k.rot_mat
    mpi.barrier()

    # determine block structure for solver
    det_blocks = None
    # load previous block_structure if possible
    if mpi.is_master_node():
        det_blocks = 'block_structure' not in archive['DMFT_input']
    det_blocks = mpi.bcast(det_blocks)

    # Previous rot_mat only not None if the rot_mat changed from load_sigma or previous run
    previous_rot_mat = None
    solver_struct_ftps = None
    # determine block structure for GF and Hyb function
    if det_blocks and not general_params['load_sigma']:
        sum_k, dm, solver_struct_ftps = _determine_block_structure(sum_k, general_params, advanced_params)
    # if load sigma we need to load everything from this h5 archive
    elif general_params['load_sigma']:
        #loading block_struc and rot_mat and deg_shells
        if mpi.is_master_node():
            with HDFArchive(general_params['path_to_sigma'], 'r') as old_calc:
                sum_k.block_structure = old_calc['DMFT_input/block_structure']
                sum_k.deg_shells = old_calc['DMFT_input/deg_shells']
                previous_rot_mat = old_calc['DMFT_input/rot_mat']
                if general_params['solver_type'] in ['ftps']:
                    solver_struct_ftps = old_calc['DMFT_input/solver_struct_ftps']

            if not np.allclose(sum_k.rot_mat, previous_rot_mat):
                print('WARNING: rot_mat in current run is different from loaded_sigma run.')
            else:
                previous_rot_mat = None

        sum_k.block_structure = mpi.bcast(sum_k.block_structure)
        sum_k.deg_shells = mpi.bcast(sum_k.deg_shells)
        previous_rot_mat = mpi.bcast(previous_rot_mat)
        solver_struct_ftps = mpi.bcast(solver_struct_ftps)

        # In a magnetic calculation, no shells are degenerate
        if not general_params['csc'] and general_params['magnetic'] and sum_k.SO == 0:
            sum_k.deg_shells = [[] for _ in range(sum_k.n_inequiv_shells)]
        dm = None
    else:
        # Master node checks if rot_mat stayed the same
        if mpi.is_master_node():
            sum_k.block_structure = archive['DMFT_input']['block_structure']
            sum_k.deg_shells = archive['DMFT_input/deg_shells']
            previous_rot_mat = archive['DMFT_input']['rot_mat']
            if not np.allclose(sum_k.rot_mat, previous_rot_mat):
                print('WARNING: rot_mat in current step is different from previous step.')
                archive['DMFT_input']['rot_mat'] = sum_k.rot_mat
            else:
                previous_rot_mat = None
            if general_params['solver_type'] in ['ftps']:
                solver_struct_ftps = archive['DMFT_input/solver_struct_ftps']

        sum_k.block_structure = mpi.bcast(sum_k.block_structure)
        sum_k.deg_shells = mpi.bcast(sum_k.deg_shells)
        previous_rot_mat = mpi.bcast(previous_rot_mat)
        solver_struct_ftps = mpi.bcast(solver_struct_ftps)
        dm = None

    # Compatibility with h5 archives from the triqs2 version
    # Sumk doesn't hold corr_to_inequiv anymore, which is in block_structure now
    if sum_k.block_structure.corr_to_inequiv is None:
        if mpi.is_master_node():
            sum_k.block_structure.corr_to_inequiv = archive['dft_input/corr_to_inequiv']
        sum_k.block_structure = mpi.bcast(sum_k.block_structure)

    # Determination of shell_multiplicity
    shell_multiplicity = [sum_k.corr_to_inequiv.count(icrsh) for icrsh in range(sum_k.n_inequiv_shells)]

    # Initializes new empty Sigma with new blockstructure for calculation of DFT density
    if general_params['solver_type'] in ['ftps']:
        zero_Sigma_w = [sum_k.block_structure.create_gf(ish=iineq, gf_function=GfReFreq,
                                                        window = general_params['w_range'],
                                                        n_points = general_params['n_w'])
                        for iineq in range(sum_k.n_inequiv_shells)]
        sum_k.put_Sigma(zero_Sigma_w)
    else:
        zero_Sigma_iw = [sum_k.block_structure.create_gf(ish=iineq,
                                                         beta=general_params['beta'],
                                                         n_points = general_params['n_iw'])
                         for iineq in range(sum_k.n_inequiv_shells)]
        sum_k.put_Sigma(zero_Sigma_iw)

    # print block structure and DFT input quantitites!
    formatter.print_block_sym(sum_k, dm, general_params)

    # extract free lattice greens function
    if general_params['solver_type'] in ['ftps']:
        G_loc_all_dft = sum_k.extract_G_loc(iw_or_w='w', broadening = general_params['eta'], with_Sigma=False,
                                            mu=dft_mu,
                                        )
    else:
        # set here with_Sigma=True because otherwise extract_G_loc assumes n_iw=1025
        G_loc_all_dft = sum_k.extract_G_loc(iw_or_w='iw', with_Sigma=True, mu=dft_mu)
    density_mat_dft = [G_loc_all_dft[iineq].density() for iineq in range(sum_k.n_inequiv_shells)]

    for iineq in range(sum_k.n_inequiv_shells):
        density_shell_dft = G_loc_all_dft[iineq].total_density()
        mpi.report('total density for imp {} from DFT: {:10.6f}'.format(iineq, np.real(density_shell_dft)))

    if not general_params['csc'] and general_params['magnetic']:
        sum_k.SP = 1

        if general_params['afm_order']:
            general_params = afm_mapping.determine(general_params, archive, sum_k.n_inequiv_shells)

    # Constructs interaction Hamiltonian and writes it to the h5 archive
    h_int = interaction_hamiltonian.construct(sum_k, general_params, advanced_params)
    if mpi.is_master_node():
        archive['DMFT_input']['h_int'] = h_int

    # If new calculation, writes input parameters and sum_k <-> solver mapping to archive
    if iteration_offset == 0:
        if mpi.is_master_node():
            archive['DMFT_input']['general_params'] = general_params
            archive['DMFT_input']['solver_params'] = solver_params
            archive['DMFT_input']['advanced_params'] = advanced_params

            archive['DMFT_input']['block_structure'] = sum_k.block_structure
            archive['DMFT_input']['deg_shells'] = sum_k.deg_shells
            archive['DMFT_input']['shell_multiplicity'] = shell_multiplicity
            if general_params['solver_type'] in ['ftps']:
                archive['DMFT_input']['solver_struct_ftps'] = solver_struct_ftps

    # Initializes the solvers
    solvers = [None] * sum_k.n_inequiv_shells
    for icrsh in range(sum_k.n_inequiv_shells):
        # Construct the Solver instances
        solvers[icrsh] = SolverStructure(general_params, solver_params, sum_k, icrsh, h_int[icrsh],
                                         iteration_offset, solver_struct_ftps)

    # store solver hash to archive
    if mpi.is_master_node():
        if 'version' not in archive['DMFT_input']:
            archive['DMFT_input'].create_group('version')
        archive['DMFT_input']['version']['solver_hash'] = solvers[0].git_hash

    # Determines initial Sigma and DC
    sum_k, solvers = initial_sigma.determine_dc_and_initial_sigma(general_params, advanced_params, sum_k,
                                                                  archive, iteration_offset, density_mat_dft, solvers)

    sum_k = manipulate_mu.set_initial_mu(general_params, sum_k, iteration_offset, archive, shell_multiplicity)

    # setup of measurement of chi(SzSz(tau) if requested
    if general_params['measure_chi'] != 'none':
        solver_params, Op_list = _chi_setup(sum_k, general_params, solver_params)
    else:
        Op_list = None

    mpi.report('\n {} DMFT cycles requested. Starting with iteration  {}.\n'.format(n_iter, iteration_offset+1))

    # Prepares observable and conv dicts
    observables = None
    conv_obs = None
    if mpi.is_master_node():
        observables = prep_observables(archive, sum_k)
        conv_obs = convergence.prep_conv_obs(archive, sum_k)
    observables = mpi.bcast(observables)
    conv_obs = mpi.bcast(conv_obs)

    if mpi.is_master_node() and iteration_offset == 0:
        write_header_to_file(general_params, sum_k)
        observables = add_dft_values_as_zeroth_iteration(observables, general_params, dft_mu, dft_energy, sum_k,
                                                         G_loc_all_dft, density_mat_dft, shell_multiplicity)
        write_obs(observables, sum_k, general_params)
        # write convergence file
        convergence.prep_conv_file(general_params, sum_k)

    # The infamous DMFT self consistency cycle
    is_converged = False
    for it in range(iteration_offset + 1, iteration_offset + n_iter + 1):
        mpi.report('#'*80)
        mpi.report('Running iteration: {} / {}'.format(it, iteration_offset + n_iter))
        (sum_k, solvers,
         observables, is_converged) = _dmft_step(sum_k, solvers, it, general_params,
                                                 solver_params, advanced_params, dft_params,
                                                 h_int, archive, shell_multiplicity, E_kin_dft,
                                                 observables, conv_obs, Op_list, dft_irred_kpt_indices, dft_energy,
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
            sum_k, solvers, observables, _ = _dmft_step(sum_k, solvers, it, general_params,
                                                        solver_params, advanced_params, dft_params,
                                                        h_int, archive, shell_multiplicity, E_kin_dft,
                                                        observables, conv_obs, Op_list, dft_irred_kpt_indices, dft_energy,
                                                        is_converged=True, is_sampling=True)

        mpi.report('** Sampling finished ***')
        mpi.report('#'*80)

    # for one-shot calculations, we can use the updated GAMMA file for postprocessing
    if not general_params['csc'] and general_params['oneshot_postproc_gamma_file']:
        # Write the density correction to file after the one-shot calculation
        sum_k.calc_density_correction(filename=os.path.join(general_params['jobname'], 'GAMMA'), dm_type='vasp')

    mpi.barrier()

    # close the h5 archive
    if mpi.is_master_node():
        del archive

    return is_converged


def _dmft_step(sum_k, solvers, it, general_params,
               solver_params, advanced_params, dft_params,
               h_int, archive, shell_multiplicity, E_kin_dft,
               observables, conv_obs, Op_list, dft_irred_kpt_indices, dft_energy,
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
    if general_params['solver_type'] in ['ftps']:
        G_loc_all = sum_k.extract_G_loc(iw_or_w='w', broadening=general_params['eta'])
    else:
        G_loc_all = sum_k.extract_G_loc(iw_or_w='iw')

    # Copies Sigma and G0 before Solver run for mixing later
    Sigma_freq_previous = [solvers[iineq].Sigma_freq.copy() for iineq in range(sum_k.n_inequiv_shells)]
    G0_freq_previous = [solvers[iineq].G0_freq.copy() for iineq in range(sum_k.n_inequiv_shells)]
    if general_params['dc'] and general_params['dc_type'] == 4:
        cpa_G_loc = gf_mixer.init_cpa(sum_k, solvers, general_params)

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

        # dyson equation to extract cpa_G0_freq
        if general_params['dc'] and general_params['dc_type'] == 4:
            cpa_G0_freq[icrsh] << inverse(solvers[icrsh].Sigma_freq + inverse(cpa_G_loc[icrsh]))

        # mixing of G0 if wanted from the second iteration on
        if it > 1:
            solvers[icrsh] = gf_mixer.mix_g0(solvers[icrsh], general_params, icrsh, archive,
                                             G0_freq_previous[icrsh], it, sum_k.deg_shells[icrsh])

        if general_params['solver_type'] in ['cthyb', 'ctint', 'hubbardI', 'inchworm']:
            solvers[icrsh].G0_freq << make_hermitian(solvers[icrsh].G0_freq)
        sum_k.symm_deg_gf(solvers[icrsh].G0_freq, ish=icrsh)

         # store solver to h5 archive
        if general_params['store_solver'] and mpi.is_master_node():
            archive['DMFT_input/solver'].create_group('it_'+str(it))
            archive['DMFT_input/solver/it_'+str(it)]['S_'+str(icrsh)] = solvers[icrsh].triqs_solver

        # store DMFT input directly in last_iter
        if mpi.is_master_node():
            archive['DMFT_results/last_iter']['G0_freq_{}'.format(icrsh)] = solvers[icrsh].G0_freq

        # setup of measurement of chi(SzSz(tau) if requested
        if general_params['measure_chi'] != 'none':
            solvers[icrsh].solver_params['measure_O_tau'] = (Op_list[icrsh], Op_list[icrsh])

        if (not general_params['csc'] and general_params['magnetic']
                and general_params['afm_order'] and general_params['afm_mapping'][icrsh][0]):
            # If we do a AFM calculation we can use the init magnetic moments to
            # copy the self energy instead of solving it explicitly
            solvers = afm_mapping.apply(general_params, solver_params, icrsh, sum_k.gf_struct_solver[icrsh], solvers)
        else:
            # Solve the impurity problem for this shell
            mpi.report('\nSolving the impurity problem for shell {} ...'.format(icrsh))
            mpi.barrier()
            start_time = timer()
            solvers[icrsh].solve()
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

        # add to cpa_G_time
        if general_params['dc'] and general_params['dc_type'] == 4:
            cpa_G_time << cpa_G_time + general_params['cpa_x'][icrsh] * solvers[icrsh].G_time

    # Done with loop over impurities

    if mpi.is_master_node():
        # Done. Now do post-processing:
        print('\n *** Post-processing the solver output ***')
        print('Total charge of all correlated shells : {:.6f}\n'.format(density_tot))

    # if CPA average Sigma over impurities before mixing
    if general_params['dc'] and general_params['dc_type'] == 4:
        solvers = gf_mixer.mix_cpa(cpa_G0_freq, sum_k.n_inequiv_shells, solvers)
    solvers = gf_mixer.mix_sigma(general_params, sum_k.n_inequiv_shells, solvers, Sigma_freq_previous)

    # calculate new DC
    if general_params['dc'] and general_params['dc_dmft']:
        sum_k = initial_sigma.calculate_double_counting(sum_k, density_mat,
                                                        general_params, advanced_params)

    # doing the dmft loop and set new sigma into sumk
    sum_k.put_Sigma([solvers[icrsh].Sigma_freq for icrsh in range(sum_k.n_inequiv_shells)])

    # saving previous mu for writing to observables file
    previous_mu = sum_k.chemical_potential
    sum_k = manipulate_mu.update_mu(general_params, sum_k, it, archive)

    # if we do a CSC calculation we need always an updated GAMMA file
    E_bandcorr = 0.0
    if general_params['csc']:
        # handling the density correction for fcsc calculations
        # TODO: keep only code in "else", the rest is redundant with handling inside
        #     sum_k.calc_density_correction. Code in "if" only implemented for
        #     compatibility with last official dft_tools release
        if dft_params['dft_code'] == 'vasp':
            if dft_irred_kpt_indices is None:
                deltaN, dens, E_bandcorr = sum_k.calc_density_correction(filename='GAMMA', dm_type='vasp')
            else:
                deltaN, dens, E_bandcorr = sum_k.calc_density_correction(filename='GAMMA', dm_type='vasp',
                                                                         kpts_to_write=dft_irred_kpt_indices)
        elif dft_params['dft_code'] == 'qe':
            deltaN, dens, E_bandcorr = sum_k.calc_density_correction(dm_type=dft_params['dft_code'])

    # for a one shot calculation we are using our own method
    if not general_params['csc'] and general_params['calc_energies']:
        E_bandcorr = calc_bandcorr_man(general_params, sum_k, E_kin_dft)

    if general_params['csc']:
        results_to_archive.write(archive, sum_k, general_params, solver_params, solvers, it,
                                 is_sampling, previous_mu, density_mat_pre, density_mat, deltaN, dens)
    else:
        results_to_archive.write(archive, sum_k, general_params, solver_params, solvers, it,
                                 is_sampling, previous_mu, density_mat_pre, density_mat)

    mpi.barrier()

    # calculate observables and write them to file
    if mpi.is_master_node():
        print('\n *** calculation of observables ***')
        observables = add_dmft_observables(observables,
                                           general_params,
                                           solver_params,
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
    if not general_params['csc'] and general_params['magnetic'] and sum_k.SO == 0:
        # if a magnetic calculation is done print out a summary of up/down occ
        formatter.print_summary_magnetic_occ(observables, sum_k.n_inequiv_shells)
    formatter.print_summary_convergence(conv_obs, general_params, sum_k.n_inequiv_shells)

    return sum_k, solvers, observables, is_converged
