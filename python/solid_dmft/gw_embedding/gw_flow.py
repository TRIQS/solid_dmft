# -*- coding: utf-8 -*-
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
# pyright: reportUnusedExpression=false
"""
Module for gw flow
"""

from timeit import default_timer as timer
import numpy as np

from h5 import HDFArchive
from triqs.utility import mpi
from triqs.gfs import (
    iOmega_n,
    inverse,
    fit_hermitian_tail,
    Gf,
    BlockGf,
    make_hermitian,
    make_gf_dlr,
    make_gf_imfreq,
    make_gf_imtime,
    make_gf_dlr_imfreq,
)
from triqs.version import git_hash as triqs_hash
from triqs.version import version as triqs_version
from triqs.mesh import MeshImFreq
from triqs.operators import c_dag, c, Operator
from triqs_dft_tools.block_structure import BlockStructure
from triqs_dft_tools.sumk_dft import SumkDFT

from solid_dmft.version import solid_dmft_hash
from solid_dmft.version import version as solid_dmft_version
from solid_dmft.dmft_tools import formatter
from solid_dmft.dmft_tools import results_to_archive
from solid_dmft.dmft_tools.solver import create_solver
from solid_dmft.dmft_tools import interaction_hamiltonian
from solid_dmft.dmft_cycle import _extract_quantity_per_inequiv, _determine_block_structure
from solid_dmft.dmft_tools import greens_functions_mixer as gf_mixer
from solid_dmft.gw_embedding.bdft_converter import convert_gw_output


class dummy_sumk(object):
    """
    create dummy sumk helper object
    """

    def __init__(self, n_inequiv_shells, n_orb_list, enforce_off_diag, use_rot, magnetic):
        self.n_inequiv_shells = n_inequiv_shells
        self.SO = 0
        self.use_rotations = use_rot
        if self.use_rotations:
            raise ValueError('rotations not implemented yet for GW embedding')
        self.gf_struct_solver = []
        self.gf_struct_sumk = []
        self.spin_block_names = []
        self.inequiv_to_corr = []
        self.corr_to_inequiv = []
        self.deg_shells = []
        self.dc_energ = [0.0 for ish in range(self.n_inequiv_shells)]
        self.sumk_to_solver = [{} for ish in range(self.n_inequiv_shells)]
        self.solver_to_sumk = [{} for ish in range(self.n_inequiv_shells)]
        self.solver_to_sumk_block = [{} for ish in range(self.n_inequiv_shells)]
        for ish in range(self.n_inequiv_shells):
            self.inequiv_to_corr.append(ish)
            self.corr_to_inequiv.append(ish)
            self.spin_block_names.append(['up', 'down'])
            self.gf_struct_sumk.append([('up', n_orb_list[ish]), ('down', n_orb_list[ish])])

            # use full off-diagonal block structure in impurity solver
            if enforce_off_diag:
                self.gf_struct_solver.append({'up_0': n_orb_list[ish], 'down_0': n_orb_list[ish]})
                if not magnetic:
                    self.deg_shells.append([['up_0', 'down_0']])
                # setup standard mapping between sumk and solver
                for block, inner_dim in self.gf_struct_sumk[ish]:
                    self.solver_to_sumk_block[ish][f'{block}_0'] = block
                    for iorb in range(inner_dim):
                        self.sumk_to_solver[ish][(block, iorb)] = (block + '_0', iorb)
                        self.solver_to_sumk[ish][(block + '_0', iorb)] = (block, iorb)
            else:
                self.gf_struct_solver.append({})
                self.deg_shells.append([])
                for block, inner_dim in self.gf_struct_sumk[ish]:
                    for iorb in range(inner_dim):
                        self.gf_struct_solver[ish][f'{block}_{iorb}'] = 1
                        if not magnetic and block == 'up':
                            self.deg_shells[ish].append([f'up_{iorb}', f'down_{iorb}'])
                        # setup standard mapping between sumk and solver
                        self.solver_to_sumk_block[ish][f'{block}_{iorb}'] = block
                        self.sumk_to_solver[ish][(block, iorb)] = (f'{block}_{iorb}', 0)
                        self.solver_to_sumk[ish][(f'{block}_{iorb}', 0)] = (block, iorb)


        self.gf_struct_solver_list = [sorted([(k, v) for k, v in list(gfs.items())], key=lambda x: x[0]) for gfs in self.gf_struct_solver]

        # creat block_structure object for solver
        self.block_structure = BlockStructure(
            gf_struct_sumk=self.gf_struct_sumk,
            gf_struct_solver=self.gf_struct_solver,
            solver_to_sumk=self.solver_to_sumk,
            sumk_to_solver=self.sumk_to_solver,
            solver_to_sumk_block=self.solver_to_sumk_block,
            deg_shells=self.deg_shells,
            corr_to_inequiv = self.corr_to_inequiv,
            transformation=None,
        )

    def symm_deg_gf(self, gf_to_symm, ish=0):
        r"""
        Averages a GF or a dict of np.ndarrays over degenerate shells.

        Degenerate shells of an inequivalent correlated shell are defined by
        `self.deg_shells`. This function enforces corresponding degeneracies
        in the input GF.

        Parameters
        ----------
        gf_to_symm : gf_struct_solver like
                     Input and output GF (i.e., it gets overwritten)
                     or dict of np.ndarrays.
        ish : int
              Index of an inequivalent shell. (default value 0)

        """

        # when reading block_structures written with older versions from
        # an h5 file, self.deg_shells might be None
        if self.deg_shells is None:
            return

        if not isinstance(gf_to_symm, BlockGf) and isinstance(gf_to_symm[list(gf_to_symm.keys())[0]], np.ndarray):
            blockgf = False
        elif isinstance(gf_to_symm, BlockGf):
            blockgf = True
        else:
            raise ValueError("gf_to_symm should be either a BlockGf or a dict of numpy arrays")

        for degsh in self.deg_shells[ish]:
            # ss will hold the averaged orbitals in the basis where the
            # blocks are all equal
            # i.e. maybe_conjugate(v^dagger gf v)
            ss = None
            n_deg = len(degsh)
            for key in degsh:
                if ss is None:
                    if blockgf:
                        ss = gf_to_symm[key].copy()
                        ss.zero()
                        helper = ss.copy()
                    else:
                        ss = np.zeros_like(gf_to_symm[key])
                        helper = np.zeros_like(gf_to_symm[key])

                # get the transformation matrix
                if isinstance(degsh, dict):
                    v, C = degsh[key]
                else:
                    # for backward compatibility, allow degsh to be a list
                    if blockgf:
                        v = np.eye(*ss.target_shape)
                    else:
                        v = np.eye(*ss.shape)
                    C = False
                # the helper is in the basis where the blocks are all equal
                if blockgf:
                    helper.from_L_G_R(v.conjugate().transpose(), gf_to_symm[key], v)
                else:
                    helper = np.dot(v.conjugate().transpose(), np.dot(gf_to_symm[key], v))

                if C:
                    helper << helper.transpose()
                # average over all shells
                ss += helper / (1.0 * n_deg)
            # now put back the averaged gf to all shells
            for key in degsh:
                if isinstance(degsh, dict):
                    v, C = degsh[key]
                else:
                    # for backward compatibility, allow degsh to be a list
                    if blockgf:
                        v = np.eye(*ss.target_shape)
                    else:
                        v = np.eye(*ss.shape)
                    C = False
                if blockgf and C:
                    gf_to_symm[key].from_L_G_R(v, ss.transpose().copy(), v.conjugate().transpose())
                elif blockgf and not C:
                    gf_to_symm[key].from_L_G_R(v, ss, v.conjugate().transpose())
                elif not blockgf and C:
                    gf_to_symm[key] = np.dot(v, np.dot(ss.transpose().copy(), v.conjugate().transpose()))
                elif not blockgf and not C:
                    gf_to_symm[key] = np.dot(v, np.dot(ss, v.conjugate().transpose()))

def embedding_driver(general_params, solver_params, gw_params, advanced_params):
    """
    Function to run the gw embedding cycle.

    Parameters
    ----------
    general_params : dict
        general parameters as a dict
    solver_params : dict
        solver parameters as a dict
    gw_params : dict
        dft parameters as a dict
    advanced_params : dict
        advanced parameters as a dict
    """

    assert gw_params['code'] == 'aimbes', 'Only AIMBES is currently supported as gw code'

    # prepare output h5 archive
    archive = general_params['jobname']+'/'+general_params['seedname']+'.h5'
    if mpi.is_master_node():
        with HDFArchive(archive, 'a') as ar:
            if 'DMFT_results' not in ar:
                ar.create_group('DMFT_results')
            if 'last_iter' not in ar['DMFT_results']:
                ar['DMFT_results'].create_group('last_iter')
            if 'DMFT_input' not in ar:
                ar.create_group('DMFT_input')
                ar['DMFT_input']['program'] = 'solid_dmft'
                ar['DMFT_input'].create_group('solver')
                ar['DMFT_input'].create_group('version')
                ar['DMFT_input']['version']['triqs_hash'] = triqs_hash
                ar['DMFT_input']['version']['triqs_version'] = triqs_version
                ar['DMFT_input']['version']['solid_dmft_hash'] = solid_dmft_hash
                ar['DMFT_input']['version']['solid_dmft_version'] = solid_dmft_version

    # make sure each iteration is saved to h5 file
    general_params['h5_save_freq'] = 1

    # lad GW input from h5 file
    if mpi.is_master_node():
        gw_data, ir_kernel = convert_gw_output(
            general_params['jobname'] + '/' + general_params['seedname'] + '.h5',
            gw_params['h5_file'],
            general_params['dlr_wmax'], general_params['dlr_eps'],
            it_1e = gw_params['it_1e'],
            it_2e = gw_params['it_2e'],
            delta_calc_type = gw_params['delta_calc_type'],
            delta_bath_fit = gw_params['delta_bath_fit']
        )
        gw_params.update(gw_data)
    mpi.barrier()
    gw_params = mpi.bcast(gw_params)
    iteration = gw_params['it_1e']

    # check if general_params[n_iw] is consistent with dlr meshes
    mesh_f_idx = np.array([iw.index for iw in gw_params['mesh_dlr_iw_f']])
    mesh_b_idx = np.array([iw.index for iw in gw_params['mesh_dlr_iw_b']])
    max_idx = max(abs(mesh_f_idx[0]), abs(mesh_f_idx[-1]), abs(mesh_b_idx[0]), abs(mesh_b_idx[-1]))
    if max_idx > general_params['n_iw']:
        mpi.report(f"\n!!! WARNING !!!!\n"
                   f"general_params['n_iw'] = {general_params['n_iw']} < maximum DLRImFreq index ({max_idx}). "
                   f"solid_dmft will automatically set general_params['n_iw'] = {max_idx+1}.\n")
        general_params['n_iw'] = max_idx + 1

    general_params['beta'] = gw_params['beta']
    sumk_mesh = MeshImFreq(beta=general_params['beta'], statistic='Fermion', n_iw=general_params['n_iw'])

    # if GW calculation was performed with spin never average spin channels
    if gw_params['number_of_spins'] == 2:
        general_params['magnetic'] = True

    sumk = SumkDFT(hdf_file=general_params['jobname'] + '/' + general_params['seedname'] + '.h5',
                   mesh=sumk_mesh, use_dft_blocks=False, h_field=general_params['h_field'])
    sumk.chemical_potential = gw_params['mu_emb']
    sumk.dc_imp = gw_params['Vhf_dc']

    for key in ['h_int_type', 'dc_type', 'enforce_off_diag']:
        general_params = _extract_quantity_per_inequiv(key, sumk.n_inequiv_shells, general_params)
    for key in ['map_solver_struct', 'pick_solver_struct', 'mapped_solver_struct_degeneracies']:
        advanced_params = _extract_quantity_per_inequiv(key, sumk.n_inequiv_shells, advanced_params)

    if len(solver_params) == 1 and solver_params[0]['idx_impurities'] is None:
        map_imp_solver = [0] * sumk.n_inequiv_shells
    else:
        all_idx_imp = [i for entry in solver_params for i in entry['idx_impurities']]
        if sorted(all_idx_imp) != list(range(sumk.n_inequiv_shells)):
            raise ValueError('All impurities must be listed exactly once in solver.idx_impurities'
                             f'but instead got {all_idx_imp}')

        map_imp_solver = []
        for iineq in range(sumk.n_inequiv_shells):
            for isolver, entry in enumerate(solver_params):
                if iineq in entry['idx_impurities']:
                    map_imp_solver.append(isolver)
                    break
    solver_type_per_imp = [solver_params[map_imp_solver[iineq]]['type'] for iineq in range(sumk.n_inequiv_shells)]
    mpi.report(f'\nSolver type per impurity: {solver_type_per_imp}')

    # determine block_structure
    dens_mat_gw = [None] * sumk.n_inequiv_shells
    for ish in range(sumk.n_inequiv_shells):
        dens_mat_gw[ish] = gw_params['Gloc_dlr'][ish].density()
    sumk, _ = _determine_block_structure(sumk, general_params, advanced_params, solver_type_per_imp, dens_mat_gw)

    # print block structure and rotation matrix based on the DFT input
    if mpi.is_master_node():
        print('\n number of ineq. correlated shells: {}'.format(sumk.n_inequiv_shells))
        print('\n block structure summary')
        for icrsh in range(sumk.n_inequiv_shells):
            shlst = [ish for ish, ineq_shell in enumerate(sumk.corr_to_inequiv) if ineq_shell == icrsh]
            print(' -- Shell type #{:3d}: '.format(icrsh) + format(shlst))
            print('  | shell multiplicity ' + str(len(shlst)))
            print('  | block struct. sum_k : ' + format(sumk.gf_struct_sumk[sumk.inequiv_to_corr[icrsh]]))
            print('  | block struct. solver: ' + format(sumk.gf_struct_solver[icrsh]))
            print('  | deg. orbitals : ' + format(sumk.deg_shells[icrsh]))

        print('\nRotation matrices')
        formatter.print_rotation_matrix(sumk)
    mpi.barrier()

    # write block structure to h5 archive
    if gw_params['it_1e'] == 1 and mpi.is_master_node():
        with HDFArchive(archive, 'a') as ar:
            if 'block_structure' not in ar['DMFT_input']:
                ar['DMFT_input']['block_structure'] = sumk.block_structure
            if 'deg_shells' not in ar['DMFT_input']:
                ar['DMFT_input']['deg_shells'] = sumk.deg_shells

    # create h_int
    h_int, gw_params = interaction_hamiltonian.construct(sumk, general_params, advanced_params, gw_params)

    # create solver objects
    solvers = [None] * sumk.n_inequiv_shells
    if mpi.is_master_node():
        Sigma_dlr = [None] * sumk.n_inequiv_shells
        Sigma_dlr_iw = [None] * sumk.n_inequiv_shells
        ir_mesh_idx = ir_kernel.wn_mesh(stats='f',ir_notation=False)
        ir_mesh = (2*ir_mesh_idx+1)*np.pi/gw_params['beta']
        Sigma_ir = np.zeros((len(ir_mesh_idx),
                             gw_params['number_of_spins'],
                             sumk.n_inequiv_shells,max(gw_params['n_orb']),max(gw_params['n_orb'])),
                            dtype=complex)
        Vhf_imp_sIab = np.zeros((gw_params['number_of_spins'],
                                 sumk.n_inequiv_shells,
                                 max(gw_params['n_orb']),max(gw_params['n_orb'])),dtype=complex)

    for icrsh in range(sumk.n_inequiv_shells):
        # Construct the Solver instances
        mpi.report("Creating solver with new interface")
        solvers[icrsh] = create_solver(general_params=general_params,
                                       solver_params=solver_params[map_imp_solver[icrsh]],
                                       gw_params=gw_params,
                                       advanced_params=advanced_params,
                                       sum_k=sumk,
                                       h_int=h_int[icrsh],
                                       icrsh=icrsh,
                                       iteration_offset=iteration,
                                       deg_orbs_ftps=None)

    # init local density matrices for observables
    density_tot = 0.0
    density_shell = np.zeros(sumk.n_inequiv_shells)
    density_mat = [None] * sumk.n_inequiv_shells
    density_mat_unsym = [None] * sumk.n_inequiv_shells
    density_shell_pre = np.zeros(sumk.n_inequiv_shells)
    density_mat_pre = [None] * sumk.n_inequiv_shells

    if sumk.SO:
        printed = ((np.real, 'real'), (np.imag, 'imaginary'))
    else:
        printed = ((np.real, 'real'),)

    for ish in range(sumk.n_inequiv_shells):
        density_shell_pre[ish] = np.real(gw_params['Gloc_dlr'][ish].total_density())
        mpi.report(
            '\n *** Correlated Shell type #{:3d} : '.format(ish)
            + 'Estimated total charge of impurity problem = {:.6f}'.format(density_shell_pre[ish])
        )
        density_mat_pre[ish] = gw_params['Gloc_dlr'][ish].density()
        mpi.report('Estimated density matrix:')
        for key, value in sorted(density_mat_pre[ish].items()):
            for func, name in printed:
                mpi.report('{}, {} part'.format(key, name))
                mpi.report(func(value))

        if not general_params['enforce_off_diag']:
            mpi.report('\n*** WARNING: off-diagonal elements are neglected in the impurity solver ***')

        # convert G0 to solver basis
        G0_dlr = sumk.block_structure.convert_gf(gw_params['G0_dlr'][ish], ish_from=ish, space_from='sumk', space_to='solver')
        # dyson equation to extract G0_freq, using Hermitian symmetry (always needed in solver postprocessing)
        solvers[ish].G0_freq << make_hermitian(make_gf_imfreq(G0_dlr, n_iw=general_params['n_iw']))

        if ((solver_type_per_imp[ish] == 'cthyb' and solvers[ish].solver_params['delta_interface'])
                or solver_type_per_imp[ish] == 'ctseg'):
            mpi.report('\n Using the delta interface for passing Delta(tau) and Hloc0 directly to the solver.\n')

            # prepare solver input
            imp_eal = sumk.block_structure.convert_matrix(gw_params['Hloc0'][ish], ish_from=ish, space_from='sumk', space_to='solver')
            Delta_dlr = sumk.block_structure.convert_gf(gw_params['delta_dlr'][ish], ish_from=ish, space_from='sumk', space_to='solver')

            for name, delta in Delta_dlr:
                if mpi.is_master_node():
                    print('H_loc0[{:2d}] block: {}'.format(ish, name))
                    fmt = '{:11.7f}' * imp_eal[name].shape[0]
                    for block in imp_eal[name]:
                        print((' '*11 + fmt).format(*block.real))

                # create now full delta(tau)
                Delta_tau = make_hermitian(make_gf_imtime(delta, n_tau=general_params['n_tau']))

                # without SOC delta_tau needs to be real
                if not sumk.SO == 1:
                    solvers[ish].Delta_time[name] << Delta_tau.real
                else:
                    solvers[ish].Delta_time[name] << Delta_tau

                if solvers[ish].solver_params['diag_delta']:
                    for o1 in range(imp_eal[name].shape[0]):
                        for o2 in range(imp_eal[name].shape[0]):
                            if o1 != o2:
                                solvers[ish].Delta_time[name].data[:, o1, o2] = 0.0 + 0.0j

            # Make non-interacting operator for Hloc0
            Hloc_0 = Operator()
            if solver_type_per_imp[ish] == 'ctseg':
                mpi.report('Discarding off-diagonals in Hloc_0 anyway since ctseg enforces it.')
                for spin, spin_block in imp_eal.items():
                    for o1 in range(spin_block.shape[0]):
                        Hloc_0 += spin_block[o1, o1].real * c_dag(spin, o1) * c(spin, o1)
            else:
                for spin, spin_block in imp_eal.items():
                    for o1 in range(spin_block.shape[0]):
                        for o2 in range(spin_block.shape[1]):
                            # check if off-diag element is larger than threshold
                            if o1 != o2 and abs(spin_block[o1, o2]) < solver_params[ish]['off_diag_threshold']:
                                continue
                            else:
                                # TODO: adapt for SOC calculations, which should keep the imag part
                                Hloc_0 += spin_block[o1, o2].real / 2 * (c_dag(spin, o1) * c(spin, o2) + c_dag(spin, o2) * c(spin, o1))
            solvers[ish].Hloc_0 = Hloc_0

        mpi.report('\nSolving the impurity problem for shell {} ...'.format(ish))
        mpi.barrier()
        start_time = timer()
        solvers[ish].solve()
        mpi.barrier()
        mpi.report('Actual time for solver: {:.2f} s'.format(timer() - start_time))

        # some printout of the obtained density matrices and some basic checks from the unsymmetrized solver output
        if solvers[ish].solver_params['type'] == 'ctseg':
            for block, occ_mat in solvers[ish].orbital_occupations.items():
                density_shell[ish] += np.trace(occ_mat)
            density_tot += density_shell[ish]
            density_mat_unsym[ish] = solvers[ish].orbital_occupations
            density_mat[ish] = density_mat_unsym[ish]
        else:
            density_shell[ish] = np.real(solvers[ish].G_freq_unsym.total_density())
            density_tot += density_shell[ish]
            density_mat_unsym[ish] = solvers[ish].G_freq_unsym.density()
            density_mat[ish] = solvers[ish].G_freq.density()
        formatter.print_local_density(density_shell[ish], density_shell_pre[ish], density_mat_unsym[ish], sumk.SO)
        mpi.report('')

        # post-processing for GW
        if mpi.is_master_node():
            if not hasattr(solvers[ish], 'Sigma_Hartree'):
                print('Moments of Sigma not measured using tail fit to extract static Hartree shift for DLR fit.')
                solvers[ish].Sigma_Hartree = {}
                for block, gf in solvers[ish].Sigma_freq:
                    tail, err = gf.fit_hermitian_tail()
                    solvers[ish].Sigma_Hartree[block] = tail[0]

            if solvers[ish].solver_params['type'] in ('cthyb', 'ctseg') and solvers[ish].solver_params['crm_dyson_solver']:
                Sigma_dlr[ish] = make_gf_dlr(solvers[ish].Sigma_dlr)
            else:
                Sigma_dlr_iw[ish] = sumk.block_structure.create_gf(ish=ish,
                                                                   gf_function=Gf,
                                                                   space='solver',
                                                                   mesh=gw_params['mesh_dlr_iw_f'])
                for w in Sigma_dlr_iw[ish].mesh:
                    for block, gf in Sigma_dlr_iw[ish]:
                        gf[w] = solvers[ish].Sigma_freq[block](w)-solvers[ish].Sigma_Hartree[block]

                sumk.symm_deg_gf(Sigma_dlr_iw[ish],ish=ish)
                Sigma_dlr[ish] = make_gf_dlr(Sigma_dlr_iw[ish])

            # mixing of impurity Sigma
            if general_params['sigma_mix'] < 1.0 and iteration > 1:
                if mpi.is_master_node():
                    print('mixing sigma with previous iteration by factor {:.3f}\n'.format(general_params['sigma_mix']))
                    with HDFArchive(general_params['jobname'] + '/' + general_params['seedname'] + '.h5', 'r') as ar:
                        Sigma_dlr_prev = ar[f'DMFT_results/it_{iteration-1}'][f'Sigma_dlr_{ish}']
                        Sigma_Hartree_prev = ar[f'DMFT_results/it_{iteration-1}'][f'Sigma_Hartree_{ish}']
                    Sigma_dlr[ish] << (general_params['sigma_mix'] * Sigma_dlr[ish]
                                                + (1-general_params['sigma_mix']) * Sigma_dlr_prev)
                    for block in solvers[ish].Sigma_Hartree.keys():
                        solvers[ish].Sigma_Hartree[block] = (general_params['sigma_mix'] * solvers[ish].Sigma_Hartree[block]
                                                    + (1-general_params['sigma_mix']) * Sigma_Hartree_prev[block])
                Sigma_dlr[ish] = mpi.bcast(Sigma_dlr[ish])
                solvers[ish].Sigma_Hartree = mpi.bcast(solvers[ish].Sigma_Hartree)


            for i, (block, gf) in enumerate(Sigma_dlr[ish]):
                # print Hartree shift
                print('Σ_HF {}'.format(block))
                fmt = '{:11.7f}' * solvers[ish].Sigma_Hartree[block].shape[0]
                for vhf in solvers[ish].Sigma_Hartree[block]:
                    print((' '*11 + fmt).format(*vhf.real))

            # average Hartree shift if not magnetic
            if not general_params['magnetic']:
                if general_params['enforce_off_diag']:
                    solvers[ish].Sigma_Hartree['up_0'] = 0.5*(solvers[ish].Sigma_Hartree['up_0']+
                                                              solvers[ish].Sigma_Hartree['down_0'])
                    solvers[ish].Sigma_Hartree['down_0'] = solvers[ish].Sigma_Hartree['up_0']
                else:
                    for iorb in range(gw_params['n_orb'][ish]):
                        solvers[ish].Sigma_Hartree[f'up_{iorb}'] = 0.5*(solvers[ish].Sigma_Hartree[f'up_{iorb}']+
                                                                      solvers[ish].Sigma_Hartree[f'down_{iorb}'])
                        solvers[ish].Sigma_Hartree[f'down_{iorb}'] = solvers[ish].Sigma_Hartree[f'up_{iorb}']

            iw_mesh = solvers[ish].Sigma_freq.mesh
            # convert Sigma to sumk basis
            Sigma_dlr_sumk = sumk.block_structure.convert_gf(Sigma_dlr[ish], ish_from=ish, space_from='solver', space_to='sumk')
            Sigma_Hartree_sumk = sumk.block_structure.convert_matrix(solvers[ish].Sigma_Hartree, ish_from=ish, space_from='solver', space_to='sumk')
            # store Sigma and V_HF in sumk basis on IR mesh
            ir_nw_half = len(ir_mesh_idx)//2
            for i, (block, gf) in enumerate(Sigma_dlr_sumk):
                Vhf_imp_sIab[i,ish] = Sigma_Hartree_sumk[block]
                # Make sure sigma_ir[iw].conj() = sigma_ir[-iw]
                for iw_idx in range(ir_nw_half):
                    iw_pos = ir_nw_half+iw_idx
                    iw_neg = ir_nw_half-1-iw_idx
                    Sigma_ir[iw_pos,i,ish] = gf(iw_mesh(ir_mesh_idx[iw_pos]))
                    Sigma_ir[iw_neg,i,ish] = gf(iw_mesh(ir_mesh_idx[iw_pos])).conj()

                if not general_params['magnetic']:
                    break

    if mpi.is_master_node():
        print("\nChecking impurity self-energy on the IR mesh...")
        ir_kernel.check_leakage(Sigma_ir, stats='f', name="impurity self-energy", w_input=True)

    # Writes results to h5 archive
    mpi.report('Writing iter {} results to h5 archives {} and {}.'.format(iteration, archive, gw_params['h5_file']))
    if mpi.is_master_node():
        results_to_archive.write(archive, sumk, general_params, solver_params, solvers,
                                 map_imp_solver, solver_type_per_imp, iteration,
                                 False, gw_params['mu_emb'], density_mat_pre, density_mat)
        with HDFArchive(archive, 'a') as ar:
            # store also IR / DLR Sigma
            ar['DMFT_results/it_{}'.format(iteration)]['ir_mesh'] = ir_mesh
            ar['DMFT_results/it_{}'.format(iteration)]['Sigma_imp_wsIab'] = Sigma_ir
            ar['DMFT_results/it_{}'.format(iteration)]['Vhf_imp_sIab'] = Vhf_imp_sIab
            for ish in range(sumk.n_inequiv_shells):
                ar['DMFT_results/it_{}'.format(iteration)][f'Sigma_dlr_{ish}'] = Sigma_dlr[ish]

        # write results to GW h5_file
        with HDFArchive(gw_params['h5_file'], 'a') as ar:
            ar[f'downfold_1e/iter{iteration}']['Sigma_imp_wsIab'] = Sigma_ir
            ar[f'downfold_1e/iter{iteration}']['Vhf_imp_sIab'] = Vhf_imp_sIab

    mpi.report('*** iteration finished ***')
    mpi.report('#'*80)
    mpi.barrier()
    return
