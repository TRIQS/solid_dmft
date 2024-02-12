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
from triqs.gf.tools import inverse
from triqs.gf import (
    Gf,
    make_hermitian,
    make_gf_dlr,
    make_gf_imfreq,
    make_gf_imtime,
    make_gf_dlr_imfreq,
)
from triqs.version import git_hash as triqs_hash
from triqs.version import version as triqs_version
from triqs.gf.meshes import MeshImFreq
from triqs.operators import c_dag, c, Operator
from triqs_dft_tools.block_structure import BlockStructure

from solid_dmft.version import solid_dmft_hash
from solid_dmft.version import version as solid_dmft_version
from solid_dmft.dmft_tools import formatter
from solid_dmft.dmft_tools import results_to_archive
from solid_dmft.dmft_tools.solver import SolverStructure
from solid_dmft.dmft_tools import interaction_hamiltonian
from solid_dmft.gw_embedding.bdft_converter import convert_gw_output


class dummy_sumk(object):
    """
    create dummy sumk helper object
    """

    def __init__(self, n_inequiv_shells, n_orb_list, use_rot, magnetic):
        self.n_inequiv_shells = n_inequiv_shells
        self.SO = 0
        self.use_rotations = use_rot
        if self.use_rotations:
            raise ValueError('rotations not implemented yet for GW embedding')
        self.gf_struct_solver = []
        self.gf_struct_sumk = []
        self.spin_block_names = []
        self.inequiv_to_corr = []
        self.deg_shells = []
        self.dc_energ = [0.0 for ish in range(self.n_inequiv_shells)]
        self.sumk_to_solver = [{} for ish in range(self.n_inequiv_shells)]
        self.solver_to_sumk = [{} for ish in range(self.n_inequiv_shells)]
        self.solver_to_sumk_block = [{} for ish in range(self.n_inequiv_shells)]
        for ish in range(self.n_inequiv_shells):
            self.gf_struct_solver.append({'up_0': n_orb_list[ish], 'down_0': n_orb_list[ish]})
            self.gf_struct_sumk.append({'up': n_orb_list[ish], 'down': n_orb_list[ish]})
            self.spin_block_names.append(['up', 'down'])
            self.inequiv_to_corr.append(ish)
            if not magnetic:
                self.deg_shells.append([['up_0', 'down_0']])
            # setup standard mapping between sumk and solver
            for block, inner_dim in self.gf_struct_sumk[ish].items():
                self.solver_to_sumk_block[ish][block] = block
                for inner in range(inner_dim):
                    self.sumk_to_solver[ish][(block, inner)] = (block + '_0', inner)
                    self.solver_to_sumk[ish][(block + '_0', inner)] = (block, inner)

        self.gf_struct_solver_list = [sorted([(k, v) for k, v in list(gfs.items())], key=lambda x: x[0]) for gfs in self.gf_struct_solver]

        # creat block_structure object for solver
        self.block_structure = BlockStructure(
            gf_struct_sumk=self.gf_struct_sumk,
            gf_struct_solver=self.gf_struct_solver,
            solver_to_sumk=self.solver_to_sumk,
            sumk_to_solver=self.sumk_to_solver,
            solver_to_sumk_block=self.solver_to_sumk_block,
            deg_shells=self.deg_shells,
            transformation=None,
        )

    def symm_deg_gf(self, gf_to_symm, ish=0):
        r"""
        Averages a GF over degenerate shells.

        Degenerate shells of an inequivalent correlated shell are defined by
        `self.deg_shells`. This function enforces corresponding degeneracies
        in the input GF.

        Parameters
        ----------
        gf_to_symm : gf_struct_solver like
                     Input and output GF (i.e., it gets overwritten)
        ish : int
              Index of an inequivalent shell. (default value 0)

        """

        # when reading block_structures written with older versions from
        # an h5 file, self.deg_shells might be None
        if self.deg_shells is None:
            return

        for degsh in self.deg_shells[ish]:
            # ss will hold the averaged orbitals in the basis where the
            # blocks are all equal
            # i.e. maybe_conjugate(v^dagger gf v)
            ss = None
            n_deg = len(degsh)
            for key in degsh:
                if ss is None:
                    ss = gf_to_symm[key].copy()
                    ss.zero()
                    helper = ss.copy()
                # get the transformation matrix
                if isinstance(degsh, dict):
                    v, C = degsh[key]
                else:
                    # for backward compatibility, allow degsh to be a list
                    v = np.eye(*ss.target_shape)
                    C = False
                # the helper is in the basis where the blocks are all equal
                helper.from_L_G_R(v.conjugate().transpose(), gf_to_symm[key], v)
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
                    v = np.eye(*ss.target_shape)
                    C = False
                if C:
                    gf_to_symm[key].from_L_G_R(v, ss.transpose().copy(), v.conjugate().transpose())
                else:
                    gf_to_symm[key].from_L_G_R(v, ss, v.conjugate().transpose())


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

    iteration = 1
    # prepare output h5 archive
    if mpi.is_master_node():
        with HDFArchive(general_params['jobname'] + '/' + general_params['seedname'] + '.h5', 'a') as ar:
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

            if 'iteration_count' in ar['DMFT_results']:
                iteration = ar['DMFT_results/iteration_count'] + 1
            else:
                iteration = 1
    iteration = mpi.bcast(iteration)

    # lad GW input from h5 file
    if mpi.is_master_node():
        gw_data, ir_kernel = convert_gw_output(
            general_params['jobname'] + '/' + general_params['seedname'] + '.h5',
            gw_params['h5_file'],
            wmax_dlr=general_params['w_max'],
        )
        gw_params.update(gw_data)
    mpi.barrier()
    gw_params = mpi.bcast(gw_params)

    # if GW calculation was performed with spin never average spin channels
    if gw_params['number_of_spins'] == 2:
        general_params['magnetic'] = True

    # dummy helper class for sumk
    sumk = dummy_sumk(gw_params['n_inequiv_shells'], gw_params['n_orb'], gw_params['use_rot'], general_params['magnetic'])
    sumk.mesh = MeshImFreq(beta=gw_params['beta'], statistic='Fermion', n_iw=general_params['n_iw'])
    sumk.chemical_potential = gw_params['mu_emb']
    sumk.dc_imp = gw_params['Vhf_dc']
    general_params['beta'] = gw_params['beta']

    # create h_int
    h_int = interaction_hamiltonian.construct(sumk, general_params, advanced_params, gw_params)

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
    for ish in range(sumk.n_inequiv_shells):
        # Construct the Solver instances
        solvers[ish] = SolverStructure(general_params, solver_params, advanced_params, sumk, ish, h_int[ish])

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

        if general_params['solver_type'] == 'cthyb' and general_params['cthyb_delta_interface']:
            mpi.report('\n Using the delta interface for cthyb passing Delta(tau) and Hloc0 directly.\n')
            # prepare solver input
            imp_eal = gw_params['Hloc0'][ish]
            # fill Delta_time from Delta_freq sumk to solver
            for name, g0 in gw_params['G0_dlr'][ish]:
                G0_dlr_iw = make_gf_dlr_imfreq(g0)
                # make non-interacting impurity Hamiltonian hermitian
                imp_eal[name] = (imp_eal[name] + imp_eal[name].T.conj())/2
                if mpi.is_master_node():
                    print('H_loc0[{:2d}] block: {}'.format(ish, name))
                    fmt = '{:11.7f}' * imp_eal[name].shape[0]
                    for block in imp_eal[name]:
                        print((' '*11 + fmt).format(*block.real))

                Delta_dlr_iw = Gf(mesh=G0_dlr_iw.mesh, target_shape=g0.target_shape)
                for iw in G0_dlr_iw.mesh:
                    Delta_dlr_iw[iw] = iw.value - inverse(G0_dlr_iw[iw]) - imp_eal[name]

                # without SOC delta_tau needs to be real
                if not sumk.SO == 1:
                    Delta_dlr = make_gf_dlr(Delta_dlr_iw)
                    Delta_tau = make_hermitian(make_gf_imtime(Delta_dlr, n_tau=general_params['n_tau']))
                    # create now full delta(tau)
                    solvers[ish].Delta_time[name] << Delta_tau.real
                else:
                    solvers[ish].Delta_time[name] << Delta_tau

                if general_params['diag_delta']:
                    for o1 in range(imp_eal[name].shape[0]):
                        for o2 in range(imp_eal[name].shape[0]):
                            if o1 != o2:
                                solvers[ish].Delta_time[name].data[:, o1, o2] = 0.0 + 0.0j

            # Make non-interacting operator for Hloc0
            Hloc_0 = Operator()
            for spin, spin_block in imp_eal.items():
                for o1 in range(spin_block.shape[0]):
                    for o2 in range(spin_block.shape[1]):
                        # check if off-diag element is larger than threshold
                        if o1 != o2 and abs(spin_block[o1, o2]) < solver_params['off_diag_threshold']:
                            continue
                        else:
                            # TODO: adapt for SOC calculations, which should keep the imag part
                            Hloc_0 += spin_block[o1, o2].real / 2 * (c_dag(spin, o1) * c(spin, o2) + c_dag(spin, o2) * c(spin, o1))
            solvers[ish].Hloc_0 = Hloc_0
        else:
            # dyson equation to extract G0_freq, using Hermitian symmetry
            solvers[ish].G0_freq << make_hermitian(make_gf_imfreq(gw_params['G0_dlr'][ish], n_iw=general_params['n_iw']))

        mpi.report('\nSolving the impurity problem for shell {} ...'.format(ish))
        mpi.barrier()
        start_time = timer()
        solvers[ish].solve()
        mpi.barrier()
        mpi.report('Actual time for solver: {:.2f} s'.format(timer() - start_time))

        # some printout of the obtained density matrices and some basic checks from the unsymmetrized solver output
        density_shell[ish] = np.real(solvers[ish].G_freq_unsym.total_density())
        density_tot += density_shell[ish]
        density_mat_unsym[ish] = solvers[ish].G_freq_unsym.density()
        density_mat[ish] = solvers[ish].G_freq.density()
        formatter.print_local_density(density_shell[ish], density_shell_pre[ish], density_mat_unsym[ish], sumk.SO)
        mpi.report('')

        # post-processing for GW
        if mpi.is_master_node():
            Sigma_dlr_iw[ish] = sumk.block_structure.create_gf(ish=ish,
                                                               gf_function=Gf,
                                                               space='solver',
                                                               mesh=gw_params['mesh_dlr_iw_f'])
            for w in Sigma_dlr_iw[ish].mesh:
                for block, gf in Sigma_dlr_iw[ish]:
                    gf[w] = solvers[ish].Sigma_freq[block](w)-solvers[ish].Sigma_Hartree[block]

            sumk.symm_deg_gf(Sigma_dlr_iw[ish],ish=ish)
            Sigma_dlr[ish] = make_gf_dlr(Sigma_dlr_iw[ish])

            iw_mesh = solvers[ish].Sigma_freq.mesh
            for block, gf in Sigma_dlr[ish]:
                # print Hartree shift
                print('Σ_HF {}'.format(block))
                fmt = '{:11.7f}' * solvers[ish].Sigma_Hartree[block].shape[0]
                for vhf in solvers[ish].Sigma_Hartree[block]:
                    print((' '*11 + fmt).format(*vhf.real))

                # if GW is performed without spin we just write one spin channel back (averaged)
                if 'up' in block:
                    i = 0
                elif 'down' in block and gw_params['number_of_spins'] == 2:
                    i = 1

                Vhf_imp_sIab[i,ish] = solvers[ish].Sigma_Hartree[block]
                for iw in range(len(ir_mesh_idx)):
                    Sigma_ir[iw,i,ish] = gf(iw_mesh(ir_mesh_idx[iw]))


            # average hartree shift for storing to be in line with averaged Sigma imp
            Vhf_imp_sIab[:,ish] = np.mean(Vhf_imp_sIab[:,ish],axis=0)


    # Writes results to h5 archive
    if mpi.is_master_node():
        with HDFArchive(general_params['jobname'] + '/' + general_params['seedname'] + '.h5', 'a') as ar:
            results_to_archive.write(ar, sumk, general_params, solver_params, solvers, iteration,
                             False, gw_params['mu_emb'], density_mat_pre, density_mat)

            # store also IR / DLR Sigma
            ar['DMFT_results/it_{}'.format(iteration)]['ir_mesh'] = ir_mesh
            ar['DMFT_results/it_{}'.format(iteration)]['Sigma_imp_wsIab'] = Sigma_ir
            ar['DMFT_results/it_{}'.format(iteration)]['Vhf_imp_sIab'] = Vhf_imp_sIab
            for ish in range(sumk.n_inequiv_shells):
                ar['DMFT_results/it_{}'.format(iteration)][f'Sigma_dlr_{ish}'] = Sigma_dlr[ish]

        # write results to GW h5_file
        with HDFArchive(gw_params['h5_file'],'a') as ar:
            ar[f'downfold_1e/iter{iteration}']['Sigma_imp_wsIab'] = Sigma_ir
            ar[f'downfold_1e/iter{iteration}']['Vhf_imp_wsIab'] = Vhf_imp_sIab


    mpi.report('*** iteration finished ***')
    mpi.report('#'*80)
    mpi.barrier()
    return
