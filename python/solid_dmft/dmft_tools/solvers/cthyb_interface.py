# %%
################################################################################
#
# solid_dmft - A versatile python wrapper to perform DFT+DMFT calculations
#              utilizing the TRIQS software library
#
# Copyright (C) 2018-2020, ETH Zurich
# Copyright (C) 2021, The Simons Foundation
#      authors: A. Carta, A. Hampel, M. Merkel, and S. Beck
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
'''
cthyb solver class for solid_dmft
'''
import numpy as np
from triqs.gfs import MeshDLRImFreq, Gf, make_gf_imfreq, make_hermitian, fit_gf_dlr, make_gf_dlr_imtime
from triqs.gfs.tools import inverse
import triqs.utility.mpi as mpi
from h5 import HDFArchive

from solid_dmft.io_tools.dict_to_h5 import prep_params_for_h5

from solid_dmft.dmft_tools import legendre_filter
from solid_dmft.dmft_tools.matheval import MathExpr


# import of the abstract class
from solid_dmft.dmft_tools.solvers.abstractdmftsolver import AbstractDMFTSolver

# import triqs solver
from triqs_cthyb.solver import Solver as cthyb_solver
from triqs_cthyb.version import triqs_cthyb_hash, version


class CTHYBInterface(AbstractDMFTSolver):
    def __init__(
        self, general_params, solver_params, sum_k, icrsh, h_int, iteration_offset, deg_orbs_ftps, gw_params=None, advanced_params=None
    ):
        # Call the base class constructor
        super().__init__(general_params, solver_params, sum_k, icrsh, h_int, iteration_offset, deg_orbs_ftps, gw_params, advanced_params)

        # sets up necessary GF objects on ImFreq
        self._init_ImFreq_objects()

        ###################################################
        # Create the cthyb solver specifics
        ###################################################

        if self.solver_params.get('random_seed') is None:
            self.random_seed_generator = None
        else:
            self.random_seed_generator = MathExpr(self.solver_params['random_seed'])

        # Separately stores all params that go into solve() call of solver
        self.triqs_solver_params = {}
        keys_to_pass = (
            'det_init_size',
            'det_n_operations_before_check',
            'det_precision_warning',
            'det_precision_error',
            'imag_threshold',
            'length_cycle',
            'max_time',
            'measure_density_matrix',
            'measure_G_l',
            'measure_pert_order',
            'move_double',
            'move_shift',
            'off_diag_threshold',
            'perform_tail_fit',
        )
        for key in keys_to_pass:
            self.triqs_solver_params[key] = self.solver_params[key]

        # Calculates number of sweeps per rank
        self.triqs_solver_params['n_cycles'] = int(self.solver_params['n_cycles_tot'] / mpi.size)
        # cast warmup cycles to int in case given in scientific notation
        self.triqs_solver_params['n_warmup_cycles'] = int(self.solver_params['n_warmup_cycles'])

        # Renames measure chi param
        self.triqs_solver_params['measure_O_tau_min_ins'] = self.solver_params['measure_chi_insertions']

        # use_norm_as_weight also required to measure the density matrix
        self.triqs_solver_params['use_norm_as_weight'] = self.triqs_solver_params['measure_density_matrix']

        if self.triqs_solver_params['perform_tail_fit']:
            for key in ('fit_max_moment', 'fit_max_n', 'fit_max_w', 'fit_min_n', 'fit_min_w'):
                self.triqs_solver_params[key] = self.solver_params[key]

        # set loc_n_min and loc_n_max
        if self.solver_params['loc_n_min'] is not None:
            self.triqs_solver_params['loc_n_min'] = self.solver_params['loc_n_min']
        if self.solver_params['loc_n_max'] is not None:
            self.triqs_solver_params['loc_n_max'] = self.solver_params['loc_n_max']

        gf_struct = self.sum_k.gf_struct_solver_list[self.icrsh]
        # Construct the triqs_solver instances
        if self.solver_params['measure_G_l']:
            self.triqs_solver = cthyb_solver(
                beta=self.general_params['beta'],
                gf_struct=gf_struct,
                n_iw=self.general_params['n_iw'],
                n_tau=self.general_params['n_tau'],
                n_l=self.solver_params['n_l'],
                delta_interface=self.solver_params['delta_interface'],
            )
        else:
            self.triqs_solver = cthyb_solver(
                beta=self.general_params['beta'],
                gf_struct=gf_struct,
                n_iw=self.general_params['n_iw'],
                n_tau=self.general_params['n_tau'],
                delta_interface=self.solver_params['delta_interface'],
            )

        # sets up metadata
        self.git_hash = triqs_cthyb_hash
        self.version = version

        return

    def solve(self, **kwargs):
        if self.random_seed_generator is not None:
            self.triqs_solver_params['random_seed'] = int(self.random_seed_generator(it=kwargs['it'], rank=mpi.rank))
        else:
            assert 'random_seed' not in self.triqs_solver_params

        # pass measure_O_tau which is prepared late in dmft cycle and only ready now
        if self.solver_params['measure_chi'] is not None:
            self.triqs_solver_params['measure_O_tau'] = self.solver_params['measure_O_tau']

        if self.solver_params['delta_interface']:
            self.triqs_solver.Delta_tau << self.Delta_time
            self.triqs_solver_params['h_loc0'] = self.Hloc_0
        else:
            # fill G0_freq from sum_k to solver
            self.triqs_solver.G0_iw << make_hermitian(self.G0_freq)

        # update solver in h5 archive one last time for debugging if solve command crashes
        if self.general_params['store_solver'] and mpi.is_master_node():
            with HDFArchive(self.general_params['jobname'] + '/' + self.general_params['seedname'] + '.h5', 'a') as archive:
                if not 'it_-1' in archive['DMFT_input/solver']:
                    archive['DMFT_input/solver'].create_group('it_-1')
                archive['DMFT_input/solver/it_-1'][f'S_{self.icrsh}'] = self.triqs_solver
                if self.solver_params['delta_interface']:
                    archive['DMFT_input/solver/it_-1'][f'Delta_time_{self.icrsh}'] = self.triqs_solver.Delta_tau
                else:
                    archive['DMFT_input/solver/it_-1'][f'G0_freq_{self.icrsh}'] = self.triqs_solver.G0_iw
                # archive['DMFT_input/solver/it_-1'][f'Delta_freq_{self.icrsh}'] = self.Delta_freq
                archive['DMFT_input/solver/it_-1'][f'solve_params_{self.icrsh}'] = prep_params_for_h5(self.solver_params)
                archive['DMFT_input/solver/it_-1'][f'triqs_solver_params_{self.icrsh}'] = prep_params_for_h5(self.triqs_solver_params)
                archive['DMFT_input/solver/it_-1']['mpi_size'] = mpi.size
        mpi.barrier()

        # Solve the impurity problem for icrsh shell
        # *************************************
        self.triqs_solver.solve(h_int=self.h_int, **self.triqs_solver_params)
        # *************************************

        # dump Delta_tau constructed internally from cthyb when delta_interface = False
        if self.general_params['store_solver'] and mpi.is_master_node():
            with HDFArchive(self.general_params['jobname'] + '/' + self.general_params['seedname'] + '.h5', 'a') as archive:
                if not self.solver_params['delta_interface']:
                    archive['DMFT_input/solver/it_-1'][f'Delta_time_{self.icrsh}'] = self.triqs_solver.Delta_tau
        mpi.barrier()

        # call postprocessing
        self.postprocess()

        return

    def postprocess(self):
        r"""
        Organize G_freq, G_time, Sigma_freq and G_l from cthyb solver
        """

        def set_Gs_from_G_l():
            # create new G_freq and G_time
            for i, g in self.G_l:
                g.enforce_discontinuity(np.identity(g.target_shape[0]))
                # set G_freq from Legendre and Fouriertransform to get G_time
                self.G_freq[i].set_from_legendre(g)
                self.G_time[i].set_from_legendre(g)

            # Symmetrize
            self.G_freq << make_hermitian(self.G_freq)
            self.G_freq_unsym << self.G_freq
            self.sum_k.symm_deg_gf(self.G_freq, ish=self.icrsh)
            self.sum_k.symm_deg_gf(self.G_time, ish=self.icrsh)
            # Dyson equation to get Sigma_freq
            self.Sigma_freq << inverse(self.G0_freq) - inverse(self.G_freq)

            return

        # get Delta_time from solver
        self.Delta_time << self.triqs_solver.Delta_tau

        # if measured in Legendre basis, get G_l from solver too
        if self.solver_params['measure_G_l']:
            # store original G_time into G_time_orig
            self.G_time_orig << self.triqs_solver.G_tau
            self.G_l << self.triqs_solver.G_l
            # get G_time, G_freq, Sigma_freq from G_l
            set_Gs_from_G_l()

        else:
            self.G_freq << make_hermitian(self.triqs_solver.G_iw)
            self.G_freq_unsym << self.G_freq
            self.sum_k.symm_deg_gf(self.G_freq, ish=self.icrsh)
            # set G_time
            self.G_time << self.triqs_solver.G_tau
            self.sum_k.symm_deg_gf(self.G_time, ish=self.icrsh)

            if self.solver_params['legendre_fit']:
                self.G_time_orig << self.triqs_solver.G_tau
                # run the filter
                self.G_l << legendre_filter.apply(self.G_time, self.solver_params['n_l'])
                # get G_time, G_freq, Sigma_freq from G_l
                set_Gs_from_G_l()
            elif self.solver_params['perform_tail_fit'] and not self.solver_params['legendre_fit']:
                # if tailfit has been used replace Sigma with the tail fitted Sigma from cthyb
                self.Sigma_freq << self.triqs_solver.Sigma_iw
            elif self.solver_params['crm_dyson_solver']:
                from triqs.gfs.dlr_crm_dyson_solver import minimize_dyson

                mpi.report('\nCRM Dyson solver to extract Σ impurity\n')
                # fit QMC G_tau to DLR
                if mpi.is_master_node():
                    if self.solver_params['crm_dlr_wmax'] is not None:
                        dlr_wmax = self.solver_params['crm_dlr_wmax']
                    else:
                        dlr_wmax = self.general_params['dlr_wmax']
                    if self.solver_params['crm_dlr_eps'] is not None:
                        dlr_eps = self.solver_params['crm_dlr_eps']
                    else:
                        dlr_eps = self.general_params['dlr_eps']
                    mpi.report(f'crm_dyson_solver with (wmax, eps) = ({dlr_wmax}, {dlr_eps}). ')
                    G_dlr = fit_gf_dlr(self.triqs_solver.G_tau, w_max=dlr_wmax, eps=dlr_eps)
                    self.G_time_dlr = make_gf_dlr_imtime(G_dlr)

                    # assume little error on G0_iw and use to get G0_dlr
                    mesh_dlr_iw = MeshDLRImFreq(G_dlr.mesh)
                    G0_dlr_iw = self.sum_k.block_structure.create_gf(ish=self.icrsh, gf_function=Gf, mesh=mesh_dlr_iw, space='solver')
                    for block, gf in G0_dlr_iw:
                        for iwn in mesh_dlr_iw:
                            gf[iwn] = self.G0_freq[block](iwn)
                    self.sum_k.symm_deg_gf(G0_dlr_iw, ish=self.icrsh)

                    # average moments
                    self.sum_k.symm_deg_gf(self.triqs_solver.Sigma_Hartree, ish=self.icrsh)
                    first_mom = {}
                    for block, mom in self.triqs_solver.Sigma_moments.items():
                        first_mom[block] = mom[1]
                    self.sum_k.symm_deg_gf(first_mom, ish=self.icrsh)

                    for block, mom in self.triqs_solver.Sigma_moments.items():
                        mom[0] = self.triqs_solver.Sigma_Hartree[block]
                        mom[1] = first_mom[block]

                    # minimize dyson for the first entry of each deg shell
                    self.Sigma_dlr = self.sum_k.block_structure.create_gf(ish=self.icrsh, gf_function=Gf, mesh=mesh_dlr_iw, space='solver')
                    # without any degenerate shells we run the minimization for all blocks
                    if self.sum_k.deg_shells[self.icrsh] == []:
                        for block, gf in self.Sigma_dlr:
                            np.random.seed(85281)
                            print('Minimizing Dyson via CRM for Σ[block]:', block)
                            gf, _, _ = minimize_dyson(
                                G0_dlr=G0_dlr_iw[block], G_dlr=G_dlr[block], Sigma_moments=self.triqs_solver.Sigma_moments[block]
                            )
                    else:
                        for deg_shell in self.sum_k.deg_shells[self.icrsh]:
                            for i, block in enumerate(deg_shell):
                                if i == 0:
                                    np.random.seed(85281)
                                    print('Minimizing Dyson via CRM for Σ[block]:', block)
                                    self.Sigma_dlr[block], _, _ = minimize_dyson(
                                        G0_dlr=G0_dlr_iw[block], G_dlr=G_dlr[block], Sigma_moments=self.triqs_solver.Sigma_moments[block]
                                    )
                                    sol_block = block
                                else:
                                    self.Sigma_dlr[block] << self.Sigma_dlr[sol_block]

                                    print(f'Copying result from first block of deg list to {block}')
                    print('\n')

                    self.Sigma_freq = make_gf_imfreq(self.Sigma_dlr, n_iw=self.general_params['n_iw'])
                    for block, gf in self.Sigma_freq:
                        gf += self.triqs_solver.Sigma_moments[block][0]

                mpi.barrier()
                self.Sigma_freq = mpi.bcast(self.Sigma_freq)
                self.Sigma_dlr = mpi.bcast(self.Sigma_dlr)
                self.G_time = mpi.bcast(self.G_time)
                self.G_time_dlr = mpi.bcast(self.G_time_dlr)
            else:
                # obtain Sigma via dyson from symmetrized G_freq
                self.Sigma_freq << inverse(self.G0_freq) - inverse(self.G_freq)

        # if density matrix is measured, get this too
        if self.solver_params['measure_density_matrix']:
            self.density_matrix = self.triqs_solver.density_matrix
            self.h_loc_diagonalization = self.triqs_solver.h_loc_diagonalization
            self.Sigma_moments = self.triqs_solver.Sigma_moments
            self.Sigma_Hartree = self.triqs_solver.Sigma_Hartree
            self.G_moments = self.triqs_solver.G_moments
            self.orbital_occupations = self.triqs_solver.orbital_occupations

        if self.solver_params['measure_pert_order']:
            self.perturbation_order = self.triqs_solver.perturbation_order
            self.perturbation_order_total = self.triqs_solver.perturbation_order_total

        if self.solver_params['measure_chi'] is not None:
            self.O_time = self.triqs_solver.O_tau

        return

