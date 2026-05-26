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
ctseg solver class for solid_dmft
'''
import numpy as np
from itertools import product

from triqs.gfs import MeshDLRImFreq, Gf, BlockGf, make_gf_imfreq, make_hermitian, make_gf_dlr, fit_gf_dlr, make_gf_dlr_imtime, make_gf_imtime
from triqs.gfs.tools import inverse, make_zero_tail
from triqs.gfs.descriptors import Fourier
from triqs.operators.util.U_matrix import reduce_4index_to_2index
import triqs.utility.mpi as mpi
from h5 import HDFArchive

from solid_dmft.io_tools.dict_to_h5 import prep_params_for_h5

from solid_dmft.dmft_tools import legendre_filter
from solid_dmft.dmft_tools.matheval import MathExpr


#  import of the abstract class
from solid_dmft.dmft_tools.solvers.abstractdmftsolver import AbstractDMFTSolver
from solid_dmft.dmft_tools import common


# import triqs solver
from triqs_ctseg import Solver as ctseg_solver
from triqs_ctseg.version import triqs_ctseg_hash, version


class CTSEGInterface(AbstractDMFTSolver):
    def __init__(
        self, general_params, solver_params, sum_k, icrsh, h_int, iteration_offset, deg_orbs_ftps, gw_params=None, advanced_params=None
    ):
        r"""
        Initialize cthyb solver instance
        """

        # Call the base class constructor
        super().__init__(general_params, solver_params, sum_k, icrsh, h_int, iteration_offset, deg_orbs_ftps, gw_params, advanced_params)

        # sets up necessary GF objects on ImFreq
        self._init_ImFreq_objects()

        if self.solver_params.get('random_seed') is None:
            self.random_seed_generator = None
        else:
            self.random_seed_generator = MathExpr(self.solver_params['random_seed'])

        # Separately stores all params that go into solve() call of solver
        self.triqs_solver_params = {}
        keys_to_pass = (
            'length_cycle',
            'max_time',
            'measure_state_hist',
            'measure_nn_tau',
            'measure_G_tau',
            'measure_pert_order',
        )
        for key in keys_to_pass:
            self.triqs_solver_params[key] = self.solver_params[key]

        # Calculates number of sweeps per rank
        self.triqs_solver_params['n_cycles'] = int(self.solver_params['n_cycles_tot'] / mpi.size)
        # cast warmup cycles to int in case given in scientific notation
        self.triqs_solver_params['n_warmup_cycles'] = int(self.solver_params['n_warmup_cycles'])

        # Makes sure measure_gw is true if improved estimators are used
        if self.solver_params['improved_estimator']:
            self.triqs_solver_params['measure_G_tau'] = True
            self.triqs_solver_params['measure_F_tau'] = True
        else:
            self.triqs_solver_params['measure_F_tau'] = False

        gf_struct = self.sum_k.gf_struct_solver_list[self.icrsh]

        # Construct the triqs_solver instances
        self.triqs_solver = ctseg_solver(
            beta=self.general_params['beta'],
            gf_struct=gf_struct,
            n_tau=self.general_params['n_tau'],
            n_tau_bosonic=int(self.solver_params['n_tau_bosonic']),
        )

        self.git_hash = triqs_ctseg_hash
        self.version = version

    def solve(self, **kwargs):
        # what does this do exactly?
        if self.random_seed_generator is not None:
            self.triqs_solver_params['random_seed'] = int(self.random_seed_generator(it=kwargs['it'], rank=mpi.rank))
        else:
            assert 'random_seed' not in self.triqs_solver_params

        # fill G0_freq from sum_k to solver
        self.triqs_solver.Delta_tau << self.Delta_time

        if self.general_params['h_int_type'][self.icrsh] == 'dyn_density_density':
            mpi.report('\nAdding dynamic interaction from AIMBES.')
            # convert 4 idx tensor to two index tensor
            Uloc_dlr = self.gw_params['Uloc_dlr'][self.icrsh]['up']
            Uloc_dlr_2idx_prime = Gf(mesh=Uloc_dlr.mesh, target_shape=[Uloc_dlr.target_shape[0], Uloc_dlr.target_shape[1]])

            for coeff in Uloc_dlr.mesh:
                Uloc_dlr_idx = Uloc_dlr[coeff]
                _, Uprime = reduce_4index_to_2index(Uloc_dlr_idx)
                Uloc_dlr_2idx_prime[coeff] = Uprime

            # extract w=0 limit for analytic Sigma_Hartree for the impurity
            Uloc_w0_2idx_prime = make_gf_imfreq(Uloc_dlr_2idx_prime, n_iw=1)
            self.Uw0_prime = Uloc_w0_2idx_prime.data[0].real
            mpi.report('Screened interaction U(iw) at iw=0 w/o the shift of V_bare (density-density only): ')
            mpi.report(self.Uw0_prime)

            # create full frequency objects
            Uloc_tau_2idx_prime = make_gf_imtime(Uloc_dlr_2idx_prime, n_tau=self.solver_params['n_tau_bosonic'])

            # fill D0_tau from Uloc_tau_2idx and Uloc_tau_2idx_prime
            ish = self.sum_k.inequiv_to_corr[self.icrsh]
            norb = Uloc_dlr.target_shape[0]
            gf_struct = self.sum_k.gf_struct_solver_list[ish]
            o1 = 0
            for name1, n1 in gf_struct:
                o2 = 0
                for name2, n2 in gf_struct:
                    # same and opposite spin are the same
                    self.triqs_solver.D0_tau[name1, name2] << Uloc_tau_2idx_prime[o1:o1 + n1, o2:o2 + n2].real
                    o2 = (o2 + n2) % norb
                o1 = (o1 + n1) % norb

            # TODO: add Jerp_Iw to the solver

            # self.triqs_solver. Jperp_iw << make_gf_imfreq(Uloc_dlr_2idx, n_iw=self.general_params['n_w_b_nn']) + V
        mpi.report('\nLocal interaction Hamiltonian is:', self.h_int)

        # update solver in h5 archive one last time for debugging if solve command crashes
        if self.general_params['store_solver'] and mpi.is_master_node():
            with HDFArchive(self.general_params['jobname'] + '/' + self.general_params['seedname'] + '.h5', 'a') as archive:
                if 'it_-1' not in archive['DMFT_input/solver']:
                    archive['DMFT_input/solver'].create_group('it_-1')
                archive['DMFT_input/solver/it_-1'][f'S_{self.icrsh}'] = self.triqs_solver
                archive['DMFT_input/solver/it_-1'][f'Delta_time_{self.icrsh}'] = self.triqs_solver.Delta_tau
                archive['DMFT_input/solver/it_-1'][f'solve_params_{self.icrsh}'] = prep_params_for_h5(self.solver_params)
                archive['DMFT_input/solver/it_-1'][f'triqs_solver_params_{self.icrsh}'] = prep_params_for_h5(self.triqs_solver_params)
                archive['DMFT_input/solver/it_-1']['mpi_size'] = mpi.size
                if self.general_params['h_int_type'][self.icrsh] == 'dyn_density_density':
                    archive['DMFT_input/solver/it_-1'][f'Uloc_dlr_2idx_prime_{self.icrsh}'] = Uloc_dlr_2idx_prime
        mpi.barrier()

        # turn of problematic move in ctseg until fixed!
        self.triqs_solver_params['move_move_segment'] = False
        # Solve the impurity problem for icrsh shell
        # *************************************
        self.triqs_solver.solve(h_int=self.h_int, h_loc0=self.Hloc_0, **self.triqs_solver_params)
        # *************************************

        # call postprocessing
        self.postprocess()

    def postprocess(self):
        r"""
        Organize G_freq, G_time, Sigma_freq and G_l from ctseg solver
        """
        from triqs.operators.util.extractors import extract_U_dict2, dict_to_matrix
        from solid_dmft.postprocessing.eval_U_cRPA_RESPACK import construct_Uijkl

        def set_Gs_from_G_l():
            if self.solver_params['improved_estimator'] and mpi.is_master_node():
                print(
                    '\n !!!!WARNING!!!! \n you enabled both improved estimators and legendre based filtering / sampling. Sigma will be overwritten by legendre result.  \n !!!!WARNING!!!!\n'
                )

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

        self.G_time << self.triqs_solver.results.G_tau
        self.sum_k.symm_deg_gf(self.G_time, ish=self.icrsh)

        # get occupation matrix
        self.orbital_occupations = {bl: np.zeros((bl_size, bl_size)) for bl, bl_size in self.sum_k.gf_struct_solver_list[self.icrsh]}
        for block, norb in self.sum_k.gf_struct_solver[self.icrsh].items():
            self.orbital_occupations[block] = np.zeros((norb, norb))
            for iorb in range(norb):
                self.orbital_occupations[block][iorb, iorb] = self.triqs_solver.results.densities[block][iorb]

        self.orbital_occupations_sumk = self.sum_k.block_structure.convert_matrix(
            self.orbital_occupations, ish_from=self.icrsh, space_from='solver', space_to='sumk'
        )
        self.Sigma_Hartree = {}
        self.Sigma_Hartree_sumk = {}
        self.Sigma_moments = {}
        if mpi.is_master_node():
            mpi.report('Evaluating static impurity self-energy analytically using interacting density from ctseg...\n'
                       '(results will be used in the subsequent tail fitting or the crm dyson solver)')
            # get density density U tensor from solver
            U_dict = extract_U_dict2(self.h_int)
            # print("sum_k is" + self.sum_k.__repr__())
            norb = common.get_n_orbitals(self.sum_k)
            norb = norb[self.icrsh]['up']
            U_dd = dict_to_matrix(U_dict, gf_struct=self.sum_k.gf_struct_solver_list[self.icrsh])
            # extract Uijij (inter- and intra-orbital Coulomb) and Uijji (Hund's coupling) terms
            # a) For static impurity problem, Us are the static screened interactions
            # b) For dynamic impurity problem, Us are the bare interactions
            Uijij = U_dd[0:norb, norb : 2 * norb]
            Uijji = Uijij - U_dd[0:norb, 0:norb]
            # and construct full Uijkl tensor for static interaction
            Uijkl = construct_Uijkl(Uijij, Uijji)

            if self.general_params['h_int_type'][self.icrsh] == 'dyn_density_density':
                # For dynamic impurity problems, separate Us are needed for the Hartree and exchange self-energy
                # a) Hartree term is evaluated using the screened interactions at w=0
                # b) Exchange term is evaluated using the bare interactions
                Uijij += self.Uw0_prime
                Uijji[:, :] = 0.0
                Uw0_ijkl = construct_Uijkl(Uijij, Uijji)
            else:
                Uw0_ijkl = Uijkl

            # now calculated Hartree shift via
            # \Sigma^0_{\alpha \beta} = \sum_{i j} n_{i j} \left( 2 Uw0_{\alpha i \beta j} - U_{\alpha i j \beta} \right)
            for block, norb in self.sum_k.gf_struct_sumk[self.icrsh]:
                self.Sigma_Hartree_sumk[block] = np.zeros((norb, norb), dtype=float)
                for iorb, jorb in product(range(norb), repeat=2):
                    for inner in range(norb):
                        # exchange diagram K
                        self.Sigma_Hartree_sumk[block][iorb, jorb] -= (
                                self.orbital_occupations_sumk[block][inner, inner].real * (Uijkl[iorb, inner, inner, jorb].real))
                        # Hartree (Coulomb) diagram J
                        self.Sigma_Hartree_sumk[block][iorb, jorb] += (
                                self.orbital_occupations_sumk[block][inner, inner].real * (2 * Uw0_ijkl[iorb, inner, jorb, inner].real))

            # convert to solver block structure
            self.Sigma_Hartree = self.sum_k.block_structure.convert_matrix(
                self.Sigma_Hartree_sumk, ish_from=self.icrsh, space_from='sumk', space_to='solver'
            )

            # use degenerate shell information to symmetrize
            self.sum_k.symm_deg_gf(self.Sigma_Hartree, ish=self.icrsh)

            # create moments array from this
            for block, hf_val in self.Sigma_Hartree.items():
                self.Sigma_moments[block] = np.array([hf_val])

        self.Sigma_Hartree = mpi.bcast(self.Sigma_Hartree)
        self.Sigma_moments = mpi.bcast(self.Sigma_moments)
        self.Sigma_Hartree_sumk = mpi.bcast(self.Sigma_Hartree_sumk)

        if mpi.is_master_node():
            # create empty moment container (list of np.arrays)
            Gf_known_moments = make_zero_tail(self.G_freq, n_moments=2)
            for i, bl in enumerate(self.G_freq.indices):
                # 0 moment is 0, dont touch it, but first moment is 1 for the Gf
                Gf_known_moments[i][1] = np.eye(self.G_freq[bl].target_shape[0])
                self.G_freq[bl] << Fourier(self.G_time[bl], Gf_known_moments[i])
        self.G_freq << mpi.bcast(self.G_freq)

        self.G_freq << make_hermitian(self.G_freq)
        self.G_freq_unsym << self.G_freq
        self.sum_k.symm_deg_gf(self.G_freq, ish=self.icrsh)

        # if measured in Legendre basis, get G_l from solver too
        if self.solver_params['legendre_fit']:
            mpi.report('Self-energy post-processing algorithm: Legendre fitting')
            self.G_time_orig << self.triqs_solver.results.G_tau
            self.G_l << legendre_filter.apply(self.G_time, self.solver_params['n_l'])
            # get G_time, G_freq, Sigma_freq from G_l
            set_Gs_from_G_l()
        elif self.solver_params['perform_tail_fit']:
            if not self.solver_params['improved_estimator']:
                mpi.report('Self-energy post-processing algorithm: tail fitting with analytic static impurity self-energy')
                self.Sigma_freq = inverse(self.G0_freq) - inverse(self.G_freq)
            else:
                mpi.report('Self-energy post-processing algorithm: '
                           'improved estimator + tail fitting with analytic static imppurity self-energy')
                self.F_freq = self.G_freq.copy()
                self.F_freq << 0.0
                self.F_time = self.G_time.copy()
                self.F_time << self.triqs_solver.results.F_tau
                F_known_moments = make_zero_tail(self.F_freq, n_moments=1)
                for i, bl in enumerate(self.F_freq.indices):
                    self.F_freq[bl] << Fourier(self.triqs_solver.results.F_tau[bl], F_known_moments[i])

                for block, fw in self.F_freq:
                    for iw in fw.mesh:
                        self.Sigma_freq[block][iw] = self.F_freq[block][iw] / self.G_freq[block][iw]

            # without any degenerate shells we run the minimization for all blocks
            self.Sigma_freq, tail = self._fit_tail_window(
                self.Sigma_freq,
                fit_min_n=self.solver_params['fit_min_n'],
                fit_max_n=self.solver_params['fit_max_n'],
                fit_min_w=self.solver_params['fit_min_w'],
                fit_max_w=self.solver_params['fit_max_w'],
                fit_max_moment=self.solver_params['fit_max_moment'],
                fit_known_moments=self.Sigma_moments,
            )

            # recompute G_freq from Sigma with fitted tail
            self.G_freq = inverse(inverse(self.G0_freq) - self.Sigma_freq)

        # if improved estimators are turned on calc Sigma from F_tau, otherwise:
        elif self.solver_params['improved_estimator'] and not self.solver_params['perform_tail_fit']:
            mpi.report('Self-energy post-processing algorithm: improved estimator')
            self.F_freq = self.G_freq.copy()
            self.F_freq << 0.0
            self.F_time = self.G_time.copy()
            self.F_time << self.triqs_solver.results.F_tau
            F_known_moments = make_zero_tail(self.F_freq, n_moments=1)
            for i, bl in enumerate(self.F_freq.indices):
                self.F_freq[bl] << Fourier(self.triqs_solver.results.F_tau[bl], F_known_moments[i])

            for block, fw in self.F_freq:
                for iw in fw.mesh:
                    self.Sigma_freq[block][iw] = self.F_freq[block][iw] / self.G_freq[block][iw]

        # should this be moved to abstract class?
        elif self.solver_params['crm_dyson_solver']:
            mpi.report('Self-energy post-processing algorithm: crm dyson solver')
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
                G_dlr = fit_gf_dlr(self.triqs_solver.results.G_tau, w_max=dlr_wmax, eps=dlr_eps)
                self.G_time_dlr = make_gf_dlr_imtime(G_dlr)
                self.G_freq = make_gf_imfreq(G_dlr, n_iw=self.general_params['n_iw'])

                # assume little error on G0_iw and use to get G0_dlr
                mesh_dlr_iw = MeshDLRImFreq(G_dlr.mesh)
                G0_dlr_iw = self.sum_k.block_structure.create_gf(ish=self.icrsh, gf_function=Gf, mesh=mesh_dlr_iw, space='solver')
                for block, gf in G0_dlr_iw:
                    for iwn in mesh_dlr_iw:
                        gf[iwn] = self.G0_freq[block](iwn)
                self.sum_k.symm_deg_gf(G0_dlr_iw, ish=self.icrsh)
                G0_dlr = make_gf_dlr(G0_dlr_iw)

                # get Hartree shift for optimizer
                G0_iw = make_gf_imfreq(G0_dlr, n_iw=self.general_params['n_iw'])
                G_iw = make_gf_imfreq(G_dlr, n_iw=self.general_params['n_iw'])
                Sigma_iw = inverse(G0_iw) - inverse(G_iw)

                # minimize dyson for the first entry of each deg shell
                self.Sigma_dlr = self.sum_k.block_structure.create_gf(ish=self.icrsh, gf_function=Gf, mesh=mesh_dlr_iw, space='solver')
                # without any degenerate shells we run the minimization for all blocks
                if self.sum_k.deg_shells[self.icrsh] == []:
                    for block, gf in self.Sigma_dlr:
                        np.random.seed(85281)
                        print('Minimizing Dyson via CRM for Σ[block]:', block)
                        gf, _, _ = minimize_dyson(G0_dlr=G0_dlr_iw[block], G_dlr=G_dlr[block], Sigma_moments=self.Sigma_moments[block])
                else:
                    for deg_shell in self.sum_k.deg_shells[self.icrsh]:
                        for i, block in enumerate(deg_shell):
                            if i == 0:
                                np.random.seed(85281)
                                print('Minimizing Dyson via CRM for Σ[block]:', block)
                                self.Sigma_dlr[block], _, _ = minimize_dyson(
                                    G0_dlr=G0_dlr_iw[block], G_dlr=G_dlr[block], Sigma_moments=self.Sigma_moments[block]
                                )
                                sol_block = block
                            else:
                                print(f'Copying result from first block of deg list to {block}')
                                self.Sigma_dlr[block] << self.Sigma_dlr[sol_block]

                            self.Sigma_freq[block] = make_gf_imfreq(self.Sigma_dlr[block], n_iw=self.general_params['n_iw'])
                            self.Sigma_freq[block] += self.Sigma_moments[block][0]

                self.G_freq = inverse(inverse(self.G0_freq) - self.Sigma_freq)
                print('\n')

            mpi.barrier()
            self.Sigma_freq = mpi.bcast(self.Sigma_freq)
            self.Sigma_dlr = mpi.bcast(self.Sigma_dlr)
            self.G_time_dlr = mpi.bcast(self.G_time_dlr)
            self.G_freq = mpi.bcast(self.G_freq)
        else:
            mpi.report('\n!!!! WARNING !!!! tail of solver output not handled! Turn on either measure_F_tau, legendre_fit\n')
            self.Sigma_freq << inverse(self.G0_freq) - inverse(self.G_freq)

        if self.solver_params['measure_state_hist']:
            self.state_histogram = self.triqs_solver.results.state_hist

        if self.solver_params['measure_pert_order']:
            self.perturbation_order_histo = self.triqs_solver.results.pert_order_Delta
            self.avg_pert_order = self.triqs_solver.results.average_order_Delta

        return
