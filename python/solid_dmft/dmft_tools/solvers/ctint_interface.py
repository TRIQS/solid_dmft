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
ctint solver class for solid_dmft
'''
from itertools import product

from triqs.gfs import make_hermitian
from triqs.gfs.tools import inverse
from triqs.gfs.descriptors import Fourier
import triqs.utility.mpi as mpi
from h5 import HDFArchive

from solid_dmft.dmft_tools.matheval import MathExpr


# import of the abstract class
from solid_dmft.dmft_tools.solvers.abstractdmftsolver import AbstractDMFTSolver

# import triqs solver
from triqs_ctint import Solver as ctint_solver
from triqs_ctint.version import triqs_ctint_hash, version


class CTINTInterface(AbstractDMFTSolver):
    def __init__(
        self, general_params, solver_params, sum_k, icrsh, h_int, iteration_offset, deg_orbs_ftps, gw_params=None, advanced_params=None
    ):
        # Call the base class constructor
        super().__init__(general_params, solver_params, sum_k, icrsh, h_int, iteration_offset, deg_orbs_ftps, gw_params, advanced_params)

        # sets up necessary GF objects on ImFreq
        self._init_ImFreq_objects()

        ###################################################
        # Create the ctint solver specifics
        ###################################################
        self._init_ImFreq_objects()
        # set up solver

        if self.solver_params.get('random_seed') is None:
            self.random_seed_generator = None
        else:
            self.random_seed_generator = MathExpr(self.solver_params['random_seed'])

        # Separately stores all params that go into solve() call of solver
        self.triqs_solver_params = {}

        # keys with same name
        keys_to_pass = ('length_cycle', 'max_time', 'n_warmup_cycles')
        for key in keys_to_pass:
            self.triqs_solver_params[key] = self.solver_params[key]

        # keys with different name
        self.triqs_solver_params['measure_histogram'] = self.solver_params.get('measure_pert_order')
        self.triqs_solver_params['use_double_insertion'] = self.solver_params.get('move_double')

        # Calculates number of sweeps per rank
        self.triqs_solver_params['n_cycles'] = int(self.solver_params['n_cycles_tot'] / mpi.size)

        gf_struct = self.sum_k.gf_struct_solver_list[self.icrsh]

        if self.general_params['h_int_type'][self.icrsh] == 'dyn_density_density':
            self.U_iw = None
            if mpi.is_master_node():
                with HDFArchive(self.general_params['jobname'] + '/' + self.general_params['seedname'] + '.h5', 'r') as archive:
                    self.U_iw = archive['dynamic_U']['U_iw']
            self.U_iw = mpi.bcast(self.U_iw)
            n_iw_dyn = self.U_iw[self.icrsh].mesh.last_index() + 1
            # Construct the triqs_solver instances
            self.triqs_solver = ctint_solver(
                beta=self.general_params['beta'],
                gf_struct=gf_struct,
                n_iw=self.general_params['n_iw'],
                n_tau=self.general_params['n_tau'],
                use_D=True,
                use_Jperp=False,
                n_iw_dynamical_interactions=n_iw_dyn,
                n_tau_dynamical_interactions=(int(n_iw_dyn * 2.5)),
            )
        else:
            # Construct the triqs_solver instances
            self.triqs_solver = ctint_solver(
                beta=self.general_params['beta'],
                gf_struct=gf_struct,
                n_iw=self.general_params['n_iw'],
                n_tau=self.general_params['n_tau'],
                use_D=False,
                use_Jperp=False,
            )

        # set up metadata
        self.git_hash = triqs_ctint_hash
        self.version = version

        return

    def solve(self, **kwargs):
        # what does this do exactly?
        if self.random_seed_generator is not None:
            self.triqs_solver_params['random_seed'] = int(self.random_seed_generator(it=kwargs['it'], rank=mpi.rank))
        else:
            assert 'random_seed' not in self.triqs_solver_params

        # fill G0_freq from sum_k to solver
        self.triqs_solver.G0_iw << self.G0_freq

        if self.general_params['h_int_type'][self.icrsh] == 'dynamic':
            for b1, b2 in product(self.sum_k.gf_struct_solver_dict[self.icrsh].keys(), repeat=2):
                self.triqs_solver.D0_iw[b1, b2] << self.U_iw[self.icrsh]

        # Solve the impurity problem for icrsh shell
        # *************************************
        self.triqs_solver.solve(h_int=self.h_int, **self.triqs_solver_params)
        # *************************************

        # call postprocessing
        self.postprocess()

        return

    def postprocess(self):
        r"""
        Organize G_freq, G_time, Sigma_freq and G_l from cthyb solver
        """
        # TODO

        # def set_Gs_from_G_l():

        #     # create new G_freq and G_time
        #     for i, g in self.G_l:
        #         g.enforce_discontinuity(np.identity(g.target_shape[0]))
        #         # set G_freq from Legendre and Fouriertransform to get G_time
        #         self.G_freq[i].set_from_legendre(g)
        #         self.G_time[i] << Fourier(self.G_freq[i])
        #     # Symmetrize
        #     self.G_freq << make_hermitian(self.G_freq)
        #     # Dyson equation to get Sigma_freq
        #     self.Sigma_freq << inverse(self.G0_freq) - inverse(self.G_freq)

        #     return

        self.G_freq << make_hermitian(self.triqs_solver.G_iw)
        self.G_freq_unsym << self.G_freq
        self.sum_k.symm_deg_gf(self.G_freq, ish=self.icrsh)
        self.Sigma_freq << inverse(self.G0_freq) - inverse(self.G_freq)
        self.G_time << Fourier(self.G_freq)

        # TODO: probably not needed/sensible
        # if self.solver_params['legendre_fit']:
        #     self.G_freq_orig << self.triqs_solver.G_iw
        #     # run the filter
        #     self.G_l << legendre_filter.apply(self.G_time, self.solver_params['n_l'])
        #     # get G_time, G_freq, Sigma_freq from G_l
        #     set_Gs_from_G_l()

        if self.solver_params['measure_pert_order']:
            self.perturbation_order = self.triqs_solver.histogram

        return
