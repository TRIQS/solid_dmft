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
hubbardI solver class for solid_dmft
'''
from triqs.gfs import MeshReFreq, Gf, make_hermitian

# import of the abstract class
from solid_dmft.dmft_tools.solvers.abstractdmftsolver import AbstractDMFTSolver

# import triqs solver
from triqs_hubbardI import Solver as hubbardI_solver
from triqs_hubbardI.version import triqs_hubbardI_hash, version

# import tail fitting from cthyb
from triqs_cthyb.tail_fit import sigma_high_frequency_moments, green_high_frequency_moments
from triqs_cthyb.util import orbital_occupations


class HubbardIInterface(AbstractDMFTSolver):
    def __init__(
        self, general_params, solver_params, sum_k, icrsh, h_int, iteration_offset, deg_orbs_ftps, gw_params=None, advanced_params=None
    ):
        # Call the base class constructor
        super().__init__(general_params, solver_params, sum_k, icrsh, h_int, iteration_offset, deg_orbs_ftps, gw_params, advanced_params)

        # Create the hartree solver specifics
        # Separately stores all params that go into solve() call of solver
        # All params need to be renamed
        self.triqs_solver_params = {}
        self.triqs_solver_params['calc_gtau'] = self.solver_params['measure_G_tau']
        self.triqs_solver_params['calc_gw'] = True
        self.triqs_solver_params['calc_gl'] = self.solver_params['measure_G_l']
        self.triqs_solver_params['calc_dm'] = self.solver_params['measure_density_matrix']

        gf_struct = self.sum_k.gf_struct_solver_list[self.icrsh]
        # Construct the triqs_solver instances
        if self.solver_params['measure_G_l']:
            self.triqs_solver = hubbardI_solver(
                beta=self.general_params['beta'],
                gf_struct=gf_struct,
                n_iw=self.general_params['n_iw'],
                n_tau=self.general_params['n_tau'],
                n_l=self.solver_params['n_l'],
                n_w=self.general_params['n_w'],
                w_min=self.general_params['w_range'][0],
                w_max=self.general_params['w_range'][1],
                idelta=self.general_params['eta'],
            )
        else:
            self.triqs_solver = hubbardI_solver(
                beta=self.general_params['beta'],
                gf_struct=gf_struct,
                n_iw=self.general_params['n_iw'],
                n_tau=self.general_params['n_tau'],
                n_w=self.general_params['n_w'],
                idelta=self.general_params['eta'],
                w_min=self.general_params['w_range'][0],
                w_max=self.general_params['w_range'][1],
            )

        # sets up necessary GF objects on ImFreq
        self._init_ImFreq_objects()
        self._init_ReFreq_hubbardI()  # definition at the end of the class

        # sets up solver
        self.git_hash = triqs_hubbardI_hash
        self.version = version
        return

    def solve(self, **kwargs):
        # fill G0_freq from sum_k to solver
        self.triqs_solver.G0_iw << self.G0_freq

        # Solve the impurity problem for icrsh shell
        # *************************************
        # this is done on every node due to very slow bcast of the AtomDiag object as of now
        self.triqs_solver.solve(h_int=self.h_int, **self.triqs_solver_params)
        # if density matrix is measured, get this too. Needs to be done here,
        # because solver property 'dm' is not initialized/broadcastable
        if self.solver_params['measure_density_matrix']:
            self.density_matrix = self.triqs_solver.dm
            self.h_loc_diagonalization = self.triqs_solver.ad
            # get moments
            self.Sigma_moments = sigma_high_frequency_moments(
                self.density_matrix, self.h_loc_diagonalization, self.sum_k.gf_struct_solver_list[self.icrsh], self.h_int
            )
            self.Sigma_Hartree = {bl: sigma_bl[0] for bl, sigma_bl in self.Sigma_moments.items()}
            self.G_moments = green_high_frequency_moments(
                self.density_matrix, self.h_loc_diagonalization, self.sum_k.gf_struct_solver_list[self.icrsh], self.h_int
            )
            self.orbital_occupations = orbital_occupations(
                self.density_matrix, self.sum_k.gf_struct_solver_list[self.icrsh], self.h_loc_diagonalization
            )

        # *************************************

        # call postprocessing
        self.postprocess()

        return

    def postprocess(self):
        r"""
        Organize G_freq, G_time, Sigma_freq and G_l from hubbardI solver
        """

        # get everything from solver
        self.Sigma_freq << self.triqs_solver.Sigma_iw
        self.G0_freq << self.triqs_solver.G0_iw
        self.G0_Refreq << self.triqs_solver.G0_w
        self.G_freq << make_hermitian(self.triqs_solver.G_iw)
        self.G_freq_unsym << self.triqs_solver.G_iw
        self.sum_k.symm_deg_gf(self.G_freq, ish=self.icrsh)
        self.G_freq << self.G_freq
        self.G_Refreq << self.triqs_solver.G_w
        self.Sigma_Refreq << self.triqs_solver.Sigma_w

        # if measured in Legendre basis, get G_l from solver too
        if self.solver_params['measure_G_l']:
            self.G_l << self.triqs_solver.G_l

        if self.solver_params['measure_G_tau']:
            self.G_time << self.triqs_solver.G_tau

        return

    def _init_ReFreq_hubbardI(self):
        r"""
        Initialize all ReFreq objects
        """

        # create all ReFreq instances
        self.n_w = self.general_params['n_w']
        self.G_Refreq = self.sum_k.block_structure.create_gf(
            ish=self.icrsh, gf_function=Gf, space='solver', mesh=MeshReFreq(n_w=self.n_w, window=self.general_params['w_range'])
        )
        # copy
        self.Sigma_Refreq = self.G_Refreq.copy()
        self.G0_Refreq = self.G_Refreq.copy()
