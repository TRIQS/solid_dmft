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
hartree solver class for solid_dmft
'''
from triqs.gfs import MeshReFreq, Gf, inverse, make_gf_dlr_imfreq, make_gf_imfreq, find_w_max
from triqs.gfs.descriptors import Fourier
import triqs.utility.mpi as mpi

# import of the abstract class
from solid_dmft.dmft_tools.solvers.abstractdmftsolver import AbstractDMFTSolver

# import triqs solver
from triqs_hartree_fock import ImpuritySolver as hartree_solver
from triqs_hartree_fock.version import triqs_hartree_fock_hash, version


class HartreeInterface(AbstractDMFTSolver):
    def __init__(
        self, general_params, solver_params, sum_k, icrsh, h_int, iteration_offset, deg_orbs_ftps, gw_params=None, advanced_params=None
    ):
        # Call the base class constructor
        super().__init__(general_params, solver_params, sum_k, icrsh, h_int, iteration_offset, deg_orbs_ftps, gw_params, advanced_params)

        # Create the hartree solver specifics

        self.triqs_solver_params = {}
        keys_to_pass = ('method', 'one_shot', 'tol', 'with_fock')
        for key in keys_to_pass:
            self.triqs_solver_params[key] = self.solver_params[key]

        # sets up necessary GF objects on ImFreq
        self._init_ImFreq_objects()
        self._init_ReFreq_hartree()  # definition at the end of the class

        # Construct the triqs_solver instances
        # Always initialize the solver with dc_U and dc_J equal to U and J and let the _interface_hartree_dc function
        # take care of changing the parameters
        gf_struct = self.sum_k.gf_struct_solver_list[self.icrsh]

        # The Hartree solver now stores its Green's functions on a DLR Matsubara mesh
        # (triqs_hartree_fock commit c948fe7). solid_dmft still works with full Matsubara
        # meshes for the lattice Green's function, so we bridge between the two with the
        # triqs functions make_gf_dlr_imfreq / make_gf_imfreq / find_w_max (commit a43d8cc):
        # the regular-mesh Weiss field G0 is sampled onto a DLR mesh in solve() and the
        # solver's DLR outputs are converted back in postprocess(). These two attributes
        # fully define that DLR mesh and are kept on the instance so they are available if
        # the lattice Green's function is itself moved onto a DLR mesh in the future.
        self.dlr_eps = self.solver_params['eps_dlr']
        self.dlr_w_max = None  # set in solve() via find_w_max from the Weiss field G0

        # w_max here is only a placeholder for the initial (empty) construction; the actual
        # DLR mesh is built in solve() once the Weiss field is known.
        self.triqs_solver = hartree_solver(
            beta=self.general_params['beta'],
            gf_struct=gf_struct,
            w_max=1.0,
            eps=self.dlr_eps,
            force_real=self.solver_params['force_real'],
            symmetries=[self._make_spin_equal],
            dc_U=self.general_params['U'][self.icrsh],
            dc_J=self.general_params['J'][self.icrsh],
        )

        # Give dc information to the solver in order to customize DC calculation
        def _interface_hartree_dc(hartree_instance, general_params, advanced_params, icrsh):
            """Modifies in-place class attributes to infercace with options in solid_dmft
                for the moment supports only DC-relevant parameters

            Parameters
            ----------
                general_params : dict
                    solid_dmft general parameter dictionary
                advanced_params : dict
                    solid_dmft advanced parameter dictionary
                icrsh : int
                    correlated shell number
            """
            setattr(hartree_instance, 'dc', general_params['dc'])
            if general_params['dc_type'][icrsh] is not None:
                setattr(hartree_instance, 'dc_type', general_params['dc_type'][icrsh])

            for key in ['dc_factor', 'dc_fixed_value']:
                if key in advanced_params and advanced_params[key] is not None:
                    setattr(hartree_instance, key, advanced_params[key])

            # list valued keys
            for key in ['dc_U', 'dc_J', 'dc_fixed_occ']:
                if key in advanced_params and advanced_params[key][icrsh] is not None:
                    setattr(hartree_instance, key, advanced_params[key][icrsh])

            # Handle special cases
            if 'dc_dmft' in general_params:
                if general_params['dc_dmft'] == False:
                    mpi.report(
                        'HARTREE SOLVER: Warning dft occupation in the DC calculations are meaningless for the hartree solver, reverting to dmft occupations'
                    )

            if hartree_instance.dc_type == 0 and not self.general_params['magnetic']:
                mpi.report(f"HARTREE SOLVER: Detected dc_type = {hartree_instance.dc_type}, changing to 'cFLL'")
                hartree_instance.dc_type = 'cFLL'
            elif hartree_instance.dc_type == 0 and self.general_params['magnetic']:
                mpi.report(f"HARTREE SOLVER: Detected dc_type = {hartree_instance.dc_type}, changing to 'sFLL'")
                hartree_instance.dc_type = 'sFLL'
            elif hartree_instance.dc_type == 1:
                mpi.report(f"HARTREE SOLVER: Detected dc_type = {hartree_instance.dc_type}, changing to 'cHeld'")
                hartree_instance.dc_type = 'cHeld'
            elif hartree_instance.dc_type == 2 and not self.general_params['magnetic']:
                mpi.report(f"HARTREE SOLVER: Detected dc_type = {hartree_instance.dc_type}, changing to 'cAMF'")
                hartree_instance.dc_type = 'cAMF'
            elif hartree_instance.dc_type == 2 and self.general_params['magnetic']:
                mpi.report(f"HARTREE SOLVER: Detected dc_type = {hartree_instance.dc_type}, changing to 'sAMF'")
                hartree_instance.dc_type = 'sAMF'

        # Give dc information to the solver in order to customize DC calculation
        _interface_hartree_dc(self.triqs_solver, self.general_params, self.advanced_params, self.icrsh)

        # set up metadata
        self.git_hash = triqs_hartree_fock_hash
        self.version = version

        return

    def _init_ReFreq_hartree(self):
        r"""
        Initialize all ReFreq objects
        """

        # create all ReFreq instances
        self.n_w = self.general_params['n_w']
        self.Sigma_Refreq = self.sum_k.block_structure.create_gf(
            ish=self.icrsh, gf_function=Gf, space='solver', mesh=MeshReFreq(n_w=self.n_w, window=self.general_params['w_range'])
        )

    def _largest_fitting_dlr_w_max(self, g):
        # Largest DLR cutoff whose Matsubara nodes still fit inside the n_iw mesh of g.
        # find_w_max(G0) is guaranteed to fit (G0 lives on the same mesh), so grow from
        # there by 1.5x until make_gf_dlr_imfreq can no longer sample g.
        w_max = find_w_max(self.G0_freq, self.dlr_eps)
        trial = w_max * 1.5
        while trial <= 200.0:
            try:
                make_gf_dlr_imfreq(g, trial, self.dlr_eps)
            except RuntimeError:
                break
            w_max = trial
            trial *= 1.5
        return w_max

    def _set_dlr_G0(self, w_max):
        # Sample the regular-mesh Weiss field onto a DLR Matsubara mesh of cutoff w_max
        # and hand it to the solver (G_iw is allocated on the same mesh for the solver).
        self.triqs_solver.w_max = w_max
        self.triqs_solver.G0_iw = make_gf_dlr_imfreq(self.G0_freq, w_max, self.dlr_eps)
        self.triqs_solver.G_iw = self.triqs_solver.G0_iw.copy()

    def solve(self, **kwargs):
        # The DLR mesh must be sized to the impurity Green's function, not just the Weiss
        # field G0: in Hartree-Fock the self energy is a constant matrix that shifts G's
        # spectral weight away from G0, so a mesh sized to G0 alone under-resolves G (and
        # hence its density). Sigma_HF is unknown before solving, and on the first iteration
        # the seeded guess is zero (the Hartree DC is computed inside the solver), so we run
        # a cheap estimate solve on a G0-sized mesh to obtain the actual (constant) Sigma_HF,
        # build a guesstimate G from it, and size the real DLR mesh from that via find_w_max.
        # The estimate solve overwrites Sigma_HF, so we restart the final solve from the
        # original guess to leave the one_shot result unchanged.
        sigma_hf_init = {bl: s.copy() for bl, s in self.triqs_solver.Sigma_HF.items()}

        self._set_dlr_G0(find_w_max(self.G0_freq, self.dlr_eps))
        self.triqs_solver.solve(h_int=self.h_int, **self.triqs_solver_params)
        sigma_est = {bl: s.copy() for bl, s in self.triqs_solver.Sigma_HF.items()}

        G_guess = self.G0_freq.copy()
        for bl, g in G_guess:
            g << inverse(inverse(self.G0_freq[bl]) - sigma_est[bl])
        try:
            self.dlr_w_max = find_w_max(G_guess, self.dlr_eps)
        except RuntimeError:
            # The impurity G shifted by the Hartree self energy cannot be represented by a
            # DLR mesh whose nodes fit inside the current n_iw range. Fall back to the
            # largest mesh-fitting cutoff (best effort) and warn that n_iw limits accuracy.
            self.dlr_w_max = self._largest_fitting_dlr_w_max(G_guess)
            mpi.report('HARTREE SOLVER: warning, n_iw is too small to fully represent the '
                       'impurity Green function for the Hartree self energy shift; falling '
                       f'back to the largest mesh-fitting DLR cutoff w_max = {self.dlr_w_max:.4f}. '
                       'Increase n_iw for higher accuracy.')

        mpi.report(f'HARTREE SOLVER: using DLR mesh with w_max = {self.dlr_w_max:.4f}, eps = {self.dlr_eps:.1e}')

        # final solve on the properly sized mesh, restarting from the original guess
        # *************************************
        # this is done on every node due to very slow bcast
        self.triqs_solver.Sigma_HF = {bl: s.copy() for bl, s in sigma_hf_init.items()}
        self._set_dlr_G0(self.dlr_w_max)
        self.triqs_solver.solve(h_int=self.h_int, **self.triqs_solver_params)

        # call postprocessing
        self.postprocess()

        return

    def postprocess(self):
        r"""
        Organize G_freq, G_time, Sigma_freq and G_l from hartree solver
        """

        # get everything from solver, converting the solver's DLR Matsubara Green's
        # functions back onto the full Matsubara mesh used throughout solid_dmft
        n_iw = self.general_params['n_iw']
        self.G0_freq << make_gf_imfreq(self.triqs_solver.G0_iw, n_iw)
        self.G_freq_unsym << make_gf_imfreq(self.triqs_solver.G_iw, n_iw)
        self.G_freq << make_gf_imfreq(self.triqs_solver.G_iw, n_iw)
        self.sum_k.symm_deg_gf(self.G_freq, ish=self.icrsh)
        for bl, gf in self.Sigma_freq:
            self.Sigma_freq[bl] << self.triqs_solver.Sigma_HF[bl]
            self.Sigma_Refreq[bl] << self.triqs_solver.Sigma_HF[bl]
        self.G_time << Fourier(self.G_freq)
        self.interaction_energy = self.triqs_solver.interaction_energy()
        self.DC_energy = self.triqs_solver.DC_energy()

        return

