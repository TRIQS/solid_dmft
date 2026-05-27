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
abstract DMFT solver class for solid_dmft
'''
import numpy as np

from triqs.gfs import MeshImTime, MeshReTime, MeshReFreq, MeshLegendre, Gf, BlockGf
from triqs.operators import Operator

# import abc
from abc import ABC, abstractmethod

# circular import of all the solver subclasses


class AbstractDMFTSolver(ABC):
    """
    Abstract base class for DMFT solvers

    This class defines the template for solvers for solving the DMFT impurity problem for the icrsh impurity.

    Parameters
    ----------
    general_paramuters: dict
                        general parameters as dict
    solver_params: dict
                        solver-specific parameters as dict
    sum_k: triqs.dft_tools.sumk object
            SumkDFT instance
    icrsh: int
            correlated shell index
    h_int: triqs.operator object
            interaction Hamiltonian of correlated shell
    iteration_offset: int
            number of iterations this run is based on
    """

    # ********************************************************************
    # General initialization common to all solvers
    # ********************************************************************

    def __init__(
        self, general_params, solver_params, sum_k, icrsh, h_int, iteration_offset, deg_orbs_ftps, gw_params=None, advanced_params=None
    ):
        self.general_params = general_params
        self.solver_params = solver_params
        self.gw_params = gw_params
        self.advanced_params = advanced_params
        self.sum_k = sum_k
        self.icrsh = icrsh
        self.h_int = h_int
        self.iteration_offset = iteration_offset
        self.deg_orbs_ftps = deg_orbs_ftps

        if self.solver_params.get('crm_dyson_solver'):
            self.G_time_dlr = None
            self.Sigma_dlr = None

        if self.solver_params.get('measure_density_matrix'):
            self.density_matrix = None
            self.h_loc_diagonalization = None

        if self.solver_params.get('measure_chi') is not None:
            self.O_time = None

        if self.solver_params.get('delta_interface'):
            self.Hloc_0 = Operator()

    # ********************************************************************
    # initialize Freq and Time objects
    # ********************************************************************

    def _init_ImFreq_objects(self):
        r"""
        Initialize all ImFreq objects
        """

        # create all ImFreq instances
        self.n_iw = self.general_params['n_iw']
        self.G_freq = self.sum_k.block_structure.create_gf(ish=self.icrsh, gf_function=Gf, space='solver', mesh=self.sum_k.mesh)
        # copy
        self.Sigma_freq = self.G_freq.copy()
        self.G0_freq = self.G_freq.copy()
        self.G_freq_unsym = self.G_freq.copy()
        self.Delta_freq = self.G_freq.copy()

        # create all ImTime instances
        self.n_tau = self.general_params['n_tau']
        self.G_time = self.sum_k.block_structure.create_gf(
            ish=self.icrsh, gf_function=Gf, space='solver', mesh=MeshImTime(beta=self.sum_k.mesh.beta, statistic='Fermion', n_tau=self.n_tau)
        )
        # copy
        self.Delta_time = self.G_time.copy()

        # create all Legendre instances
        if self.solver_params.get('measure_G_l') or self.solver_params.get('legendre_fit'):
            self.n_l = self.solver_params['n_l']
            self.G_l = self.sum_k.block_structure.create_gf(
                ish=self.icrsh,
                gf_function=Gf,
                space='solver',
                mesh=MeshLegendre(beta=self.general_params['beta'], max_n=self.n_l, statistic='Fermion'),
            )
            # move original G_freq to G_freq_orig
            self.G_time_orig = self.G_time.copy()

    def _init_ReFreq_objects(self):
        r"""
        Initialize all ReFreq objects
        """

        # create all ReFreq instances
        self.n_w = self.general_params['n_w']
        self.G_freq = self.sum_k.block_structure.create_gf(ish=self.icrsh, gf_function=Gf, space='solver', mesh=self.sum_k.mesh)
        # copy
        self.Sigma_freq = self.G_freq.copy()
        self.G0_freq = self.G_freq.copy()
        self.Delta_freq = self.G_freq.copy()
        self.G_freq_unsym = self.G_freq.copy()

        # create another Delta_freq for the solver, which uses different spin indices
        n_orb = self.sum_k.corr_shells[self.icrsh]['dim']
        n_orb = n_orb // 2 if self.sum_k.corr_shells[self.icrsh]['SO'] else n_orb
        gf = Gf(target_shape=(n_orb, n_orb), mesh=MeshReFreq(n_w=self.n_w, window=self.general_params['w_range']))

        self.Delta_freq_solver = BlockGf(name_list=tuple([block[0] for block in self.gf_struct]), block_list=(gf, gf), make_copies=True)

        # create all ReTime instances
        # FIXME: dummy G_time, since time_steps will be recalculated during run
        # time_steps = int(2 * self.solver_params['time_steps'] * self.solver_params['refine_factor']) if self.solver_params['n_bath'] != 0 else int(2 * self.solver_params['time_steps'])
        time_steps = int(2 * 1 * self.solver_params['refine_factor']) if self.solver_params['n_bath'] != 0 else int(2 * 1)
        self.G_time = self.sum_k.block_structure.create_gf(
            ish=self.icrsh,
            gf_function=Gf,
            space='solver',
            mesh=MeshReTime(n_t=time_steps + 1, window=[0, time_steps * self.solver_params['dt']]),
        )

    # define abstract solve and postprocess methods

    @abstractmethod
    def solve(self, **kwargs):
        r"""
        Solve the DMFT impurity problem
        """
        pass

    @abstractmethod
    def postprocess(self):
        r"""
        Postprocess the DMFT impurity problem
        """
        pass

    @staticmethod
    def _gf_fit_tail_fraction(Gf, fraction=0.4, replace=None, known_moments=[]):
        """
        fits the tail of Gf object by making a polynomial
        fit of the Gf on the given fraction of the Gf mesh
        and replacing that part of the Gf by the fit

        0.4 fits the last 40% of the Gf and replaces the
        part with the tail

        Parameters
        ----------
        Gf : BlockGf (Green's function) object
        fraction: float, optional default 0.4
            fraction of the Gf to fit
        replace: float, optional default fraction
            fraction of the Gf to replace
        known_moments: np.array
            known moments as numpy array
        Returns
        -------
        Gf_fit : BlockGf (Green's function) object
                fitted Gf
        """

        Gf_fit = Gf.copy()
        # if no replace factor is given use the same fraction
        if not replace:
            replace = fraction

        for i, bl in enumerate(Gf_fit.indices):
            Gf_fit[bl].mesh.set_tail_fit_parameters(tail_fraction=fraction)
            if known_moments == []:
                tail = Gf_fit[bl].fit_hermitian_tail()
            else:
                tail = Gf_fit[bl].fit_hermitian_tail(known_moments[i])
            nmax_frac = int(len(Gf_fit[bl].mesh) / 2 * (1 - replace))
            Gf_fit[bl].replace_by_tail(tail[0], n_min=nmax_frac)

        return Gf_fit

    @staticmethod
    def _fit_tail_window(
        Sigma_iw, fit_min_n=None, fit_max_n=None, fit_min_w=None, fit_max_w=None, fit_max_moment=None, fit_known_moments=None
    ):
        """
        Fit a high frequency 1/(iw)^n expansion of Sigma_iw
        and replace the high frequency part with the fitted high frequency expansion.

        Either give frequency window to fit on in terms of matsubara frequencies index
        (fit_min_n/fit_max_n) or value (fit_min_w/fit_max_w).

        Parameters
        ----------
        Sigma_iw : Gf
                Self-energy.
        fit_min_n : int, optional, default=int(0.8*len(Sigma_iw.mesh))
                    Matsubara frequency index from which tail fitting should start.
        fit_max_n : int, optional, default=int(len(Sigma_iw.mesh))
                    Matsubara frequency index at which tail fitting should end.
        fit_min_w : float, optional
                    Matsubara frequency from which tail fitting should start.
        fit_max_w : float, optional
                    Matsubara frequency at which tail fitting should end.
        fit_max_moment : int, optional
                        Highest moment to fit in the tail of Sigma_iw.
        fit_known_moments : ``ndarray.shape[order, Sigma_iw[0].target_shape]``, optional, default = None
                            Known moments of Sigma_iw, given as an numpy ndarray

        Returns
        -------
        tail_barr : dict of arr
                    fitted tail of Sigma_iw
        """
        from triqs.gfs.gf_fnt import fit_hermitian_tail_on_window, replace_by_tail

        # Define default tail quantities
        if fit_min_w is not None:
            fit_min_n = int(0.5 * (fit_min_w * Sigma_iw.mesh.beta / np.pi - 1.0))
        if fit_max_w is not None:
            fit_max_n = int(0.5 * (fit_max_w * Sigma_iw.mesh.beta / np.pi - 1.0))
        if fit_min_n is None:
            fit_min_n = int(0.8 * len(Sigma_iw.mesh) / 2)
        if fit_max_n is None:
            fit_max_n = int(len(Sigma_iw.mesh) / 2)
        if fit_max_moment is None:
            fit_max_moment = 3

        if fit_known_moments is None:
            fit_known_moments = {}
            for name, sig in Sigma_iw:
                shape = [0] + list(sig.target_shape)
                fit_known_moments[name] = np.zeros(shape, dtype=complex)  # no known moments

        # Now fit the tails of Sigma_iw and replace the high frequency part with the tail expansion
        tail_barr = {}
        Sigma_fit = Sigma_iw.copy()
        for name, sig in Sigma_fit:
            tail, err = fit_hermitian_tail_on_window(
                sig,
                n_min=fit_min_n,
                n_max=fit_max_n,
                known_moments=fit_known_moments[name],
                # set max number of pts used in fit larger than mesh size, to use all data in fit
                n_tail_max=10 * len(sig.mesh),
                expansion_order=fit_max_moment,
            )
            tail_barr[name] = tail
            replace_by_tail(sig, tail, n_min=fit_min_n)

        return Sigma_fit, tail_barr

    def _make_spin_equal(self, Sigma):
        # if not SOC than average up and down
        if not self.general_params['magnetic'] and not self.sum_k.SO == 1:
            Sigma['up_0'] = 0.5 * (Sigma['up_0'] + Sigma['down_0'])
            Sigma['down_0'] = Sigma['up_0']

        return Sigma
