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
import numpy as np
from itertools import product

from triqs.gf import GfImTime, GfReTime, GfImFreq, GfReFreq, GfLegendre, BlockGf, make_hermitian, Omega, iOmega_n, make_gf_from_fourier, fit_hermitian_tail
from triqs.gf.tools import inverse, make_zero_tail
from triqs.gf.descriptors import Fourier
from triqs.operators import c_dag, c, Operator
import triqs.utility.mpi as mpi
from h5 import HDFArchive

from . import legendre_filter

def _gf_fit_tail_fraction(Gf, fraction=0.4, replace=None, known_moments=None):
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
        if known_moments:
            tail = Gf_fit[bl].fit_hermitian_tail(known_moments[i])
        else:
            tail = Gf_fit[bl].fit_hermitian_tail()
        nmax_frac = int(len(Gf_fit[bl].mesh)/2 * (1-replace))
        Gf_fit[bl].replace_by_tail(tail[0],n_min=nmax_frac)

    return Gf_fit

class SolverStructure:

    r'''
    Handles all solid_dmft solver objects and contains TRIQS solver instance.

    Attributes
    ----------

    Methods
    -------
    solve(self)
        solve impurity problem
    '''

    def __init__(self, general_params, solver_params, sum_k, icrsh, h_int, iteration_offset, solver_struct_ftps):
        r'''
        Initialisation of the solver instance with h_int for impurity "icrsh" based on soliDMFT parameters.

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
        '''

        self.general_params = general_params
        self.solver_params = solver_params
        self.sum_k = sum_k
        self.icrsh = icrsh
        self.h_int = h_int
        self.iteration_offset = iteration_offset
        self.solver_struct_ftps = solver_struct_ftps

        # initialize solver object, options are cthyb
        if self.general_params['solver_type'] == 'cthyb':
            from triqs_cthyb.version import triqs_cthyb_hash, version

            # sets up necessary GF objects on ImFreq
            self._init_ImFreq_objects()
            # sets up solver
            self.triqs_solver = self._create_cthyb_solver()
            self.git_hash = triqs_cthyb_hash
            self.version = version

        elif self.general_params['solver_type'] == 'ctint':
            from triqs_ctint.version import triqs_ctint_hash, version

            # sets up necessary GF objects on ImFreq
            self._init_ImFreq_objects()
            # sets up solver
            self.triqs_solver = self._create_ctint_solver()
            self.solver_params['measure_histogram'] = self.solver_params.pop('measure_pert_order')
            self.solver_params['use_double_insertion'] = self.solver_params.pop('move_double')
            self.git_hash = triqs_ctint_hash
            self.version = version

        elif self.general_params['solver_type'] == 'hubbardI':
            from triqs_hubbardI.version import triqs_hubbardI_hash, version

            # sets up necessary GF objects on ImFreq
            self._init_ImFreq_objects()
            self._init_ReFreq_hubbardI()
            # sets up solver
            self.triqs_solver = self._create_hubbardI_solver()
            self.git_hash = triqs_hubbardI_hash
            self.version = version

        elif self.general_params['solver_type'] == 'ftps':
            from forktps.version import forktps_hash, version

            # additional parameters
            self.bathfit_adjusted = self.iteration_offset != 0
            self.path_to_gs_accepted = bool(self.solver_params['path_to_gs'])
            self.convert_ftps = {'up': 'up', 'down': 'dn', 'ud': 'ud', 'ud_0': 'ud_0', 'ud_1': 'ud_1'}
            self.gf_struct = self.sum_k.gf_struct_solver_list[self.icrsh]
            for ct, block in enumerate(self.gf_struct):
                spin = block[0].split('_')[0] if not self.sum_k.corr_shells[self.icrsh]['SO'] else block[0]
                # FTPS solver does not know a more complicated gf_struct list of indices, so make sure the order is correct!
                indices = block[1] if not self.sum_k.corr_shells[self.icrsh]['SO'] else list(range(3))
                self.gf_struct[ct] = (self.convert_ftps[spin], indices)
            # sets up necessary GF objects on ReFreq
            self._init_ReFreq_objects()
            self.bathfit_adjusted = self.iteration_offset != 0
            self.path_to_gs_accepted = bool(self.solver_params['path_to_gs'])
            # sets up solver
            self.triqs_solver, self.sector_params, self.dmrg_params, self.tevo_params, self.calc_me, self.calc_mapping = self._create_ftps_solver()
            self.git_hash = forktps_hash
            self.version = version

        elif self.general_params['solver_type'] == 'inchworm':

            # sets up necessary GF objects on ImFreq
            self._init_ImFreq_objects()
            # sets up solver
            self.triqs_solver = self._create_inchworm_solver()
            # self.git_hash = inchworm_hash

        elif self.general_params['solver_type'] == 'ctseg':
            from triqs_ctseg.version import triqs_ctseg_hash, version

            # sets up necessary GF objects on ImFreq
            self._init_ImFreq_objects()
            # sets up solver
            self.triqs_solver = self._create_ctseg_solver()
            self.git_hash = triqs_ctseg_hash
            self.version = version

    # ********************************************************************
    # initialize Freq and Time objects
    # ********************************************************************

    def _init_ImFreq_objects(self):
        r'''
        Initialize all ImFreq objects
        '''

        # create all ImFreq instances
        self.n_iw = self.general_params['n_iw']
        self.G_freq = self.sum_k.block_structure.create_gf(ish=self.icrsh, gf_function=GfImFreq, space='solver',
                                                           beta=self.general_params['beta'], n_points=self.n_iw)
        # copy
        self.Sigma_freq = self.G_freq.copy()
        self.G0_freq = self.G_freq.copy()
        self.G_freq_unsym = self.G_freq.copy()
        self.Delta_freq = self.G_freq.copy()

        # create all ImTime instances
        self.n_tau = self.general_params['n_tau']
        self.G_time = self.sum_k.block_structure.create_gf(ish=self.icrsh, gf_function=GfImTime, space='solver',
                                                           beta=self.general_params['beta'], n_points=self.n_tau)
        # copy
        self.Delta_time = self.G_time.copy()

        # create all Legendre instances
        if (self.general_params['solver_type'] == 'cthyb' and self.solver_params['measure_G_l']
            or self.general_params['solver_type'] == 'cthyb' and  self.general_params['legendre_fit']
            or self.general_params['solver_type'] == 'ctseg' and self.solver_params['measure_gl']
            or self.general_params['solver_type'] == 'ctseg' and  self.general_params['legendre_fit']
            or self.general_params['solver_type'] == 'hubbardI' and self.solver_params['measure_G_l']):

            self.n_l = self.general_params['n_l']
            self.G_l = self.sum_k.block_structure.create_gf(ish=self.icrsh, gf_function=GfLegendre, space='solver',
                                                            beta=self.general_params['beta'], n_points=self.n_l)
            # move original G_freq to G_freq_orig
            self.G_time_orig = self.G_time.copy()

        if self.general_params['solver_type'] in ['cthyb', 'hubbardI'] and self.solver_params['measure_density_matrix']:
            self.density_matrix = None
            self.h_loc_diagonalization = None

        if self.general_params['solver_type'] in ['cthyb'] and self.general_params['measure_chi'] != 'none':
            self.O_time = None

    def _init_ReFreq_objects(self):
        r'''
        Initialize all ReFreq objects
        '''

        # create all ReFreq instances
        self.n_w = self.general_params['n_w']
        self.G_freq = self.sum_k.block_structure.create_gf(ish=self.icrsh, gf_function=GfReFreq, space='solver',
                                                           n_points=self.n_w, window = self.general_params['w_range'])
        # copy
        self.Sigma_freq = self.G_freq.copy()
        self.G0_freq = self.G_freq.copy()
        self.Delta_freq = self.G_freq.copy()
        self.G_freq_unsym = self.G_freq.copy()

        # create another Delta_freq for the solver, which uses different spin indices
        n_orb = self.sum_k.corr_shells[self.icrsh]['dim']
        n_orb = n_orb//2 if self.sum_k.corr_shells[self.icrsh]['SO'] else n_orb
        gf = GfReFreq(target_shape = (n_orb, n_orb), n_points = self.n_w, window = self.general_params['w_range'])

        self.Delta_freq_solver = BlockGf(name_list =tuple([block[0] for block in self.gf_struct]), block_list = (gf, gf), make_copies = True)

        # create all ReTime instances
        # FIXME: dummy G_time, since time_steps will be recalculated during run
        #time_steps = int(2 * self.solver_params['time_steps'] * self.solver_params['refine_factor']) if self.solver_params['n_bath'] != 0 else int(2 * self.solver_params['time_steps'])
        time_steps = int(2 * 1 * self.solver_params['refine_factor']) if self.solver_params['n_bath'] != 0 else int(2 * 1)
        self.G_time = self.sum_k.block_structure.create_gf(ish=self.icrsh, gf_function=GfReTime, space='solver', n_points=time_steps+1,
                                                           window=[0,time_steps*self.solver_params['dt']])

    def _init_ReFreq_hubbardI(self):
        r'''
        Initialize all ReFreq objects
        '''

        # create all ReFreq instances
        self.n_w = self.general_params['n_w']
        self.G_Refreq = self.sum_k.block_structure.create_gf(ish=self.icrsh, gf_function=GfReFreq, space='solver',
                                                             n_points=self.n_w, window = self.general_params['w_range'])
        # copy
        self.Sigma_Refreq = self.G_Refreq.copy()
        self.G0_Refreq = self.G_Refreq.copy()

    # ********************************************************************
    # solver-specific solve() command
    # ********************************************************************

    def solve(self):
        r'''
        solve impurity problem with current solver
        '''

        if self.general_params['solver_type'] == 'cthyb':

            if self.general_params['cthyb_delta_interface']:
                mpi.report('\n Using the delta interface for cthyb passing Delta(tau) and Hloc0 directly.')
                 # prepare solver input
                sumk_eal = self.sum_k.eff_atomic_levels()[self.icrsh]
                solver_eal = self.sum_k.block_structure.convert_matrix(sumk_eal, space_from='sumk', ish_from=self.sum_k.inequiv_to_corr[self.icrsh])
                # fill Delta_time from Delta_freq sum_k to solver
                for name, g0 in self.G0_freq:
                    self.Delta_freq[name] << iOmega_n - inverse(g0) - solver_eal[name]
                    known_moments = make_zero_tail(self.Delta_freq[name], 1)
                    tail, err = fit_hermitian_tail(self.Delta_freq[name], known_moments)
                    # without SOC delta_tau needs to be real
                    if not self.sum_k.SO == 1:
                        self.triqs_solver.Delta_tau[name] << make_gf_from_fourier(self.Delta_freq[name], self.triqs_solver.Delta_tau.mesh, tail).real
                    else:
                        self.triqs_solver.Delta_tau[name] << make_gf_from_fourier(self.Delta_freq[name], self.triqs_solver.Delta_tau.mesh, tail)


                # Make non-interacting operator for Hloc0
                Hloc_0 = Operator()
                for spin, spin_block in solver_eal.items():
                    for o1 in range(spin_block.shape[0]):
                        for o2 in range(spin_block.shape[1]):
                            # check if off-diag element is larger than threshold
                            if o1 != o2 and abs(spin_block[o1,o2]) < self.solver_params['off_diag_threshold']:
                                continue
                            else:
                                # TODO: adapt for SOC calculations, which should keep the imag part
                                Hloc_0 += spin_block[o1,o2].real/2 * (c_dag(spin,o1) * c(spin,o2) + c_dag(spin,o2) * c(spin,o1))
                self.solver_params['h_loc0'] = Hloc_0
            else:
                # fill G0_freq from sum_k to solver
                self.triqs_solver.G0_iw << self.G0_freq

            # update solver in h5 archive one last time for debugging if solve command crashes
            if self.general_params['store_solver'] and mpi.is_master_node():
                with HDFArchive(self.general_params['jobname']+'/'+self.general_params['seedname']+'.h5', 'a') as archive:
                    if not 'it_-1' in archive['DMFT_input/solver']:
                        archive['DMFT_input/solver'].create_group('it_-1')
                    archive['DMFT_input/solver/it_-1'][f'S_{self.icrsh}'] = self.triqs_solver
                    archive['DMFT_input/solver/it_-1'][f'solve_params_{self.icrsh}'] = self.solver_params
                    archive['DMFT_input/solver/it_-1']['mpi_size'] = mpi.size

            # Solve the impurity problem for icrsh shell
            # *************************************
            self.triqs_solver.solve(h_int=self.h_int, **self.solver_params)
            # *************************************

            # call postprocessing
            self._cthyb_postprocessing()

        elif self.general_params['solver_type'] == 'ctint':
            # fill G0_freq from sum_k to solver
            self.triqs_solver.G0_iw << self.G0_freq

            if self.general_params['h_int_type'] == 'dynamic':
                for b1, b2 in product(self.sum_k.gf_struct_solver_dict[self.icrsh].keys(), repeat=2):
                    self.triqs_solver.D0_iw[b1,b2] << self.U_iw[self.icrsh]

            # Solve the impurity problem for icrsh shell
            # *************************************
            self.triqs_solver.solve(h_int=self.h_int, **self.solver_params)
            # *************************************

            # call postprocessing
            self._ctint_postprocessing()

        elif self.general_params['solver_type'] == 'hubbardI':
            # fill G0_freq from sum_k to solver
            self.triqs_solver.G0_iw << self.G0_freq

            # Solve the impurity problem for icrsh shell
            # *************************************
            # this is done on every node due to very slow bcast of the AtomDiag object as of now
            self.triqs_solver.solve(h_int=self.h_int, calc_gtau=self.solver_params['measure_G_tau'],
                                    calc_gw=True, calc_gl=self.solver_params['measure_G_l'],
                                    calc_dm=self.solver_params['measure_density_matrix'])
            # if density matrix is measured, get this too. Needs to be done here,
            # because solver property 'dm' is not initialized/broadcastable
            if self.solver_params['measure_density_matrix']:
                self.density_matrix = self.triqs_solver.dm
                self.h_loc_diagonalization = self.triqs_solver.ad
            # *************************************

            # call postprocessing
            self._hubbardI_postprocessing()

        elif self.general_params['solver_type'] == 'ftps':
            import forktps as ftps
            from forktps.DiscreteBath import DiscretizeBath, TimeStepEstimation
            from forktps.BathFitting import BathFitter
            from forktps.Helpers import MakeGFstruct
            from . import OffDiagFitter as off_fitter

            def make_positive_definite(G):
                # ensure that Delta is positive definite
                for name, gf in G:
                    for orb, w in product(range(gf.target_shape[0]), gf.mesh):
                        if gf[orb,orb][w].imag > 0.0:
                            gf[orb,orb][w] = gf[orb,orb][w].real + 0.0j
                return G

            # create h_loc solver object
            h_loc = ftps.solver_core.Hloc(MakeGFstruct(self.Delta_freq_solver), SO=bool(self.sum_k.corr_shells[self.icrsh]['SO']))
            # need eff_atomic_levels
            sumk_eal = self.sum_k.eff_atomic_levels()[self.icrsh]

            # fill Delta_time from Delta_freq sum_k to solver
            for name, g0 in self.G0_freq:
                spin = name.split('_')[0] if not self.sum_k.corr_shells[self.icrsh]['SO'] else name
                ftps_name = self.convert_ftps[spin]
                solver_eal = self.sum_k.block_structure.convert_matrix(sumk_eal, space_from='sumk', ish_from=self.sum_k.inequiv_to_corr[self.icrsh])[name]
                self.Delta_freq[name] << Omega + 1j * self.general_params['eta'] - inverse(g0) - solver_eal
                # solver Delta is symmetrized by just using 'up_0' channel
                self.Delta_freq_solver[ftps_name] << Omega + 1j * self.general_params['eta'] - inverse(g0) - solver_eal

            # ensure that Delta is positive definite
            self.Delta_freq_solver = make_positive_definite(self.Delta_freq_solver)

            # remove off-diagonal terms
            if self.general_params['diag_delta']:
                for name, delta in self.Delta_freq_solver:
                    for i_orb, j_orb in product(range(delta.target_shape[0]),range(delta.target_shape[1])):
                        if i_orb != j_orb:
                            delta[i_orb,j_orb] << 0.0 + 0.0j

            # option to increase bath sites, but run with previous eta to get increased accuracy
            if self.solver_params['n_bath'] != 0 and self.solver_params['refine_factor'] != 1:
                if not self.bathfit_adjusted or self.bathfit_adjusted and self.iteration_offset > 0:
                    mpi.report('Rescaling "n_bath" with a factor of {}'.format(self.solver_params['refine_factor']))
                    self.solver_params['n_bath'] = int(self.solver_params['refine_factor']*self.solver_params['n_bath'])

            if self.solver_params['bath_fit']:

                # bathfitter
                # FIXME: this is temporary, since off-diagonal Bathfitter is not yet integrated in FTPS
                if self.sum_k.corr_shells[self.icrsh]['SO']:
                    fitter = off_fitter.OffDiagBathFitter(Nb=self.solver_params['n_bath']) if (self.solver_params['refine_factor'] != 1 and self.solver_params['n_bath'] != 0) else off_fitter.OffDiagBathFitter(Nb=None)
                    Delta_discrete = fitter.FitBath(Delta=self.Delta_freq_solver, eta=self.general_params['eta'], ignoreWeight=self.solver_params['ignore_weight'],
                                                    SO=bool(self.sum_k.corr_shells[self.icrsh]['SO']))
                else:
                    fitter = BathFitter(Nb=self.solver_params['n_bath']) if self.solver_params['n_bath'] != 0 else BathFitter(Nb=None)
                    Delta_discrete = fitter.FitBath(Delta=self.Delta_freq_solver, eta=self.general_params['eta'], ignoreWeight=self.solver_params['ignore_weight'])
            else:
                # discretizebath
                gap_interval = self.solver_params['enforce_gap'] if self.solver_params['enforce_gap'] != 'none' else None
                Delta_discrete = DiscretizeBath(Delta=self.Delta_freq_solver, Nb=self.solver_params['n_bath'], gap=gap_interval,
                                                SO=bool(self.sum_k.corr_shells[self.icrsh]['SO']))

            # should be done only once after the first iteration
            if self.solver_params['n_bath'] != 0 and self.solver_params['refine_factor'] != 1:
                if not self.bathfit_adjusted or self.bathfit_adjusted and self.iteration_offset > 0:
                    mpi.report('Rescaling "1/eta" with a factor of {}'.format(self.solver_params['refine_factor']))
                    # rescaling eta
                    self.general_params['eta'] /= self.solver_params['refine_factor']

                    if not self.bathfit_adjusted:
                        self.bathfit_adjusted = True

            self.triqs_solver.b = Delta_discrete
            # calculate time_steps
            time_steps = TimeStepEstimation(self.triqs_solver.b, eta=self.general_params['eta'], dt=self.solver_params['dt'])
            mpi.report('TimeStepEstimation returned {} with given bath, "eta" = {} and "dt" = {}'.format(time_steps, self.general_params['eta'],
                                                                                                         self.solver_params['dt']))
            # need to update tevo_params and G_time
            self.tevo_params.time_steps = time_steps
            self.G_time = self.sum_k.block_structure.create_gf(ish=self.icrsh, gf_function=GfReTime, space='solver', n_points=2*time_steps+1,
                                                               window=[0,2*time_steps*self.solver_params['dt']])


            # fill Hloc FTPS object
            # get hloc_dft from effective atomic levels
            for name, gf in self.Delta_freq:
                solver_eal = self.sum_k.block_structure.convert_matrix(sumk_eal, space_from='sumk', ish_from=self.sum_k.inequiv_to_corr[self.icrsh])[name]
                if not self.sum_k.corr_shells[self.icrsh]['SO']:
                    name = self.convert_ftps[name.split('_')[0]]
                    solver_eal = solver_eal.real
                    # remove off-diagonal terms
                    if self.general_params['diag_delta']:
                        solver_eal = np.diag(np.diag(solver_eal))
                h_loc.Fill(name, solver_eal)

            # fill solver h_loc
            self.triqs_solver.e0 = h_loc

            # FIXME: unfortunately, in the current implementation the solver initializations aren't included yet in dmft_cycle,
            # so for debugging it is done here again
            # store solver to h5 archive
            if self.general_params['store_solver'] and mpi.is_master_node():
                archive = HDFArchive(self.general_params['jobname']+'/'+self.general_params['seedname']+'.h5', 'a')
                archive['DMFT_input/solver'].create_group('it_-1')
                archive['DMFT_input/solver/it_-1']['Delta'] = self.Delta_freq_solver
                archive['DMFT_input/solver/it_-1']['S_'+str(self.icrsh)] = self.triqs_solver

            # Solve the impurity problem for icrsh shell
            # *************************************
            path_to_gs = self.solver_params['path_to_gs'] if self.solver_params['path_to_gs'] != 'none' and self.path_to_gs_accepted else None
            # fix to make sure this is only done in iteration 1
            if self.path_to_gs_accepted:
                self.path_to_gs_accepted = False
            if path_to_gs != 'postprocess':
                self.triqs_solver.solve(h_int=self.h_int, params_GS=self.dmrg_params, params_partSector=self.sector_params,
                                        tevo=self.tevo_params, eta=self.general_params['eta'], calc_me = self.calc_me,
                                        state_storage=self.solver_params['state_storage'],path_to_gs=path_to_gs)
            else:
                self.triqs_solver.post_process(h_int=self.h_int, params_GS=self.dmrg_params, params_partSector=self.dmrg_params,
                                               tevo=self.tevo_params, eta=self.general_params['eta'], calc_me = self.calc_me,
                                               state_storage=self.solver_params['state_storage'])
            # *************************************

            # call postprocessing
            self._ftps_postprocessing()

        elif self.general_params['solver_type'] == 'inchworm':
            # fill Delta_time from Delta_freq sum_k to solver
            self.triqs_solver.Delta_tau << make_gf_from_fourier(self.Delta_freq).real

            # Solve the impurity problem for icrsh shell
            # *************************************
            self.triqs_solver.solve(h_int=self.h_int, **self.solver_params)
            # *************************************

            # call postprocessing
            self._inchworm_postprocessing()

        if self.general_params['solver_type'] == 'ctseg':
            # fill G0_freq from sum_k to solver
            self.triqs_solver.G0_iw << self.G0_freq


            if self.general_params['h_int_type'] == 'dynamic':
                for b1, b2 in product(self.sum_k.gf_struct_solver_dict[self.icrsh].keys(), repeat=2):
                    self.triqs_solver.D0_iw[b1+"|"+b2] << self.U_iw[self.icrsh]

            # Solve the impurity problem for icrsh shell
            # *************************************
            self.triqs_solver.solve(h_int=self.h_int, **self.solver_params)
            # *************************************

            # call postprocessing
            self._ctseg_postprocessing()

        return

    # ********************************************************************
    # create solvers objects
    # ********************************************************************

    def _create_cthyb_solver(self):
        r'''
        Initialize cthyb solver instance
        '''
        from triqs_cthyb.solver import Solver as cthyb_solver

        gf_struct = self.sum_k.gf_struct_solver_list[self.icrsh]
        # Construct the triqs_solver instances
        if self.solver_params['measure_G_l']:
            triqs_solver = cthyb_solver(beta=self.general_params['beta'], gf_struct=gf_struct,
                            n_iw=self.general_params['n_iw'], n_tau=self.general_params['n_tau'],
                            n_l=self.general_params['n_l'], delta_interface=self.general_params['cthyb_delta_interface'])
        else:
            triqs_solver = cthyb_solver(beta=self.general_params['beta'], gf_struct=gf_struct,
                            n_iw=self.general_params['n_iw'], n_tau=self.general_params['n_tau'],
                            delta_interface=self.general_params['cthyb_delta_interface'])

        return triqs_solver

    def _create_ctint_solver(self):
        r'''
        Initialize ctint solver instance
        '''
        from triqs_ctint import Solver as ctint_solver

        gf_struct = self.sum_k.gf_struct_solver_list[self.icrsh]

        if self.general_params['h_int_type'] == 'dynamic':
            self.U_iw = None
            if  mpi.is_master_node():
                with HDFArchive(self.general_params['jobname']+'/'+self.general_params['seedname']+'.h5', 'r') as archive:
                    self.U_iw = archive['dynamic_U']['U_iw']
            self.U_iw = mpi.bcast(self.U_iw)
            n_iw_dyn = self.U_iw[self.icrsh].mesh.last_index()+1
            # Construct the triqs_solver instances
            triqs_solver = ctint_solver(beta=self.general_params['beta'], gf_struct=gf_struct,
                            n_iw=self.general_params['n_iw'], n_tau=self.general_params['n_tau'], use_D=True, use_Jperp=False,
                            n_iw_dynamical_interactions=n_iw_dyn,n_tau_dynamical_interactions=(int(n_iw_dyn*2.5)))
        else:
            # Construct the triqs_solver instances
            triqs_solver = ctint_solver(beta=self.general_params['beta'], gf_struct=gf_struct,
                            n_iw=self.general_params['n_iw'], n_tau=self.general_params['n_tau'], use_D=False, use_Jperp=False)

        return triqs_solver

    def _create_hubbardI_solver(self):
        r'''
        Initialize hubbardI solver instance
        '''
        from triqs_hubbardI import Solver as hubbardI_solver

        gf_struct =  self.sum_k.gf_struct_solver_list[self.icrsh]
        # Construct the triqs_solver instances
        if self.solver_params['measure_G_l']:
            triqs_solver = hubbardI_solver(beta=self.general_params['beta'], gf_struct=gf_struct,
                                           n_iw=self.general_params['n_iw'], n_tau=self.general_params['n_tau'],
                                           n_l=self.general_params['n_l'], n_w=self.general_params['n_w'],
                                           w_min=self.general_params['w_range'][0], w_max=self.general_params['w_range'][1],
                                           idelta=self.general_params['eta'])
        else:
            triqs_solver = hubbardI_solver(beta=self.general_params['beta'], gf_struct=gf_struct,
                                           n_iw=self.general_params['n_iw'], n_tau=self.general_params['n_tau'],
                                           n_w=self.general_params['n_w'], idelta=self.general_params['eta'],
                                           w_min=self.general_params['w_range'][0], w_max=self.general_params['w_range'][1])

        return triqs_solver

    def _create_inchworm_solver(self):
        r'''
        Initialize inchworm solver instance
        '''

        return triqs_solver

    def _create_ctseg_solver(self):
        r'''
        Initialize cthyb solver instance
        '''
        from triqs_ctseg import Solver as ctseg_solver

        if self.general_params['h_int_type'] == 'dynamic':
            self.U_iw = None
            if  mpi.is_master_node():
                with HDFArchive(self.general_params['jobname']+'/'+self.general_params['seedname']+'.h5', 'r') as archive:
                    self.U_iw = archive['dynamic_U']['U_iw']
            self.U_iw = mpi.bcast(self.U_iw)
            n_w_b_nn = self.U_iw[self.icrsh].mesh.last_index()+1
        else:
            n_w_b_nn = 1001

        gf_struct = self.sum_k.gf_struct_solver_list[self.icrsh]
        # Construct the triqs_solver instances
        if self.solver_params['measure_gl']:
            triqs_solver = ctseg_solver(beta=self.general_params['beta'], gf_struct=gf_struct,
                            n_iw=self.general_params['n_iw'], n_tau=self.general_params['n_tau'],
                            n_legendre_g=self.general_params['n_l'], n_w_b_nn = n_w_b_nn, n_tau_k=int(n_w_b_nn*2.5), n_tau_jperp=int(n_w_b_nn*2.5))
        else:
            triqs_solver = ctseg_solver(beta=self.general_params['beta'], gf_struct=gf_struct,
                            n_iw=self.general_params['n_iw'], n_tau=self.general_params['n_tau'],
                            n_w_b_nn = n_w_b_nn, n_tau_k=int(n_w_b_nn*2.5), n_tau_jperp=int(n_w_b_nn*2.5))

        return triqs_solver

    def _create_ftps_solver(self):
        r'''
        Initialize ftps solver instance
        '''
        import forktps as ftps

        # convert self.solver_struct_ftps to mapping and solver-friendly list
        if not self.sum_k.corr_shells[self.icrsh]['SO']:
            # mapping dictionary
            calc_mapping = {self.solver_struct_ftps[self.icrsh][deg_shell][0]:
                    self.solver_struct_ftps[self.icrsh][deg_shell][1:] for deg_shell in range(len(self.solver_struct_ftps[self.icrsh]))}
            # make solver-friendly list from mapping keys
            calc_me = [[item.split('_')[0], int(item.split('_')[1])] for item in calc_mapping.keys()]
            # replace 'down' with 'dn'
            calc_me = [[item[0].replace('down','dn'),item[1]] for item in calc_me]
        else:
            # for SOC we just end up calculating everything for now
            # TODO: perhaps skip down channel
            calc_mapping = None
            calc_me = [[f'ud_{i}',j] for i,j in product(range(2), range(3))]

        # create solver
        triqs_solver = ftps.Solver(gf_struct=self.gf_struct, nw=self.general_params['n_w'],
                                   wmin=self.general_params['w_range'][0], wmax=self.general_params['w_range'][1])


        # create partSector params
        sector_params = ftps.solver.DMRGParams(maxmI=50, maxmIB=50, maxmB=50, tw=1e-10, nmax=5, sweeps=5)

        # for now prep_imagTevo, prep_method and nmax hard-coded
        # create DMRG params
        dmrg_params = ftps.solver.DMRGParams(maxmI=self.solver_params['dmrg_maxmI'], maxmIB=self.solver_params['dmrg_maxmIB'],
                                             maxmB=self.solver_params['dmrg_maxmB'], tw=self.solver_params['dmrg_tw'],
                                             prep_imagTevo=True, prep_method='TEBD', sweeps=self.solver_params['sweeps'], nmax=2,
                                             prep_time_steps=5, napph=2
                                             )

        # create TEVO params
        tevo_params = ftps.solver.TevoParams(dt=self.solver_params['dt'], time_steps=1, #dummy, will be updated during the run
                                             maxmI=self.solver_params['maxmI'], maxmIB=self.solver_params['maxmIB'],
                                             maxmB=self.solver_params['maxmB'], tw=self.solver_params['tw'])

        return triqs_solver, sector_params, dmrg_params, tevo_params, calc_me, calc_mapping

    # ********************************************************************
    # post-processing of solver output
    # ********************************************************************

    def _cthyb_postprocessing(self):
        r'''
        Organize G_freq, G_time, Sigma_freq and G_l from cthyb solver
        '''

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

            if self.general_params['legendre_fit']:
                self.G_time_orig << self.triqs_solver.G_tau
                # run the filter
                self.G_l << legendre_filter.apply(self.G_time, self.general_params['n_l'])
                # get G_time, G_freq, Sigma_freq from G_l
                set_Gs_from_G_l()
            elif self.solver_params['perform_tail_fit'] and not self.general_params['legendre_fit']:
                # if tailfit has been used replace Sigma with the tail fitted Sigma from cthyb
                self.Sigma_freq << self.triqs_solver.Sigma_iw
                self.sum_k.symm_deg_gf(self.Sigma_freq, ish=self.icrsh)
            else:
                # obtain Sigma via dyson from symmetrized G_freq
                self.Sigma_freq << inverse(self.G0_freq) - inverse(self.G_freq)

        # if density matrix is measured, get this too
        if self.solver_params['measure_density_matrix']:
            self.density_matrix = self.triqs_solver.density_matrix
            self.h_loc_diagonalization = self.triqs_solver.h_loc_diagonalization

        if self.solver_params['measure_pert_order']:
            self.perturbation_order = self.triqs_solver.perturbation_order
            self.perturbation_order_total = self.triqs_solver.perturbation_order_total

        if self.general_params['measure_chi'] != 'none':
            self.O_time = self.triqs_solver.O_tau

        return

    def _ctint_postprocessing(self):
        r'''
        Organize G_freq, G_time, Sigma_freq and G_l from cthyb solver
        '''
        #TODO

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
        # if self.general_params['legendre_fit']:
        #     self.G_freq_orig << self.triqs_solver.G_iw
        #     # run the filter
        #     self.G_l << legendre_filter.apply(self.G_time, self.general_params['n_l'])
        #     # get G_time, G_freq, Sigma_freq from G_l
        #     set_Gs_from_G_l()

        if self.solver_params['measure_histogram']:
            self.perturbation_order = self.triqs_solver.histogram

        return

    def _hubbardI_postprocessing(self):
        r'''
        Organize G_freq, G_time, Sigma_freq and G_l from hubbardI solver
        '''

        # get everything from solver
        self.Sigma_freq << self.triqs_solver.Sigma_iw
        self.G0_freq << self.triqs_solver.G0_iw
        self.G0_Refreq << self.triqs_solver.G0_w
        self.G_freq << self.triqs_solver.G_iw
        self.sum_k.symm_deg_gf(self.G_freq, ish=self.icrsh)
        self.G_freq_unsym << self.G_freq
        self.G_Refreq << self.triqs_solver.G_w
        self.Sigma_Refreq << self.triqs_solver.Sigma_w

        # if measured in Legendre basis, get G_l from solver too
        if self.solver_params['measure_G_l']:
            self.G_l << self.triqs_solver.G_l

        if self.solver_params['measure_G_tau']:
            self.G_time << self.triqs_solver.G_tau

        return

    def _inchworm_postprocessing(self):
        r'''
        Organize G_freq, G_time, Sigma_freq and G_l from inchworm solver
        '''

        return

    def _ftps_postprocessing(self):
        r'''
        Organize G_freq, G_time, Sigma_freq and G_l from ftps solver
        '''
        from forktps.DiscreteBath import SigmaDyson

        # symmetrization of reduced solver G
        def symmetrize_opt(G_in, soc):
            G = G_in.copy()
            if soc:
                def swap_2():
                    for i in range(2):
                        G['ud_1'][i,2] = -G['ud_1'][i,2]
                        G['ud_1'][2,i] = -G['ud_1'][2,i]
                swap_2()
                G['ud_0'] = 0.5*(G['ud_0'] + G['ud_1'])
                G['ud_1'] = G['ud_0']
                for name , g in G:
                    g[1,1] = 0.5*(g[1,1]+g[2,2])
                    g[2,2] = g[1,1]
                swap_2()
            else:
                switch = lambda spin: 'dn' if spin == 'down' else 'up'
                for key, mapto in self.calc_mapping.items():
                    spin, block = key.split('_')
                    for deg_item in mapto:
                        map_spin, map_block = deg_item.split('_')
                        mpi.report(f'mapping {spin}-{block} to {map_spin}-{map_block}...')
                        G[switch(map_spin)].data[:,int(map_block),int(map_block)] = G[switch(spin)].data[:,int(block),int(block)]
                # particle-hole symmetry: enforce mirror/point symmetry of G(w)
                if self.solver_params['ph_symm']:
                    for block, gf in G:
                        gf.data.real = 0.5 * ( gf.data[::1].real - gf.data[::-1].real )
                        gf.data.imag = 0.5 * ( gf.data[::1].imag + gf.data[::-1].imag )
            return G

        def symmetrize(G):
            return symmetrize_opt(G, soc=self.sum_k.corr_shells[self.icrsh]['SO'])

        def make_positive_definite(G):
            # ensure that Delta is positive definite
            for name, gf in G:
                for orb, w in product(range(gf.target_shape[0]), gf.mesh):
                    if gf[orb,orb][w].imag > 0.0:
                        gf[orb,orb][w] = gf[orb,orb][w].real + 0.0j
            return G

        G_w = symmetrize(self.triqs_solver.G_w)
        if not self.sum_k.corr_shells[self.icrsh]['SO']:
            G_w = make_positive_definite(G_w)

        # calculate Sigma_freq via Dyson
        # do not use Dyson equation directly, as G0 might have wrong eta
        Sigma_w_symm = SigmaDyson(Gret=self.triqs_solver.G_ret, bath=self.triqs_solver.b,
                                  hloc=self.triqs_solver.e0, mesh=self.Delta_freq_solver.mesh,
                                  eta=self.general_params['eta'], symmG=symmetrize)

        # convert everything to solver objects
        for block, gf in G_w:
            if not self.sum_k.corr_shells[self.icrsh]['SO']:
                reverse_convert = dict(map(reversed, self.convert_ftps.items()))
                sumk_name = reverse_convert[block.split('_')[0]] + '_0'
            else:
                sumk_name = block
            self.G_freq[sumk_name] << gf
            # in FTPS the unsym result is not calculated. Symmetries are used by construction
            self.G_freq_unsym[sumk_name] << gf
            self.Sigma_freq[sumk_name] << Sigma_w_symm[block]
            self.G_time[sumk_name] << self.triqs_solver.G_ret[block]

        return


    def _ctseg_postprocessing(self):
        r'''
        Organize G_freq, G_time, Sigma_freq and G_l from cthyb solver
        '''

        def set_Gs_from_G_l():

            if self.solver_params['measure_ft'] and mpi.is_master_node():
                print('\n !!!!WARNING!!!! \n you enabled both improved estimators and legendre based filtering / sampling. Sigma will be overwritten by legendre result.  \n !!!!WARNING!!!!\n')

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

        # first print average sign
        if mpi.is_master_node():
            print('\nAverage sign: {}'.format(self.triqs_solver.average_sign))
        # get Delta_time from solver
        self.Delta_time << self.triqs_solver.Delta_tau

        self.G_time << self.triqs_solver.G_tau
        self.sum_k.symm_deg_gf(self.G_time, ish=self.icrsh)

        if self.solver_params['measure_gw']:
            self.G_freq << self.triqs_solver.G_iw
        else:
            if mpi.is_master_node():
                # create empty moment container (list of np.arrays)
                Gf_known_moments = make_zero_tail(self.G_freq,n_moments=2)
                for i, bl in enumerate(self.G_freq.indices):
                    # 0 moment is 0, dont touch it, but first moment is 1 for the Gf
                    Gf_known_moments[i][1] = np.eye(self.G_freq[bl].target_shape[0])
                    self.G_freq[bl] << Fourier(self.G_time[bl], Gf_known_moments[i])
            self.G_freq << mpi.bcast(self.G_freq)

        self.G_freq << make_hermitian(self.G_freq)
        self.G_freq_unsym << self.G_freq
        self.sum_k.symm_deg_gf(self.G_freq, ish=self.icrsh)

        # if measured in Legendre basis, get G_l from solver too
        if self.solver_params['measure_gl']:
            # store original G_time into G_time_orig
            self.G_time_orig << self.triqs_solver.G_tau
            self.G_l << self.triqs_solver.G_l
            # get G_time, G_freq, Sigma_freq from G_l
            set_Gs_from_G_l()
        elif self.general_params['legendre_fit']:
            self.G_time_orig << self.triqs_solver.G_tau
            self.G_l << legendre_filter.apply(self.G_time, self.general_params['n_l'])
            # get G_time, G_freq, Sigma_freq from G_l
            set_Gs_from_G_l()
        # if improved estimators are turned on calc Sigma from F_tau, otherwise:
        elif self.solver_params['measure_ft']:
            self.F_freq = self.G_freq.copy()
            self.F_freq << 0.0
            self.F_time = self.G_time.copy()
            self.F_time << self.triqs_solver.F_tau
            F_known_moments = make_zero_tail(self.F_freq, n_moments=1)
            if mpi.is_master_node():
                for i, bl in enumerate(self.F_freq.indices):
                    self.F_freq[bl] << Fourier(self.triqs_solver.F_tau[bl], F_known_moments[i])
                # fit tail of improved estimator and G_freq
                self.F_freq << Gf_fit_tail_fraction(self.F_freq, fraction=0.9, replace=0.5, known_moments=F_known_moments)
                self.G_freq << Gf_fit_tail_fraction(self.G_freq ,fraction=0.9, replace=0.5, known_moments=Gf_known_moments)

            self.F_freq << mpi.bcast(self.F_freq)
            self.G_freq << mpi.bcast(self.G_freq)
            for block, fw in self.F_freq:
                for iw in fw.mesh:
                    self.Sigma_freq[block][iw] = self.F_freq[block][iw] / self.G_freq[block][iw]

        else:
            mpi.report('\n!!!! WARNING !!!! tail of solver output not handled! Turn on either measure_ft, legendre_fit, or measure_gl\n')
            self.Sigma_freq << inverse(self.G0_freq) - inverse(self.G_freq)


        if self.solver_params['measure_hist']:
            self.perturbation_order = self.triqs_solver.histogram

        return
