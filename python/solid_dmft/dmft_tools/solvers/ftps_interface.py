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
FTPS solver class for solid_dmft
'''
import numpy as np
from itertools import product

from triqs.gfs import MeshReTime, Gf, Omega
from triqs.gfs.tools import inverse
import triqs.utility.mpi as mpi
from h5 import HDFArchive

# import of the abstract class
from solid_dmft.dmft_tools.solvers.abstractdmftsolver import AbstractDMFTSolver

# import triqs solver
from forktps.version import forktps_hash, version
import forktps as ftps
from forktps.DiscreteBath import SigmaDyson, DiscretizeBath, TimeStepEstimation
from forktps.BathFitting import BathFitter
from forktps.Helpers import MakeGFstruct


class FTPSInterface(AbstractDMFTSolver):
    def __init__(
        self, general_params, solver_params, sum_k, icrsh, h_int, iteration_offset, deg_orbs_ftps, gw_params=None, advanced_params=None
    ):
        # Call the base class constructor
        super().__init__(general_params, solver_params, sum_k, icrsh, h_int, iteration_offset, deg_orbs_ftps, gw_params, advanced_params)

        # Create the  solver specifics
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

        # set up the solver
        # one needs to set up the following:
        # * triqs_solver
        # * sector_params
        # * dmrg_params
        # * tevo_params
        # * calc_me
        # * calc_mapping

        # TODO: add triqs_solver_params for ftps as well to make it analogous to other similars
        # Not necessary but nicer. For now, just keep an empty dummy dict
        self.triqs_solver_params = {}

        # convert self.deg_orbs_ftps to mapping and solver-friendly list
        if not self.sum_k.corr_shells[self.icrsh]['SO']:
            # mapping dictionary
            self.calc_mapping = {
                self.deg_orbs_ftps[self.icrsh][deg_shell][0]: self.deg_orbs_ftps[self.icrsh][deg_shell][1:]
                for deg_shell in range(len(self.deg_orbs_ftps[self.icrsh]))
            }
            # make solver-friendly list from mapping keys
            self.calc_me = [[item.split('_')[0], int(item.split('_')[1])] for item in self.calc_mapping.keys()]
            # replace 'down' with 'dn'
            self.calc_me = [[item[0].replace('down', 'dn'), item[1]] for item in self.calc_me]
        else:
            # for SOC we just end up calculating everything for now
            # TODO: perhaps skip down channel
            self.calc_mapping = None
            self.calc_me = [[f'ud_{i}', j] for i, j in product(range(2), range(3))]

        # create solver
        self.triqs_solver = ftps.Solver(
            gf_struct=self.gf_struct,
            nw=self.general_params['n_w'],
            wmin=self.general_params['w_range'][0],
            wmax=self.general_params['w_range'][1],
        )

        # create partSector params
        self.sector_params = ftps.solver.DMRGParams(maxmI=50, maxmIB=50, maxmB=50, tw=1e-10, nmax=5, sweeps=5)

        # for now prep_imagTevo, prep_method and nmax hard-coded
        # create DMRG params
        self.dmrg_params = ftps.solver.DMRGParams(
            maxmI=self.solver_params['dmrg_maxmI'],
            maxmIB=self.solver_params['dmrg_maxmIB'],
            maxmB=self.solver_params['dmrg_maxmB'],
            tw=self.solver_params['dmrg_tw'],
            prep_imagTevo=True,
            prep_method='TEBD',
            sweeps=self.solver_params['sweeps'],
            nmax=2,
            prep_time_steps=5,
            napph=2,
        )

        # create TEVO params
        self.tevo_params = ftps.solver.TevoParams(
            dt=self.solver_params['dt'],
            time_steps=1,  # dummy, will be updated during the run
            maxmI=self.solver_params['maxmI'],
            maxmIB=self.solver_params['maxmIB'],
            maxmB=self.solver_params['maxmB'],
            tw=self.solver_params['tw'],
        )

        self.git_hash = forktps_hash
        self.version = version

        return

    def solve(self, **kwargs):
        def make_positive_definite(G):
            # ensure that Delta is positive definite
            for name, gf in G:
                for orb, w in product(range(gf.target_shape[0]), gf.mesh):
                    if gf[orb, orb][w].imag > 0.0:
                        gf[orb, orb][w] = gf[orb, orb][w].real + 0.0j
            return G

        # create h_loc solver object
        h_loc = ftps.solver_core.Hloc(MakeGFstruct(self.Delta_freq_solver), SO=bool(self.sum_k.corr_shells[self.icrsh]['SO']))
        # need eff_atomic_levels
        sumk_eal = self.sum_k.eff_atomic_levels()[self.icrsh]

        # fill Delta_time from Delta_freq sum_k to solver
        for name, g0 in self.G0_freq:
            spin = name.split('_')[0] if not self.sum_k.corr_shells[self.icrsh]['SO'] else name
            ftps_name = self.convert_ftps[spin]
            solver_eal = self.sum_k.block_structure.convert_matrix(
                sumk_eal, space_from='sumk', ish_from=self.sum_k.inequiv_to_corr[self.icrsh]
            )[name]
            self.Delta_freq[name] << Omega + 1j * self.general_params['eta'] - inverse(g0) - solver_eal
            # solver Delta is symmetrized by just using 'up_0' channel
            self.Delta_freq_solver[ftps_name] << Omega + 1j * self.general_params['eta'] - inverse(g0) - solver_eal

        # ensure that Delta is positive definite
        self.Delta_freq_solver = make_positive_definite(self.Delta_freq_solver)

        if self.general_params['store_solver'] and mpi.is_master_node():
            archive = HDFArchive(self.general_params['jobname'] + '/' + self.general_params['seedname'] + '.h5', 'a')
            if not 'it_-1' in archive['DMFT_input/solver']:
                archive['DMFT_input/solver'].create_group('it_-1')
            archive['DMFT_input/solver/it_-1']['Delta_orig'] = self.Delta_freq_solver

        # remove off-diagonal terms
        if self.solver_params['diag_delta']:
            for name, delta in self.Delta_freq_solver:
                for i_orb, j_orb in product(range(delta.target_shape[0]), range(delta.target_shape[1])):
                    if i_orb != j_orb:
                        delta[i_orb, j_orb] << 0.0 + 0.0j

        # option to increase bath sites, but run with previous eta to get increased accuracy
        if self.solver_params['n_bath'] != 0 and self.solver_params['refine_factor'] != 1:
            if not self.bathfit_adjusted or self.bathfit_adjusted and self.iteration_offset > 0:
                mpi.report('Rescaling "n_bath" with a factor of {}'.format(self.solver_params['refine_factor']))
                self.solver_params['n_bath'] = int(self.solver_params['refine_factor'] * self.solver_params['n_bath'])

        if self.solver_params['bath_fit']:
            # bathfitter
            # FIXME: this is temporary, since off-diagonal Bathfitter is not yet integrated in FTPS
            if self.sum_k.corr_shells[self.icrsh]['SO']:
                fitter = (
                    off_fitter.OffDiagBathFitter(Nb=self.solver_params['n_bath'])
                    if (self.solver_params['refine_factor'] != 1 and self.solver_params['n_bath'] != 0)
                    else off_fitter.OffDiagBathFitter(Nb=None)
                )
                Delta_discrete = fitter.FitBath(
                    Delta=self.Delta_freq_solver,
                    eta=self.general_params['eta'],
                    ignoreWeight=self.solver_params['ignore_weight'],
                    SO=bool(self.sum_k.corr_shells[self.icrsh]['SO']),
                )
            else:
                fitter = BathFitter(Nb=self.solver_params['n_bath']) if self.solver_params['n_bath'] != 0 else BathFitter(Nb=None)
                Delta_discrete = fitter.FitBath(
                    Delta=self.Delta_freq_solver, eta=self.general_params['eta'], ignoreWeight=self.solver_params['ignore_weight']
                )
        else:
            # discretizebath
            gap_interval = self.solver_params['enforce_gap'] if self.solver_params['enforce_gap'] is not None else None
            Delta_discrete = DiscretizeBath(
                Delta=self.Delta_freq_solver,
                Nb=self.solver_params['n_bath'],
                gap=gap_interval,
                SO=bool(self.sum_k.corr_shells[self.icrsh]['SO']),
            )

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
        mpi.report(
            'TimeStepEstimation returned {} with given bath, "eta" = {} and "dt" = {}'.format(
                time_steps, self.general_params['eta'], self.solver_params['dt']
            )
        )
        # need to update tevo_params and G_time
        self.tevo_params.time_steps = time_steps
        self.G_time = self.sum_k.block_structure.create_gf(
            ish=self.icrsh,
            gf_function=Gf,
            space='solver',
            mesh=MeshReTime(n_t=2 * time_steps + 1, window=[0, 2 * time_steps * self.solver_params['dt']]),
        )

        # fill Hloc FTPS object
        # get hloc_dft from effective atomic levels
        for name, gf in self.Delta_freq:
            solver_eal = self.sum_k.block_structure.convert_matrix(
                sumk_eal, space_from='sumk', ish_from=self.sum_k.inequiv_to_corr[self.icrsh]
            )[name]
            if not self.sum_k.corr_shells[self.icrsh]['SO']:
                name = self.convert_ftps[name.split('_')[0]]
                solver_eal = solver_eal.real
                # remove off-diagonal terms
                if self.solver_params['diag_delta']:
                    solver_eal = np.diag(np.diag(solver_eal))
            h_loc.Fill(name, solver_eal)

        # fill solver h_loc
        self.triqs_solver.e0 = h_loc

        # FIXME: unfortunately, in the current implementation the solver initializations aren't included yet in dmft_cycle,
        # so for debugging it is done here again
        # store solver to h5 archive
        if self.general_params['store_solver'] and mpi.is_master_node():
            with HDFArchive(self.general_params['jobname'] + '/' + self.general_params['seedname'] + '.h5', 'a') as archive:
                if not 'it_-1' in archive['DMFT_input/solver']:
                    archive['DMFT_input/solver'].create_group('it_-1')
                archive['DMFT_input/solver'].create_group('it_-1')
                archive['DMFT_input/solver/it_-1']['Delta'] = self.Delta_freq_solver
                archive['DMFT_input/solver/it_-1']['S_' + str(self.icrsh)] = self.triqs_solver

        # Solve the impurity problem for icrsh shell
        # *************************************
        path_to_gs = self.solver_params['path_to_gs'] if self.solver_params['path_to_gs'] is not None and self.path_to_gs_accepted else None
        # fix to make sure this is only done in iteration 1
        if self.path_to_gs_accepted:
            self.path_to_gs_accepted = False
        if path_to_gs != 'postprocess':
            self.triqs_solver.solve(
                h_int=self.h_int,
                params_GS=self.dmrg_params,
                params_partSector=self.sector_params,
                tevo=self.tevo_params,
                eta=self.general_params['eta'],
                calc_me=self.calc_me,
                state_storage=self.solver_params['state_storage'],
                path_to_gs=path_to_gs,
            )
        else:
            self.triqs_solver.post_process(
                h_int=self.h_int,
                params_GS=self.dmrg_params,
                params_partSector=self.dmrg_params,
                tevo=self.tevo_params,
                eta=self.general_params['eta'],
                calc_me=self.calc_me,
                state_storage=self.solver_params['state_storage'],
            )
            # *************************************

        # call postprocessing
        self.postprocess()

        return

    def postprocess(self):
        r"""
        Organize G_freq, G_time, Sigma_freq and G_l from ftps solver
        """

        # symmetrization of reduced solver G
        def symmetrize_opt(G_in, soc):
            G = G_in.copy()
            if soc:

                def swap_2():
                    for i in range(2):
                        G['ud_1'][i, 2] = -G['ud_1'][i, 2]
                        G['ud_1'][2, i] = -G['ud_1'][2, i]

                swap_2()
                G['ud_0'] = 0.5 * (G['ud_0'] + G['ud_1'])
                G['ud_1'] = G['ud_0']
                for name, g in G:
                    g[1, 1] = 0.5 * (g[1, 1] + g[2, 2])
                    g[2, 2] = g[1, 1]
                swap_2()
            else:
                switch = lambda spin: 'dn' if spin == 'down' else 'up'
                for key, mapto in self.calc_mapping.items():
                    spin, block = key.split('_')
                    for deg_item in mapto:
                        map_spin, map_block = deg_item.split('_')
                        mpi.report(f'mapping {spin}-{block} to {map_spin}-{map_block}...')
                        G[switch(map_spin)].data[:, int(map_block), int(map_block)] = G[switch(spin)].data[:, int(block), int(block)]
                # particle-hole symmetry: enforce mirror/point symmetry of G(w)
                if self.solver_params['ph_symm']:
                    for block, gf in G:
                        gf.data.real = 0.5 * (gf.data[::1].real - gf.data[::-1].real)
                        gf.data.imag = 0.5 * (gf.data[::1].imag + gf.data[::-1].imag)
            return G

        def symmetrize(G):
            return symmetrize_opt(G, soc=self.sum_k.corr_shells[self.icrsh]['SO'])

        def make_positive_definite(G):
            # ensure that Delta is positive definite
            for name, gf in G:
                for orb, w in product(range(gf.target_shape[0]), gf.mesh):
                    if gf[orb, orb][w].imag > 0.0:
                        gf[orb, orb][w] = gf[orb, orb][w].real + 0.0j
            return G

        G_w = symmetrize(self.triqs_solver.G_w)
        if not self.sum_k.corr_shells[self.icrsh]['SO']:
            G_w = make_positive_definite(G_w)

        # calculate Sigma_freq via Dyson
        # do not use Dyson equation directly, as G0 might have wrong eta
        Sigma_w_symm = SigmaDyson(
            Gret=self.triqs_solver.G_ret,
            bath=self.triqs_solver.b,
            hloc=self.triqs_solver.e0,
            mesh=self.Delta_freq_solver.mesh,
            eta=self.general_params['eta'],
            symmG=symmetrize,
        )

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

