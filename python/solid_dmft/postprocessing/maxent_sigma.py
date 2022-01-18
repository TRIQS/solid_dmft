#!/usr/bin/env python3
# -*- coding: utf-8 -*-
################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2016-2018, N. Wentzell
# Copyright (C) 2018-2019, Simons Foundation
#   author: N. Wentzell
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
"""
Analytic continuation of the self-energy using maxent.

Reads Sigma(i omega) from the h5 archive and writes Sigma(omega) back.

Analytic continuation of selfenergies requires using an auxiliary Green's
function. This makes the procedure even less exact than it already is so use
the results with care!

Author: Maximilian Merkel, Materials Theory Group, ETH Zurich, 2020 - 2022
"""

import time
import sys
import numpy as np

from triqs.utility import mpi
#from triqs.gf import BlockGf
from triqs_maxent.sigma_continuator import InversionSigmaContinuator, DirectSigmaContinuator
from triqs_maxent.elementwise_maxent import PoormanMaxEnt
from triqs_maxent.omega_meshes import HyperbolicOmegaMesh
from triqs_maxent.alpha_meshes import LogAlphaMesh
from triqs_maxent.logtaker import VerbosityFlags
from h5 import HDFArchive


def _read_h5(external_path, iteration=None):
    """
    Reads the h5 archive to get the Matsubara self energy, the double counting potential
    and the chemical potential.

    Parameters:
    -----------
    external_path: string, path to h5 archive
    iteration: int, the iteration that is being read from, default is 'last_iter'

    Returns:
    --------
    sigma_iw: list, self energy as block Green's function for each impurity
    dc_potential: list, double counting for each impurity
    chemical_potential: complex, the chemical potential of the problem. Should be approximately real
    """

    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None
                                          else 'it_{}'.format(iteration))

    with HDFArchive(external_path, 'r') as archive:
        impurity_paths = [key for key in archive[h5_internal_path].keys() if 'Sigma_freq_' in key]
        # Sorts impurity paths by their indices, not sure if necessary
        impurity_indices = [int(s[s.rfind('_')+1:]) for s in impurity_paths]
        impurity_paths = [impurity_paths[i] for i in np.argsort(impurity_indices)]
        sigma_iw = [archive[h5_internal_path][p] for p in impurity_paths]

        inequiv_to_corr = archive['dft_input']['inequiv_to_corr']
        dc_potential = [archive[h5_internal_path]['DC_pot'][icrsh] for icrsh in inequiv_to_corr]
        if 'chemical_potential_post' in archive[h5_internal_path]:
            chemical_potential = archive[h5_internal_path]['chemical_potential_post']
        else:
            # Old name for chemical_potential_post
            chemical_potential = archive[h5_internal_path]['chemical_potential']

    return sigma_iw, dc_potential, chemical_potential


def _create_sigma_continuator(sigma_iw, dc_potential, chemical_potential, continuator_type):
    """
    Initializes the inversion and direct sigma continuator from Sigma(i omega) and the DC energy.
    Returns a dict of three versions of continuators, using different auxiliary Green's functions:
    * 'inversion_continuator_dc': inversion continuator, constant C = dc_potential for the impurity
    * 'inversion_continuator_sigmainf': inversion continuator, constant C = Sigma(i infinity) + chemical potential
    * 'direct_continuator': direct continuator
    """

    if continuator_type == 'inversion_dc':
        # TODO: do with sumk to take care of rotations etc
        raise NotImplementedError('Inversion continuator with DC not implemented yet.')
        continuators = [InversionSigmaContinuator(sigma_imp, dc_imp)
                        for sigma_imp, dc_imp in zip(sigma_iw, dc_potential)]
    elif continuator_type == 'inversion_sigmainf':
        shifts = [{key: sigma_block.data[-1].real + chemical_potential
                   for key, sigma_block in sigma_imp} for sigma_imp in sigma_iw]
        continuators = [InversionSigmaContinuator(sigma_imp, shift)
                        for sigma_imp, shift in zip(sigma_iw, shifts)]
    elif continuator_type == 'direct':
        # TODO: implement making input diagonal if it is not
        for sigma_imp in sigma_iw:
            for block_name, sigma_block in sigma_imp:
                if sigma_block.data.shape[1] > 1:
                    raise NotImplementedError('Continuing only diagonal elements of non-diagonal '
                                              'matrix not implemented yet')
        # TODO: test if Linalg error happens all the time
        continuators = [DirectSigmaContinuator(sigma_imp) for sigma_imp in sigma_iw]
    else:
        raise NotImplementedError

    mpi.report(f'Created sigma continuator of type "{continuator_type}"')

    return continuators


def _run_maxent(continuators, error, omega_min, omega_max, n_points_maxent,
                n_points_alpha, analyzer):
    """
    Uses maxent to continue the auxiliary Green's function obtained from the continuator.
    The default range for omega is -12 to +12, for alpha 1e-6 to 1e2. This can be changed in here.
    Returns the real-frequency auxiliary Green's function
    """

    # Finds blocks of impurities and prints summary
    mpi.report('Continuing impurities with blocks:')
    imps_blocks = []
    for i, continuator in enumerate(continuators):
        blocks = list(continuator.Gaux_iw.indices)
        mpi.report('- Imp {}: {}'.format(i, blocks))
        for block in blocks:
            imps_blocks.append((i, block))

    # Initializes arrays to save results in
    spectral_funcs = [np.zeros(1)] * len(imps_blocks)
    opt_alphas = [np.zeros(1)] * len(imps_blocks)
    omega_mesh = HyperbolicOmegaMesh(omega_min=-12, omega_max=12, n_points=n_points_maxent)

    # Runs MaxEnt while parallelizing over impurities and blocks
    imps_blocks_indices = np.arange(len(imps_blocks))
    for i in mpi.slice_array(imps_blocks_indices):
        imp, block = imps_blocks[i]
        g_aux_block = continuators[imp].Gaux_iw[block]
        solver = PoormanMaxEnt(use_complex=True)
        solver.set_G_iw(g_aux_block)
        solver.set_error(error)
        solver.omega = omega_mesh
        solver.alpha_mesh = LogAlphaMesh(alpha_min=1e-6, alpha_max=1e2, n_points=n_points_alpha)
        # Turns off MaxEnt output, it's far too messy in the parallel mode
        # For some reason, MaxEnt still prints "appending" a lot
        solver.maxent_diagonal.logtaker.verbose = VerbosityFlags.Quiet
        solver.maxent_offdiagonal.logtaker.verbose = VerbosityFlags.Quiet
        result = solver.run()

        spectral_funcs[i] = result.get_A_out(analyzer)

        opt_alphas[i] = np.full(g_aux_block.data.shape[1:], -1, dtype=int)
        for j in range(opt_alphas[i].shape[0]):
            for k in range(j+1):
                assert result.analyzer_results[k][j][1] == {}, 'Result should not be complex'
                if analyzer not in result.analyzer_results[k][j][0]:
                    print(k, j)
                    print(result.analyzer_results[k][j][0], flush=True)
                else:
                    opt_alphas[i][k, j] = result.analyzer_results[k][j][0][analyzer]['alpha_index']

    # Synchronizes information between branches
    for i in imps_blocks_indices:
        spectral_funcs[i] = mpi.all_reduce(mpi.world, spectral_funcs[i], lambda x, y: x+y)
        opt_alphas[i] = mpi.all_reduce(mpi.world, opt_alphas[i], lambda x, y: x+y)

    mpi.report('Optimal alphas: {}'.format(opt_alphas))

    # Sorts results into original order of impurities and blocks
    # and adds information from Hermitian conjugate of off-diagonal elements
    sorted_spectral_funcs = [{} for _ in range(len(continuators))]
    for (imp, block), val in zip(imps_blocks, spectral_funcs):
        for i in range(val.shape[0]):
            for j in range(i):
                val[i, j] = val[j, i]
        assert np.allclose(val.imag, 0), 'Result should not be complex'
        sorted_spectral_funcs[imp][block] = val

    return sorted_spectral_funcs, omega_mesh


def _get_sigma_omega_from_aux_spectral(continuators, aux_spectral_funcs, aux_omega_mesh,
                                       omega_min, omega_max, n_points_interp, n_points_final):
    """ Extracts the real-frequency self energy from the auxiliary Green's function. """
    for cont_imp, spec_imp in zip(continuators, aux_spectral_funcs):
        cont_imp.set_Gaux_w_from_Aaux_w(spec_imp, aux_omega_mesh, np_interp_A=n_points_interp,
                                        np_omega=n_points_final, w_min=omega_min, w_max=omega_max)

    g_aux_w = [continuator.Gaux_w for continuator in continuators]
    sigma_w = [continuator.S_w for continuator in continuators]
    return g_aux_w, sigma_w


def _write_sigma_omega_to_h5(g_aux_w, sigma_w, external_path, iteration):
    """ Writes real-frequency self energy to h5 archive. """
    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None
                                          else 'it_{}'.format(iteration))

    with HDFArchive(external_path, 'a') as archive:
        for i, (g_aux_imp, sigma_imp) in enumerate(zip(g_aux_w, sigma_w)):
            archive[h5_internal_path][f'Sigma_w_{i}'] = sigma_imp
            archive[h5_internal_path][f'G_aux_for_Sigma_w_{i}'] = g_aux_imp


def main(external_path, iteration=None, continuator_type='inversion_sigmainf', maxent_error=.02,
         omega_min=-12, omega_max=12, n_points_maxent=400, n_points_alpha=50, analyzer='LineFitAnalyzer',
         n_points_interp=2000, n_points_final=1000):
    """
    Main function that reads the Matsubara self-energy from h5, analytically continues it
    and writes the result back to the h5 archive.

    Parameters
    ----------
    external_path: string, path of the h5 archive
    iteration: int/string, optional, iteration to read from and write to

    Returns
    -------
    list of triqs.gf.BlockGf, per impurity: Sigma(omega) as a GF object
    """
    # TODO: comment on omega_min and omega_max in maxent and KK relation

    # Checks on input parameters
    if continuator_type not in ('inversion_sigmainf', 'inversion_dc', 'direct'):
        raise NotImplementedError('Unsupported type of continuator chosen')

    if analyzer not in ('LineFitAnalyzer', 'Chi2CurvatureAnalyzer', 'ClassicAnalyzer',
                        'EntropyAnalyzer', 'BryanAnalyzer'):
        raise NotImplementedError('Unsupported type of analyzer chosen')

    assert omega_min < omega_max

    # Reads in data and initializes continuator object
    start_time = time.time()
    continuators = None
    if mpi.is_master_node():
        sigma_iw, dc_potential, chemical_potential = _read_h5(external_path, iteration)
        mpi.report('Finished reading h5 archive. Found {} impurities.'.format(len(sigma_iw)))
        continuators = _create_sigma_continuator(sigma_iw, dc_potential,
                                                 chemical_potential, continuator_type)
    continuators = mpi.bcast(continuators)
    init_end_time = time.time()

    # Runs MaxEnt
    mpi.report('Starting run of maxent now.')
    aux_spectral_funcs, aux_omega_mesh = _run_maxent(continuators, maxent_error,
                                                     omega_min, omega_max,
                                                     n_points_maxent, n_points_alpha,
                                                     analyzer)
    maxent_end_time = time.time()

    # Extracts Sigma(omega)
    mpi.report(f'Extracting Σ(ω) now with {mpi.size} process(es). This might take a while.')
    g_aux_w, sigma_w = _get_sigma_omega_from_aux_spectral(continuators, aux_spectral_funcs, aux_omega_mesh,
                                                          omega_min, omega_max, n_points_interp, n_points_final)
    extraction_end_time = time.time()

    # Writes results into h5 archive
    mpi.report('Writing results to h5 archive now.')
    if mpi.is_master_node():
        _write_sigma_omega_to_h5(g_aux_w, sigma_w, external_path, iteration)
    mpi.report('Finished writing Σ(ω) to archive.')

    all_end_time = time.time()

    # Prints timing summary
    run_time_report = '\n{:<8} | {:<10}\n'.format('Task', 'Duration (s)')
    length_table = len(run_time_report) - 2
    run_time_report += '-'*length_table + '\n'
    run_time_report += '{:<8} | {:10.4f}\n'.format('Reading', init_end_time - start_time)
    run_time_report += '{:<8} | {:10.4f}\n'.format('MaxEnt', maxent_end_time - init_end_time)
    run_time_report += '{:<8} | {:10.4f}\n'.format('Extract.', extraction_end_time - maxent_end_time)
    run_time_report += '{:<8} | {:10.4f}\n'.format('Writing', all_end_time - extraction_end_time)
    run_time_report += '-'*length_table + '\n'
    run_time_report += '{:<8} | {:10.4f}\n'.format('Total', all_end_time - start_time)

    mpi.report(run_time_report)
    return sigma_w, g_aux_w


if __name__ == '__main__':
    # TODO: add help
    if '-h' in sys.argv[1:] or '--help' in sys.argv[1:]:
        pass

    # TODO: Adapt to number of valid keywords, change documentation
    if len(sys.argv) not in (2, 3):
        mpi.report('Please give the h5 name (and optionally the iteration). Exiting.')
        sys.exit(2)

    main(*sys.argv[1:], continuator_type='inversion_sigmainf', n_points_maxent=100, n_points_alpha=20,
         n_points_interp=500, n_points_final=100, analyzer='LineFitAnalyzer')
