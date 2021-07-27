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
Analytic continuation of the impurity Green's function to the impurity spectral
function using maxent.

Reads G_imp(i omega) from the h5 archive and writes A_imp(omega) back.

Author: Max Merkel, 2020
"""

import sys
import time
import glob
from multiprocessing import Pool
from functools import partial
import numpy as np

from triqs_maxent.elementwise_maxent import PoormanMaxEnt
from triqs_maxent.omega_meshes import HyperbolicOmegaMesh
from triqs_maxent.alpha_meshes import LogAlphaMesh
from h5 import HDFArchive


def _read_h5(external_path, iteration=None):
    """Reads the block Green's function G(tau) from h5 archive."""

    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None
                                          else 'it_{}'.format(iteration))

    with HDFArchive(external_path, 'r') as archive:
        impurity_paths = [key for key in archive[h5_internal_path].keys()
                          if 'Gimp_time_' in key and 'orig' not in key]
        # Sorts impurity paths by their indices, not sure if necessary
        impurity_indices = [int(s[s.rfind('_')+1:]) for s in impurity_paths]
        impurity_paths = [impurity_paths[i] for i in np.argsort(impurity_indices)]

        gf_imp_tau = [archive[h5_internal_path][p] for p in impurity_paths]
    return gf_imp_tau


def _get_nondegenerate_greens_functions(spins_degenerate, block, block_gf):
    """
    If spins_degenerate and there is a corresponding up, down pair (same names),
    it returns the averaged from up and down for the up index and None for the down
    index. This way, the gf will only be continued analytically once for each spin pair.
    """

    if spins_degenerate:
        if 'up' in block:
            degenerate_block = block.replace('up', 'down')
            if degenerate_block in block_gf.indices:
                print(' '*10 + 'Block {}: '.format(block)
                      + 'using average with degenerate block {}.'.format(degenerate_block))
                return (block_gf[block] + block_gf[degenerate_block]) / 2
        elif 'down' in block and block.replace('down', 'up') in block_gf.indices:
            print(' '*10 + 'Block {}: skipping, same as degenerate up state.'.format(block))
            return None

    return block_gf[block]


def _run_maxent(gf_imp_tau, spins_degenerate, maxent_error=.03):
    """
    Runs maxent to get the spectral function from the list of block GF.
    If spins_degenerate, pairs with the same name except up<->down switched
    will only be calculated once.
    """

    results = [{} for _ in range(len(gf_imp_tau))]
    for i, block_gf in enumerate(gf_imp_tau):
        # Prints information on the impurity and blocks found
        print('-'*50 + '\nSolving impurity {}/{}\n'.format(i+1, len(gf_imp_tau)) + '-'*50)
        print('Found blocks {}'.format(list(block_gf.indices)))
        for block in block_gf.indices:
            # Checks if gf is part of a degenerate pair
            gf = _get_nondegenerate_greens_functions(spins_degenerate, block, block_gf)
            if gf is None:
                results[i][block] = None
                continue

            # Initializes and runs the maxent solver
            solver = PoormanMaxEnt(use_complex=True)
            solver.set_G_tau(gf)
            solver.omega = HyperbolicOmegaMesh(omega_min=-20, omega_max=20, n_points=160)
            solver.alpha_mesh = LogAlphaMesh(alpha_min=1e-4, alpha_max=1e2, n_points=50)
            solver.set_error(maxent_error)
            results[i][block] = solver.run()

        # Assign up's solution to down result for degenerate calculations
        for key in results[i]:
            if results[i][key] is None:
                results[i][key] = results[i][key.replace('down', 'up')]

    return results


def _unpack_maxent_results(results):
    """
    Converts maxent result to impurity list of dict with mesh
    and spectral function from each analyzer.
    """

    mesh = [{key: np.array(r.omega) for key, r in block_res.items()} for block_res in results]
    data_linefit = [{key: r.get_A_out('LineFitAnalyzer') for key, r in block_res.items()}
                    for block_res in results]
    data_chi2 = [{key: r.get_A_out('Chi2CurvatureAnalyzer') for key, r in block_res.items()}
                 for block_res in results]

    data_per_impurity = [{'mesh': m, 'Aimp_w_line_fit': dl, 'Aimp_w_chi2_curvature': dc}
                         for m, dl, dc in zip(mesh, data_linefit, data_chi2)]
    return data_per_impurity


def _write_spectral_function_to_h5(unpacked_results, external_path, iteration=None):
    """ Writes the mesh and the maxent result for each analyzer to h5 archive. """

    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None
                                          else 'it_{}'.format(iteration))

    with HDFArchive(external_path, 'a') as archive:
        for i, res in enumerate(unpacked_results):
            impurity_path = 'Aimp_w_{}'.format(i)
            archive[h5_internal_path][impurity_path] = res


def main(external_path, iteration=None):
    """
    Main function that reads the impurity Greens function from h5, analytically continues it
    and writes the result back to the h5 archive.

    Parameters
    ----------
    external_path: string, path of the h5 archive
    iteration: int/string, optional, iteration to read from and write to

    Returns
    -------
    list of dict, per impurity: dict containing the omega mesh
        and A_imp from two different analyzers
    """
    start_time = time.time()

    gf_imp_tau = _read_h5(external_path, iteration)
    maxent_results = _run_maxent(gf_imp_tau, True)
    unpacked_results = _unpack_maxent_results(maxent_results)
    _write_spectral_function_to_h5(unpacked_results, external_path, iteration)

    total_time = time.time() - start_time
    print('-'*50 + '\nDONE')
    print('The program took {:.0f} s.'.format(total_time))
    return unpacked_results


if __name__ == '__main__':
    if len(sys.argv) not in (2, 3):
        print('Please give the h5 name (and optionally the iteration). Exiting.')
        sys.exit(2)

    files = glob.glob(sys.argv[1])
    pool = Pool(processes=min(8, len(files)))

    if len(sys.argv) == 2:
        function = main
    elif len(sys.argv) == 3:
        function = partial(main, iteration=sys.argv[2])

    pool.map(function, files)
