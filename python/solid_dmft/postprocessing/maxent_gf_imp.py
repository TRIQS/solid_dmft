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

Reads G_imp(i omega) from the h5 archive and writes A_imp(omega) back. See
the docstring of main() for more information.

Not mpi parallelized.

Author: Maximilian Merkel, Materials Theory Group, ETH Zurich, 2020 - 2022
"""

import sys
import time
import distutils.util
import numpy as np

from triqs_maxent.elementwise_maxent import PoormanMaxEnt
from triqs_maxent.omega_meshes import HyperbolicOmegaMesh
from triqs_maxent.alpha_meshes import LogAlphaMesh
from h5 import HDFArchive
from triqs.utility import mpi
from triqs.gf import BlockGf


def _read_h5(external_path, iteration):
    """
    Reads the h5 archive to get the impurity Green's functions.

    Parameters
    ----------
    external_path : string
        path to h5 archive
    iteration : int
        The iteration that is being read from, None corresponds to 'last_iter'

    Returns
    -------
    gf_imp_tau : list
        Impurity Green's function as block Green's function for each impurity
    """

    """Reads the block Green's function G(tau) from h5 archive."""

    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None
                                          else f'it_{iteration}')

    with HDFArchive(external_path, 'r') as archive:
        impurity_paths = [key for key in archive[h5_internal_path].keys()
                          if 'Gimp_time_' in key and 'orig' not in key]
        # Sorts impurity paths by their indices, not sure if necessary
        impurity_indices = [int(s[s.rfind('_')+1:]) for s in impurity_paths]
        impurity_paths = [impurity_paths[i] for i in np.argsort(impurity_indices)]

        gf_imp_tau = [archive[h5_internal_path][p] for p in impurity_paths]
    return gf_imp_tau


def _sum_greens_functions(block_gf, sum_spins):
    """
    Sums over spin channels if sum_spins. It combines "up" and "down" into one
    block "total", or for SOC, simply renames the blocks ud into "total".
    """

    if not sum_spins:
        return block_gf

    for ind in block_gf.indices:
        if ind.startswith('up_'):
            assert ind.replace('up', 'down') in block_gf.indices
        elif ind.startswith('down_'):
            assert ind.replace('down', 'up') in block_gf.indices
        elif not ind.startswith('ud_'):
            raise ValueError(f'Block {ind} in G(tau) has unknown spin type. '
                             + 'Check G(tau) or turn off sum_spins.')

    summed_gf_imp = {}

    for block_name, block in sorted(block_gf):
        if block_name.startswith('up_'):
            new_block_name = block_name.replace('up', 'total')
            opp_spin_block_name = block_name.replace('up', 'down')
            summed_gf_imp[new_block_name] = block + block_gf[opp_spin_block_name]
        elif block_name.startswith('ud_'):
            summed_gf_imp[block_name.replace('ud', 'total')] = block

    return BlockGf(name_list=summed_gf_imp.keys(), block_list=summed_gf_imp.values())


def _run_maxent(gf_imp_tau, maxent_error, n_points_maxent, n_points_alpha,
                omega_min, omega_max):
    """
    Runs maxent to get the spectral functions from the list of block GFs.
    """

    omega_mesh = HyperbolicOmegaMesh(omega_min=omega_min, omega_max=omega_max,
                                     n_points=n_points_maxent)

    if not mpi.is_master_node():
        return None, omega_mesh

    # Initializes and runs the maxent solver
    # TODO: parallelization over blocks
    results = [{} for _ in range(len(gf_imp_tau))]
    for i, block_gf in enumerate(gf_imp_tau):
        print('-'*50 + '\nSolving impurity {}/{}\n'.format(i+1, len(gf_imp_tau)) + '-'*50)
        print('Found blocks {}'.format(list(block_gf.indices)))
        for block, gf in block_gf:
            solver = PoormanMaxEnt(use_complex=True)
            solver.set_G_tau(gf)
            solver.set_error(maxent_error)
            solver.omega = omega_mesh
            solver.alpha_mesh = LogAlphaMesh(alpha_min=1e-6, alpha_max=1e2,
                                             n_points=n_points_alpha)
            results[i][block] = solver.run()

    return results, omega_mesh


def _unpack_maxent_results(results, omega_mesh):
    """
    Converts maxent result to impurity list of dict with mesh
    and spectral function from each analyzer.
    """

    data_linefit = [{key: r.get_A_out('LineFitAnalyzer') for key, r in block_res.items()}
                    for block_res in results]
    data_chi2 = [{key: r.get_A_out('Chi2CurvatureAnalyzer') for key, r in block_res.items()}
                 for block_res in results]

    data_per_impurity = [{'mesh': np.array(omega_mesh), 'Aimp_w_line_fit': dl,
                          'Aimp_w_chi2_curvature': dc}
                         for dl, dc in zip(data_linefit, data_chi2)]
    return data_per_impurity


def _write_spectral_function_to_h5(unpacked_results, external_path, iteration):
    """ Writes the mesh and the maxent result for each analyzer to h5 archive. """

    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None
                                          else f'it_{iteration}')

    with HDFArchive(external_path, 'a') as archive:
        for i, res in enumerate(unpacked_results):
            archive[h5_internal_path][f'Aimp_maxent_{i}'] = res


def main(external_path, iteration=None, sum_spins=False, maxent_error=0.02,
         n_points_maxent=200, n_points_alpha=50, omega_min=-20, omega_max=20):
    """
    Main function that reads the impurity Greens (GF) function from h5,
    analytically continues it, writes the result back to the h5 archive and
    also returns the results.

    Parameters
    ----------
    external_path : string
        Path to the h5 archive to read from and write to.
    iteration : int/string
        Iteration to read from and write to. Defaults to last_iter.
    sum_spins : bool
        Whether to sum over the spins or continue the impurity GF
        for the up and down spin separately, for example for magnetized results.
    maxent_error : float
        The error that is used for the analyzers.
    n_points_maxent : int
        Number of omega points on the hyperbolic mesh used in the continuation.
    n_points_alpha : int
        Number of points that the MaxEnt alpha parameter is varied on logarithmically.
    omega_min : float
        Lower end of range where the GF is being continued. Range has to comprise
        all features of the impurity GF for correct normalization.
    omega_max : float
        Upper end of range where the GF is being continued. See omega_min.

    Returns
    -------
    unpacked_results : list
        The omega mesh and impurity spectral function from two different analyzers
        in a dict for each impurity
    """

    start_time = time.time()

    gf_imp_tau = None
    if mpi.is_master_node():
        gf_imp_tau = _read_h5(external_path, iteration)
        for i, gf in enumerate(gf_imp_tau):
            gf_imp_tau[i] = _sum_greens_functions(gf, sum_spins)
    gf_imp_tau = mpi.bcast(gf_imp_tau)

    maxent_results, omega_mesh = _run_maxent(gf_imp_tau, maxent_error, n_points_maxent,
                                             n_points_alpha, omega_min, omega_max)

    unpacked_results = None
    if mpi.is_master_node():
        unpacked_results = _unpack_maxent_results(maxent_results, omega_mesh)
        _write_spectral_function_to_h5(unpacked_results, external_path, iteration)
    unpacked_results = mpi.bcast(unpacked_results)

    total_time = time.time() - start_time
    mpi.report('-'*80, 'DONE')
    mpi.report(f'Total run time: {total_time:.0f} s.')

    return unpacked_results


if __name__ == '__main__':
    # Casts input parameters
    if len(sys.argv) > 2:
        if sys.argv[2].lower() == 'none':
            sys.argv[2] = None
    if len(sys.argv) > 3:
        sys.argv[3] = distutils.util.strtobool(sys.argv[3])
    if len(sys.argv) > 4:
        sys.argv[4] = float(sys.argv[4])
    if len(sys.argv) > 5:
        sys.argv[5] = int(sys.argv[5])
    if len(sys.argv) > 6:
        sys.argv[6] = int(sys.argv[6])
    if len(sys.argv) > 7:
        sys.argv[7] = float(sys.argv[7])
    if len(sys.argv) > 8:
        sys.argv[8] = float(sys.argv[8])

    main(*sys.argv[1:])
