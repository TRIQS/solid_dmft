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
import numpy as np

from triqs_maxent.elementwise_maxent import PoormanMaxEnt
from triqs_maxent.omega_meshes import HyperbolicOmegaMesh
from triqs_maxent.alpha_meshes import LogAlphaMesh
from triqs_maxent.logtaker import VerbosityFlags
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


def _run_maxent(gf_imp_list, maxent_error, n_points_maxent, n_points_alpha,
                omega_min, omega_max, analyzer='LineFitAnalyzer'):
    """
    Runs maxent to get the spectral functions from the list of block GFs.
    """

    omega_mesh = HyperbolicOmegaMesh(omega_min=omega_min, omega_max=omega_max,
                                     n_points=n_points_maxent)

    mpi.report(f'Continuing impurities with blocks: ')
    imps_blocks = []
    for i, block_gf in enumerate(gf_imp_list):
        blocks = list(block_gf.indices)
        mpi.report('- Imp {}: {}'.format(i, blocks))
        for block in blocks:
            imps_blocks.append((i, block))
    mpi.report('-'*50)
    imps_blocks_indices = np.arange(len(imps_blocks))

    # Initialize collective results
    data_linefit = [0] * len(imps_blocks)
    data_chi2 =  [0] * len(imps_blocks)

    mpi.barrier()
    for i in mpi.slice_array(imps_blocks_indices):
        imp, block = imps_blocks[i]
        print(f"\nRank {mpi.rank}: solving impurity {imp+1} block '{block}'")

        solver = PoormanMaxEnt(use_complex=True)
        solver.set_G_tau(gf_imp_list[imp][block])
        solver.set_error(maxent_error)
        solver.omega = omega_mesh
        solver.alpha_mesh = LogAlphaMesh(alpha_min=1e-6, alpha_max=1e2,
                                         n_points=n_points_alpha)
        # silence output
        solver.maxent_diagonal.logtaker.verbose = VerbosityFlags.Quiet
        solver.maxent_offdiagonal.logtaker.verbose = VerbosityFlags.Quiet
        results = solver.run()

        n_orb = gf_imp_list[imp][block].target_shape[0]
        opt_alpha = np.zeros((n_orb, n_orb, 2), dtype=int)
        opt_alpha[:, :, :] = -1  # set results to -1 to distinguish them from 0
        for i_orb in range(n_orb):
            for j_orb in range(n_orb):
                for l_com in range(2):  # loop over complex numbers
                    if results.analyzer_results[i_orb][j_orb][l_com] == {}:
                        continue
                    opt_alpha[i_orb, j_orb,
                              l_com] = results.analyzer_results[i_orb][j_orb][l_com][analyzer]['alpha_index']

        print(
            f"Optimal alphas , Imp {imp+1} block '{block}': \n--- Real part ---\n", opt_alpha[:, :, 0])
        if np.any(opt_alpha[:, :, 1] != -1):
            print('Imp {i+1} block {block} Imag part ---\n', opt_alpha[:, :, 1])
        if np.any(opt_alpha[:, :, 0] == -1):
            print('(a -1 indicates that maxent did not run for this block due to symmetry)')

        # store unpacked data in flatted list / maxent res object not bcastable
        data_linefit[i] = results.get_A_out('LineFitAnalyzer')
        data_chi2[i] = results.get_A_out('Chi2CurvatureAnalyzer')

    # slow barrier to reduce CPU load of waiting ranks
    mpi.barrier(1000)
    # Synchronizes information between ranks
    for i in imps_blocks_indices:
        data_linefit[i] = mpi.all_reduce(data_linefit[i])
        data_chi2[i] = mpi.all_reduce(data_chi2[i])

    # final result list
    unpacked_results = [{'mesh': np.array(omega_mesh),
                         'Aimp_w_line_fit': {},
                         'Aimp_w_chi2_curvature': {}
                         } for _ in range(len(gf_imp_list))]

    for i in imps_blocks_indices:
        imp, block = imps_blocks[i]
        unpacked_results[imp]['Aimp_w_line_fit'][block] = data_linefit[i]
        unpacked_results[imp]['Aimp_w_chi2_curvature'][block] = data_chi2[i]

    return unpacked_results


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
    maxent_results : list
        The omega mesh and impurity spectral function from two different analyzers
        in a dict for each impurity
    """

    gf_imp_tau = []
    if mpi.is_master_node():
        start_time = time.time()

        gf_imp_tau = _read_h5(external_path, iteration)
        for i, gf in enumerate(gf_imp_tau):
            gf_imp_tau[i] = _sum_greens_functions(gf, sum_spins)
    gf_imp_tau = mpi.bcast(gf_imp_tau)

    maxent_results = _run_maxent(gf_imp_tau, maxent_error, n_points_maxent,
                                             n_points_alpha, omega_min, omega_max)

    if mpi.is_master_node():
        _write_spectral_function_to_h5(maxent_results, external_path, iteration)

        total_time = time.time() - start_time
        mpi.report('-'*80, 'DONE')
        mpi.report(f'Total run time: {total_time:.0f} s.')

    return maxent_results


def _strtobool(val):
    """Convert a string representation of truth to true (1) or false (0).
    True values are 'y', 'yes', 't', 'true', 'on', and '1'; false values
    are 'n', 'no', 'f', 'false', 'off', and '0'.  Raises ValueError if
    'val' is anything else.
    Copied from distutils.util in python 3.10.
    """
    val = val.lower()
    if val in ('y', 'yes', 't', 'true', 'on', '1'):
        return 1
    elif val in ('n', 'no', 'f', 'false', 'off', '0'):
        return 0
    else:
        raise ValueError("invalid truth value {!r}".format(val))


if __name__ == '__main__':
    # Casts input parameters
    if len(sys.argv) > 2:
        if sys.argv[2].lower() == 'none':
            sys.argv[2] = None
    if len(sys.argv) > 3:
        sys.argv[3] = _strtobool(sys.argv[3])
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
