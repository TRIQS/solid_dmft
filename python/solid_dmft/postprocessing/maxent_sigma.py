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
Analytic continuation of the self-energy using maxent on an auxiliary Green's
function.

Reads Sigma(i omega) from the h5 archive and writes Sigma(omega) back. See
the docstring of main() for more information.

mpi parallelized for the maxent routine over all blocks and for the continuator
extraction over omega points.

Author: Maximilian Merkel, Materials Theory Group, ETH Zurich, 2020 - 2022

Warnings:
    * When using this on self-energies with SOC, please check that the formalism
      is correct, in particular the Kramers-Kronig relation.
"""

import time
import sys
import numpy as np
import warnings

from triqs.utility import mpi
from triqs_maxent.sigma_continuator import InversionSigmaContinuator, DirectSigmaContinuator
from triqs_maxent.elementwise_maxent import PoormanMaxEnt
from triqs_maxent.omega_meshes import HyperbolicOmegaMesh
from triqs_maxent.alpha_meshes import LogAlphaMesh
from triqs_maxent.logtaker import VerbosityFlags
from h5 import HDFArchive


def _read_h5(external_path, iteration):
    """
    Reads the h5 archive to get the Matsubara self energy, the double-counting potential
    and the chemical potential.

    Parameters:
    -----------
    external_path : string
        path to h5 archive
    iteration : int
        The iteration that is being read from, None corresponds to 'last_iter'

    Returns:
    --------
    sigma_iw : list
        Self energy as block Green's function for each impurity
    dc_potential : list
        Double counting for each impurity
    chemical_potential : float
        The chemical potential of the problem. Should be approximately real
    chemical_potential_zero : float
        The chemical potential at 0 iteration. Should be approximately real
    """

    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None
                                          else f'it_{iteration}')

    with HDFArchive(external_path, 'r') as archive:
        impurity_paths = [key for key in archive[h5_internal_path].keys() if 'Sigma_freq_' in key]
        # Sorts impurity paths by their indices, not sure if necessary
        impurity_indices = [int(s[s.rfind('_')+1:]) for s in impurity_paths]
        impurity_paths = [impurity_paths[i] for i in np.argsort(impurity_indices)]
        sigma_iw = [archive[h5_internal_path][p] for p in impurity_paths]

        inequiv_to_corr = archive['dft_input']['inequiv_to_corr']
        dc_potential = [archive[h5_internal_path]['DC_pot'][icrsh]
                        for icrsh in inequiv_to_corr]

        if 'chemical_potential_post' in archive[h5_internal_path]:
            chemical_potential = archive[h5_internal_path]['chemical_potential_post']
        else:
            # Old name for chemical_potential_post
            chemical_potential = archive[h5_internal_path]['chemical_potential']
        chemical_potential_zero = archive['DMFT_results/observables']['mu'][0]

    return sigma_iw, dc_potential, chemical_potential, chemical_potential_zero


def _create_sigma_continuator(sigma_iw, dc_potential, chemical_potential, chemical_potential_zero, continuator_type):
    """
    Initializes the inversion and direct sigma continuator. Returns a list of
    continuators. Types of supported auxiliary Green's functions:
    * 'inversion_dc': inversion continuator, constant C = dc_potential for the impurity
    * 'inversion_sigmainf': inversion continuator, constant C = Sigma(i infinity) + chemical potential
    * 'direct': direct continuator
    """

    for sigma_imp in sigma_iw:
        for _, sigma_block in sigma_imp:
            if sigma_block.data.shape[1] > 1:
                warnings.warn('Continuation of matrix-valued selfenergies '
                              + 'with nonzero offdiagonal components can be '
                              + 'unstable since MaxEnt matrix continuation '
                              + 'does not guarantee a positive semi-definite, '
                              + 'Hermitian output.')

    n_inequiv_shells = len(sigma_iw)

    if continuator_type == 'inversion_dc':
        shifts = [None] * n_inequiv_shells
        for iineq in range(n_inequiv_shells):
            for dc_block in dc_potential[iineq].values():
                # Reads first element from matrix for shift
                if shifts[iineq] is None:
                    shifts[iineq] = dc_block[0, 0]
                # Checks that matrix for up and down is unit matrix * shift
                if not np.allclose(dc_block, np.eye(dc_block.shape[0])*shifts[iineq]):
                    raise NotImplementedError('Only scalar dc per impurity supported')

        continuators = [InversionSigmaContinuator(sigma_imp, shift)
                        for sigma_imp, shift in zip(sigma_iw, shifts)]
    elif continuator_type == 'inversion_sigmainf':
        shifts = [{key: sigma_block.data[-1].real + (chemical_potential - chemical_potential_zero)
                   for key, sigma_block in sigma_imp} for sigma_imp in sigma_iw]
        continuators = [InversionSigmaContinuator(sigma_imp, shift)
                        for sigma_imp, shift in zip(sigma_iw, shifts)]
    elif continuator_type == 'direct':
        for sigma_imp in sigma_iw:
            for _, sigma_block in sigma_imp:
                if sigma_block.data.shape[1] > 1:
                    # TODO: implement making input diagonal if it is not
                    raise NotImplementedError('Continuing only diagonal elements of non-diagonal '
                                              'matrix not implemented yet')
        continuators = [DirectSigmaContinuator(sigma_imp) for sigma_imp in sigma_iw]
    else:
        raise NotImplementedError

    mpi.report(f'Created sigma continuator of type "{continuator_type}"')

    return continuators


def _run_maxent(continuators, error, omega_min, omega_max, n_points_maxent,
                n_points_alpha, analyzer):
    """
    Uses maxent to continue the auxiliary Green's function obtained from the
    continuator. The range for alpha is set to 1e-6 to 1e2.
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
    opt_alphas = [np.zeros(1, dtype=int)] * len(imps_blocks)
    omega_mesh = HyperbolicOmegaMesh(omega_min=omega_min, omega_max=omega_max, n_points=n_points_maxent)

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
        # For some reason, MaxEnt still prints "appending"
        solver.maxent_diagonal.logtaker.verbose = VerbosityFlags.Quiet
        solver.maxent_offdiagonal.logtaker.verbose = VerbosityFlags.Quiet
        result = solver.run()

        spectral_funcs[i] = result.get_A_out(analyzer)

        opt_alphas[i] = np.full(g_aux_block.data.shape[1:] + (2, ), -1, dtype=int)
        for j in range(opt_alphas[i].shape[0]):
            for k in range(j+1):
                for l in range(2): # loop over complex numbers
                    if result.analyzer_results[k][j][l] == {}:
                        continue
                    opt_alphas[i][k, j, l] = result.analyzer_results[k][j][l][analyzer]['alpha_index']

    mpi.barrier(1000)
    # Synchronizes information between ranks
    for i in imps_blocks_indices:
        spectral_funcs[i] = mpi.all_reduce(spectral_funcs[i])
        opt_alphas[i] = mpi.all_reduce(opt_alphas[i])

    for i, block_index in enumerate(imps_blocks):
        mpi.report(f'Optimal alphas, block {block_index}:')
        mpi.report('--- Real part ---', opt_alphas[i][:, :, 0])
        if np.any(opt_alphas[i][:, :, 1] != -1):
            mpi.report('--- Imag part ---', opt_alphas[i][:, :, 1])

    # Sorts results into original order of impurities and blocks
    # and adds information from Hermitian conjugate of off-diagonal elements
    sorted_spectral_funcs = [{} for _ in range(len(continuators))]
    for (imp, block), val in zip(imps_blocks, spectral_funcs):
        for i in range(val.shape[0]):
            for j in range(i):
                val[i, j] = val[j, i].conj()
        if not np.allclose(val.imag, 0):
            mpi.report('The result is complex. This might be correct but comes '
                       + 'without guarantuee of formal correctness.')
        sorted_spectral_funcs[imp][block] = val

    return sorted_spectral_funcs, omega_mesh


def _get_sigma_omega_from_aux(continuators, aux_spectral_funcs, aux_omega_mesh,
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
                                          else f'it_{iteration}')

    with HDFArchive(external_path, 'a') as archive:
        for i, (g_aux_imp, sigma_imp) in enumerate(zip(g_aux_w, sigma_w)):
            archive[h5_internal_path][f'Sigma_maxent_{i}'] = sigma_imp
            archive[h5_internal_path][f'G_aux_for_Sigma_maxent_{i}'] = g_aux_imp


def main(external_path, iteration=None, continuator_type='inversion_sigmainf', maxent_error=.02,
         omega_min=-12., omega_max=12., n_points_maxent=400, n_points_alpha=50,
         analyzer='LineFitAnalyzer', n_points_interp=2000, n_points_final=1000):
    """
    Main function that reads the Matsubara self-energy from h5, analytically continues it,
    writes the results back to the h5 archive and also returns the results.

    Function parallelizes using MPI over impurities and blocks.

    Parameters
    ----------
    external_path : string
        Path to the h5 archive to read from and write to
    iteration : int/string
        Iteration to read from and write to. Default to last_iter
    continuator_type : string
        Type of continuator to use, one of 'inversion_sigmainf', 'inversion_dc', 'direct'
    maxent_error : float
        The error that is used for the analyzers.
    omega_min : float
        Lower end of range where Sigma is being continued. Range has to comprise
        all features of the self-energy because the real part of it comes from
        the Kramers-Kronig relation applied to the auxiliary spectral function.
        For example, if the real-frequency self-energy bends at omega_min or
        omega_max, there are neglegcted features and the range should be extended.
    omega_max : float
        Upper end of range where Sigma is being continued. See omega_min.
    n_points_maxent : int
        Number of omega points on the hyperbolic mesh used in analytically
        continuing the auxiliary GF
    n_points_alpha : int
        Number of points that the MaxEnt alpha parameter is varied on logarithmically
    analyzer : string
        Analyzer used int MaxEnt, one of 'LineFitAnalyzer', 'Chi2CurvatureAnalyzer',
        'ClassicAnalyzer', 'EntropyAnalyzer', 'BryanAnalyzer'
    n_points_interp : int
        Number of points where auxiliary GF is interpolated to integrate over
        it for the Kramers-Kronig relation
    n_points_final : int
        Number of omega points the complex auxiliary GF and therefore the
        continued self-energy has on a linear grid between omega_min and omega_max

    Returns
    -------
    sigma_w : list of triqs.gf.BlockGf
        Sigma(omega) per inequivalent shell
    g_aux_w : list of triqs.gf.BlockGf
        G_aux(omega) per inequivalent shell

    Raises
    ------
    NotImplementedError
        -- When a wrong continuator type or maxent analyzer is chosen
        -- For direct continuator: when the self energy contains blocks larger
        than 1x1 (no off-diagonal continuation possible)
        -- For inversion_dc continuator: when the DC is not a diagonal matrix with
        the same entry for all blocks of an impurity. Otherwise, issues like
        the global frame violating the block structure would come up.
    """
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
        sigma_iw, dc_potential, chemical_potential, chemical_potential_zero = _read_h5(external_path, iteration)
        mpi.report('Finished reading h5 archive. Found {} impurities.'.format(len(sigma_iw)))
        continuators = _create_sigma_continuator(sigma_iw, dc_potential,
                                                 chemical_potential, chemical_potential_zero, continuator_type)
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
    mpi.report(f'Extracting Σ(ω) now with {mpi.size} process(es).')
    g_aux_w, sigma_w = _get_sigma_omega_from_aux(continuators, aux_spectral_funcs,
                                                 aux_omega_mesh, omega_min, omega_max,
                                                 n_points_interp, n_points_final)
    extract_end_time = time.time()

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
    run_time_report += '{:<8} | {:10.4f}\n'.format('Extract.', extract_end_time - maxent_end_time)
    run_time_report += '{:<8} | {:10.4f}\n'.format('Writing', all_end_time - extract_end_time)
    run_time_report += '-'*length_table + '\n'
    run_time_report += '{:<8} | {:10.4f}\n'.format('Total', all_end_time - start_time)

    mpi.report(run_time_report)
    return sigma_w, g_aux_w


if __name__ == '__main__':
    # Casts input parameters
    if len(sys.argv) > 2:
        if sys.argv[2].lower() == 'none':
            sys.argv[2] = None
    if len(sys.argv) > 4:
        sys.argv[4] = float(sys.argv[4])
    if len(sys.argv) > 5:
        sys.argv[5] = float(sys.argv[5])
    if len(sys.argv) > 6:
        sys.argv[6] = float(sys.argv[6])
    if len(sys.argv) > 7:
        sys.argv[7] = int(sys.argv[7])
    if len(sys.argv) > 8:
        sys.argv[8] = int(sys.argv[8])
    if len(sys.argv) > 10:
        sys.argv[10] = int(sys.argv[10])
    if len(sys.argv) > 11:
        sys.argv[11] = int(sys.argv[11])

    main(*sys.argv[1:])
