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
Analytic continuation of the lattice Green's function to the lattice spectral
function using maxent.

Reads G_latt(i omega) from the h5 archive and writes A_latt(omega) back. See
the docstring of main() for more information.

mpi parallelized for the generation of the imaginary-frequency lattice GF over
k points.

Author: Maximilian Merkel, Materials Theory Group, ETH Zurich, 2020 - 2022
"""

import sys
import time
import numpy as np

from triqs_maxent.tau_maxent import TauMaxEnt
from triqs_maxent.omega_meshes import HyperbolicOmegaMesh
from triqs_maxent.alpha_meshes import LogAlphaMesh
from triqs_dft_tools.sumk_dft import SumkDFT
from h5 import HDFArchive
from triqs.utility import mpi
from triqs.gf import Gf, BlockGf


def _read_h5(external_path, iteration):
    """
    Reads the h5 archive to get the Matsubara self energy, the double-counting potential,
    the chemical potential and the block structure.

    Parameters
    ----------
    external_path : string
        path to h5 archive
    iteration : int
        The iteration that is being read from, None corresponds to 'last_iter'

    Returns
    -------
    sigma_iw : list
        Self energy as block Green's function for each impurity
    chemical_potential : float
        The chemical potential of the problem. Should be approximately real
    dc_potential : list
        Double counting for each impurity
    block_structure : triqs_dft_tools.BlockStructure
        Block structure mapping from the DMFT calculation

    """

    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None
                                          else f'it_{iteration}')

    with HDFArchive(external_path, 'r') as archive:
        impurity_paths = [key for key in archive[h5_internal_path].keys() if 'Sigma_freq_' in key]
        # Sorts impurity paths by their indices, not sure if necessary
        impurity_indices = [int(s[s.rfind('_')+1:]) for s in impurity_paths]
        impurity_paths = [impurity_paths[i] for i in np.argsort(impurity_indices)]

        sigma_iw = [archive[h5_internal_path][p] for p in impurity_paths]

        block_structure = archive['DMFT_input']['block_structure']
        # Fix for archives from triqs 2 when corr_to_inequiv was in SumkDFT, not in BlockStructure
        if block_structure.corr_to_inequiv is None:
            block_structure.corr_to_inequiv = archive['dft_input/corr_to_inequiv']

        dc_potential = archive[h5_internal_path]['DC_pot']

        if 'chemical_potential_post' in archive[h5_internal_path]:
            chemical_potential = archive[h5_internal_path]['chemical_potential_post']
        else:
            # Old name for chemical_potential_post
            chemical_potential = archive[h5_internal_path]['chemical_potential']

    return sigma_iw, chemical_potential, dc_potential, block_structure


def _generate_lattice_gf(sum_k, sum_spins):
    """
    Generates the lattice GF from the SumkDFT object. If sum_spins, it
    has one block "total". Otherwise, the block names are the spins.
    """
    # Initializes lattice GF to zero for each process
    spin_blocks = ['total'] if sum_spins else sum_k.spin_block_names[sum_k.SO]
    mesh = sum_k.Sigma_imp[0].mesh
    trace_gf_latt = {key: Gf(mesh=mesh, data=np.zeros((len(mesh), 1, 1), dtype=complex))
                     for key in spin_blocks}

    # Takes trace over orbitals (and spins). Individual entries do not make sense
    # because the KS Hamiltonian ususally has the bands sorted by energy
    for ik in mpi.slice_array(np.arange(sum_k.n_k)):
        gf_latt = sum_k.lattice_gf(ik) * sum_k.bz_weights[ik]
        if sum_spins:
            trace_gf_latt['total'].data[:] += np.trace(sum(g.data for _, g in gf_latt),
                                                       axis1=1, axis2=2).reshape(-1, 1, 1)
        else:
            for s in spin_blocks:
                trace_gf_latt[s].data[:] += np.trace(gf_latt[s].data,
                                                     axis1=1, axis2=2).reshape(-1, 1, 1)

    for s in spin_blocks:
        trace_gf_latt[s] << mpi.all_reduce(trace_gf_latt[s])

    # Lattice GF as BlockGf, required for compatibility with MaxEnt functions
    gf_lattice_iw = BlockGf(name_list=trace_gf_latt.keys(),
                            block_list=trace_gf_latt.values())
    return gf_lattice_iw


def _run_maxent(gf_lattice_iw, sum_k, error, omega_min, omega_max,
                n_points_maxent, n_points_alpha, analyzer='LineFitAnalyzer'):
    """
    Runs maxent to get the spectral function from the block GF.
    """

    # Automatic determination of energy range from hopping matrix
    if omega_max is None:
        num_ks_orbitals = sum_k.hopping.shape[2]
        hopping_diagonal = sum_k.hopping[:, :, np.arange(num_ks_orbitals), np.arange(num_ks_orbitals)]
        hopping_min = np.min(hopping_diagonal)
        hopping_max = np.max(hopping_diagonal)
        omega_min = min(-20, hopping_min - sum_k.chemical_potential)
        omega_max = max(20, hopping_max - sum_k.chemical_potential)
        mpi.report('Set omega range to {:.3f}...{:.3f} eV'.format(omega_min, omega_max))

    omega_mesh = HyperbolicOmegaMesh(omega_min=omega_min, omega_max=omega_max,
                                     n_points=n_points_maxent)

    # Prints information on the blocks found
    mpi.report('Found blocks {}'.format(list(gf_lattice_iw.indices)))

    # Initializes and runs the maxent solver
    # TODO: parallelization over blocks
    results = {}
    for block, gf in gf_lattice_iw:
        mpi.report('-'*80, f'Running MaxEnt on block "{block}" now', '-'*80)
        solver = TauMaxEnt()
        solver.set_G_iw(gf)
        solver.set_error(error)
        solver.omega = omega_mesh
        solver.alpha_mesh = LogAlphaMesh(alpha_min=1e-6, alpha_max=1e2, n_points=n_points_alpha)
        results[block] = solver.run()

        opt_alpha = results[block].analyzer_results[analyzer]['alpha_index']
        mpi.report(f'Optimal alpha, block "{block}" from {analyzer}: {opt_alpha}')

    return results, omega_mesh


def _unpack_maxent_results(results, omega_mesh):
    """
    Converts maxent result to dict with mesh and spectral function from each
    analyzer.
    """
    data_linefit = {}
    data_chi2 = {}
    for key, result in results.items():
        data_linefit[key] = result.get_A_out('LineFitAnalyzer')
        data_chi2[key] = result.get_A_out('Chi2CurvatureAnalyzer')

    data = {'mesh': np.array(omega_mesh), 'Alatt_w_line_fit': data_linefit,
            'Alatt_w_chi2_curvature': data_chi2}
    return data


def _write_spectral_function_to_h5(unpacked_results, external_path, iteration):
    """ Writes the mesh and the maxent result for each analyzer to h5 archive. """

    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None
                                          else f'it_{iteration}')

    with HDFArchive(external_path, 'a') as archive:
        archive[h5_internal_path]['Alatt_maxent'] = unpacked_results


def main(external_path, iteration=None, sum_spins=False, maxent_error=.02,
         n_points_maxent=200, n_points_alpha=50, omega_min=None, omega_max=None):
    """
    Main function that reads the lattice Green's function (GF) from h5,
    analytically continues it, writes the result back to the h5 archive and
    also returns the results.
    Only the trace can be used because the Kohn-Sham energies ("hopping") are not
    sorted by "orbital" but by energy, leading to crossovers.

    Parameters
    ----------
    external_path: string
        Path to the h5 archive to read from and write to.
    iteration: int/string
        Iteration to read from and write to. Defaults to last_iter.
    sum_spins: bool
        Whether to sum over the spins or continue the lattice GF
        for the up and down spin separately, for example for magnetized results.
    maxent_error : float
        The error that is used for the analyzers.
    n_points_maxent : int
        Number of omega points on the hyperbolic mesh used in the continuation.
    n_points_alpha : int
        Number of points that the MaxEnt alpha parameter is varied on logarithmically.
    omega_min : float
        Lower end of range where the GF is being continued. Range has to comprise
        all features of the lattice GF for correct normalization.
        If omega_min and omega_max are None, they are chosen automatically based
        on the diagonal entries in the hopping matrix but at least to -20...20 eV.
    omega_max : float
        Upper end of range where the GF is being continued. See omega_min.

    Returns
    -------
    unpacked_results : dict
        The omega mesh and lattice spectral function from two different analyzers
    """

    if (omega_max is None and omega_min is not None
            or omega_max is not None and omega_min is None):
        raise ValueError('Both or neither of omega_max and omega_min have to be None')

    start_time = time.time()

    # Sets up the SumkDFT object
    h5_content = None
    if mpi.is_master_node():
        h5_content = _read_h5(external_path, iteration)
    sigma_iw, chemical_potential, dc_potential, block_structure = mpi.bcast(h5_content)
    sum_k = SumkDFT(external_path, mesh=sigma_iw[0].mesh, use_dft_blocks=False)
    sum_k.block_structure = block_structure
    sum_k.put_Sigma(sigma_iw)
    sum_k.set_mu(chemical_potential)
    sum_k.set_dc(dc_potential, None)

    # Generates the lattice GF
    gf_lattice_iw = _generate_lattice_gf(sum_k, sum_spins)
    mpi.report('Generated the lattice GF.')

    # Runs MaxEnt
    unpacked_results = None
    if mpi.is_master_node():
        maxent_results, omega_mesh = _run_maxent(gf_lattice_iw, sum_k, maxent_error,
                                                 omega_min, omega_max, n_points_maxent,
                                                 n_points_alpha)
        unpacked_results = _unpack_maxent_results(maxent_results, omega_mesh)
        _write_spectral_function_to_h5(unpacked_results, external_path, iteration)
    unpacked_results = mpi.bcast(unpacked_results)

    total_time = time.time() - start_time
    mpi.report('-'*80, 'DONE')
    mpi.report(f'Total run time: {total_time:.0f} s.')

    return unpacked_results


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
