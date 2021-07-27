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
Analytic continuation of the lattice Green's function to the lattice spectral
function using maxent.

Reads G_latt(i omega) from the h5 archive and writes A_latt(omega) back.

Known problems:
    - when using parameter `energy_shift_orbitals` in soliDMFT, this is not
    accounted for in this routine. It requires a similar modification of the
    hopping matrix of SumkDFT as in dmft_cycle.py.

Author: Max Merkel, 2020
"""

import sys
import time
import glob
from multiprocessing import Pool
from functools import partial
import numpy as np

from triqs_maxent.tau_maxent import TauMaxEnt
from triqs_maxent.omega_meshes import HyperbolicOmegaMesh
from triqs_maxent.alpha_meshes import LogAlphaMesh
from triqs_dft_tools.sumk_dft import SumkDFT
from h5 import HDFArchive
from triqs.gf import GfImFreq


def _read_h5(external_path, iteration=None):
    """ Reads the block Green's function G(tau) from h5 archive. """

    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None
                                          else 'it_{}'.format(iteration))

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

        dc_energy = archive[h5_internal_path]['DC_energ']
        dc_potential = archive[h5_internal_path]['DC_pot']

        if 'chemical_potential_post' in archive[h5_internal_path]:
            chemical_potential = archive[h5_internal_path]['chemical_potential_post']
        else:
            # Old name for chemical_potential_post
            chemical_potential = archive[h5_internal_path]['chemical_potential']

    return sigma_iw, chemical_potential, dc_energy, dc_potential, block_structure


def _write_lattice_gf_to_h5(gf_latt_iw, external_path, iteration=None):
    """ Writes the G_latt(i omega) to the h5 archive. Mostly for faster testing. """
    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None
                                          else 'it_{}'.format(iteration))

    with HDFArchive(external_path, 'a') as archive:
        archive[h5_internal_path]['Glatt_iw'] = gf_latt_iw


def _read_lattice_gf_from_h5(external_path, iteration=None):
    """ Reads the G_latt(i omega) from the h5 archive. Mostly for faster testing. """
    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None
                                          else 'it_{}'.format(iteration))

    with HDFArchive(external_path, 'r') as archive:
        return archive[h5_internal_path]['Glatt_iw']


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


def _run_maxent(gf_latt_iw, sum_k, spins_degenerate, maxent_error=.03):
    """
    Runs maxent to get the spectral function from the list of block GF.
    If spins_degenerate, pairs with the same name except up<->down switched
    will only be calculated once.
    """

    results = {}

    # Prints information on the blocks found
    print('Found blocks {}'.format(list(gf_latt_iw.indices)))
    for block in gf_latt_iw.indices:
        # Checks if gf is part of a degenerate pair
        gf = _get_nondegenerate_greens_functions(spins_degenerate, block, gf_latt_iw)
        if gf is None:
            results[block] = None
            continue

        # Take trace of lattice GF. Individual entries do not make sense
        # because the KS Hamiltonian ususally has the bands sorted by energy
        gf = GfImFreq(mesh=gf.mesh, data=np.trace(gf.data, axis1=1, axis2=2))
        # Initializes and runs the maxent solver
        num_ks_orbitals = sum_k.hopping.shape[2]
        # Gets maximal entry on the diagonal of the hopping matrix
        hopping_max = np.max(np.abs(sum_k.hopping[:, :, np.arange(num_ks_orbitals), np.arange(num_ks_orbitals)]
                                    - sum_k.chemical_potential))
        omega_limit = max(20, hopping_max+2)
        print('Calling MaxEnt with omega in +- {:.3f} eV'.format(omega_limit))
        solver = TauMaxEnt()
        solver.set_G_iw(gf)
        solver.set_error(maxent_error)
        solver.omega = HyperbolicOmegaMesh(omega_min=-omega_limit, omega_max=omega_limit, n_points=int(8*omega_limit))
        solver.alpha_mesh = LogAlphaMesh(alpha_min=1e-4, alpha_max=1e2, n_points=50)
        results[block] = solver.run()

    # Assign up's solution to down result for degenerate calculations
    for key in results:
        if results[key] is None:
            results[key] = results[key.replace('down', 'up')]

    return results


def _unpack_maxent_results(results):
    """ Converts maxent result to impurity list of dict with mesh and spectral function from each analyzer """
    mesh = {key: np.array(r.omega) for key, r in results.items()}
    data_linefit = {}
    data_chi2 = {}
    for key, result in results.items():
        data_linefit[key] = result.get_A_out('LineFitAnalyzer')
        data_chi2[key] = result.get_A_out('Chi2CurvatureAnalyzer')

    data_per_impurity = {'mesh': mesh, 'Alatt_w_line_fit': data_linefit,
                         'Alatt_w_chi2_curvature': data_chi2}
    return data_per_impurity


def _write_spectral_function_to_h5(unpacked_results, external_path, iteration=None):
    """ Writes the mesh and the maxent result for each analyzer to h5 archive. """

    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None
                                          else 'it_{}'.format(iteration))

    with HDFArchive(external_path, 'a') as archive:
        archive[h5_internal_path]['Alatt_w'] = unpacked_results


def main(external_path, iteration=None, read_g_latt_iw=False):
    """
    Main function that reads the lattice Greens function from h5, analytically
    continues it and writes the result back to the h5 archive.
    Only the trace is used because the Kohn-Sham energies ("hopping") are not
    sorted by "orbital" but by energy, leading to crossovers which can confuse
    MaxEnt.

    Parameters
    ----------
    external_path: string, path of the h5 archive
    iteration: int/string, optional, iteration to read from and write to
    read_g_latt_iw: bool, optional, reads G_latt from the h5 instead of
        generating it from Sigma(i omega). Only works if this code has run before
    Returns
    -------
    list of dict, per impurity: dict containing the omega mesh
        and A_imp from two different analyzers
    """
    start_time = time.time()

    if read_g_latt_iw:
        gf_lattice_iw = _read_lattice_gf_from_h5(external_path, iteration)
    else:
        sum_k = SumkDFT(external_path, use_dft_blocks=False)
        sigma_iw, chemical_potential, dc_energy, dc_potential, block_structure = _read_h5(external_path, iteration)
        sum_k.block_structure = block_structure
        sum_k.put_Sigma(sigma_iw)
        sum_k.set_mu(chemical_potential)
        sum_k.set_dc(dc_potential, dc_energy)

        gf_lattice_iw = sum(sum_k.lattice_gf(i)*sum_k.bz_weights[i] for i in range(sum_k.n_k))
        _write_lattice_gf_to_h5(gf_lattice_iw, external_path, iteration)
        print('Generated the lattice GF. Starting maxent now.')

    maxent_results = _run_maxent(gf_lattice_iw, sum_k, True)
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
