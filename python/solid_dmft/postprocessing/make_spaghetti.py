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
Script to combine the information from maxent_sigma.py and the dft_bands_input
to generate the DMFT spectral function and DMFT band structure from Sigma(omega).

Known problems:
    - when using parameter `energy_shift_orbitals` in soliDMFT, this is not
    accounted for in this routine. It requires a similar modification of the
    hopping matrix of SumkDFT as in dmft_cycle.py.

Author: Max Merkel, 2020-2021
"""

import sys
import pickle
import numpy as np

from triqs.gf import GfReFreq, BlockGf
from h5 import HDFArchive
from triqs_dft_tools.sumk_dft_tools import SumkDFTTools


def _read_h5(external_path, iteration=None):
    """Reads Sigma(omega), the DC, mu and the kpoints from the h5 """

    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None else 'it_{}'.format(iteration))

    with HDFArchive(external_path, 'r') as archive:
        impurity_paths = [key for key in archive[h5_internal_path].keys() if 'Sigma_maxent_' in key]
        # Sorts impurity paths by their indices, not sure if necessary
        impurity_indices = [int(s[s.rfind('_')+1:]) for s in impurity_paths]
        impurity_paths = [impurity_paths[i] for i in np.argsort(impurity_indices)]

        sigma_w = [archive[h5_internal_path][p] for p in impurity_paths]

        dc_potential = archive[h5_internal_path]['DC_pot']
        dc_energy = archive[h5_internal_path]['DC_energ']

        block_structure = archive['DMFT_input']['block_structure']
        # Fix for archives from triqs 2 when corr_to_inequiv was in SumkDFT, not in BlockStructure
        if block_structure.corr_to_inequiv is None:
            block_structure.corr_to_inequiv = archive['dft_input/corr_to_inequiv']

        if 'chemical_potential_post' in archive[h5_internal_path]:
            chemical_potential = archive[h5_internal_path]['chemical_potential_post']
        else:
            # Old name for chemical_potential_post
            chemical_potential = archive[h5_internal_path]['chemical_potential']

    return sigma_w, dc_potential, dc_energy, chemical_potential, block_structure


def _initialize_sum_k_tools(external_path, chemical_potential, dc_potential, dc_energy, sigma_w, block_structure):
    """ Creates the SumKDFTTools objects that will calculate the spectral properties. """
    sum_k = SumkDFTTools(hdf_file=external_path, use_dft_blocks=False)
    sum_k.block_structure = block_structure
    sum_k.set_mu(chemical_potential)
    sum_k.set_dc(dc_potential, dc_energy)
    sum_k.put_Sigma(sigma_w)
    return sum_k


def _write_dos_and_spaghetti_to_h5(mesh, dos, dos_proj, dos_proj_orb, spaghetti,
                                   external_path, iteration=None):
    """ Writes different spectral functions and the spaghetti to the h5 archive. """

    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None else 'it_{}'.format(iteration))

    results = {'Alatt_w_from_Sigma_w': dos, 'Alatt_k_w_from_Sigma_w': spaghetti}
    for i, (res_total, res_per_orb) in enumerate(zip(dos_proj, dos_proj_orb)):
        results['Aimp_w_{}_from_Sigma_w'.format(i)] = {'total': res_total,
                                                       'per_orb': res_per_orb,
                                                       'mesh': mesh}

    results['Alatt_w_from_Sigma_w']['mesh'] = mesh
    results['Alatt_k_w_from_Sigma_w']['mesh'] = mesh

    # Writes results as pickle file because hdf tends to fail sometimes
    with open(external_path+'.pkl', 'wb') as file:
        pickle.dump(results, file)

    with HDFArchive(external_path, 'a') as archive:
        dos['mesh'] = mesh
        archive[h5_internal_path]['Alatt_w_from_Sigma_w'] = dos
        for i, (res_total, res_per_orb) in enumerate(zip(dos_proj, dos_proj_orb)):
            archive[h5_internal_path]['Aimp_w_{}_from_Sigma_w'.format(i)] = {'total': res_total,
                                                                             'per_orb': res_per_orb,
                                                                             'mesh': mesh}
        spaghetti['mesh'] = mesh
        archive[h5_internal_path]['Alatt_k_w_from_Sigma_w'] = spaghetti


def main(external_path, iteration=None):
    """
    Theoretically, spectral functions not needed but good as comparison.

    Parameters
    ----------
    external_path: string, path of the h5 archive
    iteration: int/string, optional, iteration to read from and write to

    Returns
    -------
    numpy array, omega mesh for all spectral functions
    dict with 'up' and 'down' with the lattice spectral function
    list of dict, per impurity: A_imp(omega)
    list of dict, per impurity: A_imp(omega), orbital resolved
    dict with 'up' and 'down' with the k resolved lattice spectral function
    """

    sigma_w, dc_potential, dc_energy, chemical_potential, block_structure = _read_h5(external_path, iteration)
    sum_k = _initialize_sum_k_tools(external_path, chemical_potential,
                                    dc_potential, dc_energy, sigma_w, block_structure)
    mesh = np.array([x.real for x in sigma_w[0].mesh])

    alatt_w, aimp_w, aimp_w_per_orb = sum_k.dos_wannier_basis(save_to_file=False, broadening=0)
    print('Calculated the spectral functions. Starting with spaghetti now.')

    alatt_k_w = sum_k.spaghettis(save_to_file=False, broadening=0)
    _write_dos_and_spaghetti_to_h5(mesh, alatt_w, aimp_w, aimp_w_per_orb, alatt_k_w,
                                   external_path, iteration)

    print('Calculated spaghetti and wrote all results to file.')
    return mesh, alatt_w, aimp_w, aimp_w_per_orb, alatt_k_w


if __name__ == '__main__':
    if len(sys.argv) == 2:
        main(sys.argv[1])
    elif len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        print('Please give the h5 name (and optionally the iteration). Exiting.')
        sys.exit(2)
