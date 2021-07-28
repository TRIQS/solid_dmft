#!/usr/bin/env python3
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
Reads the -kslice-bands.dat and the -kslice-coord.dat file (as Wannier90 writes them).
The -kslice-bands.dat contains the band energies corresponding to the slices through
k-space given in _kslice-coords.dat. The latter has the list of k points in 2D direct
coordinates.

This only works for k independent projectors as from a TB model or from Wannier90.

Writes all the information back into the h5 archive in the group 'dft_bands_input',
which is needed for plotting DMFT bands with SumkDFTTools spaghettis.

Adapted from "write_bands_to_h5.py" by Sophie Beck, 2021
"""

import sys
import numpy as np
from h5 import HDFArchive


def _read_bands(seedname):
    """ Reads the -kslice-bands.dat and the -kslice-coord.dat file. """

    print('Reading {0}-kslice-bands.dat and {0}-kslice-coord.dat'.format(seedname))

    kpoints = np.loadtxt('{}-kslice-coord.dat'.format(seedname), skiprows=0, usecols=(0, 1))
    # to avoid issues with scientific notation of decimals
    kpoints = np.around(kpoints, decimals=10)
    band_energies = np.loadtxt('{}-kslice-bands.dat'.format(seedname), skiprows=0, usecols=(0))

    # reshape to band indices
    sub_bands = band_energies.size//len(kpoints)
    band_energies = band_energies.reshape((len(kpoints), sub_bands))

    # blow up to mimic using projectors
    band_energies = np.array([np.diag(e) for e in band_energies], dtype=complex)
    # add dummy spin components
    band_energies = band_energies.reshape((kpoints.shape[0], 1, band_energies.shape[1], band_energies.shape[1]))
    kpoints = np.append(kpoints, np.zeros((kpoints.shape[0], 1)),axis=1)

    return kpoints, band_energies


def _read_h5_dft_input_proj_mat(archive_name):
    """
    Reads the projection matrix from the h5. In the following,
    it is assumed to be k independent.
    """
    with HDFArchive(archive_name, 'r') as archive:
        return archive['dft_input/proj_mat']


def _write_dft_bands_input_to_h5(archive_name, data):
    """Writes all the information back to the h5 archive. data is a dict. """
    with HDFArchive(archive_name, 'a') as archive:
        if 'dft_bands_input' in archive:
            del archive['dft_bands_input']
        archive.create_group('dft_bands_input')
        for key in data:
            archive['dft_bands_input'][key] = data[key]
    print('Written results to {}'.format(archive_name))


def main(seedname, filename_archive=None):
    """
    Executes the program on the band data from the files <seedname>_bands.dat and
    <seedname>_bands.kpt. If no seedname_archive is specified, <seedname>.h5 is used.
    """

    if filename_archive is None:
        filename_archive = seedname + '.h5'
        print('Using the archive "{}"'.format(filename_archive))

    kpoints, band_energies = _read_bands(seedname)
    dft_proj_mat = _read_h5_dft_input_proj_mat(filename_archive)

    data = {'n_k': kpoints.shape[0],
            'n_orbitals': np.ones((kpoints.shape[0], 1), dtype=int) * 3, # The 1 in here only works for SO == 0
            'proj_mat': np.broadcast_to(dft_proj_mat[0],
                                        (kpoints.shape[0], ) + dft_proj_mat.shape[1:]),
            'hopping': band_energies,
            # Quantities are not used for unprojected spaghetti
            'n_parproj': 'none',
            'proj_mat_all': 'none',
            # Quantity that SumkDFTTools does not need but that is nice for plots
            'kpoints': kpoints}

    _write_dft_bands_input_to_h5(filename_archive, data)

if __name__ == '__main__':
    if len(sys.argv) == 2:
        main(sys.argv[1])
    elif len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        print('Please give a seedname (and optionally an archive to write to). Exiting.')
        sys.exit(2)
