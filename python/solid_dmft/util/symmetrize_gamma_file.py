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

import os
import numpy as np
from h5 import HDFArchive

def _read_archive(archive_name):
    with HDFArchive(archive_name, 'r') as archive:
        n_kpts = archive['dft_input']['n_k']
        kpts = archive['dft_input/kpts']
        n_bands = archive['dft_input']['n_orbitals'][0][0]
        band_window = archive['dft_misc_input']['band_window']
    assert kpts.shape[0] == n_kpts
    assert np.all(band_window[0][0] == band_window)

    return n_kpts, kpts, n_bands, band_window[0][0]

def _read_gamma(n_bands):
    with open('GAMMA', 'r') as file:
        raw_data = np.loadtxt((line for line in file if '.' in line))

    data = raw_data[:, ::2] + 1j*raw_data[:, 1::2]
    data = data.reshape((-1, n_bands, n_bands))
    return data

def _read_symmetries_outcar(n_kpts):
    """ Only written if ISYM = 2 and LWANNIER90 = .TRUE. """
    def read_outcar(file):
        has_started_reading = False
        for line in file:
            if 'IBZKPT_HF' in line:
                has_started_reading = True
                continue

            if not has_started_reading:
                continue

            if 't-inv' in line:
                yield line
                continue

            if '-'*10 in line:
                break

    with open('OUTCAR', 'r') as file:
        outcar_data_raw = np.loadtxt(read_outcar(file), usecols=[0, 1, 2, 4])
    assert outcar_data_raw.shape[0] == n_kpts
    outcar_kpoints = outcar_data_raw[:, :3]
    outcar_indices = (outcar_data_raw[:, 3]-.5).astype(int)

    return outcar_kpoints, outcar_indices

def _create_symmetry_mapping(outcar_kpoints, outcar_indices, kpts):
    symmetry_mapping = np.full(outcar_kpoints.shape[0], -1, dtype=int)

    for i, (kpt_outcar, outcar_index) in enumerate(zip(outcar_kpoints, outcar_indices)):
        for j, kpt in enumerate(kpts):
            if np.allclose(kpt_outcar, kpt):
                # Symmetry-irreducible k points
                if i == outcar_index:
                    symmetry_mapping[j] = outcar_index
                # Symmetry-reducible
                else:
                    symmetry_mapping[j] = outcar_index
                break

        # Asserts that loop left through break, i.e. a pair was found
        assert np.allclose(kpt_outcar, kpt)

    irreds, irred_indices, irred_degcy = np.unique(symmetry_mapping, return_index=True, return_counts=True)
    assert np.all(np.diff(irreds) == 1)
    assert np.all(symmetry_mapping > -1)
    n_irreds = irreds.shape[0]

    return symmetry_mapping, irred_indices, irred_degcy, n_irreds

def _check_kpoint_order(kpts, irred_indices):
    def read_ibzkpt(file):
        has_started_reading = False
        for line in file:
            if 'Reciprocal lattice' in line:
                has_started_reading = True
                continue

            if not has_started_reading:
                continue

            if 'Tetrahedra' in line:
                break

            yield line

    with open('IBZKPT', 'r') as file:
        ibzpts_data = np.loadtxt(read_ibzkpt(file))
    assert np.allclose(ibzpts_data[:, :3], kpts[irred_indices])

def _symmetrize_gamma(gamma_data, irred_indices):
    """ Output for n_irreds k points """

    delta_n_irred = gamma_data[irred_indices]
    print('Reduced GAMMA to IBZ')

    return delta_n_irred

def _write_symmetrized_gamma_irred_only(n_irreds, band_window, n_bands, delta_n_irred):
    os.rename('GAMMA', 'GAMMA_nosym')
    assert os.path.isfile('GAMMA_nosym')

    with open('GAMMA', 'w') as file:
        file.write('{}  -1  ! Number of k-points, default number of bands\n'.format(n_irreds))
        for ik in range(n_irreds):
            ib1, ib2 = band_window
            file.write(' %i  %i  %i\n'%(ik + 1, ib1, ib2))
            for inu in range(n_bands):
                for imu in range(n_bands):
                    valre = (delta_n_irred[ik][inu, imu].real + delta_n_irred[ik][inu, imu].real) / 2.0
                    valim = (delta_n_irred[ik][inu, imu].imag + delta_n_irred[ik][inu, imu].imag) / 2.0
                    file.write(' %.14f  %.14f'%(valre, valim))
                file.write('\n')

def _write_symmetrized_gamma(n_kpts, band_window, n_bands, delta_n_irred, symmetry_mapping):
    os.rename('GAMMA', 'GAMMA_nosym')
    assert os.path.isfile('GAMMA_nosym')

    with open('GAMMA', 'w') as file:
        file.write('{}  -1  ! Number of k-points, default number of bands\n'.format(n_kpts))
        for ik in range(n_kpts):
            ib1, ib2 = band_window
            file.write(" %i  %i  %i\n"%(ik + 1, ib1, ib2))
            delta_n = delta_n_irred[symmetry_mapping[ik]]
            for inu in range(n_bands):
                for imu in range(n_bands):
                    valre = delta_n[inu, imu].real
                    valim = delta_n[inu, imu].imag
                    file.write(" %.14f  %.14f"%(valre, valim))
                file.write("\n")

def run(seedname, write_irred_only=False):
    n_kpts, kpts, n_bands, band_window = _read_archive(seedname+'.h5')
    gamma_data = _read_gamma(n_bands)
    if gamma_data.shape[0] != n_kpts:
        print('WARNING: cannot symmetrize GAMMA.')
        print('Got {} kpts, expected {}'.format(gamma_data.shape[0], n_kpts))
        return
    outcar_kpoints, outcar_indices = _read_symmetries_outcar(n_kpts)

    (symmetry_mapping, irred_indices,
     irred_degcy, n_irreds) = _create_symmetry_mapping(outcar_kpoints, outcar_indices, kpts)

    if write_irred_only:
        _check_kpoint_order(kpts, irred_indices)
    delta_n_irred = _symmetrize_gamma(gamma_data, irred_indices)

    if write_irred_only:
        _write_symmetrized_gamma_irred_only(n_irreds, band_window, n_bands, delta_n_irred)
    else:
        _write_symmetrized_gamma(n_kpts, band_window, n_bands, delta_n_irred, symmetry_mapping)

if __name__ == '__main__':
    run('wannier90', True)
