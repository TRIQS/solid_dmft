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
Updates the h5 archives with soliDMFT results for continuing old calculations.
"""
import os.path
import shutil
import sys
import numpy as np
from h5 import HDFArchive
# Import needed for transfering Green's functions in _update_results_per_iteration
from triqs.gf import BlockGf

def _backup_old_file(path):
    """ Copies file to same folder with the prefix "backup_". """
    directory, filename = os.path.split(path)
    shutil.copy2(path, os.path.join(directory, 'backup_'+filename))

def _update_results_per_iteration(archive):
    iterations = [key for key in archive['DMFT_results'] if key.startswith('it_')]
    iterations.append('last_iter')

    for it in iterations:
        keys = list(archive['DMFT_results'][it].keys())
        for key in keys:
            if '_iw' in key:
                new_key = key.replace('_iw_', '_freq_').replace('_iw', '_freq_')
                archive['DMFT_results'][it][new_key] = archive['DMFT_results'][it][key]
                del archive['DMFT_results'][it][key]
            elif '_tau' in key:
                new_key = key.replace('_tau_', '_time_').replace('_tau', '_time_')
                archive['DMFT_results'][it][new_key] = archive['DMFT_results'][it][key]
                del archive['DMFT_results'][it][key]
            # Updates chemical potential to _post
            elif key == 'chemical_potential':
                it_number = int(it[3:]) if it.startswith('it_') else -1
                archive['DMFT_results'][it]['chemical_potential_post'] = archive['DMFT_results'][it]['chemical_potential']
                archive['DMFT_results'][it]['chemical_potential_pre'] = archive['DMFT_results']['observables']['mu'][it_number]
                del archive['DMFT_results'][it]['chemical_potential']

def _update_observables(archive):
    if 'orb_Z' not in archive['DMFT_results']['observables']:
        orb_gb2 = archive['DMFT_results']['observables']['orb_gb2']
        orb_Z = []

        for i in range(len(orb_gb2)):
            orb_Z.append({})
            for spin in orb_gb2[i]:
                orb_Z[i][spin] = [np.full(np.shape(z_per_it), 'none')
                                  for z_per_it in orb_gb2[i][spin]]

        archive['DMFT_results']['observables']['orb_Z'] = orb_Z

def main(path='dmft_config.ini'):
    """ Combines methods in the full work flow for updating the h5 archive. """
    if not os.path.isfile(path):
        raise ValueError('File {} does not exist.'.format(path))

    _backup_old_file(path)

    with HDFArchive(path, 'a') as archive:
        _update_results_per_iteration(archive)
        _update_observables(archive)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    elif len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        raise TypeError('Maximally one argument supported: the archive file name')
