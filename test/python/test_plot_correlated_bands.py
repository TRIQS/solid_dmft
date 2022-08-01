# Copyright (c) 2018-2022 Simons Foundation
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You may obtain a copy of the License at
#     https:#www.gnu.org/licenses/gpl-3.0.txt
#
# Authors: Alexander Hampel

import numpy as np

from h5 import HDFArchive
from solid_dmft.postprocessing import plot_correlated_bands as pcb

import unittest


class test_convergence(unittest.TestCase):

    def setUp(self):

        self.w90_dict = {'w90_seed': 'svo', 'w90_path': './', 'mu_tb': 12.3958, 'n_orb': 3,
                         'orbital_order_w90': ['dxz', 'dyz', 'dxy'], 'add_spin': False}

        self.orbital_order_to = ['dxy', 'dxz', 'dyz']

        self.sigma_dict = {'dmft_path': './svo_example.h5', 'it': 'last_iter',
                           'orbital_order_dmft': self.orbital_order_to, 'spin': 'up',
                           'block': 0, 'eta': 0.0, 'linearize': False}

    def test_get_dmft_bands(self):
        tb_bands = {'bands_path': [('R', 'G'), ('G', 'X'), ('X', 'M'), ('M', 'G')], 'G': [0., 0., 0.],
                    'Z': np.array([0, 0, 0.5]), 'M': [0.5, 0.5, 0.], 'R': [0.5, 0.5, 0.5],
                    'X': [0.,  0.5, 0.], 'n_k': 50}

        tb_data, alatt_k_w, freq_dict = pcb.get_dmft_bands(with_sigma='calc', add_mu_tb=True,
                                                           orbital_order_to=self.orbital_order_to,
                                                           **self.w90_dict, **tb_bands, **self.sigma_dict)

        with HDFArchive('test_pcb_ref.h5', 'r') as ar:
            emat_ref = ar['tb_emat']
            Akw_ref = ar['Akw']

        assert np.allclose(tb_data['e_mat'], emat_ref)
        assert np.allclose(alatt_k_w, Akw_ref)

    def test_get_dmft_bands_proj(self):
        tb_bands = {'bands_path': [('R', 'G'), ('G', 'X'), ('X', 'M'), ('M', 'G')], 'G': [0., 0., 0.],
                    'Z': np.array([0, 0, 0.5]), 'M': [0.5, 0.5, 0.], 'R': [0.5, 0.5, 0.5],
                    'X': [0.,  0.5, 0.], 'n_k': 50}

        tb_data, alatt_k_w, freq_dict = pcb.get_dmft_bands(with_sigma='calc', add_mu_tb=True,
                                                           orbital_order_to=self.orbital_order_to,
                                                           proj_on_orb = [0,1],
                                                           **self.w90_dict, **tb_bands, **self.sigma_dict)

        with HDFArchive('test_pcb_ref.h5', 'r') as ar:
            emat_ref = ar['tb_emat_proj']
            Akw_ref = ar['Akw_proj']

        assert np.allclose(tb_data['e_mat'], emat_ref)
        assert np.allclose(alatt_k_w, Akw_ref)

    def test_get_kslice(self):

        freq_mesh_kslice = {'window': [-0.5, 0.5], 'n_w': int(1e6)}
        sigma_dict = {'dmft_path': './svo_example.h5', 'it': 'last_iter', 'w_mesh': freq_mesh_kslice,
                      'orbital_order_dmft': self.orbital_order_to, 'spin': 'up',
                      'block': 0, 'eta': 0.0, 'linearize': False}

        tb_kslice = {'bands_path': [('Y', 'G'), ('G', 'X')], 'Y': np.array([0.5, 0.0, 0]), 'G': [0., 0., 0.],
                     'Z': np.array([0, 0, 0.5]), 'M': [0.5, 0.5, 0.], 'R': [0.5, 0.5, 0.5],
                     'X': [0.,  0.5, 0.], 'n_k': 50, 'kz': 0.0}

        tb_data, alatt_k_w, freq_dict = pcb.get_dmft_bands(fermi_slice=True, with_sigma='calc', add_mu_tb=True,
                                                           orbital_order_to=self.orbital_order_to,
                                                           **self.w90_dict, **tb_kslice, **sigma_dict)

        with HDFArchive('test_pcb_ref.h5', 'r') as ar:
            emat_ref = ar['tb_emat_slice']
            Akw_ref = ar['Akw_slice']

        assert np.allclose(tb_data['e_mat'], emat_ref)
        assert np.allclose(alatt_k_w, Akw_ref)

if __name__ == '__main__':
    unittest.main()
