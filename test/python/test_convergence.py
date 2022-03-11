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
from triqs.gf import MeshImFreq, MeshImTime, iOmega_n, inverse, BlockGf, Gf, SemiCircular, make_gf_from_fourier
from triqs.utility.comparison_tests import assert_arrays_are_close
from solid_dmft.dmft_tools.convergence import max_G_diff

import unittest


class test_convergence(unittest.TestCase):

    def setUp(self):

        beta = 10.0
        n_iw = 50
        n_tau = 6*n_iw + 1
        self.iw_mesh = MeshImFreq(beta=beta, S="Fermion", n_max=n_iw)
        self.tau_mesh = MeshImTime(beta=beta, S="Fermion", n_max=n_tau)
        self.norb = 2

        self.G1_iw = Gf(mesh=self.iw_mesh, target_shape=(self.norb, self.norb))
        self.G2_iw = Gf(mesh=self.iw_mesh, target_shape=(self.norb, self.norb))
        # self.G2_iw = self.G1_iw.copy()

        self.G1_iw << inverse(iOmega_n - SemiCircular(2.0) * np.eye(self.norb))
        self.G2_iw << inverse(iOmega_n - SemiCircular(2.0) * np.eye(self.norb)
                              +0.0001j*np.diag(np.diag(np.random.normal(-1, 1, size=(self.norb, self.norb)))))

        self.G1_tau = make_gf_from_fourier(self.G1_iw)
        self.G2_tau = make_gf_from_fourier(self.G2_iw)



    def test_Gdiff_iw(self):

        B1 = BlockGf(name_list=['0', '1'], block_list=[self.G1_iw, self.G1_iw], make_copies=True)
        B2 = BlockGf(name_list=['0', '1'], block_list=[self.G2_iw, self.G2_iw], make_copies=True)

        diff = max_G_diff(B1, B2)
        print('GfImFreq diff: {:.4e}'.format(diff))

    def test_Gdiff_tau(self):

        B1 = BlockGf(name_list=['0', '1'], block_list=[self.G1_tau, self.G1_tau], make_copies=True)
        B2 = BlockGf(name_list=['0', '1'], block_list=[self.G2_tau, self.G2_tau], make_copies=True)

        diff = max_G_diff(B1, B2)
        print('GfImTime diff: {:.4e}'.format(diff))

if __name__ == '__main__':
    unittest.main()
