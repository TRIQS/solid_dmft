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
from itertools import product

from h5 import HDFArchive
from solid_dmft.postprocessing.eval_U_cRPA_RESPACK import read_interaction, fit_slater
from solid_dmft.postprocessing.eval_U_cRPA_Vasp import fit_slater_fulld as fit_slater_vasp

from triqs.operators.util.U_matrix import transform_U_matrix

import unittest


class test_respack_reader(unittest.TestCase):

    def setUp(self):
        self.path = './respack_sfo_data'
        self.seed = 'sfo'
        self.res = read_interaction(self.seed, self.path)
        return

    def test_U_J_fit(self):
        print('------------')
        print('Fitting screened interaction with fixed F4/F2 ratio')
        U_int, J_int = fit_slater(self.res.Uijij.real,
                                  self.res.Uijji.real,
                                  U_init=2, J_init=1,
                                  fixed_F4_F2=True)

        assert abs(U_int - 1.6087813069026353) < 1e-4
        assert abs(J_int - 0.6085825226244702) < 1e-4

        print('------------')
        print('Fitting screened Uijkl with Vasp routine with fixed F4/F2 ratio')
        Umat_swapped = np.moveaxis(self.res.Uijkl.real, 1, 2)
        U_int_vasp, J_int_vasp = fit_slater_vasp(Umat_swapped, 1, 2, 1)
        assert abs(U_int - U_int_vasp) < 1e-6
        assert abs(J_int - J_int_vasp) < 1e-6

        print('------------')
        print('Fitting bare interaction with fixed F4/F2 ratio')
        U_int, J_int = fit_slater(self.res.Vijij.real,
                                  self.res.Vijji.real,
                                  U_init=1.7, J_init=0.7,
                                  fixed_F4_F2=True)

        assert abs(U_int - 13.169980952573283) < 1e-4
        assert abs(J_int - 0.7852793988793887) < 1e-4

        return

    def test_U_J_fit_p(self):
        print('------------')
        print('Fitting screened interaction for p-shell with fixed F4/F2 ratio')
        print(self.res.Uijij.real[:3,:3])
        U_int, J_int = fit_slater(self.res.Uijij.real[:3,:3],
                                  self.res.Uijji.real[:3,:3],
                                  U_init=2, J_init=1,
                                  fixed_F4_F2=True)
        assert abs(U_int - 1.7803822) < 1e-4
        assert abs(J_int - 0.6158827) < 1e-4

    def test_F0_F2_F4_fit(self):
        print('------------')
        print('Fitting screened interaction to F0, F2, F4')
        U_int, J_int = fit_slater(self.res.U_R[(0, 0, 0)].real,
                                  self.res.J_R[(0, 0, 0)].real,
                                  U_init=13, J_init=1,
                                  fixed_F4_F2=False)

        assert abs(U_int - 1.6072376926602756) < 1e-4
        assert abs(J_int - 0.6166868032423679) < 1e-4

        return

    def test_construct_Uijkl(self):

        Uijkl = self.res.Uijkl

        rnb = range(self.res.n_orb)
        test_uijij = np.zeros((self.res.n_orb, self.res.n_orb), dtype=complex)
        test_uijji = np.zeros((self.res.n_orb, self.res.n_orb), dtype=complex)
        for i, j in product(rnb, rnb):
            test_uijij[i, j] = Uijkl[i, j, i, j]
            test_uijji[i, j] = Uijkl[i, j, j, i]

        assert np.allclose(test_uijij, self.res.U_R[(0, 0, 0)])
        assert np.allclose(test_uijji, self.res.J_R[(0, 0, 0)])

        Vijkl = self.res.Vijkl
        for i, j in product(rnb, rnb):
            test_uijij[i, j] = Vijkl[i, j, i, j]
            test_uijji[i, j] = Vijkl[i, j, j, i]

        assert np.allclose(test_uijij, self.res.V_R[(0, 0, 0)])
        assert np.allclose(test_uijji, self.res.X_R[(0, 0, 0)])


if __name__ == '__main__':
    unittest.main()
