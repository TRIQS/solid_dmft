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
import unittest

from triqs.gf import MeshDLRImFreq, Gf, BlockGf, make_gf_dlr, make_gf_imfreq

from solid_dmft.gw_embedding.bdft_converter import convert_screened_int, calc_Sigma_DC_gw, calc_W_from_Gloc, _get_dlr_from_IR


class test_gw_embedding(unittest.TestCase):

    def setUp(self):
        self.path = 'gw_embedding'
        self.seed = 'qpgw_svo'
        return

    def test_bdft_converter(self):

        convert_screened_int(seed=self.path+'/'+self.seed,
                             uloc_h5=self.path+'/uloc.scf.h5',
                             wloc_h5=self.path+'/wloc.scf.h5',
                             w_max=100,
                             IR_h5=self.path+'/IR_basis_1e4.1e-15.h5')

        return


    def test_calc_Sigma_DC_gw(self):

        with HDFArchive(self.path+'/'+self.seed+'.h5', 'r') as ar:
            Wloc_dlr = ar['DMFT_input']['Wloc_dlr'][0]
            Gloc_iw = ar['DMFT_input']['Gloc_iw'][0]
            Vloc = ar['DMFT_input']['Vloc'][0]

        mesh_dlr = MeshDLRImFreq(beta=Wloc_dlr.mesh.beta, statistic='Fermion', w_max=100, eps=1e-15)
        glist = []
        for block, gf in Gloc_iw:
            g_dlr_iw = Gf(mesh=mesh_dlr, target_shape=gf.target_shape)

            for iw in g_dlr_iw.mesh:
                g_dlr_iw[iw] = gf(iw)

            glist.append(make_gf_dlr(g_dlr_iw))

        Gloc_dlr = BlockGf(name_list=['up_0','down_0'], block_list=glist)

        Sig_DC_dlr, Sig_DC_hartree, Sig_DC_exchange = calc_Sigma_DC_gw(Wloc_dlr, Gloc_dlr, Vloc)

        DC_ref = np.array([[ 15.0495310-0.0j , 0.0+0.0j, 0.0+0.0j],
                           [0.0+0.0j,  15.0447285+0.0j      , 0.0+0.0j],
                           [0.0+0.0j, 0.0+0.0j, 15.0455924+0.0j]])
        for dc_mat in Sig_DC_hartree.values():
            assert np.allclose(dc_mat, DC_ref, rtol=1e-3, atol=1e-3)

    def test_calc_W(self):
        with HDFArchive(self.path+'/'+self.seed+'.h5', 'r') as ar:
            Uloc_dlr = ar['DMFT_input']['Uloc_dlr'][0]
            Gloc_iw = ar['DMFT_input']['Gloc_iw'][0]
            Vloc = ar['DMFT_input']['Vloc'][0]

        mesh_dlr = MeshDLRImFreq(beta=Uloc_dlr.mesh.beta, statistic='Fermion', w_max=100, eps=1e-15)
        glist = []
        for block, gf in Gloc_iw:
            g_dlr_iw = Gf(mesh=mesh_dlr, target_shape=gf.target_shape)

            for iw in g_dlr_iw.mesh:
                g_dlr_iw[iw] = gf(iw)

            glist.append(make_gf_dlr(g_dlr_iw))

        Gloc_dlr = BlockGf(name_list=['up_0','down_0'], block_list=glist)

        # create Coulomb tensor from Uloc_dlr
        Uloc_0 = make_gf_imfreq(Uloc_dlr['up_0'],1).data[0,:,:,:,:] + Vloc['up_0']
        U_dict = {'up_0' : Uloc_0, 'down_0' : Uloc_0}

        # test Gf | np.ndarray
        W_dlr = calc_W_from_Gloc(Gloc_dlr['up_0'], Uloc_0)
        # test BlockGf | np.ndarray
        W_dlr = calc_W_from_Gloc(Gloc_dlr, Uloc_0)
        # test BlockGf | dict of np.ndarray
        W_dlr = calc_W_from_Gloc(Gloc_dlr, U_dict)

        Sig_DC_dlr, Sig_DC_hartree, Sig_DC_exchange = calc_Sigma_DC_gw(W_dlr['up_0'], Gloc_dlr['up_0'], Uloc_0)
        DC_ref = np.array([[ 2.74430989-0.0j , 0.0+0.0j, 0.0+0.0j],
                           [0.0+0.0j,  2.74357053+0.0j      , 0.0+0.0j],
                           [0.0+0.0j, 0.0+0.0j, 2.74322484+0.0j]])
        assert np.allclose(Sig_DC_hartree, DC_ref, rtol=1e-3, atol=1e-3)




if __name__ == '__main__':
    unittest.main()
