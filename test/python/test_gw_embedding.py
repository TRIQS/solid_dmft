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
import unittest

from triqs.gf import MeshDLRImFreq, Gf, BlockGf, make_gf_dlr, make_gf_imfreq, make_gf_dlr_imfreq

from solid_dmft.gw_embedding.bdft_converter import convert_gw_output, calc_Sigma_DC_gw, calc_W_from_Gloc


class test_gw_embedding(unittest.TestCase):

    def setUp(self):
        self.path = 'gw_embedding'
        self.seed = 'qpgw_svo'
        return

    def test_bdft_converter(self):

        gw_data, ir_kernel = convert_gw_output(job_h5='emb_test.h5', gw_h5='svo_gw_emb_stat/inp.h5', wmax_dlr=10)

        return


    def test_calc_W_DC(self):
        with HDFArchive('emb_test.h5', 'r') as ar:
            Uloc_dlr = ar['DMFT_input/iter1']['Uloc_dlr'][0]
            Gloc_dlr = ar['DMFT_input/iter1']['Gloc_dlr'][0]
            Vloc = ar['DMFT_input/iter1']['Vloc'][0]

        mesh_dlr = MeshDLRImFreq(beta=Uloc_dlr.mesh.beta, statistic='Fermion', w_max=Gloc_dlr.mesh.w_max, eps=Gloc_dlr.mesh.eps)
        glist = []
        Gloc_iw = make_gf_dlr_imfreq(Gloc_dlr)
        for block, gf in Gloc_iw:
            g_dlr_iw = Gf(mesh=mesh_dlr, target_shape=gf.target_shape)

            for iw in g_dlr_iw.mesh:
                g_dlr_iw[iw] = gf[iw]

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
        DC_ref = np.array([[1.13031371e-01-0.0j , 0.0+0.0j, 0.0+0.0j],
                           [0.0+0.0j,  1.13030001e-01+0.0j      , 0.0+0.0j],
                           [0.0+0.0j, 0.0+0.0j, 1.13031969e-01+0.0j]])
        assert np.allclose(Sig_DC_hartree, DC_ref, rtol=1e-3, atol=1e-3)


if __name__ == '__main__':
    unittest.main()
