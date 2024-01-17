# -*- coding: utf-8 -*-
################################################################################
#
# solid_dmft - A versatile python wrapper to perform DFT+DMFT calculations
#              utilizing the TRIQS software library
#
# Copyright (C) 2018-2020, ETH Zurich
# Copyright (C) 2021, The Simons Foundation
#      authors: A. Hampel, M. Merkel, and S. Beck
#
# solid_dmft is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# solid_dmft is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# solid_dmft (in the file COPYING.txt in this directory). If not, see
# <http://www.gnu.org/licenses/>.
#
################################################################################
"""
converter from bdft output to edmft input for solid_dmft
"""

import numpy as np
from scipy.constants import physical_constants


from h5 import HDFArchive
from triqs.utility import mpi
from triqs.gf import Gf, BlockGf, make_gf_dlr_imtime, make_gf_dlr, make_gf_imfreq, make_gf_dlr_imfreq
from triqs.gf.meshes import MeshDLRImFreq, MeshDLRImTime
import itertools


def _get_dlr_Wloc_from_IR(Wloc_IR, IR_mesh_pos, beta, w_max, eps=1e-15):

    n_orb = Wloc_IR.shape[1]

    mesh_dlr = MeshDLRImFreq(beta=beta, statistic='Boson', w_max=w_max, eps=eps, symmetrize=False)
    mesh_dlr_arr = np.array([iw.imag for iw in mesh_dlr])
    mid_idx = int(np.where(mesh_dlr_arr == 0.0)[0])
    Wloc_dlr = Gf(mesh=mesh_dlr, target_shape=[n_orb, n_orb, n_orb, n_orb])

    for or1, or2, or3, or4 in itertools.product(range(n_orb), repeat=4):
        Wloc_dlr.data[mid_idx:, or1, or2, or3, or4] = np.interp(
            mesh_dlr_arr[mid_idx:], IR_mesh_pos, Wloc_IR[:, or1, or3, or2, or4])

        # fill negative mesh values, note that the DLR mesh is not symmetric
        Wloc_dlr.data[0:mid_idx, or1, or2, or3, or4] = np.interp(
            mesh_dlr_arr[:mid_idx], -1*IR_mesh_pos[::-1], Wloc_IR[:, or1, or3, or2, or4][::-1])

    return Wloc_dlr, mesh_dlr_arr


def calc_Sigma_DC_gw(Wloc_dlr, Gloc_dlr, Vloc, verbose=False):

    if isinstance(Gloc_dlr, BlockGf):
        Sig_DC_dlr_list = []
        Sig_DC_hartree_list = {}
        Sig_DC_exchange_list = {}
        for block, gloc in Gloc_dlr:
            res = calc_Sigma_DC_gw(Wloc_dlr[block], gloc, Vloc[block], verbose)
            Sig_DC_dlr_list.append(res[0])
            Sig_DC_hartree_list[block] = res[1]
            Sig_DC_exchange_list[block] = res[2]

        return BlockGf(name_list=list(Gloc_dlr.indices), block_list=Sig_DC_dlr_list), Sig_DC_hartree_list, Sig_DC_exchange_list

    n_orb = Gloc_dlr.target_shape[0]

    # dynamic part
    Gloc_dlr_t = make_gf_dlr_imtime(Gloc_dlr)
    Sig_dlr_t = Gf(mesh=Gloc_dlr_t.mesh, target_shape=Gloc_dlr_t.target_shape)

    Wloc_dlr_t = make_gf_dlr_imtime(Wloc_dlr)

    for tau in Gloc_dlr_t.mesh:
        # Wloc_dlr is bosonic and the mesh has a different hash, use call to get value at tau point
        Sig_dlr_t[tau] = -1*np.einsum('ijkl, jk -> li', Wloc_dlr_t[tau], Gloc_dlr_t[tau])

    Sig_DC_dlr = make_gf_dlr(Sig_dlr_t)

    # static hartree Part
    Sig_DC_hartree = np.zeros((n_orb, n_orb))
    Sig_DC_hartree = 2*np.einsum('ijkl, lj -> ik', Vloc, Gloc_dlr.density())

    if verbose:
        print('static Hartree part of DC')
        print(Sig_DC_hartree.real)
        if np.any(np.imag(Sig_DC_hartree) > 1e-3):
            print('Im:')
            print(np.imag(Sig_DC_hartree))

    # static exchange part
    Sig_DC_exchange = np.zeros((n_orb, n_orb))
    Sig_DC_exchange = -1*np.einsum('ijkl, jk -> li', Vloc, Gloc_dlr.density())

    if verbose:
        print('static exchange part of DC')
        print(Sig_DC_exchange.real)
        if np.any(np.imag(Sig_DC_exchange) > 1e-3):
            print('Im:')
            print(np.imag(Sig_DC_exchange))
    return Sig_DC_dlr, Sig_DC_hartree, Sig_DC_exchange

def calc_W_from_Gloc(Gloc_dlr: Gf | BlockGf, U: np.ndarray | dict) -> Gf | BlockGf:
    r"""

    Calculate Wijkl from given constant U tensor and Gf on DLRMesh


    triqs notation for Uijkl = phi*_i(r) phi*_j(r') U(r,r') phi_l'(r') phi_k(r)
                            = Uijkl c^+_i c^+_j' c_l' c_k

    where the ' denotes a spin index different from the other without '

    the according diagram is (left and right have same spin):

       j (phi)         k' (phi)
         \              /
          <            <
           \__________/
           /          \
          >            >
         /              \
       i (phi*)          l'

    we now have to move to a product basis form to combine two indices
    i.e. go from nb,nb,nb,nb to nb**2,nb**2 tensors:

        Uji,kl = phi*_i(r) phi_j(r) U(r,r') phi*_k(r') phi_l(r')
               = Psi*_ji(r) U(r,r') Psi_kl(r')

    So we have to transform the triqs notation of Uijkl -> Uki,jl, i.e.
    swap col/rows as (2,0,1,3) to go to the basis and the in the end
    swap W_ki,jl back in reverse.

    Then we compute pubble polarizability as

    Pi_ab,kl(tau) = -2 G_bl(tau) G_ka(beta - tau)

    So that

    [ U Pi(iwn) ]_ji,kl = sum_ab U_ji,ab Pi_ab,kl(iwn)

    i.e.
       j'              a ___
         \              /   \ k
          <            <     \
           \__________/       \
           /          \       /
          >            >     /
         /              \___/ l
       i'               b

    then the screened Coulomb interaction in product basis is:

    W_ji,kl(iwn) = [1 - U Pi(iwn) ]^-1_ji,kl Uji,kl - Uji,kl

    (subtract static shift here), and finally convert back to triqs notation.


    Parameters
    ----------
    Gloc_dlr : BlockGf or Gf with MeshDLR

    U : np.ndarray of with shape [Gloc_dlr.target_shape]*4 or dict of np.ndarray

    Returns
    -------
    W_dlr : BlockGf or Gf
        screened Coulomb interaction
    """

    if isinstance(Gloc_dlr, BlockGf):
        Wloc_list = []
        for block, gloc in Gloc_dlr:
            if isinstance(U, np.ndarray):
                Wloc_list.append(calc_W_from_Gloc(gloc, U))
            else:
                Wloc_list.append(calc_W_from_Gloc(gloc, U[block]))

        return BlockGf(name_list=list(Gloc_dlr.indices), block_list=Wloc_list)

    nb = Gloc_dlr.target_shape[0]
    Gloc_dlr_t = make_gf_dlr_imtime(Gloc_dlr)
    mesh_bos = MeshDLRImTime(beta=Gloc_dlr.mesh.beta, statistic = 'Boson', w_max=Gloc_dlr.mesh.w_max, eps=Gloc_dlr.mesh.eps)

    PI_dlr_t = Gf(mesh=mesh_bos, target_shape=[nb]*4)
    for tau in Gloc_dlr_t.mesh:
        PI_dlr_t[tau] = -2*np.einsum('bl, ka -> abkl', Gloc_dlr_t[tau], Gloc_dlr(Gloc_dlr_t.mesh.beta - tau))

    PI_dlr = make_gf_dlr(PI_dlr_t)
    PI_dlr_w = make_gf_dlr_imfreq(PI_dlr)

    # need to swap indices and go into product basis
    U_prod = np.transpose(U, (2, 0, 1, 3)).reshape(nb**2, nb**2)

    W_dlr_w = Gf(mesh=PI_dlr_w.mesh, target_shape=[nb]*4)

    ones =  np.eye(nb**2)
    for w in PI_dlr_w.mesh:
        eps = ones - U_prod@PI_dlr_w[w].reshape(nb**2, nb**2)
        # in product basis W_ji,kl
        W_dlr_w[w] = (np.linalg.inv(eps)@U_prod - U_prod).reshape(nb, nb, nb, nb)

        # swap indices back
        W_dlr_w[w] = np.transpose(W_dlr_w[w] , (1,2,0,3))
    W_dlr = make_gf_dlr(W_dlr_w)

    return W_dlr

def convert_screened_int(seed, uloc_h5, wloc_h5, w_max, IR_h5, beta=None, iter=None, eps=1e-15):
    """
    read bdft output and convert to solid_dmft readable h5 input
    """

    Hartree_eV = physical_constants['Hartree energy in eV'][0]

    mpi.report('reading input')

    with HDFArchive(seed+'.h5', 'r') as ar:
        # TODO implement multiple sites later
        shells = ar['dft_input/shells']
        n_shells = ar['dft_input/n_shells']

    with HDFArchive(uloc_h5, 'r') as ar:
        if iter is None:
            iter = ar['scf/final_iter']
        # TODO: discuss how the site index is used right now in bDFT
        Vloc_jk = ar['scf'][f'iter{iter}']['embed']['Vloc'] * Hartree_eV
        # switch inner two indices to match triqs notation
        Vloc = np.zeros(Vloc_jk.shape, dtype=complex)
        n_orb = Vloc.shape[0]
        for or1, or2, or3, or4 in itertools.product(range(n_orb), repeat=4):
            Vloc[or1,or2,or3,or4] = Vloc_jk[or1,or3,or2,or4]
        Uloc_ir = ar['scf'][f'iter{iter}']['embed']['Wloc'][:, ...] * Hartree_eV
        mu_evgw = ar['scf'][f'iter{iter}']['mu'] * Hartree_eV

        beta_uloc = ar['imaginary_fourier_transform']['beta'] / Hartree_eV
        if beta is None:
            mpi.report(f'no beta given using {beta_uloc} from bDFT calculation')
            beta = beta_uloc
        else:
            assert np.isclose(beta_uloc, beta, atol=1e-4), 'beta does not match beta in bDFT'

    # TODO should be one archive but for now one has to run 2 separate bDFT calculations
    with HDFArchive(wloc_h5, 'r') as ar:
        Wloc_ir = ar['scf'][f'iter{iter}']['embed']['Wloc'][:, ...] * Hartree_eV

    # get IR bosonic mesh (move to uloc read once avail)
    with HDFArchive(IR_h5, 'r') as ir_ar:
        iw_mesh_ir_idx = ir_ar['boson']['wn_mesh']

    # convert mesh to pos values to match Uloc / Wloc
    iw_mesh_ir_pos = []

    for indx in iw_mesh_ir_idx[np.where(iw_mesh_ir_idx == 0.0)[0][0]:]:
        iw_mesh_ir_pos.append((indx*np.pi)/beta)
    iw_mesh_ir_pos = np.array(iw_mesh_ir_pos)

    mpi.report('fitting Wloc and Uloc on DLR mesh')

    # rotate to sumk structure can be multiple sites
    U_dlr_list, W_dlr_list, V_list = [], [], []
    for ish in range(n_shells):
        # fit IR Wloc on DLR iw mesh
        wloc, _ = _get_dlr_Wloc_from_IR(Wloc_ir,
                                               iw_mesh_ir_pos,
                                               beta=beta,
                                               w_max=w_max,
                                               eps=eps)
        Wloc_dlr_iw = BlockGf(name_list=['up_0','down_0'], block_list=[wloc, wloc], make_copies=True)
        # create DLR Gf
        Wloc_dlr = make_gf_dlr(Wloc_dlr_iw)

        # fit IR Uloc on DLR iw mesh
        uloc, _ = _get_dlr_Wloc_from_IR(Uloc_ir,
                                               iw_mesh_ir_pos,
                                               beta=beta,
                                               w_max=w_max,
                                               eps=eps)
        Uloc_dlr_iw = BlockGf(name_list=['up_0','down_0'], block_list=[uloc, uloc], make_copies=True)
        # create DLR Gf
        Uloc_dlr = make_gf_dlr(Uloc_dlr_iw)

        U_dlr_list.append(Uloc_dlr)
        W_dlr_list.append(Wloc_dlr)
        V_list.append({'up_0': Vloc.copy(), 'down_0': Vloc})


    # write Uloc / Wloc back to h5 archive
    mpi.report(f'writing results in {seed}.h5/DMFT_input')
    with HDFArchive(seed+'.h5', 'a') as ar:
        if 'DMFT_input' not in ar:
            ar.create_group('DMFT_input')
        ar['DMFT_input']['Vloc'] = V_list
        ar['DMFT_input']['beta'] = beta
        ar['DMFT_input']['Uloc_dlr'] = U_dlr_list
        ar['DMFT_input']['Wloc_dlr'] = W_dlr_list
        ar['DMFT_input']['mu_evgw'] = mu_evgw

    return

