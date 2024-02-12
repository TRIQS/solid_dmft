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

from os import stat
import numpy as np
from scipy.constants import physical_constants


from h5 import HDFArchive
from triqs.utility import mpi
from triqs.gf import (
    Gf,
    BlockGf,
    make_gf_dlr_imtime,
    make_gf_dlr,
    make_gf_imfreq,
    make_gf_dlr_imfreq,
)
from triqs.gf.meshes import MeshDLRImFreq, MeshDLRImTime
import itertools

from solid_dmft.gw_embedding.iaft import IAFT

HARTREE_EV = physical_constants['Hartree energy in eV'][0]

def _get_dlr_from_IR(Gf_ir, ir_kernel, mesh_dlr_iw, dim=2):
    n_orb = Gf_ir.shape[-1]
    stats = 'f' if mesh_dlr_iw.statistic == 'Fermion' else 'b'

    if stats == 'b':
        Gf_ir_pos = Gf_ir.copy()
        Gf_ir = np.zeros([Gf_ir_pos.shape[0] * 2 - 1] + [n_orb] * dim, dtype=complex)
        Gf_ir[: Gf_ir_pos.shape[0]] = Gf_ir_pos[::-1]
        Gf_ir[Gf_ir_pos.shape[0] :] = Gf_ir_pos[1:]

    Gf_dlr_iw = Gf(mesh=mesh_dlr_iw, target_shape=[n_orb] * dim)

    # prepare idx array for spare ir
    if stats == 'f':
        mesh_dlr_iw_idx = np.array([iwn.index for iwn in mesh_dlr_iw])
    else:
        mesh_dlr_iw_idx = np.array([iwn.index for iwn in mesh_dlr_iw])

    Gf_dlr_iw.data[:] = ir_kernel.w_interpolate(Gf_ir, mesh_dlr_iw_idx, stats=stats, ir_notation=False)

    Gf_dlr = make_gf_dlr(Gf_dlr_iw)
    return Gf_dlr


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

        return (
            BlockGf(name_list=list(Gloc_dlr.indices), block_list=Sig_DC_dlr_list),
            Sig_DC_hartree_list,
            Sig_DC_exchange_list,
        )

    n_orb = Gloc_dlr.target_shape[0]

    # dynamic part
    Gloc_dlr_t = make_gf_dlr_imtime(Gloc_dlr)
    Sig_dlr_t = Gf(mesh=Gloc_dlr_t.mesh, target_shape=Gloc_dlr_t.target_shape)

    Wloc_dlr_t = make_gf_dlr_imtime(Wloc_dlr)

    for tau in Gloc_dlr_t.mesh:
        # Wloc_dlr is bosonic and the mesh has a different hash, use call to get value at tau point
        Sig_dlr_t[tau] = -1 * np.einsum('ijkl, jk -> li', Wloc_dlr_t[tau], Gloc_dlr_t[tau])

    Sig_DC_dlr = make_gf_dlr(Sig_dlr_t)

    # static hartree Part
    Sig_DC_hartree = np.zeros((n_orb, n_orb))
    Sig_DC_hartree = 2 * np.einsum('ijkl, lj -> ik', Vloc, Gloc_dlr.density())

    if verbose:
        print('static Hartree part of DC')
        print(Sig_DC_hartree.real)
        if np.any(np.imag(Sig_DC_hartree) > 1e-3):
            print('Im:')
            print(np.imag(Sig_DC_hartree))

    # static exchange part
    Sig_DC_exchange = np.zeros((n_orb, n_orb))
    Sig_DC_exchange = -1 * np.einsum('ijkl, jk -> li', Vloc, Gloc_dlr.density())

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
    mesh_bos = MeshDLRImTime(
        beta=Gloc_dlr.mesh.beta,
        statistic='Boson',
        w_max=Gloc_dlr.mesh.w_max,
        eps=Gloc_dlr.mesh.eps,
    )

    PI_dlr_t = Gf(mesh=mesh_bos, target_shape=[nb] * 4)
    for tau in Gloc_dlr_t.mesh:
        PI_dlr_t[tau] = -2 * np.einsum('bl, ka -> abkl', Gloc_dlr_t[tau], Gloc_dlr(Gloc_dlr_t.mesh.beta - tau))

    PI_dlr = make_gf_dlr(PI_dlr_t)
    PI_dlr_w = make_gf_dlr_imfreq(PI_dlr)

    # need to swap indices and go into product basis
    U_prod = np.transpose(U, (2, 0, 1, 3)).reshape(nb**2, nb**2)

    W_dlr_w = Gf(mesh=PI_dlr_w.mesh, target_shape=[nb] * 4)

    ones = np.eye(nb**2)
    for w in PI_dlr_w.mesh:
        eps = ones - U_prod @ PI_dlr_w[w].reshape(nb**2, nb**2)
        # in product basis W_ji,kl
        W_dlr_w[w] = (np.linalg.inv(eps) @ U_prod - U_prod).reshape(nb, nb, nb, nb)

        # swap indices back
        W_dlr_w[w] = np.transpose(W_dlr_w[w], (1, 2, 0, 3))
    W_dlr = make_gf_dlr(W_dlr_w)

    return W_dlr


def convert_gw_output(job_h5: str, gw_h5: str, wmax_dlr: float, iter: int = 0) -> tuple[dict, IAFT]:
    """
    read bdft output and convert to triqs Gf DLR objects

    Parameters
    ----------
    job_h5: string
        path to solid_dmft job file
    gw_h5: string
        path to GW checkpoint file for AIMBES code
    wmax_dlr: float
        DLR energy cutoff, same as Lambda / beta for the impurity problem
    iter: int, optional
        iteration to read from gw_h5 calculation, defaults to last iteration

    Returns
    -------
    gw_data: dict
        dictionary holding all read objects: mu_emb, beta, lam, w_max, prec, mesh_dlr_iw_b,
        mesh_dlr_iw_f, n_orb, G0_dlr, Gloc_dlr, Sigma_imp_dlr, Sigma_imp_DC_dlr, Uloc_dlr,
        Vloc, Hloc0, Vhf_dc, Vhf
    ir_kernel: sparse_ir kernel object
        IR kernel with AIMBES paramaters
    """

    mpi.report('reading output from bdft code')

    gw_data = {'wmax_dlr': wmax_dlr}

    with HDFArchive(gw_h5, 'r') as ar:
        if not iter:
            iter_1e = ar['downfold_1e/final_iter']
            iter_2e = ar['downfold_2e/final_iter']
        else:
            iter_1e = iter
            iter_2e = iter

        # auxilary quantities
        gw_data['mu_emb'] = ar[f'downfold_1e/iter{iter_1e}']['mu']
        gw_data['beta'] = ar['imaginary_fourier_transform']['beta']
        gw_data['lam'] = ar['imaginary_fourier_transform']['lambda']
        gw_data['w_max'] = gw_data['lam'] / gw_data['beta']
        gw_data['number_of_spins'] = ar['system/number_of_spins']
        assert gw_data['number_of_spins'] == 1, 'spin calculations not yet supported in converter'

        prec = ar['imaginary_fourier_transform']['prec']
        if prec == 'high':
            gw_data['prec'] = 1e-15
        elif prec == 'mid':
            gw_data['prec'] = 1e-10
        elif prec == 'low':
            gw_data['prec'] = 1e-6

        # 1 particle properties
        g_weiss_wsIab = ar[f'downfold_1e/iter{iter_1e}']['g_weiss_wsIab']
        Sigma_dc_wsIab = ar[f'downfold_1e/iter{iter_1e}']['Sigma_dc_wsIab']
        Gloc = ar[f'downfold_1e/iter{iter_1e}']['Gloc_wsIab']
        gw_data['n_inequiv_shells'] = Gloc.shape[1]

        # 2 particle properties
        # TODO: discuss how the site index is used right now in bDFT
        Vloc_jk = ar[f'downfold_2e/iter{iter_2e}']['Vloc_abcd']
        Uloc_ir_jk = ar[f'downfold_2e/iter{iter_2e}']['Uloc_wabcd'][:, ...]
        # switch inner two indices to match triqs notation
        Vloc = np.zeros(Vloc_jk.shape, dtype=complex)
        Uloc_ir = np.zeros(Uloc_ir_jk.shape, dtype=complex)
        n_orb = Vloc.shape[0]
        for or1, or2, or3, or4 in itertools.product(range(n_orb), repeat=4):
            Vloc[or1, or2, or3, or4] = Vloc_jk[or1, or3, or2, or4]
            for ir_w in range(Uloc_ir_jk.shape[0]):
                Uloc_ir[ir_w, or1, or2, or3, or4] = Uloc_ir_jk[ir_w, or1, or3, or2, or4]

        # get Hloc_0
        Vhf_dc_sIab = ar[f'downfold_1e/iter{iter_1e}']['Vhf_dc_sIab'][0, 0]
        Vhf_sIab = ar[f'downfold_1e/iter{iter_1e}']['Vhf_gw_sIab'][0, 0]
        H0_loc = ar[f'downfold_1e/iter{iter_1e}']['H0_sIab']

        if 'Vcorr_gw_sIab' in ar[f'downfold_1e/iter{iter_1e}']:
            mpi.report('Found Vcorr_sIab in the bdft checkpoint file, '
                       'i.e. Embedding on top of an effective QP Hamiltonian.')
            Vcorr_sIab = ar[f'downfold_1e/iter{iter_1e}/Vcorr_gw_sIab']
            Vcorr_dc_sIab = ar[f'downfold_1e/iter{iter_1e}/Vcorr_dc_sIab']
            Hloc0 = -1*(np.eye(n_orb) * gw_data['mu_emb'] - H0_loc[0,0] - Vhf_sIab - Vcorr_sIab + Vhf_dc_sIab + Vcorr_dc_sIab)
            qp_emb = True
        else:
            Sigma_wsIab = ar[f'downfold_1e/iter{iter_1e}']['Sigma_gw_wsIab']
            qp_emb = False
            Hloc0 = -1*(np.eye(n_orb) * gw_data['mu_emb'] - H0_loc[0,0] - (Vhf_sIab-Vhf_dc_sIab))

    # get IR object
    mpi.report('create IR kernel and convert to DLR')
    # create IR kernel
    ir_kernel = IAFT(beta=gw_data['beta'], lmbda=gw_data['lam'], prec=gw_data['prec'])

    gw_data['mesh_dlr_iw_b'] = MeshDLRImFreq(
        beta=gw_data['beta'],
        statistic='Boson',
        w_max=gw_data['wmax_dlr'],
        eps=gw_data['prec'],
    )
    gw_data['mesh_dlr_iw_f'] = MeshDLRImFreq(
        beta=gw_data['beta'],
        statistic='Fermion',
        w_max=gw_data['wmax_dlr'],
        eps=gw_data['prec'],
    )

    (
        U_dlr_list,
        G0_dlr_list,
        Gloc_dlr_list,
        Sigma_dlr_list,
        Sigma_DC_dlr_list,
        V_list,
        Hloc_list,
        Vhf_list,
        Vhf_dc_list,
        n_orb_list,
    ) = [], [], [], [], [], [], [], [], [], []
    for ish in range(gw_data['n_inequiv_shells']):
        # fit IR Uloc on DLR iw mesh
        temp = _get_dlr_from_IR(Uloc_ir, ir_kernel, gw_data['mesh_dlr_iw_b'], dim=4)
        Uloc_dlr = BlockGf(name_list=['up_0', 'down_0'], block_list=[temp, temp], make_copies=True)

        U_dlr_list.append(Uloc_dlr)
        V_list.append({'up_0': Vloc.copy(), 'down_0': Vloc})
        Hloc_list.append({'up_0': Hloc0.copy(), 'down_0': Hloc0})
        Vhf_list.append({'up_0': Vhf_sIab.copy(), 'down_0': Vhf_sIab})
        Vhf_dc_list.append({'up_0': Vhf_dc_sIab.copy(), 'down_0': Vhf_dc_sIab})
        n_orb_list.append(n_orb)

        temp = _get_dlr_from_IR(g_weiss_wsIab[:, 0, ish, :, :], ir_kernel, gw_data['mesh_dlr_iw_f'], dim=2)
        G0_dlr = BlockGf(name_list=['up_0', 'down_0'], block_list=[temp, temp], make_copies=True)
        G0_dlr_list.append(G0_dlr)

        temp = _get_dlr_from_IR(Gloc[:, 0, ish, :, :], ir_kernel, gw_data['mesh_dlr_iw_f'], dim=2)
        Gloc_dlr = BlockGf(name_list=['up_0', 'down_0'], block_list=[temp, temp], make_copies=True)
        Gloc_dlr_list.append(Gloc_dlr)

        # since Sigma can have a static shift we return DLR Imfreq mesh
        if not qp_emb:
            temp = _get_dlr_from_IR(Sigma_wsIab[:, 0, ish, :, :], ir_kernel, gw_data['mesh_dlr_iw_f'], dim=2)
            Sigma_dlr = BlockGf(name_list=['up_0', 'down_0'], block_list=[temp, temp], make_copies=True)
            Sigma_dlr_list.append(Sigma_dlr)

        temp = _get_dlr_from_IR(Sigma_dc_wsIab[:, 0, ish, :, :], ir_kernel, gw_data['mesh_dlr_iw_f'], dim=2)
        Sigma_DC_dlr = BlockGf(name_list=['up_0', 'down_0'], block_list=[temp, temp], make_copies=True)
        Sigma_DC_dlr_list.append(Sigma_DC_dlr)

    gw_data['G0_dlr'] = G0_dlr_list
    gw_data['Gloc_dlr'] = Gloc_dlr_list
    gw_data['Sigma_imp_dlr'] = Sigma_dlr_list
    gw_data['Sigma_imp_DC_dlr'] = Sigma_DC_dlr_list
    gw_data['Uloc_dlr'] = U_dlr_list
    gw_data['Vloc'] = V_list
    gw_data['Hloc0'] = Hloc_list
    gw_data['Vhf_dc'] = Vhf_dc_list
    gw_data['Vhf'] = Vhf_list
    gw_data['n_orb'] = n_orb_list

    # write Uloc / Wloc back to h5 archive
    mpi.report(f'writing results in {job_h5}/DMFT_input')

    with HDFArchive(job_h5, 'a') as ar:
        if 'DMFT_results' in ar and 'iteration_count' in ar['DMFT_results']:
            it = ar['DMFT_results']['iteration_count'] + 1
        else:
            it = 1
        if 'DMFT_input' not in ar:
            ar.create_group('DMFT_input')
        if f'iter{it}' not in ar['DMFT_input']:
            ar['DMFT_input'].create_group(f'iter{it}')

        for key, value in gw_data.items():
            ar[f'DMFT_input/iter{it}'][key] = value

    return gw_data, ir_kernel


