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
from triqs.gfs import (
    Gf,
    BlockGf,
    make_gf_dlr_imtime,
    make_gf_dlr,
    make_gf_dlr_imfreq,
)
from triqs.mesh import MeshDLRImFreq, MeshDLRImTime
import itertools

from solid_dmft.gw_embedding.iaft import IAFT

HARTREE_EV = physical_constants['Hartree energy in eV'][0]

def _get_dlr_from_IR(Gf_ir, ir_kernel, mesh_dlr_iw, dim=2):
    r"""
    Interpolate a given Gf from IR mesh to DLR mesh

    Parameters
    ----------
    Gf_ir : np.ndarray
        Green's function on IR mesh
    ir_kernel : sparse_ir
        IR kernel object
    mesh_dlr_iw : MeshDLRImFreq
        DLR mesh
    dim : int, optional
        dimension of the Green's function, defaults to 2

    Returns
    -------
    Gf_dlr : BlockGf or Gf
        Green's function on DLR mesh
    """

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


def check_iaft_accuracy(Aw, ir_kernel, stats,
                        beta, dlr_wmax, dlr_prec, data_name):
    mpi.report('============== DLR mesh check ==============\n')
    mpi.report(f'Estimating accuracy of the user-defined (wmax, eps) = '
               f'({dlr_wmax}, {dlr_prec}) for the DLR mesh\n')
    ir_imp_kernel = IAFT(beta=beta, lmbda=beta * dlr_wmax, prec=dlr_prec)
    Aw_imp = ir_kernel.w_interpolate(Aw, ir_imp_kernel.wn_mesh('f'), 'f')

    ir_imp_kernel.check_leakage(Aw_imp, stats, data_name, w_input=True)
    mpi.report('=================== done ===================\n')


def estimate_zero_moment(Aw, iw_mesh):
    iw_m1 = iw_mesh[-1]
    iw_m2 = iw_mesh[-2]
    t = Aw[-1].real - (Aw[-1] - Aw[-2]).real * iw_m2 ** 2 / (
           iw_m2 ** 2 - iw_m1 ** 2)
    t = t.astype(complex)

    return t


def extract_t_and_delta(g_weiss_wsIab, ir_kernel):
    iwn_mesh_imp = ir_kernel.wn_mesh('f') * np.pi / ir_kernel.beta

    ns, nImp = g_weiss_wsIab.shape[1:3]
    weiss_tmp = np.zeros(g_weiss_wsIab.shape, dtype=complex)
    for n, g in enumerate(g_weiss_wsIab):
        for s in range(ns):
            for I in range(nImp):
                weiss_tmp[n, s, I] = 1j * iwn_mesh_imp[n] - np.linalg.inv(g[s, I])

    t_sIab_estimate = estimate_zero_moment(weiss_tmp, iwn_mesh_imp)

    # construct delta: delta(iw) = iw - t - g^{-1}(iw)
    delta_estimate = np.zeros(g_weiss_wsIab.shape, dtype=complex)
    nbnd = t_sIab_estimate.shape[-1]
    for n in range(delta_estimate.shape[0]):
        for s in range(ns):
            for I in range(nImp):
                g_weiss_inv = np.linalg.inv(g_weiss_wsIab[n, 0, 0])
                delta_estimate[n, s, I] = 1j * iwn_mesh_imp[n] * np.eye(nbnd) - t_sIab_estimate[s, I] - g_weiss_inv
    ir_kernel.check_leakage(delta_estimate, 'f', 'delta_estimate', w_input=True)

    return t_sIab_estimate, delta_estimate


def read_t_and_delta(aimb_h5, it_1e=None):
    mpi.report(f"Reading the analytic H_loc0 and the corresponding hybridization from aimbes h5 {aimb_h5}.")
    with HDFArchive(aimb_h5, 'r') as ar:
        if not it_1e:
            it_1e = ar['downfold_1e/final_iter']
        iter_grp = ar[f'downfold_1e/iter{it_1e}']

        t_sIab = iter_grp['H0_sIab'] + iter_grp['Vhf_gw_sIab'] - iter_grp['Vhf_dc_sIab']
        if 'Vcorr_gw_sIab' in iter_grp:
            t_sIab += (iter_grp['Vcorr_gw_sIab'] - iter_grp['Vcorr_dc_sIab'])

        nspin, nImp, nOrbs = t_sIab.shape[:3]
        for s in np.arange(nspin):
            for I in range(nImp):
                t_sIab[s, I] -= np.eye(nOrbs) * iter_grp["mu"]

        delta_wsIab = iter_grp['delta_wsIab']

    return t_sIab, delta_wsIab


def tail_fit_g_weiss(g_weiss_wsIab, ir_kernel, gw_data, wmax_imp=None, eps_imp=None):
    mpi.report("Extracting H_loc0 and hybridization from tail fitting fermionic Weiss field g.")
    # if user-defined wmax and eps for the impurity
    if wmax_imp is not None or eps_imp is not None:
        ir_imp_kernel = IAFT(beta=gw_data['beta'],
                             lmbda=gw_data['beta'] * gw_data['gw_dlr_wmax'],
                             prec=gw_data['gw_ir_prec'], verbal=False)
        g_imp = ir_kernel.w_interpolate(g_weiss_wsIab, ir_imp_kernel.wn_mesh('f'), 'f')
        ir_imp_kernel.check_leakage(g_imp, 'f', "fermionic Weiss field on the customized imaginary meshes",
                                    w_input=True)
    else:
        ir_imp_kernel = ir_kernel
        g_imp = g_weiss_wsIab

    mpi.report("")
    # extracting H0 and Delta from g_weiss
    t_sIab, delta_wsIab = extract_t_and_delta(g_imp, ir_imp_kernel)

    return t_sIab, delta_wsIab, ir_imp_kernel


def bath_fitting(delta_wsIab, iw_mesh, Np=5):
    mpi.report(f"Bath fitting for hybridization with nbath/impurity orbital = {Np}")
    try:
        from adapol import hybfit
    except ImportError:
        raise ImportError("bath fitting requires the adapol package (https://github.com/flatironinstitute/adapol). \n"
                          "Please ensure that it is installed. ")

    nspin, nImp = delta_wsIab.shape[1:3]
    error = -1
    for s in np.arange(nspin):
        for I in range(nImp):
            bath_energy, bath_hyb, fit_error, func = hybfit(delta_wsIab[:, s, I], 1j*iw_mesh,
                                                            Np=5, solver='sdp', verbose=False)
            delta_wsIab[:, s, I] = func(1j*iw_mesh)
            error = max(error, abs(fit_error))
    mpi.report(f"Bath fitting error =  {error}")


def calc_Sigma_DC_gw(Wloc_dlr, Gloc_dlr, Vloc, verbose=False):
    r"""
    Calculate the double counting part of the self-energy from the screened Coulomb interaction

    Parameters
    ----------
    Wloc_dlr : BlockGf or Gf with MeshDLR
        screened Coulomb interaction
    Gloc_dlr : BlockGf or Gf with MeshDLR
        local Green's function
    Vloc : np.ndarray
        local Coulomb interaction
    verbose : bool, optional
        print additional information, defaults to False

    Returns
    -------
    Sig_DC_dlr : BlockGf or Gf
        double counting part of the self-energy
    Sig_DC_hartree : np.ndarray
        static Hartree part of the self-energy
    Sig_DC_exchange : np.ndarray
        static exchange part of the self-energy
    """

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

    i_tau = 0
    for tau in Gloc_dlr_t.mesh:
        # Wloc_dlr is bosonic and the mesh has a different hash, use call to get value at tau point
        Sig_dlr_t[tau] = -1 * np.einsum('ijkl, jk -> li', Wloc_dlr_t.data[i_tau,:,:], Gloc_dlr_t[tau])
        i_tau += 1

    Sig_DC_dlr = make_gf_dlr(Sig_dlr_t)

    # static hartree Part
    Sig_DC_hartree = np.zeros((n_orb, n_orb))
    Sig_DC_hartree = 2 * np.einsum('ijkl, lj -> ik', Vloc, Gloc_dlr.density())
    # symmetrize
    Sig_DC_hartree = 0.5 * (Sig_DC_hartree + Sig_DC_hartree.conj().T)

    if verbose:
        print('static Hartree part of DC')
        print(Sig_DC_hartree.real)
        if np.any(np.imag(Sig_DC_hartree) > 1e-3):
            print('Im:')
            print(np.imag(Sig_DC_hartree))

    # static exchange part
    Sig_DC_exchange = np.zeros((n_orb, n_orb))
    Sig_DC_exchange = -1 * np.einsum('ijkl, jk -> li', Vloc, Gloc_dlr.density())
    # symmetrize
    Sig_DC_exchange = 0.5 * (Sig_DC_exchange + Sig_DC_exchange.conj().T)

    if verbose:
        print('static exchange part of DC')
        print(Sig_DC_exchange.real)
        if np.any(np.imag(Sig_DC_exchange) > 1e-3):
            print('Im:')
            print(np.imag(Sig_DC_exchange))
    return Sig_DC_dlr, Sig_DC_hartree, Sig_DC_exchange


def calc_W_from_Gloc(Gloc_dlr, U):
    r"""
    Calculate Wijkl from given constant U tensor and Gf on DLRMesh
    triqs notation for Uijkl:

    phi*_i(r) phi*_j(r') U(r,r') phi_l'(r') phi_k(r) = Uijkl c^+_i c^+_j' c_l' c_k

    where the ' denotes a spin index different from the other without '

    the according diagram is (left and right have same spin)::

       j (phi)         k' (phi)
         \              /
          <            <
           \__________/
           /          \
          >            >
         /              \
       i (phi*)          l'

    we now have to move to a product basis form to combine two indices
    i.e. go from nb,nb,nb,nb to nb**2,nb**2 tensors::

        Uji,kl = phi*_i(r) phi_j(r) U(r,r') phi*_k(r') phi_l(r')
               = Psi*_ji(r) U(r,r') Psi_kl(r')

    So we have to transform the triqs notation of Uijkl -> Uki,jl, i.e.
    swap col/rows as (2,0,1,3) to go to the basis and the in the end
    swap W_ki,jl back in reverse.

    Then we compute pubble polarizability as

    Pi_ab,kl(tau) = -2 G_bl(tau) G_ka(beta - tau)

    So that::

        [ U Pi(iwn) ]_ji,kl = sum_ab U_ji,ab Pi_ab,kl(iwn)

    i.e.::

       j'              a ___
         \              /   \ k
          <            <     \
           \__________/       \
           /          \       /
          >            >     /
         /              \___/ l
       i'               b

    then the screened Coulomb interaction in product basis is::

        W_ji,kl(iwn) = [1 - U Pi(iwn) ]^-1_ji,kl Uji,kl - Uji,kl

    (subtract static shift here), and finally convert back to triqs notation.


    Parameters
    ----------
    Gloc_dlr : BlockGf or Gf with MeshDLR
        local Green's function

    U : np.ndarray of with shape [Gloc_dlr.target_shape]*4 or dict of np.ndarray
        constant U tensor

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
        symmetrize=True
    )

    PI_dlr_t = Gf(mesh=mesh_bos, target_shape=[nb] * 4)
    i_tau = 0
    for tau in Gloc_dlr_t.mesh:
        PI_dlr_t.data[i_tau,:,:] = -2 * np.einsum('bl, ka -> abkl', Gloc_dlr_t[tau], Gloc_dlr(Gloc_dlr_t.mesh.beta - tau))
        i_tau += 1

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


def convert_gw_output(job_h5, gw_h5, dlr_wmax=None, dlr_eps=None,
                      it_1e=0, it_2e=0,
                      delta_calc_type="tail_fit", delta_bath_fit=False,
                      ha_ev_conv = False):
    """
    read bdft output and convert to triqs Gf DLR objects

    Parameters
    ----------
    job_h5: string
        path to solid_dmft job file
    gw_h5: string
        path to GW checkpoint file for AIMBES code
    dlr_wmax: float
        wmax for dlr mesh, defaults to the wmax from the IR basis
    dlr_eps: float
        precision for dlr mesh, defaults to the precision from the IR basis
    it_1e: int, optional
        iteration to read from gw_h5 calculation for 1e downfolding, defaults to last iteration
    it_2e: int, optional
        iteration to read from gw_h5 calculation for 2e downfolding, defaults to last iteration
    ha_ev_conv: bool, optional
        convert energies from Hartree to eV, defaults to False

    Returns
    -------
    gw_data: dict
        dictionary holding all read objects: mu_emb, beta, lam, w_max, prec, mesh_dlr_iw_b,
        mesh_dlr_iw_f, n_orb, G0_dlr, Gloc_dlr, Sigma_imp_dlr, Sigma_imp_DC_dlr, Uloc_dlr,
        Vloc, Hloc0, Vhf_dc, Vhf
    ir_kernel: sparse_ir kernel object
        IR kernel with AIMBES paramaters
    """

    mpi.report('reading output from aimbes code')

    gw_data = {}

    if ha_ev_conv:
        conv_fac = HARTREE_EV
    else:
        conv_fac = 1.0

    with HDFArchive(gw_h5, 'r') as ar:
        if not it_1e or not it_2e:
            it_1e = ar['downfold_1e/final_iter']
            it_2e = ar['downfold_2e/final_iter']

        mpi.report(f'Reading results from downfold_1e iter {it_1e} and downfold_2e iter {it_2e} from given AIMBES chkpt file.')

        # auxilary quantities
        gw_data['it_1e'] = it_1e
        gw_data['it_2e'] = it_2e
        gw_data['mu_emb'] = ar[f'downfold_1e/iter{it_1e}']['mu']
        gw_data['beta'] = ar['imaginary_fourier_transform']['beta']
        gw_data['lam'] = ar['imaginary_fourier_transform']['lambda']
        gw_data['gw_wmax'] = gw_data['lam'] / gw_data['beta']
        gw_data['gw_dlr_wmax'] = gw_data['gw_wmax'] if dlr_wmax is None else dlr_wmax
        gw_data['number_of_spins'] = ar['system/number_of_spins']
        assert gw_data['number_of_spins'] == 1, 'spin calculations not yet supported in converter'

        prec = ar['imaginary_fourier_transform']['prec']
        if prec == 'high':
            # set to highest DLR precision possible
            gw_data['gw_ir_prec'] = 1e-15
            gw_data['gw_dlr_prec'] = 1e-13 if dlr_eps is None else dlr_eps
        elif prec == 'mid':
            gw_data['gw_ir_prec'] = 1e-10
            gw_data['gw_dlr_prec'] = 1e-10 if dlr_eps is None else dlr_eps
        elif prec == 'low':
            gw_data['gw_ir_prec'] = 1e-6
            gw_data['gw_dlr_prec'] = 1e-6 if dlr_eps is None else dlr_eps

        # 1 particle properties
        g_weiss_wsIab = ar[f'downfold_1e/iter{it_1e}']['g_weiss_wsIab']
        delta_wsIab = ar[f'downfold_1e/iter{it_1e}']['delta_wsIab']
        Sigma_dc_wsIab = ar[f'downfold_1e/iter{it_1e}']['Sigma_dc_wsIab']
        Gloc = ar[f'downfold_1e/iter{it_1e}']['Gloc_wsIab']
        gw_data['n_inequiv_shells'] = Gloc.shape[2]

        # 2 particle properties
        # TODO: discuss how the site index is used right now in bDFT
        Vloc_jk = ar[f'downfold_2e/iter{it_2e}']['Vloc_abcd']
        Uloc_ir_jk = ar[f'downfold_2e/iter{it_2e}']['Uloc_wabcd'][:, ...]
        # switch inner two indices to match triqs notation
        Vloc = np.zeros(Vloc_jk.shape, dtype=complex)
        Uloc_ir = np.zeros(Uloc_ir_jk.shape, dtype=complex)
        n_orb = Vloc.shape[0]
        for or1, or2, or3, or4 in itertools.product(range(n_orb), repeat=4):
            Vloc[or1, or2, or3, or4] = Vloc_jk[or1, or3, or2, or4]
            for ir_w in range(Uloc_ir_jk.shape[0]):
                Uloc_ir[ir_w, or1, or2, or3, or4] = Uloc_ir_jk[ir_w, or1, or3, or2, or4]

        Vhf_dc_sIab = ar[f'downfold_1e/iter{it_1e}']['Vhf_dc_sIab'][0, 0]
        Vhf_sIab = ar[f'downfold_1e/iter{it_1e}']['Vhf_gw_sIab'][0, 0]

        if 'Vcorr_gw_sIab' in ar[f'downfold_1e/iter{it_1e}']:
            mpi.report('Found Vcorr_sIab in the bdft checkpoint file, '
                       'i.e. Embedding on top of an effective QP Hamiltonian.')
            qp_emb = True
        else:
            Sigma_wsIab = ar[f'downfold_1e/iter{it_1e}']['Sigma_gw_wsIab']
            qp_emb = False
        mpi.report("")

    # get IR object
    mpi.report('Creating IR kernel and convert to DLR.')
    # create IR kernel
    mpi.report("\nReading IR representation from aimbes code...")
    ir_kernel = IAFT(beta=gw_data['beta'], lmbda=gw_data['lam'], prec=gw_data['gw_ir_prec'])

    mpi.report("Constructing DLR mesh (wmax, eps) = ({}, {})...".format(gw_data['gw_dlr_wmax'], gw_data['gw_dlr_prec']))
    gw_data['mesh_dlr_iw_b'] = MeshDLRImFreq(
        beta=gw_data['beta'] / conv_fac,
        statistic='Boson',
        w_max=gw_data['gw_dlr_wmax'] * conv_fac,
        eps=gw_data['gw_dlr_prec'],
        symmetrize=True
    )
    gw_data['mesh_dlr_iw_f'] = MeshDLRImFreq(
        beta=gw_data['beta'] / conv_fac,
        statistic='Fermion',
        w_max=gw_data['gw_dlr_wmax'] * conv_fac,
        eps=gw_data['gw_dlr_prec'],
        symmetrize=True
    )

    if delta_calc_type not in {"analytic", "tail_fit"}:
        raise ValueError("calc_type must be either \'analytic\' or \'tail_fit\'.")

    if delta_calc_type == "analytic":
        Hloc0, delta_wsIab = read_t_and_delta(gw_h5, it_1e)
        ir_imp_kernel = ir_kernel
    elif delta_calc_type == "tail_fit":
        Hloc0, delta_wsIab, ir_imp_kernel = tail_fit_g_weiss(g_weiss_wsIab, ir_kernel, gw_data,
                                                             wmax_imp=dlr_wmax, eps_imp=dlr_eps)
    if delta_bath_fit:
        bath_fitting(delta_wsIab, ir_imp_kernel.wn_mesh('f')*np.pi/ir_imp_kernel.beta)
    Hloc0 = Hloc0[0,0]

    # need to update g_weiss?

    mpi.report("")

    (
        U_dlr_list,
        G0_dlr_list,
        delta_dlr_list,
        Gloc_dlr_list,
        Sigma_dlr_list,
        Sigma_DC_dlr_list,
        V_list,
        Hloc_list,
        Vhf_list,
        Vhf_dc_list,
        n_orb_list,
    ) = [], [], [], [], [], [], [], [], [], [], []
    for ish in range(gw_data['n_inequiv_shells']):
        # fit IR Uloc on DLR iw mesh
        temp = _get_dlr_from_IR(Uloc_ir*conv_fac, ir_kernel, gw_data['mesh_dlr_iw_b'], dim=4)
        Uloc_dlr = BlockGf(name_list=['up', 'down'], block_list=[temp, temp], make_copies=True)

        U_dlr_list.append(Uloc_dlr)
        V_list.append({'up': Vloc.copy()*conv_fac, 'down': Vloc*conv_fac})
        Hloc_list.append({'up': Hloc0.copy()*conv_fac, 'down': Hloc0*conv_fac})
        Vhf_list.append({'up': Vhf_sIab.copy()*conv_fac, 'down': Vhf_sIab*conv_fac})
        Vhf_dc_list.append({'up': Vhf_dc_sIab.copy()*conv_fac, 'down': Vhf_dc_sIab*conv_fac})
        n_orb_list.append(n_orb)

        temp = _get_dlr_from_IR(g_weiss_wsIab[:, 0, ish, :, :]/conv_fac, ir_kernel, gw_data['mesh_dlr_iw_f'], dim=2)
        G0_dlr = BlockGf(name_list=['up', 'down'], block_list=[temp, temp], make_copies=True)
        G0_dlr_list.append(G0_dlr)

        # FIXME make consistent usage of ir_kernel and ir_imp_kernel
        temp = _get_dlr_from_IR(delta_wsIab[:, 0, ish, :, :]/conv_fac, ir_imp_kernel, gw_data['mesh_dlr_iw_f'], dim=2)
        delta_dlr = BlockGf(name_list=['up', 'down'], block_list=[temp, temp], make_copies=True)
        delta_dlr_list.append(delta_dlr)

        temp = _get_dlr_from_IR(Gloc[:, 0, ish, :, :]/conv_fac, ir_kernel, gw_data['mesh_dlr_iw_f'], dim=2)
        Gloc_dlr = BlockGf(name_list=['up', 'down'], block_list=[temp, temp], make_copies=True)
        Gloc_dlr_list.append(Gloc_dlr)

        # since Sigma can have a static shift we return DLR Imfreq mesh
        if not qp_emb:
            temp = _get_dlr_from_IR(Sigma_wsIab[:, 0, ish, :, :]*conv_fac, ir_kernel, gw_data['mesh_dlr_iw_f'], dim=2)
            Sigma_dlr = BlockGf(name_list=['up', 'down'], block_list=[temp, temp], make_copies=True)
            Sigma_dlr_list.append(Sigma_dlr)

        temp = _get_dlr_from_IR(Sigma_dc_wsIab[:, 0, ish, :, :]*conv_fac, ir_kernel, gw_data['mesh_dlr_iw_f'], dim=2)
        Sigma_DC_dlr = BlockGf(name_list=['up', 'down'], block_list=[temp, temp], make_copies=True)
        Sigma_DC_dlr_list.append(Sigma_DC_dlr)

    gw_data['G0_dlr'] = G0_dlr_list
    gw_data['delta_dlr'] = delta_dlr_list
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
    mpi.report(f'Writing results in {job_h5}/DMFT_input')

    with HDFArchive(job_h5, 'a') as ar:
        if 'DMFT_input' not in ar:
            ar.create_group('DMFT_input')
        if f'iter{it_1e}' not in ar['DMFT_input']:
            ar['DMFT_input'].create_group(f'iter{it_1e}')

        for key, value in gw_data.items():
            ar[f'DMFT_input/iter{it_1e}'][key] = value

    mpi.report(f'finished writing results in {job_h5}/DMFT_input')
    return gw_data, ir_kernel


