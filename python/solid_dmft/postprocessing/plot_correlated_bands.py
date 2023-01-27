#!/usr/bin/env python3
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
Reads DMFT_ouput observables such as real-frequency Sigma and a Wannier90
TB Hamiltonian to compute spectral properties. It runs in two modes,
either calculating the bandstructure or Fermi slice.
Written by Sophie Beck, 2021-2022
TODO:
- extend to multi impurity systems
- make proper use of rot_mat from DFT_Tools (atm it assumed that wannier_hr and Sigma are written in the same basis)
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import Normalize
from matplotlib import cm
from scipy.optimize import brentq
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema
import numpy as np
import itertools
import skimage.measure

from h5 import HDFArchive
from triqs.gf import BlockGf, MeshReFreq, Gf
from triqs.lattice.utils import TB_from_wannier90, k_space_path
from triqs_dft_tools.sumk_dft import SumkDFT


def _linefit(x, y, interval, spacing=50, addspace=0.0):

    def calc_Z(slope): return 1/(1-slope)

    x = np.array(x)
    lim_l, lim_r = interval
    indices = np.where(np.logical_and(x >= lim_l, x <= lim_r))
    fit = np.polyfit(x[indices], y[indices], 1)
    slope = fit[1]
    Z = calc_Z(slope)
    f_x = np.poly1d(fit)
    x_cont = np.linspace(x[indices][0] - addspace, x[indices][-1] + addspace, spacing)

    return x_cont, f_x(x_cont), fit


def lambda_matrix_w90_t2g(add_lambda):

    lambda_x, lambda_y, lambda_z = add_lambda

    lambda_matrix = np.zeros((6, 6), dtype=complex)
    lambda_matrix[0, 1] = -1j*lambda_z/2.0
    lambda_matrix[0, 5] = 1j*lambda_x/2.0
    lambda_matrix[1, 5] = -lambda_y/2.0
    lambda_matrix[2, 3] = -1j*lambda_x/2.0
    lambda_matrix[2, 4] = lambda_y/2.0
    lambda_matrix[3, 4] = 1j*lambda_z/2.0
    lambda_matrix += np.transpose(np.conjugate(lambda_matrix))

    return lambda_matrix


def change_basis(n_orb, orbital_order_to, orbital_order_from):

    change_of_basis = np.eye(n_orb)
    for ct, orb in enumerate(orbital_order_to):
        orb_idx = orbital_order_from.index(orb)
        change_of_basis[orb_idx, :] = np.roll(np.eye(n_orb, 1), ct)[:, 0]

    return change_of_basis


def print_matrix(matrix, n_orb, text):

    print('{}:'.format(text))

    if np.any(matrix.imag > 1e-4):
        fmt = '{:16.4f}' * n_orb
    else:
        fmt = '{:8.4f}' * n_orb
        matrix = matrix.real

    for row in matrix:
        print((' '*4 + fmt).format(*row))


def _sigma_from_dmft(n_orb, orbital_order, with_sigma, spin, orbital_order_dmft=None, eta=0.0, **specs):

    if orbital_order_dmft is None:
        orbital_order_dmft = orbital_order

    if with_sigma == 'calc':
        print('Setting Sigma from {}'.format(specs['dmft_path']))

        sigma_imp_list = []
        dc_imp_list = []
        with HDFArchive(specs['dmft_path'], 'r') as ar:
            for icrsh in range(ar['dft_input']['n_inequiv_shells']):
                try:
                    sigma = ar['DMFT_results'][specs['it']][f'Sigma_freq_{icrsh}']
                    assert isinstance(sigma.mesh, MeshReFreq), 'Imported Greens function must be real frequency'
                except(KeyError, AssertionError):
                    try:
                        sigma = ar['DMFT_results'][specs['it']][f'Sigma_maxent_{icrsh}']
                    except KeyError:
                        try:
                            sigma = ar['DMFT_results'][specs['it']][f'Sigma_Refreq_{icrsh}']
                        except KeyError:
                            raise KeyError('Provide either "Sigma_freq_0" in real frequency, "Sigma_Refreq_0" or "Sigma_maxent_0".')

                sigma_imp_list.append(sigma)

            # DC is stored for each correlated shell
            for ish in range(ar['dft_input']['n_corr_shells']):
                dc_imp_list.append(ar['DMFT_results'][specs['it']]['DC_pot'][ish])

            mu_dmft = ar['DMFT_results'][specs['it']]['chemical_potential_post']

            sum_k = SumkDFT(specs['dmft_path'])
            sum_k.block_structure = ar['DMFT_input/block_structure']
            sum_k.deg_shells = ar['DMFT_input/deg_shells']
            sum_k.put_Sigma(sigma_imp_list)
            sum_k.set_mu = mu_dmft
            sum_k.dc_imp = dc_imp_list

            # use add_dc function to rotate to sumk block structure and subtract the DC
            sigma_sumk = sum_k.add_dc(iw_or_w="w")

            assert np.allclose(sum_k.proj_mat[0], sum_k.proj_mat[-1]), 'upfolding works only when proj_mat is the same for all kpoints (wannier mode)'

            # from IPython import embed; embed()
            # now upfold with proj_mat to band basis, this only works for the
            # case where proj_mat is equal for all k points (wannier mode)
            sigma = Gf(mesh=sigma.mesh, target_shape=[n_orb, n_orb])
            for ish in range(ar['dft_input']['n_corr_shells']):
                sigma += sum_k.upfold(ik=0, ish=ish,
                                      bname=spin, gf_to_upfold=sigma_sumk[ish][spin],
                                      gf_inp=sigma)

        # already subtracted
        dc = 0.0

    else:
        print('Setting Sigma from memory')

        sigma = with_sigma[spin]
        dc = specs['dc'][0][spin][0, 0]
        mu_dmft = specs['mu_dmft']

    SOC = (spin == 'ud')
    w_mesh_dmft = np.linspace(sigma.mesh.omega_min, sigma.mesh.omega_max, len(sigma.mesh))
    assert sigma.target_shape[0] == n_orb, f'Number of Wannier orbitals: {n_orb} and self-energy target_shape {sigma.target_shape} does not match'

    sigma_mat = sigma.data.real - np.eye(n_orb) * dc + 1j * sigma.data.imag

    # rotate sigma from orbital_order_dmft to orbital_order
    change_of_basis = change_basis(n_orb, orbital_order, orbital_order_dmft)
    sigma_mat = np.einsum('ij, kjl -> kil', np.linalg.inv(change_of_basis), np.einsum('ijk, kl -> ijl', sigma_mat, change_of_basis))

    # set up mesh
    if 'w_mesh' in specs:
        freq_dict = specs['w_mesh']
        w_mesh = np.linspace(*freq_dict['window'], freq_dict['n_w'])
        freq_dict.update({'w_mesh': w_mesh})
    else:
        w_mesh = w_mesh_dmft
        freq_dict = {'w_mesh': w_mesh_dmft, 'n_w': len(sigma.mesh), 'window': [sigma.mesh.omega_min, sigma.mesh.omega_max]}

    sigma_interpolated = np.zeros((n_orb, n_orb, freq_dict['n_w']), dtype=complex)

    if specs['linearize']:
        print('Linearizing Sigma at zero frequency:')
        eta = eta * 1j
        iw0 = np.where(np.sign(w_mesh_dmft) is True)[0][0]-1
        if SOC:
            sigma_interpolated += np.expand_dims(sigma_mat[iw0, :, :], axis=-1)
        # linearize diagonal elements of sigma
        for ct in range(n_orb):
            _, _, fit_params = _linefit(w_mesh_dmft, sigma_mat[:, ct, ct], specs['linearize']['window'])
            zeroth_order, first_order = fit_params[::-1].real
            print('Zeroth and first order fit parameters: [{0:.4f}, {1:.4f}]'.format(zeroth_order, first_order))
            sigma_interpolated[ct, ct] = zeroth_order + freq_dict['w_mesh'] * first_order

    else:
        # eta is added on the level of the spectral function!
        eta = 0 * 1j
        # interpolate sigma
        def interpolate_sigma(w_mesh, w_mesh_dmft, orb1, orb2): return np.interp(w_mesh, w_mesh_dmft, sigma_mat[:, orb1, orb2])

        for ct1, ct2 in itertools.product(range(n_orb), range(n_orb)):
            if ct1 != ct2 and not SOC:
                continue
            sigma_interpolated[ct1, ct2] = interpolate_sigma(w_mesh, w_mesh_dmft, ct1, ct2)

    return sigma_interpolated, mu_dmft, freq_dict


def _sigma_from_model(n_orb, orbital_order, zeroth_order, first_order, eta=0.0, **w):

    print('Setting model Sigma')

    eta = eta * 1j

    # set up mesh
    freq_dict = w['w_mesh']
    w_mesh = np.linspace(*freq_dict['window'], freq_dict['n_w'])
    freq_dict.update({'w_mesh': w_mesh})

    # interpolate sigma
    sigma_interpolated = np.zeros((n_orb, n_orb, freq_dict['n_w']), dtype=complex)
    def approximate_sigma(zeroth_order, first_order, orb): return zeroth_order[orb] + freq_dict['w_mesh'] * first_order[orb]
    for ct, orb in enumerate(orbital_order):
        sigma_interpolated[ct, ct] = approximate_sigma(zeroth_order, first_order, ct)

    return sigma_interpolated, freq_dict


def _calc_alatt(n_orb, mu, eta, e_mat, sigma, qp_bands=False, e_vecs=None,
                proj_nuk=None, trace=True, **freq_dict):
    '''
    calculate slice of lattice spectral function for given TB dispersion / e_mat and self-energy
    Parameters
    ----------
    n_orb : int
          number of Wannier orbitals
    proj_nuk : optinal, 2D numpy array (n_orb, n_k)
          projections to be applied on A(k,w) in band basis. Only works when band_basis=True
    Returns
    -------
    alatt_k_w : numpy array, either (n_k, n_w) or if trace=False (n_k, n_w, n_orb)
            Lattice Green's function on specified k-path / mesh
    '''

    # adjust to system size
    def upscale(quantity, n_orb): return quantity * np.identity(n_orb)
    mu = upscale(mu, n_orb)
    eta = upscale(eta, n_orb)
    if isinstance(e_vecs, np.ndarray):
        sigma_rot = np.zeros(sigma.shape, dtype=complex)

    w_vec = np.array([upscale(freq_dict['w_mesh'][w], n_orb) for w in range(freq_dict['n_w'])])
    n_k = e_mat.shape[2]

    if not qp_bands:
        if trace:
            alatt_k_w = np.zeros((n_k, freq_dict['n_w']))
        else:
            alatt_k_w = np.zeros((n_k, freq_dict['n_w'], n_orb))

        def invert_and_trace(w, eta, mu, e_mat, sigma, trace, proj=None):
            # inversion is automatically vectorized over first axis of 3D array (omega first index now)
            Glatt = np.linalg.inv(w + eta[None, ...] + mu[None, ...] - e_mat[None, ...] - sigma.transpose(2, 0, 1))
            A_w_nu = -1.0/np.pi * np.diagonal(Glatt, axis1=1, axis2=2).imag
            if isinstance(proj, np.ndarray):
                A_w_nu = A_w_nu * proj[None, :]
            if trace:
                return np.sum(A_w_nu, axis=1)
            else:
                return A_w_nu

        for ik in range(n_k):
            # if evecs are given transform sigma into band basis
            if isinstance(e_vecs, np.ndarray):
                sigma_rot = np.einsum('ij,jkw->ikw', e_vecs[:, :, ik].conjugate().transpose(), np.einsum('ijw,jk->ikw', sigma, e_vecs[:, :, ik]))
                if isinstance(proj_nuk, np.ndarray):
                    alatt_k_w[ik, :] = invert_and_trace(w_vec, eta, mu, e_mat[:, :, ik], sigma_rot, trace, proj_nuk[:, ik])
                else:
                    alatt_k_w[ik, :] = invert_and_trace(w_vec, eta, mu, e_mat[:, :, ik], sigma_rot, trace)
            else:
                alatt_k_w[ik, :] = invert_and_trace(w_vec, eta, mu, e_mat[:, :, ik], sigma, trace)

    else:
        alatt_k_w = np.zeros((n_k, n_orb))
        kslice = np.zeros((freq_dict['n_w'], n_orb))
        def kslice_interp(orb): return interp1d(freq_dict['w_mesh'], kslice[:, orb])

        for ik in range(n_k):
            for iw, w in enumerate(freq_dict['w_mesh']):
                np.fill_diagonal(sigma[:, :, iw], np.diag(sigma[:, :, iw]).real)
                #sigma[:,:,iw] = sigma[:,:,iw].real
                kslice[iw], _ = np.linalg.eigh(upscale(w, n_orb) + eta + mu - e_mat[:, :, ik] - sigma[:, :, iw])

            for orb in range(n_orb):
                w_min, w_max = freq_dict['window']
                try:
                    x0 = brentq(kslice_interp(orb), w_min, w_max)
                    w_bin = int((x0 - w_min) / ((w_max - w_min) / freq_dict['n_w']))
                    alatt_k_w[ik, orb] = freq_dict['w_mesh'][w_bin]
                except ValueError:
                    pass

    return alatt_k_w


def _calc_kslice(n_orb, mu, eta, e_mat, sigma, qp_bands, e_vecs=None, proj_nuk=None, **freq_dict):
    '''
    calculate lattice spectral function for given TB dispersion / e_mat and self-energy
    Parameters
    ----------
    n_orb : int
          number of Wannier orbitals
    proj_nuk : optinal, 2D numpy array (n_orb, n_k)
          projections to be applied on A(k,w) in band basis. Only works when band_basis=True
    Returns
    -------
    alatt_k_w : numpy array, either (n_k, n_w) or if trace=False (n_k, n_w, n_orb)
            Lattice Green's function on specified k-path / mesh
    '''

    # adjust to system size
    def upscale(quantity, n_orb): return quantity * np.identity(n_orb)
    mu = upscale(mu, n_orb)
    eta = upscale(eta, n_orb)

    iw0 = np.where(np.sign(freq_dict['w_mesh']) == True)[0][0]-1
    print_matrix(sigma[:, :, iw0], n_orb, 'Zero-frequency Sigma')

    if isinstance(e_vecs, np.ndarray):
        sigma_rot = np.zeros(sigma.shape, dtype=complex)

    n_kx, n_ky = e_mat.shape[2:4]

    if not qp_bands:
        alatt_k_w = np.zeros((n_kx, n_ky))

        def invert_and_trace(w, eta, mu, e_mat, sigma, proj=None):
            # inversion is automatically vectorized over first axis of 3D array (omega first index now)
            Glatt = np.linalg.inv(w + eta + mu - e_mat - sigma)
            A_nu = -1.0/np.pi * np.diagonal(Glatt).imag
            if isinstance(proj, np.ndarray):
                A_nu = A_nu * proj
            return np.sum(A_nu)

        for ikx, iky in itertools.product(range(n_kx), range(n_ky)):
            if isinstance(e_vecs, np.ndarray):
                sigma_rot = np.einsum('ij,jk->ik',
                                      e_vecs[:, :, ikx, iky].conjugate().transpose(),
                                      np.einsum('ij,jk->ik', sigma[:, :, iw0], e_vecs[:, :, ikx, iky]))
            else:
                sigma_rot = sigma[:, :, iw0]

            if isinstance(proj_nuk, np.ndarray):
                alatt_k_w[ikx, iky] = invert_and_trace(upscale(freq_dict['w_mesh'][iw0], n_orb), eta, mu,
                                                       e_mat[:, :, ikx, iky], sigma_rot, proj_nuk[:, ikx, iky])
            else:
                alatt_k_w[ikx, iky] = invert_and_trace(upscale(freq_dict['w_mesh'][iw0], n_orb), eta, mu, e_mat[:, :, ikx, iky], sigma_rot)

    else:
        assert n_kx == n_ky, 'Not implemented for N_kx != N_ky'

        def search_for_extrema(data):
            # return None for no extrema, [] if ends of interval are the only extrema,
            # list of indices if local extrema are present
            answer = np.all(data > 0) or np.all(data < 0)
            if answer:
                return
            else:
                roots = []
                roots.append(list(argrelextrema(data, np.greater)[0]))
                roots.append(list(argrelextrema(data, np.less)[0]))
                roots = sorted([item for sublist in roots for item in sublist])
            return roots

        alatt_k_w = np.zeros((n_kx, n_ky, n_orb))
        # go through grid horizontally, then vertically
        for it in range(2):
            kslice = np.zeros((n_kx, n_ky, n_orb))

            for ik1 in range(n_kx):
                e_temp = e_mat[:, :, :, ik1] if it == 0 else e_mat[:, :, ik1, :]
                for ik2 in range(n_kx):
                    e_val, _ = np.linalg.eigh(eta + mu - e_temp[:, :, ik2] - sigma[:, :, iw0])
                    k1, k2 = [ik2, ik1] if it == 0 else [ik1, ik2]
                    kslice[k1, k2] = e_val

                for orb in range(n_orb):
                    temp_kslice = kslice[:,ik1,orb] if it == 0 else kslice[ik1,:,orb]
                    roots = search_for_extrema(temp_kslice)
                    # iterate through sections between extrema
                    if roots is not None:
                        idx_1 = 0
                        for root_ct in range(len(roots) + 1):
                            idx_2 = roots[root_ct] if root_ct < len(roots) else n_kx
                            root_section = temp_kslice[idx_1:idx_2+1]
                            try:
                                x0 = brentq(interp1d(np.linspace(idx_1, idx_2, len(root_section)), root_section), idx_1, idx_2)
                                k1, k2 = [int(np.floor(x0)), ik1] if it == 0 else [ik1, int(np.floor(x0))]
                                alatt_k_w[k1, k2, orb] += 1
                            except(ValueError):
                                pass
                            idx_1 = idx_2

        alatt_k_w[np.where(alatt_k_w > 1)] = 1

    return alatt_k_w


def _get_tb_bands(e_mat, proj_on_orb=[None], **specs):
    '''
    calculate eigenvalues and eigenvectors for given list of e_mat on kmesh
    Parameters
    ----------
    e_mat : numpy array of shape (n_orb, n_orb, nk) or (n_orb, n_orb, nk, nk)
    Returns
    -------
    e_val : numpy array of shape (n_orb, n_orb, nk) or (n_orb, n_orb, nk, nk)
        eigenvalues as matrix
    e_vec : numpy array of shape (n_orb, n_orb, nk) or (n_orb, n_orb, nk, nk)
        eigenvectors as matrix
    '''

    e_val = np.zeros((e_mat.shape), dtype=complex)
    e_vec = np.zeros((e_mat.shape), dtype=complex)
    n_orb = e_mat.shape[0]

    for ikx in range(e_mat.shape[2]):
        # if we have a 2d kmesh e_mat is dim=4
        if len(e_mat.shape) == 4:
            for iky in range(e_mat.shape[3]):
                e_val[range(n_orb), range(n_orb), ikx, iky], e_vec[:, :, ikx, iky] = np.linalg.eigh(e_mat[:, :, ikx, iky])
        else:
            e_val[range(n_orb), range(n_orb), ikx], e_vec[:, :, ikx] = np.linalg.eigh(e_mat[:, :, ikx])

    if proj_on_orb[0] is not None:
        print(f'calculating projection on orbitals {proj_on_orb}')
        total_proj = np.zeros(np.shape(e_vec[0]))
        for band in range(n_orb):
            for orb in proj_on_orb:
                total_proj[band] += np.real(e_vec[orb, band] * e_vec[orb, band].conjugate())
    else:
        total_proj = None

    return e_val, e_vec, total_proj


def get_tb_kslice(tb, mu_tb, **specs):

    w90_paths = list(map(lambda section: (np.array(specs[section[0]]), np.array(specs[section[1]])), specs['bands_path']))
    upper_left = np.diff(w90_paths[0][::-1], axis=0)[0]
    lower_right = np.diff(w90_paths[1], axis=0)[0]
    Z = np.array(specs['Z'])

    FS_kx_ky, band_char = get_kx_ky_FS(lower_right, upper_left, Z, tb, N_kxy=specs['n_k'], kz=specs['kz'], fermi=mu_tb)

    return FS_kx_ky, band_char


def _fract_ind_to_val(x, ind):
    ind[ind == len(x)-1] = len(x)-1-1e-6
    int_ind = [int(indi) for indi in ind]
    int_ind_p1 = [int(indi)+1 for indi in ind]
    return x[int_ind] + (x[int_ind_p1] - x[int_ind])*(np.array(ind)-np.array(int_ind))


def get_kx_ky_FS(lower_right, upper_left, Z, tb, select=None, N_kxy=10, kz=0.0, fermi=0.0):

    assert np.abs(fermi) < 1e-2, 'finite value of Fermi level not implemented. Subtract Fermi level from local Hamiltonian.'

    # create mesh
    kx = np.linspace(0, 0.5, N_kxy)
    ky = np.linspace(0, 0.5, N_kxy)

    if select is None:
        select = np.array(range(tb.NOrbitalsInUnitCell))

    # go in horizontal arrays from bottom to top
    E_FS = np.zeros((tb.NOrbitalsInUnitCell, N_kxy, N_kxy))
    for kyi in range(N_kxy):
        path_FS = [(upper_left/(N_kxy-1)*kyi + kz*Z, lower_right+upper_left/(N_kxy-1)*kyi+kz*Z)]
        k_vec, _ = k_space_path(path_FS, num=N_kxy)
        E_FS[:, :, kyi] = tb.dispersion(k_vec).transpose()

    contours = {}
    FS_kx_ky = {}
    FS_kx_ky_prim = {}
    band_char = {}
    # contour for each sheet
    for sheet in range(tb.NOrbitalsInUnitCell):
        contours[sheet] = skimage.measure.find_contours(E_FS[sheet, :, :], fermi)

    sheet_ct = 0
    for sheet in contours.keys():
        for sec_per_sheet in range(np.shape(contours[sheet])[0]):
            # once on 2D cubic mesh
            FS_kx_ky[sheet_ct] = np.vstack([_fract_ind_to_val(kx, contours[sheet][sec_per_sheet][:, 0]),
                                            _fract_ind_to_val(ky, contours[sheet][sec_per_sheet][:, 1]),
                                            kz*np.ones(len(contours[sheet][sec_per_sheet][:, 0]))]).T.reshape(-1, 3)
            # repeat on actual mesh for computing the weights
            ks_skimage = contours[sheet][sec_per_sheet]/(N_kxy-1)
            FS_kx_ky_prim[sheet_ct] = (+ np.einsum('i,j->ij', ks_skimage[:, 0], lower_right)
                                       + np.einsum('i,j->ij', ks_skimage[:, 1], upper_left)
                                       + np.einsum('i,j->ij', kz * np.ones(ks_skimage.shape[0]), Z))
            band_char[sheet_ct] = {}
            # compute the weight aka band character
            for ct_k, k_on_sheet in enumerate(FS_kx_ky_prim[sheet_ct]):
                E_mat = tb.fourier(k_on_sheet)
                e_val, e_vec = np.linalg.eigh(E_mat[select[:, np.newaxis], select])
                orb_on_FS = np.argmin(np.abs(e_val))

                band_char[sheet_ct][ct_k] = [np.round(np.real(e_vec[orb, orb_on_FS]*np.conjugate(e_vec[orb, orb_on_FS])), 4) for orb in range(len(select))]
            sheet_ct += 1

    return FS_kx_ky, band_char


def _setup_plot_bands(ax, special_k, k_points_labels, freq_dict):

    ax.axhline(y=0, c='gray', ls='--', lw=0.8, zorder=0)
    ax.set_ylabel(r'$\omega - \mu$ (eV)')
#     ax.set_ylim(*freq_dict['window'])
    for ik in special_k:
        ax.axvline(x=ik, linewidth=0.7, color='k', zorder=0.5)
    ax.set_xticks(special_k)
    ax.set_xlim(special_k[0], special_k[-1])
    k_points_labels = [r'$\Gamma$' if k == 'G' else k for k in k_points_labels]
    ax.set_xticklabels(k_points_labels)


def setup_plot_kslice(ax):

    ax.set_aspect(1)
    # ax.set_xlim(0,1)
    # ax.set_ylim(0,1)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel(r'$k_x\pi/a$')
    ax.set_ylabel(r'$k_y\pi/b$')


def plot_bands(fig, ax, alatt_k_w, tb_data, freq_dict, n_orb, tb=True, alatt=False, qp_bands=False, **plot_dict):

    proj_on_orb = tb_data['proj_on_orb']
    total_proj = tb_data['proj_nuk']

    if alatt:
        if alatt_k_w is None:
            raise ValueError('A(k,w) unknown. Specify "with_sigma = True"')
        if qp_bands:
            for orb in range(n_orb):
                ax.scatter(tb_data['k_mesh'], alatt_k_w[:, orb].T, c=np.array([eval('cm.'+plot_dict['colorscheme_qpbands'])(1.0)]), zorder=2., s=1.)
        else:
            kw_x, kw_y = np.meshgrid(tb_data['k_mesh'], freq_dict['w_mesh'])

            vmax = plot_dict['vmax'] if 'vmax' in plot_dict else np.max(alatt_k_w)
            vmin = plot_dict['vmin'] if 'vmin' in plot_dict else 0.0

            graph = ax.pcolormesh(kw_x, kw_y, alatt_k_w.T, cmap=plot_dict['colorscheme_alatt'],
                                  norm=Normalize(vmin=vmin, vmax=vmax), shading='gouraud')
            colorbar = plt.colorbar(graph)
            colorbar.set_label(r'$A(k, \omega)$')

    if tb:
        # if projection is requested, _get_tb_bands() ran already
        if proj_on_orb[0] is not None:
            eps_nuk = tb_data['e_mat']
            evec_nuk = tb_data['e_vecs']
        else:
            eps_nuk, evec_nuk, _ = _get_tb_bands(**tb_data)
        for band in range(n_orb):
            if not proj_on_orb[0] is not None:
                color = eval('cm.'+plot_dict['colorscheme_bands'])(1.0)
                ax.plot(tb_data['k_mesh'], eps_nuk[band, band].real - tb_data['mu_tb'], c=color, label=r'tight-binding', zorder=1., lw=1)
            else:
                color = eval('cm.'+plot_dict['colorscheme_bands'])(total_proj[band])
                ax.scatter(tb_data['k_mesh'], eps_nuk[band, band].real - tb_data['mu_tb'], c=color, s=1, label=r'tight-binding', zorder=1.)

    _setup_plot_bands(ax, tb_data['special_k'], tb_data['k_points_labels'], freq_dict)


def plot_kslice(fig, ax, alatt_k_w, tb_data, freq_dict, n_orb, tb_dict, tb=True, alatt=False, quarter=0, **plot_dict):

    proj_on_orb = tb_data['proj_on_orb']
    if quarter:
        assert isinstance(quarter, int) or all(isinstance(x, int) for x in quarter), 'quarter should be'\
            f'an integer or list of integers, but is {type(quarter)}.'

    if isinstance(quarter, int):
        quarter = [quarter]

    sign = [1, -1]
    quarters = np.array([sign, sign])
    four_quarters = list(itertools.product(*quarters))
    used_quarters = [four_quarters[x] for x in quarter]

    vmax = plot_dict['vmax'] if 'vmax' in plot_dict else np.max(alatt_k_w)
    vmin = plot_dict['vmin'] if 'vmin' in plot_dict else 0.0

    if alatt:
        if alatt_k_w is None:
            raise ValueError('A(k,w) unknown. Specify "with_sigma = True"')
        n_kx, n_ky = tb_data['e_mat'].shape[2:4]
        kx, ky = np.meshgrid(range(n_kx), range(n_ky))
        for (qx, qy) in used_quarters:
            if len(alatt_k_w.shape) > 2:
                for orb in range(n_orb):
                    ax.contour(qx * kx/(n_kx-1), qy * ky/(n_ky-1), alatt_k_w[:, :, orb].T,
                               colors=np.array([eval('cm.'+plot_dict['colorscheme_qpbands'])(0.7)]), levels=1, zorder=2)
            else:
                graph = ax.pcolormesh(qx * kx/(n_kx-1), qy * ky/(n_ky-1), alatt_k_w.T,
                                      cmap=plot_dict['colorscheme_kslice'],
                                      norm=Normalize(vmin=vmin, vmax=vmax),
                                      shading='gouraud')
                #colorbar = plt.colorbar(graph)
                #colorbar.set_label(r'$A(k, 0$)')

    if tb:
        FS_kx_ky, band_char = get_tb_kslice(tb_data['tb'], tb_data['mu_tb'], **tb_dict)
        for sheet in FS_kx_ky.keys():
            for k_on_sheet in range(FS_kx_ky[sheet].shape[0]):
                if not proj_on_orb[0] is not None:
                    if isinstance(plot_dict['colorscheme_kslice'], str):
                        color = eval('cm.'+plot_dict['colorscheme_kslice'])(1.0)
                    else:
                        color = plot_dict['colorscheme_kslice']
                else:
                    total_proj = 0
                    for orb in proj_on_orb:
                        total_proj += band_char[sheet][k_on_sheet][orb]
                    color = eval('cm.'+plot_dict['colorscheme_kslice'])(total_proj)
                for (qx, qy) in used_quarters:
                    ax.plot(2*qx * FS_kx_ky[sheet][k_on_sheet:k_on_sheet+2, 0], 2*qy * FS_kx_ky[sheet][k_on_sheet:k_on_sheet+2, 1], '-',
                            solid_capstyle='round', c=color, zorder=1., label=plot_dict['label'] if 'label' in plot_dict else '')

    setup_plot_kslice(ax)

    return ax


def get_dmft_bands(n_orb, w90_path, w90_seed, mu_tb, add_spin=False, add_lambda=None, add_local=None,
                   with_sigma=None, fermi_slice=False, qp_bands=False, orbital_order_to=None,
                   add_mu_tb=False, band_basis=False, proj_on_orb=None, trace=True, eta=0.0,
                   mu_shift=0.0, proj_nuk=None, **specs):
    '''
    Extract tight-binding from given w90 seed_hr.dat and seed.wout files, and then extract from
    given solid_dmft calculation the self-energy and construct the spectral function A(k,w) on
    given k-path.
    Parameters
    ----------
    n_orb : int
        Number of Wannier orbitals in seed_hr.dat
    w90_path : string
        Path to w90 files
    w90_seed : string
        Seed of wannier90 calculation, i.e. seed_hr.dat and seed.wout
    add_spin : bool, default=False
        Extend w90 Hamiltonian by spin indices
    add_lambda : float, default=None
        Add SOC term with strength add_lambda (works only for t2g shells)
    add_local : numpy array, default=None
        Add local term of dimension (n_orb x n_orb)
    with_sigma : str, or BlockGf, default=None
        Add self-energy to spectral function? Can be either directly take
        a triqs BlockGf object or can be either 'calc' or 'model'
        'calc' reads results from h5 archive (solid_dmft)
        in case 'calc' or 'model' are specified a extra kwargs dict has
        to be given sigma_dict containing information about the self-energy
    add_mu_tb : bool, default=False
        Add the TB specified chemical potential to the lattice Green function
        set to True if DMFT calculation was performed with DFT fermi subtracted.
    proj_on_orb : int or list of int, default=None
        orbital projections to be made for the spectral function and TB bands
        the integer refer to the orbitals read
    trace : bool, default=True
        Return trace over orbitals for spectral function. For special
        post-processing purposes this can be set to False giving the returned
        alatt_k_w an extra dimension n_orb
    eta : float, default=0.0
        Broadening of spectral function, finitie shift on imaginary axis
        if with_sigma=None it has to be provided !=0.0
    mu_shift : float, default=0.0
        Manual extra shift when calculating the spectral function
    proj_nuk : numpy array, default [None]
        Extra projections to be applied to the final spectral function
        per orbital and k-point. Has to match shape of final lattice Green
        function. Will be applied together with proj_on_orb if specified.
    Returns
    -------
    tb_data : dict
       tight binding dict containing the kpoint mesh, dispersion / emat, and eigenvectors
    alatt_k_w : numpy array (float) of dim n_k x n_w ( x n_orb if trace=False)
        lattice spectral function data on the kpoint mesh defined in tb_data and frequency
        mesh defined in freq_dict
    freq_dict : dict
        frequency mesh information on which alatt_k_w is evaluated
    '''

    # set default ordering
    if 'orbital_order_w90' in specs:
        orbital_order_w90 = specs['orbital_order_w90']
    else:
        orbital_order_w90 = list(range(n_orb))

    if orbital_order_to is None:
        orbital_order_to = orbital_order_w90

    # checks
    assert len(set(orbital_order_to)) == len(orbital_order_to), 'Please provide a unique identifier for each orbital.'

    assert set(orbital_order_w90) == set(orbital_order_to), f'Identifiers of orbital_order_to and orbital_order_w90'\
        f'do not match! orbital_order_to is {orbital_order_to}, but orbital_order_w90 is {orbital_order_w90}.'

    assert with_sigma or eta != 0.0, 'if no Sigma is provided eta has to be different from 0.0'

    # proj_on_orb
    assert isinstance(proj_on_orb, (int, type(None))) or all(isinstance(x, (int, type(None))) for x in proj_on_orb), 'proj_on_orb should be '\
        f'an integer or list of integers, but is {type(specs["proj_on_orb"])}.'

    if isinstance(proj_on_orb, (int, type(None))):
        proj_on_orb = [proj_on_orb]
    else:
        proj_on_orb = proj_on_orb

    # if projection is requested we have to use band_basis
    if proj_on_orb[0] is not None:
        band_basis = True

    # if proj_nuk is given we need to use the band_basis
    if isinstance(proj_nuk, np.ndarray) and not band_basis:
        band_basis = True

    # set up Wannier Hamiltonian
    n_orb_rescale = 2 * n_orb if add_spin else n_orb
    change_of_basis = change_basis(n_orb, orbital_order_to, orbital_order_w90)
    H_add_loc = np.zeros((n_orb_rescale, n_orb_rescale), dtype=complex)
    if not isinstance(add_local, type(None)):
        assert np.shape(add_local) == (n_orb_rescale, n_orb_rescale), 'add_local must have dimension (n_orb, n_orb), but has '\
                f'dimension {np.shape(add_local)}'
        H_add_loc += add_local
    if add_spin and add_lambda:
        H_add_loc += lambda_matrix_w90_t2g(add_lambda)
    eta = eta * 1j
    n_k = specs['n_k']

    tb = TB_from_wannier90(path=w90_path, seed=w90_seed, extend_to_spin=add_spin, add_local=H_add_loc)
    # print local H(R)
    h_of_r = tb.hoppings[(0, 0, 0)][2:5, 2:5] if add_spin else tb.hoppings[(0, 0, 0)]
    h_of_r = np.einsum('ij, jk -> ik', np.linalg.inv(change_of_basis), np.einsum('ij, jk -> ik', h_of_r, change_of_basis))
    if n_orb <= 12:
        print_matrix(h_of_r, n_orb, 'H(R=0)')

    # bands info
    w90_paths = list(map(lambda section: (np.array(specs[section[0]]), np.array(specs[section[1]])), specs['bands_path']))
    k_points_labels = [k[0] for k in specs['bands_path']] + [specs['bands_path'][-1][1]]

    # calculate tight-binding eigenvalues
    if not fermi_slice:
        k_vec, k_1d = k_space_path(w90_paths, bz=tb.bz, num=n_k)
        special_k = np.append(k_1d[0::n_k], k_1d[-1::])
        e_mat = tb.fourier(k_vec).transpose(1, 2, 0)
        if add_spin:
            e_mat = e_mat[2:5, 2:5]
        e_mat = np.einsum('ij, jkl -> ikl', np.linalg.inv(change_of_basis), np.einsum('ijk, jm -> imk', e_mat, change_of_basis))
    else:
        assert 'Z' in specs, 'Please provide Z point coordinate in tb_data_dict as input coordinate'
        Z = np.array(specs['Z'])

        k_vec = np.zeros((n_k*n_k, 3))
        e_mat = np.zeros((n_orb_rescale, n_orb_rescale, n_k, n_k), dtype=complex)

        upper_left = np.diff(w90_paths[0][::-1], axis=0)[0]
        lower_right = np.diff(w90_paths[1], axis=0)[0]
        for ik_y in range(n_k):
            path_along_x = [(upper_left/(n_k-1)*ik_y + specs['kz']*Z, lower_right+upper_left/(n_k-1)*ik_y+specs['kz']*Z)]
            k_vec[ik_y*n_k:ik_y*n_k+n_k, :], k_1d = k_space_path(path_along_x, bz=tb.bz, num=n_k)
            special_k = np.append(k_1d[0::n_k], k_1d[-1::])
            e_mat[:, :, :, ik_y] = tb.fourier(k_vec[ik_y*n_k:ik_y*n_k+n_k, :]).transpose(1, 2, 0)
        if add_spin:
            e_mat = e_mat[2:5, 2:5]
        e_mat = np.einsum('ij, jklm -> iklm', np.linalg.inv(change_of_basis), np.einsum('ijkl, jm -> imkl', e_mat, change_of_basis))

    if band_basis:
        e_mat, e_vecs, orb_proj = _get_tb_bands(e_mat, proj_on_orb)
    else:
        e_vecs = total_proj = orb_proj = None

    # now we merge proj_nuk and orb_proj (has reverse shape)
    if isinstance(proj_nuk, np.ndarray) and isinstance(orb_proj, np.ndarray):
        proj_nuk = proj_nuk * orb_proj
    elif not isinstance(proj_nuk, np.ndarray) and isinstance(orb_proj, np.ndarray):
        proj_nuk = orb_proj

    # dmft output
    if with_sigma:
        sigma_types = ['calc', 'model']
        if isinstance(with_sigma, str):
            if with_sigma not in sigma_types:
                raise ValueError('Invalid sigma type. Expected one of: {}'.format(sigma_types))
        elif not isinstance(with_sigma, BlockGf):
            raise ValueError('Invalid sigma type. Expected BlockGf.')

        # get sigma
        if with_sigma == 'model':
            delta_sigma, freq_dict = _sigma_from_model(n_orb, orbital_order_to, **specs)
            mu = mu_tb + mu_shift
        # else is from dmft or memory:
        else:
            delta_sigma, mu_dmft, freq_dict = _sigma_from_dmft(n_orb, orbital_order_to, with_sigma, **specs)
            mu = mu_dmft + mu_shift

        if add_mu_tb:
            print('Adding mu_tb to DMFT μ; assuming DMFT was run with subtracted dft μ.')
            mu += mu_tb

        print('μ={:2.4f} eV set for calculating A(k,ω)'.format(mu))

        assert n_orb == delta_sigma.shape[0] and n_orb == delta_sigma.shape[
            1], f'Number of orbitals n_orb={n_orb} and shape of sigma: {delta_sigma.shape} does not match'
        if isinstance(proj_nuk, np.ndarray):
            assert n_orb == proj_nuk.shape[0], f'Number of orbitals n_orb={n_orb} does not match shape of proj_nuk: {proj_nuk.shape[0]}'
            if not fermi_slice:
                assert proj_nuk.shape[-1] == e_vecs.shape[
                    2], f'Number of kpoints in proj_nuk : {proj_nuk.shape[-1]} does not match number of kpoints in e_vecs: {e_vecs.shape[2]}'
            else:
                assert proj_nuk.shape == tuple([n_orb, e_vecs.shape[2], e_vecs.shape[3]]
                                               ), f'shape of projectors {proj_nuk.shape} does not match expected shape of [{n_orb},{e_vecs.shape[2]},{e_vecs.shape[3]}]'

        # calculate alatt
        if not fermi_slice:
            alatt_k_w = _calc_alatt(n_orb, mu, eta, e_mat, delta_sigma, qp_bands, e_vecs=e_vecs,
                                    trace=trace, proj_nuk=proj_nuk, **freq_dict)
        else:
            alatt_k_w = _calc_kslice(n_orb, mu, eta, e_mat, delta_sigma, qp_bands, e_vecs=e_vecs,
                                     proj_nuk=proj_nuk, **freq_dict)
        freq_dict['sigma'] = delta_sigma
    else:
        freq_dict = {}
        freq_dict['w_mesh'] = None
        freq_dict['window'] = None
        alatt_k_w = None

    tb_data = {'k_mesh': k_1d, 'special_k': special_k, 'k_points': k_vec,
               'k_points_labels': k_points_labels, 'e_mat': e_mat,
               'e_vecs': e_vecs, 'tb': tb, 'mu_tb': mu_tb,
               'proj_on_orb': proj_on_orb, 'proj_nuk': proj_nuk}

    return tb_data, alatt_k_w, freq_dict
