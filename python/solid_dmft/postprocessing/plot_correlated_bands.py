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
import numpy as np
import itertools
import skimage.measure

from h5 import HDFArchive
from triqs.gf import BlockGf, MeshReFreq
from triqs.lattice.utils import TB_from_wannier90, k_space_path


def _linefit(x, y, interval, spacing=50, addspace=0.0):

    calc_Z = lambda slope: 1/(1-slope)

    x = np.array(x)
    lim_l, lim_r = interval
    indices = np.where(np.logical_and(x>=lim_l, x<=lim_r))
    fit = np.polyfit(x[indices], y[indices], 1)
    slope = fit[1]
    Z = calc_Z(slope)
    f_x = np.poly1d(fit)
    x_cont = np.linspace(x[indices][0] - addspace, x[indices][-1] + addspace, spacing)

    return x_cont, f_x(x_cont), fit

def lambda_matrix_w90_t2g(add_lambda):

    lambda_x, lambda_y, lambda_z = add_lambda

    lambda_matrix = np.zeros((6,6), dtype=complex)
    lambda_matrix[0,1] = -1j*lambda_z/2.0
    lambda_matrix[0,5] =  1j*lambda_x/2.0
    lambda_matrix[1,5] =    -lambda_y/2.0
    lambda_matrix[2,3] = -1j*lambda_x/2.0
    lambda_matrix[2,4] =     lambda_y/2.0
    lambda_matrix[3,4] =  1j*lambda_z/2.0
    lambda_matrix += np.transpose(np.conjugate(lambda_matrix))

    return lambda_matrix

def change_basis(n_orb, orbital_order_to, orbital_order_from):

    change_of_basis = np.eye(n_orb)
    for ct, orb in enumerate(orbital_order_to):
        orb_idx = orbital_order_from.index(orb)
        change_of_basis[orb_idx,:] = np.roll(np.eye(n_orb,1),ct)[:,0]

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

def _sigma_from_dmft(n_orb, orbital_order, with_sigma, spin, block, orbital_order_dmft, eta=0.0, **specs):

    if with_sigma == 'calc':
        print('Setting Sigma from {}'.format(specs['dmft_path']))

        with HDFArchive(specs['dmft_path'],'r') as ar:
            try:
                sigma = ar['DMFT_results'][specs['it']]['Sigma_freq_0']
                assert isinstance(sigma.mesh, MeshReFreq), 'Imported Greens function must be real frequency'
            except(KeyError, AssertionError):
                try:
                    sigma = ar['DMFT_results'][specs['it']]['Sigma_maxent_0']
                except KeyError:
                    raise KeyError('Provide either "Sigma_freq_0" in real frequency or "Sigma_maxent_0".')
            dc = ar['DMFT_results'][specs['it']]['DC_pot'][0][spin][0,0]
            mu = ar['DMFT_results'][specs['it']]['chemical_potential_post']
            dft_mu = ar['DMFT_results/observables']['mu'][0]

    else:
        print('Setting Sigma from memory')

        sigma = with_sigma
        dc = specs['dc'][0][spin][0,0]
        mu = specs['dmft_mu']
        dft_mu = 0

    block_spin = spin + '_' + str(block) if with_sigma == 'calc' else spin
    SOC = (spin == 'ud')
    w_mesh_dmft = [x.real for x in sigma[block_spin].mesh]
    sigma_mat = {block_spin: sigma[block_spin].data.real - np.eye(n_orb) * dc + 1j * sigma[block_spin].data.imag}

    # rotate sigma from orbital_order_dmft to orbital_order
    change_of_basis = change_basis(n_orb, orbital_order, orbital_order_dmft)
    sigma_mat[block_spin] = np.einsum('ij, kjl -> kil', np.linalg.inv(change_of_basis), np.einsum('ijk, kl -> ijl', sigma_mat[block_spin], change_of_basis))

    # set up mesh
    freq_dict = specs['w_mesh']
    w_mesh = np.linspace(*freq_dict['window'], freq_dict['n_w'])
    freq_dict.update({'w_mesh': w_mesh})

    sigma_interpolated = np.zeros((n_orb, n_orb, freq_dict['n_w']), dtype=complex)

    if specs['linearize']:
        print('Linearizing Sigma at zero frequency:')
        eta = eta * 1j
        iw0 = np.where(np.sign(w_mesh_dmft) == True)[0][0]-1
        if SOC: sigma_interpolated += np.expand_dims(sigma_mat[block_spin][iw0,:,:], axis=-1)
        # linearize diagonal elements of sigma
        for ct in range(n_orb):
            _, _, fit_params = _linefit(w_mesh_dmft, sigma_mat[block_spin][:,ct,ct], specs['linearize']['window'])
            zeroth_order, first_order = fit_params[::-1].real
            print('Zeroth and first order fit parameters: [{0:.4f}, {1:.4f}]'.format(zeroth_order,first_order))
            sigma_interpolated[ct,ct] = zeroth_order + freq_dict['w_mesh'] * first_order

    else:
        # eta is added on the level of the spectral function!
        eta = 0 * 1j
        # interpolate sigma
        interpolate_sigma = lambda w_mesh, w_mesh_dmft, orb1, orb2: np.interp(w_mesh, w_mesh_dmft, sigma_mat[block_spin][:, orb1, orb2])

        for ct1, ct2 in itertools.product(range(n_orb), range(n_orb)):
            if ct1 != ct2 and not SOC: continue
            sigma_interpolated[ct1,ct2] = interpolate_sigma(w_mesh, w_mesh_dmft, ct1, ct2)

    return sigma_interpolated, mu, dft_mu, freq_dict

def _sigma_from_model(n_orb, orbital_order, zeroth_order, first_order, efermi, eta=0.0, **w):

    print('Setting model Sigma')

    mu = dft_mu = efermi
    eta = eta * 1j

    # set up mesh
    freq_dict = w['w_mesh']
    w_mesh = np.linspace(*freq_dict['window'], freq_dict['n_w'])
    freq_dict.update({'w_mesh': w_mesh})

    # interpolate sigma
    sigma_interpolated = np.zeros((n_orb, n_orb, freq_dict['n_w']), dtype=complex)
    approximate_sigma = lambda zeroth_order, first_order, orb: zeroth_order[orb] + freq_dict['w_mesh'] * first_order[orb]
    for ct, orb in enumerate(orbital_order):
        sigma_interpolated[ct,ct] = approximate_sigma(zeroth_order, first_order, ct)

    return sigma_interpolated, mu, dft_mu, freq_dict

def _calc_alatt(n_orb, mu, eta, e_mat, sigma, qp_bands=False, e_vecs=np.array([None]) ,trace=True,**freq_dict):

    # adjust to system size
    upscale = lambda quantity, n_orb: quantity * np.identity(n_orb)
    mu = upscale(mu, n_orb)
    eta = upscale(eta, n_orb)
    if not e_vecs.any() == None:
        sigma_rot = np.zeros(sigma.shape, dtype=complex)

    w_vec = np.array([upscale(freq_dict['w_mesh'][w], n_orb) for w in range(freq_dict['n_w'])])
    n_k = e_mat.shape[2]

    if not qp_bands:
        if trace:
            alatt_k_w = np.zeros((n_k, freq_dict['n_w']))
        else:
            alatt_k_w = np.zeros((n_k, freq_dict['n_w'], n_orb))
        def invert_and_trace(w, eta, mu, e_mat, sigma, trace):
            # inversion is automatically vectorized over first axis of 3D array (omega first index now)
            Glatt =  np.linalg.inv(w + eta[None,...] + mu[None,...] - e_mat[None,...] - sigma.transpose(2,0,1) )
            if trace:
                return -1.0/np.pi * np.trace( Glatt ,axis1=1, axis2=2).imag
            else:
                return -1.0/np.pi * np.diagonal( Glatt ,axis1=1, axis2=2).imag

        for ik in range(n_k):
            # if evecs are given transform sigma into band basis
            if not e_vecs.any() == None:
                sigma_rot = np.einsum('ij,jkw->ikw', e_vecs[:,:,ik].conjugate().transpose(), np.einsum('ijw,jk->ikw', sigma, e_vecs[:,:,ik]))
                alatt_k_w[ik,:] = invert_and_trace(w_vec, eta, mu, e_mat[:,:,ik], sigma_rot, trace)
            else:
                alatt_k_w[ik,:] = invert_and_trace(w_vec, eta, mu, e_mat[:,:,ik], sigma, trace)

    else:
        alatt_k_w = np.zeros((n_k, n_orb))
        kslice = np.zeros((freq_dict['n_w'], n_orb))
        kslice_interp = lambda orb: interp1d(freq_dict['w_mesh'], kslice[:, orb])

        for ik in range(n_k):
            for iw, w in enumerate(freq_dict['w_mesh']):
                np.fill_diagonal(sigma[:,:,iw], np.diag(sigma[:,:,iw]).real)
                #sigma[:,:,iw] = sigma[:,:,iw].real
                kslice[iw], _ = np.linalg.eigh( upscale(w, n_orb) + eta + mu - e_mat[:,:,ik] - sigma[:,:,iw])

            for orb in range(n_orb):
                w_min, w_max = freq_dict['window']
                try:
                    x0 = brentq( kslice_interp(orb), w_min, w_max)
                    w_bin = int( (x0 - w_min) / ((w_max - w_min)/ freq_dict['n_w']) )
                    alatt_k_w[ik, orb] = freq_dict['w_mesh'][w_bin]
                except ValueError:
                    pass

    return alatt_k_w

def _calc_kslice(n_orb, mu, eta, e_mat, sigma, qp_bands, **freq_dict):

    # adjust to system size
    upscale = lambda quantity, n_orb: quantity * np.identity(n_orb)
    mu = upscale(mu, n_orb)
    eta = upscale(eta, n_orb)

    iw0 = np.where(np.sign(freq_dict['w_mesh']) == True)[0][0]-1
    print_matrix(sigma[:,:,iw0], n_orb, 'Zero-frequency Sigma')

    n_kx, n_ky = e_mat.shape[2:4]

    if not qp_bands:
        alatt_k_w = np.zeros((n_kx, n_ky))
        invert_and_trace = lambda w, eta, mu, e_mat, sigma: -1.0/np.pi * np.trace( np.linalg.inv( w + eta + mu - e_mat - sigma ).imag )

        for ikx, iky in itertools.product(range(n_kx), range(n_ky)):
            alatt_k_w[ikx, iky] = invert_and_trace(upscale(freq_dict['w_mesh'][iw0], n_orb), eta, mu, e_mat[:,:,ikx,iky], sigma[:,:,iw0])
    else:
        assert n_kx == n_ky, 'Not implemented for N_kx != N_ky'
        alatt_k_w = np.zeros((n_kx, n_ky, n_orb))
        for it in range(2):
            kslice = np.zeros((n_kx, n_ky, n_orb))
            if it == 0:
                kslice_interp = lambda ik, orb: interp1d(range(n_kx), kslice[:, ik, orb])
            else:
                kslice_interp = lambda ik, orb: interp1d(range(n_kx), kslice[ik, :, orb])

            for ik1 in range(n_kx):
                e_temp = e_mat[:,:,:,ik1] if it == 0 else e_mat[:,:,ik1,:]
                for ik2 in range(n_kx):
                    e_val, _ = np.linalg.eigh( eta + mu - e_temp[:,:,ik2] - sigma[:,:,iw0])
                    k1, k2 = [ik2, ik1] if it == 0 else [ik1, ik2]
                    kslice[k1, k2] = e_val

                for orb in range(n_orb):
                    try:
                        x0 = brentq( kslice_interp(ik1, orb), 0, n_kx - 1)
                        k1, k2 = [int(np.floor(x0)), ik1] if it == 0 else [ik1, int(np.floor(x0))]
                        alatt_k_w[k1, k2, orb] += 1
                    except ValueError:
                        pass

        alatt_k_w[np.where(alatt_k_w > 1)] = 1

    return alatt_k_w

def _get_tb_bands(k_mesh, e_mat, **specs):

    e_val = np.zeros((e_mat.shape[0], k_mesh.shape[0]), dtype=complex)
    e_vec = np.zeros(np.shape(e_mat), dtype=complex)
    for ik in range(np.shape(e_mat)[2]):
        e_val[:,ik], e_vec[:,:,ik] = np.linalg.eigh(e_mat[:,:,ik])

    return e_val, e_vec

def get_tb_kslice(tb, **specs):

    w90_paths = list(map(lambda section: (np.array(specs[section[0]]), np.array(specs[section[1]])), specs['bands_path']))
    upper_left = np.diff(w90_paths[0][::-1], axis=0)[0]
    lower_right = np.diff(w90_paths[1], axis=0)[0]
    Z = np.array(specs['Z'])

    FS_kx_ky, band_char = get_kx_ky_FS(lower_right, upper_left, Z, tb, N_kxy=specs['n_k'], kz=specs['kz'], fermi=0.0)

    return FS_kx_ky, band_char

def _fract_ind_to_val(x,ind):
    ind[ind == len(x)-1] = len(x)-1-1e-6
    int_ind = [int(indi) for indi in ind]
    int_ind_p1 = [int(indi)+1 for indi in ind]
    return x[int_ind] + (x[int_ind_p1] - x[int_ind])*(np.array(ind)-np.array(int_ind))

def get_kx_ky_FS(lower_right, upper_left, Z, tb, select=None, N_kxy=10, kz=0.0, fermi=0.0):

    assert np.abs(fermi) < 1e-2, 'finite value of Fermi level not implemented. Subtract Fermi level from local Hamiltonian.'

    # create mesh
    kx = np.linspace(0,0.5,N_kxy)
    ky = np.linspace(0,0.5,N_kxy)

    if select is None: select = np.array(range(tb.NOrbitalsInUnitCell))

    # go in horizontal arrays from bottom to top
    E_FS = np.zeros((tb.NOrbitalsInUnitCell,N_kxy,N_kxy))
    for kyi in range(N_kxy):
        path_FS = [(upper_left/(N_kxy-1)*kyi +kz*Z, lower_right+upper_left/(N_kxy-1)*kyi+kz*Z)]
        k_vec, _ = k_space_path(path_FS, num=N_kxy)
        E_FS[:,:,kyi] = tb.dispersion(k_vec).transpose()

    contours = {}
    FS_kx_ky = {}
    FS_kx_ky_prim = {}
    band_char = {}
    # contour for each sheet
    for sheet in range(tb.NOrbitalsInUnitCell):
        contours[sheet] = skimage.measure.find_contours(E_FS[sheet,:,:], fermi)

    sheet_ct = 0
    for sheet in contours.keys():
        for sec_per_sheet in range(np.shape(contours[sheet])[0]):
            # once on 2D cubic mesh
            FS_kx_ky[sheet_ct] = np.vstack([_fract_ind_to_val(kx,contours[sheet][sec_per_sheet][:,0]),
                                            _fract_ind_to_val(ky,contours[sheet][sec_per_sheet][:,1]),
                                             kz*np.ones(len(contours[sheet][sec_per_sheet][:,0]))]).T.reshape(-1,3)
            # repeat on actual mesh for computing the weights
            ks_skimage = contours[sheet][sec_per_sheet]/(N_kxy-1)
            FS_kx_ky_prim[sheet_ct] = (+ np.einsum('i,j->ij', ks_skimage[:,0], lower_right)
                                       + np.einsum('i,j->ij', ks_skimage[:,1], upper_left)
                                       + np.einsum('i,j->ij', kz * np.ones(ks_skimage.shape[0]), Z))
            band_char[sheet_ct] = {}
            # compute the weight aka band character
            for ct_k, k_on_sheet in enumerate(FS_kx_ky_prim[sheet_ct]):
                E_mat = tb.fourier(k_on_sheet)
                e_val, e_vec = np.linalg.eigh(E_mat[select[:,np.newaxis],select])
                orb_on_FS = np.argmin(np.abs(e_val))

                band_char[sheet_ct][ct_k] = [np.round(np.real(e_vec[orb,orb_on_FS]*np.conjugate(e_vec[orb,orb_on_FS])),4) for orb in range(len(select))]
            sheet_ct += 1

    return FS_kx_ky, band_char

def _setup_plot_bands(ax, special_k, k_points_labels, freq_dict):

    ax.axhline(y=0,c='gray',ls='--',lw=0.8, zorder=0)
    ax.set_ylabel(r'$\omega - \mu$ (eV)')
#     ax.set_ylim(*freq_dict['window'])
    for ik in special_k:
        ax.axvline(x=ik, linewidth=0.7, color='k',zorder=0.5)
    ax.set_xticks(special_k)
    ax.set_xlim(special_k[0], special_k[-1])
    k_points_labels = [r'$\Gamma$' if k == 'G' else k for k in k_points_labels]
    ax.set_xticklabels(k_points_labels)

def setup_plot_kslice(ax):

    ax.set_aspect(1)
    #ax.set_xlim(0,1)
    #ax.set_ylim(0,1)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel(r'$k_x\pi/a$')
    ax.set_ylabel(r'$k_y\pi/b$')

def check_and_convert_plotting(quarter=None, **specs):

    # proj_on_orb
    assert isinstance(specs['proj_on_orb'], (int, type(None))) or all(isinstance(x, (int, type(None))) for x in specs['proj_on_orb']), 'proj_on_orb should be '\
            f'an integer or list of integers, but is {type(specs["proj_on_orb"])}.'

    if isinstance(specs['proj_on_orb'], (int, type(None))):
        proj_on_orb = [specs['proj_on_orb']]
    else:
        proj_on_orb = specs['proj_on_orb']

    # quarter
    if quarter:
        assert isinstance(quarter, int) or all(isinstance(x, int) for x in quarter), 'quarter should be'\
                f'an integer or list of integers, but is {type(quarter)}.'

    if isinstance(quarter, int):
        quarter = [quarter]

    return (proj_on_orb, quarter) if quarter else proj_on_orb

def plot_bands(fig, ax, alatt_k_w, tb_data, freq_dict, n_orb, dft_mu, tb=True, alatt=False, qp_bands=False, **plot_dict):

    proj_on_orb = check_and_convert_plotting(**plot_dict, quarter=None)

    if alatt:
        if alatt_k_w is None: raise ValueError('A(k,w) unknown. Specify "with_sigma = True"')
        if qp_bands:
            for orb in range(n_orb):
                ax.scatter(tb_data['k_mesh'], alatt_k_w[:,orb].T, c=np.array([eval('cm.'+plot_dict['colorscheme_qpbands'])(1.0)]), zorder=2., s=1.)
        else:
            kw_x, kw_y = np.meshgrid(tb_data['k_mesh'], freq_dict['w_mesh'])
            graph = ax.pcolormesh(kw_x, kw_y, alatt_k_w.T, cmap=plot_dict['colorscheme_bands'],
                                  norm=Normalize(vmin=plot_dict['vmin'], vmax=np.max(alatt_k_w)), shading='gouraud')
            colorbar = plt.colorbar(graph)
            colorbar.set_label(r'$A(k, \omega)$')

    if tb:
        eps_nuk, evec_nuk = _get_tb_bands(**tb_data)
        for band in range(n_orb):
            if not proj_on_orb[0] is not None:
                color = eval('cm.'+plot_dict['colorscheme_bands'])(1.0)
                ax.plot(tb_data['k_mesh'], eps_nuk[band].real - dft_mu, c=color, label=r'tight-binding', zorder=1.)
            else:
                total_proj = np.zeros(np.shape(evec_nuk[0, band]))
                for orb in proj_on_orb:
                    total_proj += np.real(evec_nuk[orb, band] * evec_nuk[orb, band].conjugate())
                color = eval('cm.'+plot_dict['colorscheme_bands'])(total_proj)
                ax.scatter(tb_data['k_mesh'], eps_nuk[band].real - dft_mu, c=color, s=1, label=r'tight-binding', zorder=1.)

    _setup_plot_bands(ax, tb_data['special_k'], tb_data['k_points_labels'], freq_dict)

def plot_kslice(fig, ax, alatt_k_w, tb_data, freq_dict, n_orb, tb_dict, tb=True, alatt=False, quarter=0, **plot_dict):

    proj_on_orb, quarter = check_and_convert_plotting(**plot_dict, quarter=quarter)

    sign = [1,-1]
    quarters = np.array([sign,sign])
    four_quarters = list(itertools.product(*quarters))
    used_quarters = [four_quarters[x] for x in quarter]

    if alatt:
        if alatt_k_w is None: raise ValueError('A(k,w) unknown. Specify "with_sigma = True"')
        n_kx, n_ky = tb_data['e_mat'].shape[2:4]
        kx, ky = np.meshgrid(range(n_kx), range(n_ky))
        for (qx, qy) in used_quarters:
            if len(alatt_k_w.shape) > 2:
                for orb in range(n_orb):
                    ax.contour(qx * kx/(n_kx-1), qy * ky/(n_ky-1), alatt_k_w[:,:,orb].T, colors=np.array([eval('cm.'+plot_dict['colorscheme_qpbands'])(0.7)]), levels=1, zorder=2)
            else:
                graph = ax.pcolormesh(qx * kx/(n_kx-1), qy * ky/(n_ky-1), alatt_k_w.T, cmap=plot_dict['colorscheme_kslice'],
                                      norm=Normalize(vmin=plot_dict['vmin'], vmax=np.max(alatt_k_w)))
                #colorbar = plt.colorbar(graph)
                #colorbar.set_label(r'$A(k, 0$)')

    if tb:
        FS_kx_ky, band_char = get_tb_kslice(tb_data['tb'], **tb_dict)
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
                    ax.plot(2*qx * FS_kx_ky[sheet][k_on_sheet:k_on_sheet+2,0], 2*qy * FS_kx_ky[sheet][k_on_sheet:k_on_sheet+2,1], '-',
                            solid_capstyle='round', c=color, zorder=1.)

    setup_plot_kslice(ax)

    return ax

def get_dmft_bands(n_orb, w90_path, w90_seed, add_spin, mu, add_lambda, with_sigma=False, fermi_slice=False, qp_bands=False, orbital_order_to=['dxz', 'dyz', 'dxy'], band_basis=False, trace=True , eta=0.0, **specs):

    # checks
    assert len(set(orbital_order_to)) == len(orbital_order_to), 'Please provide a unique identifier for each orbital.'
    assert set(specs['orbital_order_w90']) == set(orbital_order_to), f'Identifiers of orbital_order_to and orbital_order_w90'\
            f'do not match! orbital_order_to is {orbital_order_to}, but orbital_order_w90 is {specs["orbital_order_w90"]}.'


    # set up Wannier Hamiltonian
    n_orb_rescale = 2 * n_orb if add_spin else n_orb
    change_of_basis = change_basis(n_orb, orbital_order_to, specs['orbital_order_w90'])
    H_add_loc = np.zeros((n_orb_rescale, n_orb_rescale), dtype=complex)
    H_add_loc += np.diag([-mu]*n_orb_rescale)
    if add_spin: H_add_loc += lambda_matrix_w90_t2g(add_lambda)
    eta = eta * 1j

    tb = TB_from_wannier90(path=w90_path, seed=w90_seed, extend_to_spin=add_spin, add_local=H_add_loc)
    # print local H(R)
    h_of_r = tb.hoppings[(0,0,0)][2:5,2:5] if add_spin else tb.hoppings[(0,0,0)]
    h_of_r = np.einsum('ij, jk -> ik', np.linalg.inv(change_of_basis), np.einsum('ij, jk -> ik', h_of_r, change_of_basis))
    if n_orb <=12: print_matrix(h_of_r, n_orb, 'H(R=0)')

    # bands info
    w90_paths = list(map(lambda section: (np.array(specs[section[0]]), np.array(specs[section[1]])), specs['bands_path']))
    k_points_labels = [k[0] for k in specs['bands_path']] + [specs['bands_path'][-1][1]]

    # calculate tight-binding eigenvalues
    if not fermi_slice:
        k_vec, k_1d = k_space_path(w90_paths, bz=tb.bz, num=specs['n_k'])
        special_k = np.append(k_1d[0::specs['n_k']], k_1d[-1::])
        e_mat = tb.fourier(k_vec).transpose(1, 2, 0)
        if add_spin: e_mat = e_mat[2:5,2:5]
        e_mat = np.einsum('ij, jkl -> ikl', np.linalg.inv(change_of_basis), np.einsum('ijk, jm -> imk', e_mat, change_of_basis))
        if band_basis:
            e_vals, e_vecs = _get_tb_bands(k_vec, e_mat)
            for ik in range(np.shape(e_mat)[2]):
                e_mat[:,:,ik] = np.zeros(e_mat[:,:,ik].shape)
                np.fill_diagonal(e_mat[:,:,ik],e_vals[:,ik])
        else:
            e_vecs = np.array([None])
    else:
        e_mat = np.zeros((n_orb_rescale, n_orb_rescale, specs['n_k'], specs['n_k']), dtype=complex)
        upper_left = np.diff(w90_paths[0][::-1], axis=0)[0]
        lower_right = np.diff(w90_paths[1], axis=0)[0]
        Z = np.array(specs['Z'])
        for ik_y in range(specs['n_k']):
            path_along_x = [(upper_left/(specs['n_k']-1)*ik_y +specs['kz']*Z, lower_right+upper_left/(specs['n_k']-1)*ik_y+specs['kz']*Z)]
            k_vec, k_1d = k_space_path(path_along_x, bz=tb.bz, num=specs['n_k'])
            special_k = np.append(k_1d[0::specs['n_k']], k_1d[-1::])
            e_mat[:,:,:,ik_y] = tb.fourier(k_vec).transpose(1,2,0)
        if add_spin: e_mat = e_mat[2:5,2:5]
        e_mat = np.einsum('ij, jklm -> iklm', np.linalg.inv(change_of_basis), np.einsum('ijkl, jm -> imkl', e_mat, change_of_basis))

    # dmft output
    if with_sigma:
        sigma_types = ['calc', 'model']
        if isinstance(with_sigma, str):
            if with_sigma not in sigma_types: raise ValueError('Invalid sigma type. Expected one of: {}'.format(sigma_types))
        elif not isinstance(with_sigma, BlockGf):
            raise ValueError('Invalid sigma type. Expected BlockGf.')

        # get sigma
        if with_sigma == 'model': delta_sigma, mu, dft_mu, freq_dict = _sigma_from_model(n_orb, orbital_order_to, **specs)
        # else is from dmft or memory:
        else: delta_sigma, mu, dft_mu, freq_dict = _sigma_from_dmft(n_orb, orbital_order_to, with_sigma, **specs)
        corrected_mu = mu - dft_mu

        # calculate alatt
        if not fermi_slice:
            alatt_k_w = _calc_alatt(n_orb, corrected_mu, eta, e_mat, delta_sigma, qp_bands, e_vecs=e_vecs, trace=trace, **freq_dict)
        else:
            alatt_k_w = _calc_kslice(n_orb, corrected_mu, eta, e_mat, delta_sigma, qp_bands, **freq_dict)
    else:
        dft_mu = mu
        freq_dict = {}
        freq_dict['w_mesh'] = None
        freq_dict['window'] = None
        alatt_k_w = None


    return {'k_mesh': k_1d, 'special_k': special_k, 'k_points': k_vec, 'k_points_labels': k_points_labels, 'e_mat': e_mat, 'tb': tb}, alatt_k_w, freq_dict, dft_mu

