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

# system
import numpy as np

# triqs
import triqs.utility.mpi as mpi
from triqs.gf import MeshImFreq, MeshImTime, Gf, make_gf_from_fourier
from triqs.gf.descriptors import Fourier
from triqs.gf.tools import inverse

from . import legendre_filter


def _bgf_to_vec(bgf, deg_shells=None):
    """
    flattens a BlockGf to a 1d numpy array
    if deg_shells is given, only non-deg shells are given back
    """
    if deg_shells:
        deg_vecs = []
        for deg_block in deg_shells:
            deg_vecs.append(bgf[deg_block[0]])
        vec = np.hstack(gf.data.flatten() for gf in deg_vecs)
    else:
        vec = np.hstack(bgf[bl].data.flatten() for bl in bgf.indices)
    return vec


def _vec_to_bgf(vec, bgf, deg_shells=None):
    G = bgf.copy()
    if deg_shells:
        start = 0
        for deg_block in deg_shells:
            # only the first element was stored in numpy vec
            bl = deg_block[0]
            points = len([w.value for w in G[bl].mesh]) * np.prod(list(G[bl].target_shape))
            end = start+points
            G[bl].data[:,:,:] = vec[start:end].reshape(G[bl].data.shape)
            start = end
            for other_blocks in deg_block[1:]:
                G[other_blocks] << G[bl]
    else:
        start = 0
        for i, bl in enumerate(G.indices):
            points = len([w.value for w in G[bl].mesh]) * np.prod(list(G[bl].target_shape))
            end = start+points
            G[bl].data[:,:,:] = vec[start:end].reshape(G[bl].data.shape)
            start = end
    return G


def _broyden_update(it, broyler, general_params, deg_shells, G0_freq, G0_freq_previous):
    """
    calculates the broyden update of G0 using the algorithm presented in:
    doi.org/10.1103/PhysRevB.80.125125 (Rok Zitko) &
    doi.org/10.1103/PhysRevB.38.12807 (D.D. Johnson)

    optimize function to find roots: F(G0)=dmft_step(G0) - G0 = 0
    where dmft_step gives back the normal defined G0=inverse(Gloc^-1 - Sigma)

    Parameters
    ----------
    it : int
        iteration number
    broyler: dict
        dict containing the broyden variables
    general_paramters: general params dict

    Returns
    -------
    """
    alpha = general_params['g0_mix']
    # hard code w_0
    w_0 = 0.01

    G0_broyden_update = G0_freq.copy()

    # build here leg compression
    if type(G0_freq.mesh) is MeshImFreq and 'n_l' in general_params:
        G0 = legendre_filter.apply(make_gf_from_fourier(G0_freq), order=general_params['n_l'] )
        G0_prev = legendre_filter.apply(make_gf_from_fourier(G0_freq_previous), order=general_params['n_l'] )
    else:
        G0 = G0_freq.copy()
        G0_prev = G0_freq_previous.copy()

    # if this is the first time broyden update is called (it==2)
    # store init G0 first
    if it == 2:
        # V holds the G0 as 1d numpy vectors
        broyler['V'].append(_bgf_to_vec(G0_prev, deg_shells))
        # F = dmft_step_G0 - V[-1]
        broyler['F'].append(_bgf_to_vec(G0, deg_shells)-broyler['V'][-1])
        # for the second iteration the broyden vec corresponds to simple linear mixing with alpha
        broyden_vec = np.zeros(broyler['V'][0].shape[0],dtype='complex')

    # for all other iteration calc V_mp1 = V[-1] + alpha * F[-1] - broyden_vec
    # broyden vec is sum(n=1 to it-1) gamma_m_n U^n
    # gamma_m_n = sum(k=1 to it-1)  c_k^m beta_k_n^m
    # gamma_m_n is a number for each iteration and U^n is of dim G0 for each iteration
    else:
        # get size of broyden update vector
        dim = broyler['V'][0].shape[0]
        broyden_vec = np.zeros(dim,dtype='complex')

        # define m as it-1-1 (last minus for idx starting at 0)
        m = it-2

        # update broyden input params
        broyler['F'].append(_bgf_to_vec(G0, deg_shells)-broyler['V'][-1])
        norm_F = np.linalg.norm(broyler['F'][-1]-broyler['F'][-2])
        print('broyden norm: {:.3f}'.format(norm_F))
        # dF is normed differences between last two functions
        broyler['dF'].append( (broyler['F'][-1]-broyler['F'][-2])/norm_F )
        # dV the normed difference between the two last input G0
        broyler['dV'].append( (broyler['V'][-1] - broyler['V'][-2])/norm_F )

        # now build broyden correction vector

        #only keep last m iterations
        mem_iter = general_params['broy_max_it']-1
        start = 0
        mat_dim = m
        if m > mem_iter and not general_params['broy_max_it'] == -1:
            for n in range(0, m-mem_iter):
                start += 1
                mat_dim -= 1

        # first create beta_m
        A = np.zeros((mat_dim,mat_dim),dtype='complex')
        gamma = np.zeros(mat_dim,dtype='complex')
        for k_i in range(start,m):
            for n_i in range(start,m):
                A[k_i-start,n_i-start] = np.dot(broyler['dF'][n_i].conjugate().transpose(),broyler['dF'][k_i])
                if n_i == k_i:
                    A[k_i-start,n_i-start] += w_0**2
        beta_m = np.linalg.inv(A)
        for n_i in range(start,m):
            U_n = alpha * broyler['dF'][n_i] + broyler['dV'][n_i]

            # construct gamma
            for k_i in range(start,m):
                c_k = np.dot(broyler['dF'][k_i].conjugate().transpose(),broyler['F'][-1])
                gamma[n_i-start] += c_k*beta_m[k_i-start,n_i-start]
            broyden_vec += gamma[n_i-start] * U_n

    # update input G0
    V_mp1 = broyler['V'][-1] + alpha * broyler['F'][-1] - broyden_vec
    broyler['V'].append(V_mp1)

    # convert Gl back to G0_freq
    if type(G0_freq.mesh) is MeshImFreq and 'n_l' in general_params:
        G0_l_update = _vec_to_bgf(broyler['V'][-1], G0, deg_shells)
        for block, g in G0_l_update:
            g.enforce_discontinuity(np.identity(g.target_shape[0]))
            G0_broyden_update[block].set_from_legendre(g)

    else:
        G0_broyden_update << _vec_to_bgf(broyler['V'][-1], G0, deg_shells)

    return G0_broyden_update, broyler


def mix_g0(solver, general_params, icrsh, archive, G0_freq_previous, it, deg_shell):
    if general_params['g0_mix_type'] == 'linear' and general_params['g0_mix'] < 1.0:
        mpi.report('\nlinear mixing G0 with previous iteration by factor {:.3f}\n'.format(general_params['g0_mix']))
        solver.G0_freq << (general_params['g0_mix'] * solver.G0_freq
                                   + (1-general_params['g0_mix']) * G0_freq_previous)
    elif general_params['g0_mix_type'] == 'broyden':
        mpi.report('\nbroyden mixing G0 with previous iteration with alpha {:.3f}'.format(general_params['g0_mix']))
        mpi.report('\n############\n!!!! WARNING !!!! broyden mixing is still in early testing stage ! Use with caution.\n############\n')
        # TODO implement broyden mixing for Sigma
        if mpi.is_master_node():
            broyler = archive['DMFT_results']['broyler']
            # calculate the next G0 via broyden scheme
            G0_broyden_update, broyler[icrsh] = _broyden_update(it, broyler[icrsh], general_params,
                                                                deg_shell, solver.G0_freq,
                                                                G0_freq_previous)
            # store broyden update to h5 archive
            archive['DMFT_results']['broyler'] = broyler
            solver.G0_freq << G0_broyden_update

        solver.G0_freq << mpi.bcast(solver.G0_freq)

    return solver


def mix_sigma(general_params, n_inequiv_shells, solvers, Sigma_freq_previous):
    if general_params['sigma_mix'] < 1.0 and mpi.is_master_node():
        print('mixing sigma with previous iteration by factor {:.3f}\n'.format(general_params['sigma_mix']))
        for icrsh in range(n_inequiv_shells):
            solvers[icrsh].Sigma_freq << (general_params['sigma_mix'] * solvers[icrsh].Sigma_freq
                                        + (1-general_params['sigma_mix']) * Sigma_freq_previous[icrsh])

    for icrsh in range(n_inequiv_shells):
        solvers[icrsh].Sigma_freq << mpi.bcast(solvers[icrsh].Sigma_freq)

    return solvers
