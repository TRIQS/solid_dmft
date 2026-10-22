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
from h5 import HDFArchive


def _mixed_blocks(bgf, deg_shells):
    """
    groups of block names whose first entry is stored in the Broyden vector;
    the other entries of a group are copies of the first. Blocks that are in
    no degenerate shell form a group of their own.
    """
    groups = [list(deg_block) for deg_block in deg_shells] if deg_shells else []
    in_shell = {bl for group in groups for bl in group}
    groups += [[bl] for bl in bgf.indices if bl not in in_shell]
    return groups


def _bgf_to_vec(bgf, deg_shells=None):
    """
    flattens a BlockGf to a 1d numpy array
    if deg_shells is given, only the first block of each degenerate shell is stored
    """
    return np.hstack([bgf[group[0]].data.flatten() for group in _mixed_blocks(bgf, deg_shells)])


def _vec_to_bgf(vec, bgf, deg_shells=None):
    G = bgf.copy()
    start = 0
    for group in _mixed_blocks(G, deg_shells):
        bl = group[0]
        end = start + G[bl].data.size
        G[bl].data[...] = vec[start:end].reshape(G[bl].data.shape)
        start = end
        for other_block in group[1:]:
            G[other_block] << G[bl]
    return G


def _broyden_update(broyler, alpha, broy_max_it, deg_shells, gf_out, gf_in):
    """
    modified Broyden update of a BlockGf, D. D. Johnson, PRB 38, 12807 (1988),
    see also R. Zitko, PRB 80, 125125 (2009), with w_0 = 0.01 and all weights 1.
    It solves F(V) = dmft_step(V) - V = 0 for V = G0 or Sigma.

    Parameters
    ----------
    broyler : dict
        Broyden history with the lists 'V', 'F', 'dV', 'dF' of 1d vectors;
        updated in place. An empty history starts a new sequence with a
        linear mixing step.
    alpha : float
        mixing parameter of the linear part
    broy_max_it : int
        the last broy_max_it-1 differences are used; -1 uses all of them
    deg_shells : list of list of str
        degenerate blocks, only the first block of each shell is mixed
    gf_out : BlockGf
        dmft_step(V) of the last input V
    gf_in : BlockGf
        last input V, only used to start a new history

    Returns
    -------
    gf_update : BlockGf
        next input V
    broyler : dict
        updated history
    """
    w_0 = 0.01

    if not broyler['V']:
        broyler['V'].append(_bgf_to_vec(gf_in, deg_shells))
        broyler['F'].append(_bgf_to_vec(gf_out, deg_shells) - broyler['V'][-1])
        broyden_vec = np.zeros(broyler['V'][0].shape[0], dtype=complex)
    else:
        broyden_vec = np.zeros(broyler['V'][0].shape[0], dtype=complex)

        broyler['F'].append(_bgf_to_vec(gf_out, deg_shells) - broyler['V'][-1])
        norm_F = np.linalg.norm(broyler['F'][-1] - broyler['F'][-2])
        print('broyden norm: {:.3f}'.format(norm_F))
        broyler['dF'].append((broyler['F'][-1] - broyler['F'][-2]) / norm_F)
        broyler['dV'].append((broyler['V'][-1] - broyler['V'][-2]) / norm_F)
        m = len(broyler['dF'])

        start = 0
        if broy_max_it != -1 and m > broy_max_it - 1:
            start = m - (broy_max_it - 1)
        mat_dim = m - start

        A = np.zeros((mat_dim, mat_dim), dtype=complex)
        gamma = np.zeros(mat_dim, dtype=complex)
        for k_i in range(start, m):
            for n_i in range(start, m):
                A[k_i-start, n_i-start] = np.dot(broyler['dF'][n_i].conjugate(), broyler['dF'][k_i])
                if n_i == k_i:
                    A[k_i-start, n_i-start] += w_0**2
        beta_m = np.linalg.inv(A)
        for n_i in range(start, m):
            U_n = alpha * broyler['dF'][n_i] + broyler['dV'][n_i]
            for k_i in range(start, m):
                c_k = np.dot(broyler['dF'][k_i].conjugate(), broyler['F'][-1])
                gamma[n_i-start] += c_k * beta_m[k_i-start, n_i-start]
            broyden_vec += gamma[n_i-start] * U_n

    broyler['V'].append(broyler['V'][-1] + alpha * broyler['F'][-1] - broyden_vec)

    gf_update = gf_out.copy()
    gf_update << _vec_to_bgf(broyler['V'][-1], gf_out, deg_shells)
    return gf_update, broyler


def mix_g0(solver, general_params, icrsh, archive, G0_freq_previous, deg_shell):
    if general_params['g0_mix_type'] == 'linear' and general_params['g0_mix'] < 1.0:
        mpi.report('\nlinear mixing G0 with previous iteration by factor {:.3f}\n'.format(general_params['g0_mix']))
        solver.G0_freq << (general_params['g0_mix'] * solver.G0_freq
                                   + (1-general_params['g0_mix']) * G0_freq_previous)
    elif general_params['g0_mix_type'] == 'broyden':
        mpi.report('\nbroyden mixing G0 with previous iteration with alpha {:.3f}'.format(general_params['g0_mix']))
        mpi.report('\n############\n!!!! WARNING !!!! broyden mixing is still in early testing stage ! Use with caution.\n############\n')
        # TODO implement broyden mixing for Sigma
        if mpi.is_master_node():
            with HDFArchive(archive, 'a') as ar:
                broyler = ar['DMFT_results']['broyler']
                # calculate the next G0 via broyden scheme
                G0_broyden_update, broyler[icrsh] = _broyden_update(broyler[icrsh], general_params['g0_mix'],
                                                                    general_params['broy_max_it'],
                                                                    deg_shell, solver.G0_freq,
                                                                    G0_freq_previous)
                # store broyden update to h5 archive
                ar['DMFT_results']['broyler'] = broyler
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
