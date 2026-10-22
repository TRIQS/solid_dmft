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


def init_broyden_history(archive, quantity, n_inequiv_shells):
    """
    creates an empty Broyden history DMFT_results/broyler_<quantity> with one
    entry per inequivalent shell in the h5 archive, unless it exists already.
    Call on the master node only.
    """
    key = 'broyler_' + quantity
    with HDFArchive(archive, 'a') as ar:
        if key not in ar['DMFT_results']:
            ar['DMFT_results'][key] = [{'V': [], 'dV': [], 'F': [], 'dF': []}
                                       for _ in range(n_inequiv_shells)]


def mix_gf(gf_out, gf_in, general_params, quantity, icrsh, archive, deg_shells):
    """
    mixes the output gf_out = dmft_step(gf_in) of quantity 'G0' or 'Sigma' with
    the input gf_in of the same step, according to general_params['mix_type']
    and general_params['mix_greed']. Broyden mixing reads and writes the history
    of quantity and shell icrsh in the h5 archive (see init_broyden_history).

    Returns
    -------
    BlockGf
        the next input, the same on all MPI ranks
    """
    mix_type = general_params['mix_type']
    mix_greed = general_params['mix_greed']
    mixed = gf_out.copy()

    if mix_type == 'linear':
        if mix_greed < 1.0:
            mpi.report(f'\nlinear mixing {quantity} with previous iteration by factor {mix_greed:.3f}\n')
            mixed << mix_greed * gf_out + (1 - mix_greed) * gf_in
    elif mix_type == 'broyden':
        mpi.report(f'\nbroyden mixing {quantity} with previous iteration with alpha {mix_greed:.3f}')
        mpi.report('\n############\n!!!! WARNING !!!! broyden mixing is still in early testing stage ! Use with caution.\n############\n')
        if mpi.is_master_node():
            key = 'broyler_' + quantity
            with HDFArchive(archive, 'a') as ar:
                broyler = ar['DMFT_results'][key]
                mixed, broyler[icrsh] = _broyden_update(broyler[icrsh], mix_greed, general_params['broy_max_it'],
                                                        deg_shells, gf_out, gf_in)
                ar['DMFT_results'][key] = broyler
    else:
        raise ValueError(f'Invalid mix_type "{mix_type}".')

    mixed << mpi.bcast(mixed)
    return mixed
