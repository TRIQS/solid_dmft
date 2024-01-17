#!/usr/bin/env python3
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
Contains all functions related to constructing the interaction Hamiltonian.
"""

# system
import os
import numpy as np
from itertools import product

# triqs
from h5 import HDFArchive
import triqs.utility.mpi as mpi
from triqs.gf import make_gf_imfreq
from triqs.operators import util, n, c, c_dag, Operator
from solid_dmft.dmft_tools import solver
try:
    import forktps as ftps
except ImportError:
    pass


def _extract_U_J_list(param_name, n_inequiv_shells, general_params):
    """
    Checks if param_name ('U', 'U_prime', or 'J') are a single value or
    different per inequivalent shell. If just a single value is given,
    this value is applied to each shell.
    """

    if not isinstance(param_name, str):
        formatted_param = ['none' if p == 'none' else '{:.2f}'.format(p)
                           for p in general_params[param_name]]
    else:
        formatted_param = general_params[param_name]

    if len(general_params[param_name]) == 1:
        mpi.report('Assuming {} = '.format(param_name)
                   + '{} for all correlated shells'.format(formatted_param[0]))
        general_params[param_name] *= n_inequiv_shells
    elif len(general_params[param_name]) == n_inequiv_shells:
        mpi.report('{} list for correlated shells: {}'.format(param_name, formatted_param))
    else:
        raise IndexError('Property list {} '.format(general_params[param_name])
                         + 'must have length 1 or n_inequiv_shells')

    return general_params


def _load_crpa_interaction_matrix(sum_k, general_params, filename='UIJKL'):
    """
    Loads  dynamic interaction data to use as an interaction Hamiltonian.
    """
    def _round_to_int(data):
        return (np.array(data) + .5).astype(int)

    if general_params['crpa_code'] == 'vasp':
    # Loads data from VASP cRPA file
        print('Loading Vasp cRPA matrix from file: '+str(filename))
        data = np.loadtxt(filename, unpack=True)
        u_matrix_four_indices = np.zeros(_round_to_int(np.max(data[:4], axis=1)), dtype=complex)
        for entry in data.T:
            # VASP switches the order of the indices, ijkl -> ikjl
            i, k, j, l = _round_to_int(entry[:4])-1
            u_matrix_four_indices[i, j, k, l] = entry[4] + 1j * entry[5]

        # Slices up the four index U-matrix, separating shells
        u_matrix_four_indices_per_shell = [None] * sum_k.n_inequiv_shells
        first_index_shell = 0
        for ish in range(sum_k.n_corr_shells):
            icrsh = sum_k.corr_to_inequiv[ish]
            n_orb = solver.get_n_orbitals(sum_k)[icrsh]['up']
            u_matrix_temp = u_matrix_four_indices[first_index_shell:first_index_shell+n_orb,
                                                  first_index_shell:first_index_shell+n_orb,
                                                  first_index_shell:first_index_shell+n_orb,
                                                  first_index_shell:first_index_shell+n_orb]
            # I think for now we should stick with real interactions make real
            u_matrix_temp.imag = 0.0

            if ish == icrsh:
                u_matrix_four_indices_per_shell[icrsh] = u_matrix_temp
            elif not np.allclose(u_matrix_four_indices_per_shell[icrsh], u_matrix_temp, atol=1e-6, rtol=0):
                # TODO: for some reason, some entries in the matrices differ by a sign. Check that
                # mpi.report(np.allclose(np.abs(u_matrix_four_indices_per_shell[icrsh]), np.abs(u_matrix_temp),
                # atol=1e-6, rtol=0))
                mpi.report('Warning: cRPA matrix for impurity {} '.format(icrsh)
                           + 'differs for shells {} and {}'.format(sum_k.inequiv_to_corr[icrsh], ish))

            first_index_shell += n_orb

        if not np.allclose(u_matrix_four_indices.shape, first_index_shell):
            print('Warning: different number of orbitals in cRPA matrix than in calculation.')

    if general_params['crpa_code'] == 'bdft':
        u_matrix_four_indices_per_shell = []
        for icrsh in range(sum_k.n_inequiv_shells):
            # for now we assume that up / down are equal
            if general_params['h_int_type'][icrsh] == 'crpa_density_density':
                Uloc_0 = make_gf_imfreq(general_params['Uloc_dlr'][icrsh]['up_0'],1)
                u_matrix_four_indices_per_shell.append(Uloc_0.data[0,:,:,:,:] + general_params['Vloc'][icrsh]['up_0'])
            else:
                # Uloc_0 = make_gf_imfreq(general_params['Uloc_dlr'][icrsh]['up_0'],1)
                # u_matrix_four_indices_per_shell.append(Uloc_0.data[0,:,:,:,:] + general_params['Vloc'][icrsh]['up_0'])
                u_matrix_four_indices_per_shell.append(general_params['Vloc'][icrsh]['up_0'])


    return u_matrix_four_indices_per_shell


# def _adapt_U_2index_for_SO(Umat, Upmat):
#     """
#     Changes the two-index U matrices such that for a system consisting of a
#     single block 'ud' with the entries (1, up), (1, down), (2, up), (2, down),
#     ... the matrices are consistent with the case without spin-orbit coupling.

#     Parameters
#     ----------
#     Umat : numpy array
#         The two-index interaction matrix for parallel spins without SO.
#     Upmat : numpy array
#         The two-index interaction matrix for antiparallel spins without SO.

#     Returns
#     -------
#     Umat_SO : numpy array
#         The two-index interaction matrix for parallel spins. Because in SO all
#         entries have nominal spin 'ud', this matrix now contains the original
#         Umat and Upmat.
#     Upmat_SO : numpy array
#         The two-index interaction matrix for antiparallel spins. Unused because
#         in SO, all spins have the same nominal spin 'ud'.
#     """

#     Umat_SO = np.zeros(np.array(Umat.shape)*2, dtype=Umat.dtype)
#     Umat_SO[::2, ::2] = Umat_SO[1::2, 1::2] = Umat
#     Umat_SO[::2, 1::2] = Umat_SO[1::2, ::2] = Upmat
#     Upmat_SO = None

#     return Umat_SO, Upmat_SO


def _adapt_U_4index_for_SO(Umat_full):
    """
    Changes the four-index U matrix such that for a system consisting of a
    single block 'ud' with the entries (1, up), (1, down), (2, up), (2, down),
    ... the matrix is consistent with the case without spin-orbit coupling.
    This can be derived directly from the definition of the Slater Hamiltonian.

    Parameters
    ----------
    Umat_full : numpy array
       The four-index interaction matrix without SO.

    Returns
    -------
    Umat_full_SO : numpy array
        The four-index interaction matrix with SO. For a matrix U_ijkl, the
        indices i, k correspond to spin sigma, and indices j, l to sigma'.
    """

    Umat_full_SO = np.zeros(np.array(Umat_full.shape)*2, dtype=Umat_full.dtype)
    for spin, spin_prime in ((0, 0), (0, 1), (1, 0), (1, 1)):
        Umat_full_SO[spin::2, spin_prime::2, spin::2, spin_prime::2] = Umat_full

    return Umat_full_SO


def _construct_kanamori(sum_k, general_params, icrsh):
    """
    Constructs the Kanamori interaction Hamiltonian. Only Kanamori does not
    need the full four-index matrix. Therefore, we can construct it directly
    from the parameters U and J.
    """

    n_orb = solver.get_n_orbitals(sum_k)[icrsh]['up']
    if sum_k.SO == 1:
        assert n_orb % 2 == 0
        n_orb = n_orb // 2

    if n_orb not in (2, 3):
        mpi.report('warning: are you sure you want to use the kanamori hamiltonian '
                   + 'outside the t2g or eg manifold?')

    # check if Uprime has been specified manually
    if general_params['U_prime'][icrsh] == 'U-2J':
        U_prime = general_params['U'][icrsh] - 2.0 * general_params['J'][icrsh]
    else:
        U_prime = general_params['U_prime'][icrsh]

    if general_params['solver_type'] == 'ftps':
        # 1-band modell requires J and U' equals zero
        if n_orb == 1:
            up, j = 0.0, 0.0
        else:
            up = U_prime
            j = general_params['J'][icrsh]
        h_int = ftps.solver_core.HInt(u=general_params['U'][icrsh], j=j, up=up, dd=False)
    elif sum_k.SO == 0:
        # Constructs U matrix
        Umat, Upmat = util.U_matrix_kanamori(n_orb=n_orb, U_int=general_params['U'][icrsh],
                                             J_hund=general_params['J'][icrsh],
                                             Up_int=U_prime)

        h_int = util.h_int_kanamori(sum_k.spin_block_names[sum_k.SO], n_orb,
                                    map_operator_structure=sum_k.sumk_to_solver[icrsh],
                                    U=Umat, Uprime=Upmat, J_hund=general_params['J'][icrsh],
                                    H_dump=os.path.join(general_params['jobname'], 'H.txt'))
    else:
        h_int = _construct_kanamori_soc(general_params['U'][icrsh], general_params['J'][icrsh],
                                        n_orb, sum_k.sumk_to_solver[icrsh],
                                        os.path.join(general_params['jobname'], 'H.txt'))
    return h_int


def _construct_kanamori_soc(U_int, J_hund, n_orb, map_operator_structure, H_dump=None):
    r"""
    Adapted from triqs.operators.util.hamiltonians.h_int_kanamori. Assumes
    that spin_names == ['ud'] and that map_operator_structure is given.
    """

    orb_names = list(range(n_orb))

    if H_dump:
        H_dump_file = open(H_dump, 'w')
        H_dump_file.write("Kanamori Hamiltonian:" + '\n')

    H = Operator()
    mkind = util.op_struct.get_mkind(None, map_operator_structure)

    s = 'ud'

    # density terms:
    # TODO: reformulate in terms of Umat and Upmat for consistency with triqs?
    if H_dump:
        H_dump_file.write("Density-density terms:" + '\n')
    for a1, a2 in product(orb_names, orb_names):
        if a1 == a2:  # same spin and orbital
            continue

        if a1 // 2 == a2 // 2:  # same orbital (, different spins)
            U_val = U_int
        elif a1 % 2 != a2 % 2:  # different spins (, different orbitals)
            U_val = U_int - 2*J_hund
        else:  # same spins (, different orbitals)
            U_val = U_int - 3*J_hund

        H_term = 0.5 * U_val * n(*mkind(s, a1)) * n(*mkind(s, a2))
        H += H_term

        # Dump terms of H
        if H_dump and not H_term.is_zero():
            H_dump_file.write('%s' % (mkind(s, a1), ) + '\t')
            H_dump_file.write('%s' % (mkind(s, a2), ) + '\t')
            H_dump_file.write(str(U_val) + '\n')

    # spin-flip terms:
    if H_dump:
        H_dump_file.write("Spin-flip terms:" + '\n')
    for a1, a2, a3, a4 in product(orb_names, orb_names, orb_names, orb_names):
        if a1 == a2 or a1 == a3 or a1 == a4 or a2 == a3 or a2 == a4 or a3 == a4:
            continue

        if not (a1//2 == a2//2 and a3//2 == a4//2 and a1//2 != a3//2 and a1 % 2 != a3 % 2):
            continue

        H_term = -0.5 * J_hund * c_dag(*mkind(s, a1)) * c(*mkind(s, a2)) * c_dag(*mkind(s, a3)) * c(*mkind(s, a4))
        H += H_term

        # Dump terms of H
        if H_dump and not H_term.is_zero():
            H_dump_file.write('%s' % (mkind(s, a1), ) + '\t')
            H_dump_file.write('%s' % (mkind(s, a2), ) + '\t')
            H_dump_file.write('%s' % (mkind(s, a3), ) + '\t')
            H_dump_file.write('%s' % (mkind(s, a4), ) + '\t')
            H_dump_file.write(str(-J_hund) + '\n')

    # pair-hopping terms:
    if H_dump:
        H_dump_file.write("Pair-hopping terms:" + '\n')
    for a1, a2, a3, a4 in product(orb_names, orb_names, orb_names, orb_names):
        if a1 == a2 or a1 == a3 or a1 == a4 or a2 == a3 or a2 == a4 or a3 == a4:
            continue

        if not (a1//2 == a2//2 and a3//2 == a4//2 and a1//2 != a3//2 and a1 % 2 != a3 % 2):
            continue

        H_term = 0.5 * J_hund * c_dag(*mkind(s, a1)) * c_dag(*mkind(s, a2)) * c(*mkind(s, a4)) * c(*mkind(s, a3))
        H += H_term

        # Dump terms of H
        if H_dump and not H_term.is_zero():
            H_dump_file.write('%s' % (mkind(s, a1), ) + '\t')
            H_dump_file.write('%s' % (mkind(s, a2), ) + '\t')
            H_dump_file.write('%s' % (mkind(s, a3), ) + '\t')
            H_dump_file.write('%s' % (mkind(s, a4), ) + '\t')
            H_dump_file.write(str(-J_hund) + '\n')

    return H


def _construct_dynamic(sum_k, general_params, icrsh):
    """
    Constructs the interaction Hamiltonian for a frequency-dependent interaction.
    Works only without spin-orbit coupling and only for one orbital.
    """

    mpi.report('###### Dynamic U calculation ######, load parameters from input archive.')
    U_onsite = None
    if mpi.is_master_node():
        with HDFArchive(general_params['jobname']+'/'+general_params['seedname']+'.h5', 'r') as archive:
            U_onsite = archive['dynamic_U']['U_scr']
    U_onsite = mpi.bcast(U_onsite)

    n_orb = solver.get_n_orbitals(sum_k)[icrsh]['up']
    if sum_k.SO == 1:
        raise ValueError('dynamic U not implemented for SO!=0')
    if n_orb > 1:
        raise ValueError('dynamic U not implemented for more than one orbital')

    mpi.report('onsite interaction value for imp {}: {:.3f}'.format(icrsh, U_onsite[icrsh]))
    h_int = util.h_int_density(sum_k.spin_block_names[sum_k.SO], n_orb,
                               map_operator_structure=sum_k.sumk_to_solver[icrsh],
                               U=np.array([[0]]), Uprime=np.array([[U_onsite[icrsh]]]), H_dump=os.path.join(general_params['jobname'], 'H.txt'))

    return h_int


def _generate_four_index_u_matrix(sum_k, general_params, icrsh):
    """
    Generates the four-index interaction matrix per impurity with the interaction
    parameters U and J (and ratio_F4_F2 for the d shell).
    """

    # ish points to the shell representative of the current group
    ish = sum_k.inequiv_to_corr[icrsh]
    n_orb = solver.get_n_orbitals(sum_k)[icrsh]['up']
    if sum_k.SO == 1:
        assert n_orb % 2 == 0
        n_orb = n_orb // 2

    if sum_k.corr_shells[ish]['l'] != 2:
        slater_integrals = util.U_J_to_radial_integrals(l=sum_k.corr_shells[ish]['l'],
                                                        U_int=general_params['U'][icrsh],
                                                        J_hund=general_params['J'][icrsh])
    else:
        # Implements parameter R=F4/F2. For R=0.63 equivalent to util.U_J_to_radial_integrals
        U = general_params['U'][icrsh]
        J = general_params['J'][icrsh]
        R = general_params['ratio_F4_F2'][icrsh]
        R = 0.63 if R == 'none' else R
        slater_integrals = np.array([U, 14*J/(1+R), 14*J*R/(1+R)])

    mpi.report('\nImpurity {}: The corresponding slater integrals are'.format(icrsh))
    formatted_slater_integrals = [y for x in list(zip([2*x for x in range(len(slater_integrals))], slater_integrals)) for y in x]
    mpi.report(('F{:d} = {:.2f}, '*len(slater_integrals)).format(*formatted_slater_integrals))

    # Constructs U matrix
    # construct full spherical symmetric U matrix and transform to cubic basis
    # the order for the cubic orbitals is given by the convention. The TRIQS
    # convention is as follows ("xy","yz","z^2","xz","x^2-y^2")
    # this is consistent with the order of orbitals in the VASP interface
    # but not necessarily with wannier90, qe, and wien2k!
    # This is also true for the f-shell.
    Umat_full = util.U_matrix_slater(l=sum_k.corr_shells[ish]['l'],
                              radial_integrals=slater_integrals, basis='spherical')
    Umat_full = util.transform_U_matrix(Umat_full,
                                        util.spherical_to_cubic(l=sum_k.corr_shells[ish]['l'],
                                                                convention=general_params['h_int_basis'])
                                        )

    if n_orb == 2:
        Umat_full = util.eg_submatrix(Umat_full)
        mpi.report('Using eg subspace of interaction Hamiltonian')
    elif n_orb == 3:
        Umat_full = util.t2g_submatrix(Umat_full)
        mpi.report('Using t2g subspace of interaction Hamiltonian')
    elif n_orb not in (5, 7):
        raise ValueError('Calculations for d shell only support 2, 3 or 5 orbitals'
                         + 'and for the f shell only 7 orbitals')

    return Umat_full


def _rotate_four_index_matrix(sum_k, general_params, Umat_full, icrsh):
    """ Rotates the four index matrix into the local frame. """

    ish = sum_k.inequiv_to_corr[icrsh]
    # Transposes rotation matrix here because TRIQS has a slightly different definition
    Umat_full_rotated = util.transform_U_matrix(Umat_full, sum_k.rot_mat[ish].T)

    if general_params['h_int_type'][icrsh] in ('density_density', 'crpa_density_density', 'dyn_density_density'):
        if not np.allclose(Umat_full_rotated, Umat_full):
            mpi.report('WARNING: applying a rotation matrix changes the dens-dens Hamiltonian.\n'
                       + 'This changes the definition of the ignored spin flip and pair hopping.')
    elif general_params['h_int_type'][icrsh] in ('full_slater', 'crpa', 'dyn_full'):
        if not np.allclose(Umat_full_rotated, Umat_full):
            mpi.report('WARNING: applying a rotation matrix changes the interaction Hamiltonian.\n'
                       + 'Please ensure that the rotation is correct!')

    return Umat_full_rotated


def _construct_density_density(sum_k, general_params, Umat_full_rotated, icrsh):
    """
    Constructs the density-density Slater-Hamiltonian from the four-index
    interaction matrix.
    """

    # Constructs Hamiltonian from Umat_full_rotated
    n_orb = solver.get_n_orbitals(sum_k)[icrsh]['up']

    Umat, Upmat = util.reduce_4index_to_2index(Umat_full_rotated)
    h_int = util.h_int_density(sum_k.spin_block_names[sum_k.SO], n_orb,
                               map_operator_structure=sum_k.sumk_to_solver[icrsh],
                               U=Umat, Uprime=Upmat, H_dump=os.path.join(general_params['jobname'], 'H.txt'))

    return h_int


def _construct_slater(sum_k, general_params, Umat_full_rotated, icrsh):
    """
    Constructs the full Slater-Hamiltonian from the four-index interaction
    matrix.
    """

    n_orb = solver.get_n_orbitals(sum_k)[icrsh]['up']

    h_int = util.h_int_slater(sum_k.spin_block_names[sum_k.SO], n_orb,
                              map_operator_structure=sum_k.sumk_to_solver[icrsh],
                              U_matrix=Umat_full_rotated,
                              H_dump=os.path.join(general_params['jobname'], 'H.txt'))

    return h_int

def h_int_simple_intra(spin_names,n_orb,U,off_diag=None,map_operator_structure=None,H_dump=None):
    r"""
    Create a simple intra orbital density-density Hamiltonian.
    (no inter orbital terms)

    .. math::
        H = \frac{1}{2} \sum_{i \sigma \neq \sigma')} U_{i i}^{\sigma \sigma'} n_{i \sigma} n_{i \sigma'}.

    Parameters
    ----------
    spin_names : list of strings
               Names of the spins, e.g. ['up','down'].
    n_orb : int
               Number of orbitals.
    U : float
               U value
    off_diag : boolean
               Do we have (orbital) off-diagonal elements?
               If yes, the operators and blocks are denoted by ('spin', 'orbital'),
               otherwise by ('spin_orbital',0).
    map_operator_structure : dict
               Mapping of names of GF blocks names from one convention to another,
               e.g. {('up', 0): ('up_0', 0), ('down', 0): ('down_0',0)}.
               If provided, the operators and blocks are denoted by the mapping of ``('spin', 'orbital')``.
    H_dump : string
               Name of the file to which the Hamiltonian should be written.

    Returns
    -------
    H : Operator
        The Hamiltonian.

    """
    from triqs.operators.util.op_struct import get_mkind

    if H_dump:
        H_dump_file = open(H_dump,'w')
        H_dump_file.write("Density-density Hamiltonian:" + '\n')

    H = Operator()
    mkind = get_mkind(off_diag,map_operator_structure)
    if H_dump: H_dump_file.write("Density-density terms:" + '\n')
    for s1, s2 in product(spin_names,spin_names):
        if (s1 is not s2):
            for a1 in range(n_orb):
                H_term = 0.5 * U * n(*mkind(s1,a1)) * n(*mkind(s2,a1))
                H += H_term

                # Dump terms of H
                if H_dump and not H_term.is_zero():
                    H_dump_file.write('%s'%(mkind(s1,a1),) + '\t')
                    H_dump_file.write('%s'%(mkind(s2,a1),) + '\t')
                    H_dump_file.write(str(U) + '\n')

    return H


def construct(sum_k, general_params, advanced_params):
    """
    Constructs the interaction Hamiltonian. Currently implemented are the
    Kanamori Hamiltonian (usually for 2 or 3 orbitals), the density-density and
    the full Slater Hamiltonian (for 2, 3, or 5 orbitals).
    If sum_k.rot_mat is non-identity, we have to consider rotating the interaction
    Hamiltonian: the Kanamori Hamiltonian does not change because it is invariant
    under orbital mixing but all the other Hamiltonians are at most invariant
    under rotations in space. Therefore, sum_k.rot_mat has to be correct before
    calling this method.

    The parameters U and J will be interpreted differently depending on the
    type of the interaction Hamiltonian: it is either the Kanamori parameters
    for the Kanamori Hamiltonian or the orbital-averaged parameters (consistent
    with DFT+U, https://cms.mpi.univie.ac.at/wiki/index.php/LDAUTYPE ) for all
    other Hamiltonians.

    Note also that for all Hamiltonians except Kanamori, the order of the
    orbitals matters. The correct order is specified here:
    triqs.github.io/triqs/unstable/documentation/python_api/triqs.operators.util.U_matrix.spherical_to_cubic.html
    """

    # Extracts U and J
    mpi.report('*** interaction parameters ***')

    if any(general_params['h_int_type']) in ('density_density', 'kanamori', 'full_slater', 'ntot'):
        general_params = _extract_U_J_list('h_int_type', sum_k.n_inequiv_shells, general_params)
        for param_name in ('U', 'J'):
            general_params = _extract_U_J_list(param_name, sum_k.n_inequiv_shells, general_params)
        if 'kanamori' in general_params['h_int_type']:
            general_params = _extract_U_J_list('U_prime', sum_k.n_inequiv_shells, general_params)
        for param_name in ('dc_U', 'dc_J'):
            advanced_params = _extract_U_J_list(param_name, sum_k.n_inequiv_shells, advanced_params)

    # Extracts ratio_F4_F2 if any shell uses a solver supporting it
    if 'density_density' in general_params['h_int_type'] or 'full_slater' in general_params['h_int_type']:
        general_params = _extract_U_J_list('ratio_F4_F2', sum_k.n_inequiv_shells, general_params)

    # Constructs the interaction Hamiltonian. Needs to come after setting sum_k.rot_mat
    mpi.report('\nConstructing the interaction Hamiltonians')
    h_int = [None] * sum_k.n_inequiv_shells
    for icrsh in range(sum_k.n_inequiv_shells):
        mpi.report('\nImpurity {}: constructing a {} type interaction Hamiltonian'.format(icrsh, general_params['h_int_type'][icrsh]))

        # Kanamori
        if general_params['h_int_type'][icrsh] == 'kanamori':
            h_int[icrsh] = _construct_kanamori(sum_k, general_params, icrsh)
            continue

        # for density density or full slater get full four-index U matrix
        if general_params['h_int_type'][icrsh] in ('density_density', 'full_slater'):
            mpi.report('\nNote: The input parameters U and J here are orbital-averaged parameters.')
            mpi.report('Note: The order of the orbitals is important. See also the doc string of this method.')
            # Checks that all entries are l == 2 or R == 'none'
            if (sum_k.corr_shells[sum_k.inequiv_to_corr[icrsh]]['l'] != 2
                    and general_params['ratio_F4_F2'][icrsh] != 'none'):
                raise ValueError('Ratio F4/F2 only implemented for d-shells '
                                 + 'but set in impurity {}'.format(icrsh))

            if general_params['h_int_type'][icrsh] == 'density_density' and general_params['solver_type'] == 'ftps':
                # TODO: implement
                raise NotImplementedError('\nNote: Density-density not implemented for ftps.')

            Umat_full = _generate_four_index_u_matrix(sum_k, general_params, icrsh)

            if sum_k.SO == 1:
                Umat_full = [_adapt_U_4index_for_SO(Umat_full_per_imp)
                             for Umat_full_per_imp in Umat_full]

            # Rotates the interaction matrix
            Umat_full_rotated = _rotate_four_index_matrix(sum_k, general_params, Umat_full, icrsh)

            # construct slater / density density from U tensor
            if general_params['h_int_type'][icrsh] == 'full_slater':
                h_int[icrsh] = _construct_slater(sum_k, general_params, Umat_full_rotated, icrsh)
            else:
                h_int[icrsh] = _construct_density_density(sum_k, general_params, Umat_full_rotated, icrsh)

            continue

        # simple total impurity occupation interation: U/2 (Ntot^2 - Ntot)
        if general_params['h_int_type'][icrsh] == 'ntot':
            n_tot_op = Operator()
            for block, n_orb in sum_k.gf_struct_solver[icrsh].items():
                n_tot_op += sum(n(block, orb) for orb in range(n_orb))
            h_int[icrsh] = general_params['U'][icrsh]/2.0 * (n_tot_op*n_tot_op - n_tot_op)
            continue

        if general_params['h_int_type'][icrsh] == 'simple_intra':
            h_int[icrsh] = h_int_simple_intra(sum_k.spin_block_names[sum_k.SO],
                                              solver.get_n_orbitals(sum_k)[icrsh]['up'],
                                              map_operator_structure=sum_k.sumk_to_solver[icrsh],
                                              U=general_params['U'][icrsh],
                                              H_dump=os.path.join(general_params['jobname'], 'H.txt'))
            continue


        # read from file options
        if general_params['h_int_type'][icrsh] in ('crpa', 'crpa_density_density', 'dyn_density_density', 'dyn_full'):
            Umat_full = _load_crpa_interaction_matrix(sum_k, general_params)[icrsh]

            if sum_k.SO == 1:
                Umat_full = _adapt_U_4index_for_SO(Umat_full)

            # Rotates the interaction matrix
            Umat_full_rotated = _rotate_four_index_matrix(sum_k, general_params, Umat_full, icrsh)

            # construct slater / density density from U tensor
            if general_params['h_int_type'][icrsh] in ('crpa', 'dyn_full'):
                h_int[icrsh] = _construct_slater(sum_k, general_params, Umat_full_rotated.real, icrsh)
            else:
                h_int[icrsh] = _construct_density_density(sum_k, general_params, Umat_full_rotated.real, icrsh)
            continue

        raise NotImplementedError('Error when constructing the interaction Hamiltonian.')

    return h_int
