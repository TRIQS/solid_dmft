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
Contains formatters for things that need to be printed in DMFT calculations.
"""

import numpy as np
import triqs.utility.mpi as mpi

def print_rotation_matrix(sum_k):
    """
    Prints the rotation matrix, real and imaginary part separately.
    """
    if not mpi.is_master_node():
        return

    for icrsh, rot_crsh in enumerate(sum_k.rot_mat):
        n_orb = sum_k.corr_shells[icrsh]['dim']
        print('rot_mat[{:2d}] '.format(icrsh)+'real part'.center(9*n_orb)+'  '+'imaginary part'.center(9*n_orb))
        fmt = '{:9.5f}' * n_orb
        for row in rot_crsh:
            row = np.concatenate((row.real, row.imag))
            print((' '*11 + fmt + '  ' + fmt).format(*row))
    print('\n')


def print_block_sym(sum_k, dm, general_params):
    """
    Prints a summary of block structure finder, determination of
    shell_multiplicity, local Hamiltonian, DFT density matrix.
    """
    if not mpi.is_master_node():
        return

    print('\n number of ineq. correlated shells: {}'.format(sum_k.n_inequiv_shells))
    # correlated shells and their structure
    print('\n block structure summary')
    for icrsh in range(sum_k.n_inequiv_shells):
        shlst = [ish for ish, ineq_shell in enumerate(sum_k.corr_to_inequiv) if ineq_shell == icrsh]
        print(' -- Shell type #{:3d}: '.format(icrsh) + format(shlst))
        print('  | shell multiplicity '+str(len(shlst)))
        print('  | block struct. sum_k : ' + format(sum_k.gf_struct_sumk[sum_k.inequiv_to_corr[icrsh]]))
        print('  | block struct. solver: ' + format(sum_k.gf_struct_solver[icrsh]))
        print('  | deg. orbitals : ' + format(sum_k.deg_shells[icrsh]))

    # Prints matrices
    print('\nRotation matrices')
    print_rotation_matrix(sum_k)

    # Prints dft local Hamiltonian and occupations
    print('\nEffective atomic levels')
    eal = sum_k.eff_atomic_levels()
    for icrsh, eal_crsh in enumerate(eal):
        shlst = [ish for ish, ineq_shell in enumerate(sum_k.corr_to_inequiv) if ineq_shell == icrsh]
        n_orb = sum_k.corr_shells[sum_k.inequiv_to_corr[icrsh]]['dim']
        spins = sum_k.spin_block_names[sum_k.SO]
        if sum_k.SO == 0:
            eal_blocks = {spin: eal_crsh[spin] for spin in spins}
        else:
            eal_blocks = {'ud real': eal_crsh['ud'].real, 'ud imag': eal_crsh['ud'].imag}
        print('H_loc[{:2d}] '.format(icrsh)+2*' '+' '.join([sp.center(9*n_orb) for sp in eal_blocks.keys()]))
        fmt = '{:9.5f}' * n_orb
        for block_1, block_2 in zip(*eal_blocks.values()):
            row = np.concatenate((block_1.real, block_2.real))
            print((' '*11 + fmt + '  ' + fmt).format(*row))
    print('\n')

    if dm:
        # Prints dft density matrix
        print('\nDFT density matrix')
        # Double loop to sort by impurities
        for icrsh in range(sum_k.n_inequiv_shells):
            spins = sum_k.spin_block_names[sum_k.SO]
            shlst = [ish for ish, ineq_shell in enumerate(sum_k.corr_to_inequiv) if ineq_shell == icrsh]
            for sh in shlst:
                n_orb = sum_k.corr_shells[sh]['dim']
                if sum_k.SO == 0:
                    dm_blocks = {spin: dm[sh][spin] for spin in spins}
                else:
                    dm_blocks = {'ud real': dm[sh]['ud'].real, 'ud imag': dm[sh]['ud'].imag}
                print('rho[{0:2d}] '.format(sh)+4*' '+' '.join([sp.center(9*n_orb) for sp in dm_blocks.keys()]))
                fmt = '{:9.5f}' * n_orb
                for block_1, block_2 in zip(*dm_blocks.values()):
                    row = np.concatenate((block_1.real, block_2.real))
                    print((' '*11 + fmt + '  ' + fmt).format(*row))
        print('\n')


def print_local_density(density, density_pre, density_mat, spin_orbit=False):
    if not mpi.is_master_node():
        return

    if spin_orbit:
        printed = ((np.real, 'real'), (np.imag, 'imaginary'))
    else:
        printed = ((np.real, 'real'), )

    print('\nTotal charge of impurity problem: {:7.5f}'.format(density))
    print('Total charge convergency of impurity problem: {:7.5f}'.format(density-density_pre))
    print('\nDensity matrix:')
    for key, value in sorted(density_mat.items()):
        for func, name in printed:
            print('{}, {} part'.format(key, name))
            print(func(value))
        eigenvalues = np.linalg.eigvalsh(value)
        print('eigenvalues: {}'.format(eigenvalues))
        # check for large off-diagonal elements and write out a warning
        if np.max(np.abs(value - np.diag(np.diag(value)))) >= 0.1:
            print('\n!!! WARNING !!!')
            print('!!! large off diagonal elements in density matrix detected! I hope you know what you are doing !!!')
            print('!!! WARNING !!!\n')


def print_summary_energetics(observables):
    if not mpi.is_master_node():
        return

    print('\n' + '='*60)
    print('summary of energetics:')
    print('total energy: ', observables['E_tot'][-1])
    print('DFT energy: ', observables['E_dft'][-1])
    print('correlation energy: ', observables['E_corr_en'][-1])
    print('DFT band correction: ', observables['E_bandcorr'][-1])
    print('='*60 + '\n')


def print_summary_observables(observables, n_inequiv_shells, spin_block_names):
    if not mpi.is_master_node():
        return

    print('='*60)
    print('summary of impurity observables:')
    for icrsh in range(n_inequiv_shells):
        total_occ = np.sum([observables['imp_occ'][icrsh][spin][-1] for spin in spin_block_names])
        print('total occupany of impurity {}: {:7.4f}'.format(icrsh, total_occ))
    for icrsh in range(n_inequiv_shells):
        total_gb2 = np.sum([observables['imp_gb2'][icrsh][spin][-1] for spin in spin_block_names])
        print('G(beta/2) occ of impurity {}: {:8.4f}'.format(icrsh, total_gb2))
    for icrsh in range(n_inequiv_shells):
        print('Z (simple estimate) of impurity {} per orb:'.format(icrsh))
        for spin in spin_block_names:
            Z_spin = observables['orb_Z'][icrsh][spin][-1]
            print('{:>5}: '.format(spin) + ' '.join("{:6.3f}".format(Z_orb) for Z_orb in Z_spin))
    print('='*60 + '\n')


def print_summary_magnetic_occ(observables, n_inequiv_shells):
    if not mpi.is_master_node():
        return

    occ = {'up': 0.0, 'down': 0.0}
    print('\n' + '='*60)
    print('\n *** summary of magnetic occupations: ***')
    for icrsh in range(n_inequiv_shells):
        for spin in ['up', 'down']:
            temp = observables['imp_occ'][icrsh][spin][-1]
            print('imp '+str(icrsh)+' spin '+spin+': {:6.4f}'.format(temp))
            occ[spin] += temp

    print('total spin up   occ: '+'{:6.4f}'.format(occ['up']))
    print('total spin down occ: '+'{:6.4f}'.format(occ['down']))
    print('='*60 + '\n')


def print_summary_convergence(conv_obs, general_params, n_inequiv_shells):
    if not mpi.is_master_node():
        return

    print('='*60)
    print('convergence:')
    print('δμ:      {:.4e}'.format(conv_obs['d_mu'][-1]))
    # if calc energies calc /print also the diff in Etot
    if general_params['calc_energies']:
        print('δE_tot:  {:.4e}'.format(conv_obs['d_Etot'][-1]))
        print("---")
    for icrsh in range(n_inequiv_shells):
        print('Impurity '+str(icrsh)+':')
        print('δn imp : {:.4e}'.format(conv_obs['d_imp_occ'][icrsh][-1]))
        print('δn orb : '+'  '.join("{:.4e}".format(orb) for orb in conv_obs['d_orb_occ'][icrsh][-1]))
        print('δ Gimp : {:.4e}'.format(conv_obs['d_Gimp'][icrsh][-1]))
        print('δ G0   : {:.4e}'.format(conv_obs['d_G0'][icrsh][-1]))
        print('δ Σ    : {:.4e}'.format(conv_obs['d_Sigma'][icrsh][-1]))

    print('='*60)
    print('\n')
