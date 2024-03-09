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
Contains all functions related to the observables.
"""

# system
import os.path
import numpy as np

# triqs
import triqs.utility.mpi as mpi
from triqs.gf import Gf, MeshImTime
from triqs.atom_diag import trace_rho_op
from triqs.gf.descriptors import Fourier
from solid_dmft.dmft_tools import solver

def prep_observables(h5_archive, sum_k):
    """
    prepares the observable arrays and files for the DMFT calculation

    Parameters
    ----------
    h5_archive: hdf archive instance
        hdf archive for calculation
    sum_k : SumK Object instances

    Returns
    -------
    observables : dict
        observable array for calculation
    """

    # determine number of impurities
    n_inequiv_shells = h5_archive['dft_input']['n_inequiv_shells']

    # check for previous iterations
    obs_prev = []
    if 'observables' in h5_archive['DMFT_results']:
        obs_prev = h5_archive['DMFT_results']['observables']

    # prepare observable dicts
    if len(obs_prev) > 0:
        observables = obs_prev
    else:
        observables = dict()
        observables['iteration'] = []
        observables['mu'] = []
        observables['E_tot'] = []
        observables['E_bandcorr'] = []
        observables['E_int'] = [[] for _ in range(n_inequiv_shells)]
        observables['E_corr_en'] = []
        observables['E_dft'] = []
        observables['E_DC'] = [[] for _ in range(n_inequiv_shells)]
        observables['orb_gb2'] = [{spin: [] for spin in sum_k.spin_block_names[sum_k.SO]}
                                  for _ in range(n_inequiv_shells)]
        observables['imp_gb2'] = [{spin: [] for spin in sum_k.spin_block_names[sum_k.SO]}
                                  for _ in range(n_inequiv_shells)]
        observables['orb_occ'] = [{spin: [] for spin in sum_k.spin_block_names[sum_k.SO]}
                                  for _ in range(n_inequiv_shells)]
        observables['orb_Z'] = [{spin: [] for spin in sum_k.spin_block_names[sum_k.SO]}
                                  for _ in range(n_inequiv_shells)]
        observables['imp_occ'] = [{spin: [] for spin in sum_k.spin_block_names[sum_k.SO]}
                                  for _ in range(n_inequiv_shells)]

    return observables

def _generate_header(general_params, sum_k):
    """
    Generates the headers that are used in write_header_to_file.
    Returns a dict with {file_name: header_string}
    """
    n_orb = solver.get_n_orbitals(sum_k)

    header_energy_mask = ' | {:>10} | {:>10}   {:>10}   {:>10}   {:>10}'
    header_energy = header_energy_mask.format('E_tot', 'E_DFT', 'E_bandcorr', 'E_int_imp', 'E_DC')

    headers = {}
    for iineq in range(sum_k.n_inequiv_shells):
        number_spaces = max(10*n_orb[iineq]['up'] + 3*(n_orb[iineq]['up']-1), 21)
        header_basic_mask = '{{:>3}} | {{:>10}} | {{:>{0}}} | {{:>{0}}} | {{:>17}}'.format(number_spaces)

        # If magnetic calculation is done create two obs files per imp
        if general_params['magnetic'] and sum_k.SO == 0:
            for spin in ('up', 'down'):
                file_name = 'observables_imp{}_{}.dat'.format(iineq, spin)
                headers[file_name] = header_basic_mask.format('it', 'mu', 'G(beta/2) per orbital',
                                                             'orbital occs '+spin, 'impurity occ '+spin)

                if general_params['calc_energies']:
                    headers[file_name] += header_energy
        else:
            file_name = 'observables_imp{}.dat'.format(iineq)
            headers[file_name] = header_basic_mask.format('it', 'mu', 'G(beta/2) per orbital',
                                                         'orbital occs up+down', 'impurity occ')

            if general_params['calc_energies']:
                headers[file_name] += header_energy

    return headers


def write_header_to_file(general_params, sum_k):
    """
    Writes the header to the observable files

    Parameters
    ----------
    general_params : dict
        general parameters as a dict
    n_inequiv_shells : int
        number of impurities for calculations


    Returns
    -------
    nothing
    """

    headers = _generate_header(general_params, sum_k)

    for file_name, header in headers.items():
        path = os.path.join(general_params['jobname'], file_name)
        with open(path, 'w') as obs_file:
            obs_file.write(header + '\n')


def add_dft_values_as_zeroth_iteration(observables, general_params, solver_type_per_imp, dft_mu, dft_energy,
                                       sum_k, G_loc_all_dft, shell_multiplicity):
    """
    Calculates the DFT observables that should be written as the zeroth iteration.

    Parameters
    ----------
    observables : observable arrays/dicts

    general_params : general parameters as a dict

    solver_type_per_imp : list of strings
        list of solver types for each impurity

    dft_mu : dft chemical potential

    sum_k : SumK Object instances

    G_loc_all_dft : Gloc from DFT for G(beta/2)

    shell_multiplicity : degeneracy of impurities

    Returns
    -------

    observables: list of dicts
    """
    dft_energy = 0.0 if dft_energy is None else dft_energy
    density_mat_dft = [G_loc_all_dft[iineq].density() for iineq in range(sum_k.n_inequiv_shells)]
    observables['iteration'].append(0)
    observables['mu'].append(float(dft_mu))

    if general_params['calc_energies']:
        observables['E_bandcorr'].append(0.0)
        observables['E_corr_en'].append(0.0)
        observables['E_dft'].append(dft_energy)
    else:
        observables['E_bandcorr'].append('none')
        observables['E_corr_en'].append('none')
        observables['E_dft'].append('none')

    for iineq in range(sum_k.n_inequiv_shells):
        if general_params['calc_energies']:
            observables['E_int'][iineq].append(0.0)
            double_counting_energy = shell_multiplicity[iineq]*sum_k.dc_energ[sum_k.inequiv_to_corr[iineq]] if general_params['dc'] else 0.0
            observables['E_DC'][iineq].append(double_counting_energy)
        else:
            observables['E_int'][iineq].append('none')
            observables['E_DC'][iineq].append('none')

        # Collect all occupation and G(beta/2) for spin up and down separately
        for spin in sum_k.spin_block_names[sum_k.SO]:
            g_beta_half_per_impurity = 0.0
            g_beta_half_per_orbital = []
            occupation_per_impurity = 0.0
            occupation_per_orbital = []
            Z_per_orbital = []

            # iterate over all spin channels and add the to up or down
            # we need only the keys of the blocks, but sorted! Otherwise
            # orbitals will be mixed up_0 to down_0 etc.
            for spin_channel in sorted(sum_k.gf_struct_solver[iineq].keys()):
                if not spin in spin_channel:
                    continue

                if solver_type_per_imp[iineq] == 'ftps':
                    freq_mesh = np.array([w.value for w in G_loc_all_dft[iineq][spin_channel].mesh])
                    fermi_idx = abs(freq_mesh).argmin()
                    gb2_averaged = G_loc_all_dft[iineq][spin_channel].data[fermi_idx].imag

                    # Z is not defined without Sigma, adding 1
                    Z_per_orbital.extend( [1.0] * G_loc_all_dft[iineq][spin_channel].target_shape[0]  )
                else:
                    # G(beta/2)
                    mesh = MeshImTime(beta=general_params['beta'], S="Fermion",
                                      n_max=general_params['n_tau'])
                    G_time = Gf(mesh=mesh, indices=G_loc_all_dft[iineq][spin_channel].indices)
                    G_time << Fourier(G_loc_all_dft[iineq][spin_channel])

                    # since G(tau) has always 10001 values we are sampling +-10 values
                    # hard coded around beta/2, for beta=40 this corresponds to approx +-0.05
                    mesh_mid = len(G_time.data) // 2
                    samp = 10
                    gg = G_time.data[mesh_mid-samp:mesh_mid+samp]
                    gb2_averaged = np.mean(np.real(gg), axis=0)

                    # Z is not defined without Sigma, adding 1
                    Z_per_orbital.extend( [1.0] * G_time.target_shape[0]  )

                g_beta_half_per_orbital.extend(np.diag(gb2_averaged))
                g_beta_half_per_impurity += np.trace(gb2_averaged)

                # occupation per orbital
                den_mat = np.real(density_mat_dft[iineq][spin_channel])
                occupation_per_orbital.extend(np.diag(den_mat))
                occupation_per_impurity += np.trace(den_mat)

            # adding those values to the observable object
            observables['orb_gb2'][iineq][spin].append(np.array(g_beta_half_per_orbital))
            observables['imp_gb2'][iineq][spin].append(g_beta_half_per_impurity)
            observables['orb_occ'][iineq][spin].append(np.array(occupation_per_orbital))
            observables['imp_occ'][iineq][spin].append(occupation_per_impurity)
            observables['orb_Z'][iineq][spin].append(np.array(Z_per_orbital))

    # for it 0 we just subtract E_DC from E_DFT
    if general_params['calc_energies']:
        observables['E_tot'].append(dft_energy - sum([dc_per_imp[0] for dc_per_imp in observables['E_DC']]))
    else:
        observables['E_tot'].append('none')

    return observables


def add_dmft_observables(observables, general_params, solver_params, map_imp_solver, solver_type_per_imp, dft_energy, it, solvers, h_int,
                         previous_mu, sum_k, density_mat, shell_multiplicity, E_bandcorr):
    """
    calculates the observables for given Input, I decided to calculate the observables
    not adhoc since it should be done only once by the master_node

    Parameters
    ----------
    observables : observable arrays/dicts

    general_params : general parameters as a dict

    solver_params : solver parameters as a dict

    it : iteration counter

    solvers : Solver instances

    h_int : interaction hamiltonian

    previous_mu : dmft chemical potential for which the calculation was just done

    sum_k : SumK Object instances

    density_mat : DMFT occupations

    shell_multiplicity : degeneracy of impurities

    E_bandcorr : E_kin_dmft - E_kin_dft, either calculated man or from sum_k method if CSC

    Returns
    -------
    observables: list of dicts
    """

    # init energy values
    E_corr_en = 0.0
    # Read energy from OSZICAR
    dft_energy = 0.0 if dft_energy is None else dft_energy

    # now the normal output from each iteration
    observables['iteration'].append(it)
    observables['mu'].append(float(previous_mu))
    observables['E_bandcorr'].append(E_bandcorr)
    observables['E_dft'].append(dft_energy)

    if general_params['calc_energies']:
        mpi.report('\nCalculating interaction energies')
        for icrsh in range(sum_k.n_inequiv_shells):
            if (solver_type_per_imp[icrsh] in ['cthyb', 'hubbardI']
                    and solver_params[map_imp_solver[icrsh]]["measure_density_matrix"]):
                mpi.report(f'    Imp {icrsh}: from impurity density matrix')
                # Extract accumulated density matrix
                density_matrix = solvers[icrsh].density_matrix
                # Object containing eigensystem of the local Hamiltonian
                diag_local_ham = solvers[icrsh].h_loc_diagonalization
                E_int = trace_rho_op(density_matrix, h_int[icrsh], diag_local_ham)
            elif solver_type_per_imp[icrsh] == 'hartree':
                mpi.report(f'    Imp {icrsh}: from Hartree')
                E_int = solvers[icrsh].interaction_energy
            else:
                mpi.report(f'    Imp {icrsh}: from Migdal formula. '
                           'WARNING: less stable than measuring density matrix and using trace_rho_op!')
                # calc energy for given S and G
                # dmft interaction energy with E_int = 0.5 * Tr[Sigma * G]
                E_int = 0.5 * np.real((solvers[icrsh].G_freq * solvers[icrsh].Sigma_freq).total_density())

            observables['E_int'][icrsh].append(shell_multiplicity[icrsh]*E_int.real)
            E_corr_en += shell_multiplicity[icrsh] * (E_int.real - sum_k.dc_energ[sum_k.inequiv_to_corr[icrsh]])

    observables['E_corr_en'].append(E_corr_en)

    # calc total energy
    E_tot = dft_energy + E_bandcorr + E_corr_en
    observables['E_tot'].append(E_tot)

    for icrsh in range(sum_k.n_inequiv_shells):
        if general_params['dc']:
            observables['E_DC'][icrsh].append(shell_multiplicity[icrsh]*sum_k.dc_energ[sum_k.inequiv_to_corr[icrsh]])
        else:
            observables['E_DC'][icrsh].append(0.0)

        if solver_type_per_imp[icrsh] != 'ftps':
            if solvers[icrsh].G_time:
                G_time = solvers[icrsh].G_time
            else:
                G_time = solvers[icrsh].G_time_orig

        # Collect all occupation and G(beta/2) for spin up and down separately
        for spin in sum_k.spin_block_names[sum_k.SO]:
            g_beta_half_per_impurity = 0.0
            g_beta_half_per_orbital = []
            occupation_per_impurity = 0.0
            occupation_per_orbital = []
            Z_per_orbital = []

            # iterate over all spin channels and add the to up or down
            for spin_channel in sorted(sum_k.gf_struct_solver[icrsh].keys()):
                if not spin in spin_channel:
                    continue
                if solver_type_per_imp[icrsh] == 'ftps':
                    freq_mesh = np.array([w.value for w in solvers[icrsh].G_freq[spin_channel].mesh])
                    fermi_idx = abs(freq_mesh).argmin()
                    gb2_averaged = solvers[icrsh].G_freq[spin_channel].data[fermi_idx].imag

                    # Z is not defined without Sigma, adding 1
                    Z_per_orbital.extend( [1.0] * solvers[icrsh].G_freq[spin_channel].target_shape[0]  )
                else:
                    # G(beta/2)
                    # since G(tau) has always 10001 values we are sampling +-10 values
                    # hard coded around beta/2, for beta=40 this corresponds to approx +-0.05
                    mesh_mid = len(G_time[spin_channel].data) // 2
                    samp = 10
                    gg = G_time[spin_channel].data[mesh_mid-samp:mesh_mid+samp]
                    gb2_averaged = np.mean(np.real(gg), axis=0)

                    # get Z
                    Z_per_orbital.extend( calc_Z(Sigma= solvers[icrsh].Sigma_freq[spin_channel]) )

                g_beta_half_per_orbital.extend(np.diag(gb2_averaged))
                g_beta_half_per_impurity += np.trace(gb2_averaged)

                # occupation per orbital and impurity
                den_mat = np.real(density_mat[icrsh][spin_channel])
                occupation_per_orbital.extend(np.diag(den_mat))
                occupation_per_impurity += np.trace(den_mat)

            # adding those values to the observable object
            observables['orb_gb2'][icrsh][spin].append(np.array(g_beta_half_per_orbital))
            observables['imp_gb2'][icrsh][spin].append(g_beta_half_per_impurity)
            observables['orb_occ'][icrsh][spin].append(np.array(occupation_per_orbital))
            observables['imp_occ'][icrsh][spin].append(occupation_per_impurity)
            observables['orb_Z'][icrsh][spin].append(np.array(Z_per_orbital))

    return observables

def write_obs(observables, sum_k, general_params):
    """
    writes the last entries of the observable arrays to the files

    Parameters
    ----------
    observables : list of dicts
        observable arrays/dicts

    sum_k : SumK Object instances

    general_params : dict

    Returns
    -------
    nothing

    """

    n_orb = solver.get_n_orbitals(sum_k)

    for icrsh in range(sum_k.n_inequiv_shells):
        if general_params['magnetic'] and sum_k.SO == 0:
            for spin in ('up', 'down'):
                line = '{:3d} | '.format(observables['iteration'][-1])
                line += '{:10.5f} | '.format(observables['mu'][-1])

                if n_orb[icrsh][spin] == 1:
                    line += ' '*11
                for item in observables['orb_gb2'][icrsh][spin][-1]:
                    line += '{:10.5f}   '.format(item)
                line = line[:-3] + ' | '

                if n_orb[icrsh][spin] == 1:
                    line += ' '*11
                for item in observables['orb_occ'][icrsh][spin][-1]:
                    line += '{:10.5f}   '.format(item)
                line = line[:-3] + ' | '

                line += '{:17.5f}'.format(observables['imp_occ'][icrsh][spin][-1])

                if general_params['calc_energies']:
                    line += ' | {:10.5f}'.format(observables['E_tot'][-1])
                    line += ' | {:10.5f}'.format(observables['E_dft'][-1])
                    line += '   {:10.5f}'.format(observables['E_bandcorr'][-1])
                    line += '   {:10.5f}'.format(observables['E_int'][icrsh][-1])
                    line += '   {:10.5f}'.format(observables['E_DC'][icrsh][-1])

                file_name = '{}/observables_imp{}_{}.dat'.format(general_params['jobname'], icrsh, spin)
                with open(file_name, 'a') as obs_file:
                    obs_file.write(line + '\n')
        else:
            line = '{:3d} | '.format(observables['iteration'][-1])
            line += '{:10.5f} | '.format(observables['mu'][-1])

            # Adds spaces for header to fit in properly
            if n_orb[icrsh]['up'] == 1:
                line += ' '*11
            # Adds up the spin channels
            for iorb in range(n_orb[icrsh]['up']):
                val = np.sum([observables['orb_gb2'][icrsh][spin][-1][iorb] for spin in sum_k.spin_block_names[sum_k.SO]])
                line += '{:10.5f}   '.format(val)
            line = line[:-3] + ' | '

            # Adds spaces for header to fit in properly
            if n_orb[icrsh]['up'] == 1:
                line += ' '*11
            # Adds up the spin channels
            for iorb in range(n_orb[icrsh]['up']):
                val = np.sum([observables['orb_occ'][icrsh][spin][-1][iorb] for spin in sum_k.spin_block_names[sum_k.SO]])
                line += '{:10.5f}   '.format(val)
            line = line[:-3] + ' | '

            # Adds up the spin channels
            val = np.sum([observables['imp_occ'][icrsh][spin][-1] for spin in sum_k.spin_block_names[sum_k.SO]])
            line += '{:17.5f}'.format(val)

            if general_params['calc_energies']:
                line += ' | {:10.5f}'.format(observables['E_tot'][-1])
                line += ' | {:10.5f}'.format(observables['E_dft'][-1])
                line += '   {:10.5f}'.format(observables['E_bandcorr'][-1])
                line += '   {:10.5f}'.format(observables['E_int'][icrsh][-1])
                line += '   {:10.5f}'.format(observables['E_DC'][icrsh][-1])

            file_name = '{}/observables_imp{}.dat'.format(general_params['jobname'], icrsh)
            with open(file_name, 'a') as obs_file:
                obs_file.write(line + '\n')


def calc_dft_kin_en(general_params, sum_k, dft_mu):
    """
    Calculates the kinetic energy from DFT for target states

    Parameters
    ----------
    general_params : dict
        general parameters as a dict

    sum_k : SumK Object instances

    dft_mu: float
        DFT fermi energy


    Returns
    -------
    E_kin_dft: float
        kinetic energy from DFT

    """

    H_ks = sum_k.hopping
    num_kpts = sum_k.n_k
    E_kin = 0.0
    ikarray = np.array(list(range(sum_k.n_k)))
    for ik in mpi.slice_array(ikarray):
        nb = int(sum_k.n_orbitals[ik])
        # calculate lattice greens function need here to set sigma other n_iw is assumend to be 1025!
        # TODO: implement here version for FTPS!
        G_freq_lat = sum_k.lattice_gf(ik, with_Sigma=True, mu=dft_mu).copy()
        # # calculate G(beta) via the function density, which is the same as fourier trafo G(w) and taking G(b)
        G_freq_lat_beta = G_freq_lat.density()
        for spin in sum_k.spin_block_names[sum_k.SO]:
            E_kin += np.trace(np.dot(H_ks[ik, 0, :nb, :nb], G_freq_lat_beta[spin][:, :]))
    E_kin = np.real(E_kin)
    # collect data and put into E_kin_dft
    E_kin_dft = mpi.all_reduce(E_kin)
    mpi.barrier()
    # E_kin should be divided by the number of k-points
    E_kin_dft = E_kin_dft/num_kpts

    mpi.report(f'Kinetic energy contribution dft part: {E_kin_dft:.8f}')

    return E_kin_dft

def calc_bandcorr_man(general_params, sum_k, E_kin_dft):
    """
    Calculates the correlated kinetic energy from DMFT for target states
    and then determines the band correction energy

    Parameters
    ----------
    general_params : dict
        general parameters as a dict

    sum_k : SumK Object instances

    E_kin_dft: float
        kinetic energy from DFT


    Returns
    -------
    E_bandcorr: float
        band energy correction E_kin_dmft - E_kin_dft

    """
    E_kin_dmft = 0.0j
    E_kin = 0.0j
    H_ks = sum_k.hopping
    num_kpts = sum_k.n_k

    # kinetic energy from dmft lattice Greens functions
    ikarray = np.array(list(range(sum_k.n_k)))
    for ik in mpi.slice_array(ikarray):
        nb = int(sum_k.n_orbitals[ik])
        # calculate lattice greens function
        G_freq_lat = sum_k.lattice_gf(ik, with_Sigma=True, with_dc=True).copy()
        # calculate G(beta) via the function density, which is the same as fourier trafo G(w) and taking G(b)
        G_freq_lat_beta = G_freq_lat.density()
        for spin in sum_k.spin_block_names[sum_k.SO]:
            E_kin += np.trace(np.dot(H_ks[ik, 0, :nb, :nb], G_freq_lat_beta[spin][:, :]))
    E_kin = np.real(E_kin)

    # collect data and put into E_kin_dmft
    E_kin_dmft = mpi.all_reduce(E_kin)
    mpi.barrier()
    # E_kin should be divided by the number of k-points
    E_kin_dmft = E_kin_dmft/num_kpts

    if mpi.is_master_node():
        print('Kinetic energy contribution dmft part: '+str(E_kin_dmft))

    E_bandcorr = E_kin_dmft - E_kin_dft

    return E_bandcorr

def calc_Z(Sigma):
    """
    calculates the inverse mass enhancement from the impurity
    self-energy by a simple linear fit estimate:
    [ 1 - ((Im S_iw[n_iw0+1]-S_iw[n_iw0])/(iw[n_iw0+1]-iw[n_iw0])) ]^-1

    Parameters
    ----------
    Sigma: Gf on MeshImFreq
        self-energy on Matsubara mesh


    Returns
    -------
    orb_Z: 1d numpy array
        list of Z values per orbital in Sigma

    """
    orb_Z = []

    iw = [np.imag(n) for n in Sigma.mesh]
    # find smallest iw_n
    n_iw0 = int(0.5*len(iw))

    for orb in range(0,Sigma.target_shape[0]):
        Im_S_iw = Sigma[orb,orb].data.imag
        # simple extraction from linear fit to first two Matsubara freq of Sigma
        # assuming Fermi liquid like self energy
        Z = 1/(1 - (Im_S_iw[n_iw0+1]-Im_S_iw[n_iw0]) / (iw[n_iw0+1]-iw[n_iw0]) )
        orb_Z.append(abs(Z))

    return np.array(orb_Z)
