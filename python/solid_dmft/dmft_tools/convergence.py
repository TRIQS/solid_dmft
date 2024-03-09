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
'''
contain helper functions to check convergence
'''
# system
import os.path
import numpy as np

# triqs
from triqs.gf import MeshImFreq, MeshImTime, MeshReFreq, BlockGf
from solid_dmft.dmft_tools import solver

def _generate_header(general_params, sum_k):
    """
    Generates the headers that are used in write_header_to_file.
    Returns a dict with {file_name: header_string}
    """

    n_orb = solver.get_n_orbitals(sum_k)

    header_energy_mask = '| {:>11} '
    header_energy = header_energy_mask.format('δE_tot')

    headers = {}
    for iineq in range(sum_k.n_inequiv_shells):
        number_spaces = max(13*n_orb[iineq]['up']-1, 21)
        header_basic_mask = '{{:>3}} | {{:>11}} | {{:>{0}}} | {{:>11}} | {{:>11}} | {{:>11}} | {{:>11}} '.format(number_spaces)

        file_name = 'conv_imp{}.dat'.format(iineq)
        headers[file_name] = header_basic_mask.format('it', 'δμ','δocc orb','δimp occ','δGimp', 'δG0','δΣ')

        if general_params['calc_energies']:
            headers[file_name] += header_energy


    return headers

def write_conv(conv_obs, sum_k, general_params):
    """
    writes the last entries of the conv arrays to the files

    Parameters
    ----------
    conv_obs : list of dicts
        convergence observable arrays/dicts

    sum_k : SumK Object instances

    general_params : dict

    __Returns:__

    nothing

    """

    n_orb = solver.get_n_orbitals(sum_k)

    for icrsh in range(sum_k.n_inequiv_shells):
        line = '{:3d} | '.format(conv_obs['iteration'][-1])
        line += '{:10.5e} | '.format(conv_obs['d_mu'][-1])

        # Adds spaces for header to fit in properly
        if n_orb[icrsh]['up'] == 1:
            line += ' '*11
        # Adds up the spin channels
        for iorb in range(n_orb[icrsh]['up']):
            line += '{:10.5e}   '.format(conv_obs['d_orb_occ'][icrsh][-1][iorb])
        line = line[:-3] + ' | '

        # imp occupation change
        line += '{:10.5e} | '.format(conv_obs['d_imp_occ'][icrsh][-1])
        # Gimp change
        line += '{:10.5e} | '.format(conv_obs['d_Gimp'][icrsh][-1])
        # G0 change
        line += '{:10.5e} | '.format(conv_obs['d_G0'][icrsh][-1])
        # Σ change
        line += '{:10.5e}'.format(conv_obs['d_Sigma'][icrsh][-1])

        if general_params['calc_energies']:
            line += ' | {:10.5e}'.format(conv_obs['d_Etot'][-1])

        file_name = '{}/conv_imp{}.dat'.format(general_params['jobname'], icrsh)
        with open(file_name, 'a') as obs_file:
            obs_file.write(line + '\n')

def max_G_diff(G1, G2, norm_temp = True):
    """
    calculates difference between two block Gfs
    uses numpy linalg norm on the last two indices first
    and then the norm along the mesh axis. The result is divided
    by sqrt(beta) for MeshImFreq and by sqrt(beta/#taupoints) for
    MeshImTime.

    1/ (2* sqrt(beta)) sqrt( sum_n sum_ij [abs(G1 - G2)_ij(w_n)]^2 )

    this is only done for MeshImFreq Gf objects, for all other
    meshes the weights are set to 1

    Parameters
    ----------
    G1 : Gf or BlockGf to compare

    G2 : Gf or BlockGf to compare

    norm_temp: bool, default = True
       divide by an additional sqrt(beta) to account for temperature scaling
       only correct for uniformly distributed error.

    __Returns:__

    diff : float
           difference between the two Gfs
    """

    if isinstance(G1, BlockGf):
        diff = 0.0
        for block, gf in G1:
            diff += max_G_diff(G1[block], G2[block], norm_temp)
        return diff

    assert G1.mesh == G2.mesh, 'mesh of two input Gfs does not match'
    assert G1.target_shape == G2.target_shape, 'can only compare Gfs with same shape'

    # subtract largest real value to make sure that G1-G2 falls off to 0
    if type(G1.mesh) is MeshImFreq:
        offset = np.diag(np.diag(G1.data[-1,:,:].real - G2.data[-1,:,:].real))
    else:
        offset = 0.0

    #  calculate norm over all axis but the first one which are frequencies
    norm_grid = abs(np.linalg.norm(G1.data - G2.data - offset, axis=tuple(range(1, G1.data.ndim))))
    # now calculate Frobenius norm over grid points
    norm = np.linalg.norm(norm_grid, axis=0)

    if type(G1.mesh) is MeshImFreq:
        norm = np.linalg.norm(norm_grid, axis=0) / np.sqrt(G1.mesh.beta)
    elif type(G1.mesh) is MeshImTime:
        norm = np.linalg.norm(norm_grid, axis=0) * np.sqrt(G1.mesh.beta/len(G1.mesh))
    elif type(G1.mesh) is MeshReFreq:
        norm = np.linalg.norm(norm_grid, axis=0) / np.sqrt(len(G1.mesh))
    else:
        raise ValueError('MeshReTime is not implemented')

    if type(G1.mesh) in (MeshImFreq, MeshImTime) and norm_temp:
        norm = norm / np.sqrt(G1.mesh.beta)

    return norm

def prep_conv_obs(h5_archive):
    """
    prepares the conv arrays and files for the DMFT calculation

    Parameters
    ----------
    h5_archive: hdf archive instance
        hdf archive for calculation

    __Returns:__
    conv_obs : dict
        conv array for calculation
    """

    # determine number of impurities
    n_inequiv_shells = h5_archive['dft_input']['n_inequiv_shells']

    # check for previous iterations
    conv_prev = []
    if 'convergence_obs' in h5_archive['DMFT_results']:
        conv_prev = h5_archive['DMFT_results']['convergence_obs']

    # prepare observable dicts
    if len(conv_prev) > 0:
        conv_obs = conv_prev
    else:
        conv_obs = dict()
        conv_obs['iteration'] = []
        conv_obs['d_mu'] = []
        conv_obs['d_Etot'] = []
        conv_obs['d_orb_occ'] = [[] for i in range(n_inequiv_shells)]
        conv_obs['d_imp_occ'] = [[] for i in range(n_inequiv_shells)]

        conv_obs['d_Gimp'] = [[] for i in range(n_inequiv_shells)]
        conv_obs['d_G0'] = [[] for i in range(n_inequiv_shells)]
        conv_obs['d_Sigma'] = [[] for i in range(n_inequiv_shells)]

    return conv_obs

def prep_conv_file(general_params, sum_k):
    """
    Writes the header to the conv files

    Parameters
    ----------
    general_params : dict
        general parameters as a dict
    n_inequiv_shells : int
        number of impurities for calculations


    __Returns:__
    nothing
    """

    headers = _generate_header(general_params, sum_k)

    for file_name, header in headers.items():
        path = os.path.join(general_params['jobname'], file_name)
        with open(path, 'w') as conv_file:
            conv_file.write(header + '\n')


def calc_convergence_quantities(sum_k, general_params, conv_obs, observables,
                                     solvers, G0_old, G_loc_all, Sigma_freq_previous):
    """
    Calculations convergence quantities, i.e. the difference in observables
    between the last and second to last iteration.

    Parameters
    ----------
    sum_k : SumK Object instances

    general_params : dict
        general parameters as a dict

    conv_obs : list of dicts
        convergence observable arrays

    observables : list of dicts
        observable arrays

    solvers : solver objects

    G0_old : list of block Gf object
            last G0_freq

    G_loc_all : list of block Gf objects
            G_loc extracted from before imp solver

    Sigma_freq_previous : list of block Gf objects
            previous impurity sigma to compare with

    Returns
    -------
    conv_obs : list of dicts
        updated convergence observable arrays

    """

    conv_obs['iteration'].append(observables['iteration'][-1])
    conv_obs['d_mu'].append(abs(observables['mu'][-1] - observables['mu'][-2] ))
    for icrsh in range(sum_k.n_inequiv_shells):
        if not sum_k.corr_shells[icrsh]['SO']:
            # difference in imp occupation
            conv_obs['d_imp_occ'][icrsh].append(abs((observables['imp_occ'][icrsh]['up'][-1]+
                                                 observables['imp_occ'][icrsh]['down'][-1])-
                                                (observables['imp_occ'][icrsh]['up'][-2]+
                                                 observables['imp_occ'][icrsh]['down'][-2])))
            # difference in orb occ spin absolute
            conv_obs['d_orb_occ'][icrsh].append(abs(observables['orb_occ'][icrsh]['up'][-1]-
                                                    observables['orb_occ'][icrsh]['up'][-2])+
                                                abs(observables['orb_occ'][icrsh]['down'][-1]-
                                                    observables['orb_occ'][icrsh]['down'][-2]))
        else:
            conv_obs['d_imp_occ'][icrsh].append(abs(observables['imp_occ'][icrsh]['ud'][-1]-
                                                    observables['imp_occ'][icrsh]['ud'][-2]))
            conv_obs['d_orb_occ'][icrsh].append(abs(observables['orb_occ'][icrsh]['ud'][-1]+
                                                    observables['imp_occ'][icrsh]['ud'][-2]))

        conv_obs['d_Gimp'][icrsh].append(max_G_diff(solvers[icrsh].G_freq, G_loc_all[icrsh]))
        conv_obs['d_G0'][icrsh].append(max_G_diff(solvers[icrsh].G0_freq, G0_old[icrsh]))
        conv_obs['d_Sigma'][icrsh].append(max_G_diff(solvers[icrsh].Sigma_freq, Sigma_freq_previous[icrsh]))

    if general_params['calc_energies']:
        conv_obs['d_Etot'].append(abs(observables['E_tot'][-1]-observables['E_tot'][-2]))

    return conv_obs


def check_convergence(n_inequiv_shells, general_params, conv_obs):
    """
    check last iteration for convergence

    Parameters
    ----------
    n_inequiv_shells : int
        Number of inequivalent shells as saved in SumkDFT object

    general_params : dict
        general parameters as a dict

    conv_obs : list of dicts
        convergence observable arrays

    Returns
    -------
    is_converged : bool
        true if desired accuracy is reached. None if no convergence criterion
        is set

    """

    # If no convergence criterion is set, convergence is undefined and returns None
    if (general_params['occ_conv_crit'] <= 0.0 and general_params['gimp_conv_crit'] <= 0.0
            and general_params['g0_conv_crit'] <= 0.0 and general_params['sigma_conv_crit'] <= 0.0):
        return None

    # Checks convergence criteria
    for icrsh in range(n_inequiv_shells):
        # Checks imp occ
        if (conv_obs['d_imp_occ'][icrsh][-1] > general_params['occ_conv_crit']
                and general_params['occ_conv_crit'] > 0.0):
            return False

        # Checks Gimp
        if (conv_obs['d_Gimp'][icrsh][-1] > general_params['gimp_conv_crit']
                and general_params['gimp_conv_crit'] > 0.0):
            return False

        # Checks G0
        if (conv_obs['d_G0'][icrsh][-1] > general_params['g0_conv_crit']
                and general_params['g0_conv_crit'] > 0.0):
            return False

        # Checks Sigma
        if (conv_obs['d_Sigma'][icrsh][-1] > general_params['sigma_conv_crit']
                and general_params['sigma_conv_crit'] > 0.0):
            return False

    return True
