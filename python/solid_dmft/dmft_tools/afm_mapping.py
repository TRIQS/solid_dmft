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

import numpy as np
import triqs.utility.mpi as mpi

def determine(general_params, archive, n_inequiv_shells):
    """
    Determines the symmetries that are used in AFM calculations. These
    symmetries can then be used to copy the self-energies from one impurity to
    another by exchanging up/down channels for speedup and accuracy.
    """

    afm_mapping = None
    if mpi.is_master_node():
        # Reads mapping from h5 archive if it exists already from a previous run
        if 'afm_mapping' in archive['DMFT_input']:
            afm_mapping = archive['DMFT_input']['afm_mapping']
        elif len(general_params['magmom']) == n_inequiv_shells:
            # find equal or opposite spin imps, where we use the magmom array to
            # identity those with equal numbers or opposite
            # [copy Yes/False, from where, switch up/down channel]
            afm_mapping = [None] * n_inequiv_shells
            abs_moms = np.abs(general_params['magmom'])

            for icrsh in range(n_inequiv_shells):
                # if the moment was seen before ...
                previous_occurences = np.nonzero(np.isclose(abs_moms[:icrsh], abs_moms[icrsh]))[0]
                if previous_occurences.size > 0:
                    # find the source imp to copy from
                    source = np.min(previous_occurences)
                    # determine if we need to switch up and down channel
                    switch = np.isclose(general_params['magmom'][icrsh], -general_params['magmom'][source])

                    afm_mapping[icrsh] = [True, source, switch]
                else:
                    afm_mapping[icrsh] = [False, icrsh, False]


            print('AFM calculation selected, mapping self energies as follows:')
            print('imp  [copy sigma, source imp, switch up/down]')
            print('---------------------------------------------')
            for i, elem in enumerate(afm_mapping):
                print('{}: {}'.format(i, elem))
            print('')

            archive['DMFT_input']['afm_mapping'] = afm_mapping

        # if anything did not work set afm_order false
        else:
            print('WARNING: couldn\'t determine afm mapping. No mapping used.')
            general_params['afm_order'] = False

    general_params['afm_order'] = mpi.bcast(general_params['afm_order'])
    if general_params['afm_order']:
        general_params['afm_mapping'] = mpi.bcast(afm_mapping)

    return general_params


def apply(general_params, icrsh, gf_struct_solver, solvers):
    imp_source = general_params['afm_mapping'][icrsh][1]
    invert_spin = general_params['afm_mapping'][icrsh][2]
    mpi.report('\ncopying the self-energy for shell {} from shell {}'.format(icrsh, imp_source))
    mpi.report('inverting spin channels: '+str(invert_spin))

    if invert_spin:
        for spin_channel in gf_struct_solver.keys():
            if 'up' in spin_channel:
                target_channel = spin_channel.replace('up', 'down')
            else:
                target_channel = spin_channel.replace('down', 'up')

            solvers[icrsh].Sigma_freq[spin_channel] << solvers[imp_source].Sigma_freq[target_channel]
            solvers[icrsh].G_freq[spin_channel] << solvers[imp_source].G_freq[target_channel]
            solvers[icrsh].G_freq_unsym[spin_channel] << solvers[imp_source].G_freq_unsym[target_channel]
            solvers[icrsh].G0_freq[spin_channel] << solvers[imp_source].G0_freq[target_channel]
            solvers[icrsh].G_time[spin_channel] << solvers[imp_source].G_time[target_channel]

            if solvers[icrsh].solver_params['measure_pert_order']:
                if not hasattr(solvers[icrsh], 'perturbation_order'):
                    solvers[icrsh].perturbation_order = {}
                solvers[icrsh].perturbation_order[spin_channel] = solvers[imp_source].perturbation_order[target_channel]
                solvers[icrsh].perturbation_order_total = solvers[imp_source].perturbation_order_total

    else:
        solvers[icrsh].Sigma_freq << solvers[imp_source].Sigma_freq
        solvers[icrsh].G_freq << solvers[imp_source].G_freq
        solvers[icrsh].G_freq_unsym << solvers[imp_source].G_freq_unsym
        solvers[icrsh].G0_freq << solvers[imp_source].G0_freq
        solvers[icrsh].G_time << solvers[imp_source].G_time

        if solvers[icrsh].solver_params['measure_pert_order']:
            solvers[icrsh].perturbation_order = solvers[imp_source].perturbation_order
            solvers[icrsh].perturbation_order_total = solvers[imp_source].perturbation_order_total

    if solvers[icrsh].solver_params['measure_density_matrix']:
        solvers[icrsh].density_matrix = solvers[imp_source].density_matrix
        solvers[icrsh].h_loc_diagonalization = solvers[imp_source].h_loc_diagonalization

    if 'measure_chi' in solvers[icrsh].solver_params and solvers[icrsh].solver_params['measure_chi'] is not None:
        solvers[icrsh].O_time = solvers[imp_source].O_time

    return solvers
