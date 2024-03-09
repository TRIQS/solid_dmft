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
Contains all the functions related to setting the chemical potential in the
next iteration.
"""

import numpy as np

import triqs.utility.mpi as mpi
from triqs.gf import BlockGf, GfImFreq, GfImTime, Fourier, MeshImFreq
try:
    if mpi.is_master_node():
        from solid_dmft.postprocessing import maxent_gf_latt
    imported_maxent = True
except ImportError:
    imported_maxent = False

def _mix_chemical_potential(general_params, density_tot, density_required,
                            previous_mu, predicted_mu):
    """
    Mixes the previous chemical potential and the predicted potential with linear
    mixing:
    new_mu = factor * predicted_mu + (1-factor) * previous_mu, with
    factor = mu_mix_per_occupation_offset * |density_tot - density_required| + mu_mix_const
    under the constrain of 0 <= factor <= 1.

    Parameters
    ----------
    general_params : dict
        general parameters as a dict
    density_tot : float
        total occupation of the correlated system
    density_required : float
        required density for the impurity problem
    previous_mu : float
        the chemical potential from the previous iteration
    predicted_mu : float
        the chemical potential predicted by methods like the SumkDFT dichotomy

    Returns
    -------
    new_mu : float
        the chemical potential that results from the mixing

    """
    mu_mixing = general_params['mu_mix_per_occupation_offset'] * abs(density_tot - density_required)
    mu_mixing += general_params['mu_mix_const']
    mu_mixing = max(min(mu_mixing, 1), 0)
    new_mu = mu_mixing * predicted_mu + (1-mu_mixing) * previous_mu

    mpi.report('Mixing dichotomy mu with previous iteration by factor {:.3f}'.format(mu_mixing))
    mpi.report('New chemical potential: {:.3f}'.format(new_mu))
    return new_mu



def _initialize_lattice_gf(sum_k, general_params):
    """
    Creates lattice Green's function (GF) that is averaged over orbitals,
    blocks and spins. Returns lattice GF as input for an analytical
    continuation as well as G_lattice(tau=beta/2) (proxy for the spectral
    weight) and G_lattice(beta) (proxy for the total occupation).

    Parameters
    ----------
    sum_k : SumkDFT object
        Sumk object to generate the lattice GF from.
    general_params : dict
        general parameters as dict.

    Returns
    -------
    gf_lattice_iw : triqs.gf.BlockGf
        trace of the lattice GF over all blocks, orbitals and spins in
        Matsubara frequency.
    g_betahalf : complex
        the Fourier transform of gf_lattice_iw evaluated at tau=beta/2.
    occupation : complex
        the total density from gf_lattice_iw
    """

    # Initializes lattice GF to zero for each process
    mesh = sum_k.Sigma_imp[0].mesh
    trace_gf_latt = GfImFreq(mesh=mesh, data=np.zeros((len(mesh), 1, 1), dtype=complex))
    occupation = 0

    # Takes trace over orbitals and spins
    ikarray = np.arange(sum_k.n_k)
    for ik in mpi.slice_array(ikarray):
        gf_latt = sum_k.lattice_gf(ik)*sum_k.bz_weights[ik]
        trace_gf_latt.data[:] += np.trace(sum(g.data for _, g in gf_latt), axis1=1, axis2=2).reshape(-1, 1, 1)
        occupation += gf_latt.total_density()

    trace_gf_latt << mpi.all_reduce(trace_gf_latt)
    occupation = mpi.all_reduce(occupation)

    # Lattice GF as BlockGf, required for compatibility with MaxEnt functions
    gf_lattice_iw = BlockGf(name_list=['total'], block_list=[trace_gf_latt])

    # Fourier transforms the lattice GF
    gf_tau = GfImTime(beta=general_params['beta'], n_points=general_params['n_tau'], indices=[0])
    gf_tau << Fourier(gf_lattice_iw['total'])

    tau_mesh = np.array([float(m) for m in gf_tau.mesh])
    middle_index = np.argmin(np.abs(tau_mesh-general_params['beta']/2))
    samp = 10

    # Extracts G_latt(tau) at beta/2
    g_betahalf = np.mean(gf_tau.data[middle_index-samp:middle_index+samp, 0, 0])
    mpi.report('Lattice Gf: occupation = {:.5f}'.format(occupation))
    mpi.report('            G(beta/2)  = {:.5f}'.format(g_betahalf))

    return gf_lattice_iw, g_betahalf, occupation

def _determine_band_edge(mesh, spectral_function, spectral_func_threshold, valence_band, edge_threshold=.2):
    """
    Finds the band edge of a spectral function. This is done in two steps:
    starting from the Fermi energy, looks for the first peak
    (>spectral_func_threshold and local maximum on discrete grid). Then moves
    back towards Fermi energy until the spectral function is smaller than the
    fraction edge_threshold of the peak value.

    Parameters
    ----------
    mesh : numpy.ndarray of float
        the real frequencies grid.
    spectral_function : numpy.ndarray of float
        the values of the spectral function on the grid.
    spectral_func_threshold : float
        Threshold for spectral function to cross before looking for peaks.
    valence_band : bool
        Determines if looking for valence band (i.e. the upper band edge) or
        the conduction band (i.e. the lower band edge).
    edge_threshold : float
        Fraction of the peak value that defines the band edge value.

    Returns
    -------
    float
        The frequency value of the band edge.
    """
    # Determines direction to move away from Fermi energy to find band edge
    direction = 1 if valence_band else -1

    # Starts closest to the Fermi energy
    starting_index = np.argmin(np.abs(mesh))
    print('Starting at index {} with A(omega={:.3f})={:.3f}'.format(starting_index, mesh[starting_index], spectral_function[starting_index]))
    assert spectral_function[starting_index] < spectral_func_threshold

    # Finds peak
    peak_index = None
    for i in range(starting_index+direction, mesh.shape[0] if valence_band else -1, direction):
        # If A(omega) low, go further
        if spectral_function[i] < spectral_func_threshold:
            continue

        # If spectral function still increasing, go further
        if spectral_function[i-direction] < spectral_function[i]:
            continue

        peak_index = i-direction
        break

    assert peak_index is not None, 'Band peak not found. Check frequency range of MaxEnt'
    print('Peak at index {} with A(omega={:.3f})={:.3f}'.format(peak_index, mesh[peak_index], spectral_function[peak_index]))

    # Finds band edge
    edge_index = starting_index
    for i in range(peak_index-direction, starting_index-direction, -direction):
        # If above ratio edge_threshold of peak height, go further back to starting index
        if spectral_function[i] > edge_threshold * spectral_function[peak_index]:
            continue

        edge_index = i
        break

    print('Band edge at index {} with A(omega={:.3f})={:.3f}'.format(edge_index, mesh[edge_index], spectral_function[edge_index]))
    return mesh[edge_index]

def _set_mu_to_gap_middle_with_maxent(general_params, sum_k, gf_lattice_iw, archive=None):
    """
    Bundles running maxent on the total lattice GF, analyzing the spectral
    function and determining the new chemical potential.

    Parameters
    ----------
    general_params : dict
        general parameters as dict.
    sum_k : SumkDFT object
        SumkDFT object needed for original chemical potential and frequency
        range of MaxEnt continuation.
    gf_lattice_iw : BlockGf
        trace of the lattice GF over all blocks, orbitals and spins in
        Matsubara frequency.
    archive : HDFArchive, optional
        If given, writes spectral function (i.e. MaxEnt result) to archive.

    Returns
    -------
    float
        new chemical potential located in the middle of the gap from MaxEnt.
        None if not master node or if something went wrong.
    """


    if not mpi.is_master_node():
        return None

    if not imported_maxent:
        mpi.report('WARNING: cannot find gap with MaxEnt, MaxEnt not found')
        return None

    # Runs MaxEnt using the Chi2Curvature analyzer
    maxent_results, mesh = maxent_gf_latt._run_maxent(gf_lattice_iw, sum_k, .02, None, None, 200, 30)
    mesh = np.array(mesh)
    spectral_function = maxent_results['total'].get_A_out('Chi2CurvatureAnalyzer')

    # Writes spectral function to archive
    if archive is not None:
        unpacked_results = maxent_gf_latt._unpack_maxent_results(maxent_results, mesh)
        archive['DMFT_results/last_iter']['Alatt_w'] = unpacked_results

    # Checks if spectral function at Fermi energy below threshold
    spectral_func_threshold = general_params['beta']/np.pi * general_params['mu_gap_gb2_threshold']
    if spectral_function[np.argmin(np.abs(mesh))] > spectral_func_threshold:
        mpi.report('WARNING: cannot find gap with MaxEnt, spectral function not gapped at Fermi energy')
        return None

    # Determines band edges for conduction and valence band
    edge_threshold = 0.2
    conduction_edge = _determine_band_edge(mesh, spectral_function, spectral_func_threshold, False, edge_threshold)
    valence_edge = _determine_band_edge(mesh, spectral_function, spectral_func_threshold, True, edge_threshold)

    return sum_k.chemical_potential + (valence_edge + conduction_edge) / 2

def set_initial_mu(general_params, sum_k, iteration_offset, archive, broadening):
    """
    Handles the different ways of setting the initial chemical potential mu:
    * Chemical potential set to fixed value: uses this value

    * New calculation: determines mu from dichotomy method

    * Resuming calculation and chemical potential not updated this iteration:
       loads calculation before previous iteration.

    * Resuming calculation and chemical potential is updated:
        checks if the system is gapped and potentially run MaxEnt to find gap
        middle. Otherwise, gets mu from dichotomy and applies mu mixing to result.


    Parameters
    ----------
    general_params : dict
        general parameters as dict.
    sum_k : SumkDFT object
        contains system information necessary to determine the initial mu.
    iteration_offset : int
        the number of iterations executed in previous calculations.
    archive : HDFArchive
        needed to potentially load previous results and write MaxEnt results to.

    Returns
    -------
    sum_k : SumkDFT object
        the altered SumkDFT object with the initial mu set correctly.
    """

    # Uses fixed_mu_value as chemical potential if parameter is given
    if general_params['fixed_mu_value'] is not None:
        sum_k.set_mu(general_params['fixed_mu_value'])
        mpi.report('+++ Keeping the chemical potential fixed at {:.3f} eV +++'.format(general_params['fixed_mu_value']))
        return sum_k

    # In first iteration, determines mu and returns
    if iteration_offset == 0:
        sum_k.calc_mu(precision=general_params['prec_mu'], method=general_params['calc_mu_method'],
                      broadening=broadening)
        return sum_k

    # If continuing calculation and not updating mu, loads sold value
    if iteration_offset % general_params['mu_update_freq'] != 0:
        if mpi.is_master_node():
            sum_k.chemical_potential = archive['DMFT_results/last_iter/chemical_potential_pre']
        sum_k.chemical_potential = mpi.bcast(sum_k.chemical_potential)
        mpi.report('Chemical potential not updated this step, '
                   + 'reusing loaded one of {:.3f} eV'.format(sum_k.chemical_potential))
        return sum_k

    # If continuing calculation and updating mu, reads in occupation and
    # chemical_potential_pre from the last run
    previous_mu = None
    if mpi.is_master_node():
        previous_mu = archive['DMFT_results/last_iter/chemical_potential_pre']
    previous_mu = mpi.bcast(previous_mu)

    # Runs maxent if spectral weight too low and occupation is close to desired one
    if isinstance(sum_k.mesh, MeshImFreq) and general_params['mu_gap_gb2_threshold'] is not None:
        sum_k.chemical_potential = previous_mu
        gf_lattice_iw, g_betahalf, occupation = _initialize_lattice_gf(sum_k, general_params)
        fulfills_occupation_crit = (general_params['mu_gap_occ_deviation'] is None
                                    or np.abs(occupation - sum_k.density_required) < general_params['mu_gap_occ_deviation'])

        if -np.real(g_betahalf) < general_params['mu_gap_gb2_threshold'] and fulfills_occupation_crit:
            new_mu = _set_mu_to_gap_middle_with_maxent(general_params, sum_k, gf_lattice_iw, archive)
            new_mu = mpi.bcast(new_mu)
            if new_mu is not None:
                sum_k.chemical_potential = new_mu
                mpi.report('New chemical potential in the gap: {:.3f} eV'.format(new_mu))
                return sum_k
    # Calculates occupation for mu mixing below
    elif np.isclose(general_params['mu_mix_per_occupation_offset'], 0):
        occupation = 0 # The occupation does not matter in this case
    else:
        _, _, occupation = _initialize_lattice_gf(sum_k, general_params)

    # If system not gapped, gets chemical potential from dichotomy method
    sum_k.calc_mu(precision=general_params['prec_mu'], method=general_params['calc_mu_method'],
                  broadening=broadening)

    # Applies mu mixing to dichotomy result
    sum_k.chemical_potential = _mix_chemical_potential(general_params, occupation,
                                                       sum_k.density_required,
                                                       previous_mu, sum_k.chemical_potential)

    return sum_k

def update_mu(general_params, sum_k, it, archive, broadening):
    """
    Handles the different ways of updating the chemical potential mu:
    * Chemical potential set to fixed value: uses this value

    * Chemical potential not updated this iteration: nothing happens.

    * Chemical potential is updated: checks if the system is gapped and
        potentially run MaxEnt to find gap middle. Otherwise, gets mu from
        dichotomy and applies mu mixing to result.

    Parameters
    ----------
    general_params : dict
        general parameters as dict.
    sum_k : SumkDFT object
        contains system information necessary to update mu.
    it : int
        the number of the current iteration.
    archive : HDFArchive
        needed to potentially write MaxEnt results to.

    Returns
    -------
    sum_k : SumkDFT object
        the altered SumkDFT object with the updated mu.
    """

    # Uses fixed_mu_value as chemical potential if parameter is given
    if general_params['fixed_mu_value'] is not None:
        sum_k.set_mu(general_params['fixed_mu_value'])
        mpi.report('+++ Keeping the chemical potential fixed at {:.3f} eV +++'.format(general_params['fixed_mu_value']))
        return sum_k

    # If mu won't be updated this step, don't update it...
    if it % general_params['mu_update_freq'] != 0:
        mpi.report('Chemical potential not updated this step, '
                   + 'reusing previous one of {:.3f} eV'.format(sum_k.chemical_potential))
        return sum_k

    # Runs maxent if spectral weight too low and occupation is close to desired one
    # TODO: which solvers work?
    if isinstance(sum_k.mesh, MeshImFreq) and general_params['mu_gap_gb2_threshold'] is not None:
        gf_lattice_iw, g_betahalf, occupation = _initialize_lattice_gf(sum_k, general_params)
        fulfills_occupation_crit = (general_params['mu_gap_occ_deviation'] is None
                                    or np.abs(occupation - sum_k.density_required) < general_params['mu_gap_occ_deviation'])

        if -np.real(g_betahalf) < general_params['mu_gap_gb2_threshold'] and fulfills_occupation_crit:
            new_mu = _set_mu_to_gap_middle_with_maxent(general_params, sum_k, gf_lattice_iw, archive)
            new_mu = mpi.bcast(new_mu)
            if new_mu is not None:
                sum_k.chemical_potential = new_mu
                mpi.report('New chemical potential in the gap: {:.3f} eV'.format(new_mu))
                return sum_k
    # Calculates occupation for mu mixing below
    elif np.isclose(general_params['mu_mix_per_occupation_offset'], 0):
        occupation = 0 # The occupation does not matter in this case
    else:
        _, _, occupation = _initialize_lattice_gf(sum_k, general_params)

    # If system not gapped, gets chemical potential from dichotomy method
    previous_mu = sum_k.chemical_potential
    sum_k.calc_mu(precision=general_params['prec_mu'], method=general_params['calc_mu_method'],
                  broadening=broadening)

    # Applies mu mixing to dichotomy result
    sum_k.chemical_potential = _mix_chemical_potential(general_params, occupation,
                                                       sum_k.density_required,
                                                       previous_mu, sum_k.chemical_potential)
    mpi.barrier()
    return sum_k
