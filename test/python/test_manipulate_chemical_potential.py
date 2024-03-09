#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for manipulate_chemical_potential.py.
"""

import numpy as np
from solid_dmft.dmft_tools.manipulate_chemical_potential import _determine_band_edge, _mix_chemical_potential

import unittest

class test_manipulate_chemical_potential(unittest.TestCase):
    def test_determine_band_edge(self):
        # Creates spectral function as two Gaussian peaks at 3 and -2 eV
        mesh = np.linspace(-5, 5, 1000)
        center1 = 3
        center2 = -2
        width = 0.5
        spectral_function = np.exp(-(mesh-center1)**2/2/width**2) + np.exp(-(mesh-center2)**2/2/width**2)

        # The edge position is given by the point where the spectral function is only
        # edge_threshold * peak_height high, which we can calculate for Gaussian peaks
        edge_threshold = .2
        valence_edge = _determine_band_edge(mesh, spectral_function, .05, True, edge_threshold)
        exact_valence_edge = center1 - np.sqrt(-2*width**2*np.log(edge_threshold))
        conduction_edge = _determine_band_edge(mesh, spectral_function, .05, False, edge_threshold)
        exact_conduction_edge = center2 + np.sqrt(-2*width**2*np.log(edge_threshold))

        # Error in results comes from finite mesh
        assert np.abs(valence_edge - exact_valence_edge) < 1e-2
        assert np.abs(conduction_edge - exact_conduction_edge) < 1e-2

    def test_mix_chemical_potential(self):
        general_params = {'mu_mix_const': 0, 'mu_mix_per_occupation_offset': 1}
        density_tot = 16
        density_required = 16
        previous_mu = 0.0
        predicted_mu = 1.0

        new_mu = _mix_chemical_potential(general_params, density_tot, density_required,
                                         previous_mu, predicted_mu)
        assert np.isclose(new_mu, 0)

        general_params = {'mu_mix_const': 0, 'mu_mix_per_occupation_offset': 1}
        density_tot = 15.5
        density_required = 16
        previous_mu = 0.0
        predicted_mu = 1.0

        new_mu = _mix_chemical_potential(general_params, density_tot, density_required,
                                         previous_mu, predicted_mu)
        assert np.isclose(new_mu, .5)

        general_params = {'mu_mix_const': 1, 'mu_mix_per_occupation_offset': 0}
        density_tot = 12.34
        density_required = 16
        previous_mu = 0.0
        predicted_mu = 1.0

        new_mu = _mix_chemical_potential(general_params, density_tot, density_required,
                                         previous_mu, predicted_mu)
        assert np.isclose(new_mu, 1.)

        general_params = {'mu_mix_const': .3, 'mu_mix_per_occupation_offset': 1.}
        density_tot = 15.8
        density_required = 16
        previous_mu = 0.0
        predicted_mu = 1.0

        new_mu = _mix_chemical_potential(general_params, density_tot, density_required,
                                         previous_mu, predicted_mu)
        assert np.isclose(new_mu, .5)

if __name__ == '__main__':
    unittest.main()
