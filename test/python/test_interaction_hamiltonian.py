#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains unit tests for interaction_hamiltonian
"""

import numpy as np
from solid_dmft.dmft_tools.interaction_hamiltonian import _load_crpa_interaction_matrix
from helper import Dummy

def test_load_crpa_interaction_matrix():
    sum_k = Dummy()
    sum_k.n_inequiv_shells = 2
    sum_k.n_corr_shells = 4
    sum_k.corr_shells = [{'dim': 5}, {'dim': 5}, {'dim': 5}, {'dim': 5}]
    sum_k.corr_to_inequiv = [0, 1, 1, 0]
    sum_k.inequiv_to_corr = [0, 1]
    sum_k.gf_struct_solver = [{'down_0' : 5, 'up_0' : 5},{'down_0' : 5, 'up_0' : 5}]

    general_params = {'crpa_code' : 'Vasp'}
    crpa_matrix = _load_crpa_interaction_matrix(sum_k, general_params, 'UIJKL')

    assert [c.shape for c in crpa_matrix] == [(5, 5, 5, 5), (5, 5, 5, 5)]

    # Warning: these assert were just taken from the current implementation
    #          and might be wrong if the implementation is buggy
    assert np.isclose(np.max(crpa_matrix),  2.5363308, rtol=0, atol=1e-6)
    assert np.isclose(np.min(crpa_matrix), -0.1894974, rtol=0, atol=1e-6)
    assert np.isclose(np.min(np.abs(crpa_matrix)), 7.76e-05, rtol=0, atol=1e-6)
