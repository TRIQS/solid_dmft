#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains unit tests for dmft_cycle
"""

import sys
sys.path.append('./src')

from solid_dmft.dmft_tools.afm_mapping import determine
from helper import are_iterables_equal

def test_determine_afm_mapping():
    general_params = {'magmom': [+1, -1, +1, -1], 'afm_order': True}
    archive = {'DMFT_input': {}}
    n_inequiv_shells = 4

    expected_general_params = general_params.copy()
    #                                            copy, source, switch
    expected_general_params['afm_mapping'] = [[False, 0, False], [True, 0, True],
                                                  [True, 0, False], [True, 0, True]]

    general_params = determine(general_params, archive, n_inequiv_shells)

    assert are_iterables_equal(general_params, expected_general_params)

    general_params = {'magmom': [+1, -1, +2, +2], 'afm_order': True}
    archive = {'DMFT_input': {}}
    n_inequiv_shells = 4

    expected_general_params = general_params.copy()
    #                                            copy, source, switch
    expected_general_params['afm_mapping'] = [[False, 0, False], [True, 0, True],
                                                  [False, 2, False], [True, 2, False]]

    general_params = determine(general_params, archive, n_inequiv_shells)

    assert are_iterables_equal(general_params, expected_general_params)

    # Reading in the afm_mapping from the archive
    general_params = {'magmom': [+1, -1, +2], 'afm_order': True}
    archive = {'DMFT_input': {'afm_mapping': [[False, 0, False], [False, 1, False],
                                              [False, 2, False]]}}
    n_inequiv_shells = 3

    expected_general_params = general_params.copy()
    #                                            copy, source, switch
    expected_general_params['afm_mapping'] = [[False, 0, False], [False, 1, False],
                                                  [False, 2, False]]

    general_params = determine(general_params, archive, n_inequiv_shells)

    assert are_iterables_equal(general_params, expected_general_params)

    general_params = {'magmom': [+1, -1, +2, +2], 'afm_order': True}
    archive = {'DMFT_input': {}}
    n_inequiv_shells = 3

    expected_general_params = general_params.copy()
    #                                            copy, source, switch
    expected_general_params['afm_order'] = False

    general_params = determine(general_params, archive, n_inequiv_shells)

    assert are_iterables_equal(general_params, expected_general_params)
