#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains unit tests for update_results_h5.py.
"""

import numpy as np
from solid_dmft.util.update_results_h5 import _update_results_per_iteration, _update_observables

def _generate_archive():
    dmft_results = {'it_2': {'DC_energ': 'test1',
                             'G0_iw0': 'test2',
                             'Gimp_iw_0': 'test3',
                             'chemical_potential': 12.5},
                    'it_5': {'chemical_potential': 12.9},
                    'last_iter': {},
                    'observables': {'mu': [11.9, 12.0, 12.8, 12.5, 12.1, 12.3, 12.9],
                                    'orb_gb2': [{s: [[-.05*i, -0.02/(i+1), -.1] for i in range(7)] for s in ('up', 'down')}]},
                    'iteration_count': 6,
                   }

    return {'DMFT_results': dmft_results}

def test_update_results_per_iteration():
    archive = _generate_archive()

    _update_results_per_iteration(archive)
    assert 'G0_iw0' not in archive['DMFT_results']['it_2']
    assert 'Gimp_iw_0' not in archive['DMFT_results']['it_2']
    assert archive['DMFT_results']['it_2']['G0_freq_0'] == 'test2'
    assert archive['DMFT_results']['it_2']['Gimp_freq_0'] == 'test3'
    assert np.isclose(archive['DMFT_results']['it_2']['chemical_potential_pre'], 12.8)
    assert np.isclose(archive['DMFT_results']['it_2']['chemical_potential_post'], 12.5)
    assert np.isclose(archive['DMFT_results']['it_5']['chemical_potential_pre'], 12.3)
    assert np.isclose(archive['DMFT_results']['it_5']['chemical_potential_post'], 12.9)

def test_update_observables():
    archive = _generate_archive()
    _update_observables(archive)
    assert len(archive['DMFT_results']['observables']['orb_Z']) == 1
    assert len(archive['DMFT_results']['observables']['orb_Z'][0]) == 2
    assert len(archive['DMFT_results']['observables']['orb_Z'][0]['down']) == 7
    assert len(archive['DMFT_results']['observables']['orb_Z'][0]['down'][1]) == 3
