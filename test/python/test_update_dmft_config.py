#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains unit tests for update_dmft_config.py.
"""

from configparser import ConfigParser

from solid_dmft.util.update_dmft_config import _update_section_names, _update_params

CONFIG_FILE_1 = u'''
[general]
seedname = fancy_system
jobname = out_DMFT_fancy
enforce_off_diag = True
block_threshold= 0.001
set_rot = none
solver_type = cthyb

prec_mu = -0.001

h_int_type = density_density
U = 2.0
J = 2.0

n_iter_dmft = 4

dc_type = 0
dc = True
dc_dmft = False

calc_energies = True
sigma_mix = 0.7

h5_save_freq = 5

n_iter_dmft_per = 5

[solver]
imag_threshold = 1e-5
measure_G_l = True
n_LegCoeff = 40

[advanced_parameters]
dc_factor = 0.6

[dft]
'''

def _print_config(config):
    print('Section headers:', config.sections())
    for section in config:
        print('-'*10, section)
        print([(k, v) for k,v in config[section].items()])

def test_update_legacy_section_headers():
    config = ConfigParser()
    config.read_string(CONFIG_FILE_1)

    config = _update_section_names(config)
    assert all([s in ('general', 'solver', 'advanced', 'dft') for s in config.sections()])
    assert [(k, v) for k, v in config['advanced'].items()] == [('dc_factor', '0.6')]

def test_update_params():
    config = ConfigParser()
    config.read_string(CONFIG_FILE_1)

    config = _update_section_names(config)
    config = _update_params(config)

    assert config['general']['h_int_type'] == 'density_density'
    assert config['general']['solver_type'] == 'cthyb'
    assert config['general']['n_l'] == '40'
    assert 'n_LegCoeff' not in config['solver']
