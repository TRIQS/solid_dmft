#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains unit tests for read_config.py
"""


from configparser import ConfigParser

from solid_dmft.read_config import (_config_find_default_section_entries,
                             _config_add_empty_sections, _config_remove_unused_sections,
                             _convert_params, _find_nonexistent_params,
                             _apply_default_values, _check_if_params_used,
                             _checks_validity_criterion)

from helper import are_iterables_equal

CONFIG_FILE_1 = u'''
[general]
seedname = fancy_system
jobname = out_DMFT_fancy
enforce_off_diag = True
block_threshold= 0.001
set_rot = none
solver_type = cthyb

prec_mu = -0.001

h_int_type =    density_density
U = 2.0
J = 2.0

n_iter_dmft = 4

dc_type = 0
dc = True
dc_dmft = False

calc_energies = True
g0_mix = 0.7

h5_save_freq = 5

n_iter_dmft_per = 5

[advanced]
nonexistent = 2

[weird_params]
test = 2
system_size = 1e23
'''

def test_config_file_1():
    config = ConfigParser()
    config.read_string(CONFIG_FILE_1)

    # Checks if default section is empty
    config_default_entries = _config_find_default_section_entries(config)
    assert config_default_entries == []

    # Adds empty sections if they don't exist
    config = _config_add_empty_sections(config)

    # Removes unused sections and prints a warning
    config, unused_sections = _config_remove_unused_sections(config)
    assert unused_sections == ['weird_params']

    parameters = _convert_params(config)

    nonexistent_params = _find_nonexistent_params(config)
    assert nonexistent_params == {'dft': [], 'general': [], 'advanced': ['nonexistent'], 'solver': []}

    parameters, default_values_used = _apply_default_values(parameters)

    parameters, unnecessary_params, missing_required_params = _check_if_params_used(parameters, default_values_used)
    assert unnecessary_params == {'dft': [], 'general': ['n_iter_dmft_per'], 'advanced': [], 'solver': []}
    assert missing_required_params == {'dft': [], 'general': ['beta'], 'advanced': [],
                                           'solver': ['length_cycle', 'n_warmup_cycles', 'n_cycles_tot']}

    invalid_params = _checks_validity_criterion(parameters)
    assert invalid_params == {'dft': [], 'general': ['prec_mu'], 'advanced': [], 'solver': []}

    assert are_iterables_equal(parameters, {'dft': {},
                                            'general': {'magnetic': False, 'fixed_mu_value': 'none',
                                                        'mu_update_freq': 1,
                                                        'measure_chi': 'none', 'block_threshold': 0.001,
                                                        'set_rot': u'none', 'prec_mu': -0.001,
                                                        'dft_mu': u'none', 'solver_type': 'cthyb',
                                                        'mu_mix_const': 1., 'mu_mix_per_occupation_offset': 0.,
                                                        'oneshot_postproc_gamma_file': False,
                                                        'csc': False, 'enforce_off_diag': True,
                                                        'dc_dmft': False, 'diag_delta': False,
                                                        'occ_conv_crit': -1,'g0_conv_crit': -1,'gimp_conv_crit': -1,'sigma_conv_crit': -1,
                                                        'seedname': [u'fancy_system'],
                                                        'J': [2.0], 'h5_save_freq': 5, 'ratio_F4_F2' : [u'none'],
                                                        'dc': True, 'jobname': [u'out_DMFT_fancy'],
                                                        'n_iter_dmft': 4, 'U': [2.0],
                                                        'energy_shift_orbitals': 'none', 'n_tau': 10001,
                                                        'measure_chi_insertions': 100, 'h_field': 0.0,
                                                        'calc_energies': True, 'g0_mix': 0.7, 'sigma_mix': 1.0,
                                                        'g0_mix_type' : 'linear',
                                                        'dc_type': 0, 'load_sigma': False, 'n_iw': 1025,
                                                        'noise_level_initial_sigma': 0.,
                                                        'h_int_type': 'density_density',
                                                        'mu_gap_gb2_threshold': 'none'},
                                            'advanced': {'dc_fixed_value': 'none', 'dc_fixed_occ': 'none',
                                                         'dc_nominal': False,
                                                         'dc_factor': 'none', 'dc_J': [2.0], 'dc_U': [2.0],
                                                         'map_solver_struct': 'none'},
                                            'solver': {'move_double': True, 'measure_G_l': False,
                                                       'move_shift': False, 'store_solver': False,
                                                       'measure_pert_order': False,
                                                       'measure_density_matrix': False, 'perform_tail_fit': False,
                                                       'legendre_fit': False, 'delta_interface' : False,
                                                       'off_diag_threshold' : 0.0}}
                              )



CONFIG_FILE_2 = u'''
[general]
seedname = orbital_model
jobname = out_60M_afm
enforce_off_diag = True
block_threshold= 0.001
set_rot = none
solver_type = cthyb

prec_mu = 0.001

h_int_type = kanamori
U = 5.5
J = 1.0
n_l = 35

beta = 40
energy_shift_orbitals = 2*%(U)s, 0, 3.8*%(J)s + %(U)s, 2.1

n_iter_dmft = 6

path_to_sigma = orbital_2site_model_Sigma_eg_swapped.h5

dc_type = 0
dc = True
dc_dmft = False

magnetic = True
magmom = 1, -1
afm_order = True

calc_energies = False
sigma_mix = 0.6

h5_save_freq = 2

[solver]
length_cycle = 140
n_warmup_cycles = 10000
n_cycles_tot = 60e+6
imag_threshold = 1e-5
measure_G_l = True
delta_interface = True
off_diag_threshold = 0.1

measure_density_matrix = False

perform_tail_fit = True

[advanced]
map_solver_struct = {('ud_0', 0): ('up_0', 0)}
'''

def test_config_file_2():
    config = ConfigParser()
    config.read_string(CONFIG_FILE_2)

    # Checks if default section is empty
    config_default_entries = _config_find_default_section_entries(config)
    assert config_default_entries == []

    # Adds empty sections if they don't exist
    config = _config_add_empty_sections(config)

    # Removes unused sections and prints a warning
    config, unused_sections = _config_remove_unused_sections(config)
    assert unused_sections == []

    parameters = _convert_params(config)

    nonexistent_params = _find_nonexistent_params(config)
    assert nonexistent_params == {'dft': [], 'general': [], 'advanced': [], 'solver': []}

    parameters, default_values_used = _apply_default_values(parameters)

    parameters, unnecessary_params, missing_required_params = _check_if_params_used(parameters, default_values_used)
    assert unnecessary_params == {'dft': [], 'general': ['path_to_sigma'],
                                      'advanced': [], 'solver': ['perform_tail_fit']}
    assert missing_required_params == {'dft': [], 'general': [], 'advanced': [], 'solver': []}


    invalid_params = _checks_validity_criterion(parameters)
    assert invalid_params == {'dft': [], 'general': [], 'advanced': [], 'solver': []}

    print(parameters)
    assert are_iterables_equal(parameters, {'dft': {},
                                            'general': {'afm_order': True, 'magnetic': True,
                                                        'measure_chi': 'none',
                                                        'block_threshold': 0.001, 'set_rot': u'none',
                                                        'prec_mu': 0.001, 'dft_mu': u'none',
                                                        'mu_mix_const': 1., 'mu_mix_per_occupation_offset': 0.,
                                                        'oneshot_postproc_gamma_file': False, 'csc': False,
                                                        'enforce_off_diag': True, 'fixed_mu_value': 'none',
                                                        'mu_update_freq': 1, 'solver_type': 'cthyb',
                                                        'seedname': [u'orbital_model'], 'dc_dmft': False,
                                                        'occ_conv_crit': -1,'g0_conv_crit': -1,'gimp_conv_crit': -1,'sigma_conv_crit': -1,
                                                        'J': [1.0], 'h5_save_freq': 2,
                                                        'dc': True, 'jobname': [u'out_60M_afm'],
                                                        'beta': 40.0, 'U': [5.5], 'diag_delta': False,
                                                        'energy_shift_orbitals': [11, 0, 9.3, 2.1],
                                                        'measure_chi_insertions': 100, 'h_field': 0.0,
                                                        'calc_energies': False, 'sigma_mix': 0.6, 'g0_mix': 1.0,
                                                        'g0_mix_type': 'linear',
                                                        'magmom': [1.0, -1.0], 'dc_type': 0, 'n_tau': 10001,
                                                        'load_sigma': False, 'n_l': 35, 'n_iw': 1025,
                                                        'noise_level_initial_sigma': 0., 'n_iter_dmft': 6,
                                                        'h_int_type': 'kanamori',
                                                        'mu_gap_gb2_threshold': 'none'},
                                            'advanced': {'dc_fixed_value': 'none', 'dc_fixed_occ': 'none',
                                                         'dc_nominal': False,
                                                         'dc_factor': 'none', 'dc_J': [1.0], 'dc_U': [5.5],
                                                         'map_solver_struct': {('ud_0', 0): ('up_0', 0)},
                                                         'mapped_solver_struct_degeneracies': 'none', },
                                            'solver': {'imag_threshold': 1e-05, 'n_warmup_cycles': 10000,
                                                       'length_cycle': 140, 'measure_G_l': True,
                                                       'n_cycles_tot': 60000000, 'store_solver': False,
                                                       'move_double': True, 'measure_pert_order': False,
                                                       'move_shift': False, 'legendre_fit' : False,
                                                       'measure_density_matrix': False, 'delta_interface': True,
                                                       'off_diag_threshold' : 0.1}}

                             )
