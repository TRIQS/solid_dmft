#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains unit tests for observables.py
"""

from triqs.gf import BlockGf, Gf, SemiCircular, MeshImFreq, MeshImTime
from triqs.gf.descriptors import Fourier

from solid_dmft.dmft_tools.observables import add_dft_values_as_zeroth_iteration, add_dmft_observables, _generate_header
from helper import are_iterables_equal, Dummy


import unittest

def set_up_observables(n_inequiv_shells):
    observables = {}
    observables['iteration'] = []
    observables['mu'] = []
    observables['E_tot'] = []
    observables['E_bandcorr'] = []
    observables['E_int'] = [[] for _ in range(n_inequiv_shells)]
    observables['E_corr_en'] = []
    observables['E_dft'] = []
    observables['E_DC'] = [[] for _ in range(n_inequiv_shells)]
    observables['orb_gb2'] = [{'up': [], 'down': []} for _ in range(n_inequiv_shells)]
    observables['imp_gb2'] = [{'up': [], 'down': []} for _ in range(n_inequiv_shells)]
    observables['orb_occ'] = [{'up': [], 'down': []} for _ in range(n_inequiv_shells)]
    observables['imp_occ'] = [{'up': [], 'down': []} for _ in range(n_inequiv_shells)]
    observables['orb_Z'] =  [{'up': [], 'down': []} for _ in range(n_inequiv_shells)]

    return observables

class test_observables(unittest.TestCase):

    # ---------- add_dft_values_as_zeroth_iteration ----------
    def test_add_dft_values_one_impurity_one_band(self):
        sum_k = Dummy()
        sum_k.n_inequiv_shells = 1
        sum_k.dc_energ = [0.0]
        sum_k.inequiv_to_corr = [0]
        sum_k.gf_struct_solver = [{'up_0': None, 'down_0': None}]
        sum_k.spin_block_names = [['up', 'down'], ['ud']]
        sum_k.SO = 0

        general_params = {}
        general_params['calc_energies'] = False
        general_params['csc'] = False
        general_params['dc'] = False
        general_params['beta'] = 40.
        general_params['solver_type'] = 'cthyb'
        general_params['n_tau'] = 10001

        mesh_iw = MeshImFreq(beta=general_params['beta'], S='Fermion', n_iw = 1025)
        gf_up = Gf(mesh=mesh_iw, target_shape=[1,1])
        gf_down = gf_up.copy()
        gf_up << SemiCircular(2, 0)
        gf_down << SemiCircular(1, 0)
        G_loc_all_dft = [BlockGf(name_list=('up_0', 'down_0'), block_list=(gf_up, gf_down), make_copies=True)]

        solver_type_per_imp = ['cthyb']
        dft_mu = 12.
        shell_multiplicity = [4]

        observables = set_up_observables(sum_k.n_inequiv_shells)

        observables = add_dft_values_as_zeroth_iteration(observables, general_params, solver_type_per_imp, dft_mu, None, sum_k,
                                                         G_loc_all_dft, shell_multiplicity)

        expected_observables = {'E_bandcorr': ['none'], 'E_tot': ['none'], 'E_dft': ['none'],
                                'E_DC': [['none']], 'E_int': [['none']],
                                'orb_occ': [{'down': [[0.5]], 'up': [[0.5]]}],
                                'iteration': [0], 'E_corr_en': ['none'], 'mu': [12.0],
                                'orb_gb2': [{'down': [[-0.0498]], 'up': [[-0.0250]]}],
                                'imp_gb2': [{'down': [-0.0498], 'up': [-0.0250]}],
                                'imp_occ': [{'down': [0.5], 'up': [0.5]}],
                                'orb_Z': [{'down': [[1.0]], 'up': [[1.0]]}]}

        assert are_iterables_equal(observables, expected_observables)


    def test_add_dft_values_two_impurites_two_bands(self):
        sum_k = Dummy()
        sum_k.n_inequiv_shells = 2
        # mapping corr to inequiv = [0, 0, 1, 0]
        sum_k.dc_energ = [0.1, 0.1, 5.0, 0.1]
        sum_k.inequiv_to_corr = [0, 2]
        sum_k.gf_struct_solver = [{'up_0': None, 'down_0': None, 'up_1': None, 'down_1': None},
                                  {'up_0': None, 'down_0': None}]
        sum_k.spin_block_names = [['up', 'down'], ['ud']]
        sum_k.SO = 0

        general_params = {}
        general_params['calc_energies'] = True
        general_params['csc'] = False
        general_params['dc'] = True
        general_params['beta'] = 40.
        general_params['solver_type'] = 'cthyb'
        general_params['n_tau'] = 10001

        mesh_iw = MeshImFreq(beta=general_params['beta'], S='Fermion', n_iw = 1025)
        gf_up_one_band = Gf(mesh=mesh_iw, target_shape=[1,1])
        gf_down_one_band = gf_up_one_band.copy()
        gf_up_one_band << SemiCircular(2, 1)
        gf_down_one_band << SemiCircular(1, 2)

        gf_up_two_bands = Gf(mesh=mesh_iw, target_shape=[2, 2])
        gf_down_two_bands = gf_up_two_bands.copy()
        gf_up_two_bands << SemiCircular(2, 0)
        gf_down_two_bands << SemiCircular(1, 0)

        G_loc_all_dft = [BlockGf(name_list=('up_0', 'down_1', 'down_0', 'up_1'),
                                 block_list=(gf_up_one_band, gf_down_one_band, gf_down_one_band, gf_up_one_band),
                                 make_copies=True),
                         BlockGf(name_list=('up_0', 'down_0'), block_list=(gf_up_two_bands, gf_down_two_bands), make_copies=True)]

        solver_type_per_imp = ['cthyb', 'cthyb']
        dft_mu = 2.
        dft_energy = None
        shell_multiplicity = [3, 1]

        observables = set_up_observables(sum_k.n_inequiv_shells)

        observables = add_dft_values_as_zeroth_iteration(observables, general_params, solver_type_per_imp, dft_mu, dft_energy, sum_k,
                                                         G_loc_all_dft, shell_multiplicity)
        expected_observables = {'E_bandcorr': [0.0], 'E_tot': [-5.3], 'E_dft': [0.0],
                                'E_DC': [[0.3], [5.0]], 'E_int': [[0.0], [0.0]],
                                'orb_occ': [{'down': [[1., 1.]], 'up': [[0.8044, 0.8044]]}, {'down': [[0.5, 0.5]], 'up': [[0.5, 0.5]]}],
                                'iteration': [0], 'E_corr_en': [0.0], 'mu': [2.0],
                                'orb_gb2': [{'down': [[0, 0]], 'up': [[-0.0216, -0.0216]]}, {'down': [[-0.0498, -0.0498]], 'up': [[-0.0250, -0.0250]]}],
                                'imp_gb2': [{'down': [0], 'up': [-0.0432]}, {'down': [-0.0996], 'up': [-0.0500]}],
                                'imp_occ': [{'down': [2.], 'up': [1.6088]}, {'down': [1.], 'up': [1.]}],
                                'orb_Z': [{'down': [[1., 1.]], 'up': [[1., 1.]]}, {'down': [[1., 1.]], 'up': [[1., 1.]]}],
                               }

        assert are_iterables_equal(observables, expected_observables)


    # ---------- add_dmft_observables ----------
    def test_add_dmft_observables_one_impurity_one_band(self):
        sum_k = Dummy()
        sum_k.n_inequiv_shells = 1
        sum_k.dc_energ = [0.0]
        sum_k.inequiv_to_corr = [0]
        sum_k.gf_struct_solver = [{'up_0': None, 'down_0': None}]
        sum_k.spin_block_names = [['up', 'down'], ['ud']]
        sum_k.SO = 0

        general_params = {}
        general_params['calc_energies'] = False
        general_params['csc'] = False
        general_params['dc'] = False
        general_params['beta'] = 40.
        general_params['solver_type'] = 'cthyb'
        general_params['n_tau'] = 10001

        solver_params = [{'measure_density_matrix': False}]

        observables = set_up_observables(sum_k.n_inequiv_shells)

        it = 1
        h_int = None
        previous_mu = 1.56
        shell_multiplicity = [4]
        E_bandcorr = 10.43

        mesh_iw = MeshImFreq(beta=general_params['beta'], S='Fermion', n_iw = 1025)
        gf_up = Gf(mesh=mesh_iw, target_shape=[1,1])
        gf_down = gf_up.copy()
        gf_up << SemiCircular(2, 0)
        gf_down << SemiCircular(1, 0)
        gf_iw = [BlockGf(name_list=('up_0', 'down_0'), block_list=(gf_up, gf_down), make_copies=True)]

        mesh_tau = MeshImTime(beta=general_params['beta'], S='Fermion', n_tau = 10001)
        gf_imtime = Gf(mesh=mesh_tau, target_shape=[1,1])
        gf_tau = [BlockGf(name_list=('up_0', 'down_0'), block_list=(gf_imtime, gf_imtime), make_copies=True)]
        gf_tau[0] << Fourier(gf_iw[0])

        solvers = [Dummy()]
        solvers[0].G_time = gf_tau[0]
        # misusing here G_freq to get some comparable dummy values for Z
        solvers[0].Sigma_freq = gf_iw[0]

        map_imp_solver = [0]
        solver_type_per_imp = ['cthyb']
        density_mat = [gf_iw[iineq].density() for iineq in range(sum_k.n_inequiv_shells)]

        observables = add_dmft_observables(observables, general_params, solver_params, map_imp_solver, solver_type_per_imp, None, it, solvers, h_int,
                                           previous_mu, sum_k, density_mat, shell_multiplicity, E_bandcorr)
        print(observables['orb_Z'])
        expected_observables = {'iteration': [1], 'mu': [1.56], 'E_tot': [10.43], 'E_bandcorr': [10.43],
                                'E_int': [[]], 'E_corr_en': [0.0], 'E_dft': [0.0], 'E_DC': [[0.0]],
                                'orb_gb2': [{'up': [[-0.0250]], 'down': [[-0.0498]]}], 'imp_gb2': [{'up': [-0.0250], 'down': [-0.0498]}],
                                'orb_occ': [{'up': [[0.5]], 'down': [[0.5]]}],
                                'imp_occ': [{'up': [0.5], 'down': [0.5]}],
                                'orb_Z':   [{'up': [[1.85487611]], 'down': [[1.4481126]]}]
                                }

        assert are_iterables_equal(observables, expected_observables)


    # ---------- _generate_header ----------
    def test_generate_header(self):
        general_params = {}
        general_params['magnetic'] = False
        general_params['calc_energies'] = False
        general_params['csc'] = False
        sum_k = Dummy()
        sum_k.n_inequiv_shells = 1
        sum_k.inequiv_to_corr = [0]
        sum_k.gf_struct_solver = [{'up_0': 3, 'down_0': 3}]
        sum_k.corr_shells = [{'dim': 3}]
        sum_k.SO = 0

        headers = _generate_header(general_params, sum_k)
        expected_headers = {'observables_imp0.dat':
                                ' it |         mu |                G(beta/2) per orbital |'
                                + '                 orbital occs up+down |      impurity occ'}
        assert headers == expected_headers


        general_params = {}
        general_params['magnetic'] = False
        general_params['calc_energies'] = True
        general_params['csc'] = False
        sum_k = Dummy()
        sum_k.n_inequiv_shells = 2
        sum_k.inequiv_to_corr = [0, 2]
        sum_k.gf_struct_solver = [{'up_0': 3, 'down_0': 3},
                                  {'up_0': 1, 'down_0': 1}]
        sum_k.corr_shells = [{'dim': 3}, {'dim': 3}, {'dim': 1}, {'dim': 1}]
        sum_k.SO = 0

        headers = _generate_header(general_params, sum_k)
        expected_headers = {'observables_imp0.dat':
                                ' it |         mu |                G(beta/2) per orbital |'
                                + '                 orbital occs up+down |      impurity occ |'
                                + '      E_tot |      E_DFT   E_bandcorr    E_int_imp         E_DC',
                            'observables_imp1.dat':
                                ' it |         mu | G(beta/2) per orbital |'
                                + '  orbital occs up+down |      impurity occ |'
                                + '      E_tot |      E_DFT   E_bandcorr    E_int_imp         E_DC'}
        assert headers == expected_headers

        general_params = {}
        general_params['magnetic'] = True
        general_params['calc_energies'] = True
        general_params['csc'] = False
        sum_k = Dummy()
        sum_k.n_inequiv_shells = 2
        sum_k.inequiv_to_corr = [0, 1]
        sum_k.gf_struct_solver = [{'up_0': 2, 'down_0': 2}, {'up_0': 1, 'down_0': 1}]
        sum_k.corr_shells = [{'dim': 2}, {'dim': 1}]
        sum_k.SO = 0

        headers = _generate_header(general_params, sum_k)
        expected_headers = {'observables_imp0_up.dat':
                                ' it |         mu |   G(beta/2) per orbital |'
                                + '         orbital occs up |   impurity occ up |'
                                + '      E_tot |      E_DFT   E_bandcorr    E_int_imp         E_DC',
                            'observables_imp0_down.dat':
                                ' it |         mu |   G(beta/2) per orbital |'
                                + '       orbital occs down | impurity occ down |'
                                + '      E_tot |      E_DFT   E_bandcorr    E_int_imp         E_DC',
                            'observables_imp1_up.dat':
                                ' it |         mu | G(beta/2) per orbital |'
                                + '       orbital occs up |   impurity occ up |'
                                + '      E_tot |      E_DFT   E_bandcorr    E_int_imp         E_DC',
                            'observables_imp1_down.dat':
                                ' it |         mu | G(beta/2) per orbital |'
                                + '     orbital occs down | impurity occ down |'
                                + '      E_tot |      E_DFT   E_bandcorr    E_int_imp         E_DC'}
        assert headers == expected_headers

if __name__ == '__main__':
    unittest.main()
