from solid_dmft.io_tools.postproc_toml_dict import merge_config_with_default

import unittest

class test_postproc_toml_dict_integration(unittest.TestCase):

    def setUp(self):
        self.default_simple = {'general': {'U': 4.0, 'J': 0.9, 'beta': 40.0},
                          'solver': {'type': 'cthyb', 'n_cycles_tot': 1000000}}

        self.default = {'general': {'U': '<no default>', 'J': '<no default>',
                               'beta': '<no default>', 'seedname': '<no default>',
                               'jobname': '<general.seedname>', 'n_tau': 100001, 'w_range': [-10, 10]},
                   'solver': {'type': 'cthyb', 'idx_imp': '<none>',
                              'n_cycles_tot': '<no default>', 'off_diag_threshold': 0.0}}

        self.default_solver_matching = {'general': {'U': '<no default>', 'J': '<no default>',
                                               'beta': '<no default>', 'seedname': '<no default>',
                                               'jobname': '<general.seedname>', 'n_tau': 100001},
                                   'solver': [{'type': 'cthyb', 'idx_imp': '<none>',
                                               'n_cycles_tot': '<no default>', 'off_diag_threshold': 0.0},
                                              {'type': 'ctint', 'idx_imp': '<none>',
                                               'n_cycles_tot': '<no default>', 'random_seed': '<none>'}]}

    def test_simple_empty_config(self):
        config = {}
        assert merge_config_with_default(config, self.default_simple) == self.default_simple

    def test_simple_complete_config(self):
        config = {'general': {'U': 5.0, 'J': 0, 'beta': 40.0},
                  'solver': {'type': 'ctint', 'n_cycles_tot': 2000000}}
        assert merge_config_with_default(config, self.default_simple) == config

    def test_simple_partial_config(self):
        config = {'general': {'U': 5.0},
                  'solver': {'type': 'ctint'}}
        expected_output = {'general': {'U': 5.0, 'J': 0.9, 'beta': 40.0},
                           'solver': {'type': 'ctint', 'n_cycles_tot': 1000000}}
        assert merge_config_with_default(config, self.default_simple) == expected_output

    def test_minimal_input(self):
        config = {'general': {'U': 4.0, 'J': 0.9, 'beta': 40.0, 'seedname': 'test'},
                  'solver': {'n_cycles_tot': 1000000}}

        expected_output = {'general': {'U': 4.0, 'J': 0.9,
                                       'beta': 40.0, 'seedname': 'test',
                                       'jobname': 'test', 'n_tau': 100001, 'w_range': [-10, 10]},
                           'solver': {'type': 'cthyb', 'idx_imp': None,
                                      'n_cycles_tot': 1000000, 'off_diag_threshold': 0.0}}

        assert merge_config_with_default(config, self.default) == expected_output

    def test_missing_mandatory_field(self):
        config = {'general': {'U': 4.0, 'J': 0.9, 'beta': 40.0},
                  'solver': {'n_cycles_tot': 1000000}}

        with self.assertRaises(ValueError):
            merge_config_with_default(config, self.default)

    def test_minimal_input_solver_matching(self):
        config = {'general': {'U': 4.0, 'J': 0.9, 'beta': 40.0, 'seedname': 'test'},
                  'solver': {'type': 'ctint', 'n_cycles_tot': 1000000, 'random_seed': 1234}}

        expected_output = {'general': {'U': 4.0, 'J': 0.9,
                                       'beta': 40.0, 'seedname': 'test',
                                       'jobname': 'test', 'n_tau': 100001},
                           'solver': [{'type': 'ctint', 'idx_imp': None,
                                       'n_cycles_tot': 1000000, 'random_seed': 1234}]}

        assert merge_config_with_default(config, self.default_solver_matching,
                                         {'solver': 'type'}) == expected_output

    def test_multiple_solvers_matching(self):
        config = {'general': {'U': 4.0, 'J': 0.9, 'beta': 40.0, 'seedname': 'test'},
                  'solver': [{'type': 'ctint', 'n_cycles_tot': 1000000, 'random_seed': 1234},
                             {'type': 'cthyb', 'n_cycles_tot': 1000000}]}

        expected_output = {'general': {'U': 4.0, 'J': 0.9,
                                       'beta': 40.0, 'seedname': 'test',
                                       'jobname': 'test', 'n_tau': 100001},
                           'solver': [{'type': 'ctint', 'idx_imp': None,
                                       'n_cycles_tot': 1000000, 'random_seed': 1234},
                                      {'type': 'cthyb', 'idx_imp': None,
                                       'n_cycles_tot': 1000000, 'off_diag_threshold': 0.0}]}

        assert merge_config_with_default(config, self.default_solver_matching,
                                         {'solver': 'type'}) == expected_output

    def test_multiple_same_solvers_matching(self):
        config = {'general': {'U': 4, 'J': 0.9, 'beta': 40.0, 'seedname': 'test'},
                  'solver': [{'type': 'ctint', 'n_cycles_tot': 1000000, 'random_seed': '<general.U>'},
                             {'type': 'ctint', 'n_cycles_tot': 2000000}]}

        expected_output = {'general': {'U': 4, 'J': 0.9,
                                       'beta': 40.0, 'seedname': 'test',
                                       'jobname': 'test', 'n_tau': 100001},
                           'solver': [{'type': 'ctint', 'idx_imp': None,
                                       'n_cycles_tot': 1000000, 'random_seed': 4},
                                      {'type': 'ctint', 'idx_imp': None,
                                       'n_cycles_tot': 2000000, 'random_seed': None}]}

        assert merge_config_with_default(config, self.default_solver_matching,
                                         {'solver': 'type'}) == expected_output

    def test_unmatched_section(self):
        config = {'general': {'U': 4.0, 'J': 0.9, 'beta': 40.0, 'seedname': 'test'},
                  'solver': {'type': 'ctint', 'n_cycles_tot': 1000000, 'random_seed': 1234}}

        with self.assertRaises(ValueError):
            merge_config_with_default(config, self.default_solver_matching, {'solver': 'type', 'unknown': 'unknown'})

if __name__ == '__main__':
    unittest.main()
