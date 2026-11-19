import contextlib
import io
import json
import os
import tempfile
import textwrap
import unittest

from solid_dmft.io_tools import params_cli

VALID_CONFIG = textwrap.dedent('''\
    [general]
    seedname = "inp"
    jobname = "out"
    beta = 40
    n_iter_dmft = 2
    h_int_type = "kanamori"
    U = 4.0
    J = 0.6

    [solver]
    type = "ctseg"
    length_cycle = 50
    n_warmup_cycles = 1000
    n_cycles_tot = 10000
    ''')


def run(args):
    """Returns exit code, stdout and stderr of params_cli.main(args)."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = params_cli.main(args)
    return code, out.getvalue(), err.getvalue()


class test_params_cli(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def write_config(self, text):
        path = os.path.join(self.tmp.name, 'dmft_config.toml')
        with open(path, 'w') as file:
            file.write(text)
        return path

    def test_is_cli_call(self):
        self.assertFalse(params_cli.is_cli_call(['solid_dmft']))
        self.assertFalse(params_cli.is_cli_call(['solid_dmft', 'my_config.toml']))
        self.assertTrue(params_cli.is_cli_call(['solid_dmft', '--params']))
        self.assertTrue(params_cli.is_cli_call(['solid_dmft', '--check', 'my_config.toml']))

    def test_params_solver_type_json(self):
        code, out, _ = run(['--params', 'solver', '--solver-type', 'ctseg', '--json'])
        self.assertEqual(code, 0)
        params = json.loads(out)
        self.assertTrue(all(p['section'] == 'solver' and p['solver_type'] == 'ctseg' for p in params))
        n_cycles, = [p for p in params if p['key'] == 'n_cycles_tot']
        self.assertEqual(set(n_cycles), {'section', 'solver_type', 'key', 'type', 'mandatory',
                                         'default', 'default_ref', 'description'})
        self.assertTrue(n_cycles['mandatory'])

    def test_params_key_text(self):
        code, out, _ = run(['--params', 'general', 'n_iw'])
        self.assertEqual(code, 0)
        self.assertIn('[general] n_iw', out)
        self.assertIn('default: 1025', out)
        self.assertIn('Matsubara', out)

    def test_params_key_without_section(self):
        code, out, _ = run(['--params', 'dc_U', '--json'])
        self.assertEqual(code, 0)
        self.assertEqual([p['section'] for p in json.loads(out)], ['advanced'])

    def test_params_unknown_key_suggests(self):
        code, _, err = run(['--params', 'solver', 'n_cycle_tot'])
        self.assertEqual(code, 1)
        self.assertIn('did you mean "n_cycles_tot"', err)

    def test_params_unknown_section(self):
        code, _, err = run(['--params', 'genral', 'n_iw'])
        self.assertEqual(code, 1)
        self.assertIn('did you mean "general"', err)

    def test_check_valid(self):
        code, out, _ = run(['--check', self.write_config(VALID_CONFIG)])
        self.assertEqual(code, 0)
        self.assertIn('config is valid', out)
        self.assertRegex(out, r'n_iw +1025  \(default\)')

    def test_check_valid_json(self):
        code, out, _ = run(['--check', '--json', self.write_config(VALID_CONFIG)])
        self.assertEqual(code, 0)
        result = json.loads(out)
        self.assertTrue(result['valid'])
        self.assertEqual(result['config']['general']['n_iw'], 1025)
        self.assertIs(result['config']['solver'][0]['measure_pert_order'], True)

    def test_check_unknown_keys(self):
        config = VALID_CONFIG.replace('beta = 40', 'betta = 40').replace('n_cycles_tot', 'n_cycle_tot')
        config += '\n[dft]\nU = 3.0\n'
        code, out, _ = run(['--check', '--json', self.write_config(config)])
        self.assertEqual(code, 1)
        errors = json.loads(out)['errors']
        self.assertEqual(len(errors), 3, errors)
        self.assertIn('did you mean "beta"', errors[0])
        self.assertIn('did you mean "n_cycles_tot"', errors[1])
        self.assertIn('"U" is valid in [general]', errors[2])

    def test_check_invalid_value(self):
        code, out, _ = run(['--check', self.write_config(VALID_CONFIG.replace('"kanamori"', '"kanamori2"'))])
        self.assertEqual(code, 1)
        self.assertIn('h_int_type', out)

    def test_check_missing_mandatory(self):
        code, out, _ = run(['--check', self.write_config(VALID_CONFIG.replace('length_cycle = 50\n', ''))])
        self.assertEqual(code, 1)
        self.assertIn('"length_cycle" in section "solver" is mandatory', out)

    def test_check_does_not_write(self):
        path = self.write_config(VALID_CONFIG)
        run(['--check', path])
        self.assertEqual(os.listdir(self.tmp.name), ['dmft_config.toml'])


if __name__ == '__main__':
    unittest.main()
