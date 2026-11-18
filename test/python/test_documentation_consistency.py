import os
import tempfile
import textwrap
import unittest

from solid_dmft.io_tools import documentation


class test_documentation_consistency(unittest.TestCase):
    """documentation.txt is the source of the input documentation and of `solid_dmft --params`,
    default.toml is what the code uses. Both have to list the same keys with the same defaults."""

    def test_documentation_matches_default_toml(self):
        problems = documentation.check_consistency()
        self.assertEqual(problems, [], 'documentation.txt and default.toml disagree:\n' + '\n'.join(problems))

    def test_list_params(self):
        params = documentation.list_params(section='solver', solver_type='ctseg')
        self.assertTrue(params)
        self.assertTrue(all(p['solver_type'] == 'ctseg' for p in params))
        by_key = {p['key']: p for p in params}
        # per-solver default, differs from cthyb
        self.assertIs(by_key['measure_pert_order']['default'], True)
        self.assertTrue(by_key['n_cycles_tot']['mandatory'])
        self.assertNotIn('move_shift', by_key)
        # documented once for all solver types
        self.assertIn('list of impurities', by_key['idx_impurities']['description'])

        dc_u, = documentation.list_params(section='advanced', key='dc_U')
        self.assertIsNone(dc_u['default'])
        self.assertEqual(dc_u['default_ref'], 'general.U')

        seedname, = documentation.list_params(key='seedname')
        self.assertTrue(seedname['mandatory'])
        self.assertEqual(seedname['type'], 'str')


class test_check_consistency_detects_drift(unittest.TestCase):
    DOC = textwrap.dedent('''\
        [  general  ]
        -------------

        beta : float, default = None
                inverse temperature
        n_iw : int, default = 1025
                number of Matsubara frequencies
        seedname : str, mandatory
                seedname
        jobname : str, default = general.seedname
                output directory

        [  solver  ]
        ------------
        type : str, mandatory
                solver type
        idx_impurities : list of int, default = None
                impurities of this solver

        cthyb
        =====
        n_cycles_tot : int, mandatory
                total number of sweeps
        move_double : bool, default = True
                double moves
        ''')

    TOML = textwrap.dedent('''\
        [general]
        beta = "<none>"
        n_iw = 1025
        seedname = "<no default>"
        jobname = "<general.seedname>"

        [[solver]]
        type = "cthyb"
        idx_impurities = "<none>"
        n_cycles_tot = "<no default>"
        move_double = true
        ''')

    def _check(self, doc, toml):
        with tempfile.TemporaryDirectory() as tmp:
            doc_path = os.path.join(tmp, 'documentation.txt')
            toml_path = os.path.join(tmp, 'default.toml')
            with open(doc_path, 'w') as file:
                file.write(doc)
            with open(toml_path, 'w') as file:
                file.write(toml)
            return documentation.check_consistency(doc_path, toml_path)

    def assertOneProblem(self, problems, *fragments):
        self.assertEqual(len(problems), 1, problems)
        for fragment in fragments:
            self.assertIn(fragment, problems[0])

    def test_consistent(self):
        self.assertEqual(self._check(self.DOC, self.TOML), [])

    def test_undocumented_key(self):
        problems = self._check(self.DOC, self.TOML + 'max_time = -1\n')
        self.assertOneProblem(problems, '(cthyb) max_time', 'not in documentation.txt')

    def test_key_only_documented(self):
        problems = self._check(self.DOC, self.TOML.replace('n_iw = 1025\n', ''))
        self.assertOneProblem(problems, '[general] n_iw', 'not in default.toml')

    def test_wrong_default(self):
        problems = self._check(self.DOC, self.TOML.replace('n_iw = 1025', 'n_iw = 501'))
        self.assertOneProblem(problems, '[general] n_iw', 'differs')

    def test_mandatory_mismatch(self):
        problems = self._check(self.DOC, self.TOML.replace('seedname = "<no default>"', 'seedname = "<none>"'))
        self.assertOneProblem(problems, '[general] seedname', 'mandatory')

    def test_wrong_reference(self):
        problems = self._check(self.DOC, self.TOML.replace('"<general.seedname>"', '"<general.beta>"'))
        self.assertOneProblem(problems, '[general] jobname', 'differs')

    def test_common_key_missing_for_a_solver(self):
        problems = self._check(self.DOC, self.TOML.replace('idx_impurities = "<none>"\n', ''))
        self.assertOneProblem(problems, '(cthyb) idx_impurities', 'not in default.toml')

    def test_unknown_solver_subsection(self):
        doc = self.DOC + 'ctint\n=====\nmove_double : bool, default = True\n        double moves\n'
        problems = self._check(doc, self.TOML)
        self.assertOneProblem(problems, '(ctint)', 'no such type')

    def test_duplicate_entry(self):
        doc = self.DOC + 'move_double : bool, default = True\n        double moves\n'
        problems = self._check(doc, self.TOML)
        self.assertOneProblem(problems, 'move_double', 'more than once')

    def test_missing_default(self):
        problems = self._check(self.DOC.replace('n_iw : int, default = 1025', 'n_iw : int'), self.TOML)
        self.assertOneProblem(problems, '[general] n_iw', 'neither a default nor "mandatory"')

    def test_unparsable_entry(self):
        problems = self._check(self.DOC.replace('n_iw : int, default = 1025', 'n_iw : int, default 1025'),
                               self.TOML)
        self.assertIn('cannot parse', problems[0])

    def test_type_mismatch(self):
        problems = self._check(self.DOC.replace('move_double : bool', 'move_double : int'), self.TOML)
        self.assertOneProblem(problems, 'move_double', 'does not match documented type')

    def test_unknown_type(self):
        problems = self._check(self.DOC.replace('move_double : bool', 'move_double : boold'), self.TOML)
        self.assertOneProblem(problems, 'move_double', 'unknown type')


if __name__ == '__main__':
    unittest.main()
