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
        self.assertOneProblem(problems, '(ctint)', 'such a type')

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



class test_rst_pages(unittest.TestCase):
    """The web pages are generated from the same parser as the consistency check."""

    def test_write_rst(self):
        sections, problems = documentation.parse_documentation()
        self.assertEqual(problems, [])
        with tempfile.TemporaryDirectory() as tmp:
            documentation.write_rst(tmp)
            with open(os.path.join(tmp, 'input.rst')) as file:
                index_page = file.read()
            for section in sections:
                with open(os.path.join(tmp, f'{section.name}.rst')) as file:
                    page = file.read()
                self.assertIn(f'    {section.name}\n', index_page)
                for item in section.items:
                    if isinstance(item, documentation.ParamDoc):
                        self.assertIn(f'.. admonition:: {item.name}\n', page)
                        self.assertIn(item.name, index_page)

    def test_description_line_with_colon_stays_description(self):
        # a description line containing " : " must not become an entry of its own
        doc = '[  general  ]\n-------------\n\nmagmom : list of float, default = None\n        -1 : spin down\n'
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'documentation.txt')
            with open(path, 'w') as file:
                file.write(doc)
            (section,), problems = documentation.parse_documentation(path)
        self.assertEqual(problems, [])
        self.assertEqual(len(section.items), 1)
        self.assertEqual(section.items[0].description, '-1 : spin down')


if __name__ == '__main__':
    unittest.main()
