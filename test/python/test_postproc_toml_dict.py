from solid_dmft.io_tools.postproc_toml_dict import merge_config_with_default, _apply_default_values, _verify_restrictions_on_default_and_config, _resolve_references

import unittest

class test_postproc_toml_dict(unittest.TestCase):
    def test_verify_restrictions_invalid_section(self):
        config = {
            'section1': {'key1': 'newval11', 'key2': 'newval12'},
            'section4': {'key1': 'newval41'},
        }
        default = {
            'section1': {'key1': 'defval11', 'key2': 'defval12', 'key3': 'defval13'},
            'section2': {'key1': 'defval21', 'key2': 'defval22', 'key3': 'defval23'},
            'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
        }
        with self.assertRaises(ValueError):
            _verify_restrictions_on_default_and_config(config, default, {})

    def test_verify_restrictions_key_outside_section_default(self):
        config = {
            'section1': {'key1': 'newval11', 'key2': 'newval12'},
            'section2': {'key1': 'newval31'},
        }
        default = {
            'section1': {'key1': 'defval11', 'key2': 'defval12', 'key3': 'defval13'},
            'section2': {'key1': 'defval21', 'key2': 'defval22', 'key3': 'defval23'},
            'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
            'invalid_key': 'value',
        }
        with self.assertRaises(ValueError):
            _verify_restrictions_on_default_and_config(config, default, {})

    def test_verify_restrictions_key_outside_section_config(self):
        config = {
            'section1': {'key1': 'newval11', 'key2': 'newval12'},
            'section2': 'value',
        }
        default = {
            'section1': {'key1': 'defval11', 'key2': 'defval12', 'key3': 'defval13'},
            'section2': {'key1': 'defval21', 'key2': 'defval22', 'key3': 'defval23'},
            'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
        }
        with self.assertRaises(ValueError):
            _verify_restrictions_on_default_and_config(config, default, {})

    def test_verify_restriction_missing_listed_section(self):
        config = {}
        default = {
            'section1': [{'match': 'defval11', 'key2': 'defval12', 'key3': 'defval13'}],
            'section2': {'key1': 'defval21', 'key2': 'defval22', 'key3': 'defval23'},
            'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
        }

        with self.assertRaises(ValueError):
            _verify_restrictions_on_default_and_config(default, config, {'section1': 'match'})

    def test_merge_config_nonexistent_listed_section(self):
        config = {
            'section1': {'key1': 'newval11', 'key2': 'newval12'},
            'section2': [{'key1': 'newval31'}],
        }
        default = {
            'section1': {'key1': 'defval11', 'key2': 'defval12', 'key3': 'defval13'},
            'section2': {'key1': 'defval21', 'key2': 'defval22', 'key3': 'defval23'},
            'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
        }
        with self.assertRaises(ValueError):
            merge_config_with_default(config, default, {})

    def test_resolve_references_simple(self):
        config = {
            'section1': {'key1': '<section2.key3>', 'key2': 'defval12', 'key3': 'defval13'},
            'section2': {'key1': 'defval21', 'key2': 'defval22', 'key3': 'defval23'},
            'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
        }
        expected_output = {
            'section1': {'key1': 'defval23', 'key2': 'defval12', 'key3': 'defval13'},
            'section2': {'key1': 'defval21', 'key2': 'defval22', 'key3': 'defval23'},
            'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
        }

        _resolve_references(config['section1'], 'section1', config)
        assert config == expected_output

    def test_resolve_references_chained_reference(self):
        config = {
            'section1': {'key1': '<section2.key3>', 'key2': 'defval12', 'key3': 'defval13'},
            'section2': {'key1': 'defval21', 'key2': 'defval22', 'key3': '<section3.key3>'},
            'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
        }
        with self.assertRaises(ValueError):
            _resolve_references(config['section1'], 'section1', config)

    def test_apply_default_values_partial_config(self):
        config = {
            'section1': {'key1': 'newval11', 'key2': 'newval12'},
            'section2': {'key1': 'newval21', 'key2': 'newval22', 'key3': 'newval23'},
            'section3': {'key1': 'newval31'},
        }
        default = {
            'section1': {'key1': 'defval11', 'key2': 'defval12', 'key3': 'defval13'},
            'section2': {'key1': 'defval21', 'key2': 'defval22', 'key3': 'defval23'},
            'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
        }
        expected_output = {
            'section1': {'key1': 'newval11', 'key2': 'newval12', 'key3': 'defval13'},
            'section2': {'key1': 'newval21', 'key2': 'newval22', 'key3': 'newval23'},
            'section3': {'key1': 'newval31', 'key2': 'defval32', 'key3': 'defval33'},
        }
        assert _apply_default_values(config, default, {}) == expected_output

    def test_apply_default_values_complete_listed_config(self):
        config = {
            'section1': {'key1': 'newval11', 'key2': 'newval12', 'key3': 'newval13'},
            'section2': [{'match': 'matchval', 'key2': 'newval22', 'key3': 'newval23'}],
            'section3': {'key1': 'newval31', 'key2': 'newval32', 'key3': 'newval33'},
        }
        default = {
            'section1': {'key1': 'defval11', 'key2': 'defval12', 'key3': 'defval13'},
            'section2': [{'match': 'matchval', 'key2': 'defval22', 'key3': 'defval23'}],
            'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
        }
        assert _apply_default_values(config, default, {'section2': 'match'}) == config

    def test_apply_default_values_partial_listed_config(self):
        config = {
            'section1': {'key1': 'newval11', 'key2': 'newval12'},
            'section2': [{'match': 'matchval', 'key2': 'newval22', 'key3': 'newval23'}],
            'section3': {'key1': 'newval31'},
        }
        default = {
            'section1': {'key1': 'defval11', 'key2': 'defval12', 'key3': 'defval13'},
            'section2': [{'match': 'matchval', 'key2': 'defval22', 'key3': 'defval23'}],
            'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
        }
        expected_output = {
            'section1': {'key1': 'newval11', 'key2': 'newval12', 'key3': 'defval13'},
            'section2': [{'match': 'matchval', 'key2': 'newval22', 'key3': 'newval23'}],
            'section3': {'key1': 'newval31', 'key2': 'defval32', 'key3': 'defval33'},
        }
        assert _apply_default_values(config, default, {'section2': 'match'}) == expected_output

    def test_apply_default_values_invalid_key(self):
        config = {
            'section1': {'key1': 'newval11', 'key2': 'newval12'},
            'section2': {'key4': 'newval24'},
        }
        default = {
            'section1': {'key1': 'defval11', 'key2': 'defval12', 'key3': 'defval13'},
            'section2': {'key1': 'defval21', 'key2': 'defval22', 'key3': 'defval23'},
            'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
        }
        with self.assertRaises(ValueError):
            _apply_default_values(config, default, {})

    def test_apply_default_values_invalid_listed_key(self):
        config = {
            'section1': {'key1': 'newval11', 'key2': 'newval12'},
            'section2': [{'match': 'matchval', 'key4': 'newval24'}],
        }
        default = {
            'section1': {'key1': 'defval11', 'key2': 'defval12', 'key3': 'defval13'},
            'section2': [{'match': 'matchval', 'key2': 'defval22', 'key3': 'defval23'}],
            'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
        }
        with self.assertRaises(ValueError):
            _apply_default_values(config, default, {'section2': 'match'})

if __name__ == '__main__':
    unittest.main()
