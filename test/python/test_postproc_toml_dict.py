import pytest
from solid_dmft.io_tools.postproc_toml_dict import _apply_default_values, _verify_restrictions_on_default_and_config, _resolve_references

def test_verify_restrictions_invalid_key():
    config = {
        'section1': {'key1': 'newval11', 'key2': 'newval12'},
        'section4': {'key1': 'newval31'},
    }
    default = {
        'section1': {'key1': 'defval11', 'key2': 'defval12', 'key3': 'defval13'},
        'section2': {'key1': 'defval21', 'key2': 'defval22', 'key3': 'defval23'},
        'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
    }
    with pytest.raises(ValueError):
        _verify_restrictions_on_default_and_config(config, default, {})

def test_verify_restrictions_key_outside_section_default():
    config = {
        'section1': {'key1': 'newval11', 'key2': 'newval12'},
        'section4': {'key1': 'newval31'},
    }
    default = {
        'section1': {'key1': 'defval11', 'key2': 'defval12', 'key3': 'defval13'},
        'section2': {'key1': 'defval21', 'key2': 'defval22', 'key3': 'defval23'},
        'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
        'invalid_key': 'value',
    }
    with pytest.raises(ValueError):
        _verify_restrictions_on_default_and_config(config, default, {})

def test_verify_restrictions_key_outside_section_config():
    config = {
        'section1': {'key1': 'newval11', 'key2': 'newval12'},
        'section2': 'value',
    }
    default = {
        'section1': {'key1': 'defval11', 'key2': 'defval12', 'key3': 'defval13'},
        'section2': {'key1': 'defval21', 'key2': 'defval22', 'key3': 'defval23'},
        'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
    }
    with pytest.raises(ValueError):
        _verify_restrictions_on_default_and_config(config, default, {})

def test_verify_restriction_missing_listed_section():
    config = {}
    default = {
        'section1': [{'match': 'defval11', 'key2': 'defval12', 'key3': 'defval13'}],
        'section2': {'key1': 'defval21', 'key2': 'defval22', 'key3': 'defval23'},
        'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
    }

    with pytest.raises(ValueError):
        _verify_restrictions_on_default_and_config(default, config, {'section1': 'match'})

def test_verify_restrictions_nonexistent_listed_section():
    config = {
        'section1': {'key1': 'newval11', 'key2': 'newval12'},
        'section2': [{'key1': 'newval31'}],
    }
    default = {
        'section1': {'key1': 'defval11', 'key2': 'defval12', 'key3': 'defval13'},
        'section2': {'key1': 'defval21', 'key2': 'defval22', 'key3': 'defval23'},
        'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
    }
    with pytest.raises(ValueError):
        _verify_restrictions_on_default_and_config(config, default, {})

def test_resolve_references_simple():
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

def test_resolve_references_chained_reference():
    config = {
        'section1': {'key1': '<section2.key3>', 'key2': 'defval12', 'key3': 'defval13'},
        'section2': {'key1': 'defval21', 'key2': 'defval22', 'key3': '<section3.key3>'},
        'section3': {'key1': 'defval31', 'key2': 'defval32', 'key3': 'defval33'},
    }
    with pytest.raises(ValueError):
        _resolve_references(config['section1'], 'section1', config)

def test_apply_default_values_partial_config():
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

def test_apply_default_values_complete_listed_config():
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

def test_apply_default_values_partial_listed_config():
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
