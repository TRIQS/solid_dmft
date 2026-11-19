from typing import Any, Dict, List, Union
import copy
import difflib
import os

try:
    import tomllib
except ImportError:
    import tomli as tomllib

ParamDict = Dict[str, Any]
FullConfig = Dict[str, Union[ParamDict, List[ParamDict]]]

DEFAULT_CONFIG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'default.toml')
# Sections that are lists of tables in default.toml, with the key that selects the entry
MATCH_KEY = {'solver': 'type'}
MANDATORY = '<no default>'
NONE = '<none>'


def load_default_config(path=DEFAULT_CONFIG_FILE):
    """ Reads the default config of solid_dmft (default.toml). """
    with open(path, 'rb') as file:
        return tomllib.load(file)


def is_reference(value):
    """ True for a default of the form "<section.key>", taken from another key. """
    return (isinstance(value, str) and value.startswith('<') and value.endswith('>')
            and value not in (MANDATORY, NONE))


def iter_default_keys(cfg_def, match_key=MATCH_KEY):
    """ Yields (section, match value, key, default) for every key of the default config;
    the match value is the entry's match_key value (the solver type) in listed sections, else None. """
    for section_name, section in cfg_def.items():
        if section_name in match_key:
            for entry in section:
                for key, value in entry.items():
                    yield section_name, entry[match_key[section_name]], key, value
        else:
            for key, value in section.items():
                yield section_name, None, key, value


def did_you_mean(word, candidates):
    matches = difflib.get_close_matches(str(word), [str(c) for c in candidates], n=2, cutoff=0.6)
    return ', did you mean ' + ' or '.join(f'"{m}"' for m in matches) + '?' if matches else ''


def find_unknown_entries(cfg_inp, cfg_def, match_key=MATCH_KEY):
    """
    All unknown sections, listed-section entries (e.g. solver types) and keys of the config,
    with suggestions for close matches and the sections where an unknown key is valid.

    Returns a list of messages, empty if the config only uses known names.
    """
    errors = []
    keys_per_section = {name: set(sec) for name, sec in cfg_def.items() if name not in match_key}
    keys_per_entry = {name: {entry[key]: set(entry) for entry in cfg_def[name]} for name, key in match_key.items()
                      if name in cfg_def}

    def valid_elsewhere(key, own):
        places = [f'[{name}]' for name, keys in keys_per_section.items() if key in keys and name != own]
        places += [f'[{name}] {match_key[name]} = "{value}"' for name, entries in keys_per_entry.items()
                   for value, keys in entries.items() if key in keys and value != own]
        return f' ("{key}" is valid in {", ".join(places)})' if places else ''

    for section_name, section in cfg_inp.items():
        if section_name not in cfg_def:
            errors.append(f'unknown section [{section_name}]' + did_you_mean(section_name, cfg_def))
            continue
        if section_name in match_key:
            key_name = match_key[section_name]
            for entry in (section if isinstance(section, list) else [section]):
                if not isinstance(entry, dict) or key_name not in entry:
                    errors.append(f'[{section_name}] "{key_name}" is mandatory, one of '
                                  + ', '.join(f'"{v}"' for v in keys_per_entry[section_name]))
                    continue
                value = entry[key_name]
                if value not in keys_per_entry[section_name]:
                    errors.append(f'[{section_name}] unknown {key_name} = "{value}"'
                                  + did_you_mean(value, keys_per_entry[section_name]))
                    continue
                valid = keys_per_entry[section_name][value]
                errors += [f'[{section_name}] unknown key "{key}" for {key_name} = "{value}"'
                           + did_you_mean(key, valid) + valid_elsewhere(key, value)
                           for key in entry if key not in valid]
        elif isinstance(section, dict):
            valid = keys_per_section[section_name]
            errors += [f'[{section_name}] unknown key "{key}"' + did_you_mean(key, valid)
                       + valid_elsewhere(key, section_name) for key in section if key not in valid]
    return errors

def _verify_dict_is_param_dict(d: Any) -> None:
    """ Checks that the input is of type ParamDict. """
    if not isinstance(d, dict):
        raise ValueError(f'Expected a dict, but got {d} of type {type(d)}.')
    for key in d:
        if not isinstance(key, str):
            raise ValueError(f'Expected a string as key, but got {key} of type {type(key)}.')

def _verify_dict_is_full_config(d: Dict[str, Any]) -> None:
    """ Checks that dict is of type FullConfig. """
    # Checks that no keys outside a section
    # This is, that d is of type Dict[str, Union[Dict, List]]
    for section_name, section in d.items():
        if not isinstance(section, (dict, list)):
            raise ValueError(f'Key "{section_name}" found outside of a section.')

    # Checks that entries are ParamDicts or List[ParamDict]
    for section in d.values():
        if isinstance(section, list):
            for entry in section:
                _verify_dict_is_param_dict(entry)
        else:
            _verify_dict_is_param_dict(section)

def _verify_restrictions_on_default_and_config(cfg_inp: Dict[str, Any], cfg_def: Dict[str, Any], match_key: Dict[str, str]) -> None:
    """ Checks that the restrictions described in the docstring of merge_config_with_default are met. """
    # Checks that type of cfg_def dict is FullConfig
    _verify_dict_is_full_config(cfg_def)

    # Checks that keys listed in match_key are lists and all other are dicts
    for section_name, section in cfg_def.items():
        if section_name in match_key and not isinstance(section, list):
            raise ValueError(f'"{section_name}" is in match_key so it should be a list in the default config.')
        if section_name not in match_key and not isinstance(section, dict):
            raise ValueError(f'"{section_name}" is not in match_key so it should be a dict in the default config.')

    # Checks that no sections in config that are not in default
    unknown_sections = set(cfg_inp.keys()) - set(cfg_def.keys())
    if unknown_sections:
        raise ValueError('Unknown sections were found in the config file: ' + str(unknown_sections)
                         + '. Please refer to the default config file for valid keys.')

    # Checks that all sections listed in match_key are in default and config
    unmatched_sections = set(match_key.keys()) - set(cfg_def.keys())
    if unmatched_sections:
        raise ValueError('Sections ' + str(unmatched_sections) + ' found in match_key '
                            + 'that are not in the default config.')

    unmatched_sections = set(match_key.keys()) - set(cfg_inp.keys())
    if unmatched_sections:
        raise ValueError('Sections ' + str(unmatched_sections) + ' found in match_key '
                            + 'that are not in the config.')

    # Checks type of config dict
    _verify_dict_is_full_config(cfg_inp)

def _apply_default_values(cfg_inp: FullConfig, cfg_def: FullConfig, match_key: Dict[str, str]) -> FullConfig:
    """ Fills in the default values where the input config does not specify a value. """
    output: FullConfig = {}
    for section_name, section in cfg_def.items():
        if isinstance(section, list):
            key = match_key[section_name]
            output[section_name] = []
            for entry in cfg_inp[section_name]:
                # Finds matching section through match_key in cfg_def
                for default_entry in section:
                    if default_entry[key] == entry[key]:
                        output[section_name].append(copy.deepcopy(default_entry))
                        break
                else:
                    raise ValueError(f'No matching section with same "{section_name}.{key}"="{entry[key]}" found in defaults.')
                # Updates config values in output
                unknown_keys = set(entry.keys()) - set(output[section_name][-1].keys())
                if unknown_keys:
                    raise ValueError(f'Unknown keys {unknown_keys} found in section "{section_name}". '
                                     'All valid keys have to be in the default config.')
                output[section_name][-1].update(entry)
        else:
            entry = cfg_inp.get(section_name, {})
            output[section_name] = copy.deepcopy(section)
            unknown_keys = set(entry.keys()) - set(output[section_name].keys())
            if unknown_keys:
                raise ValueError(f'Unknown keys {unknown_keys} found in section "{section_name}". '
                                 'All valid keys have to be in the default config.')
            output[section_name].update(entry)

    return output

def _replace_none(d: ParamDict) -> None:
    """ Replace '<none>' by None in a ParamDict. This also works inside lists. """
    for key, value in d.items():
        if value == NONE:
            d[key] = None
        elif isinstance(value, list):
            for i, v in enumerate(value):
                if v == NONE:
                    value[i] = None

def _verify_all_mandatory_fields_present(d: ParamDict, section_name: str) -> None:
    """ Verifies that all fields with "<no default>" have been replaced after reading in the config. """
    for key, value in d.items():
        if value == MANDATORY:
            raise ValueError(f'"{key}" in section "{section_name}" is mandatory and was left empty.')

def _resolve_references(d: ParamDict, section_name: str, output: FullConfig) -> None:
    """ Resolve all references of type "<section.key>" in a ParamDict. """
    for key, value in d.items():
        if is_reference(value):
            ref_key = value[1:-1].split('.')
            if len(ref_key) != 2:
                raise ValueError(f'Invalid reference "{value}" in section "{section_name}".')
            if isinstance(output[ref_key[0]], list):
                raise ValueError(f'Invalid reference "{value}" to listed section "{section_name}".')

            referenced_val = output[ref_key[0]][ref_key[1]]
            if is_reference(referenced_val):
                raise ValueError(f'"{ref_key[1]}" in section "{ref_key[0]}" is a reference itself.')
            d[key] = referenced_val

# type hints currently not supported by sphinx autodoc
# def merge_config_with_default(cfg_inp: Dict[str, Any], cfg_def: Dict[str, Any],
#                               match_key: Dict[str, str] = {}) -> FullConfig:
def merge_config_with_default(cfg_inp, cfg_def, match_key={}):
    """
    Merge a TOML config dict with a default TOML dict.
    The default dict dictates the structure of the input:

    - Only sections and keys in the default are allowed in the input
    - All sections listed in match_key must be lists of dicts in the default
      and can be lists of dicts or dicts in the config

    The dicts allows for the following extensions:

    - Mandatory inputs for all calculations indicated by "<no default>"
    - None indicated by "<none>". Also works inside lists
    - References within the dictionary indicated by "<section.key>"

    Parameters
    ----------
    cfg_inp : dict
        The input config dict
    cfg_def : dict
        The default config dict
    match_key : dict, optional
        A dictionary that contains section/key pairs to map entries in listed sections
        between the input and default config.

    Returns
    -------
    dict
        The merged config dict
    """

    # Reports all unknown sections and keys at once, with suggestions for typos
    unknown = find_unknown_entries(cfg_inp, cfg_def, match_key)
    if unknown:
        raise ValueError('Invalid config:\n  ' + '\n  '.join(unknown)
                         + '\nAll valid keys have to be in the default config, see "solid_dmft --params".')

    # Check restrictions and makes sure that config and default are of type FullConfig
    _verify_restrictions_on_default_and_config(cfg_inp, cfg_def, match_key)

    # Checks that keys not listed in match_key are dicts
    # The others can be lists or dicts. This differs from cfg_def
    # to allow users to use multiple sections or not
    for section_name, section in cfg_inp.items():
        if section_name in match_key and not isinstance(section, list):
            cfg_inp[section_name] = [section]
        if section_name not in match_key and not isinstance(section, dict):
            raise ValueError(f'"{section_name}" should be a dict and not a list in the config.')

    # Merges config with default
    output = _apply_default_values(cfg_inp, cfg_def, match_key)

    # Converts "<none>" to None, checks that no mandatory fields were left empty
    # and resolves referencing defaults
    for section_name, section in output.items():
        if isinstance(section, dict):
            _replace_none(section)
            _verify_all_mandatory_fields_present(section, section_name)
            _resolve_references(section, section_name, output)
        else:
            for entry in section:
                _replace_none(entry)
                _verify_all_mandatory_fields_present(entry, section_name)
                _resolve_references(entry, section_name, output)

    return output
