"""
Command line modes of solid_dmft that only read input and never start a
calculation:

- ``solid_dmft --params [section] [key] [--solver-type TYPE] [--json]`` lists the
  input parameters with type, code default and description
- ``solid_dmft --check [config.toml] [--json]`` reads and verifies a config file
  and prints it with all defaults filled in

This module must not import triqs: these modes run without MPI, and importing
triqs.utility.mpi writes to stdout, which would break the --json output.
"""

import argparse
import copy
import json
import sys

try:
    import tomllib
except ImportError:
    import tomli as tomllib

from solid_dmft.io_tools import documentation, postproc_toml_dict as pp, verify_input_params


def is_cli_call(argv):
    """True if argv asks for one of the modes here. A plain ``solid_dmft [config]`` runs the calculation."""
    return any(arg.startswith('-') for arg in argv[1:])


def _format_value(value):
    """TOML spelling of a value; None is not representable in TOML and printed as None."""
    if value is None:
        return 'None'
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, list):
        return '[' + ', '.join(_format_value(v) for v in value) + ']'
    if isinstance(value, dict):
        return '{' + ', '.join(f'{k} = {_format_value(v)}' for k, v in value.items()) + '}'
    return str(value)


def _format_default(param):
    if param['mandatory']:
        return 'mandatory'
    if param['default_ref'] is not None:
        return f'value of {param["default_ref"]}'
    return _format_value(param['default'])


def _print_params(params):
    for param in params:
        header = f'[{param["section"]}]'
        if 'solver_type' in param:
            header += f' (type = "{param["solver_type"]}")'
        print(f'{header} {param["key"]}')
        print(f'    type:    {param["type"]}')
        print(f'    default: {_format_default(param)}')
        for line in (param['description'] or '').splitlines():
            print(f'    {line}' if line.strip() else '')
        print('')


def params_mode(selection, solver_type=None, as_json=False):
    """Prints the parameters selected by [section] [key]. Returns the exit code."""
    defaults = pp.load_default_config()
    sections = list(defaults)
    solver_types = [entry[pp.MATCH_KEY['solver']] for entry in defaults['solver']]

    if len(selection) > 2:
        print('error: --params takes at most a section and a key', file=sys.stderr)
        return 2
    section = key = None
    if len(selection) == 2:
        section, key = selection
    elif len(selection) == 1:
        # a single argument is a section if it names one, otherwise a key searched in all sections
        if selection[0] in sections:
            section = selection[0]
        else:
            key = selection[0]

    if section is not None and section not in sections:
        print(f'error: unknown section "{section}"' + pp.did_you_mean(section, sections), file=sys.stderr)
        return 1
    if solver_type is not None and solver_type not in solver_types:
        print(f'error: unknown solver type "{solver_type}"' + pp.did_you_mean(solver_type, solver_types),
              file=sys.stderr)
        return 1

    params = documentation.list_params(section=section, key=key, solver_type=solver_type)
    if not params:
        candidates = {p['key'] for p in documentation.list_params(section=section, solver_type=solver_type)}
        where = f' in section "{section}"' if section else ''
        if solver_type is not None:
            where += f' for solver type "{solver_type}"'
        elsewhere = sorted({p['section'] for p in documentation.list_params(key=key)})
        hint = f' ("{key}" is in {", ".join(f"[{s}]" for s in elsewhere)})' if elsewhere else ''
        print(f'error: unknown parameter "{key}"{where}' + pp.did_you_mean(key, candidates) + hint, file=sys.stderr)
        return 1

    if as_json:
        print(json.dumps(params, indent=2, default=str))
    else:
        _print_params(params)
    return 0


def check_config(config_file):
    """
    Reads config_file, merges it with the defaults and runs verify_before_dmft_cycle.

    Returns
    -------
    full_params : dict or None
        the config with all defaults filled in, None if it is invalid
    errors : list of str
    """
    try:
        with open(config_file, 'rb') as file:
            params = tomllib.load(file)
    except (OSError, tomllib.TOMLDecodeError) as error:
        return None, [f'cannot read {config_file}: {error}']

    defaults = pp.load_default_config()
    errors = pp.find_unknown_entries(params, defaults)
    if errors:
        return None, errors

    try:
        full_params = pp.merge_config_with_default(copy.deepcopy(params), defaults, pp.MATCH_KEY)
        verify_input_params.verify_before_dmft_cycle(full_params)
    except (ValueError, TypeError, KeyError, IndexError) as error:
        return None, [str(error)]
    return full_params, []


def _print_config(full_params, params):
    """Prints the resolved config; keys not set in the config file are marked as defaults."""
    for section_name, section in full_params.items():
        given = params.get(section_name, {})
        entries = section if isinstance(section, list) else [section]
        given_entries = given if isinstance(given, list) else [given]
        for i, entry in enumerate(entries):
            print(f'[{section_name}]' if not isinstance(section, list) else f'[[{section_name}]]')
            set_keys = given_entries[i] if i < len(given_entries) else {}
            for key, value in entry.items():
                marker = '' if key in set_keys else '  (default)'
                print(f'    {key: <34} {_format_value(value)}{marker}')
            print('')


def check_mode(config_file, as_json=False):
    """Prints the result of check_config. Returns the exit code: 0 if valid, 1 otherwise."""
    full_params, errors = check_config(config_file)
    if as_json:
        print(json.dumps({'config_file': config_file, 'valid': not errors,
                          'errors': errors, 'config': full_params}, indent=2, default=str))
        return 1 if errors else 0

    if errors:
        print(f'{config_file}: invalid config')
        for error in errors:
            print(f'  error: {error}')
        return 1

    with open(config_file, 'rb') as file:
        params = tomllib.load(file)
    _print_config(full_params, params)
    print(f'{config_file}: config is valid')
    return 0


def main(args):
    """Entry point for the arguments after the program name. Returns the exit code."""
    parser = argparse.ArgumentParser(
        prog='solid_dmft',
        description='Without options, "solid_dmft [config]" runs the calculation defined in config '
                    '(default: dmft_config.toml). The options below only read input.')
    parser.add_argument('config', nargs='?', default='dmft_config.toml',
                        help='config file (default: dmft_config.toml)')
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--params', nargs='*', metavar='SECTION_OR_KEY',
                      help='list input parameters with type, default and description: '
                           '[section] [key], or a key searched in all sections')
    mode.add_argument('--check', action='store_true',
                      help='read and verify the config file without running a calculation, '
                           'and print it with the defaults filled in')
    parser.add_argument('--solver-type', help='with --params: only list solver keys valid for this solver type')
    parser.add_argument('--json', action='store_true', help='machine-readable output')
    args = parser.parse_args(args)

    if args.params is not None:
        return params_mode(args.params, solver_type=args.solver_type, as_json=args.json)
    if args.solver_type is not None:
        parser.error('--solver-type is only valid with --params')
    return check_mode(args.config, as_json=args.json)
