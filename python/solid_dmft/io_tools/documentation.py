"""
Reads the input parameter documentation (documentation.txt) and the code
defaults (default.toml), and checks that both describe the same parameters.

documentation.txt is also turned into the web documentation by
doc/input_output/DMFT_input/generate_doc_from_comments.sh, so its format is
fixed by that script:

- a section starts with ``[  name  ]`` followed by a line of dashes
- in listed sections (solver), a subsection for one solver type starts with
  a line whose first word is the type, followed by a line of ``=``; entries
  before the first subsection apply to all solver types
- an entry starts in the first column as ``name : type, default = value`` or
  ``name : type, mandatory``, followed by the indented description

default.toml uses "<no default>" for mandatory keys, "<none>" for None and
"<section.key>" for a default taken from another key.
"""

import ast
import math
import os
import re
import textwrap
from dataclasses import dataclass, field
from typing import List, Optional

try:
    import tomllib
except ImportError:
    import tomli as tomllib

_IO_TOOLS_DIR = os.path.dirname(os.path.abspath(__file__))
DOCUMENTATION_FILE = os.path.join(_IO_TOOLS_DIR, 'documentation.txt')
DEFAULT_CONFIG_FILE = os.path.join(_IO_TOOLS_DIR, 'default.toml')

# Sections that are lists of tables in default.toml, with the key selecting the entry
LISTED_SECTIONS = {'solver': 'type'}

MANDATORY = '<no default>'
NONE = '<none>'

_SECTION_RE = re.compile(r'^\[\s+(\w+)\s+\]\s*$')
# the sed script only recognizes entries with a blank before the colon
_ENTRY_RE = re.compile(r'^(\w+) +:\s*(.*?)\s*$')
_SPEC_RE = re.compile(r'^(?P<type>[^,]*?)\s*'
                      r'(?:,\s*(?:(?P<mandatory>mandatory)|default\s*=\s*(?P<default>.+?)))?\s*$')
_REFERENCE_RE = re.compile(r'^[A-Za-z_]\w*\.\w+$')

# Python types accepted for the first word of each alternative in a documented type
_TYPE_WORDS = {
    'bool': (bool,),
    'int': (int,),
    'float': (float, int),
    'str': (str,),
    'string': (str,),
    'list': (list,),
    'dict': (dict,),
    'table': (dict,),
}


@dataclass
class ParamDoc:
    """One entry of documentation.txt. solver_type is None for entries valid in all solver types."""
    section: str
    name: str
    type: str
    mandatory: bool
    default: Optional[str]
    solver_type: Optional[str]
    line: int
    description_lines: List[str] = field(default_factory=list, repr=False)

    @property
    def description(self):
        return textwrap.dedent('\n'.join(self.description_lines)).strip()


def parse_documentation(path=DOCUMENTATION_FILE):
    """
    Parse documentation.txt.

    Returns
    -------
    entries : list of ParamDoc
    problems : list of str
        lines that could not be parsed
    """
    with open(path, 'r', encoding='utf-8') as file:
        lines = file.read().splitlines()

    entries, problems = [], []
    section = solver_type = current = None
    for i, line in enumerate(lines):
        next_line = lines[i + 1] if i + 1 < len(lines) else ''
        match = _SECTION_RE.match(line)
        if match:
            section, solver_type, current = match.group(1), None, None
            continue
        if re.fullmatch(r'[-=]+\s*', line):
            continue
        if line.strip() and not line[0].isspace() and re.fullmatch(r'=+\s*', next_line):
            solver_type, current = line.split()[0], None
            continue
        if line.strip() and not line[0].isspace():
            current = None
            match = _ENTRY_RE.match(line)
            spec = _SPEC_RE.match(match.group(2)) if match else None
            if section is None or spec is None:
                problems.append(f'documentation.txt line {i+1}: cannot parse "{line}"')
                continue
            current = ParamDoc(section=section, name=match.group(1), type=spec.group('type'),
                               mandatory=spec.group('mandatory') is not None,
                               default=spec.group('default'), solver_type=solver_type, line=i + 1)
            entries.append(current)
            continue
        if current is not None:
            current.description_lines.append(line)

    return entries, problems


def load_defaults(path=DEFAULT_CONFIG_FILE):
    with open(path, 'rb') as file:
        return tomllib.load(file)


def parse_documented_default(text):
    """
    Convert the default written in documentation.txt to a value comparable with default.toml.
    A reference to another key ("general.J") is returned as "<general.J>".
    """
    if _REFERENCE_RE.match(text):
        return f'<{text}>'
    value = ast.literal_eval(text)
    return NONE if value is None else value


def _values_equal(a, b):
    if isinstance(a, bool) or isinstance(b, bool):
        return type(a) is type(b) and a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return math.isclose(a, b, rel_tol=1e-12, abs_tol=0)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_values_equal(x, y) for x, y in zip(a, b))
    return a == b


def _type_alternatives(doc_type):
    """Python types for each alternative of e.g. "int or list of int"; None for an unknown word."""
    alternatives = []
    for alternative in re.split(r'\s+or\s+', doc_type.strip()):
        words = alternative.split()
        alternatives.append(_TYPE_WORDS.get(words[0]) if words else None)
    return alternatives


def _type_matches(doc_type, value):
    for types in _type_alternatives(doc_type):
        # bool is a subclass of int, but an int/float parameter must not accept true/false
        if isinstance(value, bool) and bool not in types:
            continue
        if isinstance(value, types):
            return True
    return False


def _is_reference(value):
    return isinstance(value, str) and value.startswith('<') and value.endswith('>') and value not in (MANDATORY, NONE)


def _doc_index(entries):
    """Maps (section, solver_type, name) to ParamDoc; solver_type is None outside listed sections."""
    index, duplicates = {}, []
    for entry in entries:
        key = (entry.section, entry.solver_type, entry.name)
        if key in index:
            duplicates.append(entry)
        else:
            index[key] = entry
    return index, duplicates


def _code_keys(defaults):
    """Yields (section, solver_type, key, default value) for every key the code accepts."""
    for section_name, section in defaults.items():
        if section_name in LISTED_SECTIONS:
            match_key = LISTED_SECTIONS[section_name]
            for entry in section:
                for key, value in entry.items():
                    yield section_name, entry[match_key], key, value
        else:
            for key, value in section.items():
                yield section_name, None, key, value


def find_doc(index, section, solver_type, key):
    """Entry specific to the solver type if there is one, else the entry common to all types."""
    return index.get((section, solver_type, key)) or index.get((section, None, key))


def check_consistency(doc_path=DOCUMENTATION_FILE, default_path=DEFAULT_CONFIG_FILE):
    """
    Compare documentation.txt with default.toml.

    Returns a list of problems, empty if both agree: every key in default.toml
    has an entry in documentation.txt and vice versa, every entry states a
    default or "mandatory" that matches default.toml, and every documented type
    is known and fits the code default.
    """
    entries, problems = parse_documentation(doc_path)
    defaults = load_defaults(default_path)
    index, duplicates = _doc_index(entries)

    for entry in duplicates:
        where = f'[{entry.section}]' + (f' ({entry.solver_type})' if entry.solver_type else '')
        problems.append(f'{where} {entry.name}: documented more than once (line {entry.line})')

    doc_sections = {entry.section for entry in entries}
    for section in sorted(set(defaults) - doc_sections):
        problems.append(f'[{section}]: section in default.toml but not in documentation.txt')
    for section in sorted(doc_sections - set(defaults)):
        problems.append(f'[{section}]: section in documentation.txt but not in default.toml')

    # Keys present in the code but not documented
    code_index = {}
    for section, solver_type, key, value in _code_keys(defaults):
        code_index[(section, solver_type, key)] = value
        where = f'[{section}]' + (f' ({solver_type})' if solver_type else '')
        if section in doc_sections and find_doc(index, section, solver_type, key) is None:
            problems.append(f'{where} {key}: in default.toml but not in documentation.txt')

    # Keys documented but not present in the code
    code_solver_types = {section: {entry[match_key] for entry in defaults.get(section, [])}
                         for section, match_key in LISTED_SECTIONS.items()}
    for (section, solver_type, key), entry in index.items():
        if section not in defaults:
            continue
        if section in LISTED_SECTIONS:
            if solver_type is None:
                types = sorted(code_solver_types[section])
                if any((section, t, key) in index for t in types):
                    problems.append(f'[{section}] {key}: documented for all types and for a single type')
                missing = [t for t in types if (section, t, key) not in code_index]
            elif solver_type not in code_solver_types[section]:
                problems.append(f'[{section}] ({solver_type}): subsection in documentation.txt '
                                f'but no such type in default.toml (line {entry.line})')
                continue
            else:
                missing = [solver_type] if (section, solver_type, key) not in code_index else []
            for t in missing:
                problems.append(f'[{section}] ({t}) {key}: in documentation.txt but not in default.toml')
        elif (section, None, key) not in code_index:
            problems.append(f'[{section}] {key}: in documentation.txt but not in default.toml')

    # Defaults and types
    for entry in entries:
        if None in _type_alternatives(entry.type):
            problems.append(f'[{entry.section}] {entry.name}: unknown type "{entry.type}" (line {entry.line})')
    for (section, solver_type, key), code_value in code_index.items():
        doc = find_doc(index, section, solver_type, key)
        if doc is None:
            continue
        where = f'[{section}]' + (f' ({solver_type})' if solver_type else '') + f' {key}'
        # the match key holds the solver type itself, there is no default to compare
        if LISTED_SECTIONS.get(section) == key:
            continue
        if doc.mandatory:
            if code_value != MANDATORY:
                problems.append(f'{where}: documented as mandatory, but default.toml has {code_value!r}')
            continue
        if doc.default is None:
            problems.append(f'{where}: documentation.txt states neither a default nor "mandatory" (line {doc.line})')
            continue
        try:
            doc_value = parse_documented_default(doc.default)
        except (ValueError, SyntaxError):
            problems.append(f'{where}: cannot parse documented default "{doc.default}" (line {doc.line})')
            continue
        if not _values_equal(doc_value, code_value):
            problems.append(f'{where}: documented default {doc.default} differs from default.toml {code_value!r}')
        elif code_value not in (MANDATORY, NONE) and not _is_reference(code_value) \
                and None not in _type_alternatives(doc.type) and not _type_matches(doc.type, code_value):
            problems.append(f'{where}: default {code_value!r} does not match documented type "{doc.type}"')

    return problems


def list_params(section=None, key=None, solver_type=None,
                doc_path=DOCUMENTATION_FILE, default_path=DEFAULT_CONFIG_FILE):
    """
    All parameters accepted by the code with their documentation, in the order of default.toml.

    Parameters
    ----------
    section, key : str, optional
        restrict to this section / key
    solver_type : str, optional
        restrict listed sections (solver) to this type

    Returns
    -------
    list of dict
        keys: section, solver_type (listed sections only), key, type, mandatory,
        default (None if mandatory or a reference), default_ref ("section.key"
        if the default is taken from another key, else None), description
    """
    entries, _ = parse_documentation(doc_path)
    index, _ = _doc_index(entries)
    defaults = load_defaults(default_path)

    params = []
    for sec, stype, name, value in _code_keys(defaults):
        if (section is not None and sec != section) or (key is not None and name != key):
            continue
        if solver_type is not None and stype is not None and stype != solver_type:
            continue
        doc = find_doc(index, sec, stype, name)
        mandatory = value == MANDATORY or (LISTED_SECTIONS.get(sec) == name and doc is not None and doc.mandatory)
        param = {'section': sec}
        if stype is not None:
            param['solver_type'] = stype
        param.update({
            'key': name,
            'type': doc.type if doc else None,
            'mandatory': mandatory,
            'default': None if mandatory or value == NONE or _is_reference(value) else value,
            'default_ref': value[1:-1] if _is_reference(value) else None,
            'description': doc.description if doc else None,
        })
        params.append(param)
    return params
