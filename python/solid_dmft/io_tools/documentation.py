"""
Reads the input parameter documentation (documentation.txt), checks it against
the code defaults (default.toml, via postproc_toml_dict) and renders it as the
reStructuredText pages of the web documentation.

Format of documentation.txt:

- a section starts with ``[  name  ]`` followed by a line of dashes
- in listed sections (solver), a subsection for one solver type starts with
  a line whose first word is the type, followed by a line of ``=``; entries
  before the first subsection apply to all solver types
- an entry starts in the first column as ``name : type, default = value`` or
  ``name : type, mandatory``, followed by the indented description

Run ``python -m solid_dmft.io_tools.documentation <output_dir>`` to write the
rst pages (done by the doc build).
"""

import ast
import math
import os
import re
import sys
import textwrap
from dataclasses import dataclass, field
from typing import List, Optional

from solid_dmft.io_tools import postproc_toml_dict as pp

DOCUMENTATION_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'documentation.txt')

_SECTION_RE = re.compile(r'^\[\s+(\w+)\s+\]\s*$')
_ENTRY_RE = re.compile(r'^(\w+) +:\s*(.*?)\s*$')
_SPEC_RE = re.compile(r'^(?P<type>[^,]*?)\s*'
                      r'(?:,\s*(?:(?P<mandatory>mandatory)|default\s*=\s*(?P<default>.+?)))?\s*$')
_REFERENCE_RE = re.compile(r'^[A-Za-z_]\w*\.\w+$')


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


@dataclass
class Section:
    """A section of documentation.txt: its entries and, for listed sections, the subsection titles."""
    name: str
    items: list = field(default_factory=list)   # ParamDoc or ('subsection', title)


def parse_documentation(path=DOCUMENTATION_FILE):
    """
    Parse documentation.txt.

    Returns
    -------
    sections : list of Section, in file order
    problems : list of str, lines that could not be parsed
    """
    with open(path, 'r', encoding='utf-8') as file:
        lines = file.read().splitlines()

    sections, problems = [], []
    solver_type = current = None
    for i, line in enumerate(lines):
        next_line = lines[i + 1] if i + 1 < len(lines) else ''
        match = _SECTION_RE.match(line)
        if match:
            sections.append(Section(match.group(1)))
            solver_type = current = None
            continue
        if re.fullmatch(r'[-=]+\s*', line):
            continue
        if line.strip() and not line[0].isspace() and re.fullmatch(r'=+\s*', next_line):
            solver_type, current = line.split()[0], None
            if sections:
                sections[-1].items.append(('subsection', line.strip()))
            continue
        if line.strip() and not line[0].isspace():
            current = None
            match = _ENTRY_RE.match(line)
            spec = _SPEC_RE.match(match.group(2)) if match else None
            if not sections or spec is None:
                problems.append(f'documentation.txt line {i+1}: cannot parse "{line}"')
                continue
            current = ParamDoc(section=sections[-1].name, name=match.group(1), type=spec.group('type'),
                               mandatory=spec.group('mandatory') is not None,
                               default=spec.group('default'), solver_type=solver_type, line=i + 1)
            sections[-1].items.append(current)
            continue
        if current is not None:
            current.description_lines.append(line)

    return sections, problems


def _entries(sections):
    return [item for section in sections for item in section.items if isinstance(item, ParamDoc)]


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


def find_doc(index, section, solver_type, key):
    """Entry specific to the solver type if there is one, else the entry common to all types."""
    return index.get((section, solver_type, key)) or index.get((section, None, key))


def _documented_value(text):
    """The documented default as a value comparable with default.toml ("general.J" -> "<general.J>")."""
    if _REFERENCE_RE.match(text):
        return f'<{text}>'
    value = ast.literal_eval(text)
    return pp.NONE if value is None else value


def _values_equal(a, b):
    if isinstance(a, bool) or isinstance(b, bool):
        return type(a) is type(b) and a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return math.isclose(a, b, rel_tol=1e-12, abs_tol=0)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_values_equal(x, y) for x, y in zip(a, b))
    return a == b


def check_consistency(doc_path=DOCUMENTATION_FILE, default_path=pp.DEFAULT_CONFIG_FILE):
    """
    Compare documentation.txt with default.toml.

    Returns a list of problems, empty if both agree: every key of default.toml has an
    entry in documentation.txt and vice versa, and every entry states a default or
    "mandatory" that matches default.toml.
    """
    sections, problems = parse_documentation(doc_path)
    defaults = pp.load_default_config(default_path)
    index, duplicates = _doc_index(_entries(sections))

    for entry in duplicates:
        where = f'[{entry.section}]' + (f' ({entry.solver_type})' if entry.solver_type else '')
        problems.append(f'{where} {entry.name}: documented more than once (line {entry.line})')

    doc_sections = {section.name for section in sections}
    for section in sorted(set(defaults) ^ doc_sections):
        where = 'default.toml' if section in defaults else 'documentation.txt'
        problems.append(f'[{section}]: section only in {where}')

    # every key of the code is documented, with the same default
    code_index = {}
    for section, solver_type, key, code_value in pp.iter_default_keys(defaults):
        code_index[(section, solver_type, key)] = code_value
        where = f'[{section}]' + (f' ({solver_type})' if solver_type else '') + f' {key}'
        doc = find_doc(index, section, solver_type, key)
        if doc is None:
            if section in doc_sections:
                problems.append(f'{where}: in default.toml but not in documentation.txt')
            continue
        # the match key holds the solver type itself, there is no default to compare
        if pp.MATCH_KEY.get(section) == key:
            continue
        if doc.mandatory or code_value == pp.MANDATORY:
            if not (doc.mandatory and code_value == pp.MANDATORY):
                problems.append(f'{where}: mandatory in only one of documentation.txt and default.toml')
            continue
        if doc.default is None:
            problems.append(f'{where}: documentation.txt states neither a default nor "mandatory" (line {doc.line})')
            continue
        try:
            doc_value = _documented_value(doc.default)
        except (ValueError, SyntaxError):
            problems.append(f'{where}: cannot parse documented default "{doc.default}" (line {doc.line})')
            continue
        if not _values_equal(doc_value, code_value):
            problems.append(f'{where}: documented default {doc.default} differs from default.toml {code_value!r}')

    # every documented key exists in the code
    solver_types = {section: {entry[key] for entry in defaults.get(section, [])}
                    for section, key in pp.MATCH_KEY.items()}
    for (section, solver_type, key), entry in index.items():
        if section not in defaults:
            continue
        if section not in pp.MATCH_KEY:
            types = [None]
        elif solver_type is None:
            types = sorted(solver_types[section])
        elif solver_type in solver_types[section]:
            types = [solver_type]
        else:
            problems.append(f'[{section}] ({solver_type}): subsection without such a type in default.toml '
                            f'(line {entry.line})')
            continue
        for t in types:
            if (section, t, key) not in code_index:
                problems.append(f'[{section}]' + (f' ({t})' if t else '') + f' {key}: '
                                'in documentation.txt but not in default.toml')

    return problems


def list_params(section=None, key=None, solver_type=None,
                doc_path=DOCUMENTATION_FILE, default_path=pp.DEFAULT_CONFIG_FILE):
    """
    All parameters accepted by the code with their documentation, in the order of default.toml.

    Returns a list of dicts with keys section, solver_type (listed sections only), key, type,
    mandatory, default (None if mandatory or a reference), default_ref ("section.key" if the
    default is taken from another key, else None) and description.
    """
    sections, _ = parse_documentation(doc_path)
    index, _ = _doc_index(_entries(sections))

    params = []
    for sec, stype, name, value in pp.iter_default_keys(pp.load_default_config(default_path)):
        if (section is not None and sec != section) or (key is not None and name != key):
            continue
        if solver_type is not None and stype is not None and stype != solver_type:
            continue
        doc = find_doc(index, sec, stype, name)
        mandatory = value == pp.MANDATORY or (pp.MATCH_KEY.get(sec) == name and doc is not None and doc.mandatory)
        param = {'section': sec}
        if stype is not None:
            param['solver_type'] = stype
        param.update({
            'key': name,
            'type': doc.type if doc else None,
            'mandatory': mandatory,
            'default': None if mandatory or value == pp.NONE or pp.is_reference(value) else value,
            'default_ref': value[1:-1] if pp.is_reference(value) else None,
            'description': doc.description if doc else None,
        })
        params.append(param)
    return params


# ---------------------------------------------------------------------------
# reStructuredText pages of the web documentation

_PAGE_HEADERS = {
    'general': ('[general]: General parameters',
                'Frequently used parameters that apply to the whole simulation.'),
    'solver': ('[solver]: solver specific parameters',
               'Here are the parameters that are uniquely dependent on the solver chosen. Some parameters '
               'are used within solid_dmft and some are passed directly into the triqs solver.\n'
               'To see which parameters were passed to the solver for a given calculation, look at the '
               'triqs_solver_params in DMFT_input/solver in the h5 archive.\n'
               'Solver-specific parameters are listed in the respective sections.'),
    'dft': ('[dft]: DFT related inputs',
            'List of parameters for the DFT calculation in charge-self-consistent calculations.\n'
            'The parameters are ignored in one-shot calculations.'),
    'gw': ('[GW]: GW embedding inputs',
           'List of parameters for the GW embedding calculation. The parameters are ignored unless '
           'gw_embedding=true, or interactions are read from cRPA.'),
    'advanced': ('[advanced]: Advanced inputs',
                 'Advanced parameters, do not modify the default value unless you know what you are doing.'),
}


def _rst_entry(entry):
    spec = f'**type** = {entry.type}'
    if entry.mandatory:
        spec += ';  **mandatory**'
    elif entry.default is not None:
        spec += f';  **default** = {entry.default}'
    lines = [f'.. admonition:: {entry.name}', '    :class: intag', '', f'    {spec}', '']
    lines += ['    ' + line if line.strip() else '' for line in entry.description.splitlines()]
    return '\n'.join(lines) + '\n'


def _rst_section_page(section):
    title, intro = _PAGE_HEADERS.get(section.name, (f'[{section.name}]', ''))
    out = [title, '-' * len(title), '', intro, '']
    for item in section.items:
        if isinstance(item, ParamDoc):
            out.append(_rst_entry(item))
        else:
            out += [item[1], '=' * len(item[1]), '']
    return '\n'.join(out) + '\n'


def _rst_input_page(sections):
    out = ['Input', '-------', '',
           'The aim of this section is to provide a comprehensive listing of all the input flags available '
           'for the `dmft_config.toml` input file. We begin by listing the possible sections and follow with '
           'the input parameters.', '',
           '.. toctree::', '    :maxdepth: 1', '']
    out += [f'    {section.name}' for section in sections]
    out += ['', 'Below an exhaustive list containing all the parameters marked by section.', '']
    for section in sections:
        names = [item.name for item in section.items if isinstance(item, ParamDoc)]
        out += [f'**[{section.name}]**', '', '; '.join(names), '']
    return '\n'.join(out) + '\n'


def write_rst(output_dir, doc_path=DOCUMENTATION_FILE):
    """Writes input.rst and one page per section of documentation.txt to output_dir."""
    sections, problems = parse_documentation(doc_path)
    if problems:
        raise ValueError('\n'.join(problems))
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, 'input.rst'), 'w', encoding='utf-8') as file:
        file.write(_rst_input_page(sections))
    for section in sections:
        with open(os.path.join(output_dir, f'{section.name}.rst'), 'w', encoding='utf-8') as file:
            file.write(_rst_section_page(section))


if __name__ == '__main__':
    write_rst(sys.argv[1] if len(sys.argv) > 1 else '.')
