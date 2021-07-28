# -*- coding: utf-8 -*-
################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2016-2018, N. Wentzell
# Copyright (C) 2018-2019, Simons Foundation
#   author: N. Wentzell
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
"""
Updates the config file (usually dmft_config.ini) for continuing old calculations.
"""

import os.path
import shutil
import sys
from configparser import ConfigParser

def _backup_old_file(path):
    """ Copies file to same folder with the prefix "backup_". """
    directory, filename = os.path.split(path)
    shutil.copy2(path, os.path.join(directory, 'backup_'+filename))

def _load_config_file(path):
    """ Loads file with the configparser module. Returns ConfigParser object. """
    config = ConfigParser()
    config.read(path)
    return config

def _update_section_names(config):
    """
    Applies the mapping between legacy names of sections and new names. The
    mapping is saved in LEGACY_SECTION_NAME_MAPPING.

    Parameters
    ----------
    config : ConfigParser
        A configparser instance that has read the config file already.

    Returns
    -------
    config : ConfigParser
        The configparser with correct section names.
    """

    LEGACY_SECTION_NAME_MAPPING = {'solver': 'solver_parameters', 'advanced': 'advanced_parameters'}

    for new_name, legacy_name in LEGACY_SECTION_NAME_MAPPING.items():
        # Only new name in there, everything is okay
        if legacy_name not in config.keys():
            continue

        # Only legacy name exists, creates updated section
        if new_name not in config.keys():
            config.add_section(new_name)

        # Transfers parameters in legacy section to updated section
        for param_name, param_value in config[legacy_name].items():
            config.set(new_name, param_name, param_value)
        config.remove_section(legacy_name)

    return config

def _update_params(config):
    """
    Updates parameter names/values and adds new required parameters.

    Parameters
    ----------
    config : ConfigParser
        The config parser with correct section names.

    Returns
    -------
    config : ConfigParser
        The completely updated config parser.
    """

    # Updates h_int_type
    if config['general']['h_int_type'] in ('1', '2', '3'):
        config['general']['h_int_type'] = {'1': 'density_density',
                                           '2': 'kanamori',
                                           '3': 'full_slater'
                                          }[config['general']['h_int_type']]

    # ---new params
    # Updates solver_type - if not existent, uses previous default cthyb
    if 'solver_type' not in config['general']:
        config['general']['solver_type'] = 'cthyb'

    if 'dft_code' not in config['dft']:
        config['dft']['dft_code'] = 'vasp'

    # ---updated params
    # Updates legendre coefficients
    if 'n_LegCoeff' in config['solver']:
        config['general']['n_l'] = config['solver']['n_LegCoeff']
        del config['solver']['n_LegCoeff']

    # Updates dft_executable
    if 'executable' in config['dft']:
        config['dft']['dft_exec'] = config['dft']['executable']
        del config['dft']['executable']

    # Updates dft_executable
    if 'wannier90_exec' in config['dft']:
        config['dft']['w90_exec'] = config['dft']['wannier90_exec']
        del config['dft']['wannier90_exec']

    return config

def _write_config_file(config, path):
    """ Writes config parser content to a file. """
    with open(path, 'w') as file:
        config.write(file)

def main(path='dmft_config.ini'):
    """ Combines methods in the full work flow for updating the config file. """
    if not os.path.isfile(path):
        raise ValueError('File {} does not exist.'.format(path))

    _backup_old_file(path)

    config = _load_config_file(path)
    config = _update_section_names(config)
    config = _update_params(config)
    _write_config_file(config, path)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    elif len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        raise TypeError('Maximally one argument supported: the config file name')
