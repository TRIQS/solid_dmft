#!@TRIQS_PYTHON_EXECUTABLE@
################################################################################
#
# solid_dmft - A versatile python wrapper to perform DFT+DMFT calculations
#              utilizing the TRIQS software library
#
# Copyright (C) 2018-2020, ETH Zurich
# Copyright (C) 2021, The Simons Foundation
#      authors: A. Hampel, M. Merkel, and S. Beck
#
# solid_dmft is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# solid_dmft is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# solid_dmft (in the file COPYING.txt in this directory). If not, see
# <http://www.gnu.org/licenses/>.
#
################################################################################
"""
This python script allows one to perform DFT+DMFT calculations with VASP
or with a pre-defined h5 archive (only one-shot) for
multiband/many-correlated-shells systems using the TRIQS package,
in combination with the CThyb solver and SumkDFT from DFT-tools.
triqs version 3.0 or higher is required
"""

# system
import os
import sys
import shutil
from timeit import default_timer as timer

try:
    import tomllib
except ImportError:
    import tomli as tomllib

# triqs
import triqs.utility.mpi as mpi

# own modules
from solid_dmft.dmft_cycle import dmft_cycle
from solid_dmft.csc_flow import csc_flow_control
from solid_dmft.io_tools import postproc_toml_dict, verify_input_params

def run_dmft(params, config_file_name=None):
    # timing information
    if mpi.is_master_node():
        global_start = timer()

    # Reads in default params file and merges with params
    full_params = None
    default_config_name = os.path.join(os.path.dirname(__file__),
                                       'io_tools', 'default.toml')
    if mpi.is_master_node():
        with open(default_config_name, 'rb') as file:
            default_params = tomllib.load(file)
        full_params = postproc_toml_dict.merge_config_with_default(params, default_params,
                                                                   {'solver': 'type'})
        verify_input_params.verify_before_dmft_cycle(full_params)
    full_params = mpi.bcast(full_params)

    # Prints parameters
    if mpi.is_master_node():
        # Prints only sections that are in the input file
        for section_name in params:
            section = full_params[section_name]
            print(f'\n{section_name} parameters')
            if isinstance(section, dict):
                for key, value in section.items():
                    print(f'    {key: <30} {str(value)}')
            else:
                for i, entry in enumerate(section):
                    print(f'entry {i+1}')
                    for key, value in entry.items():
                        print(f'    {key: <30} {str(value)}')
        print('')
    mpi.barrier()

    general_params = full_params['general']
    solver_params = full_params['solver']
    dft_params = full_params['dft']
    advanced_params = full_params['advanced']

    if general_params['csc']:
        # Start CSC calculation, always in same folder as dmft_config
        general_params['jobname'] = '.'
        csc_flow_control(general_params, solver_params, dft_params, advanced_params)
    else:
        # Sets up one-shot calculation
        mpi.report('', '#'*80)
        mpi.report(f'Using input file "{general_params["seedname"]}.h5" '
                   + f'and running in folder "{general_params["jobname"]}".')

        if mpi.is_master_node():
            # Checks for h5 file
            if not os.path.exists(general_params['seedname']+'.h5'):
                raise FileNotFoundError('Input h5 file not found')

            # Creates output directory if it does not exist
            if not os.path.exists(general_params['jobname']):
                os.makedirs(general_params['jobname'])

            # Copies h5 archive to subfolder if it is not there
            h5_name = general_params['seedname']+'.h5'
            if not os.path.isfile(general_params['jobname']+'/'+os.path.basename(h5_name)):
                shutil.copyfile(h5_name, general_params['jobname']+'/'+os.path.basename(h5_name))

            # Copies config file to subfolder
            if config_file_name is not None:
                shutil.copyfile(config_file_name, general_params['jobname']+'/'+os.path.basename(config_file_name))
        mpi.barrier()

        # Runs dmft_cycle
        dmft_cycle(general_params, solver_params, advanced_params,
                   dft_params, general_params['n_iter_dmft'])

    mpi.barrier()
    if mpi.is_master_node():
        global_end = timer()
        print('-'*80)
        print('overall elapsed time: %10.4f seconds'%(global_end-global_start))

def main(argv=sys.argv):
    """The main function for one-shot and charge-self-consistent calculations"""
    config_file_name = str(argv[1]) if len(argv) > 1 else 'dmft_config.toml'

    params = None
    if mpi.is_master_node():
        print('Reading the config file', config_file_name)
        with open(config_file_name, 'rb') as file:
            params = tomllib.load(file)
    params = mpi.bcast(params)

    run_dmft(params, config_file_name)

if __name__ == '__main__':
    main(sys.argv)
