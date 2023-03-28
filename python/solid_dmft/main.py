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

# triqs
import triqs.utility.mpi as mpi

# own modules
from solid_dmft.read_config import read_config
from solid_dmft.dmft_cycle import dmft_cycle
from solid_dmft.csc_flow import csc_flow_control


def main(argv=sys.argv):
    """The main function for one-shot and charge-self-consistent calculations"""
    # timing information
    if mpi.is_master_node():
        global_start = timer()

    # reading configuration for calculation
    general_params = None
    solver_params = None
    dft_params = None
    advanced_params = None
    if len(argv) > 1:
        config_file_name = str(argv[1])
    else:
        config_file_name = 'dmft_config.ini'
    if not os.path.isfile(config_file_name):
        raise FileNotFoundError(f'Could not find config file {config_file_name}.')

    if mpi.is_master_node():
        print('Reading the config file ' + config_file_name)
        general_params, solver_params, dft_params, advanced_params = read_config(config_file_name)
        general_params['config_file'] = config_file_name

        print('-'*25 + '\nGeneral parameters:')
        for key, value in general_params.items():
            print('{0: <20} {1: <4}'.format(key, str(value)))
        print('-'*25 + '\nSolver parameters:')
        for key, value in solver_params.items():
            print('{0: <20} {1: <4}'.format(key, str(value)))
        print('-'*25 + '\nDFT parameters:')
        for key, value in dft_params.items():
            print('{0: <20} {1: <4}'.format(key, str(value)))
        print('-'*25 + '\nAdvanced parameters, don\'t change them unless you know what you are doing:')
        for key, value in advanced_params.items():
            print('{0: <20} {1: <4}'.format(key, str(value)))

    general_params = mpi.bcast(general_params)
    solver_params = mpi.bcast(solver_params)
    dft_params = mpi.bcast(dft_params)
    advanced_params = mpi.bcast(advanced_params)

    if general_params['csc']:
        # Start CSC calculation, always in same folder as dmft_config
        general_params['jobname'] = '.'
        csc_flow_control(general_params, solver_params, dft_params, advanced_params)
    else:
        # Sets up one-shot calculation
        mpi.report('', '#'*80)
        mpi.report(f'Using input file {general_params["seedname"]}.h5 '
                   + f'and running in folder {general_params["jobname"]}')

        if mpi.is_master_node():
            # Checks for h5 file
            if not os.path.exists(general_params['seedname']+'.h5'):
                raise FileNotFoundError('Input h5 file not found')

            # Creates output directory if it does not exist
            if not os.path.exists(general_params['jobname']):
                os.makedirs(general_params['jobname'])

            # Copies h5 archive and config file to subfolder if are not there
            for file in (general_params['seedname']+'.h5',
                         general_params['config_file']):
                if not os.path.isfile(general_params['jobname']+'/'+file):
                    shutil.copyfile(file, general_params['jobname']+'/'+file)
        mpi.barrier()

        # Runs dmft_cycle
        dmft_cycle(general_params, solver_params, advanced_params,
                   dft_params, general_params['n_iter_dmft'])

    mpi.barrier()
    if mpi.is_master_node():
        global_end = timer()
        print('-------------------------------')
        print('overall elapsed time: %10.4f seconds'%(global_end-global_start))


if __name__ == '__main__':
    main(sys.argv)
