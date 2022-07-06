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

    # start CSC calculation if csc is set to true
    if general_params['csc']:

        # check if seedname is only one Value
        if len(general_params['seedname']) > 1:
            mpi.report('!!! WARNING !!!')
            mpi.report('CSC calculations can only be done for one set of file at a time')

        # some basic setup that needs to be done for CSC calculations
        general_params['seedname'] = general_params['seedname'][0]
        general_params['jobname'] = '.'
        general_params['previous_file'] = 'none'

        # run the whole machinery
        csc_flow_control(general_params, solver_params, dft_params, advanced_params)

    # do a one-shot calculation with given h5 archive
    else:
        # extract filenames and do a dmft iteration for every h5 archive given
        number_calculations = len(general_params['seedname'])
        filenames = general_params['seedname']
        foldernames = general_params['jobname']
        mpi.report('{} DMFT calculation will be made for the following files: {}'.format(number_calculations, filenames))

        # check for h5 file(s)
        if mpi.is_master_node():
            for file in filenames:
                if not os.path.exists(file+'.h5'):
                    mpi.report('*** Input h5 file(s) not found! I was looking for '+file+'.h5 ***')
                    mpi.MPI.COMM_WORLD.Abort(1)

        for i, file in enumerate(foldernames):
            general_params['seedname'] = filenames[i]
            general_params['jobname'] = foldernames[i]
            if i == 0:
                general_params['previous_file'] = 'none'
            else:
                previous_file = filenames[i-1]
                previous_folder = foldernames[i-1]
                general_params['previous_file'] = previous_folder+'/'+previous_file+'.h5'

            if mpi.is_master_node():
                # create output directory
                print('calculation is performed in subfolder: '+general_params['jobname'])
                if not os.path.exists(general_params['jobname']):
                    os.makedirs(general_params['jobname'])

                    # copy h5 archive and config file to created folder
                    shutil.copyfile(general_params['seedname']+'.h5',
                                    general_params['jobname']+'/'+general_params['seedname']+'.h5')
                    shutil.copyfile(general_params['config_file'],
                                    general_params['jobname']+'/'+general_params['config_file'])
                else:
                    print('#'*80+'\n WARNING! specified job folder already exists continuing previous job! \n'+'#'*80+'\n')

            mpi.report('#'*80)
            mpi.report('starting the DMFT calculation for '+str(general_params['seedname']))
            mpi.report('#'*80)

            ############################################################
            # run the dmft_cycle
            dmft_cycle(general_params, solver_params, advanced_params,
                       dft_params, general_params['n_iter_dmft'])
            ############################################################

    if mpi.is_master_node():
        global_end = timer()
        print('-------------------------------')
        print('overall elapsed time: %10.4f seconds'%(global_end-global_start))


if __name__ == '__main__':
    main(sys.argv)
