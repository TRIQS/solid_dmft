#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
Contains the handling of the QE process. It can start QE, reactivate it,
check if the lock file is there and finally kill QE. Needed for CSC calculations.
"""

import os
import subprocess
import sys

import triqs.utility.mpi as mpi

from solid_dmft.dft_managers import mpi_helpers


def _fork_and_start_qe(mpi_exe, arguments, env_vars, seedname):
    """
    Forks a process from the master process that then calls mpi to start vasp.
    The child process running VASP never leaves this function whereas the main
    process returns the child's process id and continues.

    Parameters
    ----------
    mpi_exe: string, mpi command
    arguments: list of string, arguments to start mpi with
    env_vars: dict of string, environment variables containing PATH
    seedname: string, QE input file

    Returns
    -------
    int: id of the VASP child process
    """

    mpi_arguments, qe_file_ext = arguments[0], arguments[1]

    if qe_file_ext in ['scf', 'nscf', 'pw2wan', 'mod_scf', 'bnd', 'bands', 'proj']:

        inp = open(f'{seedname}.{qe_file_ext}.in', 'r')
        out = open(f'{seedname}.{qe_file_ext}.out', 'w')
        err = open(f'{seedname}.{qe_file_ext}.err', 'w')

        print('  solid_dmft: Starting {} calculation...'.format(qe_file_ext))

        # start subprocess
        qe_process_id = subprocess.run(mpi_arguments, stdin=inp, env=env_vars, capture_output=True,
                                       text=True, shell=False)

        # write output and error file
        output = qe_process_id.stdout
        error = qe_process_id.stderr
        out.writelines(output)
        err.writelines(error)

        if qe_process_id.returncode != 0:
            mpi.report('QE calculation failed. Exiting programm.')
            sys.exit(1)

    elif 'win' in qe_file_ext:

        print('  solid_dmft: Starting Wannier90 {}...'.format(qe_file_ext))

        # don't need any piping for Wannier90
        qe_process_id = subprocess.check_call(mpi_arguments + [seedname], env=env_vars, shell=False)

    return qe_process_id


def start(number_cores, qe_file_ext, qe_exec, mpi_profile, seedname):
    """
    Starts the VASP child process. Takes care of initializing a clean
    environment for the child process. This is needed so that VASP does not
    get confused with all the standard slurm environment variables.

    Parameters
    ----------
    number_cores: int, the number of cores that vasp runs on
    qe_file_ext: string, qe executable
    qe_exec: string, path to qe executables
    mpi_profile: string, name of the cluster so that settings can be tailored to it
    """

    # get MPI env
    qe_process_id = 0

    hostfile = mpi_helpers.create_hostfile(number_cores, mpi_profile)
    qe_exec_path = qe_exec.strip(qe_exec.rsplit('/')[-1])
    qe_exec = qe_exec_path

    if mpi.is_master_node():
        # clean environment
        env_vars = {}
        for var_name in ['PATH', 'LD_LIBRARY_PATH', 'SHELL', 'PWD', 'HOME',
                         'OMP_NUM_THREADS', 'OMPI_MCA_btl_vader_single_copy_mechanism']:
            var = os.getenv(var_name)
            if var:
                env_vars[var_name] = var

        # assuming that mpirun points to the correct mpi env
        mpi_exe = mpi_helpers.find_path_to_mpi_command(env_vars, 'mpirun')

        if qe_file_ext in ['scf', 'nscf', 'mod_scf', 'bnd']:
            qe_exec += f'pw.x -nk {number_cores}'
        elif qe_file_ext in ['pw2wan']:
            qe_exec += 'pw2wannier90.x -nk 1 -pd .true.'
        elif qe_file_ext in ['bands']:
            qe_exec += f'bands.x -nk {number_cores}'
        elif qe_file_ext in ['proj']:
            qe_exec += f'projwfc.x -nk {number_cores}'
        elif qe_file_ext in ['win_pp']:
            qe_exec += 'wannier90.x -pp'
        elif qe_file_ext in ['win']:
            qe_exec += 'wannier90.x'

        arguments = mpi_helpers.get_mpi_arguments(mpi_profile, mpi_exe, number_cores, qe_exec, hostfile)
        qe_process_id = _fork_and_start_qe(mpi_exe, (arguments, qe_file_ext), env_vars, seedname)

    mpi_helpers.poll_barrier(mpi.MPI.COMM_WORLD)
    qe_process_id = mpi.bcast(qe_process_id)

    return qe_process_id
