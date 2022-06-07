
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
import socket
from collections import defaultdict
import shlex
import time

import triqs.utility.mpi as mpi


def create_hostfile(number_cores, cluster_name):
    """
    Writes a host file for the mpirun. This tells mpi which nodes to ssh into
    and start VASP on. The format of the hist file depends on the type of MPI
    that is used.

    Parameters
    ----------
    number_cores: int, the number of cores that vasp runs on
    cluster_name: string, the name of the server

    Returns
    -------
    string: name of the hostfile if not run locally and if called by master node
    """

    if cluster_name == 'default':
        return None

    hostnames = mpi.world.gather(socket.gethostname(), root=0)
    if mpi.is_master_node():
        # create hostfile based on first number_cores ranks
        hosts = defaultdict(int)
        for hostname in hostnames[:number_cores]:
            hosts[hostname] += 1

        mask_hostfile = {'openmpi': '{} slots={}',  # OpenMPI format
                         'openmpi-intra': '{} slots={}',  # OpenMPI format
                         'mpich': '{}:{}',  # MPICH format
                         }[cluster_name]

        hostfile = 'dft.hostfile'
        with open(hostfile, 'w') as file:
            file.write('\n'.join(mask_hostfile.format(*i) for i in hosts.items()))
        return hostfile

    return None


def find_path_to_mpi_command(env_vars, mpi_exe):
    """
    Finds the complete path for the mpi executable by scanning the directories
    of $PATH.

    Parameters
    ----------
    env_vars: dict of string, environment variables containing PATH
    mpi_exe: string, mpi command

    Returns
    -------
    string: absolute path to mpi command
    """

    for path_directory in env_vars.get('PATH').split(os.pathsep):
        if path_directory:
            potential_path = os.path.join(path_directory, mpi_exe)
            if os.access(potential_path, os.F_OK | os.X_OK):
                return potential_path

    return None


def get_mpi_arguments(mpi_profile, mpi_exe, number_cores, dft_exe, hostfile):
    """
    Depending on the settings of the cluster and the type of MPI used,
    the arguments to the mpi call have to be different. The most technical part
    of the vasp handler.

    Parameters
    ----------
    cluster_name: string, name of the cluster so that settings can be tailored to it
    mpi_exe: string, mpi command
    number_cores: int, the number of cores that vasp runs on
    dft_exe: string, the command to start the DFT code
    hostfile: string, name of the hostfile

    Returns
    -------
    list of string: arguments to start mpi with
    """

    if mpi_profile == 'default':
        return [mpi_exe, '-np', str(number_cores)] + shlex.split(dft_exe)

    # For the second node, mpirun starts DFT by using ssh
    # Therefore we need to handover the env variables with -x
    if mpi_profile == 'openmpi':
        return [mpi_exe, '-hostfile', hostfile, '-np', str(number_cores),
                '-mca', 'mtl', '^psm,psm2,ofi',
                '-mca', 'btl', '^vader,openib,usnix',
                '-x', 'LD_LIBRARY_PATH',
                '-x', 'PATH', '-x', 'OMP_NUM_THREADS'] + shlex.split(dft_exe)
    # Run mpi with intra-node communication among ranks (on a single node).
    if mpi_profile == 'openmpi-intra':
        return [mpi_exe, '-np', str(number_cores),
                '--mca', 'pml', 'ob1', '--mca', 'btl', 'self,vader',
                '-x', 'LD_LIBRARY_PATH',
                '-x', 'PATH', '-x', 'OMP_NUM_THREADS'] + shlex.split(dft_exe)

    if mpi_profile == 'mpich':
        return [mpi_exe, '-launcher', 'ssh', '-hostfile', hostfile,
                '-np', str(number_cores), '-envlist', 'PATH'] + shlex.split(dft_exe)

    return None


def poll_barrier(comm, poll_interval=0.1):
    """
    Use asynchronous synchronization, otherwise mpi.barrier uses up all the CPU time during
    the run of subprocess.

    Parameters
    ----------
    comm: MPI communicator
    poll_interval: float, time step for pinging the status of the sleeping ranks
    """

    req = comm.Ibarrier()
    while not req.Test():
        time.sleep(poll_interval)
