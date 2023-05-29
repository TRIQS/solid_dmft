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
Contains the function to run a QuantumEspresso iteration. Needed for CSC calculations.
"""

import os
import socket
from collections import defaultdict
import subprocess
import shlex
import sys
import time

import triqs.utility.mpi as mpi


def _create_hostfile(number_cores, cluster_name):
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

    if cluster_name == 'local':
        return None

    hostnames = mpi.world.gather(socket.gethostname(), root=0)
    if mpi.is_master_node():
        # create hostfile based on first number_cores ranks
        hosts = defaultdict(int)
        for hostname in hostnames[:number_cores]:
            hosts[hostname] += 1

        mask_hostfile = {'rusty': '{} slots={}', # OpenMPI format
                         'rusty-intra': '{} slots={}', #OpenMPI format
                         'rusty-ompi2': '{} slots={}', # OpenMPI format
                         'daint': '{}:{}', # MPICH format
                        }[cluster_name]

        hostfile = 'dft.hostfile'
        with open(hostfile, 'w') as file:
            file.write('\n'.join(mask_hostfile.format(*i) for i in hosts.items()))
        return hostfile

    return None

def _find_path_to_mpi_command(env_vars, mpi_exe):
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

def _get_mpi_arguments(cluster_name, mpi_exe, number_cores, qe_exec, hostfile):
    """
    Depending on the settings of the cluster and the type of MPI used,
    the arguments to the mpi call have to be different. The most technical part
    of the vasp handler.

    Parameters
    ----------
    cluster_name: string, name of the cluster so that settings can be tailored to it
    mpi_exe: string, mpi command
    number_cores: int, the number of cores that vasp runs on
    vasp_command: string, the command to start vasp
    hostfile: string, name of the hostfile

    Returns
    -------
    list of string: arguments to start mpi with
    """

    if cluster_name == 'local':
        return [mpi_exe, '-np', str(number_cores)] + shlex.split(qe_exec)

    # For the second node, mpirun starts VASP by using ssh
    # Therefore we need to handover the env variables with -x
    if cluster_name == 'rusty':
        return [mpi_exe, '-hostfile', hostfile, '-np', str(number_cores),
                '-mca', 'mtl', '^psm,psm2,ofi',
                '-mca', 'btl', '^vader,openib,usnix' ,
                '-x', 'LD_LIBRARY_PATH',
                '-x', 'PATH', '-x', 'OMP_NUM_THREADS'] + shlex.split(qe_exec)
    # Run mpi with intra-node communication among ranks (on a single node).
    if cluster_name == 'rusty-intra':
        return [mpi_exe, '-np', str(number_cores),
                '--mca', 'pml', 'ob1', '--mca', 'btl', 'self,vader',
                '-x', 'LD_LIBRARY_PATH',
                '-x', 'PATH', '-x', 'OMP_NUM_THREADS'] + shlex.split(qe_exec)
    if cluster_name == 'rusty-ompi2':
        return [mpi_exe, '-hostfile', hostfile, '-np', str(number_cores),
                '-mca', 'mtl', '^psm2,ofi', '-mca', 'btl', '^vader' ,'-x', 'LD_LIBRARY_PATH',
                '-x', 'PATH', '-x', 'OMP_NUM_THREADS'] + shlex.split(qe_exec)

    if cluster_name == 'daint':
        return [mpi_exe, '-launcher', 'ssh', '-hostfile', hostfile,
                '-np', str(number_cores), '-envlist', 'PATH'] + shlex.split(qe_exec)

    return None

def _poll_barrier(comm, poll_interval = 0.1):
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

def _start_with_piping(mpi_exe, mpi_arguments, qe_file_ext, env_vars, seedname):
    """
    Handles the piping of the output when starting QE.

    Parameters
    ----------
    mpi_exe: string, mpi command
    mpi_arguments: list of string, arguments to start mpi with
    qe_file_ext : string, file name for QE
    env_vars: dict of string, environment variables containing PATH
    seedname: string, QE input file

    Returns
    -------
    int: id of the VASP child process
    """

    if qe_file_ext in ['scf', 'nscf', 'pw2wan', 'mod_scf', 'bnd', 'bands', 'proj']:

        inp = open(f'{seedname}.{qe_file_ext}.in', 'r')
        out = open(f'{seedname}.{qe_file_ext}.out','w')
        err = open(f'{seedname}.{qe_file_ext}.err','w')

        print('  solid_dmft: Starting {} calculation...'.format(qe_file_ext))

        # start subprocess
        qe_result = subprocess.run(mpi_arguments, stdin=inp, env=env_vars, capture_output=True,
                                   text=True, shell=False)

        # write output and error file
        output = qe_result.stdout
        error = qe_result.stderr
        out.writelines(output)
        err.writelines(error)

        if qe_result.returncode != 0:
            mpi.report('QE calculation failed. Exiting programm.')
            sys.exit(1)

    elif 'win' in qe_file_ext:
        print('  solid_dmft: Starting Wannier90 {}...'.format(qe_file_ext))
        # don't need any piping for Wannier90
        subprocess.check_call(mpi_arguments + [seedname], env=env_vars, shell=False)


def run(number_cores, qe_file_ext, qe_exec, mpi_profile, seedname):
    """
    Starts the VASP child process. Takes care of initializing a clean
    environment for the child process. This is needed so that VASP does not
    get confused with all the standard slurm environment variables.

    Parameters
    ----------
    number_cores: int, the number of cores that vasp runs on
    vasp_command: string, the command to start vasp
    cluster_name: string, name of the cluster so that settings can be tailored to it
    """

    # get MPI env
    hostfile = _create_hostfile(number_cores, mpi_profile)
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
        mpi_exe = _find_path_to_mpi_command(env_vars, 'mpirun')

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

        arguments = _get_mpi_arguments(mpi_profile, mpi_exe, number_cores, qe_exec, hostfile)
        _start_with_piping(mpi_exe, arguments, qe_file_ext, env_vars, seedname)

    _poll_barrier(mpi.MPI.COMM_WORLD)

def read_dft_energy(seedname, iter_dmft):
    """
    Reads DFT energy from quantum espresso's out files

    1. At the first iteration, the DFT energy is read from the scf file.

    2. After the first iteration the band energy computed in the mod_scf calculation is wrong,
       and needs to be subtracted from the reported total energy. The correct band energy
       is computed in the nscf calculation.

    """
    dft_energy = 0.0
    RYDBERG = 13.605693123 # eV

    if iter_dmft == 1:
        with open(f'{seedname}.scf.out', 'r') as file:
            dft_output = file.readlines()
        for line in dft_output:
            if '!' in line:
                print("\nReading total energy from the scf calculation \n")
                dft_energy = float(line.split()[-2]) * RYDBERG
                print(f"The DFT energy is: {dft_energy} eV")
                break
            if  line =="":
                raise EOFError("Did not find scf total energy")
    else:
        with open(f'{seedname}.mod_scf.out', 'r') as file:
            dft_output = file.readlines()
        for line in dft_output:
            #if 'eband, Ef (eV)' in line:
            if "(sum(wg*et))" in line:
                print("\nReading band energy from the mod_scf calculation \n")
                #band_energy = float(line.split())
                band_energy_modscf = float(line.split()[-2])*RYDBERG
                print(f"The mod_scf band energy is: {band_energy_modscf} eV")
            if 'total energy' in line:
                print("\nReading total energy from the mod_scf calculation \n")
                dft_energy = float(line.split()[-2]) * RYDBERG
                print(f"The uncorrected DFT energy is: {dft_energy} eV")
        dft_energy -= band_energy_modscf
        print(f"The DFT energy without kinetic part is: {dft_energy} eV")

        with open(f'{seedname}.nscf.out', 'r') as file:
            dft_output = file.readlines()
        for line in dft_output:
            if 'The nscf band energy' in line:
                print("\nReading band energy from the nscf calculation\n")
                band_energy_nscf = float(line.split()[-2]) * RYDBERG
                dft_energy += band_energy_nscf
                print(f"The nscf band energy is: {band_energy_nscf} eV")
                print(f"The corrected DFT energy is: {dft_energy} eV")
                break
    return dft_energy
