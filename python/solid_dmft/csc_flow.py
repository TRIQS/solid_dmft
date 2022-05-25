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
contains the charge self-consistency flow control functions
"""

import time
from timeit import default_timer as timer
import subprocess
import shlex
import os
import numpy as np

# triqs
from h5 import HDFArchive
import triqs.utility.mpi as mpi
from triqs_dft_tools.converters.wannier90 import Wannier90Converter
try:
    from triqs_dft_tools.converters.vasp import VaspConverter
    from triqs_dft_tools.converters.plovasp.vaspio import VaspData
    import triqs_dft_tools.converters.plovasp.converter as plo_converter
except ImportError:
    pass

from solid_dmft.dmft_cycle import dmft_cycle
from solid_dmft.dft_managers import vasp_manager as vasp
from solid_dmft.dft_managers import qe_manager as qe

def _run_plo_converter(general_params):
    if not mpi.is_master_node():
        return

    # Checks for plo file for projectors
    if not os.path.exists(general_params['plo_cfg']):
        print('*** Input PLO config file not found! '
              + 'I was looking for {} ***'.format(general_params['plo_cfg']))
        mpi.MPI.COMM_WORLD.Abort(1)

    # Runs plo converter
    plo_converter.generate_and_output_as_text(general_params['plo_cfg'], vasp_dir='./')
    # Writes new H(k) to h5 archive
    converter = VaspConverter(filename=general_params['seedname'])
    converter.convert_dft_input()

def _run_wannier90(general_params, dft_params):
    if not mpi.is_master_node():
        return

    if not os.path.exists(general_params['seedname'] + '.win'):
        print('*** Wannier input file not found! '
              + 'I was looking for {0}.win ***'.format(general_params['seedname']))
        mpi.MPI.COMM_WORLD.Abort(1)

    # Runs wannier90 twice:
    # First preprocessing to write nnkp file, then normal run
    command = shlex.split(dft_params['w90_exec'])
    subprocess.check_call(command + ['-pp', general_params['seedname']], shell=False)
    subprocess.check_call(command + [general_params['seedname']], shell=False)

def _run_w90converter(seedname, tolerance):
    if (not os.path.exists(seedname + '.win')
        or not os.path.exists(seedname + '.inp')):
        print('*** Wannier input/converter config file not found! '
              + 'I was looking for {0}.win and {0}.inp ***'.format(seedname))
        mpi.MPI.COMM_WORLD.Abort(1)

    #TODO: choose rot_mat_type with general_params['set_rot']
    converter = Wannier90Converter(seedname, rot_mat_type='hloc_diag', bloch_basis=True, w90zero=tolerance)
    converter.convert_dft_input()
    mpi.barrier()

    # Checks if creating of rot_mat succeeded
    if mpi.is_master_node():
        with HDFArchive(seedname+'.h5', 'r') as archive:
            assert archive['dft_input']['use_rotations'], 'Creation of rot_mat failed in W90 converter'
    mpi.barrier()

def _run_qe(general_params, dft_params, iter_dmft, iteration_offset):
    # runs a full iteration of DFT
    start_qe = lambda n_cores, calc_type: qe.start(n_cores, calc_type, dft_params['dft_exec'],
                                                   dft_params['mpi_env'], general_params['seedname'])

    if iter_dmft == 1: # scf
        qe_scf = start_qe(dft_params['n_cores'], 'scf')
    else: # use modified scf

        # if calculation is restarted, need to check in first iteration if DFT step needs to be skipped
        iter_one_shot = (iter_dmft - 1 - general_params['n_iter_dmft_first'])%general_params['n_iter_dmft_per']
        if iteration_offset > 0 and iter_dmft == iteration_offset + 1 and iter_one_shot > 0:
            mpi.report('  solid_dmft: ...skipping DFT step')
            return

        qe_scf = start_qe(dft_params['n_cores'], 'mod_scf')
    # optionally do bnd, bands, proj if files are present
    for nscf in ['bnd', 'bands', 'proj']:
        if os.path.isfile(f'{general_params["seedname"]}.{nscf}.in'):
            qe_nscf = start_qe(dft_params['n_cores'], nscf)
    # nscf
    qe_nscf = start_qe(dft_params['n_cores'], 'nscf')
    # w90 parts
    qe_w90 = start_qe(dft_params['n_cores'], 'win_pp')
    qe_pw2wan = start_qe(dft_params['n_cores'], 'pw2wan')
    qe_w90 = start_qe(dft_params['n_cores'], 'win')

    # launch Wannier90Converter
    _run_w90converter(general_params['seedname'], dft_params['w90_tolerance'])

def read_dft_energy_vasp():
    """
    Reads DFT energy from the last line of Vasp's OSZICAR.
    """
    try:
        with open('OSZICAR', 'r') as file:
            nextline = file.readline()
            while nextline.strip():
                line = nextline
                nextline = file.readline()
        dft_energy = float(line.split()[2])
    except FileNotFoundError:
        print('OSZICAR not found, cannot read DFT energy')
        return None
    except ValueError:
        print('Failed to read DFT energy from OSZICAR')
        return None

    print('DFT energy read from OSZICAR')
    return dft_energy

def read_dft_energy_qe(seedname):
    """
    Reads DFT energy from quantum espresso's seedname.scf.out.
    """
    RYDBERG = 13.605698066 # eV

    for mode in ['mod_scf', 'scf']:
        try:
            with open(f'{seedname}.{mode}.out', 'r') as file:
                dft_output = file.readlines()
            for line in dft_output:
                if 'total energy' in line:
                    dft_energy = float(line.split()[-2]) * RYDBERG
            break
        except FileNotFoundError:
            if mode == 'scf':
                print(f'{seedname}.{mode}.out not found, cannot read DFT energy.')
        except ValueError:
            pass

    print(f'DFT energy read from {seedname}.{mode}.out')
    return dft_energy


def _set_projections_suppressed(suppressed):
    if mpi.is_master_node():
        if suppressed:
            print('  solid_dmft: Writing suppress projectors file', flush=True)
            open('./vasp.suppress_projs', 'w').close()
        elif os.path.isfile('./vasp.suppress_projs'):
            print('  solid_dmft: Removing suppress projectors file', flush=True)
            os.remove('./vasp.suppress_projs')
    mpi.barrier()

def _read_irred_kpoints_from_vasp(general_params):
    irred_indices = None
    if mpi.is_master_node():
        with HDFArchive(general_params['seedname'] + '.h5', 'r') as archive:
            kpts = archive['dft_input/kpts']

        def read_outcar(file):
            has_started_reading = False
            for line in file:
                if 'IBZKPT_HF' in line:
                    has_started_reading = True
                    continue

                if not has_started_reading:
                    continue

                if 't-inv' in line:
                    yield line
                    continue

                if '-'*10 in line:
                    break

        with open('OUTCAR', 'r') as file:
            outcar_data_raw = np.loadtxt(read_outcar(file), usecols=[0, 1, 2, 4])
        outcar_kpoints = outcar_data_raw[:, :3]
        outcar_indices = (outcar_data_raw[:, 3]-.5).astype(int)

        symmetry_mapping = np.full(outcar_kpoints.shape[0], -1, dtype=int)

        for i, (kpt_outcar, outcar_index) in enumerate(zip(outcar_kpoints, outcar_indices)):
            for j, kpt in enumerate(kpts):
                if np.allclose(kpt_outcar, kpt):
                    # Symmetry-irreducible k points
                    if i == outcar_index:
                        symmetry_mapping[j] = outcar_index
                    # Symmetry-reducible
                    else:
                        symmetry_mapping[j] = outcar_index
                    break

            # Asserts that loop left through break, i.e. a pair was found
            assert np.allclose(kpt_outcar, kpt)

        irreds, irred_indices = np.unique(symmetry_mapping, return_index=True)
        assert np.all(np.diff(irreds) == 1)
        assert np.all(symmetry_mapping >= 0)

    return mpi.bcast(irred_indices)

def _store_dft_eigvals(path_to_h5, iteration, projector_type):
    """
    save the eigenvalues from LOCPROJ/wannier90 file to h5 archive
    """
    with HDFArchive(path_to_h5, 'a') as archive:
        if 'dft_eigvals' not in archive:
            archive.create_group('dft_eigvals')

        if projector_type == 'plo':
            vasp_data = VaspData('./')
            eigenvals = np.array(vasp_data.plocar.eigs[:, :, 0]) - vasp_data.plocar.efermi
        elif projector_type == 'w90':
            with open('LOCPROJ') as locproj_file:
                fermi_energy = float(locproj_file.readline().split()[4])
            n_k = archive['dft_input']['n_k']
            num_ks_bands = archive['dft_input']['n_orbitals'][0, 0]
            eigenvals = np.loadtxt('wannier90.eig', usecols=2)
            eigenvals = eigenvals.reshape((n_k, num_ks_bands)) - fermi_energy

        archive['dft_eigvals']['it_'+str(iteration)] = eigenvals

# Main CSC flow method
def csc_flow_control(general_params, solver_params, dft_params, advanced_params):
    """
    function to run the csc cycle. It writes and removes the vasp.lock file to
    start and stop Vasp, run the converter, run the dmft cycle and abort the job
    if all iterations are finished.

    Parameters
    ----------
    general_params : dict
        general parameters as a dict
    solver_params : dict
        solver parameters as a dict
    dft_params : dict
        dft parameters as a dict
    advanced_params : dict
        advanced parameters as a dict

    Returns
    -------
    nothing

    """

    iter_dft = 0
    irred_indices = None
    if dft_params['dft_code'] == 'vasp':
        # if GAMMA file already exists, load it by doing extra DFT iterations
        if os.path.exists('GAMMA'):
            iter_dft = -dft_params['n_iter']
            mpi.barrier()

        vasp_process_id = vasp.start(dft_params['n_cores'], dft_params['dft_exec'],
                                     dft_params['mpi_env'])
        # Removes projection suppression file if present
        _set_projections_suppressed(False)
        mpi.report('  solid_dmft: Waiting for VASP to start (lock appears)...')
        while not vasp.is_lock_file_present():
            time.sleep(1)

    # Reads in iteration offset if restarting
    iteration_offset = 0
    if mpi.is_master_node() and os.path.isfile(general_params['seedname']+'.h5'):
        with HDFArchive(general_params['seedname']+'.h5', 'r') as archive:
            if 'DMFT_results' in archive and 'iteration_count' in archive['DMFT_results']:
                iteration_offset = archive['DMFT_results']['iteration_count']
    iteration_offset = mpi.bcast(iteration_offset)

    iter_dmft = iteration_offset+1
    start_time_dft = timer()
    while iter_dmft <= general_params['n_iter_dmft'] + iteration_offset:
        mpi.report('  solid_dmft: Running {}...'.format(dft_params['dft_code'].upper()))
        mpi.barrier()

        # vasp dft run
        if dft_params['dft_code'] == 'vasp':
            while vasp.is_lock_file_present():
                time.sleep(1)
        # qe dft run
        elif dft_params['dft_code'] == 'qe':
            _run_qe(general_params, dft_params, iter_dmft, iteration_offset)

        # check if we should do another DFT iteration or go on with DMFT
        iter_dft += 1
        if dft_params['dft_code'] == 'vasp' and ((iter_dft-1) % dft_params['n_iter'] != 0 or iter_dft < 0):
            suppressed = iter_dft%dft_params['n_iter'] != 0 or iter_dft+1 < 0
            _set_projections_suppressed(suppressed)
            vasp.reactivate()
            continue

        end_time_dft = timer()
        mpi.report('  solid_dmft: DFT cycle took {:10.4f} seconds'.format(end_time_dft-start_time_dft))

        if dft_params['dft_code'] == 'vasp':
            # Runs the converter
            if dft_params['projector_type'] == 'plo':
                _run_plo_converter(general_params)
            elif dft_params['projector_type'] == 'w90':
                _run_wannier90(general_params, dft_params)
                mpi.barrier()
                _run_w90converter(general_params['seedname'], dft_params['w90_tolerance'])
                mpi.barrier()
                irred_indices = _read_irred_kpoints_from_vasp(general_params)
                mpi.barrier()

        # Writes eigenvals to archive if requested
        dft_energy = None
        if mpi.is_master_node():
            # TODO: not yet compatible with dft_code = qe
            if dft_params['store_eigenvals']:
                _store_dft_eigvals(path_to_h5=general_params['seedname']+'.h5',
                                   iteration=iter_dmft,
                                   projector_type=dft_params['projector_type'])

            if dft_params['dft_code'] == 'vasp':
                dft_energy = read_dft_energy_vasp()
            elif dft_params['dft_code'] == 'qe':
                dft_energy = read_dft_energy_qe(general_params['seedname'])
        dft_energy = mpi.bcast(dft_energy)

        if mpi.is_master_node():
            print('\n' + '='*80)
            print('Calling dmft_cycle')
            print('DMFT iteration {} / {}'.format(iter_dmft,
                                                  general_params['n_iter_dmft']+iteration_offset))
            print('='*80 + '\n')

        if mpi.is_master_node():
            start_time_dmft = timer()

        # Determines number of steps
        if iter_dmft == 1:
            iter_one_shot = general_params['n_iter_dmft_first']
        elif iteration_offset > 0 and iter_dmft == iteration_offset + 1:
            iter_one_shot = general_params['n_iter_dmft_per'] - (iter_dmft - 1
                            - general_params['n_iter_dmft_first'])%general_params['n_iter_dmft_per']
        else:
            iter_one_shot = general_params['n_iter_dmft_per']
        # Maximum total number of iterations is n_iter_dmft+iteration_offset
        iter_one_shot = min(iter_one_shot,
                            general_params['n_iter_dmft'] + iteration_offset - iter_dmft + 1)

        ############################################################
        # run the dmft_cycle
        is_converged = dmft_cycle(general_params, solver_params, advanced_params,
                                  dft_params, iter_one_shot, irred_indices, dft_energy)
        ############################################################

        iter_dmft += iter_one_shot

        if mpi.is_master_node():
            end_time_dmft = timer()
            print('\n' + '='*80)
            print('DMFT cycle took {:10.4f} seconds'.format(end_time_dmft-start_time_dmft))
            print('='*80 + '\n')

        if is_converged:
            break

        # Reactivates Vasp if we do another run
        if iter_dmft <= general_params['n_iter_dmft'] + iteration_offset:
            mpi.barrier()

            start_time_dft = timer()
            if dft_params['dft_code'] == 'vasp':
                _set_projections_suppressed(dft_params['n_iter'] > 1)
                mpi.report('  solid_dmft: Reactivating VASP')
                vasp.reactivate()

    # Stops after maximum number of dmft iterations or convergence
    if mpi.is_master_node():
        print('  solid_dmft: Stopping {}\n'.format(dft_params['dft_code'].upper()), flush=True)
        if dft_params['dft_code'] == 'vasp':
            vasp.kill(vasp_process_id)
