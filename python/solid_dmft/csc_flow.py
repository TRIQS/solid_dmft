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

from timeit import default_timer as timer
import subprocess
import shlex
import os
import numpy as np

# triqs
from h5 import HDFArchive
import triqs.utility.mpi as mpi

from triqs_dft_tools.converters.wannier90 import Wannier90Converter
from triqs_dft_tools.converters.vasp import VaspConverter
from triqs_dft_tools.converters.plovasp.vaspio import VaspData
import triqs_dft_tools.converters.plovasp.converter as plo_converter

from solid_dmft.dmft_cycle import dmft_cycle
from solid_dmft.dft_managers import vasp_manager as vasp
from solid_dmft.dft_managers import qe_manager as qe

def _run_plo_converter(general_params, dft_params):
    if not mpi.is_master_node():
        return

    # Checks for plo file for projectors
    if not os.path.exists(dft_params['plo_cfg']):
        print('*** Input PLO config file not found! '
              + 'I was looking for {} ***'.format(dft_params['plo_cfg']))
        mpi.MPI.COMM_WORLD.Abort(1)

    # Runs plo converter
    plo_converter.generate_and_output_as_text(dft_params['plo_cfg'], vasp_dir='./')
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

def _full_qe_run(seedname, dft_params, mode):
    assert mode in ('initial', 'restart', 'update')

    # runs a full iteration of DFT
    qe_wrapper = lambda calc_type: qe.run(dft_params['n_cores'], calc_type, dft_params['dft_exec'],
                                          dft_params['mpi_env'], seedname)

    # Initially run an scf calculation
    if mode == 'initial':
        qe_wrapper('scf')
    # For charge update, use mode scf
    elif mode == 'update':
        qe_wrapper('mod_scf')

    # Rest is executed regardless of mode
    # Optionally does bnd, bands, proj if files are present
    for nscf in ['bnd', 'bands', 'proj']:
        if os.path.isfile(f'{seedname}.{nscf}.in'):
            qe_wrapper(nscf)

    # nscf
    qe_wrapper('nscf')
    # w90 parts
    qe_wrapper('win_pp')
    qe_wrapper('pw2wan')
    qe_wrapper('win')
    _run_w90converter(seedname, dft_params['w90_tolerance'])


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

def _full_vasp_run(general_params, dft_params, initial_run, n_iter_dft=1, sum_k=None):
    """
    Performs a complete DFT cycle in Vasp and the correct converter. If
    initial_run, Vasp is starting and performing a normal scf calculation
    followed by a converter run. Otherwise, it performs n_iter_dft runs of DFT,
    generating the projectors with the converter, and recalculating the charge
    density correction with the new projectors.

    Parameters
    ----------
    general_params : dict
        general parameters as a dict
    dft_params : dict
        dft parameters as a dict
    initial_run : bool
        True when VASP is called for the first time. initial_run = True requires
        n_iter_dft = 1.
    n_iter_dft : int, optional
        Number of DFT iterations to perform. The default is 1.
    sum_k : SumkDFT, optional
        The SumkDFT object required to recalculate the charge-density correction
        if n_iter_dft > 1. The default is None.

    Returns
    -------
    vasp_process_id : int
        The process ID of the forked VASP process.
    irred_indices : np.array
        Integer indices of kpts in the irreducible Brillouin zone. Only needed
        for Wannier projectors, which are normally run with symmetries.
    """


    if initial_run:
        assert n_iter_dft == 1
    else:
        assert n_iter_dft == 1 or sum_k is not None, 'Sumk object needed to run multiple DFT iterations'

    for i in range(n_iter_dft):
        if initial_run:
            vasp_process_id = vasp.run_initial_scf(dft_params['n_cores'], dft_params['dft_exec'],
                                                   dft_params['mpi_env'])
        else:
            vasp_process_id = None
            vasp.run_charge_update()

        if dft_params['projector_type'] == 'plo':
            _run_plo_converter(general_params, dft_params)
            irred_indices = None
        elif dft_params['projector_type'] == 'w90':
            _run_wannier90(general_params, dft_params)
            mpi.barrier()
            _run_w90converter(general_params['seedname'], dft_params['w90_tolerance'])
            mpi.barrier()
            kpts = None
            if mpi.is_master_node():
                with HDFArchive(general_params['seedname']+'.h5', 'r') as archive:
                    kpts = archive['dft_input/kpts']
            irred_indices = vasp.read_irred_kpoints(kpts)

        # No need for recalculation of density correction if we run DMFT next
        if i == n_iter_dft - 1:
            break

        # Recalculates the density correction
        # Reads in new projectors and hopping and updates chemical potential
        # rot_mat is not updated since it's more closely related to the local problem than DFT
        # New fermi weights are directly read in calc_density_correction
        mpi.barrier()
        if mpi.is_master_node():
            with HDFArchive(general_params['seedname']+'.h5', 'r') as archive:
                sum_k.proj_mat = archive['dft_input/proj_mat']
                sum_k.hopping = archive['dft_input/hopping']
        sum_k.proj_mat = mpi.bcast(sum_k.proj_mat)
        sum_k.hopping = mpi.bcast(sum_k.hopping)
        sum_k.calc_mu(precision=general_params['prec_mu'])

        # Writes out GAMMA file
        sum_k.calc_density_correction(dm_type='vasp',  kpts_to_write=irred_indices)

    return vasp_process_id, irred_indices


# Main CSC flow method
def csc_flow_control(general_params, solver_params, dft_params, advanced_params):
    """
    Function to run the csc cycle. It writes and removes the vasp.lock file to
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
    """

    # Removes legacy file vasp.suppress_projs if present
    vasp.remove_legacy_projections_suppressed()

    # if GAMMA file already exists, load it by doing extra DFT iterations
    if dft_params['dft_code'] == 'vasp' and os.path.exists('GAMMA'):
        # TODO: implement
        raise NotImplementedError('GAMMA file found but restarting from updated '
                                  + 'charge density not yet implemented for Vasp.')

    # Reads in iteration offset if restarting
    iteration_offset = 0
    if mpi.is_master_node() and os.path.isfile(general_params['seedname']+'.h5'):
        with HDFArchive(general_params['seedname']+'.h5', 'r') as archive:
            if 'DMFT_results' in archive and 'iteration_count' in archive['DMFT_results']:
                iteration_offset = archive['DMFT_results']['iteration_count']
    iteration_offset = mpi.bcast(iteration_offset)

    iter_dmft = iteration_offset+1

    # Runs DFT once and converter
    mpi.barrier()
    irred_indices = None
    start_time_dft = timer()
    mpi.report('  solid_dmft: Running {}...'.format(dft_params['dft_code'].upper()))

    if dft_params['dft_code'] == 'qe':
        if iteration_offset == 0:
            _full_qe_run(general_params['seedname'], dft_params, 'initial')
        else:
            _full_qe_run(general_params['seedname'], dft_params, 'restart')
    elif dft_params['dft_code'] == 'vasp':
        vasp_process_id, irred_indices = _full_vasp_run(general_params, dft_params, True)

    mpi.barrier()
    end_time_dft = timer()
    mpi.report('  solid_dmft: DFT cycle took {:10.3f} seconds'.format(end_time_dft-start_time_dft))

    # Now that everything is ready, starts DFT+DMFT loop
    while True:
        dft_energy = None
        if mpi.is_master_node():
            # Writes eigenvals to archive if requested
            if dft_params['store_eigenvals']:
                if dft_params['dft_code'] == 'qe':
                    # TODO: implement
                    raise NotImplementedError('store_eigenvals not yet compatible with dft_code = qe')
                _store_dft_eigvals(path_to_h5=general_params['seedname']+'.h5',
                                   iteration=iter_dmft,
                                   projector_type=dft_params['projector_type'])

            # Reads the DFT energy
            if dft_params['dft_code'] == 'vasp':
                dft_energy = vasp.read_dft_energy()
            elif dft_params['dft_code'] == 'qe':
                dft_energy = qe.read_dft_energy(general_params['seedname'], iter_dmft)
        dft_energy = mpi.bcast(dft_energy)

        mpi.report('', '#'*80, 'Calling dmft_cycle')

        if mpi.is_master_node():
            start_time_dmft = timer()

        # Determines number of DMFT steps
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
        is_converged, sum_k = dmft_cycle(general_params, solver_params, advanced_params,
                                         dft_params, iter_one_shot, irred_indices, dft_energy)
        ############################################################

        iter_dmft += iter_one_shot

        if mpi.is_master_node():
            end_time_dmft = timer()
            print('\n' + '='*80)
            print('DMFT cycle took {:10.3f} seconds'.format(end_time_dmft-start_time_dmft))
            print('='*80 + '\n')

        # If all steps are executed or calculation is converged, finish DFT+DMFT loop
        if is_converged or iter_dmft > general_params['n_iter_dmft'] + iteration_offset:
            break

        # Restarts DFT
        mpi.barrier()
        start_time_dft = timer()
        mpi.report('  solid_dmft: Running {}...'.format(dft_params['dft_code'].upper()))

        # Runs DFT and converter
        if dft_params['dft_code'] == 'qe':
            _full_qe_run(general_params['seedname'], dft_params, 'update')
        elif dft_params['dft_code'] == 'vasp':
            # Determines number of DFT steps
            if iter_dmft == general_params['n_iter_dmft_first'] + 1:
                n_iter_dft = dft_params['n_iter_first']
            else:
                n_iter_dft = dft_params['n_iter']
            _, irred_indices = _full_vasp_run(general_params, dft_params, False, n_iter_dft, sum_k)

        mpi.barrier()
        end_time_dft = timer()
        mpi.report('  solid_dmft: DFT cycle took {:10.3f} seconds'.format(end_time_dft-start_time_dft))

    # Kills background VASP process for clean end
    if mpi.is_master_node() and dft_params['dft_code'] == 'vasp':
        print('  solid_dmft: Stopping VASP\n', flush=True)
        vasp.kill(vasp_process_id)
