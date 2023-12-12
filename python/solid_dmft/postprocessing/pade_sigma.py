# pyright: reportUnusedExpression=false

import numpy as np

from triqs.utility import mpi
from h5 import HDFArchive
from triqs.gf import Gf, MeshReFreq, BlockGf

from solid_dmft.postprocessing.maxent_sigma import _read_h5

def _write_sigma_omega_to_h5(sigma_w, external_path, iteration):
    """ Writes real-frequency self energy to h5 archive. """
    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None
                                          else f'it_{iteration}')

    with HDFArchive(external_path, 'a') as archive:
        for i, sigma_imp in enumerate(sigma_w):
            archive[h5_internal_path][f'Sigma_Refreq_{i}'] = sigma_imp

def _run_pade(sigma_iw_list, n_w, w_min, w_max, n_iw, eta):
    """
    Run pade in parallel. Call via main function.
    """
    mpi.report('Continuing impurities with blocks:')

    imps_blocks = []
    sigma_iw_flat_list = []

    # create flattened list of self-energies
    for i, sigma_iw in enumerate(sigma_iw_list):
        blocks = list(sigma_iw.indices)
        mpi.report('- Imp {}: {}'.format(i, blocks))
        for block in blocks:
            imps_blocks.append((i, block))
            sigma_iw_flat_list.append(sigma_iw[block])

    sigma_w_flat_list = []
    wmesh = MeshReFreq(w_min=w_min,w_max=w_max,n_w=n_w)
    imps_blocks_indices = np.arange(len(imps_blocks))
    for i in imps_blocks_indices:
        sigma_w_flat_list.append(Gf(mesh=wmesh, target_shape=sigma_iw_flat_list[i].target_shape))

    # Runs Pade while parallelizing over impurities and blocks
    for i in mpi.slice_array(imps_blocks_indices):
        print(f'Rank {mpi.rank} continuing Σ {i}/{len(imps_blocks)}')
        sigma_w_flat_list[i].set_from_pade(sigma_iw_flat_list[i],n_points=n_iw, freq_offset=eta)

    # sync Pade data
    for i in imps_blocks_indices:
        sigma_w_flat_list[i] = mpi.all_reduce(sigma_w_flat_list[i])

    # Create list of BlockGf
    sigma_w_list = []
    for i, sigma_iw in enumerate(sigma_iw_list):
        block_list = []
        for block in sigma_iw.indices:
            block_list.append(sigma_w_flat_list.pop(0))
        sigma_w_list.append(BlockGf(name_list=list(sigma_iw.indices), block_list=block_list, make_copies=True))

    return sigma_w_list

def main(external_path, n_w, w_min, w_max, n_iw, iteration=None, eta=0.0):
    """
    Main function that reads the Matsubara self-energy from h5, analytically continues it,
    writes the results back to the h5 archive and also returns the results.

    Function parallelizes using MPI over impurities and blocks.

    Parameters
    ----------
    external_path : string
        Path to the h5 archive to read from and write to
    n_w : int
        number of real frequencies of the final self-energies returned
    w_min : float
        Lower end of range where Sigma is being continued.
    w_max : float
        Upper end of range where Sigma is being continued.
    n_iw : int
        number of Matsubara frequencies to consider for the Pade approximant
    iteration : int/string
        Iteration to read from and write to. Default to last_iter
    eta : float
        frequency offset within Pade

    Returns
    -------
    sigma_w : list of triqs.gf.BlockGf
        Sigma(omega) per inequivalent shell
    """

    sigma_iw = None
    if mpi.is_master_node():
        sigma_iw, _, _, _ = _read_h5(external_path, iteration)
    sigma_iw = mpi.bcast(sigma_iw)

    # run pade in parallel
    sigma_w = _run_pade(sigma_iw, n_w, w_min, w_max, n_iw, eta)

    mpi.report('Writing results to h5 archive now.')
    if mpi.is_master_node():
        _write_sigma_omega_to_h5(sigma_w, external_path, iteration)
    mpi.report('Finished writing Σ(ω) to archive.')

    return sigma_w

