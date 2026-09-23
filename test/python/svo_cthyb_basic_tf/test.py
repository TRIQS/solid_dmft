import sys
import shutil

import numpy as np
from triqs.gfs import *
from triqs.utility.comparison_tests import assert_block_gfs_are_close, assert_arrays_are_close
from h5 import HDFArchive
import triqs.utility.mpi as mpi

import solid_dmft.main as solid


if mpi.is_master_node():
    shutil.rmtree('out', ignore_errors=True)

mpi.barrier()

solid.main([None, 'dmft_config.toml'])

# dc_orb_shift breaks the t2g degeneracy. If the shifted orbitals were still
# symmetrized, the impurity levels passed via the delta interface would not
# match G0 and Delta(tau) would pick up spikes of the wrong sign at tau=0, beta
# of order 10. The tolerance leaves room for the QMC noise of the last iteration
if mpi.is_master_node():
    with HDFArchive('out/inp.h5', 'r') as ar:
        for block, delta in ar['DMFT_results/last_iter/Delta_time_0']:
            diag = np.einsum('tii->ti', delta.data).real
            assert np.all(diag < 0.05), f'Delta(tau) of block {block} is not negative, max = {diag.max()}'

# with HDFArchive('out/inp.h5','r') as ar:
#     G_iw = ar['DMFT_results/last_iter/Gimp_freq_0']
#     G_iw_direct = ar['DMFT_results/last_iter/Gimp_freq_direct_0']
#
#
# # compare direct measured G_iw
# iw0 = len(G_iw.mesh)//2
# n_iw = len(G_iw_direct.mesh) // 2
#
# # works only if statistics are improved to higher precision
# for block, gf in G_iw_direct:
#     assert_arrays_are_close(G_iw[block].data[iw0-n_iw:iw0+n_iw, :, :], G_iw_direct[block].data, precision=5e-2)
