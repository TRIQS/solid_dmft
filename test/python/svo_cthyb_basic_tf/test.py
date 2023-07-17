import sys
import shutil

from triqs.gf import *
from triqs.utility.comparison_tests import assert_block_gfs_are_close, assert_arrays_are_close
from h5 import HDFArchive
import triqs.utility.mpi as mpi

import solid_dmft.main as solid


if mpi.is_master_node():
    shutil.rmtree('out', ignore_errors=True)

mpi.barrier()

solid.main([None, 'dmft_config.ini'])

with HDFArchive('out/inp.h5','r') as ar:
    G_iw = ar['DMFT_results/last_iter/Gimp_freq_0']
    G_iw_direct = ar['DMFT_results/last_iter/Gimp_freq_direct_0']


# compare direct measured G_iw
iw0 = len(G_iw.mesh)//2
n_iw = len(G_iw_direct.mesh) // 2

# works only if statistics are improved to higher precision
for block, gf in G_iw_direct:
    assert_arrays_are_close(G_iw[block].data[iw0-n_iw:iw0+n_iw, :, :], G_iw_direct[block].data, precision=5e-2)
