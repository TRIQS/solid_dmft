import sys
import shutil

from triqs.gfs import *
from triqs.utility.comparison_tests import assert_block_gfs_are_close
from h5 import HDFArchive
import triqs.utility.mpi as mpi

import solid_dmft.main as solid


if mpi.is_master_node():
    shutil.rmtree('out', ignore_errors=True)

mpi.barrier()

solid.main([None, 'dmft_config.toml'])

mpi.barrier()

if mpi.is_master_node():
    out = 'out/inp.h5'
    ref = 'ref.h5'

    out =  HDFArchive(out, 'r')['DMFT_results']['last_iter']
    ref =  HDFArchive(ref, 'r')['DMFT_results']['last_iter']

