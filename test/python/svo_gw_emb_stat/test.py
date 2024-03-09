import sys
import shutil

from triqs.gf import *
from triqs.utility.comparison_tests import assert_block_gfs_are_close, assert_arrays_are_close
from h5 import HDFArchive
import triqs.utility.mpi as mpi

import solid_dmft.main as solid

try:
    from triqs_cthyb import Solver
except ImportError:
    print('ctseg solver not installed skipping')
    sys.exit()

if mpi.is_master_node():
    shutil.rmtree('out', ignore_errors=True)

mpi.barrier()

solid.main([None, 'dmft_config.toml'])
