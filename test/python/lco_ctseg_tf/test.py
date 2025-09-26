import sys
import shutil
import importlib.util

import triqs.utility.mpi as mpi

import solid_dmft.main as solid

# try triqs_ctseg import
ctseg = importlib.util.find_spec("triqs_ctseg") is not None
if not ctseg:
    mpi.report('ImportWarning: ctseg needs to be installed to run this test')
    sys.exit()

if mpi.is_master_node():
    shutil.rmtree('out', ignore_errors=True)

mpi.barrier()

solid.main([None, 'dmft_config.toml'])
