import sys
import shutil

from triqs.gf import *
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

    for key in ['Delta_time_0', 'G0_Refreq_0', 'G0_freq_0', 'Gimp_Refreq_0', 'Gimp_freq_0', 'Gimp_time_0', 'Sigma_Refreq_0', 'Sigma_freq_0',
                'Delta_time_1', 'G0_Refreq_1', 'G0_freq_1', 'Gimp_Refreq_1', 'Gimp_freq_1', 'Gimp_time_1', 'Sigma_Refreq_1', 'Sigma_freq_1']:
        print(key)
        assert_block_gfs_are_close(out[key],ref[key], precision= 1e-5)

    for key in ['chemical_potential_pre', 'chemical_potential_post']:
        print(key)
        assert abs(out[key]-ref[key]) < 0.001, "chemical potential mismatch"
