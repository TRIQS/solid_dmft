import shutil

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
    with HDFArchive('out/inp.h5', 'r')['DMFT_results']['last_iter'] as out, \
            HDFArchive('ref.h5', 'r')['DMFT_results']['last_iter'] as ref:
        for key in ['G0_freq_0', 'Gimp_freq_0', 'Gimp_time_0', 'Sigma_Refreq_0', 'Sigma_freq_0']:
            print(key)
            assert_block_gfs_are_close(out[key],ref[key],precision=1e-05)

        for key in ['chemical_potential_pre', 'chemical_potential_post']:
            print(key)
            assert abs(out[key]-ref[key]) < 0.001, "chemical potential mismatch"
