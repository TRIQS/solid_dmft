import shutil
import triqs.utility.mpi as mpi
import solid_dmft.main as solid

if mpi.is_master_node():
    shutil.rmtree('out', ignore_errors=True)
mpi.barrier()

solid.main([None, 'dmft_config.toml'])
