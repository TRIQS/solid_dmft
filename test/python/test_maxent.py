from h5 import HDFArchive
import triqs.utility.mpi as mpi

from helper import are_iterables_equal

from solid_dmft.postprocessing import maxent_gf_imp, maxent_gf_latt

if mpi.is_master_node():
    maxent_gf_imp.main('svo_hubbardI_basic/out/inp.h5', sum_spins=True, n_points_maxent=100, n_points_alpha=25)
mpi.barrier()

maxent_gf_latt.main('svo_hubbardI_basic/out/inp.h5', sum_spins=True, n_points_maxent=100, n_points_alpha=25,
                    omega_min=-20, omega_max=20)
mpi.barrier()

if mpi.is_master_node():
    with HDFArchive('svo_hubbardI_basic/out/inp.h5', 'r')['DMFT_results']['last_iter'] as out, \
            HDFArchive('svo_hubbardI_basic/ref.h5', 'r')['DMFT_results']['last_iter'] as ref:
        for key in ['Aimp_maxent_0', 'Alatt_maxent']:
            are_iterables_equal(out[key], ref[key])
