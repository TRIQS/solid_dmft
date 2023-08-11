from h5 import HDFArchive
import triqs.utility.mpi as mpi

from helper import are_iterables_equal

from solid_dmft.postprocessing import maxent_gf_imp, maxent_gf_latt

# Runs maxent on lattice Green function and compares afterwards
maxent_gf_latt.main('svo_hubbardI_basic/out/inp.h5', sum_spins=True, n_points_maxent=100, n_points_alpha=25, omega_min=-20, omega_max=20)

# Runs maxent on the impurity Green function
# No comparison to reference because the spectral function would be too "spiky"
# because of HubbardI so that numerical differences can dominate the comparison
maxent_gf_imp.main('svo_hubbardI_basic/out/inp.h5', sum_spins=True,
                   n_points_maxent=50, n_points_alpha=20)

if mpi.is_master_node():
    print('Comparing Alatt_maxent')
    with HDFArchive('svo_hubbardI_basic/out/inp.h5', 'r')['DMFT_results']['last_iter'] as out, HDFArchive('svo_hubbardI_basic/ref.h5', 'r')['DMFT_results']['last_iter'] as ref:
        assert are_iterables_equal(out['Alatt_maxent'], ref['Alatt_maxent'])
