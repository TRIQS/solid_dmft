from h5 import HDFArchive
import triqs.utility.mpi as mpi

from helper import are_iterables_equal

from solid_dmft.postprocessing import maxent_gf_imp, maxent_gf_latt, maxent_sigma, pade_sigma

# Note: no unit tests here because lack of MPI support

# Runs maxent on lattice Green function and compares afterwards
mpi.report('#########\nTesting lattice Gf Maxent\n#########')
maxent_gf_latt.main('svo_hubbardI_basic/out/inp.h5',
                    sum_spins=True,
                    n_points_maxent=100,
                    n_points_alpha=25,
                    omega_min=-20,
                    omega_max=20)

if mpi.is_master_node():
    print('Comparing Alatt_maxent')
    with HDFArchive('svo_hubbardI_basic/out/inp.h5', 'r') as out, HDFArchive('svo_hubbardI_basic/ref.h5', 'r') as ref:
        assert are_iterables_equal(out['DMFT_results']['last_iter']['Alatt_maxent'], ref['DMFT_results']['last_iter']['Alatt_maxent'])

# Runs maxent on the impurity Green function
# No comparison to reference because the spectral function would be too "spiky"
# because of HubbardI so that numerical differences can dominate the comparison
mpi.report('#########\nTesting impurity Gf Maxent\n#########')
maxent_gf_imp.main('svo_hubbardI_basic/out/inp.h5', sum_spins=True,
                   n_points_maxent=40, n_points_alpha=6, maxent_error=0.01)

# Run sigma maxent
mpi.report('#########\nTesting Sigma Maxent\n#########')
maxent_sigma.main(external_path='svo_hubbardI_basic/out/inp.h5',
                  omega_min=-12, omega_max=12,
                  maxent_error=0.002, iteration=None,
                  n_points_maxent=40,
                  n_points_alpha=6,
                  analyzer='LineFitAnalyzer',
                  n_points_interp=501,
                  n_points_final=501,
                  continuator_type='inversion_dc')[0]


# Run sigma pade
mpi.report('#########\nTesting Sigma Pade\n#########')
pade_sigma.main(external_path='svo_hubbardI_basic/out/inp.h5',
                n_w = 4001,
                w_min=-4.5,
                w_max=4.5,
                n_iw=100,
                eta=0.0
                )

