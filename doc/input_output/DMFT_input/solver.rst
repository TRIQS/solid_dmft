[solver]: solver specific parameters
------

Here are the parameters that are uniquely dependent on the solver chosen. Below a list of the supported solvers:

===================
.. toctree::
    :maxdepth: 1



List of possible entries:

seedname; jobname; csc; plo_cfg; h_int_type; U; J; ratio_F4_F2; beta; n_iter_dmft_first; n_iter_dmft_per; n_iter_dmft; dc_type; prec_mu; dc_dmft; cpa_zeta; cpa_x; solver_type; n_iw; n_tau; n_l; n_w; w_range; eta; diag_delta; h5_save_freq; magnetic; magmom; enforce_off_diag; h_field; energy_shift_orbitals; sigma_mix; g0_mix; g0_mix_type; broy_max_it; dc; calc_energies; block_threshold; block_suppress_orbital_symm; load_sigma; path_to_sigma; load_sigma_iter; noise_level_initial_sigma; occ_conv_crit; gimp_conv_crit; g0_conv_crit; sigma_conv_crit; sampling_iterations; sampling_h5_save_freq; fixed_mu_value; mu_update_freq; dft_mu; mu_mix_const; mu_mix_per_occupation_offset; afm_order; set_rot; oneshot_postproc_gamma_file; measure_chi_SzSz; measure_chi_insertions; mu_gap_gb2_threshold; mu_gap_occ_deviation; 


.. admonition:: store_solver: 
 
            **type=** bool;  **optional** default= False

            store the whole solver object under DMFT_input in h5 archive

cthyb parameters
================

.. admonition:: length_cycle: 
 
            **type=** int

            length of each cycle; number of sweeps before measurement is taken

.. admonition:: n_warmup_cycles: 
 
            **type=** int

            number of warmup cycles before real measurement sets in

.. admonition:: n_cycles_tot: 
 
            **type=** int

            total number of sweeps

.. admonition:: measure_G_l: 
 
            **type=** bool

            measure Legendre Greens function

.. admonition:: measure_G_tau: 
 
            **type=** bool;  **optional**;  **default=** True

            should the solver measure G(tau)?

.. admonition:: measure_G_iw: 
 
            **type=** bool;  **optional**;  **default=** False

            should the solver measure G(iw)?

.. admonition:: measure_density_matrix: 
 
            **type=** bool;  **optional**;  **default=** False

            measures the impurity density matrix and sets also
            use_norm_as_weight to true

.. admonition:: measure_pert_order: 
 
            **type=** bool;  **optional**;  **default=** False

            measure perturbation order histograms: triqs.github.io/cthyb/latest/guide/perturbation_order_notebook.html

            The result is stored in the h5 archive under 'DMFT_results' at every iteration
            in the subgroups 'pert_order_imp_X' and 'pert_order_total_imp_X'

.. admonition:: max_time: 
 
            **type=** int;  **optional**;  **default=** -1

            maximum amount the solver is allowed to spend in each iteration

.. admonition:: imag_threshold: 
 
            **type=** float;  **optional**;  **default=**  10e-15

            threshold for imag part of G0_tau. be warned if symmetries are off in projection scheme imag parts can occur in G0_tau

.. admonition:: off_diag_threshold: 
 
            **type=** float;  **optional**

            threshold for off-diag elements in Hloc0

.. admonition:: delta_interface: 
 
            **type=** bool;  **optional**;  **default=** False

            use new delta interface in cthyb instead of input G0

.. admonition:: move_double: 
 
            **type=** bool;  **optional**;  **default=** True

            double moves in solver

.. admonition:: perform_tail_fit: 
 
            **type=** bool;  **optional**;  **default=** False

            tail fitting if legendre is off?

.. admonition:: fit_max_moment: 
 
            **type=** int;  **optional**

            max moment to be fitted

.. admonition:: fit_min_n: 
 
            **type=** int;  **optional**

            number of start matsubara frequency to start with

.. admonition:: fit_max_n: 
 
            **type=** int;  **optional**

            number of highest matsubara frequency to fit

.. admonition:: fit_min_w: 
 
            **type=** float;  **optional**

            start matsubara frequency to start with

.. admonition:: fit_max_w: 
 
            **type=** float;  **optional**

            highest matsubara frequency to fit

.. admonition:: random_seed: 
 
            **type=** int;  **optional** default by triqs

            if specified the int will be used for random seeds! Careful, this will give the same random
            numbers on all mpi ranks

.. admonition:: legendre_fit: 
 
            **type=** bool;  **optional** default= False

            filter noise of G(tau) with G_l, cutoff is taken from n_l

ftps parameters
===============

.. admonition:: n_bath: 
 
            **type=** int

            number of bath sites

.. admonition:: bath_fit: 
 
            **type=** bool;  **default=** False

            DiscretizeBath vs BathFitter

.. admonition:: refine_factor: 
 
            **type=** int;  **optional**;  **default=** 1

            rerun ftps cycle with increased accuracy

.. admonition:: ph_symm: 
 
            **type=** bool;  **optional**;  **default=** False

            particle-hole symmetric problem

.. admonition:: calc_me: 
 
            **type=** bool;  **optional**;  **default=** True

            calculate only symmetry-inequivalent spins/orbitals, symmetrized afterwards

.. admonition:: enforce_gap: 
 
            **type=** list of floats;  **optional**;  **default=** 'none'

            enforce gap in DiscretizeBath between interval

.. admonition:: ignore_weight: 
 
            **type=** float;  **optional**;  **default=** 0.0

            ignore weight of peaks for bath fitter

.. admonition:: dt: 
 
            **type=** float

            time step

.. admonition:: state_storage: 
 
            **type=** string;  **default=**  './'

            location of large MPS states

.. admonition:: path_to_gs: 
 
            **type=** string;  **default=**  'none'

            location of GS if already present. Use 'postprocess' to skip solver and go directly to post-processing
            of previously terminated time-evolved state

.. admonition:: sweeps: 
 
            **type=** int;  **optional**;  **default=**  10

            Number of DMRG sweeps

.. admonition:: maxmI: 
 
            **type=** int;  **optional**;  **default=**  100

            maximal imp-imp bond dimensions

.. admonition:: maxmIB: 
 
            **type=** int;  **optional**;  **default=**  100

            maximal imp-bath bond dimensions

.. admonition:: maxmB: 
 
            **type=** int;  **optional**;  **default=**  100

            maximal bath-bath bond dimensions

.. admonition:: tw: 
 
            **type=** float, default 1E-9

            truncated weight for every link

.. admonition:: dmrg_maxmI: 
 
            **type=** int;  **optional**;  **default=**  100

            maximal imp-imp bond dimensions

.. admonition:: dmrg_maxmIB: 
 
            **type=** int;  **optional**;  **default=**  100

            maximal imp-bath bond dimensions

.. admonition:: dmrg_maxmB: 
 
            **type=** int;  **optional**;  **default=**  100

            maximal bath-bath bond dimensions

.. admonition:: dmrg_tw: 
 
            **type=** float, default 1E-9

            truncated weight for every link

ctseg parameters
================

.. admonition:: measure_hist: 
 
            **type=** bool;  **optional**;  **default=** False

               measure perturbation_order histograms

.. admonition:: improved_estimator : 
 
            **type=** bool;  **optional**;  **default=** False

              measure improved estimators
              Sigma_iw will automatically be calculated via
              http://dx.doi.org/10.1103/PhysRevB.85.205106
