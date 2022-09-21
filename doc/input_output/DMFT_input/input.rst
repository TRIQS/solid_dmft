DMFT input
------------------------

The aim of this section is to provide a comprehensive listing of all the input flags available for the `dmft_config.ini` input file. We begin by listing the possible sections and follow with the input parameters.

Input/Output
===================
.. toctree::
    :maxdepth: 1

    general
    solver
    dft
    advanced

Below an exhaustive list containing all the parameters marked by section.

**[general]**

seedname; jobname; csc; plo_cfg; h_int_type; U; J; ratio_F4_F2; beta; n_iter_dmft_first; n_iter_dmft_per; n_iter_dmft; dc_type; prec_mu; dc_dmft; cpa_zeta; cpa_x; solver_type; n_iw; n_tau; n_l; n_w; w_range; eta; diag_delta; h5_save_freq; magnetic; magmom; enforce_off_diag; h_field; energy_shift_orbitals; sigma_mix; g0_mix; g0_mix_type; broy_max_it; dc; calc_energies; block_threshold; block_suppress_orbital_symm; load_sigma; path_to_sigma; load_sigma_iter; noise_level_initial_sigma; occ_conv_crit; gimp_conv_crit; g0_conv_crit; sigma_conv_crit; sampling_iterations; sampling_h5_save_freq; fixed_mu_value; mu_update_freq; dft_mu; mu_mix_const; mu_mix_per_occupation_offset; afm_order; set_rot; oneshot_postproc_gamma_file; measure_chi_SzSz; measure_chi_insertions; mu_gap_gb2_threshold; mu_gap_occ_deviation; 

**[solver]**


store_solver; length_cycle; n_warmup_cycles; n_cycles_tot; measure_G_l; measure_G_tau; measure_G_iw; measure_density_matrix; measure_pert_order; max_time; imag_threshold; off_diag_threshold; delta_interface; move_double; perform_tail_fit; fit_max_moment; fit_min_n; fit_max_n; fit_min_w; fit_max_w; random_seed; legendre_fit; n_bath; bath_fit; refine_factor; ph_symm; calc_me; enforce_gap; ignore_weight; dt; state_storage; path_to_gs; sweeps; maxmI; maxmIB; maxmB; tw; dmrg_maxmI; dmrg_maxmIB; dmrg_maxmB; dmrg_tw; measure_hist; improved_estimator ; 

**[dft]**

dft_code; n_cores; n_iter; dft_exec; store_eigenvals; mpi_env; projector_type; w90_exec; w90_tolerance; 

**[advanced]**

dc_factor; dc_fixed_value; dc_fixed_occ; dc_U; dc_J; map_solver_struct; mapped_solver_struct_degeneracies; 