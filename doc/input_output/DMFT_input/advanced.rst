[advanced]: Advanced inputs
------

Advanced parameters, do not modify default value unless you know what you are doing

List of possible entries:

seedname; jobname; csc; plo_cfg; h_int_type; U; J; ratio_F4_F2; beta; n_iter_dmft_first; n_iter_dmft_per; n_iter_dmft; dc_type; prec_mu; dc_dmft; cpa_zeta; cpa_x; solver_type; n_iw; n_tau; n_l; n_w; w_range; eta; diag_delta; h5_save_freq; magnetic; magmom; enforce_off_diag; h_field; energy_shift_orbitals; sigma_mix; g0_mix; g0_mix_type; broy_max_it; dc; calc_energies; block_threshold; block_suppress_orbital_symm; load_sigma; path_to_sigma; load_sigma_iter; noise_level_initial_sigma; occ_conv_crit; gimp_conv_crit; g0_conv_crit; sigma_conv_crit; sampling_iterations; sampling_h5_save_freq; fixed_mu_value; mu_update_freq; dft_mu; mu_mix_const; mu_mix_per_occupation_offset; afm_order; set_rot; oneshot_postproc_gamma_file; measure_chi_SzSz; measure_chi_insertions; mu_gap_gb2_threshold; mu_gap_occ_deviation; 


.. admonition:: dc_factor: 
 
            **type=** float;  **optional**;  **default=**  'none' (corresponds to 1)

            If given, scales the dc energy by multiplying with this factor, usually < 1

.. admonition:: dc_fixed_value: 
 
            **type=** float;  **optional**;  **default=**  'none'

            If given, it sets the DC (energy/imp) to this fixed value. Overwrites EVERY other DC configuration parameter if DC is turned on

.. admonition:: dc_fixed_occ: 
 
            **type=** list of float;  **optional**;  **default=**  'none'

            If given, the occupation for the DC for each impurity is set to the provided value.
            Still uses the same kind of DC!

.. admonition:: dc_U: 
 
            **type=** float or comma seperated list of floats;  **optional**;  **default=**  general_params['U']

            U values for DC determination if only one value is given, the same U is assumed for all impurities

.. admonition:: dc_J: 
 
            **type=** float or comma seperated list of floats;  **optional**;  **default=**  general_params['J']

            J values for DC determination if only one value is given, the same J is assumed for all impurities

.. admonition:: map_solver_struct: 
 
            **type=** dict;  **optional**;  **default=** no additional mapping

            Additional manual mapping of the solver block structure, applied
            after the block structure finder to all impurities.

.. admonition:: mapped_solver_struct_degeneracies: 
 
            **type=** list;  **optional**;  **default=** none

