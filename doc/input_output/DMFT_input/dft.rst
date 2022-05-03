
[dft]: DFT related inputs
------

List of parameters that relate to the DFT calculation, useful mostly when doing CSC.

List of possible entries:


seedname; jobname; csc; plo_cfg; h_int_type; U; J; ratio_F4_F2; beta; n_iter_dmft_first; n_iter_dmft_per; n_iter_dmft; dc_type; prec_mu; dc_dmft; cpa_zeta; cpa_x; solver_type; n_iw; n_tau; n_l; n_w; w_range; eta; diag_delta; h5_save_freq; magnetic; magmom; enforce_off_diag; h_field; energy_shift_orbitals; sigma_mix; g0_mix; g0_mix_type; broy_max_it; dc; calc_energies; block_threshold; block_suppress_orbital_symm; load_sigma; path_to_sigma; load_sigma_iter; noise_level_initial_sigma; occ_conv_crit; gimp_conv_crit; g0_conv_crit; sigma_conv_crit; sampling_iterations; sampling_h5_save_freq; fixed_mu_value; mu_update_freq; dft_mu; mu_mix_const; mu_mix_per_occupation_offset; afm_order; set_rot; oneshot_postproc_gamma_file; measure_chi_SzSz; measure_chi_insertions; mu_gap_gb2_threshold; mu_gap_occ_deviation; 


.. admonition:: dft_code: 
 
            **type=** string

            Choose the DFT code interface, for now Quantum Espresso and Vasp are available.

            Possible values:

            * 'vasp'
            * 'qe'

.. admonition:: n_cores: 
 
            **type=** int

            number of cores for the DFT code (VASP)

.. admonition:: n_iter: 
 
            **type=** int;  **optional**;  **default=**  6

            number of dft iterations per cycle

.. admonition:: dft_exec: 
 
            **type=** string;  **default=**  'vasp_std'

            command for the DFT / VASP executable

.. admonition:: store_eigenvals: 
 
            **type=** bool;  **optional**;  **default=**  False

            stores the dft eigenvals from LOCPROJ (projector_type=plo) or
            wannier90.eig (projector_type=w90) file in h5 archive

.. admonition:: mpi_env: 
 
            **type=** string;  **default=**  'local'

            selection for mpi env for DFT / VASP in default this will only call VASP as mpirun -np n_cores_dft dft_exec

.. admonition:: projector_type: 
 
            **type=** string;  **optional**;  **default=**  'plo'

            plo: uses VASP's PLO formalism, requires LOCPROJ in the INCAR
            w90: uses Wannier90, requires LWANNIER90 in the INCAR

.. admonition:: w90_exec: 
 
            **type=** string;  **default=** 'wannier90.x'

            the command to start a single-core wannier run

.. admonition:: w90_tolerance: 
 
            **type=** float;  **default=** 1e-6

            threshold for mapping of shells and checks of the Hamiltonian
