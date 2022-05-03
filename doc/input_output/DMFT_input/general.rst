[general]: General parameters
------

Includes the majority of the parameters

List of possible entries:

seedname; jobname; csc; plo_cfg; h_int_type; U; J; ratio_F4_F2; beta; n_iter_dmft_first; n_iter_dmft_per; n_iter_dmft; dc_type; prec_mu; dc_dmft; cpa_zeta; cpa_x; solver_type; n_iw; n_tau; n_l; n_w; w_range; eta; diag_delta; h5_save_freq; magnetic; magmom; enforce_off_diag; h_field; energy_shift_orbitals; sigma_mix; g0_mix; g0_mix_type; broy_max_it; dc; calc_energies; block_threshold; block_suppress_orbital_symm; load_sigma; path_to_sigma; load_sigma_iter; noise_level_initial_sigma; occ_conv_crit; gimp_conv_crit; g0_conv_crit; sigma_conv_crit; sampling_iterations; sampling_h5_save_freq; fixed_mu_value; mu_update_freq; dft_mu; mu_mix_const; mu_mix_per_occupation_offset; afm_order; set_rot; oneshot_postproc_gamma_file; measure_chi_SzSz; measure_chi_insertions; mu_gap_gb2_threshold; mu_gap_occ_deviation; 



.. admonition:: seedname: 
 
            **type=** str or list of str

            seedname for h5 archive or for multiple if calculations should be connected

.. admonition:: jobname: 
 
            **type=** str or list of str;  **optional**;  **default=** seedname

            one or multiple jobnames specifying the output directories

.. admonition:: csc: 
 
            **type=** bool;  **optional**;  **default=** False

            are we doing a CSC calculation?

.. admonition:: plo_cfg: 
 
            **type=** str;  **optional**;  **default=** 'plo.cfg'

            config file for PLOs for the converter

.. admonition:: h_int_type: 
 
            **type=** string

            interaction type:

            * density_density: used for full d-shell or eg- or t2g-subset
            * kanamori: only physical for the t2g or the eg subset
            * full_slater: used for full d-shell or eg- or t2g-subset
            * crpa: use the cRPA matrix as interaction Hamiltonian
            * crpa_density_density: use the density-density terms of the cRPA matrix
            * dynamic: use dynamic U from h5 archive
            
            Needs to be stored as Matsubara Gf under dynamic_U/U_iw in the input h5

.. admonition:: U: 
 
            **type=** float or comma separated list of floats

            U values for impurities if only one value is given, the same U is assumed for all impurities

.. admonition:: J: 
 
            **type=** float or comma separated list of floats

            J values for impurities if only one value is given, the same J is assumed for all impurities

.. admonition:: ratio_F4_F2: 
 
            **type=** float or comma separated list of floats;  **optional**;  **default=** 'none'

            Ratio between the Slater integrals  F_4 and F_2. Only used for the
            interaction Hamiltonians 'density_density' and 'full_slater' and
            only for d-shell impurities, where the default is 0.63.

.. admonition:: beta: 
 
            **type=** float, only used if solver ImFreq

            inverse temperature for Greens function etc

.. admonition:: n_iter_dmft_first: 
 
            **type=** int;  **optional**;  **default=**  10

            number of iterations in first dmft cycle to converge dmft solution

.. admonition:: n_iter_dmft_per: 
 
            **type=** int;  **optional**;  **default=**  2

            number of iterations per dmft step in CSC calculations

.. admonition:: n_iter_dmft: 
 
            **type=** int

            number of iterations per dmft cycle after first cycle

.. admonition:: dc_type: 
 
            **type=** int

            Type of double counting correction considered:
            * 0: FLL
            * 1: held formula, needs to be used with slater-kanamori h_int_type=2
            * 2: AMF
            * 3: FLL for eg orbitals only with U,J for Kanamori

.. admonition:: prec_mu: 
 
            **type=** float

            general precision for determining the chemical potential at any time calc_mu is called

.. admonition:: dc_dmft: 
 
            **type=** bool

           Whether to use DMFT or DFT occupations:

           * DC with DMFT occupation in each iteration -> True
           * DC with DFT occupations after each DFT cycle -> False

.. admonition:: cpa_zeta: 
 
            **type=** float or comma separated list of floats

            shift of local levels per impurity in CPA

.. admonition:: cpa_x: 
 
            **type=** float or comma separated list of floats

            probability distribution for summing G(tau) in CPA

.. admonition:: solver_type: 
 
            **type=** str

            type of solver chosen for the calculation, currently supports:

            * 'cthyb'
            * 'ctint'
            * 'ftps'
            * 'hubbardI'
            * 'ctseg'


.. admonition:: n_iw: 
 
            **type=** int;  **optional**;  **default=** 1025

            number of Matsubara frequencies

.. admonition:: n_tau: 
 
            **type=** int;  **optional**;  **default=** 10001

            number of imaginary time points

.. admonition:: n_l: 
 
            **type=** int, needed if measure_G_l=True or legendre_fit=True

            number of Legendre coefficients

.. admonition:: n_w: 
 
            **type=** int;  **optional**;  **default=** 5001

            number of real frequency points

.. admonition:: w_range: 
 
            **type=** tuple;  **optional**;  **default=** (-10, 10)

            w_min and w_max, example: w_range = -10, 10

.. admonition:: eta: 
 
            **type=** float, only used if solver ReFreq

            broadening of Green's function

.. admonition:: diag_delta: 
 
            **type=** bool;  **optional**;  **default=** False

            option to remove off-diagonal terms in the hybridization function



.. admonition:: h5_save_freq: 
 
            **type=** int;  **optional**;  **default=** 5

            how often is the output saved to the h5 archive

.. admonition:: magnetic: 
 
            **type=** bool;  **optional**;  **default=** False

            are we doing a magnetic calculations? If yes put magnetic to True.
            Not implemented for CSC calculations

.. admonition:: magmom: 
 
            **type=** list of float seperated by comma;  **optional** default=[]

            initialize magnetic moments if magnetic is on. length must be #imps.
            This will be used as factor for each imp in the initial self
            energy, with up (or ud for spin-orbit coupling) (1+fac)*sigma, and
            with down (1-fac)*sigma

.. admonition:: enforce_off_diag: 
 
            **type=** bool;  **optional**;  **default=** False

            enforce off diagonal elements in block structure finder

.. admonition:: h_field: 
 
            **type=** float;  **optional**;  **default=** 0.0

            magnetic field

.. admonition:: energy_shift_orbitals: 
 
            **type=** list of floats;  **optional**;  **default=**  'none'

            orbitals will be shifted by this energy
            The entries can be python code, to be combined with configparser's interpolation

.. admonition:: sigma_mix: 
 
            **type=** float;  **optional**;  **default=** 1.0

            careful: Sigma mixing can break orbital symmetries, use G0 mixing
            mixing sigma with previous iteration sigma for better convergency. 1.0 means no mixing

.. admonition:: g0_mix: 
 
            **type=** float;  **optional**;  **default=** 1.0

            mixing the weiss field G0 with previous iteration G0 for better convergency. 1.0 means no mixing

.. admonition:: g0_mix_type: 
 
            **type=** string;  **optional**;  **default=** 'linear'

            which type of mixing is used. Possible values are:
            linear: linear mixing
            broyden: broyden mixing

.. admonition:: broy_max_it: 
 
            **type=** int;  **optional**;  **default=** 1

            maximum number of iteration to be considered for broyden mixing
            1 corresponds to simple linear mixing

.. admonition:: dc: 
 
            **type=** bool;  **optional**;  **default=** True

            dc correction on yes or no?

.. admonition:: calc_energies: 
 
            **type=** bool;  **optional**;  **default=** False, not compatible with 'ftps' solver

            calc energies explicitly within the dmft loop

.. admonition:: block_threshold: 
 
            **type=** float;  **optional**;  **default=** 1e-05

            threshold for finding block structures in the input data (off-diag yes or no)

.. admonition:: block_suppress_orbital_symm: 
 
            **type=** bool;  **optional**;  **default=** False

            should blocks be checked if symmetry-equiv. between orbitals?
            Does not affect spin symmetries.

.. admonition:: load_sigma: 
 
            **type=** bool;  **optional**;  **default=** False

            load a old sigma from h5 file

.. admonition:: path_to_sigma: 
 
            **type=** str, needed if load_sigma is true

            path to h5 file from which the sigma should be loaded

.. admonition:: load_sigma_iter: 
 
            **type=** int;  **optional**;  **default=**  last iteration

            load the sigma from a specific iteration if wanted

.. admonition:: noise_level_initial_sigma: 
 
            **type=** float;  **optional**;  **default=** 0.0

            spread of Gaussian noise applied to the initial Sigma

.. admonition:: occ_conv_crit: 
 
            **type=** float;  **optional**;  **default=**  -1

            stop the calculation if a certain threshold for the imp occ change is reached

.. admonition:: gimp_conv_crit: 
 
            **type=** float;  **optional**;  **default=**  -1

            stop the calculation if  sum_w 1/(w^0.6) ||Gimp-Gloc|| is smaller than threshold

.. admonition:: g0_conv_crit: 
 
            **type=** float;  **optional**;  **default=**  -1

            stop the calculation if sum_w 1/(w^0.6) ||G0-G0_prev|| is smaller than threshold

.. admonition:: sigma_conv_crit: 
 
            **type=** float;  **optional**;  **default=**  -1

            stop the calculation if sum_w 1/(w^0.6) ||Sigma-Sigma_prev|| is smaller than threshold

.. admonition:: sampling_iterations: 
 
            **type=** int;  **optional**;  **default=**  0

            for how many iterations should the solution sampled after the CSC loop is converged

.. admonition:: sampling_h5_save_freq: 
 
            **type=** int;  **optional**;  **default=**  5

            overwrites h5_save_freq when sampling has started

.. admonition:: fixed_mu_value: 
 
            **type=** float;  **optional**;  **default=**  'none'

            If given, the chemical potential remains fixed in calculations

.. admonition:: mu_update_freq: 
 
            **type=** int;  **optional**;  **default=**  1

            The chemical potential will be updated every # iteration

.. admonition:: dft_mu: 
 
            **type=** float;  **optional**;  **default=**  'none'

            The chemical potential of the DFT calculation.
            If not given, mu will be calculated from the DFT bands

.. admonition:: mu_mix_const: 
 
            **type=** float;  **optional**;  **default=**  1.0

            Constant term of the mixing of the chemical potential. See mu_mix_per_occupation_offset.

.. admonition:: mu_mix_per_occupation_offset: 
 
            **type=** float;  **optional**;  **default=**  0.0

            Mu mixing proportional to the occupation offset.
            Mixing between the dichotomy result and the previous mui,

            mu_next = factor * mu_dichotomy + (1-factor) * mu_previous, with
            factor = mu_mix_per_occupation_offset * abs(n - n\_target) + mu_mix_const.

            The program ensures that 0 <= factor <= 1.
            mu_mix_const = 1.0 and mu_mix_per_occupation_offset = 0.0 means no mixing.

.. admonition:: afm_order: 
 
            **type=** bool;  **optional**;  **default=** False

            copy self energies instead of solving explicitly for afm order

.. admonition:: set_rot: 
 
            **type=** string;  **optional**;  **default=** 'none'

            use density_mat_dft to diagonalize occupations = 'den'
            use hloc_dft to diagonalize occupations = 'hloc'

.. admonition:: oneshot_postproc_gamma_file: 
 
            **type=** bool;  **optional**;  **default=** False

            write the GAMMA file for vasp after completed one-shot calculations

.. admonition:: measure_chi_SzSz: 
 
            **type=** bool;  **optional**;  **default=** False

            measure the dynamic spin suszeptibility chi(sz,sz(tau))
            triqs.github.io/cthyb/unstable/guide/dynamic_susceptibility_notebook.html

.. admonition:: measure_chi_insertions: 
 
            **type=** int;  **optional**;  **default=** 100

            number of insertation for measurement of chi

.. admonition:: mu_gap_gb2_threshold: 
 
            **type=** float;  **optional**;  **default=** none

            Threshold of the absolute of the lattice GF at tau=beta/2 for use
            of MaxEnt's lattice spectral function to put the chemical potential
            into the middle of the gap. Does not work if system completely full
            or empty, mu mixing is not applied to it. Recommended value 0.01.

.. admonition:: mu_gap_occ_deviation: 
 
            **type=** float;  **optional**;  **default=** none

            Only used if mu_gap_gb2_threshold != none. Sets additional criterion
            for finding the middle of the gap through occupation deviation to
            avoid getting stuck in an insulating state with wrong occupation.
