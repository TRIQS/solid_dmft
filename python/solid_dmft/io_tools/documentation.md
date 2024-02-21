List of all parameters, sorted by sections:

[  general  ]
-------------

seedname : str
            seedname for h5 archive with DMFT input and output
jobname : str, optional, default='dmft_dir'
            the output directory for one-shot calculations
csc : bool, optional, default=False
            are we doing a CSC calculation?
plo_cfg : str, optional, default='plo.cfg'
            config file for PLOs for the converter
h_int_type : string
            interaction type:

            * density_density: used for full d-shell or eg- or t2g-subset
            * kanamori: only physical for the t2g or the eg subset
            * full_slater: used for full d-shell or eg- or t2g-subset
            * ntot: U/2 (Ntot^2 - Ntot) interaction
            * simple_intra: density-density like but only intra orbital with given U value (no rotations applied)
            * crpa: use the cRPA matrix as interaction Hamiltonian
            * crpa_density_density: use the density-density terms of the cRPA matrix
            * dynamic: use dynamic U from h5 archive

            Needs to be stored as Matsubara Gf under dynamic_U/U_iw in the input h5
h_int_basis : string
              cubic basis convention to compute the interaction U matrix
              * 'triqs'
              * 'vasp' (equivalent to 'triqs')
              * 'wien2k'
              * 'wannier90'
              * 'qe' (equivalent to 'wannier90')
U :  float or comma separated list of floats
            U values for impurities if only one value is given, the same U is assumed for all impurities
U_prime :  float or comma separated list of floats
            U prime values for impurities if only one value is given, the same U prime is assumed for all impurities
            only used if h_int_type is kanamori
J :  float or comma separated list of floats
            J values for impurities if only one value is given, the same J is assumed for all impurities
ratio_F4_F2 : float or comma separated list of floats, optional, default='none'
            Ratio between the Slater integrals  F_4 and F_2. Only used for the
            interaction Hamiltonians 'density_density' and 'full_slater' and
            only for d-shell impurities, where the default is 0.63.
beta : float, only used if solver ImFreq
            inverse temperature for Greens function etc
n_iter_dmft_first : int, optional, default= 10
            number of iterations in first dmft cycle to converge dmft solution
n_iter_dmft_per : int, optional, default= 2
            number of iterations per dmft step in CSC calculations
n_iter_dmft : int
            number of iterations per dmft cycle after first cycle
dc_type : int
            Type of double counting correction considered:
            * 0: FLL
            * 1: held formula, needs to be used with slater-kanamori h_int_type=2
            * 2: AMF
            * 3: FLL for eg orbitals only with U,J for Kanamori
dc_dmft : bool
           Whether to use DMFT or DFT occupations:

           * DC with DMFT occupation in each iteration -> True
           * DC with DFT occupations after each DFT cycle -> False
cpa_zeta : float or comma separated list of floats
            shift of local levels per impurity in CPA
cpa_x : float or comma separated list of floats
            probability distribution for summing G(tau) in CPA
solver_type : str
            type of solver chosen for the calculation, currently supports:

            * 'cthyb'
            * 'ctint'
            * 'ftps'
            * 'hubbardI'
            * 'hartree'
            * 'ctseg'

n_iw : int, optional, default=1025
            number of Matsubara frequencies
n_tau : int, optional, default=10001
            number of imaginary time points
n_l : int, needed if measure_G_l=True or legendre_fit=True
            number of Legendre coefficients
n_w : int, optional, default=5001
            number of real frequency points
w_range : tuple, optional, default=(-10, 10)
            w_min and w_max, example: w_range = -10, 10
eta : float, only used if solver ReFreq
            broadening of Green's function
diag_delta : bool, optional, default=False
            option to remove off-diagonal terms in the hybridization function


h5_save_freq : int, optional, default=5
            how often is the output saved to the h5 archive
magnetic : bool, optional, default=False
            are we doing a magnetic calculations? If yes put magnetic to True.
            Not implemented for CSC calculations
magmom : list of float seperated by comma, optional default=[]
            Initialize magnetic moments if magnetic is on. length must be #imps.
            List composed of energetic shifts written in electronvolts.
            This will initialize the spin blocks of the sigma with a diagonal shift
            With -shift for the up block, and +shift for the down block
            (positive shift favours the up spin component, not compatible with spin-orbit coupling)
enforce_off_diag : bool, optional, default=False
            enforce off diagonal elements in block structure finder
h_field : float, optional, default=0.0
            magnetic field
h_field_it : int, optional, default=0
            number of iterations the magnetic field is kept on
sigma_mix : float, optional, default=1.0
            careful: Sigma mixing can break orbital symmetries, use G0 mixing
            mixing sigma with previous iteration sigma for better convergency. 1.0 means no mixing
g0_mix : float, optional, default=1.0
            Mixing the weiss field G0 with previous iteration G0 for better convergency. 1.0 means no mixing.
            Setting g0_mix to 0.0 with linear mixing can be used for statistic sampling when
            restarting a calculation
g0_mix_type : string, optional, default='linear'
            which type of mixing is used. Possible values are:
            linear: linear mixing
            broyden: broyden mixing
broy_max_it : int, optional, default=1
            maximum number of iteration to be considered for broyden mixing
            1 corresponds to simple linear mixing
dc : bool, optional, default=True
            dc correction on yes or no?
calc_energies : bool, optional, default=False, not compatible with 'ftps' solver
            calc energies explicitly within the dmft loop
block_threshold : float, optional, default=1e-05
            threshold for finding block structures in the input data (off-diag yes or no)
block_suppress_orbital_symm : bool, optional, default=False
            should blocks be checked if symmetry-equiv. between orbitals?
            Does not affect spin symmetries.
load_sigma : bool, optional, default=False
            load a old sigma from h5 file
path_to_sigma : str, needed if load_sigma is true
            path to h5 file from which the sigma should be loaded
load_sigma_iter : int, optional, default= last iteration
            load the sigma from a specific iteration if wanted
noise_level_initial_sigma : float, optional, default=0.0
            spread of Gaussian noise applied to the initial Sigma
occ_conv_crit : float, optional, default= -1
            stop the calculation if a certain threshold for the imp occ change is reached
gimp_conv_crit : float, optional, default= -1
            stop the calculation if  sum_w 1/(w^0.6) ||Gimp-Gloc|| is smaller than threshold
g0_conv_crit : float, optional, default= -1
            stop the calculation if sum_w 1/(w^0.6) ||G0-G0_prev|| is smaller than threshold
sigma_conv_crit : float, optional, default= -1
            stop the calculation if sum_w 1/(w^0.6) ||Sigma-Sigma_prev|| is smaller than threshold
sampling_iterations : int, optional, default= 0
            for how many iterations should the solution sampled after the CSC loop is converged
sampling_h5_save_freq : int, optional, default= 5
            overwrites h5_save_freq when sampling has started
calc_mu_method : string, optional, default = 'dichotomy'
            optimization method used for finding the chemical potential:

            * 'dichotomy': usual method from TRIQS, should always converge but may be slow
            * 'newton': scipy Newton root finder, much faster but might be unstable
            * 'brent': scipy hyperbolic Brent root finder preconditioned with dichotomy to find edge, a compromise between speed and stability
prec_mu : float
            general precision for determining the chemical potential at any time calc_mu is called
fixed_mu_value : float, optional, default= 'none'
            If given, the chemical potential remains fixed in calculations
mu_update_freq : int, optional, default= 1
            The chemical potential will be updated every # iteration
mu_initial_guess : float, optional, default= 'none'
            The chemical potential of the DFT calculation.
            If not given, mu will be calculated from the DFT bands
mu_mix_const : float, optional, default= 1.0
            Constant term of the mixing of the chemical potential. See mu_mix_per_occupation_offset.
mu_mix_per_occupation_offset : float, optional, default= 0.0
            Mu mixing proportional to the occupation offset.
            Mixing between the dichotomy result and the previous mui,

            mu_next = factor * mu_dichotomy + (1-factor) * mu_previous, with
            factor = mu_mix_per_occupation_offset * abs(n - n\_target) + mu_mix_const.

            The program ensures that 0 <= factor <= 1.
            mu_mix_const = 1.0 and mu_mix_per_occupation_offset = 0.0 means no mixing.
afm_order : bool, optional, default=False
            copy self energies instead of solving explicitly for afm order
set_rot : string, optional, default='none'
            use density_mat_dft to diagonalize occupations = 'den'
            use hloc_dft to diagonalize occupations = 'hloc'
measure_chi_SzSz : bool, optional, default=False
            measure the dynamic spin suszeptibility chi(sz,sz(tau))
            triqs.github.io/cthyb/unstable/guide/dynamic_susceptibility_notebook.html
measure_chi_insertions : int, optional, default=100
            number of insertation for measurement of chi
mu_gap_gb2_threshold : float, optional, default=none
            Threshold of the absolute of the lattice GF at tau=beta/2 for use
            of MaxEnt's lattice spectral function to put the chemical potential
            into the middle of the gap. Does not work if system completely full
            or empty, mu mixing is not applied to it. Recommended value 0.01.
mu_gap_occ_deviation : float, optional, default=none
            Only used if mu_gap_gb2_threshold != none. Sets additional criterion
            for finding the middle of the gap through occupation deviation to
            avoid getting stuck in an insulating state with wrong occupation.

[  solver  ]
------------
store_solver : bool, optional default= False
            store the whole solver object under DMFT_input in h5 archive

cthyb parameters
================
length_cycle : int
            length of each cycle; number of sweeps before measurement is taken
n_warmup_cycles : int
            number of warmup cycles before real measurement sets in
n_cycles_tot : int
            total number of sweeps
measure_G_l : bool
            measure Legendre Greens function
measure_G_tau : bool,optional, default=True
            should the solver measure G(tau)?
measure_G_iw : bool,optional, default=False
            should the solver measure G(iw)?
measure_density_matrix : bool, optional, default=False
            measures the impurity density matrix and sets also
            use_norm_as_weight to true
measure_pert_order : bool, optional, default=False
            measure perturbation order histograms: triqs.github.io/cthyb/latest/guide/perturbation_order_notebook.html

            The result is stored in the h5 archive under 'DMFT_results' at every iteration
            in the subgroups 'pert_order_imp_X' and 'pert_order_total_imp_X'
max_time : int, optional, default=-1
            maximum amount the solver is allowed to spend in each iteration
imag_threshold : float, optional, default= 10e-15
            threshold for imag part of G0_tau. be warned if symmetries are off in projection scheme imag parts can occur in G0_tau
off_diag_threshold : float, optional
            threshold for off-diag elements in Hloc0
delta_interface : bool, optional, default=False
            use new delta interface in cthyb instead of input G0
move_double : bool, optional, default=True
            double moves in solver
perform_tail_fit : bool, optional, default=False
            tail fitting if legendre is off?
fit_max_moment : int, optional
            max moment to be fitted
fit_min_n : int, optional
            number of start matsubara frequency to start with
fit_max_n : int, optional
            number of highest matsubara frequency to fit
fit_min_w : float, optional
            start matsubara frequency to start with
fit_max_w : float, optional
            highest matsubara frequency to fit
random_seed : str, optional default by triqs
            if specified the int will be used for random seeds! Careful, this will give the same random
            numbers on all mpi ranks
            You can also pass a string that will convert the keywords it or rank on runtime, e.g.
            34788 * it + 928374 * rank will convert each iteration the variables it and rank for the random
            seed
legendre_fit : bool, optional default= False
            filter noise of G(tau) with G_l, cutoff is taken from n_l
loc_n_min : int, optional
            Restrict local Hilbert space to states with at least this number of particles
loc_n_max : int, optional
            Restrict local Hilbert space to states with at most this number of particles

ftps parameters
===============
n_bath :     int
            number of bath sites
bath_fit :   bool, default=False
            DiscretizeBath vs BathFitter
refine_factor : int, optional, default=1
            rerun ftps cycle with increased accuracy
ph_symm :    bool, optional, default=False
            particle-hole symmetric problem
calc_me :    bool, optional, default=True
            calculate only symmetry-inequivalent spins/orbitals, symmetrized afterwards
enforce_gap : list of floats, optional, default='none'
            enforce gap in DiscretizeBath between interval
ignore_weight : float, optional, default=0.0
            ignore weight of peaks for bath fitter
dt :         float
            time step
state_storage : string, default= './'
            location of large MPS states
path_to_gs : string, default= 'none'
            location of GS if already present. Use 'postprocess' to skip solver and go directly to post-processing
            of previously terminated time-evolved state
sweeps :     int, optional, default= 10
            Number of DMRG sweeps
maxmI :      int, optional, default= 100
            maximal imp-imp bond dimensions
maxmIB :     int, optional, default= 100
            maximal imp-bath bond dimensions
maxmB :      int, optional, default= 100
            maximal bath-bath bond dimensions
tw :         float, default 1E-9
            truncated weight for every link
dmrg_maxmI : int, optional, default= 100
            maximal imp-imp bond dimensions
dmrg_maxmIB : int, optional, default= 100
            maximal imp-bath bond dimensions
dmrg_maxmB : int, optional, default= 100
            maximal bath-bath bond dimensions
dmrg_tw :    float, default 1E-9
            truncated weight for every link

ctseg parameters
================
measure_hist : bool, optional, default=False
               measure perturbation_order histograms
improved_estimator  : bool, optional, default=False
              measure improved estimators
              Sigma_iw will automatically be calculated via
              http://dx.doi.org/10.1103/PhysRevB.85.205106

hartree parameters
================
with_fock : bool, optional, default=False
        include Fock exchange terms in the self-energy
force_real : bool, optional, default=True
        force the self energy from Hartree fock to be real
one_shot : bool, optional, default=True
        Perform a one-shot or self-consitent root finding in each DMFT step of the Hartree solver.
method : bool, optional, default=True
        method for root finder. Only used if one_shot=False, see scipy.optimize.root for options.
tol : float, optional, default=1e-5
        tolerance for root finder if one_shot=False.

[  dft  ]
---------
dft_code : string
            Choose the DFT code interface, for now Quantum Espresso and Vasp are available.

            Possible values:

            * 'vasp'
            * 'qe'
n_cores : int
            number of cores for the DFT code (VASP)
n_iter : int, optional, default= 6
            only needed for VASP. Number of DFT iterations to feed the DMFT
            charge density into DFT, which generally takes multiple Davidson steps.
            For every DFT iterations, the charge-density correction is recalculated
            using newly generated projectors and hoppings from the previous DFT run
n_iter_first : int, optional, default= dft/n_iter
            number of DFT iterations in the first charge correction because this
            first charge correction usually changes the DFT wave functions the most.
dft_exec :  string, default= 'vasp_std'
            command for the DFT executable
store_eigenvals : bool, optional, default= False
            stores the dft eigenvals from LOCPROJ (projector_type=plo) or
            wannier90.eig (projector_type=w90) file in h5 archive
mpi_env : string, default= 'local'
            selection for mpi env for DFT / VASP in default this will only call VASP as mpirun -np n_cores_dft dft_exec
projector_type : string, optional, default= 'w90'
            plo: uses VASP's PLO formalism, requires LOCPROJ in the INCAR
            w90: uses Wannier90 (for VASP and QuantumEspresso)
w90_exec :  string, default='wannier90.x'
            the command to start a single-core wannier run
w90_tolerance :  float, default=1e-6
            threshold for mapping of shells and checks of the Hamiltonian

[  advanced  ]
--------------
dc_factor : float, optional, default= 'none' (corresponds to 1)
            If given, scales the dc energy by multiplying with this factor, usually < 1
dc_fixed_value : float, optional, default= 'none'
            If given, it sets the DC (energy/imp) to this fixed value. Overwrites EVERY other DC configuration parameter if DC is turned on
dc_fixed_occ : list of float, optional, default= 'none'
            If given, the occupation for the DC for each impurity is set to the provided value.
            Still uses the same kind of DC!
dc_orb_shift : list of float, optional, default= 'none'
            extra potential shift per orbital per impurity added to the DC
dc_U :  float or comma seperated list of floats, optional, default= general_params['U']
            U values for DC determination if only one value is given, the same U is assumed for all impurities
dc_J :  float or comma seperated list of floats, optional, default= general_params['J']
            J values for DC determination if only one value is given, the same J is assumed for all impurities
map_solver_struct : list of dict, optional, default=no additional mapping
            Additional manual mapping of the solver block structure, applied
            after the block structure finder for each impurity.
            Give exactly one dict per ineq impurity.
            see also triqs.github.io/dft_tools/latest/_python_api/triqs_dft_tools.block_structure.BlockStructure.map_gf_struct_solver.html
mapped_solver_struct_degeneracies : list, optional, default=none
            Degeneracies applied when using map_solver_struct, for each impurity.
            If not given and map_solver_struct is used, no symmetrization will happen.
pick_solver_struct : list of dict, optional, default=no additional picking
            input a solver dictionary for each ineq impurity to reduce dimensionality of
            solver block structure. Similar to to map_solver_struct, but with simpler syntax.
            Not listed blocks / orbitals will be not treated in impurity solver.
            Keeps degenerate shells.
