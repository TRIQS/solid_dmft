################################################################################
#
# solid_dmft - A versatile python wrapper to perform DFT+DMFT calculations
#              utilizing the TRIQS software library
#
# Copyright (C) 2018-2020, ETH Zurich
# Copyright (C) 2021, The Simons Foundation
#      authors: A. Hampel, M. Merkel, and S. Beck
#
# solid_dmft is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# solid_dmft is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# solid_dmft (in the file COPYING.txt in this directory). If not, see
# <http://www.gnu.org/licenses/>.
#
################################################################################
"""
Provides the read_config function to read the config file

Reads the config file (default dmft_config.ini) with python's configparser
module. It consists of at least the section 'general', and optionally of the
sections 'solver', 'dft' and 'advanced'.
Comments are in general possible with with the delimiters ';' or
'#'. However, this is only possible in the beginning of a line not within
the line! For default values, the string 'none' is used. NoneType cannot be
saved in an h5 archive (in the framework that we are using).

List of all parameters, sorted by sections:

---XXX---start
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
            * 'w2dyn_cthyb'

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

w2dyn_cthyb parameters
======================
measure_G2_iw_ph : bool, optional, default=False
                   Measure two-particle Green's function in particle-hole
                   frequency convention
measure_G2_n_bosonic : int, optional, default=30
                       Number of bosonic Matsubara frequencies for two-particle
                       Green's function measurement
measure_G2_n_fermionic : int, optional, default=30
                         Number of fermionic Matsubara frequencies for
                         two-particle Green's function measurement
worm_components : iterable of int, optional, default=None
                  Overrides worm components to measure with compound indices
                  (see ``w2dyn.auxiliaries.compound_index``) in iterable if not
                  None (Advanced setting: do not change if unsure)
move_global_prob : float, optional, default=0.005
                          Overall probability of the global moves
flavourchange_moves : int, optional, default=0
                      Use flavourchange-moves instead of 4-operator insertions;
                      they seemed to be inferior to 4 operator moves...
statesampling : int, optional, default=0
                Activate state-sampling algorithm

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


---XXX---end

"""

from configparser import ConfigParser
import triqs.utility.mpi as mpi
import numpy as np

# Workaround to get the default configparser boolean converter
BOOL_PARSER = lambda b: ConfigParser()._convert_to_boolean(b)

# TODO: it might be nicer to not have optional parameters at all and instead use
#       explicit default values

# Dictionary for the parameters. Contains the four sections general, dft, solver
# and advanced. Inside, all parameters are listed with their properties:
#   - converter: converter applied on the string value of the parameter
#   - valid for: a criterion for validity. If not fulfilled, the program crashes.
#                Always of form lambda x, params: ..., with x being the current parameter
#   - used: determines if parameter is used (and if not given, set to default value)
#           or unused and ignored. If 'used' and no default given, the program crashes.
#           If 'used' and default=None, this is an optional parameter
#   - default: default value for parameter. Can be a function of params but can only
#              use values that have NO default value. If it is None but 'used'
#              is True, the parameter becomes an optional parameter
PROPERTIES_PARAMS = {'general': {'seedname': {'used': True},

                                 'h_int_type': {'valid for': lambda x, _: all(hint in ('density_density',
                                                                                       'kanamori',
                                                                                       'full_slater',
                                                                                       'crpa',
                                                                                       'crpa_density_density',
                                                                                       'dynamic',
                                                                                       'ntot') for hint in x),
                                                'converter': lambda s: list(map(str, s.replace(" ", "").split(','))),
                                                'used': True},

                                 'U': {'converter': lambda s: list(map(float, s.split(','))), 'used': True},

                                 'U_prime': {'converter': lambda s: list(map(float, s.split(','))),
                                             'default': ['U-2J'],
                                             'valid for': lambda x, params: all(r == 'U-2J' or hint in ('kanamori')
                                                                                    for r, hint in zip(x, params['general']['h_int_type'])),
                                             'used': True},

                                 'J': {'converter': lambda s: list(map(float, s.split(','))), 'used': True},

                                 'ratio_F4_F2': {'converter': lambda s: list(map(float, s.split(','))),
                                                 'default': ['none'],
                                                 'valid for': lambda x, params: all(r == 'none' or hint in ('density_density','full_slater')
                                                                                    for r, hint in zip(x, params['general']['h_int_type'])),
                                                 'used': True},

                                 'beta': {'converter': float, 'valid for': lambda x, _: x > 0,
                                          'used': lambda params: params['general']['solver_type'] in ['cthyb', 'ctint', 'inchworm', 'hubbardI', 'ctseg', 'w2dyn_cthyb','hartree']},

                                 'n_iter_dmft': {'converter': int, 'valid for': lambda x, _: x >= 0, 'used': True},

                                 'dc': {'converter': BOOL_PARSER, 'used': True, 'default': True},

                                 'dc_type': {'converter': int, 'valid for': lambda x, _: x in (0, 1, 2, 3, 4),
                                             'used': lambda params: params['general']['dc']},

                                 'prec_mu': {'converter': float, 'valid for': lambda x, _: x > 0, 'used': True},

                                 'dc_dmft': {'converter': BOOL_PARSER,
                                             'used': lambda params: params['general']['dc']},

                                 'cpa_zeta': {'converter': lambda s: list(map(float, s.split(','))),
                                              'used': lambda params: params['general']['dc'] and params['general']['dc_type'] == 4},

                                 'cpa_x': {'converter': lambda s: list(map(float, s.split(','))),
                                           'used': lambda params: params['general']['dc'] and params['general']['dc_type'] == 4},

                                 'solver_type': {'valid for': lambda x, _: x in ['cthyb', 'ctint', 'ftps', 'hubbardI', 'ctseg', 'w2dyn_cthyb', 'hartree'],
                                                 'used': True},
                                 

                                 'n_l': {'converter': int, 'valid for': lambda x, _: x > 0,
                                         'used': lambda params: params['general']['solver_type'] in ['cthyb', 'inchworm', 'hubbardI', 'ctseg', 'w2dyn_cthyb']
                                         and (params['solver']['measure_G_l'] or params['solver']['legendre_fit'])},

                                 'n_iw': {'converter': int, 'valid for': lambda x, _: x > 0,
                                          'used': lambda params: params['general']['solver_type'] in ['cthyb', 'ctint', 'inchworm', 'hubbardI', 'ctseg', 'w2dyn_cthyb','hartree'], 'default': 1025},

                                 'n_tau': {'converter': int, 'valid for': lambda x, _: x > 0,
                                           'used': lambda params: params['general']['solver_type'] in ['cthyb', 'ctint', 'inchworm', 'hubbardI', 'ctseg', 'w2dyn_cthyb','hartree'], 'default': 10001},

                                 'n_w': {'converter': int, 'valid for': lambda x, _: x > 0,
                                         'used': lambda params: params['general']['solver_type'] in ['ftps', 'hubbardI', 'hartree'], 'default': 5001},

                                 'w_range': {'converter': lambda s: tuple(map(float, s.split(','))),
                                             'valid for': lambda x, _: x[0] < x[1],
                                             'used': lambda params: params['general']['solver_type'] in ['ftps', 'hubbardI', 'hartree'], 'default': (-10, 10)},

                                 'eta': {'converter': float, 'valid for': lambda x, _: x >= 0,
                                         'used': lambda params: params['general']['solver_type'] in ['ftps', 'hubbardI', 'hartree']},

                                 'diag_delta': {'converter': BOOL_PARSER, 'used': True, 'default': False},

                                 'csc': {'converter': BOOL_PARSER,
                                         'used': True, 'default': False},

                                 'n_iter_dmft_first': {'converter': int, 'valid for': lambda x, _: x > 0,
                                                       'used': lambda params: params['general']['csc'], 'default': 10},

                                 'n_iter_dmft_per': {'converter': int, 'valid for': lambda x, _: x > 0,
                                                     'used': lambda params: params['general']['csc'], 'default': 2},

                                 'plo_cfg': {'used': lambda params: (params['general']['csc']
                                                                     and params['dft']['projector_type'] == 'plo'),
                                             'default': 'plo.cfg'},

                                 'jobname': {'used': lambda params: not params['general']['csc'], 'default': lambda params: 'dmft_dir'},

                                 'h5_save_freq': {'converter': int, 'valid for': lambda x, _: x > 0,
                                                  'used': True, 'default': 5},

                                 'magnetic': {'converter': BOOL_PARSER,
                                              'used': True, 'default': False},

                                 'magmom': {'converter': lambda s: list(map(float, s.split(','))),
                                            'used': lambda params: params['general']['magnetic'],
                                            'default': []},

                                 'h_field': {'converter': float, 'used': True, 'default': 0.0},

                                 'h_field_it': {'converter': int, 'used': True, 'default': 0},

                                 'afm_order': {'converter': BOOL_PARSER,
                                               'used': lambda params: params['general']['magnetic'],
                                               'default': False},

                                 'sigma_mix': {'converter': float,
                                               'valid for': lambda x, params: x >= 0 and (np.isclose(params['general']['g0_mix'], 1)
                                                                                         or np.isclose(x, 1)),
                                               'used': True, 'default': 1.0},

                                 'g0_mix': {'converter': float, 'valid for': lambda x, _: x >= 0,
                                               'used': True, 'default': 1.0},

                                 'g0_mix_type': {'valid for': lambda x, _: x in ('linear', 'broyden'),
                                                'used': True, 'default': 'linear'},

                                 'broy_max_it': {'converter': int, 'valid for': lambda x, _: x >= 1 or x==-1 ,
                                                'used': lambda params: params['general']['g0_mix_type'] == 'broyden', 'default': -1},

                                 'calc_energies': {'converter': BOOL_PARSER, 'used': True, 'default': False},

                                 'block_threshold': {'converter': float, 'valid for': lambda x, _: x > 0,
                                                     'used': True, 'default': 1e-5},

                                 'block_suppress_orbital_symm': {'converter': BOOL_PARSER,
                                                                 'used': lambda params: not params['general']['enforce_off_diag'], 'default': False},

                                 'enforce_off_diag': {'converter': BOOL_PARSER, 'used': True, 'default': False},

                                 'load_sigma': {'converter': BOOL_PARSER, 'used': True, 'default': False},

                                 'path_to_sigma': {'used': lambda params: params['general']['load_sigma']},

                                 'load_sigma_iter': {'converter': int,
                                                     'used': lambda params: params['general']['load_sigma'], 'default': -1},

                                 'noise_level_initial_sigma': {'converter': float,
                                                               'valid for': lambda x, _: x > 0 or np.isclose(x, 0),
                                                               'used': True, 'default': 0.},

                                 # TODO: change default to 'none'
                                 'occ_conv_crit': {'converter': float, 'used': True, 'default': -1},
                                 'gimp_conv_crit': {'converter': float, 'used': True, 'default': -1},
                                 'g0_conv_crit': {'converter': float, 'used': True, 'default': -1},
                                 'sigma_conv_crit': {'converter': float, 'used': True, 'default': -1},

                                 'sampling_iterations': {'converter': int, 'valid for': lambda x, _: x >= 0,
                                                         'used': lambda params: params['general']['occ_conv_crit'] > 0 or
                                                                                params['general']['gimp_conv_crit'] > 0 or
                                                                                params['general']['g0_conv_crit'] > 0,
                                                         'default': 0},

                                 'sampling_h5_save_freq': {'converter': int, 'valid for': lambda x, _: x > 0,
                                                           'used': lambda params: (params['general']['occ_conv_crit'] > 0
                                                                                   and params['general']['sampling_iterations'] > 0),
                                                           'default': 5},

                                 'fixed_mu_value': {'converter': float, 'used': True, 'default': 'none'},

                                 'calc_mu_method': {'valid for': lambda x, _: x in ['dichotomy', 'newton', 'brent'],
                                                 'used': True,
                                                 'default': 'dichotomy',
                                                 },

                                 'mu_update_freq': {'converter': int, 'valid for': lambda x, _: x > 0,
                                                    'used': lambda params: params['general']['fixed_mu_value'] == 'none',
                                                    'default': 1},

                                 'mu_initial_guess': {'converter': float, 'used': True, 'default': 'none'},

                                 'mu_mix_const': {'converter': float,
                                                  'valid for': lambda x, _: x > 0 or np.isclose(x, 0),
                                                  'used': lambda params: params['general']['fixed_mu_value'] == 'none',
                                                  'default': 1.},

                                 'mu_mix_per_occupation_offset': {'converter': float,
                                                                  'valid for': lambda x, _: x > 0 or np.isclose(x, 0),
                                                                  'used': lambda params: params['general']['fixed_mu_value'] == 'none',
                                                                  'default': 0.},

                                 'set_rot': {'valid for': lambda x, _: x in ('none', 'den', 'hloc'),
                                             'used': True, 'default': 'none'},

                                 'measure_chi': {'valid for': lambda x, _: x in ('SzSz', 'NN', 'none'), 'used': True, 'default': 'none'},

                                 'measure_chi_insertions': {'converter': int, 'used': True, 'default': 100},

                                 # TODO: used for which solvers? Generalize to real freq. solvers without maxent?
                                 'mu_gap_gb2_threshold': {'converter': float,
                                                          'valid for': lambda x, _: x == 'none' or x > 0 or np.isclose(x, 0),
                                                          'used': lambda params: params['general']['solver_type'] in ['cthyb', 'ctint', 'ctseg', 'w2dyn_cthyb'],
                                                          'default': 'none'},

                                 'mu_gap_occ_deviation': {'converter': float,
                                                          'valid for': lambda x, _: x == 'none' or x > 0 or np.isclose(x, 0),
                                                          'used': lambda params: (params['general']['solver_type'] in ['cthyb', 'ctint', 'ctseg', 'w2dyn_cthyb']
                                                                                  and params['general']['mu_gap_gb2_threshold'] != 'none'),
                                                          'default': 'none'},

                                 'h_int_basis' : {'valid for': lambda x, _: x in ('triqs', 'wien2k', 'wannier90', 'qe', 'vasp'),
                                            'used': True, 'default' : 'triqs'}

                                },
                     'dft': {'dft_code': {'used': lambda params: params['general']['csc'],
                                          'valid for': lambda x, _: x in ('vasp', 'qe')},

                             'n_cores': {'converter': int, 'valid for': lambda x, _: x > 0,
                                         'used': lambda params: params['general']['csc']},

                             'n_iter': {'converter': int, 'valid for': lambda x, _: x > 0,
                                        'used': lambda params: params['general']['csc'] and params['dft']['dft_code'] == 'vasp',
                                        'default': 4},

                             'n_iter_first': {'converter': int, 'valid for': lambda x, _: x > 0,
                                              'used': lambda params: params['general']['csc'] and params['dft']['dft_code'] == 'vasp',
                                              'default': lambda params: params['dft']['n_iter']},

                             'dft_exec': {'used': lambda params: params['general']['csc'], 'default': 'vasp_std'},

                             'store_eigenvals': {'converter': BOOL_PARSER,
                                                 'used': lambda params: params['general']['csc'],
                                                 'default': False},

                             'mpi_env': {'valid for': lambda x, _: x in ('default', 'openmpi', 'openmpi-intra', 'mpich'),
                                         'used': lambda params: params['general']['csc'], 'default': 'default'},

                             'projector_type': {'valid for': lambda x, params: x == 'w90' or x == 'plo' and params['dft']['dft_code'] == 'vasp',
                                                'used': lambda params: params['general']['csc'], 'default': 'w90'},

                             'w90_exec': {'used': lambda params: (params['general']['csc']
                                                                  and params['dft']['dft_code'] == 'vasp' and params['dft']['projector_type'] == 'w90'),
                                                'default': 'wannier90.x'},

                             'w90_tolerance': {'converter': lambda s: float(s),
                                                'used': lambda params: params['general']['csc'] and params['dft']['projector_type'] == 'w90',
                                                'default': 1e-6},
                            },
                     'solver': {
                                #
                                'store_solver': {'converter': BOOL_PARSER, 'used': True, 'default': False},
                                #
                                # cthyb parameters
                                #
                                'length_cycle': {'converter': int, 'valid for': lambda x, _: x > 0,
                                                 'used': lambda params: params['general']['solver_type'] in ['cthyb', 'ctint', 'ctseg', 'w2dyn_cthyb']},

                                'n_warmup_cycles': {'converter': lambda s: int(float(s)),
                                                    'valid for': lambda x, _: x > 0,
                                                    'used': lambda params: params['general']['solver_type'] in ['cthyb', 'ctint', 'ctseg', 'w2dyn_cthyb']},

                                'n_cycles_tot': {'converter': lambda s: int(float(s)),
                                                 'valid for': lambda x, _: x >= 0,
                                                 'used': lambda params: params['general']['solver_type'] in ['cthyb', 'ctint', 'ctseg', 'w2dyn_cthyb']},

                                'max_time': {'converter': int, 'valid for': lambda x, _: x >= 0,
                                             'default': None,
                                             'used': lambda params: params['general']['solver_type'] in ['cthyb', 'ctint', 'ctseg', 'w2dyn_cthyb']},

                                'imag_threshold': {'converter': float, 'default': None,
                                                   'used': lambda params: params['general']['solver_type'] in ['cthyb']},

                                'off_diag_threshold': {'converter': float, 'default': 0.0,
                                                   'used': lambda params: params['general']['solver_type'] in ['cthyb']},

                                'delta_interface': {'converter': BOOL_PARSER, 'default': False,
                                                  'used': lambda params: params['general']['solver_type'] in ['cthyb', 'wd2yn_cthyb']},

                                'measure_G_tau': {'converter': BOOL_PARSER, 'default': True,
                                                  'used': lambda params: params['general']['solver_type'] in ['hubbardI', 'ctseg', 'w2dyn_cthyb']},

                                'measure_G_iw': {'converter': BOOL_PARSER, 'default': False,
                                                  'used': lambda params: params['general']['solver_type'] in ['ctseg']},

                                'measure_G_l': {'converter': BOOL_PARSER, 'default': False,
                                                'used': lambda params: params['general']['solver_type'] in ['cthyb', 'hubbardI', 'ctseg', 'w2dyn_cthyb']},

                                'measure_density_matrix': {'converter': BOOL_PARSER, 'default': False,
                                                           'used': lambda params: params['general']['solver_type'] in ['cthyb', 'hubbardI']},

                                'move_double': {'converter': BOOL_PARSER, 'default': True,
                                                'used': lambda params: params['general']['solver_type'] in ['cthyb', 'ctint']},

                                'measure_pert_order': {'converter': BOOL_PARSER, 'default': False,
                                                       'used': lambda params: params['general']['solver_type'] in ['cthyb', 'ctint', 'ctseg', 'w2dyn_cthyb']},

                                'move_shift': {'converter': BOOL_PARSER, 'default': False,
                                                'used': lambda params: params['general']['solver_type'] in ['cthyb']},

                                'random_seed': {'converter': str, 'default': None,
                                                'used': lambda params: params['general']['solver_type'] in ['cthyb', 'ctint', 'ctseg', 'w2dyn_cthyb']},

                                'perform_tail_fit': {'converter': BOOL_PARSER,
                                                     'used': lambda params: params['general']['solver_type'] in ['cthyb']
                                                             and not params['solver']['measure_G_l'],
                                                     'default': False},

                                'fit_max_moment': {'converter': int, 'valid for': lambda x, _: x >= 0,
                                                   'used': lambda params: 'perform_tail_fit' in params['solver']
                                                           and params['solver']['perform_tail_fit']
                                                           and params['general']['solver_type'] in ['cthyb'],
                                                   'default': None},

                                'fit_min_n': {'converter': int, 'valid for': lambda x, _: x >= 0,
                                              'used': lambda params: 'perform_tail_fit' in params['solver']
                                                      and params['solver']['perform_tail_fit']
                                                      and params['general']['solver_type'] in ['cthyb'],
                                              'default': None},

                                'fit_max_n': {'converter': int, 'valid for': lambda x, params: x >= params['solver']['fit_min_n'],
                                              'used': lambda params: 'perform_tail_fit' in params['solver']
                                                      and params['solver']['perform_tail_fit']
                                                      and params['general']['solver_type'] in ['cthyb'],
                                              'default': None},

                                'fit_min_w': {'converter': float, 'valid for': lambda x, _: x >= 0,
                                              'used': lambda params: 'perform_tail_fit' in params['solver']
                                                      and params['solver']['perform_tail_fit']
                                                      and params['general']['solver_type'] in ['cthyb'],
                                              'default': None},

                                'fit_max_w': {'converter': float, 'valid for': lambda x, params: x >= params['solver']['fit_min_w'],
                                              'used': lambda params: 'perform_tail_fit' in params['solver']
                                                      and params['solver']['perform_tail_fit']
                                                      and params['general']['solver_type'] in ['cthyb'],
                                              'default': None},

                                'legendre_fit': {'converter': BOOL_PARSER,
                                                 'used': lambda params: params['general']['solver_type'] in ['cthyb', 'ctseg', 'w2dyn_cthyb'],
                                                 'default': False},

                                'loc_n_min': {'converter': int, 'valid for': lambda x, _: x >= 0,
                                              'used': lambda params: params['general']['solver_type'] in ['cthyb'],
                                              'default': None},

                                'loc_n_max': {'converter': int, 'valid for': lambda x, _: x >= 0,
                                              'used': lambda params: params['general']['solver_type'] in ['cthyb'],
                                              'default': None},


                                #
                                # extra ctseg params
                                #
                                'improved_estimator': {'converter': BOOL_PARSER,
                                                 'used': lambda params: params['general']['solver_type'] in ['ctseg'],
                                                 'default': False},

                                #
                                # extra hartree params
                                #
                                'with_fock': {'converter': BOOL_PARSER,
                                                 'used': lambda params: params['general']['solver_type'] in ['hartree'],
                                                 'default': False},
                                'one_shot': {'converter': BOOL_PARSER,
                                                 'used': lambda params: params['general']['solver_type'] in ['hartree'],
                                                 'default': False},
                                'force_real': {'converter': BOOL_PARSER,
                                                 'used': lambda params: params['general']['solver_type'] in ['hartree'],
                                                 'default': True},
                                'method': {'valid for': lambda x, _: x in ['krylov', 'broyden1', 'broyden2', 'hybr', 'linearmixing'],
                                                 'used': lambda params: params['general']['solver_type'] in ['hartree'],
                                                 'default': 'krylov'},
                                'tol': {'converter': float, 'valid for': lambda x, _: x >= 0,
                                              'used': lambda params: params['general']['solver_type'] in ['hartree'],
                                              'default': 1e-5},

                                #
                                # ftps parameters
                                #
                                'n_bath': {'converter': int, 'valid for': lambda x, _: x >= 0,
                                           'used': lambda params: params['general']['solver_type'] in ['ftps'],
                                           'default': 0},
                                'bath_fit': {'converter': BOOL_PARSER,
                                             'used': lambda params: params['general']['solver_type'] in ['ftps']},
                                'refine_factor': {'converter': int, 'valid for': lambda x, _: x >= 0,
                                                  'used': lambda params: params['general']['solver_type'] in ['ftps'],
                                                  'default': 1},
                                'ph_symm': {'converter': BOOL_PARSER, 'default': False,
                                            'used': lambda params: params['general']['solver_type'] in ['ftps']},
                                'calc_me': {'converter': BOOL_PARSER, 'default': True,
                                            'used': lambda params: params['general']['solver_type'] in ['ftps']},
                                'enforce_gap': {'converter': lambda s: [float(eval(x)) for x in s.split(',')],
                                                'used': lambda params: params['general']['solver_type'] in ['ftps'],
                                                'default': 'none'},
                                'ignore_weight': {'converter': float, 'valid for': lambda x, _: x >= 0,
                                                  'used': lambda params: params['general']['solver_type'] in ['ftps'],
                                                  'default': 0.0},
                                'dt': {'converter': float, 'valid for': lambda x, _: x >= 0,
                                       'used': lambda params: params['general']['solver_type'] in ['ftps']},
                                'state_storage': {'used': lambda params: params['general']['solver_type'] in ['ftps'],
                                                  'default': './'},
                                'path_to_gs': {'used': lambda params: params['general']['solver_type'] in ['ftps'],
                                               'default': 'none'},
                                'sweeps': {'converter': int, 'valid for': lambda x, _: x >= 0,
                                           'used': lambda params: params['general']['solver_type'] in ['ftps'],
                                           'default': 10},
                                'maxm': {'converter': int, 'valid for': lambda x, _: x > 0, 'default': 100,
                                         'used': lambda params: params['general']['solver_type'] in ['ftps']},
                                'maxmI': {'converter': int, 'valid for': lambda x, _: x > 0, 'default': 100,
                                          'used': lambda params: params['general']['solver_type'] in ['ftps']},
                                'maxmIB': {'converter': int, 'valid for': lambda x, _: x > 0, 'default': 100,
                                           'used': lambda params: params['general']['solver_type'] in ['ftps']},
                                'maxmB': {'converter': int, 'valid for': lambda x, _: x > 0, 'default': 100,
                                          'used': lambda params: params['general']['solver_type'] in ['ftps']},
                                'tw': {'converter': float, 'valid for': lambda x, _: x > 0, 'default': 1e-9,
                                       'used': lambda params: params['general']['solver_type'] in ['ftps']},
                                'dmrg_maxm': {'converter': int, 'valid for': lambda x, _: x > 0, 'default': 100,
                                              'used': lambda params: params['general']['solver_type'] in ['ftps']},
                                'dmrg_maxmI': {'converter': int, 'valid for': lambda x, _: x > 0, 'default': 100,
                                               'used': lambda params: params['general']['solver_type'] in ['ftps']},
                                'dmrg_maxmIB': {'converter': int, 'valid for': lambda x, _: x > 0, 'default': 100,
                                                'used': lambda params: params['general']['solver_type'] in ['ftps']},
                                'dmrg_maxmB': {'converter': int, 'valid for': lambda x, _: x > 0, 'default': 100,
                                               'used': lambda params: params['general']['solver_type'] in ['ftps']},
                                'dmrg_tw': {'converter': float, 'valid for': lambda x, _: x > 0, 'default': 1e-9,
                                            'used': lambda params: params['general']['solver_type'] in ['ftps']},
                                #
                                # w2dyn_cthyb parameters
                                #
                                'measure_G2_iw_ph': {'converter': BOOL_PARSER, 'default': False,
                                                     'used': lambda params: params['general']['solver_type'] in ['w2dyn_cthyb']},
                                'measure_G2_n_bosonic': {'converter': int, 'valid for': lambda x, _: x > 0, 'default': 30,
                                                     'used': lambda params: params['general']['solver_type'] in ['w2dyn_cthyb']},
                                'measure_G2_n_fermionic': {'converter': int, 'valid for': lambda x, _: x > 0, 'default': 30,
                                                     'used': lambda params: params['general']['solver_type'] in ['w2dyn_cthyb']},
                                'worm_components': {'converter': lambda s: list(map(int, s.split(','))), 'default': None,
                                                     'used': lambda params: params['general']['solver_type'] in ['w2dyn_cthyb']},
                                'move_global_prob': {'converter': float, 'default': 0.005,
                                                     'used': lambda params: params['general']['solver_type'] in ['w2dyn_cthyb']},
                                'flavourchange_moves': {'converter': int, 'default': 0,
                                                        'used': lambda params: params['general']['solver_type'] in ['w2dyn_cthyb']},
                                'statesampling': {'converter': int, 'default': 0,
                                                  'used': lambda params: params['general']['solver_type'] in ['w2dyn_cthyb']},

                               },
                     'advanced': {'dc_factor': {'converter': float, 'used': True, 'default': 'none'},

                                  'dc_fixed_value': {'converter': float, 'used': True, 'default': 'none'},

                                  'dc_fixed_occ': {'converter': lambda s: list(map(float, s.split(','))),
                                                   'used': True, 'default': 'none'},

                                  'dc_nominal': {'converter': BOOL_PARSER, 'used': True, 'default': False},

                                  'dc_U': {'converter': lambda s: list(map(float, s.split(','))),
                                           'used': True, 'default': lambda params: params['general']['U']},

                                  'dc_J': {'converter': lambda s: list(map(float, s.split(','))),
                                           'used': True, 'default': lambda params: params['general']['J']},

                                  'map_solver_struct': {'converter': eval,
                                                        'valid for': lambda x, _: x =='none' or isinstance(x, dict) or isinstance(x, list),
                                                        'used': True, 'default': 'none'},

                                  'mapped_solver_struct_degeneracies': {'converter': eval,
                                            'valid for': lambda x, _: x=='none' or isinstance(x, list),
                                            'used': lambda params: params['advanced']['map_solver_struct'] != 'none',
                                            'default': 'none'},
                                  'pick_solver_struct': {'converter': eval,
                                                        'valid for': lambda x, _: x =='none' or isinstance(x, dict) or isinstance(x, list),
                                                        'used': True, 'default': 'none'},
                                 }
                    }


# -------------------------- config section cleanup --------------------------
def _config_find_default_section_entries(config):
    """
    Returns all items in the default section.

    Parameters
    ----------
    config : ConfigParser
        A configparser instance that has read the config file already.

    Returns
    -------
    list
        All entries in the default section.
    """
    return list(config['DEFAULT'].keys())

def _config_add_empty_sections(config):
    """
    Adds empty sections if no parameters in the whole section were given.

    Parameters
    ----------
    config : ConfigParser
        A configparser instance that has read the config file already.

    Returns
    -------
    config : ConfigParser
        The config parser with all required sections.
    """
    for section_name in PROPERTIES_PARAMS:
        if section_name not in config:
            config.add_section(section_name)

    return config

def _config_remove_unused_sections(config):
    """
    Removes sections that are not supported by this program.

    Parameters
    ----------
    config : ConfigParser
        A configparser instance that has read the config file already.

    Returns
    -------
    config : ConfigParser
        The config parser without unnecessary sections.
    unused_sections : list
        All sections that are not supported.
    """
    unused_sections = []
    for section_name in list(config.keys()):
        if section_name != 'DEFAULT' and section_name not in PROPERTIES_PARAMS.keys():
            unused_sections.append(section_name)
            config.remove_section(section_name)

    return config, unused_sections

# -------------------------- parameter reading --------------------------
def _convert_params(config):
    """
    Applies the converter given in the PROPERTIES_PARAMS to the config. If no
    converter is given, a default string conversion is used.

    Parameters
    ----------
    config : ConfigParser
        A configparser instance that has passed through the above clean-up
        methods.

    Returns
    -------
    parameters : dict
        Contains dicts for each section. These dicts contain all parameter
        names and their respective value that are in the configparser and in
        the PROPERTIES_PARAMS.
    """
    parameters = {name: {} for name in PROPERTIES_PARAMS}

    for section_name, section_params in parameters.items():
        for param_name, param_props in PROPERTIES_PARAMS[section_name].items():
            if param_name not in config[section_name]:
                continue

            # Uses converter for parameters
            if 'converter' in param_props:
                section_params[param_name] = param_props['converter'](str(config[section_name][param_name]))
            else:
                section_params[param_name] = str(config[section_name][param_name])

    return parameters


def _find_nonexistent_params(config):
    """
    Returns all parameters that are in the config but not in the
    PROPERTIES_PARAMS and are therefore not recognized by the program.

    Parameters
    ----------
    config : ConfigParser
        A configparser instance that has passed through the above clean-up
        methods.

    Returns
    -------
    nonexistent_params : dict
        Contains a list for each section, which contains all unused parameter
        names.
    """
    nonexistent_params = {section_name: [] for section_name in PROPERTIES_PARAMS}

    for section_name, section_params in PROPERTIES_PARAMS.items():
        for param_name in config[section_name]:
            if param_name not in (key.lower() for key in section_params):
                nonexistent_params[section_name].append(param_name)

    return nonexistent_params


def _apply_default_values(parameters):
    """
    Applies default values to all parameters that were not given in the config.

    Parameters
    ----------
    parameters : dict
        Contains dicts for each section, which contain the parameter names and
        their values as read from the config file.

    Returns
    -------
    parameters : dict
        The parameters dict including the default values.
    default_values_used : dict
        Contains a list for each section, which contains all parameters that
        were set to their default values. Used to find out later which
        unnecessary parameters were actually given in the config file.
    """
    default_values_used = {section_name: [] for section_name in PROPERTIES_PARAMS}

    for section_name, section_params in PROPERTIES_PARAMS.items():
        for param_name, param_props in section_params.items():
            if 'default' not in param_props:
                continue

            if param_name in parameters[section_name]:
                continue

            default_values_used[section_name].append(param_name)
            if callable(param_props['default']):
                parameters[section_name][param_name] = param_props['default'](parameters)
            else:
                parameters[section_name][param_name] = param_props['default']

    return parameters, default_values_used


def _check_if_params_used(parameters, default_values_used):
    """
    Checks if the parameters in the config file are used or unnecessary.

    Parameters
    ----------
    parameters : dict
        Contains dicts for each section, which contain the parameter names and
        their values as read from the config file or otherwise set to default.
    default_values_used : dict
        Contains a list for each section, which contains all parameters that
        were set to their default values.

    Returns
    -------
    parameters : dict
        The parameters dict where all unnecessary parameters were removed.
    unnecessary_params : dict
        Contains a list for each section, which contains all parameters that
        were given in the config file but are unnecessary.
    missing_required_params : dict
        Contains a list for each section, which contains all parameters that
        are required, have no default and are missing from the config file.
    """
    unnecessary_params = {section_name: [] for section_name in PROPERTIES_PARAMS}
    missing_required_params = {section_name: [] for section_name in PROPERTIES_PARAMS}

    for section_name, section_params in PROPERTIES_PARAMS.items():
        for param_name, param_props in section_params.items():
            # 'used' could be bool or function returning bool
            if callable(param_props['used']):
                required = param_props['used'](parameters)
            else:
                required = param_props['used']

            if required:
                if param_name not in parameters[section_name]:
                    missing_required_params[section_name].append(param_name)
                elif parameters[section_name][param_name] is None:
                    del parameters[section_name][param_name]
                continue

            if param_name in parameters[section_name]:
                del parameters[section_name][param_name]

                if param_name not in default_values_used[section_name]:
                    unnecessary_params[section_name].append(param_name)

    return parameters, unnecessary_params, missing_required_params

def _checks_validity_criterion(parameters):
    """
    Checks the validity criterion from the PROPERTIES_PARAMS.

    Parameters
    ----------
    parameters : dict
        Contains dicts for each section, which contain the parameter names and
        their values as read from the config file, set to default if required
        or removed if they are unnecessary.

    Returns
    -------
    invalid_params : dict
        Contains a list for each section, which contains all parameters that
        do not fulfill their validity criterion.
    """
    invalid_params = {section_name: [] for section_name in PROPERTIES_PARAMS}

    for section_name, section_params in PROPERTIES_PARAMS.items():
        for name, value in parameters[section_name].items():
            if 'valid for' not in section_params[name]:
                continue

            condition = section_params[name]['valid for']
            if not condition(value, parameters):
                invalid_params[section_name].append(name)

    return invalid_params


def read_config(config_file):
    """
    Reads in the config file, checks its sections and parameters and returns
    the parameters sorted by their categories.

    Parameters
    ----------
    config_file : string
        File name of the config file usable for configparser.

    Raises
    ------
    ValueError
        Required parameters are missing or parameters do not fulfill their
        validity criterion.

    Returns
    -------
    general_params : dict

    solver_params : dict

    dft_params : dict

    advanced_params : dict
    """
    config = ConfigParser()
    config.read(config_file)

    # Checks if default section is empty
    config_default_entries = _config_find_default_section_entries(config)
    if config_default_entries:
        print('Warning: the following parameters are not in any section and will be ignored:')
        print(', '.join(config_default_entries))

    # Adds empty sections if they don't exist
    config = _config_add_empty_sections(config)

    # Removes unused sections and prints a warning
    config, unused_sections = _config_remove_unused_sections(config)
    if unused_sections:
        print('Warning: ignoring parameters in following unexpected sections:')
        print(', '.join(unused_sections))

    # Reads and converts all valid parameters
    parameters = _convert_params(config)

    # Checks for unused parameters given in the config file
    nonexistent_params = _find_nonexistent_params(config)
    if any(nonexistent_params.values()):
        print('Warning: the following parameters are not supported:')
        for section_name, section_params in nonexistent_params.items():
            if section_params:
                print('- Section "{}": '.format(section_name) + ', '.join(section_params))

    # Applies default values
    parameters, default_values_used = _apply_default_values(parameters)

    # Prints warning if unnecessary parameters are given
    parameters, unnecessary_params, missing_required_params = _check_if_params_used(parameters, default_values_used)
    if any(unnecessary_params.values()):
        print('Warning: the following parameters are given but not used in this calculation:')
        for section_name, section_params in unnecessary_params.items():
            if section_params:
                print('- Section "{}": '.format(section_name) + ', '.join(section_params))

    # Raises error if required parameters are not given
    if any(missing_required_params.values()):
        required_error_string = ''
        for section_name, section_params in missing_required_params.items():
            if section_params:
                required_error_string += '\n- Section "{}": '.format(section_name) + ', '.join(section_params)
        raise ValueError('The following parameters are required but not given:'
                         + required_error_string)

    # Raises error if parameters invalid
    invalid_params = _checks_validity_criterion(parameters)
    if any(invalid_params.values()):
        invalid_error_string = ''
        for section_name, section_params in invalid_params.items():
            if section_params:
                invalid_error_string += '\n- Section "{}": '.format(section_name) + ', '.join(section_params)
        raise ValueError('The following parameters are not valid:'
                         + invalid_error_string)

    # warning if sigma mixing is used, remove in future versions
    if parameters['general']['sigma_mix'] < 1.0:
        if parameters['general']['g0_mix'] < 1.0:
            raise ValueError('You shall not use Sigma and G0 mixing together!')

    # Workarounds for some parameters

    # make sure that pick_solver_struct and map_solver_struct are a list of dict
    if isinstance(parameters['advanced']['pick_solver_struct'], dict):
        parameters['advanced']['pick_solver_struct'] = [parameters['advanced']['pick_solver_struct']]

    if isinstance(parameters['advanced']['map_solver_struct'], dict):
        parameters['advanced']['map_solver_struct'] = [parameters['advanced']['map_solver_struct']]

    if parameters['general']['solver_type'] in ['cthyb', 'ctint', 'ctseg', 'w2dyn_cthyb']:
        parameters['solver']['n_cycles'] = parameters['solver']['n_cycles_tot'] // mpi.size
        del parameters['solver']['n_cycles_tot']

    if parameters['general']['solver_type'] in ['cthyb']:
        parameters['general']['cthyb_delta_interface'] = parameters['solver']['delta_interface']
        del parameters['solver']['delta_interface']

    if parameters['general']['solver_type'] in ['ctseg']:
        # some parameters have different names for ctseg
        parameters['solver']['measure_gt'] = parameters['solver']['measure_G_tau']
        del parameters['solver']['measure_G_tau']

        parameters['solver']['measure_gw'] = parameters['solver']['measure_G_iw']
        del parameters['solver']['measure_G_iw']

        # make sure measure_gw is true if improved estimators are used
        if parameters['solver']['improved_estimator']:
            parameters['solver']['measure_gt'] = True
            parameters['solver']['measure_ft'] = True
        else:
            parameters['solver']['measure_ft'] = False
        del parameters['solver']['improved_estimator']

        parameters['solver']['measure_gl'] = parameters['solver']['measure_G_l']
        del parameters['solver']['measure_G_l']

        parameters['solver']['measure_hist'] = parameters['solver']['measure_pert_order']
        del parameters['solver']['measure_pert_order']

    if parameters['general']['solver_type'] in ['cthyb'] and parameters['solver']['measure_density_matrix']:
        # also required to measure the density matrix
        parameters['solver']['use_norm_as_weight'] = True

    if parameters['general']['solver_type'] in ['ftps'] and parameters['general']['calc_energies']:
        raise ValueError('"calc_energies" is not valid for solver_type = "ftps"')

    if parameters['general']['dc'] and parameters['general']['dc_type'] == 4:
        assert sum(parameters['general']['cpa_x']) == 1., 'Probability distribution for CPA must equal 1.'

    # little workaround since #leg coefficients is not directly a solver parameter
    if 'legendre_fit' in parameters['solver']:
        parameters['general']['legendre_fit'] = parameters['solver']['legendre_fit']
        del parameters['solver']['legendre_fit']

    parameters['general']['store_solver'] = parameters['solver']['store_solver']
    del parameters['solver']['store_solver']

    return parameters['general'], parameters['solver'], parameters['dft'], parameters['advanced']
