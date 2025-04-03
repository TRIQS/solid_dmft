(changelog)=

# Changelog

## Version 3.3.2

solid_dmft version 3.3.2 is a minor patch release that fixes a few small bugs and improves documentation:

### doc
* update Ce2O3 tutorials input file

### fix
* wrong error message for tail treatment in cthyb
* pcb projection on orb feat for Fermi surfaces in TB only mode was not giving correct results
* cthyb measure_chi was not passed into solve params since 3.3.x
* small bug in SOC calculation

We thank all contributors: Thomas Hahn, Nils Wentzell, Ina Park, Alexander Hampel

## Version 3.3.1

solid_dmft version 3.3.1 is a minor release that fixes bugs in the new triqs/ctseg interface and improves the interface to GW:

### fix
* small bug in dc initialization
* dynamic U tensor corrections
* initial self-energy from dc for magnetic calculation
* align ctseg interface with triqs/ctseg release
* ctseg hartree_shift rename in triqs

### feat
* add print of avg order for ctseg
* use symmetric meshes for DLR
* fix delta_interface in ctseg

### build
* add missing default.toml to pypi packaging

We thank all contributors: Alexander Hampel, Nils Wentzell, Chia-Nan Yeh

## Version 3.3.0

solid_dmft 3.3.0 is a major release, compatible with TRIQS 3.3, updated to the latest app4triqs skeleton, and bringing major changes to the code:

* the input parser is switched to a general toml parser, i.e. strings have to be passed in quotes, boolean parameters given without capitalization, and lists passed with brackets. Check below separate section for detailed changes.
* the new input parser allows to define for each impurity problem a different solver if wanted, i.e. d-shell with cthyb and p-shell with Hartree-Fock. See new text NiO-cthyb
* docker images are automatically build on each push for all major releases to ghcr.io
* switch from old ctseg to new ctseg_j solver
* allow CRM Dyson solver for both cthyb and ctseg to obtain Sigma_imp
  from G_tau:   "crm_dyson_solver=true" and dlr_wmax and dlr_eps (see https://triqs.github.io/triqs/unstable/documentation/python_api/triqs.gf.dlr_crm_dyson_solver.html#module-triqs.gf.dlr_crm_dyson_solver for details)
* add new DC schemes 'crpa_static', 'crpa_static_qp', 'crpa_dynamic'
* use cRPA calculated Uijkl as interaction via 'crpa',
  'crpa_density_density', 'dyn_density_density', 'dyn_full' hint types
* read interaction tensor from AIMBES h5
* new experimental gw_embedding module. See gw_embedding/gw_flow.py for details allowing to run solid_dmft on top of AIMBES
* allow to use Pade for AC in post-processing

We thank all contributors: Sophie Beck, Thomas Hahn, Alexander Hampel, Henri Menke, Maximilian Merkel, and Nils Wentzell

Find below an itemized list of changes in this release.

### General
* merge dev GW embedding (includes other fixes as well) (#78)
* pass gw params to all methods
* multiple solvers and toml input parser (#74)
* added toml to docker images
* Restore Python 3.8 compatibility for dictionary merge (#63)
* Allow mathematical expression to be passed for random_seed (#61)
* allow PCB to read from TRIQS TB object
* respack fit slater for p shell
* add Pade Sigma analytic continuation and refine tests
* add simple_intra interaction, for intro orbital only interaction
* add dc_orb_shift param to allow orbital dependent shift in impurity levels
* allow 0.0 mixing to perform stat sampling
* switch all pytest to unit tests

### new toml input parser
* The following input parameters can now be a list per impurity:
    * `general_params`: U, J, U_prime, ratio_F4_F2, h_int_type, enforce_off_diag, dc_type
    * `advanced_params`: dc_U, dc_J, dc_fixed_occ, map_solver_struct, pick_solver_struct, mapped_solver_struct_degeneracies
* Multiple solvers can be used, which only solve the impurity problems specified in `idx_impurities`
    * general parameter `solver_type` moved to solver section and renamed to `type`
    * general parameter `n_l` moved to solver section
    * general parameter `measure_chi` moved to solver section
    * general parameter `delta_interface` moved to solver section
* All possible input parameters are defined in the `python/solid_dmft/io_tools/default.toml`
* according to toml format the config file is now called .toml (instead of .ini), and boolean are not capitalized, strings are given with quotes and lists are given with brackets.
* Documentation of the input is now generated from `python/solid_dmft/io_tools/documentation.txt`
* For an example, refer to the new integration test (see below)
* Updated interface to python scripts wrapping solid_dmft: new routine `main.run_dmft` that expects the params as python dictionaries, which are then supplemented with the defaults etc equivalent to what happens when reading in a toml file

* the existence of the parameter `general_params['beta']` now determines if a imaginary- or real-frequency grid is used within solid_dmft
* Bug fix: Slater interaction for p orbitals can now be constructed
* Renaming of solver parameters for the different solvers is now moved to `solver.py`. The idea is that every other part of solid_dmft should care as little as possible what solvers are used, with the details abstracted by the SolverStructure class
    * In `solver.py`, all solver parameters that are passed to the triqs solver are transferred to a dict `triqs_solver_params`. When adding new triqs solver parameter to solid_dmft in the future, they also need to be added within solver.py.
* In the determination of the block structure, the largely unused parameter `general_params['block_suppress_orbital_symm']` removed. Its behavior can be replaced by using `advanced_params['mapped_solver_struct_degeneracies']`
* Integration tests: previously existing tests updated, new tests added. One with ftps solver (requires installation of ftps, otherwise just passes without doing anything) and one with a combination of CT-HYB and Hartree solver
* Unit tests: added test for toml-related functionality
* `read_config.py` removed and the functionality for dealing with the dicts from reading a toml file moved to `postproc_toml_dict.py`
* `io_tools/verify_input_params.py` contains all checks of the input params that the code performs before starting the DMFT calculations
* Updated the documentation of the input parameters

### doc
* add comment that proj in postprocessing is only correct for diag A(k,w)
* update NNO magnetic tutorial
* fix Vasp CSC tutorial for PNO after CSC fixes

### build
* add new tests for CRM Dyson solver (requires triqs 3.3)
* add new GW embedding tests that run optionally with -DTest_GW_embedding=ON
* modify basic SVO test to do crm test instead of gl
* add useful apt packages to openmpi image
* use ghcr.io images when testing PR
* ci: build and cache base image separately (#70)
* use new auto build ghcr.io docker images
* add GitHub Actions workflow for Docker images (#66)
* simplify dockerfile for github ci
* trigger pypi build on tags
* add pypi workflow
* update Vasp patches for ver 6.4
* Cleaned up VASP diff files for CSC
* use cmake variable to determine max number of mpi ranks during testing

### other fixes
* Added warning for matrix-valued selfenergy continuation
* draw colorbar only once in kslices
* PCB bug aprx Sigma as diagonal if interpolation is used
* broken FS: np.shape -> len
* fix small FTPS problems and introduce a different eta for FTPS
* maxent test precision fix and test dependency
* use of origin in Fermi surface
* fix calculation of Akw for off-diag Sigma


## Version 3.2.3
solid_dmft version 3.2.3 and 3.2.2 are minor releases that fixes bugs in the post-processing routines and brings small new improvements:

* allow 0.0 mixing to perform stat sampling
* allow mathematical expression to be passed for random_seed
* fix broken FS plot in PCB: np.shape -> len
* fix PCB bug aprx Sigma as diagonal if interpolation is used
* fix PCB to draw coloarbar only once in kslices

We thank all contributors: Alexander Hampel, Henri Menke

## Version 3.2.1
solid_dmft version 3.2.1 is a minor release that automatizes the pypi packaging release

## Version 3.2.0

solid_dmft version 3.2.0 is a release that
* adds Jenkins CI support via flatiron-jenkins
* includes several fixes to match the latest triqs 3.2.x release
* changes the Z estimate to a correct linear fit of the first two Matsubara frequencies
* fixes for QE and Vasp CSC
* add option to add a magnetic field in DMFT
* add solid_dmft JOSS paper reference (doi.org/10.21105/joss.04623)
* add simple Ntot interaction
* allow Uprime!=U-2J in Kanamori
* updates the tutorials
* introduces input output documentation
* add support for the TRIQS Hartree Solver
* add RESPACK support

We thank all contributors: Sophie Beck, Alberto Carta, Alexander Hampel, Max Merkel, Harrison LaBollita, Nils Wentzell

Find below an itemized list of changes in this release.

General
-------
* fix SzSz measurement in triqs unstable
* Updated mpich VASP5 docker file to include HF solver
* add hartree solver
* feat: add regular kmesh option to pcb postproc
* Fix to charge-self-consistency with Vasp (#48)
* removed QE fix files which are now in official release
* Modified dockerfile to add pmi support for cray supercomputing environments
* add RESPACK postprocessing routines (#38)
* Added correction to energy calculation
* add triqs logos to skeleton and include ico in install directive of doc
* change name of dft_mu to mu_initial_guess
* support different DFT cubic basis conventions (#36)
* allow magnetic calculation for CSC (output den correction is always averaged)
* fix sym bug in hubbardI postprocessing
* always calculate dft_mu at start of calculation
* add h_field_it to remove magnetic field after x iterations
* Write solid_dmft hash to h5
* fix delta interface of cthyb for multiple sites with different block structures
* correctly use tail fitted Sigma from cthyb not via double dyson equation
* add paper ref to toml
* minor addition of post-processing script: add_local Hamiltonian, separate from add_lambda. We might remove add_lambda
* update doc with JOSS references
* Bug fix for changes in sumk mesh definition in maxent_gf_latt
* adapt vasp patch files for ver6.3.2
* function to det n_orb_solver, fix test
* apply block picker before block mapping
* fix header writing for obs file
* add pick solver struct option to select specific blocks for the impurity problem
* fix print for failing comparison test
* allow different interaction Hamiltonians per impurity
* enforce PEP standard in interaction Hamiltonian
* print optimal alpha in other maxent scripts
* final corrections for PCB functions
* add proj_on_orb functionality to Akw
* fix bug in max_G_diff function ignoring norm_temp
* change Sigma_imp_iw / _w to Sigma_imp (DFTTools unstable)
* fix load Sigma with new gf_struct in triqs 3.1.x
* adapt to sumk mesh changes in dfttools
* Made the way mesh is stored in maxent_gf_latt consistent with maxent_gf_imp

fix
---
* fix deg shells in magnetic calculations
* fix parameter n_orb in hint construction
* doc strings of cRPA avering for Slater
* critical bug in hubbardI interface
* PCB fermi surface plot
* updates from triqs unstable
* simple Z estimate as linear fit
* PCB: removing "linearize" function, changing the model
* delta_interface with SOC and store solver options
* convert warmup cycles to int automatically
* problem with ish vs icrsh in PCB Thanks @HenryScottx for reporting!
* h_int uses now n_orb instead of orb_names

build
-----
* adapt jenkins CI files
* simplify docker image
* update openmpi docker file with clang-15
* update CI dockerfile
* Updated docker file to ubuntu 22

feat
----
* enable MPI for maxent_gf_imp post-processing routines
* add possibility to specify Uprime in Kanamori interaction
* add loc_n_min / max arg for cthyb
* add additional support for hartree when computing DC from the solver
* add Ntot interaction

doc
---
* Added observables documentation for DMFT output
* Updated tutorial svo one-shot

test
----
* fix tests after Hartree additions
* add Hartree Solver test
* Integration test for maxent gf imp and latt, bug fixes to both scripts (#30)
* add new test for pcb get_dmft_bands function


## Version 3.1.5

solid_dmft version 3.1.5 is a patch-release that improves / fixes the following issues:

* fix to charge-self-consistency with Vasp and QE
* feat add loc_n_min / max arg for cthyb
* fix simple Z estimate as linear fit
* adapt docker images for ubuntu 22.04

Contributors: Sophie Beck, Alberto Carta, Alexander Hampel, Max Merkel:

## Version 3.1.4

solid_dmft version 3.1.4 is a patch-release that improves / fixes the following issues:

* fix and improve rootfinder in PCB for quasiparticle dispersion
* fix pypi package version.py module

Contributors: Sophie Beck, Alberto Carta, Alexander Hampel, Max Merkel:

## Version 3.1.3

solid_dmft version 3.1.3 is a patch-release that improves / fixes the following issues:

* fix delta interface of cthyb for multiple sites with different block structures
* correctly use tail fitted Sigma from cthyb not via double dyson equation
* magnetic param not available in CSC crash PM calc
* improve PCB script from unstable branch
* convert warmup cycles to int automatically
* fix function calls in gap finder
* fix delta_interface with SOC and store solver options
* fix: update svo example for PCB test from unstable

Contributors: Sophie Beck, Alberto Carta, Alexander Hampel, Max Merkel

## Version 3.1.2

solid_dmft version 3.1.1 is a patch-release that improves / fixes the following issues:

* fix deg shells in magnetic calculations
* fix bug in max_G_diff function ignoring norm_temp
* fix load Sigma with new gf_struct in triqs 3.1.x
* Made the way mesh is stored in maxent_gf_latt consistent with maxent_gf_imp
* adapt vasp patch files for ver6.3.2
* update README.md for Joss publication
* print optimal alpha in other maxent scripts
* update postprocessing routines for plotting spectral functions
* add new test for pcb get_dmft_bands function
* DOC: extend install instructions & improve readme for #21 #22
* DOC: update support & contribute section, bump ver to 3.1.1
* add proj_on_orb functionality to Akw
* Added observables documentation for DMFT output
* Added input documentation
* Added ETH logo to website, small fixes to documentation
* rename examples to debbuging_examples
* pip package build files

Contributors: Sophie Beck, Alberto Carta, Alexander Hampel, Max Merkel


## Version 3.1.1

solid_dmft version 3.1.1 is a patch-release that improves / fixes the following issues:

* delete obsolete make_spaghetti.py
* SOC self energies can be continued in maxent
* run hubbardI solver on all nodes due to slow bcast performance of atomdiag object
* fix DFT energy read when running CSC QE
* updated documentation, small fixes to tutorials
* exposed params of maxent_gf_imp
* fix the way dft_mu is loaded in PCB
* fix executable in SVO tutorial
* fix shift in sigma continuator to remove dft_mu
* fix chemical potential in plot Akw and minor fixes
* correct plotlabels in postprocessing
* tiny modification of printing H_loc in postprocessing

Contributors: Sophie Beck, Alberto Carta, Max Merkel

## Version 3.1.0

solid_dmft version 3.1.0 is a major release that provides tutorials in the documentation, changes to app4triqs skeleton, allows CSC calculations with QE, improves postprocessing routines, and add functionality for SOC calculations.

* all new tutorials
* generalize measure_chi functionality
* CSC with Vasp 6.3.0 works, examples updated
* fix two bugs in w90 interface in vasp
* Renamed files
* fix Fermi level print in mlwf.F LPRJ_WRITE call
* Automatic patching of vasp 6.3.0 with Docker
* Updated tutorial
* Added check on all mpi ranks if dmft_config exists at beginning of run
* fix small bug in convergence.py thanks @merkelm
* Rework convergence metrics
* remove gf_struct_flatten from solver in accordance with latest dfttools version
* Renaming to solid_dmft
* Update of maxent_gf_latt.py: more parameters exposed and spin averaging is not default anymore
* fix bug in afm calculation when measuring density matrix
* Add w90_tolerance flag for CSC
* use sphinx autosummary for module reference
* small changes in IO, additional mpi barriers in csc flow for better stability
* With SOC now program prints real and imag part of matrices
* Fixed creation of Kanamori Hamiltonian with SOC
* Improvements in plot_correlated_bands.py and updated tutorial
* change output name of MaxEnt Sigma to Sigma_maxent
* change to develop version of w90 because of mpi bug in openmpi dockerfile
* bugfix in plot_correlated_bands and cleaning up
* update OpenMPI Dockerfile to latest Ubuntu
* Tutorial to explore correlated bands using the postprocessing script
* check in CSC with QE if optional files are presesnt, otherwise skip calculation
* Updated maxent_sigma: mpi parallelization, continuator types, bug fixes, parameters exposed
* update installation instructions
* add workflow and code structure images
* Updated maxent sigma script
* W90 runs in parallel
* Fixing a bug related to measure_pert_order and measure_chi_SzSz for afm_order
* add vasp crpa scripts and tutorials
* add delta interface for cthyb
* fix get_dmft_bands and pass eta to alatt_k_w correctly
* allows to recompute rotation matrix even if W90 is used
* bugfix in initial_self_energies.py in case dc = False
* flatten gf_struct for triqs solvers to remove depracted warning
* add example files for SVO and LNO
* bump triqs and package version to 3.1

Contributors: Sophie Beck, Alberto Carta, Max Merkel

## Version 3.0.0

solid_dmft version 3.0.0 is a compatibility
release for TRIQS version 3.0.0 that
* introduces compatibility with Python 3 (Python 2 no longer supported)
* adds a cmake-based dependency management
* fixes several application issues

