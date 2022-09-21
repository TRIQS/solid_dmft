(changelog)=

# Changelog

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

