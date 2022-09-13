---
title: 'solid\_dmft: gray-boxing DFT+DMFT materials simulations with TRIQS'
tags:
  - Python
  - electronic structure theory
  - solid-state
  - computational materials science
  - correlated electrons
  - dynamical mean-field theory
authors:
  - name: Maximilian E. Merkel
    orcid: 0000-0003-2589-9625
    affiliation: 1
  - name: Alberto Carta
    orcid: 0000-0003-0705-0281
    affiliation: 1
  - name: Sophie Beck
    orcid: 0000-0002-9336-6065
    affiliation: 2
  - name: Alexander Hampel^[Corresponding author]
    orcid: 0000-0003-1041-8614
    affiliation: 2
affiliations:
 - name: Materials Theory, ETH Zürich, Wolfgang-Pauli-Strasse 27, 8093 Zürich, Switzerland
   index: 1
 - name: Center for Computational Quantum Physics, Flatiron Institute, 162 5th Avenue, New York, NY 10010, USA
   index: 2
date: 24 June 2022
bibliography: paper.bib

---


# Summary

Strongly correlated systems are a class of materials whose electronic structure is heavily influenced by the effect of electron-electron interactions.
In these systems, an effective single-particle description may not capture the many-body effects accurately.
Although density functional theory (DFT) plus dynamical mean-field theory (DMFT) has proven successful in describing strongly correlated electron systems for over two decades, it has only been very recently that ready-to-use software packages began to become available, with most scientific research carried out by self-written codes developed and used in research groups.
Given the complexity of the method, there is also the question of whether users should implement the formalism themselves for each problem or whether ready-to-use black-box software, as is the case with many DFT software packages, is beneficial to the community.

The goal of solid\_dmft is to find a middle ground, i.e., a *gray-box* tool as a ready-to-use implementation.
This means that while the code contains all the functionality needed for many standard DMFT calculations, it is highly modular, based on open-source and community-developed software, and therefore can be easily adapted to specific applications and needs.
Hence, this project is targeted towards researchers aiming to apply DMFT methods on top of DFT simulations to describe the physics of strongly correlated electron systems.
While our approach allows one to fully perform these computations using standardized input flags without need for coding, the final user can easily extend the functionalities by modifying the corresponding modules in the code.

The package is MPI-parallelized and written in Python 3, utilizing the publicly available TRIQS software library [@triqs:2015] and the applications based on TRIQS, such as different solvers or DFT interfaces.
The philosophy of the package is to increase reproducibility of DFT+DMFT calculations, provide clearer convergence metrics, and allow one to run calculations for a large variety of systems without adapting the code manually, i.e., on a level similar to widely available DFT simulation packages.

## Design Principles

solid_dmft uses the state-of-the-art implementations provided by the TRIQS ecosystem.
This allows the user to easily run ab-initio calculations for strongly correlated materials, as well as implement and test new features of TRIQS and benchmark new TRIQS solvers against existing ones.
solid\_dmft either manages the DFT run itself in charge self-consistency or simply uses the DFT output for one-shot calculations, then creates the downfolded Hamiltonian, solves the resulting Hubbard-like Hamiltonian via DMFT, postprocesses the data to calculate physical observables, and allows feeding back the corrected charge via charge self-consistency.
The full DFT+DMFT cycle is presented in \autoref{fig:downfolding}.

![Fully charge self-consistent DFT+DMFT cycle. Starting from a DFT calculation (top left), a downfolded Hamiltonian and projector functions are created via optimized projections on a local basis set (top right). By adding a specified interaction Hamiltonian $H_\text{int}$, a full interacting electron problem is created, to be solved via the DMFT equations in TRIQS (bottom). After convergence in DMFT is reached, physical observables are calculated (bottom left). For fully charge self-consistent calculations, the DFT cycle is restarted with a DMFT-corrected charge density.\label{fig:downfolding}](downfolding.png){ width=100% }

The code is designed to run on the output of a DFT calculation or a tight-binding model, which provide the low-energy (downfolded) description of a periodic system.
This DFT output is converted into a standardized HDF5 archive by TRIQS/DFTTools [@dfttools:2016] to be used by solid\_dmft so that all DFT codes compatible with TRIQS/DFTTools are supported.
The input for the DMFT calculation can be provided either as a Hamiltonian in reciprocal $\mathbf{k}$-space in a localized basis set or in terms of the overlap between the localized basis set and the Kohn-Sham wavefunctions (so-called projector functions) and their respective eigenvalues.

The code follows the same modular philosophy as the TRIQS software package and relies on TRIQS functionalities to perform basic operations on Green functions.
Each part of the simulation is split into separate stand-alone functions to limit statefulness to a minimum and allows easily extending the functionalities.
The modularity of the program also allows users to run, for example, the DMFT loop via a call of a single pure function with well-defined input and output, i.e., without running solid\_dmft as a monolithic code.
An abstracted `solver` class implements the various impurity solvers available in TRIQS.
solid\_dmft allows one to seamlessly switch between impurity solvers with the change of a simple input flag and by adjusting the solver parameters.
Fully charge self-consistent interfaces are implemented for Quantum ESPRESSO [@Giannozzi_et_al:2009] and the Vienna ab-initio simulation package (VASP) [@kresse_ab_1993; @kresse_efficient_1996].
solid\_dmft also allows users to perform inhomogeneous DMFT calculations, i.e., the treatment of multiple correlated and uncorrelated shells (impurity problems).
Postprocessing scripts are available to perform analytic continuation of imaginary Green functions or self-energies, and to calculate spectral functions.

As of now, solid\_dmft has been successfully used in various peer-reviewed research studies [@Hampel:2019; @Hampel:2020; @Hampel:2021; @merkel_charge_2021; @Beck:2022; @zhang_training_2022].
We provide releases matching those of the TRIQS library, as well as up-to-date documentation with automatic API documentation and tutorials.
Examples and benchmark calculations can be found in the tutorials section of the documentation.
Furthermore, we utilize a continuous-integration workflow on GitHub to test every pull request and commit.

# Statement of need

As of now, only few ready-to-use DFT+DMFT codes are available, all of them released rather recently.
Most of these codes adopt a black-box approach, where the complexity of the DMFT part is abstracted away from the final user (as in EDMFT [@Haule:2010], Amulet [@amulet] or the DMFT implementation in Abinit [@Aldo:2020]) and therefore reduces the number of free parameters to tune. However, this approach may limit the flexibility of the implementation.
solid\_dmft is designed as a more modular, and open source implementation, similar to other software packages like DFTwDMFT [@Singh:2021] and DCORE [@Shinaoka:2021], and acts as a flagship program connecting the TRIQS functionality.
The benefits of this approach are twofold: on the one hand, developers of the TRIQS ecosystem are able to benchmark their applications in a well-tested framework.
On the other hand, users can benefit from a standardized input-output structure compatible with the TRIQS framework, fundamentally increasing robustness and reproducibility.
solid_dmft is developed in the spirit of a community code and supports external contributions that advance the capabilities of the software.

# Acknowledgements

This research was supported by ETH Zurich and the NCCR MARVEL, a National Centre of Competence in Research, funded by the Swiss National Science Foundation (grant number 182892). The Flatiron Institute is a division of the Simons Foundation.

# References
