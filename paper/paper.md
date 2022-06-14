---
title: 'solid\_dmft: gray-boxing ab initio DFT+DMFT utilizing TRIQS'
tags:
  - Python
  - solid state physics
  - quantum materials
  - correlated electrons
  - dynamical mean-field theory
  - Anderson impurity problem
  - DMFT
authors:
  - name: Maximilian E. Merkel
    orcid: 0000-0003-2589-9625
    affiliation: 1
  - name: Alberto Carta
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
date: 18 April 2022
bibliography: paper.bib

---


# Summary

Although density functional theory (DFT) plus dynamical mean-field theory (DMFT) has proven successful in describing correlated electron systems for over two decades, it has only been very recently that ready-to-use software packages began to become available, with most scientific research carried out by self-written codes developed and used in research groups.
Given the complexity of the method, there is also the question of whether a black-box approach is beneficial to the community or whether users should be able to implement the formalism themselves.
The goal of solid\_dmft is to find a middle ground, i.e. a gray-box tool as a ready-to-use implementation.
The package is an MPI-parallelized scientific simulation code written in Python 3, allowing to perform ab initio DFT+DMFT calculations.
solid\_dmft utilizes the TRIQS software library [@triqs:2015], which handles most numerical operations.
The philosophy of the package is to increase reproducibility of DFT+DMFT calculations, provide clearer convergence metrics, and allow to run calculations for a large variety of systems without adapting the code manually, i.e. on this level similar to widely available DFT simulation packages.
Hence, the targeted user group of this software project is researchers that aim to perform DMFT calculations on top of their DFT simulations to describe the physics of strongly correlated electron systems, without the need of elaborate coding, but rather using a standardized input file to control the calculation.
On the other hand, the code is highly modular, based on open-source and community-developed software and therefore easily adjustable for specific applications and needs.

## Design Principles


The idea is to provide the full functionality of a full a DFT+DMFT calculation by merging the state-of-the-art implementations provided by the TRIQS library and its applications.
This allows to easily run ab initio calculations for strongly correlated materials, as well as implement and test new features of TRIQS or benchmark new solvers against existing ones.
solid\_dmft manages the calls of the necessary routines to run the DFT calculations, to create the downfolded Hamiltonian, solve the resulting Hubbard-like Hamiltonian via DMFT, postprocess the data to calculate physically meaningful observables, and allow for charge-corrected feedback via charge self-consistency.
The full DFT+DMFT cycle is presented in \autoref{fig:downfolding}.

![Full CSC DFT+DFMT cycle. Starting from a DFT calculation (top left), a downfolded Hamiltonian and projector functions are created via optimized projections on a local basis set (top right). By adding a specified interaction Hamiltonian $H_\text{int}$, a full interacting electron problem is created, to be solved via the DMFT equations (bottom) in TRIQS. After convergence in DMFT is reached, physical observables are calculated (bottom left). For fully CSC calculations a charge density correction is added to DFT and the cycle is restarted.\label{fig:downfolding}](downfolding.png){ width=100% }

The code is designed to run on top of a DFT calculation or model system providing a low-energy (downfolded) description of the periodic solid system.
The DFT calculations can be performed with any code that is compatible with TRIQS/DFTTools [@dfttools:2016].
The input for the DMFT calculation can be either provided directly as a Hamiltonian in reciprocal $\mathbf{k}$-space in a localized basis set, or in terms of the overlap between the localized basis set and the Kohn-Sham wavefunctions (so-called projector functions), and their respective eigenvalues.
The DFT output is converted by TRIQS/DFTTools into an HDF archive in a standardized structure to be read or called by solid\_dmft.
The code is designed to be modular in the same philosophy as the TRIQS software package, relying on TRIQS functionalities to perform basic operations.
Therefore, we split each part of the simulation into separate stand-alone functions, to limit statefulness to a minimum and allow for an easy extension to include new features.
The modularity of the program also allows to run, for example, the DMFT loop only via a call of a single function with well-defined input and output, i.e. without running solid\_dmft as a monolithic code.
This ensures that the code can be used in other projects.
A abstracted `solver` class implements the various impurity solvers available in triqs: cthyb, HubbardI, ForkTPS, ctint, and ctseg.
Even though these solvers operate differently, solid\_dmft allows to seamlessly switch impurity solvers, with a simple input flag and adjusting the solver paramters.
A fully charge self-consistent (CSC) interface is implemented for Quantum ESPRESSO and the Vienna ab-initio simulation package (VASP).
solid\_dmft allows also to perform inhomogenous DMFT calculations, i.e. the treatment of multiple correlated and uncorrelated shells (impurity problems) while converging the full lattice self-energy.
After self-consistency is reached, either via full CSC or just within the DMFT cycle, postprocessing scripts are available to perform analytic continuation of imaginary Green's functions, and to calculate spectral functions.

As of now, solid\_dmft has been successfully used in various peer-reviewed research studies [@Beck:2022; @Hampel:2019; @Hampel:2020; @Hampel:2021], and provides stable releases matching the releases of the TRIQS library.
We provide a full documentation including several tutorials and a reference manual.
Examples and benchmark calculations can be found in the tutorials section of the documentation.
Furthermore, we utilize an extensive CI workflow on github to test every pull request and commit.

# Statement of need

The number of ready-to-use DFT+DMFT codes is small, and all codes have been developed rather recently.
Some of these operate in a black-box way, such as, for example, EDMFT, Amulet[@amulet] and the DMFT implementation included in Abinit[@Aldo:2020].
Other software packages like DFTwDMFT[@Singh:2021] and DCORE[@Shinaoka:2021] follow a very similar strategy as solid\_dmft but have not implemented different impurity solvers as of now.
solid\_dmft provides a flagship implementation of the TRIQS functionality to perform DFT+DMFT calculation and is easily extended to beyond-DMFT formalisms.
This is beneficial both to developers of TRIQS applications, by enabling them to benchmark their applications in a well-tested framework, and to users, who benefit from the most recent features of TRIQS.
Via the gray-box approach the software provides a robust and flexible implementation of the DFT+DMFT method, controlled via a single input file.
It is developed in the spirit of a community code and supports external contributions that advance the capabilities of the program.

# Acknowledgements

The Flatiron Institute is a division of the Simons Foundation.

# References
