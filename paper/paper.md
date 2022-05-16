---
title: 'solid\_dmft: grayboxing ab initio DFT+DMFT utilizing TRIQS'
tags:
  - Python
  - solid state physics
  - quantum materials
  - correlated electrons
  - dynamical mean field theory
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

* target user group
* has the functionality be verified: list papers that used solid\_dmft 
* describe state of the field
  * EDMFT
  * DFTwDMFT
  * Abinit
  * DCORE
  * Amulet
* how does our software compare / our design philosophy
* TRIQS flagship implementation
* flexible solver choice
* add plot for DMFT loop
* outline functionality in general terms

solid\_dmft is a MPI parallelized scientific simulation code written in Python 3, allowing to perform ab-initio density functional theory (DFT) plus dynamical mean-field theory (DMFT) calculations.
The software is utilizing the TRIQS software library [@triqs:2015],handling most numerical operations. 
Although, DMFT has been very successfully in describing correlated electron systems for over two decades, ready to use software packages are only available very recently, with most scientific research carried out by self-written codes in research groups.   
The aim of solid\_dmft is to provide such ready to use implementation, to increase reproducibility of results, provide clearer convergence metrics, and being able run DMFT calculations for all kind of systems without adapting the code manually, very similar to widely available DFT simulation packages. 
Hence, the targeted user group are researchers that aim to perform DMFT calculations on top of their DFT simulations to describe the physics of strongly correlated electron systems, without the need of elaborate coding, but rather using a standardized input file to control the calculation. 

DFT calculations are performed with one of the TRIQS/DFTTools [@dfttools:2016] compatible codes, with a fully charge self-consistent (CSC) interface implemented for Quantum Espresso and the Vienna ab-initio simulation package (VASP). The DFT output is converted by TRIQS/DFTTools into a HDF archive in a standardized structure to be read by solid_dmft. 


![Caption for example figure.\label{fig:downfolding}](downfolding.png){ width=100% }

The code is mend to run with any input DFT calculation providing a low energy (downfolded) description of the periodic solid system.
This can be either provided directly as Hamiltonian in reciprocal $\mathbf{k}$ space or a overlap between localized basis functions and the DFT wave function. 
Furthermore, the treatment of multiple correlated and uncorrelated shells (impurity problems) is possible. 
To solve the occurring impurity problem in the DMFT loop solid_dmft takes advantage of various impurity solvers available in triqs: cthyb, HubbardI, ForkTPS, ctint, and ctseg.
Even though, these solvers operate very differently solid\_dmft allows to seamlessly switch impurity solvers, with a simple input flag.
After self-consistency is reached, either via full CSC, or just within the DMFT cycle, postprocessing scripts are available to perform analytic continuation of imaginary Green's functions and to calculate spectral functions. 
The full schematic cycle of a DFT+DMFT loop is presented in \autoref{fig:downfolding}.


## Design Principles


# Statement of need

* homebrew codes and a few openly available code
* there only completely black box solutions (EDMFT) or one has to write DFT+DMFT on their own


# Acknowledgements

The Flatiron Institute is a division of the Simons Foundation.

# References
