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

This program allows to perform DFT+DMFT ‘’one-shot’’ and charge self-consistent (CSC) calculations from h5 archives or VASP/Quantum Espresso input files for multiband systems using the TRIQS software library, and the DFT code interface TRIQS/DFTTools. Works with triqs >3.x.x. solid_dmft takes advantage of various impurity solvers available in triqs: cthyb, HubbardI, ForkTPS, ctint, and ctseg. Postprocessing scripts are available to perform analytic continuation and calculate spectral functions.

## Design Principles

![Caption for example figure.\label{fig:downfolding}](downfolding.png){ width=100% }

# Statement of need

* homebrew codes and a few openly available code
* there only completely black box solutions (EDMFT) or one has to write DFT+DMFT on their own


# Acknowledgements

The Flatiron Institute is a division of the Simons Foundation.

# References
