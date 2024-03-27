.. _tutorials:

Tutorials
==========

These tutorials provide an overview about typical workflows to perform DFT+DMFT calculations with solid_dmft. The tutorials are sorted by complexity and introduce one after another more available features. 

.. note::
   The tutorials are run with the 3.1.x branch of triqs. Please use the 3.1.x branch for triqs and all applications to reproduce the results shown here. 

Short description of the tutorials linked below: 

1. Typical one-shot (OS) DMFT calculation based on prepared hdf5 archive for SrVO3
2. Full charge self-consistent (CSC) DFT+DMFT calculation using the PLO formalism with Vasp for PrNiO3
3. Full CSC DFT+DMFT calculation using w90 in combination with Quantum Espresso utilizing the lighter HubbardI solver
4. OS magnetic DMFT calculation for NdNiO2 in a large energy window for 5 d orbitals
5. Postprocessing: plot the spectral function after a DFT+DMFT calculation

----

.. toctree::
    :maxdepth: 2

    tutorials/SVO_os_qe/tutorial
    tutorials/PrNiO3_csc_vasp_plo_cthyb/tutorial
    tutorials/Ce2O3_csc_w90/tutorial
    tutorials/NNO_os_plo_mag/tutorial
    tutorials/correlated_bandstructure/plot_correlated_bands
