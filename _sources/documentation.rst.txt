.. _documentation:

***************
Documentation
***************

Code structure:
===============

.. image:: _static/code_structure.png
   :width: 100%
   :align: center

more details in the reference manual below.

To get started with the code after a successful :ref:`installation`, take a look at the :ref:`tutorials` section. Here we provide further special information and a reference manual for all available functions.

DFT interface notes
===================

.. toctree::
    :maxdepth: 1

    md_notes/w90_interface.md
    md_notes/vasp_csc.md
    cRPA_VASP/README.md

Further details for running
===========================
   
.. toctree::
    :maxdepth: 1
    
    md_notes/docker.md
    md_notes/run_locally.md
    md_notes/run_cluster.md

Module reference manual
=======================

.. autosummary::
    :toctree: _autosummary
    :template: autosummary_module_template.rst
    :recursive:

    csc_flow
    dft_managers
    dmft_cycle
    dmft_tools
    postprocessing
    read_config
    util
    

   

