.. _documentation:

***************
Documentation
***************

Code structure
==============

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

Input/Output
===================
.. toctree::
    :maxdepth: 1

    input_output/DMFT_input/input
    input_output/DMFT_output/results

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
    :toctree: _ref
    :template: autosummary_module_template.rst
    :recursive:

    csc_flow
    dft_managers
    dmft_cycle
    dmft_tools
    gw_embedding
    io_tools
    postprocessing
    util




