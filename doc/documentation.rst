.. _documentation:

***************
documentation
***************

code structure:
===============

.. image:: _static/code_structure.png
   :width: 100%
   :align: center

more details in the reference manual below.

DFT code interfaces
===================

.. toctree::
   md_notes/w90_interface.md
   md_notes/vasp_csc.md
   cRPA_VASP/README.md

run solid_dmft
===================
   
   .. toctree::
      md_notes/docker.md
      md_notes/run_locally.md
      md_notes/run_cluster.md

module reference manual
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
    

   

