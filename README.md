![logo_soliDMFT](doc/logos/logo_solid_dmft.png)

This program allows to perform DFT+DMFT ''one-shot'' and CSC
calculations from h5 archives, VASP, ro Quantum Espresso for multiband systems using
the TRIQS package, in combination with various impurity solvers and SumkDFT from
DFT-tools. 

Runs with triqs 3.x.x

Copyright (C) 2021: A. Hampel, M. Merkel, and S. Beck
(see LICENSE.txt for details)


## Source code files and their use

- __bin/solid_dmft:__ main file that runs the DMFT calculation and starts a CSC flow 
- __main.py:__ main function that invokes `csc_flow_control` or a one shot 
  calculation directly by invoking `dmft_cycle` on a given h5 archives
- __read_config.py:__ contains the functions to read the dmft config file. Take a
  look in `read_config_doc.md` for a detailed list of parameters
- __dmft_cycle.py:__ contains the `dmft_cycle` function that run a predefined
  number of DMFT iterations
- __csc_flow.py:__ contains the `csc_flow_control` function to steer the CSC
  calculation and call then ones per DFT+DMFT cycle the `dmft_cycle` function
- __/dmft_tools__
    - __convergence.py:__ contains functions to calculate convergence properties and
      defining definitions for Gf differences. Results will be stored in a dictionary:
      `convergence_obs`
    - __formatter.py:__
    - __interaction_hamiltonian.py:__
    - __legendre_filter.py:__
    - __manipulate_chemical_potential.py:__
    - __observables.py:__ contains all functions necessary to calculate and write the
      observables during the run, which are stored in a general dictionary: `observables`
    - __solver.py:__
- __/dft_managers.py__
    - __vasp_manager.py:__ contains all functions to communicate with vasp,
    which in a CSC calculation in running the whole time (in the background,
    if not needed)
    - __qe_manager.py:__ contains all functions to start quantum espresso
- __/util:__ non-essential scripts that are independent of the rest of the code
- __/postprocessing:__ different scripts to perform postprocessing steps like analytical 
  continuations with triqs's MaxEnt


## Getting started

To start take a look in the `example` directory. There one finds several
examples to run. Best start with the svo-one-shot example. The
`dmft_config.ini` file contains the configuration for the DMFT run, which is
explained in the read\_config method in the main script. The `svo.h5` is the DMFT
input data, which is obtained from projection on localized Wannier functions
(see folder `svo-one-shot/converter`).

Furthermore, there is a `read_config_doc.md` file containing the docstrings from
the main script in a readable format. If one wishes to do CSC calculations the
docker container must contain also a installed VASP version >5.4.4 that
understands the ICHARG=5 flag.

To test triqs, you can use the official image from https://hub.docker.com/r/flatironinstitute/triqs/.
For a full installation including Vasp run our own Docker image (see folder `/docker`).
