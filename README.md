![logo_soliDMFT](https://raw.githubusercontent.com/triqs/solid_dmft/unstable/doc/logos/logo_solid_dmft.png)

![CI](https://github.com/triqs/solid_dmft/actions/workflows/build.yml/badge.svg)
[![PyPI version](https://badge.fury.io/py/solid_dmft.svg)](https://badge.fury.io/py/solid_dmft)
[![status](https://joss.theoj.org/papers/48eb529b08c6bb464b235ba919d78922/status.svg)](https://joss.theoj.org/papers/48eb529b08c6bb464b235ba919d78922)


---

This program allows to perform DFT+DMFT one-shot and charge self-consistent (CSC) calculations from h5 archives or VASP/Quantum Espresso input files for multiband systems using the [TRIQS](https://triqs.github.io/triqs/latest/) software library, and the DFT code interface [TRIQS/DFTTools](https://triqs.github.io/dft_tools/latest/). solid_dmft takes advantage of various impurity solvers available in [TRIQS](https://triqs.github.io/triqs/unstable/applications.html#impurity-solvers). Postprocessing scripts are available to perform analytic continuation and calculate spectral functions.

### Documentation & tutorials

To learn how to use solid_dmft, take a look at the [online documentation](https://triqs.github.io/solid_dmft/). There you can find:
* [input / output documentation](https://triqs.github.io/solid_dmft/documentation.html#input-output)
* [reference manual of functions](https://triqs.github.io/solid_dmft/documentation.html#module-reference-manual)
* [code structure](https://triqs.github.io/solid_dmft/documentation.html#code-structure)
* [tutorials](https://triqs.github.io/solid_dmft/tutorials.html)

Check also the [solid_dmft publication](https://doi.org/10.21105/joss.04623) in the JOSS journal for more information and further references.

### Installation

You can install the latest solid_dmft release simply via pip (PyPi):
```
pip install solid_dmft
```
However, please make sure that you have a valid TRIQS and TRIQS/DFTTools installation matching the version of solid_dmft. Furthermore, you need at least one of the supported DMFT impurity solvers installed to use solid_dmft. 

A more thorough installation can be performed manually via `cmake`, which will also check if you have a working and matching TRIQS installation.

Please check the [installation page](https://triqs.github.io/solid_dmft/install.html) on the online documentation for more detailed instructions.

---

Copyright (C) 2018-2020, ETH Zurich
Copyright (C) 2021, The Simons Foundation 
  authors: A. Hampel, S. Beck, M. Merkel, and A. Carta
(see LICENSE.txt for details)

Developed by the Materials Theory Group, ETH Zurich
and the Center for Computational Quantum Physics, Flatiron Institute.

If you are using this code for your research, please cite it with this
[bib file](https://github.com/TRIQS/solid_dmft/blob/unstable/cite_solid_dmft.bib).

