.. highlight:: bash
.. _installation:

Installation 
#############

Prerequisites
-------------

#. The :ref:`TRIQS <triqslibs:welcome>` library, see `TRIQS installation instruction <https://triqs.github.io/triqs/latest/install.html>`_.
   In the following, we assume that :ref:`TRIQS <triqslibs:welcome>`, `triqs/dft_tools <https://triqs.github.io/dft_tools>`_, and at least one of the impurity solvers `available in TRIQS <https://triqs.github.io/triqs/latest/applications.html>`_, e.g. cthyb, HubbardI, ctseg, FTPS, or ctint is installed in the directory ``path_to_triqs``.

#. Make sure to install besides the triqs requirements also the python packages::

     $ pip3 install --user scipy argparse pytest

#. To build the documentation the following extra python packages are needed::

     $ pip3 install --user sphinx sphinx-autobuild pandoc nbsphinx linkify-it-py sphinx_rtd_theme myst-parser


Installation via pip
--------------------

You can install the latest solid_dmft release simply via pip (PyPi):
```
pip install solid_dmft
```
However, please make sure that you have a valid TRIQS and TRIQS/DFTTools installation matching the version of solid_dmft. Furthermore, you need at least one of the supported DMFT impurity solvers installed to use solid_dmft. 

Manual installation via CMake
-----------------------------

We provide hereafter the build instructions in the form of a documented bash script. Please change the variable INSTALL_PREFIX to point to your TRIQS installation directory::
    
    INSTALL_PREFIX=/path/to/triqs
    # source the triqsvars.sh file from your TRIQS installation to load the TRIQS environment
    source $(INSTALL_PREFIX)/share/triqs/triqsvars.sh

    # clone the flatironinstitute/solid_dmft repository from GitHub
    git clone https://github.com/flatironinstitute/solid_dmft solid_dmft.src

    # checkout the branch of solid_dmft matching your triqs version. 
    # For example if you use the 3.1.x branch of triqs, dfttools. and cthyb
    git checkout 3.1.x

    # Create and move to a new directory where you will compile the code
    mkdir solid_dmft.build && cd solid_dmft.build

    # In the build directory call cmake, including any additional custom CMake options, see below
    cmake ../solid_dmft.src

    # Compile the code, run the tests, and install the application
    make test
    make install

This installs solid_dmft into your TRIQS installation folder.

To build ``solid_dmft`` with documentation you should run:: 

     $ cmake path/to/solid_dmft.src -DBuild_Documentation=ON 
     $ make 
     $ sphinx-autobuild path/to/solid_dmft.src/doc ./doc/html -c ./doc/

The last line will automatically search for changes in your src dir, rebuild the documentation, 
and serve it locally as under `127.0.0.1:8000`. 

Docker files & images
---------------------

We `provide docker files <https://github.com/TRIQS/solid_dmft/tree/3.1.x/Docker>`_ to build solid_dmft inside a docker container with all dependencies and instructions on how to integrate the connected DFT codes as well. Additionally, we host a most recent unstable version of the docker image used for the github CI `on dockerhub <https://hub.docker.com/r/materialstheory/solid_dmft_ci>`_. To use this version, which includes the cthyb solver, the hubbardI solver, dfttools, and the maxent package, pull the following image::

    $ docker pull materialstheory/solid_dmft_ci


Version compatibility
---------------------

Keep in mind that the version of ``solid_dmft`` must be compatible with your TRIQS library version,
see :ref:`TRIQS website <triqslibs:versions>`.
In particular the Major Version numbers have to be the same.
To use a particular version, go into the directory with the sources, and look at all available branches::

     $ cd solid_dmft.src && git branch -vv

Checkout the version of the code that you want::

     $ git checkout 3.1.x

and follow steps 3 to 6 above to compile the code.

Custom CMake options
--------------------

The compilation of ``solid_dmft`` can be configured using CMake-options::

    cmake ../solid_dmft.src -DOPTION1=value1 -DOPTION2=value2 ...

+-----------------------------------------------------------------+-----------------------------------------------+
| Options                                                         | Syntax                                        |
+=================================================================+===============================================+
| Specify an installation path other than path_to_triqs           | -DCMAKE_INSTALL_PREFIX=path_to_solid_dmft     |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build in Debugging Mode                                         | -DCMAKE_BUILD_TYPE=Debug                      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Disable testing (not recommended)                               | -DBuild_Tests=OFF                             |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build the documentation                                         | -DBuild_Documentation=ON                      |
+-----------------------------------------------------------------+-----------------------------------------------+
