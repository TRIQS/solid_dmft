.. highlight:: bash
.. _installation:

Installation 
#############

Prerequisites
-------------

#. The :ref:`TRIQS <triqslibs:welcome>` library, see :ref:`TRIQS installation instruction <triqslibs:installation>`.
   In the following, we assume that TRIQS is installed in the directory ``path_to_triqs``.

#. make sure to install besides the triqs requirements also the python packages::

     $ pip3 install --user scipy argparse

Installation steps
------------------

#. Download the source code of the latest stable version by cloning the ``flatironinstitute/solid_dmft`` repository from GitHub::

     $ git clone https://github.com/flatironinstitute/solid_dmft solid_dmft.src

#. Create and move to a new directory where you will compile the code::

     $ mkdir solid_dmft.build && cd solid_dmft.build

#. Ensure that your shell contains the TRIQS environment variables by sourcing the ``triqsvars.sh`` file from your TRIQS installation::

     $ source path_to_triqs/share/triqs/triqsvars.sh

#. In the build directory call cmake, including any additional custom CMake options, see below::

     $ cmake ../solid_dmft.src

#. Compile the code, run the tests and install the application::

     $ make test
     $ make install


to build ``solid_dmft`` with documentation you should run:: 

     $ cmake path/to/solid_dmft.src -DBuild_Documentation=ON 
     $ make 
     $ sphinx-autobuild path/to/solid_dmft.src/doc ./doc/html -c ./doc/

The last line will automatically search for changes in your src dir, rebuild the documentation, 
and serve it locally as under `127.0.0.1:8000`. To build the documentation the following extra 
python packages are needed::

     $ pip3 install --user sphinx sphinx-autobuild pandoc nbsphinx linkify-it-py sphinx_rtd_theme myst-parser

Version compatibility
---------------------

Keep in mind that the version of ``solid_dmft`` must be compatible with your TRIQS library version,
see :ref:`TRIQS website <triqslibs:versions>`.
In particular the Major and Minor Version numbers have to be the same.
To use a particular version, go into the directory with the sources, and look at all available versions::

     $ cd solid_dmft.src && git tag

Checkout the version of the code that you want::

     $ git checkout 3.0.x

and follow steps 2 to 4 above to compile the code.

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
