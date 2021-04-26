.. highlight:: bash

.. _install:

Compiling solid_dmft from source
===============================


Prerequisites
-------------

#. The :ref:`TRIQS <triqslibs:welcome>` library, see :ref:`TRIQS installation instruction <triqslibs:installation>`.
   In the following, we assume that TRIQS is installed in the directory ``path_to_triqs``.

Installation steps
------------------

#. Download the source code of the latest stable version by cloning the ``TRIQS/solid_dmft`` repository from GitHub::

     $ git clone https://github.com/TRIQS/solid_dmft solid_dmft.src

#. Create and move to a new directory where you will compile the code::

     $ mkdir solid_dmft.build && cd solid_dmft.build

#. Ensure that your shell contains the TRIQS environment variables by sourcing the ``triqsvars.sh`` file from your TRIQS installation::

     $ source path_to_triqs/share/triqs/triqsvars.sh

#. In the build directory call cmake, including any additional custom CMake options, see below::

     $ cmake ../solid_dmft.src

#. Compile the code, run the tests and install the application::

     $ make test
     $ make install

Version compatibility
---------------------

Keep in mind that the version of ``solid_dmft`` must be compatible with your TRIQS library version,
see :ref:`TRIQS website <triqslibs:versions>`.
In particular the Major and Minor Version numbers have to be the same.
To use a particular version, go into the directory with the sources, and look at all available versions::

     $ cd solid_dmft.src && git tag

Checkout the version of the code that you want::

     $ git checkout 2.1.0

and follow steps 2 to 4 above to compile the code.

Custom CMake options
--------------------

The compilation of ``solid_dmft`` can be configured using CMake-options::

    cmake ../solid_dmft.src -DOPTION1=value1 -DOPTION2=value2 ...

+-----------------------------------------------------------------+-----------------------------------------------+
| Options                                                         | Syntax                                        |
+=================================================================+===============================================+
| Specify an installation path other than path_to_triqs           | -DCMAKE_INSTALL_PREFIX=path_to_solid_dmft      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build in Debugging Mode                                         | -DCMAKE_BUILD_TYPE=Debug                      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Disable testing (not recommended)                               | -DBuild_Tests=OFF                             |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build the documentation                                         | -DBuild_Documentation=ON                      |
+-----------------------------------------------------------------+-----------------------------------------------+
