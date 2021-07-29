#!@TRIQS_PYTHON_EXECUTABLE@
################################################################################
#
# solid_dmft - A versatile python wrapper to perform DFT+DMFT calculations
#              utilizing the TRIQS software library
#
# Copyright (C) 2018-2020, ETH Zurich
# Copyright (C) 2021, The Simons Foundation
#      authors: A. Hampel, M. Merkel, and S. Beck
#
# solid_dmft is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# solid_dmft is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# solid_dmft (in the file COPYING.txt in this directory). If not, see
# <http://www.gnu.org/licenses/>.
#
################################################################################

import sys
import time
# triqs
import triqs.utility.mpi as mpi

import solid_dmft.main as run_solid_dmft

# Needed for clean kill of mpi job
def mpiabort_excepthook(type, value, traceback):
    # branch here allows the master node to print the error
    if mpi.is_master_node():
        sys.__excepthook__(type, value, traceback)
    else:
        # wait before calling mpi abort to let the error be printed
        time.sleep(2)
        mpi.MPI.COMM_WORLD.Abort(1)

# this sometimes weirdly suppresses error output completely
sys.excepthook = mpiabort_excepthook
run_solid_dmft.main(sys.argv)


