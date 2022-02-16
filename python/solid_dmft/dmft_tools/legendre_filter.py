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
import numpy as np

# triqs
from triqs.gf import BlockGf
from triqs.gf.tools import fit_legendre


def apply(G_tau, order=100, G_l_cut=1e-19):
    """ Filter binned imaginary time Green's function
    using a Legendre filter of given order and coefficient threshold.

    Parameters
    ----------
    G_tau : TRIQS imaginary time Block Green's function
    auto : determines automatically the cut-off nl
    order : int
        Legendre expansion order in the filter
    G_l_cut : float
        Legendre coefficient cut-off

    Returns
    -------
    G_l : TRIQS Legendre Block Green's function
        Fitted Green's function on a Legendre mesh
    """

    l_g_l = []

    for _, g in G_tau:

        g_l = fit_legendre(g, order=order)
        g_l.data[:] *= (np.abs(g_l.data) > G_l_cut)
        g_l.enforce_discontinuity(np.identity(g.target_shape[0]))

        l_g_l.append(g_l)

    G_l = BlockGf(name_list=list(G_tau.indices), block_list=l_g_l)

    return G_l
