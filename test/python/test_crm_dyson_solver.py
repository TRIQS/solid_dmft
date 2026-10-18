#!/usr/bin/env python
# CRM Dyson solver of the solver interfaces: every block of Sigma must be
# obtained, with and without degenerate shells, and blocks outside any
# degenerate shell must not be skipped. The input is an exact model with
# Sigma(iw) = Sigma_0 + V V^+ / (iw - e) per block.

import unittest
import numpy as np

from triqs.gf import Gf, BlockGf, MeshDLRImFreq, make_gf_dlr, inverse
from solid_dmft.dmft_tools.solver import _crm_dyson_solver

beta, w_max, eps = 10.0, 5.0, 1e-8


def model(norb, seed):
    rng = np.random.default_rng(seed)
    E = rng.normal(size=(norb, norb))
    E = 0.5 * (E + E.T)
    S0 = rng.normal(size=(norb, norb))
    S0 = 0.5 * (S0 + S0.T)
    V = rng.normal(size=(norb, norb))
    e = 0.5 * rng.normal()
    return E, S0, V, e


def make_block(mesh, norb, seed):
    E, S0, V, e = model(norb, seed)
    G0 = Gf(mesh=mesh, target_shape=[norb, norb])
    Sig = Gf(mesh=mesh, target_shape=[norb, norb])
    for iw in mesh:
        z = complex(iw)
        G0[iw] = np.linalg.inv(z * np.eye(norb) - E)
        Sig[iw] = V @ V.T / (z - e)
    G = inverse(inverse(G0) - Sig - S0)
    return G0, G, Sig, np.array([S0, V @ V.T], dtype=complex)


class TestCRMDysonSolver(unittest.TestCase):
    def setUp(self):
        self.mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)

    def run_case(self, block_seeds, deg_shells):
        names = list(block_seeds)
        blocks = {name: make_block(self.mesh, norb, seed) for name, (norb, seed) in block_seeds.items()}
        G0_dlr_iw = BlockGf(name_list=names, block_list=[blocks[n][0] for n in names])
        G_dlr = BlockGf(name_list=names, block_list=[make_gf_dlr(blocks[n][1]) for n in names])
        moments = {n: blocks[n][3] for n in names}

        Sigma_dlr = _crm_dyson_solver(G0_dlr_iw, G_dlr, moments, deg_shells)

        for n in names:
            exact = blocks[n][2].data
            self.assertGreater(np.abs(exact).max(), 1e-2)
            np.testing.assert_allclose(Sigma_dlr[n].data, exact, atol=1e-4, err_msg=f'block {n}')

    def test_no_degenerate_shells(self):
        self.run_case({'up_0': (2, 1), 'down_0': (2, 2), 'up_1': (1, 3), 'down_1': (1, 4)}, deg_shells=[])

    def test_partially_degenerate_shells(self):
        # up_0/down_0 degenerate, up_1/down_1 in no shell
        self.run_case(
            {'up_0': (2, 1), 'down_0': (2, 1), 'up_1': (1, 3), 'down_1': (1, 4)}, deg_shells=[['up_0', 'down_0']]
        )

    def test_degenerate_shells(self):
        self.run_case(
            {'up_0': (2, 1), 'down_0': (2, 1), 'up_1': (1, 3), 'down_1': (1, 3)},
            deg_shells=[['up_0', 'down_0'], ['up_1', 'down_1']],
        )


if __name__ == '__main__':
    unittest.main()
