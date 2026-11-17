#!/usr/bin/env python
# Mixing of the hybridization for solvers with Delta(tau) input
# (greens_functions_mixer.mix_delta). Two Weiss fields with causal discrete
# baths are built with local levels that differ by a chemical-potential shift,
# as in two DMFT iterations. Mixing G0 and subtracting the current levels keeps
# a static part in Delta(iw), and Delta(tau) changes sign; mixing Delta keeps
# Delta(tau) <= 0.

import contextlib
import io
import unittest
from types import SimpleNamespace

import numpy as np

from triqs.gfs import Gf, BlockGf, MeshImFreq, MeshImTime, iOmega_n, inverse, make_gf_from_fourier, fit_hermitian_tail
from triqs.gfs.tools import make_zero_tail
from solid_dmft.dmft_tools import greens_functions_mixer as gf_mixer

beta = 20.0
mesh = MeshImFreq(beta=beta, statistic='Fermion', n_iw=1025)
tau_mesh = MeshImTime(beta=beta, statistic='Fermion', n_tau=4001)
blocks = ['up_0', 'down_0']


def bath_delta(eps_bath, v2):
    D = BlockGf(name_list=blocks, block_list=[Gf(mesh=mesh, target_shape=[1, 1]) for _ in blocks])
    for _, d in D:
        d << sum(v * inverse(iOmega_n - e) for e, v in zip(eps_bath, v2))
    return D


def weiss_field(eps, Delta):
    G0 = Delta.copy()
    for bl, d in Delta:
        G0[bl] << inverse(iOmega_n - eps - d)
    return G0


def delta_tau(Delta_freq):
    """Fourier transform as in dmft_cycle: tail fit with vanishing static part"""
    d = Delta_freq['up_0']
    tail, _ = fit_hermitian_tail(d, make_zero_tail(d, 1))
    return make_gf_from_fourier(d, tau_mesh, tail).data[:, 0, 0].real


eps_prev, eps_new = -1.0, -1.3          # local levels of the two iterations (mu shifted by 0.3)
Delta_prev = bath_delta([-2.0, 1.5], [0.8, 0.5])
Delta_new = bath_delta([-2.5, 1.0], [1.0, 0.4])
G0_prev = weiss_field(eps_prev, Delta_prev)
params = {'g0_mix': 0.5, 'g0_mix_type': 'linear'}


class TestDeltaMixing(unittest.TestCase):

    def test_delta_from_g0(self):
        Delta = gf_mixer._delta_from_g0(G0_prev)
        self.assertLess(np.abs(Delta['up_0'].data - Delta_prev['up_0'].data).max(), 1e-5)

    def test_linear_delta_mixing(self):
        solver = SimpleNamespace(Delta_freq=Delta_new.copy())
        with contextlib.redirect_stdout(io.StringIO()):
            gf_mixer.mix_delta(solver, params, 0, None, G0_prev, [])
        expected = 0.5 * (Delta_new['up_0'].data + Delta_prev['up_0'].data)
        self.assertLess(np.abs(solver.Delta_freq['up_0'].data - expected).max(), 1e-5)
        self.assertTrue(np.all(delta_tau(solver.Delta_freq) < 0))

    def test_g0_mixing_breaks_delta(self):
        # the previous behaviour: mix G0, then Delta = iw - G0^-1 - eps_new
        G0_new = weiss_field(eps_new, Delta_new)
        G0_mix = G0_new.copy()
        G0_mix << 0.5 * G0_new + 0.5 * G0_prev
        Delta = G0_mix.copy()
        for bl, g in G0_mix:
            Delta[bl] << iOmega_n - inverse(g) - eps_new
        self.assertTrue(np.any(delta_tau(Delta) > 0))

    def test_no_mixing(self):
        solver = SimpleNamespace(Delta_freq=Delta_new.copy())
        gf_mixer.mix_delta(solver, {'g0_mix': 1.0, 'g0_mix_type': 'linear'}, 0, None, G0_prev, [])
        self.assertTrue(np.array_equal(solver.Delta_freq['up_0'].data, Delta_new['up_0'].data))


if __name__ == '__main__':
    unittest.main()
