#!/usr/bin/env python
# Broyden mixing (greens_functions_mixer._broyden_update and mix_gf) on a
# model fixed-point problem G0 = (iw - H - lam G0)^-1, without any solver. Every
# update is compared with an independent implementation of the modified
# Broyden scheme of D. D. Johnson, PRB 38, 12807 (1988), with w_0 = 0.01
# and all weights w_n = 1.

import contextlib
import io
import os
import tempfile
import unittest

import numpy as np

from h5 import HDFArchive
from triqs.gfs import Gf, BlockGf, MeshImFreq
from solid_dmft.dmft_tools import greens_functions_mixer as gf_mixer

mesh = MeshImFreq(beta=10.0, statistic='Fermion', n_iw=30)
H = {'a1': np.array([[0.3, 0.2], [0.2, -0.1]]),
     'a2': np.array([[0.3, 0.2], [0.2, -0.1]]),
     'b': np.array([[0.5]])}
lam = 1.6
w_0 = 0.01


def new_bgf(value=None):
    G = BlockGf(name_list=list(H), block_list=[Gf(mesh=mesh, target_shape=H[bl].shape) for bl in H])
    if value is not None:
        for _, g in G:
            g.data[:] = value
    return G


def dmft_step(G0):
    out = G0.copy()
    for bl, g in out:
        n = H[bl].shape[0]
        for iw in mesh:
            g[iw] = np.linalg.inv(iw.value * np.eye(n) - H[bl] - lam * G0[bl][iw])
    return out


def empty_history():
    return {'mu': [], 'V': [], 'dV': [], 'F': [], 'dF': []}


def broyden_update(*args):
    with contextlib.redirect_stdout(io.StringIO()):
        return gf_mixer._broyden_update(*args)


def johnson_step(V, F, dF, dV, alpha, n_hist):
    """V_{m+1} from the last n_hist normalized differences"""
    if n_hist == 0:
        return V[-1] + alpha * F[-1]
    DF = np.array(dF[-n_hist:])
    DV = np.array(dV[-n_hist:])
    a = DF.conj() @ DF.T
    gamma = np.linalg.solve(a + w_0**2 * np.eye(n_hist), DF.conj() @ F[-1])
    return V[-1] + alpha * F[-1] - gamma @ (alpha * DF + DV)


class TestBroydenG0Mixing(unittest.TestCase):

    def run_broyden(self, n_it, alpha, broy_max_it, deg_shells):
        broyler = empty_history()
        G0 = new_bgf()
        for bl, g in G0:
            for iw in mesh:
                g[iw] = np.linalg.inv(iw.value * np.eye(H[bl].shape[0]))
        residuals = []
        for it in range(1, n_it + 1):
            G0_out = dmft_step(G0)
            residuals.append(max(np.abs(G0_out[bl].data - G0[bl].data).max() for bl in H))
            if it == 1:
                G0 = G0_out
                continue
            n_before = len(broyler['V'])
            G0, broyler = broyden_update(broyler, alpha, broy_max_it, deg_shells, G0_out, G0)
            # V[-1] of this step from the history the step had available
            V, F, dF, dV = (broyler[k] for k in ('V', 'F', 'dF', 'dV'))
            n_avail = len(dF) if broy_max_it == -1 else min(len(dF), broy_max_it - 1)
            ref = johnson_step(V[:-1], F, dF, dV, alpha, n_avail)
            self.assertEqual(len(V), max(n_before, 1) + 1)
            np.testing.assert_allclose(V[-1], ref, rtol=0, atol=1e-12)
            # the vector is put back into every block, degenerate copies included
            np.testing.assert_array_equal(gf_mixer._bgf_to_vec(G0, deg_shells), V[-1])
            if deg_shells:
                np.testing.assert_array_equal(G0['a1'].data, G0['a2'].data)
        return np.array(residuals)

    def test_full_history(self):
        res = self.run_broyden(n_it=25, alpha=0.5, broy_max_it=-1, deg_shells=None)
        self.assertLess(res[-1], 1e-12)

    def test_truncated_history_with_deg_shells(self):
        res = self.run_broyden(n_it=25, alpha=0.5, broy_max_it=3, deg_shells=[['a1', 'a2']])
        self.assertLess(res[-1], 1e-12)

    def test_block_outside_deg_shells_is_mixed(self):
        # first step is linear mixing: 0.5 * 3 + 0.5 * 1 = 2 in every block
        G0, _ = broyden_update(empty_history(), 0.5, -1,
                               [['a1', 'a2']], new_bgf(3.0), new_bgf(1.0))
        for bl in H:
            np.testing.assert_allclose(G0[bl].data, 2.0, atol=1e-15, err_msg=bl)

    def test_empty_history_starts_with_linear_step(self):
        # an archive without history, at any iteration, starts a new sequence
        broyler = empty_history()
        G0, broyler = broyden_update(broyler, 0.3, -1,
                                     None, new_bgf(2.0), new_bgf(1.0))
        for bl in H:
            np.testing.assert_allclose(G0[bl].data, 1.3, atol=1e-15, err_msg=bl)
        self.assertEqual((len(broyler['V']), len(broyler['F']), len(broyler['dF'])), (2, 1, 0))
        G0, broyler = broyden_update(broyler, 0.3, -1,
                                     None, dmft_step(G0), G0)
        self.assertEqual((len(broyler['V']), len(broyler['F']), len(broyler['dF'])), (3, 2, 1))


class TestMixGf(unittest.TestCase):
    """mix_gf, the entry point used in the DMFT cycle for G0 and Sigma"""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.archive = os.path.join(self.tmp.name, 'mix.h5')
        with HDFArchive(self.archive, 'w') as ar:
            ar.create_group('DMFT_results')

    def tearDown(self):
        self.tmp.cleanup()

    def mix(self, gf_out, gf_in, quantity, mix_type, mix_greed, icrsh=0, deg_shells=None):
        params = {'mix_type': mix_type, 'mix_greed': mix_greed, 'broy_max_it': -1}
        with contextlib.redirect_stdout(io.StringIO()):
            return gf_mixer.mix_gf(gf_out, gf_in, params, quantity, icrsh, self.archive, deg_shells)

    def test_linear(self):
        for quantity in ('G0', 'Sigma'):
            mixed = self.mix(new_bgf(3.0), new_bgf(1.0), quantity, 'linear', 0.25)
            for bl in H:
                np.testing.assert_allclose(mixed[bl].data, 1.5, atol=1e-15)
            # mix_greed = 1 returns the output unchanged
            mixed = self.mix(new_bgf(3.0), new_bgf(1.0), quantity, 'linear', 1.0)
            for bl in H:
                np.testing.assert_array_equal(mixed[bl].data, 3.0)

    def test_broyden_histories_per_quantity_and_shell(self):
        # G0 of two shells and Sigma, interleaved in one archive, must each give
        # the same sequence as a standalone Broyden run
        for quantity in ('G0', 'Sigma'):
            gf_mixer.init_broyden_history(self.archive, quantity, n_inequiv_shells=2)
        runs = {('G0', 0): 0.5, ('G0', 1): 0.3, ('Sigma', 0): 0.7}
        V_in = {key: new_bgf(0.1 * (i + 1)) for i, key in enumerate(runs)}
        ref_hist = {key: empty_history() for key in runs}
        ref_in = {key: V_in[key].copy() for key in runs}
        for _ in range(6):
            for (quantity, icrsh), alpha in runs.items():
                key = (quantity, icrsh)
                V_in[key] = self.mix(dmft_step(V_in[key]), V_in[key], quantity, 'broyden', alpha, icrsh=icrsh)
                ref_in[key], ref_hist[key] = broyden_update(ref_hist[key], alpha, -1, None,
                                                            dmft_step(ref_in[key]), ref_in[key])
                for bl in H:
                    np.testing.assert_allclose(V_in[key][bl].data, ref_in[key][bl].data, rtol=0, atol=1e-14)
        with HDFArchive(self.archive, 'r') as ar:
            self.assertEqual([len(h['V']) for h in ar['DMFT_results']['broyler_G0']], [7, 7])
            self.assertEqual(len(ar['DMFT_results']['broyler_Sigma'][0]['V']), 7)
            self.assertEqual(len(ar['DMFT_results']['broyler_Sigma'][1]['V']), 0)

    def test_init_keeps_existing_history(self):
        gf_mixer.init_broyden_history(self.archive, 'Sigma', n_inequiv_shells=1)
        self.mix(new_bgf(3.0), new_bgf(1.0), 'Sigma', 'broyden', 0.5)
        gf_mixer.init_broyden_history(self.archive, 'Sigma', n_inequiv_shells=1)
        with HDFArchive(self.archive, 'r') as ar:
            self.assertEqual(len(ar['DMFT_results']['broyler_Sigma'][0]['V']), 2)


if __name__ == '__main__':
    unittest.main()
