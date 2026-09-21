# Copyright (c) 2018-2022 Simons Foundation
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You may obtain a copy of the License at
#     https:#www.gnu.org/licenses/gpl-3.0.txt
#
# Authors: Alexander Hampel

import shutil
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from h5 import HDFArchive
from solid_dmft.postprocessing import plot_correlated_bands as pcb

import unittest


def _write_bloch_basis_archive(src, dst, n_extra_bands=0, scramble_phases=True, seed=1234):
    """
    Rewrites a Wannier-basis archive as the equivalent Bloch-basis one.

    For an isolated set of bands the Bloch-basis representation of a wannier90
    run is fully determined by the Wannier one: the KS eigenvalues of the
    Wannier bands are the eigenvalues of H_W(k) and the gauge matrix is the
    corresponding eigenvector matrix V(k), so that P H P^dag = V diag(e) V^dag
    = H_W(k). The per-band phase is left free by wannier90 as well and is
    randomised here, which makes P(k) vary discontinuously along k.

    With n_extra_bands > 0 additional decoupled bands far outside the window
    are added, so that the number of bands exceeds the number of Wannier
    functions as it does for a run with band disentanglement.
    """

    shutil.copyfile(src, dst)
    rng = np.random.default_rng(seed)
    with HDFArchive(dst, 'a') as ar:
        dft_input = ar['dft_input']
        hopping = dft_input['hopping']
        n_k, n_sp, n_orb, _ = hopping.shape
        n_bands = n_orb + n_extra_bands

        hopping_bloch = np.zeros((n_k, n_sp, n_bands, n_bands), dtype=complex)
        proj_mat_bloch = np.zeros((n_k, n_sp, dft_input['n_corr_shells'], n_orb, n_bands), dtype=complex)

        for ik in range(n_k):
            for isp in range(n_sp):
                evals, evecs = np.linalg.eigh(hopping[ik, isp])
                if scramble_phases:
                    evecs = evecs * np.exp(1j * rng.uniform(0, 2 * np.pi, n_orb))[None, :]
                hopping_bloch[ik, isp, :n_orb, :n_orb] = np.diag(evals)
                proj_mat_bloch[ik, isp, 0, :, :n_orb] = evecs
                if n_extra_bands > 0:
                    hopping_bloch[ik, isp, n_orb:, n_orb:] = np.diag(
                        20.0 + np.arange(n_extra_bands) + 0.1 * rng.standard_normal(n_extra_bands))

                assert np.allclose(proj_mat_bloch[ik, isp, 0] @ hopping_bloch[ik, isp]
                                   @ proj_mat_bloch[ik, isp, 0].conj().T, hopping[ik, isp])

        dft_input['hopping'] = hopping_bloch
        dft_input['proj_mat'] = proj_mat_bloch
        dft_input['n_orbitals'] = np.full((n_k, n_sp), n_bands)
        dft_input['dft_code'] = 'w90'

    return dst


class _FakeSumk:
    """Minimal stand-in for SumkDFT for the embedding checks."""

    def __init__(self, dims, hopping, proj_mat, dft_code='w90'):
        self.corr_shells = [{'dim': dim} for dim in dims]
        self.n_corr_shells = len(dims)
        self.hopping = hopping
        self.proj_mat = proj_mat
        self.n_k = hopping.shape[0]
        self.n_orbitals = np.full((self.n_k, 1), hopping.shape[2])
        self.bz_weights = np.full(self.n_k, 1.0 / self.n_k)
        self.dft_code = dft_code


class test_convergence(unittest.TestCase):

    def setUp(self):

        self.w90_dict = {'w90_seed': 'svo', 'w90_path': './', 'mu_tb': 12.3958, 'n_orb': 3,
                         'orbital_order_w90': ['dxz', 'dyz', 'dxy'], 'add_spin': False}

        self.orbital_order_to = ['dxy', 'dxz', 'dyz']

        self.sigma_dict = {'dmft_path': './svo_example.h5', 'it': 'last_iter',
                           'orbital_order_dmft': self.orbital_order_to, 'spin': 'up',
                           'block': 0, 'eta': 0.0, 'linearize': False}

    def test_get_dmft_bands(self):
        tb_bands = {'bands_path': [('R', 'G'), ('G', 'X'), ('X', 'M'), ('M', 'G')], 'G': [0., 0., 0.],
                    'Z': np.array([0, 0, 0.5]), 'M': [0.5, 0.5, 0.], 'R': [0.5, 0.5, 0.5],
                    'X': [0.,  0.5, 0.], 'n_k': 50}

        tb_data, alatt_k_w, freq_dict = pcb.get_dmft_bands(with_sigma='calc', add_mu_tb=True,
                                                           orbital_order_to=self.orbital_order_to,
                                                           **self.w90_dict, **tb_bands, **self.sigma_dict)

        with HDFArchive('test_pcb_ref.h5', 'r') as ar:
            emat_ref = ar['tb_emat']
            Akw_ref = ar['Akw']

        assert np.allclose(tb_data['e_mat'], emat_ref)
        assert np.allclose(alatt_k_w, Akw_ref)

    def test_get_dmft_bands_proj(self):
        tb_bands = {'bands_path': [('R', 'G'), ('G', 'X'), ('X', 'M'), ('M', 'G')], 'G': [0., 0., 0.],
                    'Z': np.array([0, 0, 0.5]), 'M': [0.5, 0.5, 0.], 'R': [0.5, 0.5, 0.5],
                    'X': [0.,  0.5, 0.], 'n_k': 50}

        tb_data, alatt_k_w, freq_dict = pcb.get_dmft_bands(with_sigma='calc', add_mu_tb=True,
                                                           orbital_order_to=self.orbital_order_to,
                                                           proj_on_orb = [0,1],
                                                           **self.w90_dict, **tb_bands, **self.sigma_dict)

        with HDFArchive('test_pcb_ref.h5', 'r') as ar:
            emat_ref = ar['tb_emat_proj']
            Akw_ref = ar['Akw_proj']

        assert np.allclose(tb_data['e_mat'], emat_ref)
        assert np.allclose(alatt_k_w, Akw_ref)

    def test_get_kslice(self):

        freq_mesh_kslice = {'window': [-0.5, 0.5], 'n_w': int(1e6)}
        sigma_dict = {'dmft_path': './svo_example.h5', 'it': 'last_iter', 'w_mesh': freq_mesh_kslice,
                      'orbital_order_dmft': self.orbital_order_to, 'spin': 'up',
                      'block': 0, 'eta': 0.0, 'linearize': False}

        # use kz to create non kz=0.0 slice first and then explicity coordinates
        tb_kslice = {'bands_path': [('Y', 'G'), ('G', 'X')],
                     'Y': np.array([0.5, 0.0, 0]), 'G': [0., 0., 0.],
                     'Z': np.array([0, 0, 1.0]), 'X': [0., 0.5, 0.],
                     'n_k': 50, 'kz': 0.34}

        tb_data_kz, alatt_k_w_kz, _ = pcb.get_dmft_bands(fermi_slice=True, with_sigma='calc', add_mu_tb=True,
                                                         orbital_order_to=self.orbital_order_to,
                                                         **self.w90_dict, **tb_kslice, **sigma_dict)

        # now explicittly define the coordinates
        tb_kslice = {'bands_path': [('YZ', 'GZ'), ('GZ', 'XZ')],
                     'YZ': np.array([0.5, 0.0, 0.34]), 'GZ': [0., 0., 0.34], 'XZ': [0., 0.5, 0.34],
                     'n_k': 50}

        tb_data_exp, alatt_k_w_exp, _ = pcb.get_dmft_bands(fermi_slice=True, with_sigma='calc', add_mu_tb=True,
                                                           orbital_order_to=self.orbital_order_to,
                                                           **self.w90_dict, **tb_kslice, **sigma_dict)


        assert np.allclose(tb_data_kz['e_mat'], tb_data_exp['e_mat'])
        assert np.allclose(alatt_k_w_kz, alatt_k_w_exp)

    def test_get_kslice_nokz(self):

        freq_mesh_kslice = {'window': [-0.5, 0.5], 'n_w': int(1e6)}
        sigma_dict = {'dmft_path': './svo_example.h5', 'it': 'last_iter', 'w_mesh': freq_mesh_kslice,
                      'orbital_order_dmft': self.orbital_order_to, 'spin': 'up',
                      'block': 0, 'eta': 0.0, 'linearize': False}

        tb_kslice = {'bands_path': [('Y', 'G'), ('G', 'X')], 'Y': np.array([0.5, 0.0, 0]), 'G': [0., 0., 0.],
                     'M': [0.5, 0.5, 0.], 'R': [0.5, 0.5, 0.5],
                     'X': [0.,  0.5, 0.], 'n_k': 50}

        tb_data, alatt_k_w, freq_dict = pcb.get_dmft_bands(fermi_slice=True, with_sigma='calc', add_mu_tb=True,
                                                           orbital_order_to=self.orbital_order_to,
                                                           **self.w90_dict, **tb_kslice, **sigma_dict)

        with HDFArchive('test_pcb_ref.h5', 'r') as ar:
            emat_ref = ar['tb_emat_slice']
            Akw_ref = ar['Akw_slice']

        assert np.allclose(tb_data['e_mat'], emat_ref)
        assert np.allclose(alatt_k_w, Akw_ref)

    def test_get_dmft_bands_reg_mesh_read_TB_obj(self):
        tb_bands = {'kmesh': 'regular', 'n_k': 7}

        tb_data, alatt_k_w, freq_dict = pcb.get_dmft_bands(with_sigma='calc', add_mu_tb=True,
                                                           orbital_order_to=self.orbital_order_to,
                                                           **self.w90_dict, **tb_bands, **self.sigma_dict)

        with HDFArchive('test_pcb_ref.h5', 'a') as ar:
            emat_ref = ar['tb_emat_reg_mesh']
            Akw_ref = ar['Akw_reg_mesh']

        assert np.allclose(tb_data['e_mat'], emat_ref)
        assert np.allclose(alatt_k_w, Akw_ref)

        # read now from TB_obj
        w90_dict = {'TB_obj': tb_data['tb'], 'mu_tb': 12.3958, 'n_orb': 3,
                         'orbital_order_w90': ['dxz', 'dyz', 'dxy']}

        tb_data_obj, alatt_k_w_obj, freq_dict_obj = pcb.get_dmft_bands(with_sigma='calc', add_mu_tb=True,
                                                           orbital_order_to=self.orbital_order_to,
                                                           **w90_dict, **tb_bands, **self.sigma_dict)

        assert np.allclose(tb_data_obj['e_mat'], emat_ref)
        assert np.allclose(alatt_k_w_obj, Akw_ref)


    def test_plot_bands_with_projection(self):
        """
        plot_bands colours the tight-binding bands with the orbital character it
        reads from tb_data['proj_nuk'], and draws them from tb_data['e_mat'],
        which has to be the eigenvalue matrix. Both have to survive an orbital
        projection, also when no self-energy is given at all, in which case the
        projection affects nothing else.
        """

        tb_bands = {'bands_path': [('R', 'G'), ('G', 'X')], 'G': [0., 0., 0.],
                    'Z': np.array([0, 0, 0.5]), 'R': [0.5, 0.5, 0.5],
                    'X': [0.,  0.5, 0.], 'n_k': 30}

        # no self-energy: get_dmft_bands returns alatt_k_w = None and the
        # projection only feeds the band colouring
        tb_data, alatt_k_w, freq_dict = pcb.get_dmft_bands(orbital_order_to=self.orbital_order_to,
                                                           proj_on_orb=[0, 1], eta=0.01,
                                                           **self.w90_dict, **tb_bands)
        assert alatt_k_w is None
        assert isinstance(tb_data['proj_nuk'], np.ndarray), 'orbital character missing from tb_data'
        assert tb_data['proj_nuk'].shape == (3, tb_data['e_mat'].shape[2])
        # orbital character of a complete set of bands sums to one per band
        full = pcb.get_dmft_bands(orbital_order_to=self.orbital_order_to, proj_on_orb=[0, 1, 2],
                                  eta=0.01, **self.w90_dict, **tb_bands)[0]
        assert np.allclose(full['proj_nuk'], 1.0)

        fig, ax = plt.subplots()
        try:
            pcb.plot_bands(fig, ax, alatt_k_w, tb_data, freq_dict, n_orb=3, tb=True,
                           alatt=False, colorscheme_bands='coolwarm')
        finally:
            plt.close(fig)

        # and the same with a self-energy, where alatt is plotted on top
        sigma_dict = dict(self.sigma_dict)
        tb_data, alatt_k_w, freq_dict = pcb.get_dmft_bands(with_sigma='calc', add_mu_tb=True,
                                                           orbital_order_to=self.orbital_order_to,
                                                           proj_on_orb=[0, 1],
                                                           **self.w90_dict, **tb_bands, **sigma_dict)
        assert isinstance(tb_data['proj_nuk'], np.ndarray)
        fig, ax = plt.subplots()
        try:
            pcb.plot_bands(fig, ax, alatt_k_w, tb_data, freq_dict, n_orb=3, tb=True,
                           alatt=True, colorscheme_bands='coolwarm',
                           colorscheme_alatt='magma')
        finally:
            plt.close(fig)

    def test_orbital_projection_is_exact(self):
        """
        An orbital projection is a mask on the diagonal of the lattice Green
        function in the Wannier basis, so summing the projections over all
        orbitals has to reproduce the orbital-resolved spectral function exactly.

        Going through the band basis instead weights the band-resolved spectral
        function with the orbital character |<orb|band>|^2, which drops the
        off-diagonal band components of G. That is exact only for an orbital
        independent self-energy, which is why the SVO archive used by the other
        tests cannot detect the difference: its Sigma is isotropic to 1e-14.
        Both conventions give the same trace.
        """

        tb_bands = {'bands_path': [('R', 'G'), ('G', 'X')], 'G': [0., 0., 0.],
                    'Z': np.array([0, 0, 0.5]), 'R': [0.5, 0.5, 0.5],
                    'X': [0.,  0.5, 0.], 'n_k': 30}
        model = {'with_sigma': 'model', 'Sigma_0': [0., 0., 0.], 'Sigma_Z': [0.9, 0.5, 0.3],
                 'mu_dmft': 0.0, 'eta': 0.03}

        def run(**extra):
            return pcb.get_dmft_bands(orbital_order_to=self.orbital_order_to, **self.w90_dict,
                                      **tb_bands, **model, w_mesh={'window': [-1., 1.], 'n_w': 301},
                                      **extra)[1]

        alatt_orb = run(trace=False)     # -1/pi Im G_orb,orb in the Wannier basis
        alatt_trace = run()

        projected = np.stack([run(proj_on_orb=[orb]) for orb in range(3)], axis=-1)
        assert np.allclose(projected, alatt_orb), 'orbital projection is not the orbital-resolved A'
        assert np.allclose(run(proj_on_orb=[0, 1]), alatt_orb[:, :, [0, 1]].sum(-1))

        # the band-character weighting is still reachable and still sums to the trace,
        # but differs substantially for an orbital dependent self-energy
        band_char = np.stack([run(proj_on_orb=[orb], band_basis=True) for orb in range(3)], axis=-1)
        assert np.allclose(band_char.sum(-1), alatt_trace)
        assert np.allclose(projected.sum(-1), alatt_trace)
        assert np.max(np.abs(band_char - alatt_orb)) > 1.0, 'expected the two conventions to differ here'

        # same for the Fermi slice, which goes through _calc_kslice. A static
        # orbital dependent shift is needed here, because the Fermi liquid
        # self-energy vanishes at w=0 when Sigma_0 is zero.
        kslice = {'bands_path': [('Y', 'G'), ('G', 'X')], 'Y': np.array([0.5, 0., 0.]),
                  'G': [0., 0., 0.], 'X': [0., 0.5, 0.], 'n_k': 20}
        slice_model = dict(model, Sigma_0=[0.4, -0.3, 0.1], eta=0.05)

        def run_slice(**extra):
            return pcb.get_dmft_bands(fermi_slice=True, orbital_order_to=self.orbital_order_to,
                                      **self.w90_dict, **kslice, **slice_model,
                                      w_mesh={'window': [-0.5, 0.5], 'n_w': 201}, **extra)[1]

        slice_total = run_slice()
        slice_proj = [run_slice(proj_on_orb=[orb]) for orb in range(3)]
        slice_band = [run_slice(proj_on_orb=[orb], band_basis=True) for orb in range(3)]
        assert np.allclose(sum(slice_proj), slice_total)
        assert np.allclose(sum(slice_band), slice_total)
        assert max(np.max(np.abs(a - b)) for a, b in zip(slice_proj, slice_band)) > 0.1

        # fat-band data for plot_bands must still come back in the band basis
        tb_data, _, _ = pcb.get_dmft_bands(orbital_order_to=self.orbital_order_to, **self.w90_dict,
                                           **tb_bands, **model, w_mesh={'window': [-1., 1.], 'n_w': 301},
                                           proj_on_orb=[0, 1])
        e_mat = tb_data['e_mat']
        assert tb_data['e_vecs'] is not None
        offdiag = e_mat - np.stack([np.diag(np.diag(e_mat[:, :, ik])) for ik in range(e_mat.shape[2])], axis=-1)
        assert np.allclose(offdiag, 0.0), 'returned e_mat is not diagonal, plot_bands needs eigenvalues'

    def test_bloch_basis_gives_identical_bands(self):
        """
        The self-energy is written into the Wannier basis of seed_hr.dat, which
        is a k-independent index operation. Whether the archive stores H_W(k)
        (bloch_basis=False) or the KS eigenvalues plus k-dependent projectors
        (bloch_basis=True) only fixes the gauge the archive is written in. As
        long as the projectors are square and unitary the two are related by a
        unitary at every k, so A(k,w) must come out identical no matter how
        strongly proj_mat varies along k.
        """

        _write_bloch_basis_archive('./svo_example.h5', './svo_bloch.h5')

        with HDFArchive('./svo_bloch.h5', 'r') as ar:
            proj_mat = ar['dft_input']['proj_mat'][:, 0, 0]
        # the projectors really are k-dependent and unitary
        assert not np.allclose(proj_mat[0], proj_mat[-1])
        assert np.mean(np.abs(np.diff(proj_mat, axis=0))) > 0.1
        assert np.allclose(np.einsum('kab,kcb->kac', proj_mat, proj_mat.conj()), np.eye(3))

        tb_bands = {'bands_path': [('R', 'G'), ('G', 'X'), ('X', 'M'), ('M', 'G')], 'G': [0., 0., 0.],
                    'Z': np.array([0, 0, 0.5]), 'M': [0.5, 0.5, 0.], 'R': [0.5, 0.5, 0.5],
                    'X': [0.,  0.5, 0.], 'n_k': 50}

        sigma_dict = dict(self.sigma_dict, dmft_path='./svo_bloch.h5')
        tb_data, alatt_k_w, _ = pcb.get_dmft_bands(with_sigma='calc', add_mu_tb=True,
                                                   orbital_order_to=self.orbital_order_to,
                                                   **self.w90_dict, **tb_bands, **sigma_dict)

        with HDFArchive('test_pcb_ref.h5', 'r') as ar:
            emat_ref = ar['tb_emat']
            Akw_ref = ar['Akw']

        assert np.allclose(tb_data['e_mat'], emat_ref)
        assert np.allclose(alatt_k_w, Akw_ref)

    def test_bloch_basis_projected_spectral_function(self):
        """
        The orbital projection and the band-resolved output are built from the
        eigenvectors of the tight-binding H(k) and from Sigma in the Wannier
        basis, neither of which depends on the gauge the archive is written in.
        Unlike the trace, an orbital-resolved quantity is not gauge invariant in
        general, so this is checked explicitly rather than inferred.
        """

        _write_bloch_basis_archive('./svo_example.h5', './svo_bloch.h5')

        tb_bands = {'bands_path': [('R', 'G'), ('G', 'X')], 'G': [0., 0., 0.],
                    'Z': np.array([0, 0, 0.5]), 'R': [0.5, 0.5, 0.5],
                    'X': [0.,  0.5, 0.], 'n_k': 30}

        for extra in ({'proj_on_orb': [0, 1]}, {'proj_on_orb': [2]},
                      {'trace': False}, {'band_basis': True}):
            results = []
            for archive in ('./svo_example.h5', './svo_bloch.h5'):
                sigma_dict = dict(self.sigma_dict, dmft_path=archive)
                results.append(pcb.get_dmft_bands(with_sigma='calc', add_mu_tb=True,
                                                  orbital_order_to=self.orbital_order_to,
                                                  **self.w90_dict, **tb_bands, **sigma_dict,
                                                  **extra))
            (tb_wan, alatt_wan, _), (tb_blo, alatt_blo, _) = results
            assert np.allclose(alatt_blo, alatt_wan), f'A(k,w) differs for {extra}'
            assert np.allclose(tb_blo['e_mat'], tb_wan['e_mat']), f'e_mat differs for {extra}'

    def test_bloch_basis_disentanglement_warns(self):
        """
        With more bands than Wannier functions the projectors are an isometry,
        A(k,w) is the spectral function of the Wannier model rather than of the
        DMFT lattice problem and mu was converged for a different electron
        count. That is a caveat for the user, not a reason to refuse.
        """

        _write_bloch_basis_archive('./svo_example.h5', './svo_bloch_dis.h5', n_extra_bands=2)

        tb_bands = {'bands_path': [('G', 'X')], 'G': [0., 0., 0.], 'X': [0., 0.5, 0.], 'n_k': 10}
        sigma_dict = dict(self.sigma_dict, dmft_path='./svo_bloch_dis.h5')

        with warnings.catch_warnings(record=True) as raised:
            warnings.simplefilter('always')
            pcb.get_dmft_bands(with_sigma='calc', add_mu_tb=True,
                               orbital_order_to=self.orbital_order_to,
                               **self.w90_dict, **tb_bands, **sigma_dict)

        messages = [str(w.message) for w in raised]
        assert any('disentanglement' in m for m in messages), messages

    def test_sigma_embedding_default_and_explicit(self):

        hopping = np.zeros((4, 1, 8, 8), dtype=complex)
        # two well separated shells: a 5-fold one around -2 eV, a 3-fold one around +3 eV
        hopping[:, 0, :5, :5] = np.diag([-2.4, -2.2, -2.0, -1.8, -1.6])
        hopping[:, 0, 5:, 5:] = np.diag([2.6, 3.0, 3.4])
        proj_mat = np.zeros((4, 1, 2, 5, 8), dtype=complex)
        proj_mat[:, 0, 0, :5, :5] = np.eye(5)
        proj_mat[:, 0, 1, :3, 5:] = np.eye(3)
        sum_k = _FakeSumk([5, 3], hopping, proj_mat)
        tb_hloc = np.einsum('k,kab->ab', np.full(4, 0.25), hopping[:, 0])

        # default follows the wannier90 convention: correlated shells come first, in order
        embedding = pcb._get_sigma_embedding(sum_k, 8)
        assert [idx.tolist() for idx in embedding] == [[0, 1, 2, 3, 4], [5, 6, 7]]

        with warnings.catch_warnings(record=True) as raised:
            warnings.simplefilter('always')
            pcb._check_sigma_embedding(sum_k, 8, embedding, tb_hloc=tb_hloc)
        assert len(raised) == 0, [str(w.message) for w in raised]

        # an explicit mapping that swaps the two shells must be caught by
        # comparing the local Hamiltonians
        swapped = pcb._get_sigma_embedding(sum_k, 8, [[3, 4, 5, 6, 7], [0, 1, 2]])
        with warnings.catch_warnings(record=True) as raised:
            warnings.simplefilter('always')
            pcb._check_sigma_embedding(sum_k, 8, swapped, tb_hloc=tb_hloc, explicit_embedding=True)
        assert any('Local Hamiltonian' in str(w.message) for w in raised)

        # malformed mappings are rejected outright
        with self.assertRaises(AssertionError):
            pcb._get_sigma_embedding(sum_k, 8, [[0, 1], [5, 6, 7]])
        with self.assertRaises(AssertionError):
            pcb._get_sigma_embedding(sum_k, 8, [[0, 1, 2, 3, 4], [4, 5, 6]])
        with self.assertRaises(AssertionError):
            pcb._get_sigma_embedding(sum_k, 6)

    def test_non_w90_archive_warns_about_orbital_order(self):

        hopping = np.zeros((4, 1, 3, 3), dtype=complex)
        proj_mat = np.zeros((4, 1, 1, 3, 3), dtype=complex)
        proj_mat[:, 0, 0] = np.eye(3)
        sum_k = _FakeSumk([3], hopping, proj_mat, dft_code='vasp')

        with warnings.catch_warnings(record=True) as raised:
            warnings.simplefilter('always')
            pcb._check_sigma_embedding(sum_k, 3, pcb._get_sigma_embedding(sum_k, 3))
        assert any('vasp' in str(w.message) for w in raised), [str(w.message) for w in raised]


if __name__ == '__main__':
    unittest.main()
