#!/usr/bin/env python3
# -*- coding: utf-8 -*-
################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2016-2018, N. Wentzell
# Copyright (C) 2018-2019, Simons Foundation
#   author: N. Wentzell
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
"""
Analytic continuation of the self energy using maxent.

Reads Sigma(i omega) from the h5 archive and writes Sigma(omega) back.

Analytic continuation of selfenergies requires using an auxiliary Green's
function. This makes the procedure even less exact than it already is so use
the results with care!

Author: Max Merkel, 2020
"""

import time
import sys
import glob
from multiprocessing import Pool
from functools import partial
import numpy as np

from triqs_maxent.sigma_continuator import InversionSigmaContinuator, DirectSigmaContinuator
from triqs_maxent.elementwise_maxent import PoormanMaxEnt
from triqs_maxent.omega_meshes import HyperbolicOmegaMesh
from triqs_maxent.alpha_meshes import LogAlphaMesh
from triqs.gf import BlockGf
from h5 import HDFArchive


def _read_h5(external_path, iteration=None):
    """
    Reads the h5 archive to get the Matsubara self energy, the double counting energy
    and the chemical potential.

    Parameters:
    -----------
    external_path: string, path to h5 archive
    iteration: int, the iteration that is being read from, default is 'last_iter'

    Returns:
    --------
    sigma_iw: list, self energy as block Green's function for each impurity
    dc_energy: list, double counting for each impurity
    chemical_potential: complex, the chemical potential of the problem. Should be approximately real
    """

    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None
                                          else 'it_{}'.format(iteration))

    with HDFArchive(external_path, 'r') as archive:
        impurity_paths = [key for key in archive[h5_internal_path].keys() if 'Sigma_freq_' in key]
        # Sorts impurity paths by their indices, not sure if necessary
        impurity_indices = [int(s[s.rfind('_')+1:]) for s in impurity_paths]
        impurity_paths = [impurity_paths[i] for i in np.argsort(impurity_indices)]
        sigma_iw = [archive[h5_internal_path][p] for p in impurity_paths]

        dc_energy = archive[h5_internal_path]['DC_energ']
        if 'chemical_potential_post' in archive[h5_internal_path]:
            chemical_potential = archive[h5_internal_path]['chemical_potential_post']
        else:
            # Old name for chemical_potential_post
            chemical_potential = archive[h5_internal_path]['chemical_potential']

    return sigma_iw, dc_energy, chemical_potential


def _create_sigma_continuators(sigma_iw, dc_energy, chemical_potential):
    """
    Initializes the inversion and direct sigma continuator from Sigma(i omega) and the DC energy.
    Returns a dict of three versions of continuators, using different auxiliary Green's functions:
    * 'inversion_continuator_dc': inversion continuator, constant C = dc_energy for the impurity
    * 'inversion_continuator_mu': inversion continuator, constant C = Sigma(i infinity) + chemical potential
          (NOT TESTED YET)
    * 'direct_continuator': direct continuator, not implemented because it makes
    for non-diagonal Sigma crash even if not used
    """
    inversion_continuator_dc = [{key: InversionSigmaContinuator(block_GF, imp_dc)
                                 for key, block_GF in imp_GF}
                                for imp_GF, imp_dc in zip(sigma_iw, dc_energy)]
    # TODO: implement inversion sigma continuator using chemical potential and Sigma(i infty) as constant
    # Implementation might or might not be correct
    inversion_continuator_mu = [{key: InversionSigmaContinuator(block_GF, block_GF.data[-1].real
                                                                + chemical_potential)
                                 for key, block_GF in imp_GF} for imp_GF in sigma_iw]

    return {'inversion_continuator_dc': inversion_continuator_dc,
            'inversion_continuator_mu': inversion_continuator_mu}


def _run_maxent(sigma_continuator, error=.03):
    """
    Uses maxent to continue the auxiliary Green's function obtained from the continuator.
    The default range for omega is -20 to +20, for alpha 1e-4 to 1e2. This can be changed in here.
    Returns the real-frequency auxiliary Green's function
    """
    results = [{} for _ in range(len(sigma_continuator))]
    for i in range(len(sigma_continuator)):
        print('-'*50 + '\nSolving impurity {}/{}\n'.format(i+1, len(sigma_continuator)) + '-'*50)
        print('Found blocks {}'.format(list(sigma_continuator[i].keys())))
        for key, continuator in sigma_continuator[i].items():
            solver = PoormanMaxEnt(use_complex=True)
            solver.set_G_iw(continuator.Gaux_iw)
            solver.set_error(error)
            solver.omega = HyperbolicOmegaMesh(omega_min=-12, omega_max=12, n_points=400)
            solver.alpha_mesh = LogAlphaMesh(alpha_min=1e-6, alpha_max=1e2, n_points=50)
            results[i][key] = solver.run()
            print(results[i][key].analyzer_results[0][0][0]['LineFitAnalyzer']['info'])
    return results


def _get_sigma_omega_from_aux_spectral(sigma_continuator, results):
    """ Extracts the real-frequency self energy from the auxiliary Green's function. """
    for imp_sc, imp_res in zip(sigma_continuator, results):
        for block_sc, block_res in zip(imp_sc.values(), imp_res.values()):
            # Uses that matrix is Hermitian
            data = block_res.get_A_out('LineFitAnalyzer')[:]
            for i in range(data.shape[0]):
                for j in range(i):
                    data[i, j] = data[j, i].conjugate()
            block_sc.set_Gaux_w_from_Aaux_w(data, block_res.omega, np_interp_A=2000,
                                            np_omega=5000, w_min=-12, w_max=12)

    sigma_w = [BlockGf(name_list=continuator.keys(),
                       block_list=[c.S_w for c in continuator.values()])
               for continuator in sigma_continuator]
    return sigma_w


def _write_sigma_omega_to_h5(sigma_w, results, external_path, iteration=None):
    """ Writes real-frequency self energy to h5 archive. """
    h5_internal_path = 'DMFT_results/' + ('last_iter' if iteration is None
                                          else 'it_{}'.format(iteration))

    with HDFArchive(external_path, 'a') as archive:
        for i, sigma in enumerate(sigma_w):
            impurity_path = 'Sigma_w_{}'.format(i)
            archive[h5_internal_path][impurity_path] = sigma
            # result_path = 'Sigma_w_ME_res_{}'.format(i)
            # res = {}
            # for block_name, block_res in results[i].iteritems():
                # res[block_name] = block_res.data

            # archive[h5_internal_path][result_path] = res


def main(external_path, iteration=None):
    """
    Main function that reads the Matsubara self-energy from h5, analytically continues it
    and writes the result back to the h5 archive.

    Parameters
    ----------
    external_path: string, path of the h5 archive
    iteration: int/string, optional, iteration to read from and write to

    Returns
    -------
    list of triqs.gf.GfReFreq, per impurity: Sigma(omega) as a GF object
    """
    start_time = time.time()

    sigma_iw, dc_energy, chemical_potential = _read_h5(external_path, iteration)
    print('Finished reading h5 archive. Found {} impurities.'.format(len(sigma_iw)))
    sigma_continuator = _create_sigma_continuators(sigma_iw, dc_energy,
                                                   chemical_potential)['inversion_continuator_dc']
    print('Starting run of maxent now.')
    results = _run_maxent(sigma_continuator, error=0.02)
    print(u'Extracting Σ(ω) now. This might take a while.')
    sigma_w = _get_sigma_omega_from_aux_spectral(sigma_continuator, results)
    print('Writing results to h5 archive now.')
    _write_sigma_omega_to_h5(sigma_w, results, external_path, iteration)
    print('Finished writing Σ(ω) to archive.')

    total_time = time.time() - start_time
    print('Run time: {:.0f} s.'.format(total_time))
    return sigma_w


if __name__ == '__main__':
    if len(sys.argv) not in (2, 3):
        print('Please give the h5 name (and optionally the iteration). Exiting.')
        sys.exit(2)

    files = glob.glob(sys.argv[1])
    pool = Pool(processes=min(8, len(files)))

    if len(sys.argv) == 2:
        function = main
    elif len(sys.argv) == 3:
        function = partial(main, iteration=sys.argv[2])

    pool.map(function, files)
