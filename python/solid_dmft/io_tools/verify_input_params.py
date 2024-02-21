from typing import Any, Dict, List, Union
import numpy as np
import triqs.utility.mpi as mpi

ParamDict = Dict[str, Any]
FullConfig = Dict[str, Union[ParamDict, List[ParamDict]]]

def verify_input_params_general(params: FullConfig) -> None:
    # warning if sigma mixing is used, remove in future versions
    if params['general']['sigma_mix'] < 1.0 and params['general']['g0_mix'] < 1.0:
        raise ValueError('You shall not use Sigma and G0 mixing together!')

    if params['general']['calc_energies'] and any(entry['type'] == 'ftps' for entry in params['solver']):
        raise ValueError('"calc_energies" is not valid for solver of type = "ftps"')

    if (params['general']['dc'] and params['general']['dc_type'] == 4
            and not np.isclose(sum(params['general']['cpa_x']), 1)):
        raise ValueError('Probability distribution for CPA must equal 1.')

def verify_input_params_dft(params: FullConfig) -> None:
    pass

def verify_input_params_solver(params: FullConfig) -> None:
    pass

def verify_input_params_advanced(params: FullConfig) -> None:
    pass

def manual_changes_input_params(params: FullConfig) -> None:
    """ Necessary workarounds for some of the parameters. """

    # Makes sure that pick_solver_struct and map_solver_struct are a list of dict
    if isinstance(params['advanced']['pick_solver_struct'], dict):
        params['advanced']['pick_solver_struct'] = [params['advanced']['pick_solver_struct']]
    if isinstance(params['advanced']['map_solver_struct'], dict):
        params['advanced']['map_solver_struct'] = [params['advanced']['map_solver_struct']]

    for entry in params['solver']:
        # Calculates the number of solver cycles per rank
        if entry['type'] in ('cthyb', 'ctint', 'ctseg'):
            entry['n_cycles'] = entry['n_cycles_tot'] // mpi.size
            del entry['n_cycles_tot']

        # Some parameters have different names for ctseg
        if entry['type'] == 'ctseg':
            entry['measure_gt'] = entry['measure_G_tau']
            del entry['measure_G_tau']

            entry['measure_gw'] = entry['measure_G_iw']
            del entry['measure_G_iw']

            # Makes sure measure_gw is true if improved estimators are used
            if entry['improved_estimator']:
                entry['measure_gt'] = True
                entry['measure_ft'] = True
            else:
                entry['measure_ft'] = False
            del entry['improved_estimator']

            entry['measure_gl'] = entry['measure_G_l']
            del entry['measure_G_l']

            entry['measure_hist'] = entry['measure_pert_order']
            del entry['measure_pert_order']

        # use_norm_as_weight also required to measure the density matrix
        if entry['type'] == 'cthyb' and entry['measure_density_matrix']:
            entry['use_norm_as_weight'] = True


    return

    # TODO: treat the following parameters in the solver.py class?
    if params['general']['solver_type'] in ['cthyb']:
        params['general']['cthyb_delta_interface'] = params['solver']['delta_interface']
        del params['solver']['delta_interface']

    # little workaround since #leg coefficients is not directly a solver parameter
    if 'legendre_fit' in params['solver']:
        params['general']['legendre_fit'] = params['solver']['legendre_fit']
        del params['solver']['legendre_fit']

    params['general']['store_solver'] = params['solver']['store_solver']
    del params['solver']['store_solver']