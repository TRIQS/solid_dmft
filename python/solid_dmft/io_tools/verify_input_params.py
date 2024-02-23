from typing import Any, Dict, List, Union
import numpy as np
import triqs.utility.mpi as mpi

ParamDict = Dict[str, Any]
FullConfig = Dict[str, Union[ParamDict, List[ParamDict]]]

def _verify_input_params_general(params: FullConfig) -> None:
    # Checks that grid parameters are specified completely
    if params['general']['beta'] is not None and (params['general']['n_iw'] is None
                                                  or params['general']['n_tau'] is None):
        raise ValueError('Imaginary-frequency grid chosen, specify "n_iw" and "n_tau".')

    if params['general']['beta'] is None and (params['general']['eta'] is None or params['general']['n_w'] is None
                                              or params['general']['w_range'] is None):
        raise ValueError('Real-frequency grid chosen, specify "eta", "n_w", and "w_range".')

    # warning if sigma mixing is used, remove in future versions
    if params['general']['sigma_mix'] < 1.0 and params['general']['g0_mix'] < 1.0:
        raise ValueError('You shall not use Sigma and G0 mixing together!')

    if params['general']['calc_energies'] and any(entry['type'] == 'ftps' for entry in params['solver']):
        raise ValueError('"calc_energies" is not valid for solver of type = "ftps"')

    if (params['general']['dc'] and params['general']['dc_type'] == 4
            and not np.isclose(sum(params['general']['cpa_x']), 1)):
        raise ValueError('Probability distribution for CPA must equal 1.')

def _verify_input_params_dft(params: FullConfig) -> None:
    pass

def _verify_input_params_solver(params: FullConfig) -> None:
    solver_params = params['solver']

    # Checks that the solver impurities index are all lists if there are multiple solvers
    if len(solver_params) > 1 or solver_params[0]['idx_impurities'] is not None:
        if any(not isinstance(entry['idx_impurities'], list) for entry in solver_params):
            raise ValueError('All "solver.idx_impurities" need to specify the impurities '
                             'as a list of ints when there are multiple solvers.')
        for entry in solver_params:
            if any(not isinstance(imp, int) for imp in entry['idx_impurities']):
                raise ValueError('All "solver.idx_impurities" need to specify the impurities '
                                 'as a list of ints when there are multiple solvers.')

    # Checks that all solvers support the specified grid
    # TODO: add real-frequency support for solvers that do both (e.g., hartree)
    supported_grids = {'real': ['ftps'],
                       'imag': ['cthyb', 'ctint', 'ctseg', 'hubbardI', 'hartree']}
    if params['general']['beta'] is not None:
        for entry in solver_params:
            if entry['type'] not in supported_grids['imag']:
                raise ValueError(f'Solver {entry["type"]} does not support real-frequency grid.')
    else:
        for entry in solver_params:
            if entry['type'] not in supported_grids['real']:
                raise ValueError(f'Solver {entry["type"]} does not support imaginary-frequency grid.')

def _verify_input_params_advanced(params: FullConfig) -> None:
    pass

def verify_input_params(params: FullConfig) -> None:
    _verify_input_params_general(params)
    _verify_input_params_dft(params)
    _verify_input_params_solver(params)
    _verify_input_params_advanced(params)

def manual_changes_input_params(params: FullConfig) -> None:
    """ Necessary workarounds for some of the parameters. """

    return
