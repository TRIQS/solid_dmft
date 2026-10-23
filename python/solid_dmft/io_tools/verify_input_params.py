from typing import Any, Dict, List, Union
import warnings
import numpy as np

ParamDict = Dict[str, Any]
FullConfig = Dict[str, Union[ParamDict, List[ParamDict]]]


def _resolve_mixing_params(params: FullConfig) -> None:
    """
    Maps the deprecated keys g0_mix, g0_mix_type and sigma_mix to mix_quantity,
    mix_type and mix_greed, and fills in the defaults mix_type = 'linear' and
    mix_greed = 1.0. Modifies params in place.
    """
    general = params['general']
    old_keys = [key for key in ('g0_mix', 'g0_mix_type', 'sigma_mix') if general[key] is not None]
    new_keys = [key for key in ('mix_quantity', 'mix_type', 'mix_greed') if general[key] is not None]

    if old_keys:
        if new_keys:
            raise ValueError(f'Mixing is specified both with the deprecated keys {old_keys} and with {new_keys}. '
                             'Use only "mix_quantity", "mix_type" and "mix_greed".')
        g0_active = general['g0_mix_type'] == 'broyden' or (general['g0_mix'] is not None and general['g0_mix'] < 1.0)
        sigma_active = general['sigma_mix'] is not None and general['sigma_mix'] < 1.0
        if g0_active and sigma_active:
            raise ValueError('You shall not use Sigma and G0 mixing together!')
        if g0_active:
            general['mix_quantity'] = 'G0'
            general['mix_type'] = general['g0_mix_type'] if general['g0_mix_type'] is not None else 'linear'
            general['mix_greed'] = general['g0_mix'] if general['g0_mix'] is not None else 1.0
        elif sigma_active:
            general['mix_quantity'] = 'Sigma'
            general['mix_type'] = 'linear'
            general['mix_greed'] = general['sigma_mix']
        replacement = ', '.join(f'{key} = {general[key]!r}' for key in ('mix_quantity', 'mix_type', 'mix_greed')
                                if general[key] is not None) or 'nothing (no mixing)'
        warnings.warn(f'The general parameters {old_keys} are deprecated and will be removed. '
                      f'They are interpreted as {replacement}.', FutureWarning)
        for key in ('g0_mix', 'g0_mix_type', 'sigma_mix'):
            general[key] = None

    if general['mix_type'] is None:
        general['mix_type'] = 'linear'
    if general['mix_greed'] is None:
        general['mix_greed'] = 1.0


def _verify_mixing_params(params: FullConfig) -> None:
    general = params['general']
    if general['mix_quantity'] not in (None, 'G0', 'Sigma'):
        raise ValueError(f'Invalid "mix_quantity" = {general["mix_quantity"]}, allowed are "G0" and "Sigma".')
    if general['mix_type'] not in ('linear', 'broyden'):
        raise ValueError(f'Invalid "mix_type" = {general["mix_type"]}, allowed are "linear" and "broyden".')
    if general['mix_type'] == 'linear' and not 0.0 <= general['mix_greed'] <= 1.0:
        raise ValueError(f'"mix_greed" = {general["mix_greed"]} must be in [0, 1] for linear mixing.')
    if general['mix_type'] == 'broyden' and not general['mix_greed'] > 0.0:
        raise ValueError(f'"mix_greed" = {general["mix_greed"]} must be positive for Broyden mixing.')
    if general['mix_quantity'] is None and (general['mix_type'] == 'broyden' or general['mix_greed'] < 1.0):
        raise ValueError('"mix_type" or "mix_greed" is set, but "mix_quantity" is not. '
                         'Set "mix_quantity" to "G0" or "Sigma".')
    if general['broy_max_it'] != -1 and general['broy_max_it'] < 1:
        raise ValueError(f'"broy_max_it" = {general["broy_max_it"]} must be -1 or at least 1.')
    if general['gw_embedding'] and (general['mix_quantity'] == 'G0' or general['mix_type'] == 'broyden'):
        raise ValueError('GW embedding supports only linear mixing of Sigma.')


def _verify_input_params_general(params: FullConfig) -> None:
    # Checks that grid parameters are specified completely
    if (
        params['general']['beta'] is not None
        and not params['general']['gw_embedding']
        and (params['general']['n_iw'] is None or params['general']['n_tau'] is None)
    ):
        raise ValueError('Imaginary-frequency grid chosen, specify "n_iw" and "n_tau".')

    if (
        params['general']['beta'] is None
        and not params['general']['gw_embedding']
        and (params['general']['eta'] is None or params['general']['n_w'] is None or params['general']['w_range'] is None)
    ):
        raise ValueError('Real-frequency grid chosen, specify "eta", "n_w", and "w_range".')

    if params['general']['calc_energies'] and any(entry['type'] == 'ftps' for entry in params['solver']):
        raise ValueError('"calc_energies" is not valid for solver of type = "ftps"')

    if params['general']['n_iter_dmft'] is None:
        raise ValueError('"n_iter_dmft" must be specified.')
    elif params['general']['n_iter_dmft'] < 0:
        raise ValueError('"n_iter_dmft" must be at least 0.')

    # Checks validity of other general params
    h_int_type_options = (
        'density_density',
        'kanamori',
        'kanamori_den_den',
        'full_slater',
        'crpa',
        'crpa_density_density',
        'dyn_density_density',
        'dyn_full',
        'ntot',
        'simple_intra',
    )
    if isinstance(params['general']['h_int_type'], str):
        if params['general']['h_int_type'] not in h_int_type_options:
            raise ValueError(f'Invalid "h_int_type" = {params["general"]["h_int_type"]}.')
    elif isinstance(params['general']['h_int_type'], list):
        if any(entry not in h_int_type_options for entry in params['general']['h_int_type']):
            raise ValueError('Invalid "h_int_type" in input list.')
    else:
        raise ValueError('Invalid "h_int_type" input type. String or list expected.')

    _verify_mixing_params(params)

    if params['general']['calc_mu_method'] not in ('dichotomy', 'newton', 'brent'):
        raise ValueError(f'Invalid "calc_mu_method" = {params["general"]["calc_mu_method"]}.')

    if params['general']['set_rot'] not in (None, 'den', 'hloc'):
        raise ValueError(f'Invalid "set_rot" = {params["general"]["set_rot"]}.')

    if params['general']['h_int_basis'] not in ('triqs', 'wien2k', 'wannier90', 'qe', 'vasp'):
        raise ValueError(f'Invalid "h_int_basis" = {params["general"]["h_int_basis"]}.')


def _verify_input_params_dft(params: FullConfig) -> None:
    if params['dft']['dft_code'] not in ('vasp', 'qe', None):
        raise ValueError(f'Invalid "dft.dft_code" = {params["dft"]["dft_code"]}.')

    if params['dft']['mpi_env'] not in ('default', 'openmpi', 'openmpi-intra', 'mpich', 'slurm'):
        raise ValueError(f'Invalid "dft.mpi_env" = {params["dft"]["mpi_env"]}.')

    if params['dft']['projector_type'] not in ('w90', 'plo'):
        raise ValueError(f'Invalid "dft.projector_type" = {params["dft"]["projector_type"]}.')


def _verify_input_params_solver(params: FullConfig) -> None:
    solver_params = params['solver']

    # Checks that the solver impurities index are all lists if there are multiple solvers
    if len(solver_params) > 1 or solver_params[0]['idx_impurities'] is not None:
        if any(not isinstance(entry['idx_impurities'], list) for entry in solver_params):
            raise ValueError(
                'All "solver.idx_impurities" need to specify the impurities ' 'as a list of ints when there are multiple solvers.'
            )
        for entry in solver_params:
            if any(not isinstance(imp, int) for imp in entry['idx_impurities']):
                raise ValueError(
                    'All "solver.idx_impurities" need to specify the impurities ' 'as a list of ints when there are multiple solvers.'
                )

    # Checks that all solvers support the specified grid
    # TODO: add real-frequency support for solvers that do both (e.g., hartree)
    supported_grids = {'real': ['ftps'], 'imag': ['cthyb', 'ctint', 'ctseg', 'hubbardI', 'hartree']}
    if params['general']['beta'] is not None:
        for entry in solver_params:
            if entry['type'] not in supported_grids['imag']:
                raise ValueError(f'Solver {entry["type"]} does not support real-frequency grid.')
    elif not params['general']['gw_embedding']:
        for entry in solver_params:
            if entry['type'] not in supported_grids['real']:
                raise ValueError(f'Solver {entry["type"]} does not support imaginary-frequency grid.')

    for entry in solver_params:
        if entry['type'] == 'cthyb' and entry['crm_dyson_solver'] and not entry['measure_density_matrix']:
            raise ValueError(
                'Solver "cthyb" when using "crm_dyson_solver" must also measure the density matrix: "measure_density_matrix" = True'
            )
        if entry['type'] == 'ctseg':
            tail_op = [entry['crm_dyson_solver'],
                       entry['legendre_fit'],
                       entry['improved_estimator'],
                       entry['perform_tail_fit']]
            if sum(tail_op) > 1:
                # allow improved_estimator + tail_fit
                if sum([entry['improved_estimator'], entry['perform_tail_fit']]) < 2:
                    raise ValueError(
                        'Except "improved_estimator" + "perform_tail_fit", '
                        'only one of the options "crm_dyson_solver", "legendre_fit", '
                        '"improved_estimator", and "perform_tail_fit" can be set to True.'
                    )
        if entry['type'] == 'cthyb':
            tail_op = [entry['crm_dyson_solver'],
                       entry['legendre_fit'],
                       entry['measure_G_l'],
                       entry['perform_tail_fit']]
            if sum(tail_op) > 1:
                raise ValueError('Only one of the options "crm_dyson_solver", "legendre_fit", "measure_G_l", or "perform_tail_fit" can be set to True.')


def _verify_input_params_gw(params: FullConfig) -> None:
    pass

def _verify_input_params_advanced(params: FullConfig) -> None:
    pass


def verify_before_dmft_cycle(params):
    """Maps deprecated keys and runs all checks of params that can be done before dmft_cycle is called."""
    _resolve_mixing_params(params)
    _verify_input_params_general(params)
    _verify_input_params_dft(params)
    _verify_input_params_solver(params)
    _verify_input_params_gw(params)
    _verify_input_params_advanced(params)


def verify_h5_dependent(sum_k, solver_type_per_imp, general_params):
    """All checks of params that depend on the h5 file content that is stored in the SumkDFT object."""
    # Incompatabilities for SO coupling
    if sum_k.SO == 1 and general_params['magnetic'] and general_params['afm_order']:
        raise ValueError('AFM order not supported with SO coupling')

    # Checks that enforce_off_diag true for ftps and hartree
    if any(s in ['ftps', 'hartree'] and not e for s, e in zip(solver_type_per_imp, general_params['enforce_off_diag'])):
        raise ValueError('enforce_off_diag must be True for a impurities solver by ftps or hartree solvers')

    # Checks that the interaction Hamiltonian and the parameters match
    if any(
        h not in ['density_density', 'slater'] and r is not None
        for h, r in zip(general_params['h_int_type'], general_params['ratio_F4_F2'])
    ):
        raise ValueError(
            '"ratio_F4_F2" only considered for interaction Hamiltonians "density_density" and "slater". '
            'Please set to None for all other Hamiltonians.'
        )
    if any(h not in ('kanamori', 'kanmori_den_den') and up is not None for h, up in zip(general_params['h_int_type'], general_params['U_prime'])):
        raise ValueError(
            '"U_prime" only considered for interaction Hamiltonian "kanamori". ' 'Please set to None for all other Hamiltonians.'
        )
