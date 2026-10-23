#!/usr/bin/env python
# Mixing parameters of the general section: mix_quantity, mix_type, mix_greed,
# the mapping of the deprecated g0_mix, g0_mix_type and sigma_mix to them,
# and the checks on their values.

import os
import tomllib
import unittest
import warnings

import solid_dmft.io_tools as io_tools
from solid_dmft.io_tools.postproc_toml_dict import merge_config_with_default
from solid_dmft.io_tools.verify_input_params import _resolve_mixing_params, _verify_mixing_params

with open(os.path.join(os.path.dirname(io_tools.__file__), 'default.toml'), 'rb') as file:
    DEFAULT = tomllib.load(file)


def params(**general):
    config = {'general': {'seedname': 'inp', 'h_int_type': 'kanamori', **general},
              'solver': {'type': 'hubbardI'}}
    return merge_config_with_default(config, DEFAULT, {'solver': 'type'})


def resolve(**general):
    return resolve_params(params(**general))


def resolve_params(p):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        _resolve_mixing_params(p)
    _verify_mixing_params(p)
    g = p['general']
    deprecated = any(issubclass(w.category, FutureWarning) for w in caught)
    return (g['mix_quantity'], g['mix_type'], g['mix_greed']), deprecated


class TestMixingParams(unittest.TestCase):

    def test_default_is_no_mixing(self):
        self.assertEqual(resolve(), ((None, 'linear', 1.0), False))

    def test_new_keys(self):
        self.assertEqual(resolve(mix_quantity='G0', mix_type='broyden', mix_greed=0.4),
                         (('G0', 'broyden', 0.4), False))
        self.assertEqual(resolve(mix_quantity='Sigma', mix_greed=0.7), (('Sigma', 'linear', 0.7), False))
        self.assertEqual(resolve(mix_quantity='Sigma', mix_type='broyden'), (('Sigma', 'broyden', 1.0), False))
        self.assertEqual(resolve(mix_quantity='G0', mix_greed=0.0), (('G0', 'linear', 0.0), False))

    def test_deprecated_keys(self):
        self.assertEqual(resolve(g0_mix=0.8), (('G0', 'linear', 0.8), True))
        self.assertEqual(resolve(g0_mix_type='broyden'), (('G0', 'broyden', 1.0), True))
        self.assertEqual(resolve(g0_mix_type='broyden', g0_mix=0.5), (('G0', 'broyden', 0.5), True))
        self.assertEqual(resolve(sigma_mix=0.9), (('Sigma', 'linear', 0.9), True))
        self.assertEqual(resolve(sigma_mix=0.9, g0_mix=1.0), (('Sigma', 'linear', 0.9), True))
        # tutorials set both to 1.0, which is no mixing
        self.assertEqual(resolve(sigma_mix=1.0, g0_mix=1.0), ((None, 'linear', 1.0), True))

    def test_deprecated_keys_are_cleared(self):
        p = params(g0_mix=0.8)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            _resolve_mixing_params(p)
        self.assertEqual([p['general'][k] for k in ('g0_mix', 'g0_mix_type', 'sigma_mix')], [None] * 3)

    def test_invalid(self):
        invalid = [
            dict(g0_mix=0.5, mix_greed=0.5),
            dict(g0_mix=0.5, sigma_mix=0.5),
            dict(mix_quantity='Gimp'),
            dict(mix_quantity='G0', mix_type='anderson'),
            dict(mix_greed=0.5),
            dict(mix_type='broyden'),
            dict(mix_quantity='G0', mix_greed=1.5),
            dict(mix_quantity='G0', mix_type='broyden', mix_greed=0.0),
            dict(mix_quantity='G0', mix_type='broyden', broy_max_it=0),
            dict(mix_quantity='Sigma', mix_type='broyden', gw_embedding=True),
            dict(mix_quantity='G0', mix_greed=0.5, gw_embedding=True),
        ]
        for general in invalid:
            p = params(**general)
            with self.subTest(**general), self.assertRaises(ValueError):
                resolve_params(p)


if __name__ == '__main__':
    unittest.main()
