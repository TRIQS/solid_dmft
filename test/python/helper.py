# -*- coding: utf-8 -*-

"""
Contains routine used for checking equality between nested iterables, using
numpy.isclose for float comparison.
"""

import numpy as np
from collections.abc import Iterable

class Dummy():
    """Needed to create test objects resembling e.g. SumkDFT or Solver. """
    pass

def are_iterables_equal(obs1, obs2):
    """
    Compares two iterables.

    Dict: compare all keys and values
    String: equality with ==
    Other iterables: compare all values
    Float: numpy.isclose
    Else: Equality with ==
    """

    if isinstance(obs1, dict):
        if sorted(obs1.keys()) != sorted(obs2.keys()):
            print(set(obs1.keys()) ^ set(obs2.keys()))
            return False

        for key, value in obs2.items():
            if not are_iterables_equal(obs1[key], value):
                print('Key {} contains unequal elements'.format(key))
                return False
        return True

    if isinstance(obs1, str) or isinstance(obs1, type(u'')):
        return obs1 == obs2

    if isinstance(obs1, Iterable):
        if len(obs1) != len(obs2):
            return False

        for value1, value2 in zip(obs1, obs2):
            if not are_iterables_equal(value1, value2):
                return False
        return True

    if isinstance(obs1, float):
        return np.isclose(obs1, obs2, atol=1e-4, rtol=1e-4)

    return obs1 == obs2
