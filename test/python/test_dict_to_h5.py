from solid_dmft.io_tools import dict_to_h5

def test_prep_params_for_h5():
    inp = {'a': None, 'b': {'c': None, 'd': 'e'}, 'f': [None, 'g']}
    expected = {'a': 'none', 'b': {'c': 'none', 'd': 'e'}, 'f': ['none', 'g']}

    assert dict_to_h5.prep_params_for_h5(inp) == expected

def test_prep_params_from_h5():
    inp = {'a': 'none', 'b': {'c': 'none', 'd': 'e'}, 'f': ['none', 'g']}
    expected = {'a': None, 'b': {'c': None, 'd': 'e'}, 'f': [None, 'g']}

    assert dict_to_h5.prep_params_from_h5(inp) == expected
