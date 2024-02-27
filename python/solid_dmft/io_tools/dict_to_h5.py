from copy import deepcopy

def _iteratively_replace_none(to_write, replace_from, replace_with):
    """ Limitation: can only replace None with a string, or a string with None. """
    # First two checks needed because comparison to triqs many-body operator fails
    if (isinstance(to_write, str) or to_write is None) and to_write == replace_from:
        return replace_with

    if isinstance(to_write, dict):
        for key, value in to_write.items():
            to_write[key] = _iteratively_replace_none(value, replace_from, replace_with)
    elif isinstance(to_write, list):
        for i, value in enumerate(to_write):
            to_write[i] = _iteratively_replace_none(value, replace_from, replace_with)

    return to_write

def prep_params_for_h5(dict_to_write):
    """ Replace all NoneType with a string 'none' to be able to write to h5. """
    return _iteratively_replace_none(deepcopy(dict_to_write), None, 'none')

# Not sure if the reverse route is actually needed
def prep_params_from_h5(dict_to_read):
    """ Replace all 'none' strings with NoneType to parse the dict coming from h5. """
    return _iteratively_replace_none(deepcopy(dict_to_read), 'none', None)
