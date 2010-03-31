import numpy as np

from gensei.base import ordered_iteritems

def get_random(ranges):
    """Get an array of random numbers 0 <= a_i < ranges[i]."""
    ranges = np.atleast_1d(ranges)
    return ranges * np.random.random(len(ranges))

def get_suffix(n):
    """Get suffix format string given a number of files.
    
    Examples:

        n = 5 -> '%01d'
        n = 15 -> '%02d'
        n = 1005 -> '%04d'
    """
    if n > 1:
        n_digit = int(np.log10(n - 1) + 1)
        suffix = '%%0%dd' % n_digit
    else:
        suffix = '%d'
    return suffix

def format_dict(d, raw=None, indent=2):
    """Format a dictionary for printing.

    Parameters:

        d : dict
            The dictionary.

        raw : dict
            The raw (unadjusted) dictionary to compare with.

        indent : int
            The indentation level.

    Return:

        msg : string
           The string with dictionary's key : val formatted in two columns.
    """
    if raw is None:
        raw = d
        
    msg = ''
    for key, val in ordered_iteritems(d):
        if val == raw[key]:
            msg += (' ' * indent) + ('%s : %s\n' % (key, val))
        else:
            msg += (' ' * indent) + ('%s : %s (%s)\n' % (key, val, raw[key]))

    msg = msg[:-1]

    return msg
