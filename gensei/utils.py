import numpy as np

from gensei.base import is_sequence, ordered_iteritems

def get_random(ranges0, ranges1):
    """
    Get an array of random numbers `a` with the uniform distribution so that
    `ranges0[i] <= a[i] < ranges1[i]`.

    Returns
    -------
    a : array
        The array of random numbers.
    """
    ranges0 = np.atleast_1d(ranges0)
    return (ranges1 - ranges0) * np.random.random(len(ranges0)) + ranges0

def get_random_direction(dim):
    """
    Get a uniform random direction in unit ball of dimension `dim`. 
    """
    while 1:
        coor = np.random.normal(0.0, 1.0, size=dim)
        r = np.linalg.norm(coor)
        if (r <= 1.0) and (r > 1e-5) :
            break
    out = coor / r

    return out

def evaluate(val, shape=None):
    """
    Evaluate a single value or a sequence of values.

    The values can be have a uniform (keyword 'random') or normal (keyword
    'normal') random distribution.

    Parameters
    ----------
    val : config value
        The value to be evaluated. It can be one of:

        - list: evaluate() is called recursively for each item
        - tuple:
            - ('random', [`min. values`], [`max. values`]) (uniform)
            - ('random', `min. value`, `max. value`)
            - ('normal', `mean`, `standard deviation`)
            - ('random direction', `dim`) : uniform distribution of unit vectors
              in dimension `dim`
        - 'random' : random value(s) in [0, 1)
        - 'random direction' : uniform distribution of unit vectors
              in dimension 3
        - other values are returned as they are.
    shape : tuple, optional
        Shape of the data when `val` is 'random'.

    Returns
    -------
    out : data
        The evaluated data.
    """
    if isinstance(val, list):
        out = []
        for sub_val in val:
            out.append(evaluate(sub_val))
        out = np.array(out, dtype=np.float64)

    elif isinstance(val, tuple):
        if val[0] == 'normal':
            out = np.random.normal(val[1], val[2])

        elif val[0] == 'random':
            if is_sequence(val[1]):
                out = get_random(val[1], val[2])

            else:
                out = (val[2] - val[1]) * np.random.random() + val[1]

        elif val[0] == 'random direction':
            out = get_random_direction(val[1])

        else:
            raise ValueError('unsupported value type! (%s)' % val[0])

    elif val == 'random':
        out = np.random.random(shape)

    elif val == 'random direction':
        out = get_random_direction(3)
        
    else:
        out = val

    return out


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
