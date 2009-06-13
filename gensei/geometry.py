import numpy as np

def make_axis_rotation_matrix(direction, angle):
    """
    Create a rotation matrix corresponding to the rotation around a general
    axis by a specified angle.
    
    R = dd^T + cos(a) (I - dd^T) + sin(a) skew(d)
    
    Parameters:
    
        angle : float a
        direction : array d
    """
    d = np.array(direction, dtype=np.float64)
    d /= np.linalg.norm(d)
    
    eye = np.eye(3, dtype=np.float64)
    ddt = np.outer(d, d)
    skew = np.array([[    0,  d[2],  -d[1]],
                     [-d[2],     0,  d[0]],
                     [d[1], -d[0],    0]], dtype=np.float64)

    mtx = ddt + np.cos(angle) * (eye - ddt) + np.sin(angle) * skew
    return mtx

def make_rotation_matrix_hc(mtx):
    """Extend a rotation matrix to homogenous coordinates.

    R_hc = [[R, 0],
            [0, 1]]

    Parameters:

        mtx : 3 x 3 array R
    """
    zz = np.zeros((3,), dtype=np.float64)
    mtx = np.r_[np.c_[mtx, zz], [np.r_[zz, 1]]]
    return mtx

def make_translation_matrix_hc(point):
    """Create a matrix whose application corresponds to translation to the
    point in homogenous coordinates.

    T_hc = [[1, 0, 0, x],
            [0, 1, 0, y],
            [0, 0, 1, z],
            [0, 0, 0, 1]]

    Parameters:

        point : array (x, y, z)
    """
    mtx = np.eye(4, dtype=np.float64)
    mtx[:-1,3] = point
    return mtx

def get_average_semiaxes(volume, length_to_width):
    """Get semiaxes of an ellipsoid given its volume and length-to-width
    ratio."""
    b = c = np.power(volume / (4.0/3.0 * np.pi * length_to_width),
                     1.0/3.0)
    a = length_to_width * b
    return a, b, c

def transform_to_pixels(coors, max_coors, resolution):
    """
    Transform real coordinates (in [0, max_coors]) to pixel coordinates, given
    the figure resolution and max. coordinates (block dimensions).
    """
    pc = np.asarray(resolution)[np.newaxis,:] * coors / max_coors[np.newaxis,:]
    return pc.astype(np.int32)
