#!/usr/bin/env python
from optparse import OptionParser
import numpy as np
from scipy.linalg import eigh, inv
import pylab as pl

def make_axis_rotation_matrix(direction, angle):
    """
    angle : a
    direction : d
    
    R = dd^T + cos(a) (I - dd^T) + sin(a) skew(d)
    """
    d = np.array(direction, dtype=np.float64)
    d /= np.linalg.norm(d)
    
    eye = np.eye(3, dtype=np.float64)
    ddt = np.outer(d, d)
    skew = np.array([[    0,  d[2],  d[1]],
                     [-d[2],     0,  d[0]],
                     [-d[1], -d[0],    0]], dtype=np.float64)

    mtx = ddt + np.cos(angle) * (eye - ddt) + np.sin(angle) * skew
    return mtx

def mtx_to_homogenized( mtx ):
    zz = np.zeros((3,), dtype=np.float64)
    return np.r_[np.c_[mtx, zz], [np.r_[zz, 1]]]

class Object(object):
    def __contains__(self, point):
        return False


class Ellipsoid(object):
    def __init__(self, semiaxes, centre, rot_axis, rot_angle):
        """
        Parameters:

            semiaxes : (float, float, float)
                The semiaxes a, b, c of the ellipsoid.

            centre : (float, float, float)
                The position of the ellipsoid's centre in space.

            rot_axis : (float, float, float)
                The direction vector of rotation axis defining the orientation
                in space.

            rot_angle : float
                The rotation angle around the rotation axis.
            
        """
        self.semiaxes = np.array(semiaxes, dtype=np.float64)
        self.centre = np.array(centre, dtype=np.float64)
        self.rot_axis = np.array(rot_axis, dtype=np.float64)
        self.rot_angle = rot_angle

        self.volume = 4.0 / 3.0 * np.pi * np.prod(self.semiaxes)
        self.mtx0 = np.diag( 1.0 / (self.semiaxes**2) )

        self.rot_mtx = make_axis_rotation_matrix(self.rot_axis, self.rot_angle)
        self.mtx = np.dot(self.rot_mtx.T, np.dot(self.mtx0, self.rot_mtx))
        
    def __contains__( self, point ):
        """
        Point x in ellipsoid A <=> x^T A x <= 1.
        """
        x = point - self.centre
        aux = np.dot(x, np.dot(self.mtx, x))
        return aux <= 1.0

    def contains( self, points ):
        """
        Point x in ellipsoid A <=> x^T A x <= 1.
        Works for array of points.

        Parameters:
            points : (n_point, 3) array
        """
        points = np.array(points, ndmin=2, dtype=np.float64)
        x = points.T - self.centre[:,np.newaxis]
        aux = np.sum(x * np.dot(self.mtx, x), axis = 0)
        return np.where(aux <= 1.0, True, False)

    def intersects( self, other ):
        A = mtx_to_homogenized( self.mtx )
        B = mtx_to_homogenized( other.mtx )
        roots = eigh(np.dot(-inv(A), B), eigvals_only=True)
        print roots
        if roots[2] > 0:
            if roots[2] != root[3]:
                return 2
            else:
                return 1
        else:
            return 0

usage = """%prog [options] filename_in"""

help = {
    'filename' :
    'basename of output file(s) [default: %default]',
    'output_format' :
    'output file format, one of: {png, jpg} [default: %default]',
    'n_slice' :
    'number of slices to generate [default: %default]',
    'dims' :
    'dimensions of specimen in milimetres [default: %default]',
    'resolution' :
    'figure resolution [default: %default]',
    'fraction' :
    'volume fraction of objects [default: %default]',
    'length_to_width' :
    'length-to-width ratio of objects [default: %default]',
    'n_object' :
    'number of objects [default: %default]',
}

def main():
    parser = OptionParser(usage=usage, version="%prog ")
    parser.add_option("-o", "", metavar='filename',
                      action="store", dest="output_filename_trunk",
                      default='./slices/slice', help=help['filename'])
    parser.add_option("-f", "--format", metavar='format',
                      action="store", dest="output_format",
                      default="png", help=help['output_format'])
    parser.add_option("-n", "--n-slice", type=int, metavar='int',
                      action="store", dest="n_slice",
                      default=100, help=help['n_slice'])
    parser.add_option("-d", "--dims", metavar='dims',
                      action="store", dest="dims",
                      default='(10, 10, 10)', help=help['n_slice'])
    parser.add_option("-r", "--resolution", metavar='resolution',
                      action="store", dest="resolution",
                      default='600x600', help=help['resolution'])
    parser.add_option("", "--fraction", type=float, metavar='float',
                      action="store", dest="fraction",
                      default='0.3', help=help['fraction'])
    parser.add_option("", "--length-to-width", type=float, metavar='float',
                      action="store", dest="length_to_width",
                      default='8.0', help=help['length_to_width'])
    parser.add_option("", "--n-object", type=int, metavar='int',
                      action="store", dest="n_object",
                      default='5', help=help['n_object'])
    options, args = parser.parse_args()

    options.dims = eval(options.dims)
    options.resolution = [int( r ) for r in  options.resolution.split('x')]
    print options

    total_volume = np.prod(options.dims)
    total_object_volume = options.fraction * total_volume
    average_object_volume = total_object_volume / options.n_object

    print 'total volume [mm^3]: %.2f' % total_volume
    print 'total object volume [mm^3]: %.2f' % total_object_volume
    print 'average object volume [mm^3]: %.2f' % average_object_volume

    object_volume = 0.0
    els = []
    for ii in xrange( options.n_object ):
        print ii
        volume = average_object_volume
        b = c = np.power(volume / (4.0/3.0 * np.pi * options.length_to_width),
                         1.0/3.0)
        a = options.length_to_width * b
        el = Ellipsoid((a, b, c), (0, 0, 0), (1, 1, 1), np.pi/4)
        for prev in els:
            print prev.intersects(el)
        print 'el:', el.semiaxes, el.volume
        els.append(el)
    pause()
##     el = Ellipsoid((1, 1, 1), (0, 0, 0), (0, 0, 1), 0)
##     print el.volume
##     print (0, 0, 0) in el
##     print (1, 0, 0) in el
##     print (1.1, 0, 0) in el

    el = Ellipsoid((1.5, 2.3, 1), (-1, 1, 0), (1, 1, 1), np.pi/4)
    print el.volume
    print (0, 0, 0) in el
    print (1, 0, 0) in el
    print (1.1, 0, 0) in el
    print el.contains([(0, 0, 0), (1, 0, 0), (1.1, 0, 0), (2, 0, 0)])

    x0 = (0, 0, 0)
    v1 = (1, 0, 0)
    v2 = (0, 1, 0)

    xr = np.linspace(-5, 5, 1001)
    yr = xr
    x, y = np.meshgrid(xr, yr)
    z = np.zeros_like( x.flat ) + 0.03
#    import pdb; pdb.set_trace()
    points = np.c_[x.ravel(), y.ravel(), z]
    print points
    mask = el.contains(points)
    print mask
    pl.figure(1)
    pl.clf()
    pl.imshow(mask.reshape(1001, 1001))
    pl.show()
    
if __name__ == "__main__":
    main()
