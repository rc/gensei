#!/usr/bin/env python
import os, glob
from optparse import OptionParser
import numpy as np
from scipy.linalg import eig, inv
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
    skew = np.array([[    0,  d[2],  -d[1]],
                     [-d[2],     0,  d[0]],
                     [d[1], -d[0],    0]], dtype=np.float64)

    mtx = ddt + np.cos(angle) * (eye - ddt) + np.sin(angle) * skew
    return mtx

def make_rotation_matrix_hc(mtx):
    zz = np.zeros((3,), dtype=np.float64)
    mtx = np.r_[np.c_[mtx, zz], [np.r_[zz, 1]]]
    return mtx

def make_translation_matrix_hc(point):
    mtx = np.eye(4, dtype=np.float64)
    mtx[:-1,3] = point
    return mtx

def get_random(ranges):
    ranges = np.atleast_1d(ranges)
    return ranges * np.random.random(len(ranges))

def get_average_semiaxes(volume, length_to_width):
    b = c = np.power(volume / (4.0/3.0 * np.pi * length_to_width),
                     1.0/3.0)
    a = length_to_width * b
    return a, b, c

def get_suffix( n ):
    if n > 1:
        n_digit = int( np.log10( n - 1 ) + 1 )
        suffix = '%%0%dd' % n_digit
    else:
        suffix = '%d'
    return suffix

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

        self.rot_mtx_hc = make_rotation_matrix_hc(self.rot_mtx)
        self.mtx_hc = self._get_matrix_hc()

    def _get_matrix_hc( self ):
        """Return:
               mtx_hc : 4 x 4 array
                   The matrix describing the ellipsoid in homogenous
                   coordinates.
        """
        M0 = np.zeros((4, 4), dtype=np.float64)
        M0[:3,:3] = self.mtx0
        M0[3, 3] = -1

        M1 = np.dot(self.rot_mtx_hc.T, np.dot(M0, self.rot_mtx_hc))

        T = make_translation_matrix_hc(-self.centre)

        mtx_hc = np.dot(T.T, np.dot(M1, T))

##         print M0
##         print M1
##         print T
##         print mtx_hc
        
        return mtx_hc

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
        mask = np.where(aux <= 1.0, True, False)
##         x2 = np.r_[points.T, np.ones((1,points.shape[0]))]
##         aux2 = np.sum(x2 * np.dot(self.mtx_hc, x2), axis = 0)
##         mask2 = np.where(aux2 <= 0.0, True, False)
##         print np.alltrue(mask == mask2)
        
        return mask

    def intersects( self, other ):
        A, B = self.mtx_hc, other.mtx_hc
        eigs = eig(np.dot(-inv(A), B), left=False, right=False).real
        roots = np.sort(eigs)
#        print roots
        if roots[2] > 0:
            if roots[2] != roots[3]:
                return 0
            else:
                return 1
        else:
            return 2

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
                      default=101, help=help['n_slice'])
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
    ok = 0
    while 1:
        total_object_volume = options.fraction * total_volume
        average_object_volume = total_object_volume / options.n_object

        semiaxes = get_average_semiaxes(average_object_volume,
                                        options.length_to_width)
        print 'fraction: %s, major semiaxis: %s' % (options.fraction,
                                                    semiaxes[0])
        box_dims = np.array(options.dims, dtype=np.float64)
        rbox = box_dims - 2.0 * semiaxes[0]
        if ok: # postpone one iteration
            break
        if np.alltrue(rbox > 0):
            ok = 1
        options.fraction *= 0.5

    print 'total volume [mm^3]: %.2f' % total_volume
    print 'total object volume [mm^3]: %.2f' % total_object_volume
    print 'average object volume [mm^3]: %.2f' % average_object_volume
    raw_input()

    object_volume = 0.0
    els = []
    for ii in xrange( options.n_object ):
        print ('\n*** %d ' % ii) + 70*'*' + '\n'

        while 1:
            # Ensure the whole ellipsoid in the box. 
            centre = get_random(rbox) + semiaxes[0]
            axis = get_random((1.0, 1.0, 1.0))
            angle = get_random(np.pi/2)
            print 'centre:', centre
            print 'rot. axis, angle:', axis, angle
            el = Ellipsoid(semiaxes, centre, axis, angle)
            for ip, prev in enumerate(els):
                bad = prev.intersects(el)
                print '%d. intersects: %d' % (ip, bad)
                if bad:
                    break
            else:
                print 'ok'
                break

        print 'accepted el:', el.semiaxes, el.volume
        els.append(el)

        object_volume += el.volume
    print 'object volume error:', abs( total_object_volume - object_volume )

    output_dir = os.path.dirname(options.output_filename_trunk)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        files = glob.glob(os.path.join(output_dir, '*'))
        for name in files:
            os.remove(name)
    suffix = get_suffix(options.n_slice)

    # All points in the box.
    xb = np.linspace(0, box_dims[0], options.resolution[0])
    yb = np.linspace(0, box_dims[1], options.resolution[1])
    zb = np.linspace(0, box_dims[2], options.n_slice)
    x, y = np.meshgrid(xb, yb)
    x = x.ravel()
    y = y.ravel()
    z = np.empty_like( x.flat )

    for iz, zb1 in enumerate(zb):
        print iz, zb1
        z.fill(zb1)
        points = np.c_[x, y, z]

        mask = np.zeros(points.shape[0], dtype=np.int32)
        for el in els:
            mask += el.contains(points)

        pl.figure(1)
        pl.clf()
        pl.imshow(mask.reshape(options.resolution), origin='upper')
        pl.axis('off')

        zb_name = ('%05.2f' % zb1).replace('.', '_')
        filename = '.'.join((options.output_filename_trunk,
                             suffix % iz, zb_name, options.output_format))
        pl.savefig(filename, format=options.output_format)
##        pl.show()
    
if __name__ == "__main__":
    main()
