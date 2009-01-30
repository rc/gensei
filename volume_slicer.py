#!/usr/bin/env python
import os, glob, copy, time
from optparse import OptionParser
import numpy as np
from scipy.linalg import eig, inv
import pylab as pl

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

def get_random(ranges):
    """Get an array of random numbers 0 <= a_i < ranges[i]."""
    ranges = np.atleast_1d(ranges)
    return ranges * np.random.random(len(ranges))

def get_average_semiaxes(volume, length_to_width):
    """Get semiaxes of an ellipsoid given its volume and length-to-width
    ratio."""
    b = c = np.power(volume / (4.0/3.0 * np.pi * length_to_width),
                     1.0/3.0)
    a = length_to_width * b
    return a, b, c

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
    for key, val in d.iteritems():
        if val == raw[key]:
            msg += (' ' * indent) + ('%s : %s\n' % (key, val))
        else:
            msg += (' ' * indent) + ('%s : %s (%s)\n' % (key, val, raw[key]))
    return msg

def transform_to_pixels(coors, max_coors, resolution):
    """
    Transform real coordinates (in [0, max_coors]) to pixel coordinates, given
    the figure resolution and max. coordinates (block dimensions).
    """
    pc = np.asarray(resolution)[np.newaxis,:] * coors / max_coors[np.newaxis,:]
    return pc.astype(np.int32)

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
        self.rot_axis = np.array(rot_axis, dtype=np.float64)
        self.rot_angle = rot_angle

        self.volume = 4.0 / 3.0 * np.pi * np.prod(self.semiaxes)
        self.mtx0 = np.diag(1.0 / (self.semiaxes**2))

        self.rot_mtx = make_axis_rotation_matrix(self.rot_axis, self.rot_angle)
        self.mtx = np.dot(self.rot_mtx.T, np.dot(self.mtx0, self.rot_mtx))

        self.rot_mtx_hc = make_rotation_matrix_hc(self.rot_mtx)

        self.set_centre(centre)

    def set_centre(self, centre):
        """Set the ellipsoid's centre and update its description matrix in
        homogenous coordinates."""
        self.centre = np.array(centre, dtype=np.float64)
        self.mtx_hc = self._get_matrix_hc()

    def __str__(self):
        msg = '%s\n' % object.__str__(self)
        msg += '    semiaxes: %s\n' % self.semiaxes
        msg += '    centre: %s\n' % self.centre
        msg += '    rot. axis: %s\n' % self.rot_axis
        msg += '    rot. angle: %s' % (self.rot_angle * 180.0 / np.pi)
        return msg

    def _get_matrix_hc(self):
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

    def get_origin_bounding_box(self):
        """Return:
               bbox : 3 x 2 array
                   The bounding box of the ellipsoid placed in the origin.
        """
        aux = np.sqrt(np.diag(inv(self.mtx)))[:,np.newaxis]
        return np.c_[-aux, aux]

    def __contains__(self, point):
        """
        Point x in ellipsoid A <=> x^T A x <= 1.
        """
        x = point - self.centre
        aux = np.dot(x, np.dot(self.mtx, x))
        return aux <= 1.0

    def contains(self, points):
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

    def intersects(self, other):
        """Test if two ellipsoids self and other intersect.

        Return:
            value : int
                0 -> the ellipsoids are disjoint
                1 -> touch in a single surface point
                2 -> have common inner points
        """
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
    'output file format (supported by the matplotlib backend used) '\
    '[default: %default]',
    'n_slice' :
    'number of slices to generate [default: %default]',
    'dims' :
    'dimensions of specimen in units given by --units [default: %default]',
    'units' :
    'length units to use [default: %default]',
    'resolution' :
    'figure resolution [default: %default]',
    'fraction' :
    'volume fraction of objects [default: %default]',
    'fraction_reduction' :
    'volume fraction reduction factor [default: %default]',
    'length_to_width' :
    'length-to-width ratio of objects [default: %default]',
    'n_object' :
    'number of objects [default: %default]',
    'timeout' :
    'timeout in seconds for attempts to place more ellipsiods into '\
    'the block [default: %default]',
}

def main():
    time_start = time.time()

    parser = OptionParser(usage=usage, version="%prog ")
    parser.add_option("-o", "", metavar='filename',
                      action="store", dest="output_filename_trunk",
                      default='./slices/slice', help=help['filename'])
    parser.add_option("-f", "--format", metavar='format',
                      action="store", dest="output_format",
                      default="png", help=help['output_format'])
    parser.add_option("-n", "--n-slice", type=int, metavar='int',
                      action="store", dest="n_slice",
                      default=21, help=help['n_slice'])
    parser.add_option("-d", "--dims", metavar='dims',
                      action="store", dest="dims",
                      default='(10, 10, 10)', help=help['dims'])
    parser.add_option("-u", "--units", metavar='units',
                      action="store", dest="units",
                      default='mm', help=help['units'])
    parser.add_option("-r", "--resolution", metavar='resolution',
                      action="store", dest="resolution",
                      default='600x600', help=help['resolution'])
    parser.add_option("", "--fraction", type=float, metavar='float',
                      action="store", dest="fraction",
                      default='0.1', help=help['fraction'])
    parser.add_option("", "--fraction-reduction", type=float, metavar='float',
                      action="store", dest="fraction_reduction",
                      default='0.9', help=help['fraction_reduction'])
    parser.add_option("", "--length-to-width", type=float, metavar='float',
                      action="store", dest="length_to_width",
                      default='8.0', help=help['length_to_width'])
    parser.add_option("", "--n-object", type=int, metavar='int',
                      action="store", dest="n_object",
                      default='10', help=help['n_object'])
    parser.add_option("-t", "--timeout", type=float, metavar='float',
                      action="store", dest="timeout",
                      default='5.0', help=help['timeout'])
    options, args = parser.parse_args()

    options.dims = eval(options.dims)
    options.resolution = [int(r) for r in  options.resolution.split('x')]
    print options

    orig_options = copy.deepcopy(options.__dict__)

    # Adjust the volume fraction of ellipsoids so that they fit in the
    # specimen's block.
    total_volume = np.prod(options.dims)
    while 1:
        total_object_volume = options.fraction * total_volume
        average_object_volume = total_object_volume / options.n_object

        semiaxes = get_average_semiaxes(average_object_volume,
                                        options.length_to_width)
        print 'fraction: %s, major semiaxis: %s' % (options.fraction,
                                                    semiaxes[0])
        box_dims = np.array(options.dims, dtype=np.float64)
        rbox_conservative = box_dims - 2.0 * semiaxes[0]
        if np.alltrue(rbox_conservative > 0):
            break
        options.fraction *= options.fraction_reduction

    print 'total volume [(%s)^3]: %.2f' % (options.units, total_volume)
    print 'total object volume [(%s)^3]: %.2f' % (options.units,
                                                  total_object_volume)
    print 'average object volume [(%s)^3]: %.2f' % (options.units,
                                                    average_object_volume)
    raw_input(""">>> press <Enter> to generate objects
if it takes too long, press <Ctrl-C> and retry with different parameters""")

    # Generate non-intersecting ellipsoids fully contained in the specimen's
    # block.
    object_volume = 0.0
    els = []
    for ii in xrange(options.n_object):
        print ('\n*** %d ' % ii) + 70*'*' + '\n'

        t0 = time.clock()
        ok = True
        while 1:
            if (time.clock() - t0) > options.timeout:
                print 'timeout!'
                print '-> try reducing --fraction'
                print '   or adjusting --length-to-width, --n-object,'\
                      ' --timeout options'
                ok = False
                break
            # Ensure the whole ellipsoid in the box. 
            axis = get_random((1.0, 1.0, 1.0))
            angle = get_random(np.pi)
            el = Ellipsoid(semiaxes, (0.0, 0.0, 0.0), axis, angle)
            bbox = el.get_origin_bounding_box()
##             print bbox
            rbox = box_dims - 2 * bbox[:,1]
            centre = get_random(rbox) + bbox[:,1]
            el.set_centre(centre)
##             print 'centre:', centre
##             print 'rot. axis, angle:', axis, angle
#            el = Ellipsoid(semiaxes, centre, (0,0,1), 135*np.pi/180)
            for ip, prev in enumerate(els):
                bad = prev.intersects(el)
##                 print '%d. intersects: %d' % (ip, bad)
                if bad:
                    break
            else:
##                 print 'ok'
                break

        if ok:
            print 'accepted:', el
            els.append(el)
            object_volume += el.volume
        else:
            break

    total_volume_error = abs(total_object_volume - object_volume)
    print 'total volume error:', total_volume_error

    output_dir = os.path.dirname(options.output_filename_trunk)
    raw_input(""">>> press <Enter> to save slices in '%s'
all files in that directory will be deleted""" % output_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        files = glob.glob(os.path.join(output_dir, '*'))
        for name in files:
            os.remove(name)
    suffix = get_suffix(options.n_slice)

    # All points in the block.
    xb = np.linspace(0, box_dims[0], options.resolution[0])
    yb = np.linspace(0, box_dims[1], options.resolution[1])
    zb = np.linspace(0, box_dims[2], options.n_slice)
    if options.n_slice > 1:
        dz = zb[1] - zb[0]
    else:
        dz = 0.0
    x, y = np.meshgrid(xb, yb)
    x = x.ravel()
    y = y.ravel()
    z = np.empty_like(x.flat)

    # Save images of the specimen slices along the z axis of the block. Each
    # image displays a planar cut plane of the block intersecting the
    # ellipsoids.
    imshape = options.resolution[::-1]
    aspect = float(options.resolution[1]) / options.resolution[0]
    figsize = pl.figaspect(aspect)
    dpi = options.resolution[0] / figsize[0]

    pl.figure(1, figsize=figsize, dpi=dpi)
    ax = pl.gcf().add_axes([0, 0, 1, 1])
    for iz, zb1 in enumerate(zb):
        zb_name = ('%05.2f' % zb1).replace('.', '_')
        filename = '.'.join((options.output_filename_trunk,
                             suffix % iz, zb_name, options.output_format))
        print iz, zb1, filename
        print 'computing...'
        z.fill(zb1)
        points = np.c_[x, y, z]

        mask = np.zeros(points.shape[0], dtype=np.int32)
        for el in els:
            mask += el.contains(points)

        print 'drawing...'
        ax.cla()
        ax.set_axis_off()
        ax.imshow(mask.reshape(imshape), origin='upper')

        print 'saving...'
        pl.savefig(filename, format=options.output_format, dpi=dpi)
        print 'done.'
##        pl.show()

    time_end = time.time()

    # Save the statistics to a text file.
    reportname = options.output_filename_trunk + '_info.txt'
    print 'saving report to %s...' % reportname
    fd = open(reportname, 'w')
    fd.write('started: %s\n' % time.ctime(time_start))
    fd.write('elapsed: %.1f [s]\n' % (time_end - time_start))
    fd.write('-'*50 + '\n')
    fd.write('dimensions of specimen [%s]: (%f, %f, %f)\n' %\
             ((options.units,) + options.dims))
    fd.write('volume of specimen [(%s)^3]: %f\n' % (options.units, total_volume))
    fd.write('number of slices: %d\n' % options.n_slice)
    fd.write('slice distance [%s]: %f\n' % (options.units, dz))
    fd.write('%d (required: %d) objects (ellipsiods):\n' % (len(els),
                                                            options.n_object))
    fd.write('  volume fraction: %f\n' % options.fraction)
    fd.write('  total volume [(%s)^3]: %f\n' % (options.units,
                                                total_object_volume))
    fd.write('  total volume error [(%s)^3]: %f\n' % (options.units,
                                                      total_volume_error))
    fd.write('  average volume [(%s)^3]: %f\n' % (options.units,
                                                  average_object_volume))
    fd.write('  length-to-width ratio: %f\n' % options.length_to_width)
    fd.write('  semiaxes [%s]: (%f, %f, %f)\n' % ((options.units,) + semiaxes))
    fd.write('-'*50 + '\n')
    fd.write('run with adjusted (raw) options:\n')
    fd.write(format_dict(options.__dict__, raw=orig_options))
    fd.close()
    print 'done.'
    
if __name__ == "__main__":
    main()
