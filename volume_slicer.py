#!/usr/bin/env python
from optparse import OptionParser
import numpy as np
import pylab as pl

class Object( object ):
    def __contains__( self, point ):
        return False


class Ellipsoid( object ):
    def __init__( self, semiaxes, centre, rot_axis, rot_angle ):
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
        # TODO
        self.mtx = self.mtx0
        
    def __contains__( self, point ):
        """
        Point x in ellipsoid A <=> x^T A x <= 0.
        """
        # TODO: vectorize properly.
        x = point - self.centre
        aux = np.dot(x, np.dot(self.mtx, x))
        return np.where(aux <= 0)[0]

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
    parser.add_option("-n", "--n-slice", type=int, metavar='n_slice',
                      action="store", dest="n_slice",
                      default=100, help=help['n_slice'])
    parser.add_option("-d", "--dims", metavar='dims',
                      action="store", dest="dims",
                      default='(10, 10, 10)', help=help['n_slice'])
    parser.add_option("-r", "--resolution", metavar='resolution',
                       action="store", dest="resolution",
                       default='600x600', help=help['resolution'])
    parser.add_option("", "--fraction", type=float, metavar='fraction',
                       action="store", dest="fraction",
                       default='0.3', help=help['fraction'])
    parser.add_option("", "--n-object", type=int, metavar='n_object',
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
    for ii in xrange( options.n_object ):
        print ii


    u = np.linspace(0, 2.0 * np.pi, 100)
    v = u
    a = 10.0
    b = 3.0
    c = 2.0
    
    x = a * np.sin( u ) * np.cos( v )
    y = b * np.cos( u ) * np.cos( v )
    z = c * np.sin( v )

    
    
if __name__ == "__main__":
    main()
