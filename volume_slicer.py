#!/usr/bin/env python
from optparse import OptionParser
import numpy as np

class Object( object ):
    pass

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

if __name__ == "__main__":
    main()
