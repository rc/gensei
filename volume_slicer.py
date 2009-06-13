#!/usr/bin/env python
import os, glob, copy, time
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt

from gensei.base import *
from gensei import Ellipsoid
from gensei.utils import get_random, get_suffix, format_dict
from gensei.geometry import get_average_semiaxes

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
    output(options)

    orig_options = copy.deepcopy(options.__dict__)

    # Adjust the volume fraction of ellipsoids so that they fit in the
    # specimen's block.
    total_volume = np.prod(options.dims)
    while 1:
        total_object_volume = options.fraction * total_volume
        average_object_volume = total_object_volume / options.n_object

        semiaxes = get_average_semiaxes(average_object_volume,
                                        options.length_to_width)
        output('fraction: %s, major semiaxis: %s' % (options.fraction,
                                                     semiaxes[0]))
        box_dims = np.array(options.dims, dtype=np.float64)
        rbox_conservative = box_dims - 2.0 * semiaxes[0]
        if np.alltrue(rbox_conservative > 0):
            break
        options.fraction *= options.fraction_reduction

    output('total volume [(%s)^3]: %.2f' % (options.units, total_volume))
    output('total object volume [(%s)^3]: %.2f' % (options.units,
                                                  total_object_volume))
    output('average object volume [(%s)^3]: %.2f' % (options.units,
                                                    average_object_volume))
    raw_input(""">>> press <Enter> to generate objects
if it takes too long, press <Ctrl-C> and retry with different parameters""")

    # Generate non-intersecting ellipsoids fully contained in the specimen's
    # block.
    object_volume = 0.0
    els = []
    for ii in xrange(options.n_object):
        output(('\n*** %d ' % ii) + 70*'*' + '\n')

        t0 = time.clock()
        ok = True
        while 1:
            if (time.clock() - t0) > options.timeout:
                output('timeout!')
                output('-> try reducing --fraction')
                output('   or adjusting --length-to-width, --n-object,'\
                      ' --timeout options')
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
            output('accepted:', el)
            els.append(el)
            object_volume += el.volume
        else:
            break

    total_volume_error = abs(total_object_volume - object_volume)
    output('total volume error:', total_volume_error)

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
    figsize = plt.figaspect(aspect)
    dpi = options.resolution[0] / figsize[0]

    plt.figure(1, figsize=figsize, dpi=dpi)
    ax = plt.gcf().add_axes([0, 0, 1, 1])
    for iz, zb1 in enumerate(zb):
        zb_name = ('%05.2f' % zb1).replace('.', '_')
        filename = '.'.join((options.output_filename_trunk,
                             suffix % iz, zb_name, options.output_format))
        output(iz, zb1, filename)
        output('computing')
        z.fill(zb1)
        points = np.c_[x, y, z]

##         mask = np.zeros(points.shape[0], dtype=np.int32)
##         for el in els:
## #            tt = time.clock()
##             mask += el.contains(points)
## #            print time.clock() - tt
## #        mask2 = mask

        mask = np.zeros(points.shape[0], dtype=np.int32)
        for el in els:
##             tt = time.clock()
##             bbox = el.get_bounding_box()
##             ii = np.where((x > bbox[0,0]) & (x < bbox[0,1])
##                           & (y > bbox[1,0]) & (y < bbox[1,1]))[0]
##             ii2 = ii
##             print '*', time.clock() - tt
##             tt = time.clock()
            bbox = el.get_bounding_box()
            ix = np.where((xb > bbox[0,0]) & (xb < bbox[0,1]))[0]
            iy = np.where((yb > bbox[1,0]) & (yb < bbox[1,1]))[0]
            a, b = np.meshgrid(ix, options.resolution[0]*iy)
            ii = (a + b).ravel()
##             print time.clock() - tt

##             tt = time.clock()
            mask[ii] += el.contains(points[ii])
##             print time.clock() - tt
#        print np.alltrue(mask==mask2)

        output('drawing')
        ax.cla()
        ax.set_axis_off()
        ax.imshow(mask.reshape(imshape), origin='upper')

        output('saving')
        plt.savefig(filename, format=options.output_format, dpi=dpi)
        output('done')
##        plt.show()

    time_end = time.time()

    # Save the statistics to a text file.
    reportname = options.output_filename_trunk + '_info.txt'
    output('saving report to %s' % reportname)
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
    output('done')
    
if __name__ == "__main__":
    main()
