#!/usr/bin/env python
import os, glob, copy, time
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import colorConverter

from gensei.base import *
from gensei import Objects, Ellipsoid, Box
from gensei.utils import get_random, get_suffix
from gensei.geometry import get_average_semiaxes

axis_map = {'x' : [1, 2, 0], 'y' : [2, 0, 1], 'z' : [0, 1, 2]}

def get_points(box):
    # All points in the block.
    n_slice = box.n_slice
    if not isinstance(n_slice, dict):
        n_slice = {'z' : int(n_slice)}

    for axis, num in n_slice.iteritems():
        am = axis_map[axis]

        shape = np.array((box.resolution[0], box.resolution[1], num))
        pb = np.zeros((3,), dtype=np.object)
        for ii in range(3):
            pb[am[ii]] = np.linspace(0, box.dims[ii], shape[ii])

        if num > 1:
            delta = pb[am[2]][1] - pb[am[2]][0]
        else:
            delta = 0.0
        x1, x2 = np.meshgrid(pb[am[0]], pb[am[1]])
        x1 = x1.ravel()
        x2 = x2.ravel()

        points = np.empty((x1.shape[0], 3), dtype=np.float64)
        points[:,am[0]] = x1
        points[:,am[1]] = x2

        yield pb, points, delta, num, axis, am

def generate_slices(objects, box, options, output_filename_trunk):
    # Save images of the specimen slices along the z axis of the block. Each
    # image displays a planar cut plane of the block intersecting the
    # ellipsoids.
    resolution = box.resolution

    imshape = resolution[::-1] + (3,)
    aspect = float(resolution[1]) / resolution[0]
    figsize = plt.figaspect(aspect)
    dpi = resolution[0] / figsize[0]

    for pb, points, delta, n_slice, axis, am in get_points(box):
        suffix = get_suffix(n_slice)

        # dpi=dpi in plt.figure() messes with figsize... ???
        fig = plt.figure(1, figsize=figsize, dpi=dpi)
        fig.set_figwidth(figsize[0])
        fig.set_figheight(figsize[1])
        ax = fig.add_axes([0, 0, 1, 1])

        x1b, x2b = pb[am[0]], pb[am[1]]
        for islice, x3b in enumerate(pb[am[2]]):
            x3b_name = ('%05.2f' % x3b).replace('.', '_')
            filename = '.'.join((output_filename_trunk, axis,
                                 suffix % islice, x3b_name,
                                 options.output_format))
            output(islice, x3b, filename, '...')
            output('computing')
            points[:,am[2]] = x3b

            mask = np.zeros(points.shape[0], dtype=np.int8)
            cmask = np.zeros((points.shape[0], 3), dtype=np.float64)
            for obj in objects.itervalues():
                color = np.array(colorConverter.to_rgb(obj.conf.color))

                bbox = obj.get_bounding_box()[am]
                
                ix1 = np.where((x1b > bbox[0,0]) & (x1b < bbox[0,1]))[0]
                ix2 = np.where((x2b > bbox[1,0]) & (x2b < bbox[1,1]))[0]
                a, b = np.meshgrid(ix1, resolution[0]*ix2)
                ii = (a + b).ravel()

                _mask = obj.contains(points[ii])
                mask[ii] += _mask
                cmask[ii[_mask]] = color

            assert_(np.alltrue(mask <= 1))
            output('drawing')
            ax.cla()
            ax.set_axis_off()
            ax.imshow(cmask.reshape(imshape), origin='upper')

            output('saving')
            plt.savefig(filename, format=options.output_format, dpi=dpi)
            output('...done')
##        plt.show()

usage = """%prog [options] [filename]

If an input file is given, the object class options have no effect.

Default option values do _not_ override the input file options.
"""

defaults = {
    'fraction' : 0.1,
    'fraction_reduction' : 0.9,
    'length_to_width' : 8.0,

    'n_slice' : 21,
    'dims' : '(10, 10, 10)',
    'units' : 'mm',
    'resolution' : '600x600',
    'n_object' : 10,
    'output_format' : 'png',
    'timeout' : 5.0,
}

default_objects = {
    'class 1' : {
        'kind' : 'ellipsoid',
        'color' : 'r',
        'fraction' : defaults['fraction'],
        'length_to_width' : defaults['length_to_width'],
        'reduce_to_fit' : {'fraction' : defaults['fraction_reduction']},
        'centre' : 'random',
        'rot_axis' : 'random',
        'rot_angle': 'random',
    },
}

default_box = {
    'dims' : defaults['dims'],
    'units' : defaults['units'],
    'resolution' : defaults['resolution'],
    'n_object' : defaults['n_object'],
    'n_slice' : defaults['n_slice'],
}

default_options = {
    'output_format' : defaults['output_format'],
    'timeout' : defaults['timeout'],
}

help = {
    'filename' :
    'basename of output file(s) [default: %default]',
    'output_format' :
    'output file format (supported by the matplotlib backend used) '\
    '[default: %s]' % defaults['output_format'],
    'n_slice' :
    'number of slices to generate [default: %s]' % defaults['n_slice'],
    'dims' :
    'dimensions of specimen in units given by --units [default: %s]' \
     % defaults['dims'],
    'units' :
    'length units to use [default: %s]' % defaults['units'],
    'resolution' :
    'figure resolution [default: %s]' % defaults['resolution'],
    'fraction' :
    'volume fraction of objects [default: %s]' % defaults['fraction'],
    'fraction_reduction' :
    'volume fraction reduction factor [default: %s]' \
    % defaults['fraction_reduction'],
    'length_to_width' :
    'length-to-width ratio of objects [default: %s]' \
    % defaults['length_to_width'],
    'n_object' :
    'number of objects [default: %s]' % defaults['n_object'],
    'timeout' :
    'timeout in seconds for attempts to place more ellipsiods into '\
    'the block [default: %s]' % defaults['timeout'],
}

def main():
    time_start = time.time()

    parser = OptionParser(usage=usage, version="%prog ")
    parser.add_option("-o", "", metavar='filename',
                      action="store", dest="output_filename_trunk",
                      default='./slices/slice', help=help['filename'])
    parser.add_option("-f", "--format", metavar='format',
                      action="store", dest="output_format",
                      default=None, help=help['output_format'])
    parser.add_option("-n", "--n-slice", type=int, metavar='int',
                      action="store", dest="n_slice",
                      default=None, help=help['n_slice'])
    parser.add_option("-d", "--dims", metavar='dims',
                      action="store", dest="dims",
                      default=None, help=help['dims'])
    parser.add_option("-u", "--units", metavar='units',
                      action="store", dest="units",
                      default=None, help=help['units'])
    parser.add_option("-r", "--resolution", metavar='resolution',
                      action="store", dest="resolution",
                      default=None, help=help['resolution'])
    parser.add_option("", "--fraction", type=float, metavar='float',
                      action="store", dest="fraction",
                      default=None, help=help['fraction'])
    parser.add_option("", "--fraction-reduction", type=float, metavar='float',
                      action="store", dest="fraction_reduction",
                      default=None, help=help['fraction_reduction'])
    parser.add_option("", "--length-to-width", type=float, metavar='float',
                      action="store", dest="length_to_width",
                      default=None, help=help['length_to_width'])
    parser.add_option("", "--n-object", type=int, metavar='int',
                      action="store", dest="n_object",
                      default=None, help=help['n_object'])
    parser.add_option("-t", "--timeout", type=float, metavar='float',
                      action="store", dest="timeout",
                      default=None, help=help['timeout'])
    cmdl_options, args = parser.parse_args()

    can_override = set()
    for key, default in defaults.iteritems():
        val = getattr(cmdl_options, key)
        if val is None:
            setattr(cmdl_options, key, default)
        else:
            can_override.add(key)

    if len(args) == 1:
        filename = args[0]
        config = Config.from_file(filename, required=['objects', 'box'],
                                  optional=['options'])
    else:
        conf = {'objects' : default_objects,
                'box' : default_box,
                'options' : default_options}
        config = Config.from_conf(conf, required=['objects', 'box'],
                                  optional=['options'])

    config.override(cmdl_options, can_override)

    if isinstance(config.box['dims'], str):
        config.box['dims'] = eval(config.box['dims'])
    if isinstance(config.box['resolution'], str):
        aux = tuple([int(r) for r in  config.box['resolution'].split('x')])
        config.box['resolution'] = aux

##     print config

    box = Box(**config.box)
    options = Object(name='options', **config.options)

    output(box)
    output(options)
    
    object_classes = Objects.from_conf(config.objects, box)
    print object_classes

    output('total volume [(%s)^3]: %.2f' % (box.units, box.volume))
##     output('total object volume [(%s)^3]: %.2f' % (options.units,
##                                                   total_object_volume))
##     output('average object volume [(%s)^3]: %.2f' % (options.units,
##                                                     average_object_volume))
    spause(""">>> press a key to generate objects
if it takes too long, press <Ctrl-C> and retry with different parameters""")

    objects = object_classes.place_objects(box, options)
    print objects

    output_dir = os.path.dirname(cmdl_options.output_filename_trunk)
    spause(""">>> press a key to save slices in '%s'
all files in that directory will be deleted""" % output_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        files = glob.glob(os.path.join(output_dir, '*'))
        for name in files:
            os.remove(name)

    generate_slices(objects, box, options, cmdl_options.output_filename_trunk)

    time_end = time.time()

    # Save the statistics to a text file.
    reportname = cmdl_options.output_filename_trunk + '_info.txt'
    output('saving report to %s' % reportname)
    fd = open(reportname, 'w')

    fd.write('started: %s\n' % time.ctime(time_start))
    fd.write('elapsed: %.1f [s]\n' % (time_end - time_start))

    box.report(fd)
    options.report(fd)

    fd.write(objects.format_statistics()+'\n')

    object_classes.report(fd)
    objects.report(fd)

    fd.close()

    output('all done.')

if __name__ == "__main__":
    main()
