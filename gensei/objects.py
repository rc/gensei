import time

from gensei.base import np, output, Object, pause, assert_, _dashes, \
     ordered_iteritems
from gensei.utils import get_random
from gensei.ellipsoid import Ellipsoid

def reduce_to_fit(cls, obj_conf, box):
    """
    Adjusts obj_conf parameters of object of class cls so that the circumscribed
    sphere radius fits in the box.
    """
    rtf = obj_conf.reduce_to_fit

    while 1:
        obj = cls.from_conf(obj_conf, box)

        r = obj.get_radius()
        rbox = 0.5 * box.dims.min()
        if rbox > r:
            break

        for attr, ratio in rtf.iteritems():
            val0 = getattr(obj_conf, attr)
            val = ratio * val0
            setattr(obj_conf, attr, val)
            output('reducing %s: %s -> %s' % (attr, val0, val))

    return obj

class Objects(Object, dict):
    traits = {
        'names' : None,
        'objects' : lambda self: dict.__repr__(self),
    }

    @staticmethod
    def from_conf(conf, box):

##         objs = Objects(conf.keys(), conf.values())
##         print objs
        objs = Objects(is_placed=False)

        for key, val in conf.iteritems():
            output(('*** %s ' % key) + 50*'*')

            try:
                cls = eval(val['kind'].capitalize())
            except NameError:
                raise ValueError('unknown object kind! (%s)' % val['kind'])

            obj_conf = Object.objects_from_dict(val, name='obj_conf', flag=(1,))
            obj_conf0 = obj_conf.copy(deep=True)

            obj = reduce_to_fit(cls, obj_conf, box)
            obj.set_conf(obj_conf, obj_conf0)
            
            output(obj._format(mode='set_only'))
##             print obj.conf
##             print obj.requested_conf
            objs[key] = obj
        
        return objs
        
    def __init__(self, names=None, objs=None, is_placed=False):
        if names is None:
            names, objs = [], []

        dict.__init__(self, zip(names, objs))
        Object.__init__(self, is_placed=is_placed)

        self.names = sorted(self.keys())

    def __setitem__(self, key, item):
        super(Objects, self).__setitem__(key, item)

        self.names = sorted(self.keys())
        self.n_class = len(self.names)

    def init_counts(self, count):
        """If count is a number, split it evenly among object classes."""
        if isinstance(count, dict):
            n_object = count
        else:
            n_object = {}
            n_per_class = int(count / self.n_class)
            for name in self.names:
                n_object[name] = n_per_class

            n_extra = count % self.n_class
            for ii in xrange(n_extra):
                ic = ii % self.n_class
                n_object[self.names[ic]] += 1

            assert_(count == sum(n_object.values()))

        return n_object

    def place_objects(self, box, options):
        """Generate non-intersecting objects fully contained in the specimen's
        block.
        """
        objs = Objects(is_placed=True)
        objs.box = box

        n_object = self.init_counts(box.n_object)
        print n_object

        objs.init_trait('n_object', {}.fromkeys(n_object.keys(), 0))
        
        stats_per_class = {}
        for key in self.names:
            obj_class = self[key]

            stats = Object(volume=0.0, surface=0.0)
            stats_per_class[key] = stats
            
            for ii in xrange(n_object[key]):
                output(('*** %s: %d ' % (key, ii)) + 50*'*')

                obj = obj_class.copy(deep=True)
                obj.init_trait('obj_class', key)
                
                t0 = time.clock()
                ok = True

                while 1:
                    if (time.clock() - t0) > options.timeout:
                        output('timeout -> try reducing object size!')
                        ok = False
                        break

                    obj.setup_orientation()
                    bbox = obj.get_origin_bounding_box()

                    rbox = box.dims - 2 * bbox[:,1]
                    centre = get_random(rbox) + bbox[:,1]
                    obj.set_centre(centre)
                    obj.is_placed = True
##                     print obj
                    
                    for ip, prev in enumerate(objs.itervalues()):
                        bad = prev.intersects(obj)
##                         print '%d. intersects: %d' % (ip, bad)
                        if bad:
                            break
                    else:
##                         print 'ok'
                        break

                if ok:
                    output('accepted:', obj)
                    objs[key + ('_%d' % ii)] = obj

                    stats.volume += obj.volume
                    stats.surface += obj.surface

                    objs.n_object[key] += 1
                else:
                    break

        object_volume = object_surface = 0.0
        for stats in stats_per_class.itervalues():
            object_volume += stats.volume
            object_surface += stats.surface

        objs.stats_per_class = stats_per_class
        objs.init_trait('total_object_volume', object_volume)
        objs.init_trait('total_object_volume_fraction',
                        object_volume / box.volume)
        objs.init_trait('total_object_surface', object_surface)
        objs.init_trait('total_object_surface_fraction',
                        object_surface / box.volume)
        objs.section_volumes = {}
        objs.section_surfaces = {}

        return objs

    def init_section_based_data(self, axis):
        """The section based data are stored per axis and object class as
        follows: {axis : {obj_class : value, ...}, ...}
        """
        names = list(set(obj.obj_class for obj in self.itervalues()))
        self.section_volumes[axis] = {}.fromkeys(names, 0.0)
        self.section_surfaces[axis] = {}.fromkeys(names, 0.0)
        for obj in self.itervalues():
            obj.init_intersection_counters(axis)


    def update_section_based_data(self, mask, shape, axis, delta, obj_class):
        """
        Estimate volume and surface of object classes from their sections.
        WARNING: Surface estimation is very crude!

        Parameters
        ----------

        mask : bool array
            Slice mask, True where the object inside is.
        shape : tuple
            Real 2D shape of the mask.
        axis : 'x', 'y' or 'z'
            Axis perpendicular to the slices.
        delta : float
            Slice distance.
        obj_class : str
            Geometrical object class.
        """
        pixel_area = self.box.get_pixel_area(axis)

        volumes = self.section_volumes[axis]
        volumes[obj_class] += pixel_area * np.sum(mask) * delta

        mask.shape = shape
        pixel_sizes = self.box.get_pixel_sizes(axis)
        grad0, grad1 = np.gradient(mask)
        val0 = len(grad0[np.where(grad0)]) * pixel_sizes[0]
        val1 = len(grad1[np.where(grad1)]) * pixel_sizes[1]
##         if np.sum(mask):
##             import pylab as pll
##             pll.spy(grad0+grad1)
##             pll.show()
##             import pdb; pdb.set_trace()
        surfaces = self.section_surfaces[axis]
        surfaces[obj_class] += 0.5 * (val0 + val1) * delta

    def format_intersection_statistics(self, is_output=False):
        if is_output:
            space = ''
        else:
            space = '  '

        msg = []
        for axis, num in ordered_iteritems(self.box.n_slice):
            msg.append(space + 'axis %s (%d sections):' % (axis, num))
            ok = True
            for key, obj in ordered_iteritems(self):
                if not obj.has_intersection(axis):
                    msg.append(space + '  warning: object %s is not intersected'
                               % key)
                    ok = False
            if ok:
                msg.append(space
                           + '  all objects intersected by at least one section')

        if is_output:
            output('\n'.join(msg))

        return msg

    def format_statistics(self):
        units = self.box.units

        msg = [_dashes, 'statistics', _dashes]

        msg.append('  number of objects per class:')
        for key, num in ordered_iteritems(self.n_object):
            msg.append('    %s: %d' % (key, num))
        msg.append(_dashes)

        msg.append('  total object volume fraction:           %f' \
                   % self.total_object_volume_fraction)
        msg.append('  total object volume [(%s)^3]:           %f' \
                   % (units, self.total_object_volume))
        msg.append('  total object surface fraction [1/(%s)]: %f' \
                   % (units, self.total_object_surface_fraction))
        msg.append('  total object surface [(%s)^2]:          %f' \
                   % (units, self.total_object_surface))
        msg.append(_dashes)

        msg.append('  missed objects per axis:')
        missed = {}.fromkeys(self.box.n_slice.keys(), 0)
        for obj in self.itervalues():
            for axis, ints in obj.intersection_counters.iteritems():
                missed[axis] += len(ints) == 0
        for axis, val in ordered_iteritems(missed):
            msg.append('    %s: %d' % (axis, val))
        msg.extend(self.format_intersection_statistics())
        msg.append(_dashes)

        ipac = {}.fromkeys(self.box.n_slice.keys())
        for axis in ipac.iterkeys():
            volumes = self.section_volumes[axis]
            ipac[axis] = {}.fromkeys(volumes.keys(), 0)

        for key, obj in self.iteritems():
            for axis, ints in obj.intersection_counters.iteritems():
                val = len(ints)
                ipac[axis][obj.obj_class] += val

        msg.append('  intersections per class:')
        for axis, iac in ordered_iteritems(ipac):
            msg.append('    axis: %s' % axis)
            for key, ints in ordered_iteritems(iac):
                msg.append('      %s: % 3d average: %5.2f for %d objects' \
                           % (key, ints,
                              float(ints) / self.n_object[key],
                              self.n_object[key]))
        msg.append(_dashes)

        msg.append('  volumes per class:')
        for axis, volumes in ordered_iteritems(self.section_volumes):
            msg.append('    axis: %s' % axis)
            for key, val in ordered_iteritems(volumes):
                stats = self.stats_per_class[key]
                msg.append('      class: %s' % key)
                msg.append('        true volume V1 [(%s)^3]:          %f' \
                           % (units, stats.volume))
                msg.append('        section-based volume V2 [(%s)^3]: %f' \
                           % (units, val))
                msg.append('        true volume fraction:             %f' \
                           % (stats.volume / self.box.volume))
                msg.append('        section-based volume fraction:    %f' \
                           % (val / self.box.volume))
                msg.append('        ratio V2/V1 (accuracy estimate):  %f' \
                           % (val / stats.volume))
        msg.append(_dashes)

        msg.append('  surfaces per class:')
        for axis, surfaces in ordered_iteritems(self.section_surfaces):
            msg.append('    axis: %s' % axis)
            for key, val in ordered_iteritems(surfaces):
                stats = self.stats_per_class[key]
                msg.append('      class: %s' % key)
                msg.append('        true (approximate) surface S1 [(%s)^2]: %f' \
                           % (units, stats.surface))
                msg.append('        section-based surface S2 [(%s)^2]:      %f' \
                           % (units, val))
                msg.append('        true (approximate) surface fraction:    %f' \
                           % (stats.surface / self.box.volume))
                msg.append('        section surface fraction [[1/(%s)]:     %f' \
                           % (units, val / self.box.volume))
                msg.append('        ratio S2/S1 (accuracy estimate):        %f' \
                           % (val / stats.surface))

        return msg

    def report(self, filename):
        fd = self.fd_open(filename)

        if self.is_placed:
            msg = '\n'.join([_dashes, 'objects', _dashes])
            fd.write(msg+'\n')
        else:
            msg = '\n'.join([_dashes, 'object classes', _dashes])
            fd.write(msg+'\n')

        for key in self.names:
            obj_class = self[key]
            msg = '\n'.join([_dashes, key, _dashes])
            fd.write(msg+'\n')

            obj_class.report(fd)

        self.fd_close()
        
