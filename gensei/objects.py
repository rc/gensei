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
        
        object_volume = 0.0
        stats_per_class = {}
        for key in self.names:
            obj_class = self[key]

            stats = Object(volume=0.0)
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
                    object_volume += obj.volume
                    stats.volume += obj.volume
                else:
                    break

        objs.stats_per_class = stats_per_class
        objs.init_trait('total_object_volume', object_volume)
        objs.init_trait('total_object_fraction', object_volume / box.volume)
        objs.rough_volumes = {}

        return objs

    def init_rough_volumes(self, axis):
        names = list(set(obj.obj_class for obj in self.itervalues()))
        self.rough_volumes[axis] = {}.fromkeys(names, 0.0)

    def update_rough_volumes(self, mask, axis, delta, obj_class):
        """
        Parameters
        ----------

        mask : bool array
            Slice mask, True where the object inside is.
        axis : 'x', 'y' or 'z'
            Axis perpendicular to the slices.
        delta : float
            Slice distance.
        obj_class : str
            Geometrical object class.
        """
        pixel_area = self.box.get_pixel_area(axis)

        rvs = self.rough_volumes[axis]
        rvs[obj_class] += pixel_area * np.sum(mask) * delta

    def format_statistics(self):
        units = self.box.units

        msg = [_dashes, 'statistics', _dashes]

        msg.append('  total object volume fraction: %f' \
                   % self.total_object_fraction)
        msg.append('  total object volume [(%s)^3]: %f' \
                   % (units, self.total_object_volume))
        msg.append(_dashes)
        msg.append('  per class:')
        for axis, rvs in ordered_iteritems(self.rough_volumes):
            msg.append('    axis: %s' % axis)
            for key, val in ordered_iteritems(rvs):
                stats = self.stats_per_class[key]
                msg.append('      class: %s' % key)
                msg.append('        true volume [(%s)^3]:  %f' \
                           % (units, stats.volume))
                msg.append('        rough volume [(%s)^3]: %f' \
                           % (units, val))
                msg.append('        true volume fraction:  %f' \
                           % (stats.volume / self.box.volume))
                msg.append('        rough volume fraction: %f' \
                           % (val / self.box.volume))
            
        return '\n'.join(msg)

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
        
