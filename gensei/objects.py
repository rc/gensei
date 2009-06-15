import time

from gensei.base import np, output, Object, pause, assert_
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
        
        objs = Objects()

        for key, val in conf.iteritems():
            print key

            try:
                cls = eval(val['kind'].capitalize())
            except NameError:
                raise ValueError('unknown object kind! (%s)' % val['kind'])

            obj_conf = Object.objects_from_dict(val, name='obj_conf', flag=(1,))
            obj_conf0 = obj_conf.copy(deep=True)

            obj = reduce_to_fit(cls, obj_conf, box)
            obj.set_conf(obj_conf, obj_conf0)
            
            print obj
##             print obj.conf
##             print obj.requested_conf
            objs[key] = obj
        
        return objs
        
    def __init__(self, names=None, objs=None):
        if names is None:
            names, objs = [], []
        dict.__init__(self, zip(names, objs))
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
        objs = Objects()

        n_object = self.init_counts(box.n_object)
        print n_object
        
        object_volume = 0.0
        for key in self.names:
            obj_class = self[key]
            
            for ii in xrange(n_object[key]):
                output(('*** %s: %d ' % (key, ii)) + 50*'*' + '\n')

                obj = obj_class.copy(deep=True)

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
                else:
                    break

        objs.init_trait('total_object_volume', object_volume)
        objs.init_trait('total_object_fraction', object_volume / box.volume)

        return objs
