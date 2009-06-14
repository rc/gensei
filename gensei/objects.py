from gensei.base import np, output, Object, pause
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
            obj = reduce_to_fit(cls, obj_conf, box)
            obj.set_conf(obj_conf, val)
            
            print obj
            print obj.conf
            print obj.requested_conf
            pause()
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
        

##         if setup_position:
##             self.setup_position(centre, rot_axis, rot_angle)
