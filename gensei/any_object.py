from gensei.base import np, Object, pause

class AnyObject(Object):
    """
    Base class of all geometrical objects to be sliced.
    """

    @staticmethod
    def from_conf(conf, box):
        from gensei.ellipsoid import Ellipsoid

        try:
            cls = eval(conf.kind.capitalize())
        except NameError:
            raise ValueError('unknown object kind! (%s)' % val['kind'])
        
        obj = cls.from_conf(conf, box)
        return obj
