from gensei.base import np, Object, pause, ordered_iteritems
from gensei.utils import format_dict

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

    def set_conf(self, conf, requested_conf):
        self.conf = conf
        self.conf.n_object = 0
        self.requested_conf = requested_conf

    def accepted(self):
        """
        Called when an object's placement is accepted.

        Just adds one to the actual object class count.
        """
        self.conf.n_object += 1

    def init_intersection_counters(self, axis):
        """Initialize for using store_intersection()."""
        self.intersection_counters[axis] = []

    def store_intersection(self, mask, axis, coor):
        """
        Store intersection if it occurred.
        
        Parameters
        ----------

        mask : bool array
            Slice mask, True where the object inside is.
        axis : 'x', 'y' or 'z'
            Axis perpendicular to the slices.
        coor: float
            Coordinate along the axis, where intersection might occur.
        """
        if mask.any():
            self.intersection_counters[axis].append(coor)

    def has_intersection(self, axis):
        return len(self.intersection_counters[axis]) > 0

    def report(self, filename):
        fd = self.fd_open(filename)

        if self.is_placed:
            Object.report(self, fd, header=False)
            fd.write('intersections per axis:\n')
            for axis, ints in ordered_iteritems(self.intersection_counters):
                fd.write('  %s (%d): %s\n' % (axis, len(ints), ints))
        else:
            fd.write(format_dict(self.conf.get_dict(),
                                 raw=self.requested_conf.get_dict()) + '\n')
            
        self.fd_close()

    def intersects(self, other):
        """
        Test whether two objects intersect. Some objects may provide only an
        approximate answer using a kind of bounding box. In that case -1 is
        returned in place of 1 or 2, and the objects may not actually intersect
        (but should be mutually very close).
        
        Returns
        -------

        flag : int
             0 -> the objects are disjoint
             1 -> touch in a single surface point
             2 -> have common inner points
            -1 -> may have common surface or inner points (inexact computation)
            Some objects may return 2 instead of 1.
        """
        o1 = self.get_intersector()
        o2 = other.get_intersector()

        flag = o1.intesects_fast(o2)
        if flag:
            val = 0
            for segment1 in o1.iter_segments():
                for segment2 in o1.iter_segments():
                    val = segment1.intersects(segment2)
                    if val:
                        break
                if val:
                    break
            flag = val

        return flag
