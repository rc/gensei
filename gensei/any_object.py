import gensei.geometry as gm
from gensei.base import np, Object, pause, ordered_iteritems, is_sequence
from gensei.utils import evaluate, format_dict

class AnyObject(Object):
    """
    Base class of all geometrical objects to be sliced.
    """

    @staticmethod
    def from_conf(conf, box):
        """
        Return an object class instance, based on `conf.kind`.
        The `conf.kind` should correspond to a module in the 'gensei'
        package directory. That module should contain a class named as the
        capitalized `conf.kind`, for example::
            conf.kind = 'ellipsoid'
        means the module `gensei.ellipsoid` with the `Ellipsoid` class.
        """
        try:
            mod = __import__('gensei.' + conf.kind, fromlist=conf.kind)
            cls = getattr(mod, conf.kind.capitalize())
        except (ImportError, AttributeError):
            raise ValueError('unknown object kind! (%s)' % conf.kind)
        
        obj = cls.from_conf(conf, box)
        return obj

    def set_conf(self, conf, requested_conf):
        self.conf = conf
        self.conf.n_object = 0
        self.requested_conf = requested_conf

    def setup_orientation(self):
        """
        If `direction` (orientation of the long axis of an object in space) is
        set in `self.conf`, compute the corresponding `rot_axis` (the direction
        vector of rotation axis) and `rot_angle` (the rotation angle around the
        rotation axis), that map the unrotated object (with `direction0`
        orientation) to the rotated one.

        From `rot_axis` and `rot_angle` form the rotation matrix
        `rot_mtx`, so that:

          - `direction = dot(rot_mtx.T, direction0)`
          - `direction0 = dot(rot_mtx, direction)`

        If `direction` is not set in `self.conf`, use `rot_axis` and
        `rot_angle` from `self.conf`, and compute the `direction` using the
        above relation.
        """
        if hasattr(self.conf, 'direction'):
            nn = np.linalg.norm
            self.direction = evaluate(self.conf.direction, shape=(3,))

            rd = np.dot(self.direction0, self.direction)
            
            self.rot_angle = np.arccos(rd / (nn(self.direction)
                                             * nn(self.direction0)))

            direction = self.direction
            while 1:
                self.rot_axis = np.cross(self.direction0, direction)
                ra = np.linalg.norm(self.rot_axis)
                if ra > 1e-5:
                    break
                direction = evaluate('random', shape=(3,))

            self.rot_axis /= ra

        else:
            self.rot_axis = evaluate(self.conf.rot_axis, shape=(3,))
            self.rot_angle = evaluate(self.conf.rot_angle)

        self.rot_mtx = gm.make_axis_rotation_matrix(self.rot_axis,
                                                    self.rot_angle)

        if not hasattr(self.conf, 'direction'):
            self.direction = np.dot(self.rot_mtx.T, self.direction0)

    def accepted(self):
        """
        Called when an object's placement is accepted.

        Just adds one to the actual object class count.
        """
        self.conf.n_object += 1

    def update_stats(self, stats):
        """
        Add the object volume, surface and length to global statistics.
        """
        stats.volume += self.volume
        stats.surface += self.surface
        stats.length += self.length

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

    def get_intersector(self):
        return self.intersector

    def get_origin_bounding_box(self):
        """
        Get the objects's axes-aligned bounding box as if centered at the
        origin.
        
        Return:
            bbox : 3 x 2 array
                The bounding box.
        """
        return self.intersector.get_origin_bounding_box()

    def get_aligned_bounding_box(self):
        """
        Get the objects's axes-aligned bounding box.
        
        Return:
            bbox : 3 x 2 array
                The bounding box.
        """
        return self.intersector.get_aligned_bounding_box()

    def intersects(self, other):
        """
        Test whether two objects intersect. Some objects may provide only an
        approximate answer using a kind of bounding box (an intersector). In
        that case -1 is returned in place of 1 or 2, and the objects may not
        actually intersect (but should be mutually very close).
        
        Returns
        -------
        flag : int
            -  0 -> the objects are disjoint
            -  1 -> touch in a single surface point
            -  2 -> have common inner points
            - -1 -> may have common surface or inner points (inexact
              computation)
            Some objects may return 2 instead of 1.
        """
        o1 = self.get_intersector()
        o2 = other.get_intersector()

        flag = o1.intersects_fast(o2)

        if flag:
            val = 0
            for segment1 in o1.iter_segments():
                for segment2 in o2.iter_segments():
                    val = segment1.intersects(segment2)
                    if val:
                        break
                if val:
                    break
            flag = val

        ## print flag

        return flag
