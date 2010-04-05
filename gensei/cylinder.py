import gensei.geometry as gm
from gensei.base import np, Object, pause
from gensei.single_object import SingleObject
from gensei.geometry import get_cylinder_shape
from gensei.intersectors import EllipsoidIntersector

class Cylinder(SingleObject):
    traits = {
        'radius' : '%s',
        'height' : '%s',
        'centre' : '%s',
        'rot_axis' : '%s',
        'rot_angle' : ('%.2f', lambda x: x * 180.0 / np.pi),
    }

    @staticmethod
    def from_conf(conf, box):
        volume = conf.fraction * box.volume / conf.n_object
        radius, height = get_cylinder_shape(volume, conf.length_to_width)

        obj = Cylinder(radius, height)
        return obj

    def __init__(self, radius, height):
        """
        Parameters
        ----------
        radius : float
            The radius of the cylinder.
        height : float
            The height of the cylinder.
        """
        self.radius = radius
        self.height = height
        self.length_to_width = 0.5 * height / radius
        self.volume = height * np.pi * (radius**2)
        self.surface = 2.0 * np.pi * radius * height
        self.mtx0 = self.get_circumscribed_ellipsoid()
        self.is_placed = False
        self.intersection_counters = {}
        self.intersector = EllipsoidIntersector()

    def get_circumscribed_ellipsoid(self):
        """
        Return a circumscribed ellipsoid of the cylinder.
        """
        semiaxes = np.array(3 * [np.sqrt(2.0) * self.radius], dtype=np.float64)
        semiaxes[2] *= self.length_to_width
        mtx = np.diag(1.0 / (semiaxes**2))

        return mtx
        
    def get_radius(self):
        r = np.sqrt((0.5 * self.height)**2 + self.radius**2)
        return r

    def __contains__(self, point):
        """
        Point x in cylinder.
        """
        x = point - self.centre
        aux = np.dot(self.rot_mtx, x)
        r = np.sqrt(aux[0]**2 + aux[1]**2)
        val = (np.abs(aux[2]) <= self.height) & (r <= self.radius)

        return val

    def contains(self, points):
        """
        Point x in cylinder. Works for array of points.

        Parameters
        ----------
        points : (n_point, 3) array
            The points to be tested for inclusion.
        """
        points = np.array(points, ndmin=2, dtype=np.float64)
        x = points.T - self.centre[:,np.newaxis]
        aux = np.dot(self.rot_mtx, x)
        r = np.sqrt(aux[0,:]**2 + aux[1,:]**2)
        mask = (np.abs(aux[2,:]) <= (0.5 * self.height)) & (r <= self.radius)

        return mask
