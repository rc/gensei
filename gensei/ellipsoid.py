import gensei.geometry as gm
from gensei.base import np, Object, pause
from gensei.single_object import SingleObject
from gensei.geometry import get_average_semiaxes
from gensei.intersectors import EllipsoidIntersector

class Ellipsoid(SingleObject):
    traits = {
        'semiaxes' : '%s',
        'centre' : '%s',
        'rot_axis' : '%s',
        'rot_angle' : ('%.2f', lambda x: x * 180.0 / np.pi),
    }

    @staticmethod
    def from_conf(conf, box):
        volume = conf.fraction * box.volume / conf.n_object
        semiaxes = get_average_semiaxes(volume, conf.length_to_width)

        obj = Ellipsoid(semiaxes)
        return obj

    def __init__(self, semiaxes):
        """
        Parameters:

            semiaxes : (float, float, float)
                The semiaxes a, b, c of the ellipsoid.
        """
        self.semiaxes = np.array(semiaxes, dtype=np.float64)
        self.volume = 4.0 / 3.0 * np.pi * np.prod(self.semiaxes)
        self.surface = self.compute_approximate_surface()
        self.mtx0 = np.diag(1.0 / (self.semiaxes**2))
        self.is_placed = False
        self.intersection_counters = {}
        self.intersector = EllipsoidIntersector()

    def compute_approximate_surface(self):
        """Approximate Knud Thomsen's formula, where p about 1.6075 yields a
        relative error of at most 1.061%."""
        p = 1.6075
        a, b, c = self.semiaxes
        val = ((np.power(a, p) * np.power(b, p))
                + (np.power(a, p) * np.power(c, p))
                + (np.power(b, p) * np.power(c, p))) / 3.0
        return 4.0 * np.pi * np.power(val, 1.0/p)

    def get_radius(self):
        return self.semiaxes.max()

    def __contains__(self, point):
        """
        Point x in ellipsoid A <=> x^T A x <= 1.
        """
        x = point - self.centre
        aux = np.dot(x, np.dot(self.mtx, x))
        return aux <= 1.0

    def contains(self, points):
        """
        Point x in ellipsoid A <=> x^T A x <= 1.
        Works for array of points.

        Parameters:
            points : (n_point, 3) array
        """
        points = np.array(points, ndmin=2, dtype=np.float64)
        x = points.T - self.centre[:,np.newaxis]
        aux = np.sum(x * np.dot(self.mtx, x), axis = 0)
        mask = np.where(aux <= 1.0, True, False)
##         x2 = np.r_[points.T, np.ones((1,points.shape[0]))]
##         aux2 = np.sum(x2 * np.dot(self.mtx_hc, x2), axis = 0)
##         mask2 = np.where(aux2 <= 0.0, True, False)
##         print np.alltrue(mask == mask2)
        
        return mask

    def get_intersector(self):
        return self.intersector
