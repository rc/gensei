from scipy.linalg import eig, inv

import gensei.geometry as gm
from gensei.base import np, Object, pause
from gensei.geometry import get_average_semiaxes

class Ellipsoid(Object):
    traits = {
        'semiaxes' : '%s',
        'centre' : '%s',
        'rot_axis' : '%s',
        'rot_angle' : ('%.2f', lambda x: x * 180.0 / np.pi),
    }

    @staticmethod
    def from_conf(conf, box):
        volume = conf.fraction * box.volume
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

    def set_conf(self, conf, requested_conf):
        self.conf = conf
        self.requested_conf = requested_conf
            
    def setup_orientation(self):
        """
        Sets rot_axis, the direction vector of rotation axis defining the
        orientation in space, and rot_angle, the rotation angle around the
        rotation axis according to self.conf.
        """

        self.volume = 4.0 / 3.0 * np.pi * np.prod(self.semiaxes)
        self.mtx0 = np.diag(1.0 / (self.semiaxes**2))

        self.rot_axis = np.array(rot_axis, dtype=np.float64)
        self.rot_angle = rot_angle

        self.rot_mtx = gm.make_axis_rotation_matrix(self.rot_axis,
                                                    self.rot_angle)
        self.mtx = np.dot(self.rot_mtx.T, np.dot(self.mtx0, self.rot_mtx))

        self.rot_mtx_hc = gm.make_rotation_matrix_hc(self.rot_mtx)

    def set_centre(self, centre):
        """Set the ellipsoid's centre and update its description matrix in
        homogenous coordinates."""
        self.centre = np.array(centre, dtype=np.float64)
        self.mtx_hc = self._get_matrix_hc()

    def _get_matrix_hc(self):
        """
        Get the matrix describing the ellipsoid in homogenous coordinates.
        It incorporates both the rotation and translation.
        
        Return:
               mtx_hc : 4 x 4 array
                   The matrix describing the ellipsoid in homogenous
                   coordinates.
        """
        M0 = np.zeros((4, 4), dtype=np.float64)
        M0[:3,:3] = self.mtx0
        M0[3, 3] = -1

        M1 = np.dot(self.rot_mtx_hc.T, np.dot(M0, self.rot_mtx_hc))

        T = gm.make_translation_matrix_hc(-self.centre)

        mtx_hc = np.dot(T.T, np.dot(M1, T))

##         print M0
##         print M1
##         print T
##         print mtx_hc
        
        return mtx_hc

    def get_radius(self):
        return self.semiaxes.max()

    def get_origin_bounding_box(self):
        """
        Get the ellipsoid's bounding box as if centered at the origin.
        
        Return:
            bbox : 3 x 2 array
                The bounding box.
        """
        aux = np.sqrt(np.diag(inv(self.mtx)))[:,np.newaxis]
        return np.c_[-aux, aux]

    def get_bounding_box(self):
        """
        Get the ellipsoid's bounding box.
        
        Return:
            bbox : 3 x 2 array
                The bounding box.
        """
        obb = self.get_origin_bounding_box()
        return obb + self.centre[:,np.newaxis]

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

    def intersects(self, other):
        """Test if two ellipsoids self and other intersect.

        Return:
            value : int
                0 -> the ellipsoids are disjoint
                1 -> touch in a single surface point
                2 -> have common inner points
        """
        A, B = self.mtx_hc, other.mtx_hc
        eigs = eig(np.dot(-inv(A), B), left=False, right=False).real
        roots = np.sort(eigs)
#        print roots
        if roots[2] > 0:
            if roots[2] != roots[3]:
                return 0
            else:
                return 1
        else:
            return 2
