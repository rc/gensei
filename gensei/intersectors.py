from scipy.linalg import eig, inv

from gensei.base import np, Object, pause

class Intersector(Object):
    """
    Base class defining the intersector interface.

    An intersector is a generalization of the concept of bounding box. Its
    purpose is to quickly decide that two objects do not intersect. False
    positives are, on the other hand possible. The possible results of an
    intersector test are thus:

    - objects do not intersect (no common inner or surface point)
    - objects might intersect
    """
    def __init__(self, **kwargs):

        self.ibox = [[0, 0, 0],
                     [1, 0, 0],
                     [0, 1, 0],
                     [1, 1, 0],
                     [0, 0, 1],
                     [1, 0, 1],
                     [0, 1, 1],
                     [1, 1, 1]]

        self.segments = []

        self.set_data(**kwargs)

    def set_centre(self, centre):
        """
        Set the intersector's centre.
        """
        self.centre = np.array(centre, dtype=np.float64)

    def intersects_fast(self, other):
        """
        Test intersection of axes-aligned bounding boxes.
        """
        sbox = self.get_aligned_bounding_box()
        obox = other.get_aligned_bounding_box()

        val = np.empty((8, 3), dtype=np.float64)
        for ii, ib in enumerate(self.ibox):
            val[ii,0] = sbox[0, ib[0]]
            val[ii,1] = sbox[1, ib[1]]
            val[ii,2] = sbox[2, ib[2]]

        flag = np.any((obox[:,0] <= val) & (val <= obox[:,1]))

        return flag

    def iter_segments(self):
        return iter(self.segments)

class EllipsoidIntersector(Intersector):
    """
    The intersector corresponding to a bounding ellipsoid.
    """
    traits = {
        'mtx' : '%s',
        'mtx_hc' : '%s',
        'centre' : '%s',
    }

    def __init__(self, mtx=None, mtx_hc=None, centre=None):
        Intersector.__init__(self, mtx=mtx, mtx_hc=mtx_hc, centre=centre)

        self.segments = [self]

    def set_data(self, mtx=None, mtx_hc=None, centre=None):
        if mtx is not None:
            self.mtx = mtx

        if mtx_hc is not None:
            self.mtx_hc = mtx_hc

        if centre is not None:
            self.centre = centre

    def get_origin_bounding_box(self):
        """
        Get the ellipsoid's axes-aligned bounding box as if centered at the
        origin.
        
        Return:
            bbox : 3 x 2 array
                The bounding box.
        """
        aux = np.sqrt(np.diag(inv(self.mtx)))[:,np.newaxis]
        return np.c_[-aux, aux]

    def get_aligned_bounding_box(self):
        """
        Get the ellipsoid's axes-aligned bounding box.
        
        Return:
            bbox : 3 x 2 array
                The bounding box.
        """
        obb = self.get_origin_bounding_box()
        return obb + self.centre[:,np.newaxis]

    def intersects(self, other):
        """Test if two ellipsoids self and other intersect.

        Returns
        -------
        flag : int
            - 0 -> the ellipsoids are disjoint
            - 1 -> touch in a single surface point
            - 2 -> have common inner points
        """
        A, B = self.mtx_hc, other.mtx_hc
        eigs = eig(np.dot(-inv(A), B), left=False, right=False).real
        roots = np.sort(eigs)

        ## print A, B, roots
        if roots[2] > 0:
            if roots[2] != roots[3]:
                return 0
            else:
                return 1
        else:
            return 2
