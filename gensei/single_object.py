import gensei.geometry as gm
from gensei.base import np, Object, pause, ordered_iteritems, is_sequence
from gensei.any_object import AnyObject

class SingleObject(AnyObject):
    """
    Base class of all geometrical objects to be sliced that are not compound
    (multi-object).
    """
        
    def setup_orientation(self):
        """
        Sets rot_axis, the direction vector of rotation axis defining the
        orientation in space, and rot_angle, the rotation angle around the
        rotation axis according to self.conf. Also update the intersector.
        """
        AnyObject.setup_orientation(self)

        self.mtx = np.dot(self.rot_mtx.T, np.dot(self.mtx0, self.rot_mtx))

        self.rot_mtx_hc = gm.make_rotation_matrix_hc(self.rot_mtx)

        self.intersector.set_data(mtx=self.mtx)

    def set_centre(self, centre):
        """
        Set the objects's centre and update its description matrix in
        homogenous coordinates. Also update the intersector. 
        """
        self.centre = np.array(centre, dtype=np.float64)
        self.mtx_hc = self._get_matrix_hc()

        self.intersector.set_data(mtx_hc=self.mtx_hc, centre=self.centre)

    def _get_matrix_hc(self):
        """
        Get the matrix describing the circumscribed ellipsoid in homogenous
        coordinates. It incorporates both the rotation and translation.
        
        Returns
        -------
        mtx_hc : 4 x 4 array
            The matrix describing the circumscribed ellipsoid in homogenous
            coordinates.
        """
        M0 = np.zeros((4, 4), dtype=np.float64)
        M0[:3,:3] = self.mtx0
        M0[3, 3] = -1

        M1 = np.dot(self.rot_mtx_hc.T, np.dot(M0, self.rot_mtx_hc))

        T = gm.make_translation_matrix_hc(-self.centre)

        mtx_hc = np.dot(T.T, np.dot(M1, T))
        
        return mtx_hc
