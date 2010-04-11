from gensei.base import np, Object, pause, ordered_iteritems

class Box(Object):

    def __init__(self, dims=None, units=None, resolution=None,
                 n_object=None, n_slice=None):
        Object.__init__(self, name='box',
                        dims=np.array(dims, dtype=np.float64),
                        units=units, resolution=resolution,
                        n_object=n_object, n_slice=n_slice)

        if not isinstance(self.n_slice, dict):
            self.n_slice = {'z' : int(self.n_slice)}

        self.init_trait('volume', np.prod(self.dims))

        self._axis_map = {'x' : [1, 2, 0], 'y' : [2, 0, 1], 'z' : [0, 1, 2]}

    def get_sizes(self, axis):
        """
        Get sizes (width, height) of a slice for the given axis.
        """
        sizes = self.dims[self._axis_map[axis][:-1]]
        return sizes

    def get_area(self, axis):
        """
        Get the area of a slice for the given axis.
        """
        return np.prod(self.get_sizes(axis))

    def get_pixel_sizes(self, axis):
        """
        Get sizes (width, height) in pixels of a slice for the given axis.
        """
        pixel_sizes = self.get_sizes(axis) / self.resolution
        return pixel_sizes
    
    def get_pixel_area(self, axis):
        """
        Get the area in pixels squared of a slice for the given axis.
        """
        return np.prod(self.get_pixel_sizes(axis))

    def get_points(self):
        """
        All points in the block.
        """
        for axis, num in ordered_iteritems(self.n_slice):
            am = self._axis_map[axis]

            shape = np.array((self.resolution[0], self.resolution[1], num))
            pb = np.zeros((3,), dtype=np.object)
            for ii in range(3):
                pb[am[ii]] = np.linspace(0, self.dims[ii], shape[ii])

            if num > 1:
                delta = pb[am[2]][1] - pb[am[2]][0]
            else:
                delta = 0.0
            x1, x2 = np.meshgrid(pb[am[0]], pb[am[1]])
            x1 = x1.ravel()
            x2 = x2.ravel()

            points = np.empty((x1.shape[0], 3), dtype=np.float64)
            points[:,am[0]] = x1
            points[:,am[1]] = x2

            yield pb, points, delta, num, axis, am
