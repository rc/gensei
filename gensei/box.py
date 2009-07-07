from gensei.base import np, Object, pause

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

    def get_pixel_sizes(self, axis):
        indx = {'x' : [1, 2], 'y' : [0, 2], 'z' : [0, 1]}
        pixel_sizes = self.dims[indx[axis]] / self.resolution
        return pixel_sizes
    
    def get_pixel_area(self, axis):
        return np.prod(self.get_pixel_sizes(axis))
