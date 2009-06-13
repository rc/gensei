from gensei.base import np, Object, pause

class Box(Object):

    def __init__(self, dims=None, units=None, resolution=None,
                 n_object=None, n_slice=None):
        Object.__init__(self, dims=np.array(dims, dtype=np.float64),
                        units=units, resolution=resolution,
                        n_object=n_object, n_slice=n_slice)

        self.init_trait('volume', np.prod(self.dims))
