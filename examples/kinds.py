"""
Demonstration of using several kinds of objects.
"""
import numpy as np

#--------Settings of the properties of objects-----------
objects = {
    'class 1' : {
        'kind' : 'ellipsoid',
        'color' : 'r',
        'fraction' : 0.1,
        'length_to_width' : 2.0,
        'reduce_to_fit' : {'length_to_width' : 0.9},
        'centre' : 'random',
        # Ellipsoids should be approximately aligned with the x axis.
        'rot_axis': [('normal', 0.0, 0.1), 1.0, ('random', -0.2, 0.2)],
        'rot_angle': np.pi/2,
    },
    'class 2' : {
        'kind' : 'cylinder',
        'color' : (0.1, 0.2, 0.7),
        'fraction' : 0.01,
        'length_to_width' : 40.0,
        'reduce_to_fit' : {'fraction' : 0.9},
        'centre' : 'random',
        # Cylinders should have uniform random distribution of directions. 
        'direction' : 'random direction',
    },
}
#--------End of settings of the properties of objects-----------

#--------Settings of the properties of the box-----------
box = {
    # dimensions of the box
    'dims' : (10, 10, 10),
    # arbitrary units 
    'units' : 'mm',
    # size of the generated image in pixels
    'resolution' : (300, 300),
    # number of objects to be generated within the box, either a number, or a
    # per class dictionary, for example:
    'n_object' : {'class 1' : 10, 'class 2' : 50},
    # number of slices generated - a dictionary of numbers for the
    # directions perpendicular to the Z, X, and Y axis, i.e. in the XY, YZ, and
    # XY planes, or a single integer.
    'n_slice' : {'z' : 21, 'x' : 21, 'y' : 21},
    }
#--------End of settings of the properties of objects-----------

#--------General options-----------
options = {
    # output file format supported by the matplotlib backend used
    'output_format' : 'png',
    # timeout in seconds to place more ellipsoids in to the box
    'timeout' : 5.0,
}
