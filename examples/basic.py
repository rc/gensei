# This is an example file with advanced settings. It can be called in ipython
# shell as: %run volume_slicer.py examples/myexample.py

#--------Settings of the properties of objects-----------
objects = {
    # objects of the first class
    'class 1' : {
        # geometry of the objects
        'kind' : 'ellipsoid',
        # color in RGB space, or one of
        # 'b'         blue
        # 'g'         green
        # 'r'         red
        # 'c'         cyan
        # 'm'         magenta
        # 'y'         yellow
        # 'w'         white
        'color' : 'r',
        # volume fraction of the objects of this class
        'fraction' : 0.02,
        # length-to-width ratio
        'length_to_width' : 8.0,
        # (0,1) allows to change the geometry of the objects when they can not
        # fit the size of the box in order to keep the desirable volume fraction
        # and number of objects;
        # another option is 'reduce_to_fit' : {'fraction' : 0.9} which tries to
        # change the volume fraction in order to keep the number and geometry
        # of objects
        'reduce_to_fit' : {'length_to_width' : 0.9},
        'centre' : 'random',
        'rot_axis' : 'random',
        'rot_angle': 'random',
    },
    'class 2' : {
        'kind' : 'ellipsoid',
        'color' : (0.1, 0.2, 0.7),
        'fraction' : 0.01,
        'length_to_width' : 1.0,
        'reduce_to_fit' : {'fraction' : 0.9},
        'centre' : 'random',
        'rot_axis' : 'random',
        'rot_angle': 'random',
    },
    'class 3' : {
        'kind' : 'ellipsoid',
        'color' : (0.7, 1.0, 0.0),
        'fraction' : 0.001,
        'length_to_width' : 30.0,
        'reduce_to_fit' : {'fraction' : 0.9},
        'centre' : 'random',
        'rot_axis' : 'random',
        'rot_angle': 'random',
    },
    'class 4' : {
        'kind' : 'ellipsoid',
        'color' : (0.7, 1.0, 1.0),
        'fraction' : 0.0001,
        'length_to_width' : 1.0,
        'reduce_to_fit' : {'fraction' : 0.9},
        'centre' : 'random',
        'rot_axis' : 'random',
        'rot_angle': 'random',
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
    # 'n_object' : {'class 1' : 5, 'class 2' : 3, 'class 3' : 2, 'class 4' : 1},
    'n_object' : 30,
    # number of slices generated - a dictionary of numbers for the
    # directions perpendicular to the Z, X, and Y axis, i.e. in the XY, YZ, and
    # XY planes, or a single integer:
    # 'n_slice' : 10, is equivalent to
    # 'n_slice' : {'z' : 10},
    'n_slice' : {'z' : 21, 'x' : 11, 'y' : 5},
    }
#--------End of settings of the properties of objects-----------

#--------General options-----------
options = {
    # output file format supported by the matplotlib backend used
    'output_format' : 'png',
    # timeout in seconds to place more ellipsoids in to the box
    'timeout' : 5.0,
}
