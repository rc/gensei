objects = {
    'class 1' : {
        'kind' : 'ellipsoid',
        'color' : 'r',
        'fraction' : 0.02,
        'length_to_width' : 8.0,
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

box = {
    'dims' : (10, 10, 10),
    'units' : 'mm',
    'resolution' : (300, 300),
    'n_object' : 30,
    'n_slice' : {'z' : 21, 'x' : 11, 'y' : 5},
    }

options = {
    'output_format' : 'png',
    'timeout' : 5.0,
}
