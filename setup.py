#!/usr/bin/env python
"""gensei (generate serial images)

The script volume_slicer.py generates images corresponding to slices of a block
specimen filled with ellipsoids. The ellipsoids have known properties such as
the volume fraction within the block, dimensions (semiaxes a, b, c), volumes of
individual ellipsoids, a length-to-width ratio and number.

Although some basic parameters can be set using the command line options, it is
preferable to describe the objects and generation rules in an input file - see
'examples/basic.py', where four classes of objects are defined, and each
configuration item is described.

Its purpose is to verify results, obtained by various software for
reconstructing 3D data (e.g. a living tissue microstructure) given a set of 2D
slices using usually stereological methods, on a dataset with well-defined and
known properties.
"""

DOCLINES = __doc__.split("\n")

import os
import sys

VERSION = '2009.2'

CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: BSD License
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft :: Windows
"""

DOWNLOAD_URL = "http://code.google.com/p/gensei/wiki/Downloads?tm=2"

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'): os.remove('MANIFEST')

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('gensei')

    config.add_data_files(('gensei', ('README', 'LICENSE')))
    config.add_data_files(('../../../share/gensei/examples', 'examples/*'))

    config.get_version(version_variable=VERSION) # sets config.version
    print config
    return config

def setup_package():
    from numpy.distutils.core import setup
    from numpy.distutils.misc_util import Configuration

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    os.chdir(local_path)
    sys.path.insert(0, local_path)

    main_scripts = [
        'volume_slicer.py',
    ]

    try:
        setup(name = 'gensei',
              maintainer = "Robert Cimrman",
              maintainer_email = "cimrman3@ntc.zcu.cz",
              description = DOCLINES[0],
              long_description = "\n".join(DOCLINES[2:]),
              url = "http://http://code.google.com/p/gensei",
              download_url = DOWNLOAD_URL,
              license = 'BSD',
              classifiers = filter(None, CLASSIFIERS.split('\n')),
              platforms = ["Linux", "Mac OS-X", 'Windows'],
              scripts = main_scripts,
#              cmdclass = {'install_scripts' : install_scripts},
              configuration = configuration)
    finally:
        del sys.path[0]
        os.chdir(old_path)

    return

if __name__ == '__main__':
    setup_package()
