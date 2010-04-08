.. gensei documentation master file, created by
   sphinx-quickstart on Tue Jun 16 18:45:54 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to gensei's documentation!
==================================
*gensei* = generator of serial images, or a karate style (gen = mysterious, sei =
control, system), meaning that the methods of quantitative histology (or
stereology) are often a mysterious process for us :)

web: http://code.google.com/p/gensei/

code: http://github.com/rc/gensei/

Get the latest sources::

  $ git clone git://github.com/rc/gensei.git
  $ cd gensei

and read the ``README`` for  installation instructions.

Abstract
--------

The software is used for generating series of two-dimensional sections through
a volume filled with ellipsoids, cylinders and fibres with a priori known
statistical properties. Arbitrary series of images are produced in order to
represent models of histological structures, such as cells and fibres of the
extracellular matrix forming biological tissues. These images represent either
serial consecutive sections or sections sampled in a systematic uniform
manner. The user can set the number and types of classes of three-dimensional
objects (ellipsoids, cylinders), dimensions of the reference volume, size of
the generated images, the number of objects and sections to be generated,
length-to-width ratio and volume fraction of the objects. Three perpendicular
series of sections through the same reference volume and the same objects may
be generated, thus testing the isotropy of the objects. The output of the
simulation is stored as three series of image files together with a detailed
report containing the statistical properties of the generated objects,
i.e. number of objects per class, volume fraction of the objects within the
reference volume, surface density of the objects, number of intersections
between the objects and the sections, and rotation of the objects around
uniquely determined spatial vectors. By comparing the true and section-based
estimates of the volume, surface and length, the settings of various
stereological grids used for estimation can be tested. The user can
e.g. simulate the effect of section thickness on the stereological estimate of
volume, surface and length. The resulting images have a high contrast so the
objects can be segmented.

Contents:

.. toctree::
   :maxdepth: 2

   users_guide
   developer_guide

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

