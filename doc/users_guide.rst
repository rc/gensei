User's Guide
============

Gensei is a software for generating serial sections through objects with known
statistical properties and simulating stereological quantification in 3D.

Introduction
------------

The purpose of gensei is to have a simple, open source tool to generate images
of slices of a volume filled with objects (e.g. ellipsoids) with known
statistical and other properties (e.g. total volume, volume fractions, number
of objects), to verify information obtained by the stereological methods. 

Stereology
----------

Modern stereological methods are used for quantitative and qualitative
description of objects, such as metals, stones, biological tissues, and their
components. The stereology was defined by Weibel (1981) as "a body of
mathematical methods relating three-dimensional measurements obtainable on
sections of the structure". The principles of stereology have one common
feature: They make use of random tissue sections of negligible thickness, on
which a test system is placed at random. This test system is made up of probes
in the form of points, lines, or planar areas, depending on the information
sought, see Weibel (1966).

The earlier decade of application of stereological methods on biological
specimens (liver cells, muscles, and other systems) is described in Weibel
(1966, 1981). From its beginning the stereology reached a great advancement and
presently belongs to the most used methods for unbiased quantification of
number, length, surface area and volume of specimen attributes of samples of
various sizes and structures in the fields of biology, metallography, and
petrography and may reveal important information about the function and
organization of the parts being studied (organs, tissues, grains, etc.).

The stereological techniques are in detail described in following works:
Gundersen and Jensen (1987), Gundersen et al. (1988a), Gundersen et
al. (1988b), Howard and Reed (1998), Mandarim-de-Lacerda (2003), Mouton (2001),
Nyengaard (1999), Nyengaard and Gundersen (2006), Russ and Dehoff (1999).

References
^^^^^^^^^^
Gundersen and Jensen (1987): Gundersen, H. J. G. and Jensen, E. B. (1987). The
efficiency of systematic sampling in stereology and its prediction. Journal of
Microscopy, 147:229-263. (http://www.ncbi.nlm.nih.gov/pubmed/3430576)

Gundersen et al. (1988a): Gundersen, H. J. G., Bendtsen, T. F., Korbo, L.,
Marcussen, N., Møller, A., Nielsen, K., Nyengaard, J. R., Pakkenberg, B.,
Sørensen, F. B., Vesterby, A., and West, M. J. (1988b). Some new, simple and
efficient stereological methods and their use in pathological research and
diagnosis. Acta Pathologica, Microbiologica et Immunologica Scandinavica,
96:379-394. (http://www.ncbi.nlm.nih.gov/pubmed/3288247)

Gundersen et al. (1988b): Gundersen, H. J. G., Bagger, P., Bendtsen, T. F.,
Evans, S. M., Korbo, L., Marcussen, N., Møller, A., Nielsen, K., Nyengaard,
J. R., Pakkenberg, B., Sørensen, F. B., Vesterby, A., and West,
M. J. (1988a). The new stereological tools: Disector, fractionator, nucleator
and point sampled intercepts and their use in pathological research and
diagnosis. Acta Pathologica, Microbiologica et Immunologica Scandinavica,
96:857-881. (http://www.ncbi.nlm.nih.gov/pubmed/3056461)

Howard and Reed (1998): Howard, C. V. and Reed, M. G. (1998). Unbiased Stereology: Three Dimensional Measurement in Microscopy. Royal Microscopical Society, Microscopy Handbook Series No. 41. Springer, New York.

Mandarim-de-Lacerda (2003): Mandarim-de-Lacerda, C. A. (2003). Stereological
tools in biomedical research. Annals of the Brazilian Academy of Sciences,
75(4):469-486.(http://www.ncbi.nlm.nih.gov/pubmed/14605681)

Mouton (2001): Mouton, P. R. (2001). Principles and Practices of Unbiased
Stereology: An Introduction for Bioscientists. The Johns Hopkins University
Press, USA, Baltimore.

Nyengaard (1999): Nyengaard, J. R. (1999). Stereologic methods and their
application in kidney research. Journal of the American Society of Nephrology,
10:1100-1123. (http://jasn.asnjournals.org/cgi/content/full/10/5/1100)

Nyengaard and Gundersen (2006): Nyengaard, J. R. and Gundersen,
H. J. G. (2006). Sampling for stereology in lungs. European Respiratory Review,
15(101):107-114. (http://err.ersjournals.com/cgi/content/abstract/15/101/107)

Russ and Dehoff (1999): Russ, J. C. and Dehoff, R. T. (1999). Practical
Stereology. Plenum Press, New York.

Weibel (1966): Weibel, E. R., Kistler, G. S., and
Scherle,W. F. (1966). Practical stereological methods for morphometric
cytology. The Journal of Cell Biology,
30:23-38. (http://www.ncbi.nlm.nih.gov/pubmed/5338131)

Weibel (1981): Weibel, E. R. (1981). Stereological methods in cell biology:
Where are we - where are we going? The Journal of Histochemistry and
Cytochemistry,
29(9):1043-1052. (http://www.jhc.org/cgi/reprint/29/9/1043)

Basic usage
-----------

The script *volume_slicer.py* can be used to place objects of known (given)
geometrical and statistical properties into a (virtual) 3D box and then
generate the slices. It can be called as::

    Usage: volume_slicer.py [options] [filename]                                 

    If an input file is given, the object class options have no effect.

    Default option values do _not_ override the input file options.


    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit       
      -o filename           basename of output file(s) [default: ./slices/slice]
      -f format, --format=format
                            output file format (supported by the matplotlib
                            backend used) [default: png]
      -n int, --n-slice=int
                            number of slices to generate [default: 21]
      -d dims, --dims=dims  dimensions of specimen in units given by --units
                            [default: (10, 10, 10)]
      -u units, --units=units
                            length units to use [default: mm]
      -r resolution, --resolution=resolution
                            figure resolution [default: 600x600]
      --fraction=float      volume fraction of objects [default: 0.1]
      --fraction-reduction=float
                            volume fraction reduction factor [default: 0.9]
      --length-to-width=float
                            length-to-width ratio of objects [default: 8.0]
      --n-object=int        number of objects [default: 10]
      -t float, --timeout=float
                            timeout in seconds for attempts to place more
                            ellipsoids into the block [default: 5.0]
      --no-pauses           do not wait for a key press between fitting,
                            generation and slicing phases (= may overwrite
                            previous slices without telling!)

Further options can be set in an input file, see below.

Example input file
^^^^^^^^^^^^^^^^^^

This file defines four classes of ellipsoids, the box and other necessary
options.

.. literalinclude:: ../examples/basic.py

Examples
--------

One class of objects
^^^^^^^^^^^^^^^^^^^^

.. figure:: images/examples/slice.05.00_70.png
    :width: 50%
    :align: center
    :alt: One of generated slices.

    One of generated slices.

We assessed the volume fractions of objects (in our case red ellipsoids
obtained by gensei) :math:`V_V(objects,ref)` using one of the stereological
methods made on two-dimensional serial sections and points test probe:

.. math::

    V_V(objects,ref) = V(objects)/V(ref) = \sum{i=1}^m P_i/P_T,

where objects represents ellipsoids, :math:`V(ref)` is the volume of reference
space, in this case the volume of whole stack of serial
sections. :math:`V(objects)` and :math:`V(ref)` were estimated using the
Cavalieri principle according to equation:

.. math::

    estV(c) = T a/p \sum{i=1}^m P_i,

where c corresponds to the ellipsoids and the reference space, respectively,
a/p is the known area associated with each point of the test system, :math:`T=`
0.140845 :math:`\mu{\rm m}` is the distance between two subsequent sections and
:math:`P_i` is the number of points landing within the relevant component
transect on the i-th section. Because the parameters of test point
probe of objects and reference space were similar, the equation
:math:`V_V(objects,ref)` could be simplified to:

.. math::

    V_V(objects,ref) = \sum{i=1}^m P_i/P_T,

where :math:`P_i` was a number of points of randomly translated point grid
landing within the object transect on the i-th section, :math:`P_T` was
a number of all points of the test system.

The examples of point test system and counting the intersections using the
module PointGrid of the software Ellipse3D (http://www.ellipse.sk/index.htm,
ViDiTo Systems, KoĹĄice, Slovak Republic) follow.

.. figure:: images/examples/slice.05.00_70_point_grid.png
    :width: 50%
    :align: center
    :alt: Point grid.

    Point grid.

.. figure:: images/examples/slice.05.00_70_points.png
    :width: 50%
    :align: center
    :alt: Points.

    Points.

The color channel serial images obtained by gensei was changed to gray scale
using the software IrfanView (http://www.irfanview.com, Irfan Skiljan, Vienna
University of Technology, Vienna, Austria). The software Amira
(http://www.tgs.com, TGS Inc., Massachusetts, USA) was then used for the
three-dimensional reconstruction. The stack of the grey-scale images was loaded
in the Amira and the threshold of objects was determined by using histogram
curves and the slices were then joined together. The resultant non-smoothing
model of individual objects is shown in following figure.

.. figure:: images/examples/snapshot.png
    :width: 50%
    :align: center
    :alt: Reconstructed 3D structure (not smoothed).

    Reconstructed 3D structure (not smoothed).

Multiple classes of objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Currently, gensei (git version) supports generating images of more classes. For
each class, the geometrical parameters, volume fraction, number of objects, and
colour can be set separately. Moreover, the user can choose to generate series
of sections in three perpendicular cutting planes. The following image contains
objects of three classes (red, green, blue).

.. figure:: images/examples/three_classes/slice.x.14.35_90.png
    :width: 50%
    :align: center
    :alt: Objects of three classes.
    
    Objects of three classes.

Such a contrast image can be easily segmented. The following segmentation was
performed with the Thresh module of the software Ellipse 3D (Vidito, KoĹĄice,
Slovak Republic):

.. figure:: images/examples/three_classes/s_014.png
    :width: 50%
    :align: center
    :alt: Segmentation.

    Segmentation.

The circumferences of all the object profiles can be saved and their
geometrical properties (e.g., circularity, circumference, etc.) can be
assessed. Moreover, the objects can be reconstructed in three dimensions after
calibration and the contours can be visualized as follows (module Contours,
software Ellipse 3D):

.. figure:: images/examples/three_classes/snap2.png 
    :width: 50%
    :align: center
    :alt: Contours.

    Contours.

The surface of the objects can be reconstructed. The reconstruction gives
usually best results when the objects were cut along their long axis. The
following visualization was obtained with the Surface module of the software
Ellipse 3D:

.. figure:: images/examples/three_classes/snap1.png
    :width: 50%
    :align: center
    :alt: Surface.

    Surface.

The volume fraction and the surface of the objects is known and saved together
with the image series. Therefore it can be used e.g. for training in
histology - the trainee can compare the results of estimating the volume
fraction (using the Cavalieri principle) either with the true volume fraction
of the whole simulated objects or with the volume fraction based on the series
of sections (the latter is usually lower as not all parts of the objects appear
on the sections). Another application of the images is testing the settings of
various stereological grids used for estimation of volume. When adjusting the
geometrical parameters of stereological grids, it is usually desirable to find
the lowest number of intersections between the objects and the grid, which
still yields an acceptable amount of error. This minimum number of
intersections depends on geometry of the objects, and the desired error (for
details, see the nomogram published by Gundersen and Jensen, 1987). The user
can easily check, whether the current settings give a good estimate of the true
(known) volume. In the following example, the volume of red objects was
quantified with a rectangular point grid (counted intersections are yellow,
module LineSystem?, software Ellipse 3D):

.. figure:: images/examples/three_classes/snap3.png
    :width: 50%
    :align: center
    :alt: Red objects were quantified with a rectangular point grid.

    Red objects were quantified with a rectangular point grid.
