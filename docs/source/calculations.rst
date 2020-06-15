Calculations
========================================

In MagTense, a magnetostatic and a :ref:`micromagnetism <micromagnetism>` calculation 
framework is available. The basic pipeline for a magnetostatic 
calculation consists of three parts:

* MagTiles
* Evaluation points
* State function

As result, the three-dimensional H-field vector in 
the evaluation points is returned.

Examples of how to calculate magnetostatic and micromagnetic 
problems can be found for `Matlab <https://github.com/cmt-dtu-energy/MagTense/tree/master/matlab/examples>`_ 
and `python <https://github.com/cmt-dtu-energy/MagTense/tree/master/python/examples>`_.


========================================
MagTile
========================================

::
    
    type(MagTile),dimension(n_tiles) :: tiles

MagTiles are the basic structure of MagTense, where **n_tiles** 
is the number of given magnetic tiles.
They can be specified with several parameters corresponding
to the following properties:

.. toctree::
   :maxdepth: 2

   geometry
   magnetization
   other_parameters

========================================
Evaluation points
========================================

::
    
    real,dimension(n_ele,3) :: pts

The three-dimensional evaluation points are defined with respect 
to the global coordinate system, where **n_ele** is the number
of given points.

========================================
State function
========================================

::

    type(MagStateFunction),dimension(n_stf) :: stateFunction 

The state functions (hysteresis loops) for different materials 
and temperatures can be given and used if the given MagTile 
is representing a soft magnet, where **n_stf** is the number of
given state functions.

========================================
H-field
========================================

::

    real,dimension(n_ele,3) :: H

The three-dimensional H-field is returned in the unit :math:`[A/m]`
and its vector is given with respect to the global coordinate system.