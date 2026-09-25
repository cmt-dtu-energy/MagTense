Calculations
========================================

In MagTense, a magnetostatic and a :ref:`micromagnetism <micromagnetism>` calculation 
framework is available. The basic pipeline for a magnetostatic 
calculation consists of three parts:

* MagTiles
* Evaluation points
* State function

As output, the three-dimensional H-field vector in 
the evaluation points is returned.

Examples of how to calculate magnetostatic and micromagnetic
problems can be found for `Matlab <https://github.com/cmt-dtu-energy/MagTense/tree/master/matlab/examples>`_
and `python <https://github.com/cmt-dtu-energy/MagTense/tree/master/python/examples>`_,
and the micromagnetic ones are described in
:ref:`Micromagnetic examples and validation`.

The parameters of a micromagnetic problem are of a different kind than the
MagTile properties below, and are documented separately in
:ref:`Micromagnetic parameter reference`.


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
Applied field
========================================

An external field that is the same at every point, such as the field of a
large electromagnet or a Helmholtz coil, is entered as a tile of type 102.
This tile is **not a geometry**. It has no size, position or orientation; it
is a field source, and its ``M`` vector holds the applied field
:math:`\mathbf{H}_\mathrm{app}` in A/m. Like the planar coil it is never
updated by the iteration. It acts in two places:

* it is part of the field that magnetizes every other tile in the iteration,
  which is what lets a soft magnet alone in space acquire a magnetization, and
  it enters before the self-consistent field of a soft tile is formed, so the
  permeability of the tile sees the total field;
* it is added to the H-field returned at the evaluation points, so the
  returned field is the total field, applied plus that of the tiles.

The value is an H field. The source is by definition outside all tiles, where
:math:`\mathbf{B}/\mu_0` and :math:`\mathbf{H}` coincide, so a field known as
:math:`\mathbf{B}` in tesla is divided by :math:`\mu_0` before it is entered.
A field evaluated inside a magnetized body is not an applied field: for the
cylindrical tile MagTense computes :math:`\mathbf{B}/\mu_0` inside the tile
and subtracts :math:`\mathbf{M}` to get :math:`\mathbf{H}` there, and that
correction belongs to the tile's own field, not to this one. The field of other
MagTense tiles is not entered this way either: put those tiles in the same
problem, with ``includeInIteration`` set to zero if they are permanent magnets
whose magnetization should stay fixed, and the iteration uses their field
directly, correction included.

In Python the source is appended with ``tiles.add_uniform_field(H_app)``, in
Matlab with ``tile.setMagTileType('Uniformfield')`` and ``tile.M = H_app``, and
in the standalone input file with tile type 102 and the field in the ``M``
line. A soft sphere with constant permeability :math:`\mu_r` in such a field is
the standard check: its magnetization is
:math:`\mathbf{M} = 3 (\mu_r - 1)/(\mu_r + 2)\,\mathbf{H}_\mathrm{app}`, and
`soft_sphere_in_uniform_field.py <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/examples/magnetostatics/soft_sphere_in_uniform_field.py>`_
and its Matlab counterpart
``Example_005_soft_sphere_in_uniform_field`` reproduce it.

A field that differs from tile to tile cannot be entered this way; that is a
separate, per-tile applied field, which is not available yet.

========================================
H-field
========================================

::

    real,dimension(n_ele,3) :: H

The three-dimensional H-field is returned in the unit :math:`[A/m]`
and its vector is given with respect to the global coordinate system. It is
the total field: the field of the tiles plus any :ref:`Applied field`.
