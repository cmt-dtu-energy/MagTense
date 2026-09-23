Grids and meshes
========================================

The micromagnetic model supports three grid types. The grid type is selected
with ``problem.setMicroMagGridType(...)`` in Matlab and with the ``grid_type``
argument in Python, and it decides which of the geometry parameters are read.

.. list-table::
   :widths: 22 8 70
   :header-rows: 1

   * - Grid type
     - Value
     - Required geometry input
   * - ``uniform``
     - 1
     - ``grid_n`` (``res`` in Python) and ``grid_L``. The grid is generated
       internally.
   * - ``tetrahedron``
     - 2
     - ``grid_pts``, ``grid_nod``, ``grid_ele`` and ``grid_nnod``.
   * - ``unstructuredPrisms``
     - 3
     - ``grid_pts`` and ``grid_abc``.

----------------------------------------
Uniform grid
----------------------------------------

This is the default grid type, and by far the simplest. The domain is a
rectangular box of size ``grid_L`` :math:`= (L_x, L_y, L_z)` divided into
``grid_n`` :math:`= (n_x, n_y, n_z)` identical rectangular prism cells, so that
:math:`\Delta x = L_x/n_x` and similarly for :math:`y` and :math:`z`. The box is
centred on the origin, i.e. it spans :math:`[-L_x/2, L_x/2]` and so on, and the
cell centres are placed at the centre of each cell. A direction with a single
cell is placed at zero and given the full length of the domain as its cell
size.

The cells are ordered with :math:`x` running fastest and :math:`z` slowest,

.. math::

    \mathrm{ind} = i + (j-1)\,n_x + (k-1)\,n_x n_y .

The same ordering applies to every array of single-cell properties (``Ms``, ``K0``, ``u_ea``,
``m0``, ``temperature``, ...) and to the returned magnetization.

For a uniform grid the exchange operator is built as a finite-difference
Laplacian, see :ref:`Exchange interaction`.

----------------------------------------
Unstructured prisms
----------------------------------------

The unstructuredPrisms grid is a mesh of rectangular prisms of varying size. 
The prisms are all aligned in the sense their faces are perpendicular to the global xyz-axes, 
so every face is either parallel or perpendicular to every other.
The mesh is specified through

* ``grid_pts`` - an ``(ntot,3)`` array with the centre of every prism;
* ``grid_abc`` - an ``(ntot,3)`` array with the full side lengths
  :math:`(a,b,c)` of every prism.

The problem resolution must be set to ``[ntot,1,1]``, since the number of cells
is what the solver uses, not a rectangular subdivision:

.. code-block:: matlab

    resolution = [length(mesh.pos_out) 1 1];
    problem = DefaultMicroMagProblem( resolution(1), resolution(2), resolution(3) );
    problem = problem.setMicroMagGridType('unstructuredPrisms');
    problem.grid_pts = mesh.pos_out;
    problem.grid_abc = mesh.dims_out;

.. code-block:: python

    problem = MicromagProblem(
        res=(len(grid_pts), 1, 1),
        grid_type="unstructuredPrisms",
        grid_pts=grid_pts,
        grid_abc=grid_abc,
    )

For this grid type the mesh is analysed by
``CartesianUnstructuredMeshAnalysis``, which finds the faces, face normals,
face areas, cell volumes and the interpolation stencils that the exchange
operator is built from. The result is returned to Matlab as the ``GridInfo``
struct. The mesh analysis is also where periodic exchange boundaries are
resolved, see :ref:`Periodic exchange boundaries`.

An unstructured mesh may be plotted from Matlab with
``cartesianUnstructuredMeshPlot``.

----------------------------------------
Tetrahedral grid
----------------------------------------

A tetrahedral mesh is specified through the node coordinates ``grid_nod``, the
number of nodes ``grid_nnod``, the element definitions ``grid_ele`` (the four
node indices of each element) and the element centres ``grid_pts``. As for the
unstructured prisms, the resolution is ``[ntot,1,1]``.

Giving the mesh is all that is required. MagTense analyses it and builds the
exchange operator itself, exactly as it does for a grid of unstructured prisms
given its centres and sizes. In MATLAB the whole grid is set up by

.. code-block:: matlab

   problem = problem.setMicroMagGridTetrahedron(model.Mesh.Nodes, model.Mesh.Elements);

which fills in ``grid_nod``, ``grid_ele``, ``grid_nnod`` and the element centres
``grid_pts``, and in Python by passing ``grid_nod`` and ``grid_ele`` to
``MicromagProblem``, where ``grid_nnod`` and ``grid_pts`` likewise default to the
mesh. A quadratic mesh is accepted and analysed as its linear counterpart.

The mesh analysis is in ``TetrahedralMeshAnalysis.f90``. It works from the
connectivity rather than from the geometry: two tetrahedra are neighbours when
they share three nodes, and the face they share is the one opposite the node
that they do not.

Periodic exchange boundaries are requested through ``exchPBC``. They are linked
either by identifying the nodes on the two boundary planes, when the mesh is
periodic along that direction, or by mortar coupling the boundary faces when it
is not, with the choice made per direction and a warning printed when the
fallback is used. See :ref:`Periodic boundaries on a tetrahedral mesh` for what
each mechanism does and what the fallback costs.

The demagnetization tensor is computed from the analytical tetrahedron tensor.
Note the following limitations of this grid type:

* Averaging of the demagnetization tensor over the receiving cell (``N_ave``)
  is not supported and a warning is printed if it is requested.
* The shape correction is disabled, because a tetrahedral mesh carries no
  per-element size array from which the occupied volume fraction can be found.

----------------------------------------------
Passing an externally computed exchange matrix
----------------------------------------------

For any grid type the exchange operator may be computed outside MagTense and
passed in directly, which bypasses the mesh analysis entirely. This is
described in :ref:`Supplying the exchange matrix directly`.
