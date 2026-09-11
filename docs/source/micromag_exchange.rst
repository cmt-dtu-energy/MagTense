Exchange interaction
========================================

The exchange field is computed as

.. math::

    \mathbf{H}_\mathrm{exc} = \frac{2}{\mu_0 M_s}\, \nabla\cdot\left( A_0 \nabla \mathbf{m} \right),

which the solver evaluates as a sparse matrix - the *exchange operator* -
acting on each Cartesian component of :math:`\mathbf{m}`. How that matrix is
built depends on the grid type.

----------------------------------------
Exchange on a uniform grid
----------------------------------------

On a uniform grid the operator is the finite-difference Laplacian
:math:`\partial^2/\partial x^2 + \partial^2/\partial y^2 + \partial^2/\partial z^2`
with Neumann (free) boundary conditions, assembled directly in CSR format. A
direction with a single cell is dropped from the operator.

A spatially varying exchange stiffness is handled through the modified stencil
of Heistracher *et al.*, `Proposal for a micromagnetic standard problem: domain
wall pinning at phase boundaries
<https://doi.org/10.1016/j.jmmm.2021.168875>`_, so that a phase boundary
between two materials with different :math:`A_0` is treated correctly. This is
what the standard problem 6 example tests.

Periodic boundaries along any subset of the three directions are available on
the uniform grid, see :ref:`Periodic exchange boundaries`.

----------------------------------------
Exchange on unstructured meshes
----------------------------------------

For the ``unstructuredPrisms`` and ``tetrahedron`` grids the operator is built
from the mesh by ``computeDifferentialOperatorsFromMesh_DirectLap``, which is
given the faces, face normals, areas, volumes and interpolation stencils that a
mesh analysis has worked out beforehand. There is one analysis per grid type:
``UnstructuredMeshAnalysis.f90`` for the prisms works from the cell
centres and sizes and finds the neighbours geometrically.
``TetrahedralMeshAnalysis.f90`` for the tetrahedra works from the
connectivity - two tetrahedra are neighbours when they share three nodes, and
the face they share is the one opposite the node they do not. Both are run by
MagTense itself, so giving the mesh is all that is needed; see
:ref:`Tetrahedral grid` for how a tetrahedral mesh is passed in.

The method used to turn that into the operator is described in

`Poulsen, E. B., Insinga, A. R. and Bjørk, R.: Direct exchange calculation for
unstructured micromagnetic meshes, Journal of Magnetism and Magnetic Materials,
551, 169093, 2022 <https://www.sciencedirect.com/science/article/pii/S0304885322000671>`_

and in the ``DifferentialOperatorMeshes`` note in the ``documentation`` folder
of the repository. Three parameters control it.

**Exchange method** - ``exch_meth`` (Matlab
``setMicroMagExchMethod``, Python ``exch_meth``):

.. list-table::
   :widths: 30 8 62
   :header-rows: 1

   * - Name
     - Value
     - Description
   * - ``DirectLaplacianNeumann``
     - 1
     - Default. The face interpolation produces the *gradient* on each face and
       the second step assembles :math:`\nabla\cdot(A\nabla\phi)` directly,
       taking the local exchange stiffness into account.
   * - ``GGNeumann``
     - 2
     - Green-Gauss. The face interpolation produces :math:`\phi` on each face
       and the operator is formed as the product of two first-order operators.

**Interpolation stencil** - ``exch_intpn`` (Matlab
``setMicroMagExchInterpn``, Python ``exch_intpn``):

.. list-table::
   :widths: 30 8 62
   :header-rows: 1

   * - Name
     - Value
     - Description
   * - ``Extended``
     - 1
     - Default. The face value is interpolated from all elements that share at
       least a vertex with the face, which gives a wide and robust stencil on
       an irregular mesh.
   * - ``Compact``
     - 2
     - Only the elements immediately adjacent to the face are used. The solver
       prints ``Warning: untested method: compact`` when this is selected.

**Weighting** - ``exch_weigh`` (default ``8.0``) is the exponent of the
inverse-distance weighting used in the face interpolation. A large exponent
makes the interpolation approach a nearest-neighbour average.

----------------------------------------
Exchange across a material interface
----------------------------------------

Where two cells with different exchange stiffness share a face, the exchange
that couples them is by default the harmonic mean of the two values,

.. math::

    A_\mathrm{face} = \frac{2 A_1 A_2}{A_1 + A_2} ,

on both grid types. The harmonic mean is the natural choice when the two
materials are simply different regions of one exchange-coupled body, but it is
a modelling assumption rather than a measured quantity: the exchange across a
real phase boundary depends on the interface itself and is often weaker than
either bulk value. It can therefore be specified directly.

The interface value is given per **pair of materials**, not as a single number.
Each cell is labelled with a material index, and a table gives the exchange to
use for each pair of labels:

.. list-table::
   :widths: 22 22 12 44
   :header-rows: 1

   * - Matlab
     - Python
     - Unit
     - Description
   * - ``n_phase``
     - ``n_phase``
     - \-
     - Number of distinct materials. The default of 1 disables the feature.
   * - ``phase_id``
     - ``phase_id``
     - \-
     - Material index of every cell, an ``(ntot,1)`` array with entries in
       ``1..n_phase``.
   * - ``A_int``
     - ``A_int``
     - J/m
     - Symmetric ``(n_phase, n_phase)`` table of interface exchange values. A
       negative entry means *use the harmonic mean* for that pair.

Keying the value on the pair of materials keeps it well defined for any number
of materials. Every internal face lies between exactly two cells, so even where
three or more materials meet, each individual face is still an unambiguous
two-material interface.

Only pairs of *different* materials are consulted: two cells of the same
material always use the harmonic mean, so the diagonal of ``A_int`` is never
read. Because a negative entry falls back to the harmonic mean, a table can
override only the pairs of interest and leave the rest at the default.

In Matlab:

.. code-block:: matlab

    problem.n_phase  = int32(2);
    problem.phase_id = phase_id;            % 1 or 2 for every tile
    problem.A_int    = [   -1, 3e-12; ...   % 3 pJ/m across the 1-2 interface
                        3e-12,    -1];      % diagonal unused, negative = default

and in Python, either as constructor arguments or afterwards:

.. code-block:: python

    problem = MicromagProblem(..., phase_id=phase_id, A_int=A_int)
    # or
    problem.set_interface_exchange(phase_id, A_int)

The table must be symmetric. An asymmetric table would make the exchange across
a face depend on which of the two cells it is asked from, which is not a
physical operator, so it is rejected with an error rather than silently
symmetrised. Out-of-range or wrongly sized ``phase_id`` arrays are rejected in
the same way.

.. note::
   Setting ``A_int`` to the harmonic mean of the two materials reproduces the
   default operator exactly, to the last bit, on both grid types. That makes it
   easy to confirm a model has been set up as intended before changing the
   interface value.

.. note::
   A cell with :math:`A_0 = 0` is not magnetic and carries no exchange, so a
   face touching one stays uncoupled whatever the table says. This is what the
   harmonic mean does by itself, since it vanishes as soon as either side is
   zero, and it holds on both grid types. Such a cell also takes no part in the
   face interpolation on an unstructured mesh, so nothing in the mesh depends on
   its magnetisation.

----------------------------------------
Supplying the exchange matrix directly
----------------------------------------

The exchange operator may be computed outside MagTense and passed in, which
skips the mesh analysis entirely. The matrix is given in coordinate (COO)
format. In Matlab:

.. code-block:: matlab

    problem = problem.setExchangeMatrixCOO( nrows, ncols, rows, cols, values );

which sets ``passExch = 1`` together with ``exch_nrow``, ``exch_ncol``,
``exch_nval``, ``exch_rows``, ``exch_cols`` and ``exch_val``. In Python the
same fields are constructor arguments (``exch_val``, ``exch_rows``,
``exch_cols``, ``exch_nval``, ``exch_nrow``, ``exch_ncols``) together with
``passexch=1``.

A matrix passed in this way carries its own boundary conditions, so the
``exchPBC`` setting is not applied to it.

.. note::
   Row and column indices are 1-based, as in Fortran and Matlab.

----------------------------------------
Returning the exchange matrix
----------------------------------------

The assembled exchange matrix is returned to the caller, which is useful both
for inspection and for reusing the same operator in a later run. In Matlab it
is part of the ``GridInfo`` output as ``ExchMat_r``, ``ExchMat_c``,
``ExchMat_v``, ``ExchMat_nr`` and ``ExchMat_nc``. In Python the same arrays are
elements 7-12 of the result list returned by ``run_simulation``.

Because Python has to preallocate the arrays that the values are copied into,
the parameter ``exch_presize`` sets their length as ``exch_presize * ntot``.
The default of 12 is enough for a uniform grid; unstructured meshes with wide
stencils need more, and the solver prints the required value if the array is
too small:

.. code-block:: python

    problem.exch_presize = 28
