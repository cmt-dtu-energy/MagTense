Periodic boundaries and macrogeometry
========================================

A micromagnetic simulation resolves a domain of sub-mm size, while real
samples are macroscopic. MagTense bridges this gap in length scale with three
separate methods, which can be combined:

#. a **macrogeometry** that repeats the simulated domain when the
   demagnetization tensor is built;
#. a **shape correction** that replaces the demagnetization field from the average magnetisation 
   with that of the real sample shape;
#. **periodic exchange boundaries** that couple the cells at opposite ends of
   the domain.

The periodic demgnetization field is verified against an analytical result from 
Durhuus *et al.*, "Exact demagnetisation field for periodic
one-dimensional array of rectangular prisms" (see :ref:`Publications`). The
implementation is validated by the ``macrogeometry_PBC_test``,
``macrogeometry_supercell_test`` and ``shape_correction_test`` examples that
ship with MagTense in both Matlab and Python.

That average magnetisation is enough to fully account for shape effects is proven 
in Durhuus *et al.*, "Including sample shape in micromagnetics with 3D periodic 
boundary conditions" (see :ref:`Publications`). 
The idea is that the field from the average magnetisation :math:`H_\mathrm{avg}` 
is longer-ranged than the field from the non-uniform remainder, :math:`H_\mathrm{rest}`. 
Consequently, :math:`H_\mathrm{avg}` depends on the macroscopic distribution of magnetic 
material (sample-shape), while :math:`H_\mathrm{rest}` just depends on the local environment.
Thus for simulating the bulk of a sample, i.e. a simulation domain far from the surface,
one can use the macrogeometry method for the local environment, and the shape correction
field for the distant sample regions, including sample-shape effects.

----------------------------------------
Macrogeometry
----------------------------------------

The macrogeometry is a regular array of copies of the simulated domain. Its
effect is included when the demagnetization tensor is built: for every copy the
tensor is evaluated with the evaluation points shifted by the copy offset, and
the contributions are summed. The magnetization of every copy is by
construction identical to that of the simulated domain, so the cost is one
extra tensor evaluation per copy and no extra memory. Since every copy
carries the same material properties and magnetisation, the simulated domain
represents the bulk sample, and it is chosen large enough to contain the
features of interest, e.g. long domain walls.

.. list-table::
   :widths: 22 22 56
   :header-rows: 1

   * - Matlab
     - Python
     - Description
   * - ``n_macro``
     - ``n_macro``
     - Number of copies on **each** side of the simulated domain along x, y and
       z. The macrogeometry therefore contains
       :math:`(2n_x+1)(2n_y+1)(2n_z+1)` copies in total. The default
       ``[0 0 0]`` means no macrogeometry.
   * - ``shiftVec``
     - ``shiftVec``
     - Distance in metres between neighbouring copies along x, y and z.

A one-dimensional array of period :math:`A` along x, with 50 copies on either
side, is therefore

.. code-block:: python

    n_macro = np.zeros(3)
    n_macro[0] = 50
    shift_vec = np.zeros(3)
    shift_vec[0] = A

    problem = MicromagProblem(res=res, n_macro=n_macro, shiftVec=shift_vec, ...)

.. code-block:: matlab

    problem.n_macro  = int32([50 0 0]);
    problem.shiftVec = [A 0 0];

Note that ``shiftVec`` is the *spacing between copies*, not the size of the
domain: setting it larger than the domain leaves a gap between the copies,
which is exactly how a periodic array of separated particles is modelled.

For a gapless periodic medium, set ``shiftVec`` equal to ``grid_L`` along
every periodic direction, i.e. :math:`n_x \Delta x` for :math:`n_x` cells of
size :math:`\Delta x`. Do not use the distance between the centres of the two
end cells, :math:`(n_x-1)\Delta x`. The copies would then overlap by one layer
of cells, and for a small domain the field would be wrong by up to 100 %.
MagTense stops with an error when ``shiftVec`` is smaller than ``grid_L`` along
a direction with copies. The ``macrogeometry_supercell_test`` example checks the
gapless case: a periodic domain has to give exactly the field of an explicitly
replicated domain.

.. note::
   Some analytical benchmarks for the periodic field are derived for the
   demagnetization tensor evaluated at the cell centres rather than averaged
   over the cell. Those cases need ``useAvgN`` switched off, see
   :ref:`Cell-averaged versus point-evaluated tensor`.

----------------------------------------
Sample shape correction
----------------------------------------

The macrogeometry represents wires and thin films directly. For a
three-dimensional sample the shape correction goes one step further and adds
the difference between the demagnetization field of the **real sample** and
that of the **macrogeometry**, both approximated by rectangular prisms centred
on the mesh:

.. math::

    \mathbf{H}_\mathrm{shape} = \left( \mathrm{N}_\mathrm{sample}
        - \mathrm{N}_\mathrm{macro} \right) \langle \mathbf{M} \rangle.

Here :math:`\langle \mathbf{M} \rangle` is the average magnetization of the
simulated cells multiplied by the volume fraction of the domain that is
actually occupied by magnetic material. That fraction is computed from the cell
volumes and from ``macroShape``, so voids and non-magnetic regions are
accounted for. :math:`\mathrm{N}_\mathrm{sample}` is the 
demagnetisation tensor for the sample and :math:`\mathrm{N}_\mathrm{macro}` for 
the macrogeometry. Note that the absolute size of the sample matters for the spatial
distribution of :math:`\mathrm{N}_\mathrm{sample}`. For a sample orders of magnitude
larger than the simulation domain, :math:`\mathrm{N}_\mathrm{sample}` is spatially uniform.

.. list-table::
   :widths: 22 22 56
   :header-rows: 1

   * - Matlab
     - Python
     - Description
   * - ``macroShape``
     - ``macroShape``
     - Side lengths of the prism that approximates the whole macrogeometry, in
       metres. Default ``[1 1 1]``.
   * - ``sampleShape``
     - ``sampleShape``
     - Side lengths of the prism that approximates the physical sample, in
       metres. Default ``[1 1 1]``.

Leaving both at their defaults makes the two tensors identical, so the
correction vanishes; a zero side length in ``macroShape`` switches it off
explicitly.

The evaluation points are shifted into a frame centred on the bounding box of
the mesh before the two tensors are evaluated, so an externally supplied mesh
that runs from 0 to L is handled the same way as the internally generated
uniform grid, which is centred on the origin. For an unstructured prism mesh
the bounding box is spanned by the outer faces of the cells, taken from
``grid_abc``. This holds even when the mesh is refined at one end, so the end
cells differ in size.

.. note::
   The correction is applied on both the dense and the FMM demagnetization
   paths, for the uniform and the unstructured prism grids, where the occupied
   volume fraction follows from the cell sizes.

----------------------------------------
Periodic exchange boundaries
----------------------------------------

``exchPBC`` is a three-element integer array of 0/1 flags that makes the
exchange coupling periodic along x, y and z. It is independent of the
macrogeometry: a domain can be exchange-periodic without a macrogeometry and
vice versa.

.. code-block:: python

    problem = MicromagProblem(
        res=res, exchPBC=np.array([1, 1, 1], dtype=np.int32), ...
    )

.. code-block:: matlab

    problem.exchPBC = int32([1 1 1]);

How it is realised depends on the grid:

* **Uniform grid** - the finite-difference stencil wraps around, so the cell at
  one end of a periodic direction gets the cell at the other end as its
  neighbour instead of a free boundary.
* **Unstructured prisms** - the mesh analysis links the cells at the two ends of
  each periodic direction, i.e. they are treated as sharing a face and enter
  each other's interpolation stencils. The exchange operator is then assembled
  exactly as it is without periodic boundaries. The period is taken as the
  extent of the mesh along that direction.
* **Tetrahedral mesh** - the mesh analysis links the two boundary planes of each
  periodic direction, by one of two mechanisms chosen per direction depending on
  whether the mesh itself is periodic. Both are described in
  :ref:`Periodic boundaries on a tetrahedral mesh` below. The period is taken as
  the extent of the mesh along that direction, as for the prisms.
* An exchange matrix **passed in from outside** carries its own boundary
  conditions, periodic or otherwise.

The mesh has to be thick enough along a periodic direction for the wrap-around
neighbour of a boundary cell to be distinct from its ordinary neighbour: at
least three cells on a uniform grid, and on an unstructured mesh an extent
larger than twice the largest cell size along that direction. The solver
checks both. A single cell along a direction is simply dropped from the
operator, which is the correct periodic answer for one cell.

The ``periodic_exchange_test`` example verifies that exchange-coupled moments
across a periodic boundary end up identical, on a uniform grid, on an
unstructured mesh of prisms and on a tetrahedral mesh.

------------------------------------------
Periodic boundaries on a tetrahedral mesh
------------------------------------------

A prism mesh is analysed geometrically, so linking a periodic boundary is a
matter of finding which faces at one end of the domain overlap which faces at
the other. A tetrahedral mesh is analysed from its connectivity instead - two
tetrahedra are neighbours when they share three nodes - and that changes what
the periodic boundary needs. ``TetrahedralMeshAnalysis.f90`` therefore has two
mechanisms, and picks between them separately for each periodic direction.

Periodic mesh: node identification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A mesh is *periodic along a direction* when the surface triangulation on the two
boundary planes normal to it is a translated copy of itself: every node on one
plane has a partner on the other, and the triangles match one for one. This is
what a mesh generator produces when it is asked for a periodic mesh - in gmsh
through ``Periodic Surface``, in COMSOL through a Copy Face mesh operation.

When that holds, each pair of boundary nodes is merged into a single node and
the whole analysis runs on the resulting connectivity. Two tetrahedra on
opposite sides of the domain then share three nodes and are found to be
neighbours in exactly the same way as any interior pair: they share one face,
with one area and one normal, and they enter each other's interpolation
stencils. Nothing downstream distinguishes them from an interior pair, apart
from measuring their separation through the boundary rather than across the
domain.

This mechanism involves no geometry beyond pairing the nodes, and it is exact.
When it is used, the analysis reports::

    Periodic exchange along x: mesh is periodic, 16 node pairs merged

Nodes on an edge or a corner of the domain have more than one partner - three
with two periodic directions, seven with three - and the merging is transitive,
so they collapse to a single node rather than to a set of pairs.

Non-periodic mesh: mortar coupling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A mesh generated without asking for periodicity, such as one from Matlab's
``generateMesh`` used by ``CreateTetraMesh``, triangulates the two boundary
planes independently, so the triangles do not match one for one. MagTense
handles such a mesh as well: the boundary faces are **mortar coupled**. Every boundary triangle
on one plane is clipped against the triangles it overlaps on the other, and each
intersection polygon becomes a *sub-face* shared by the two elements on either
side, carrying its own area and centroid. A triangle that straddles three
triangles on the opposite plane therefore contributes three sub-faces rather
than one face.

This is well behaved because a periodic boundary is a special case of a
non-conforming interface: once the period is subtracted the two planes are
*exactly coplanar*, so the clipping is a plain two-dimensional intersection of
two triangles in the plane of the boundary, with no projection error. The
sub-faces therefore tile both parent triangles exactly - no gaps and no overlaps
- and the flux leaving an element through its boundary face is the sum over its
sub-faces and equals what enters on the other side. **Flux conservation is
exact.**

The interpolation stencil is where the two mechanisms differ. A sub-face has no
nodes of its own, so its stencil is taken as the union of the stencils of its
two parent triangles, a superset of what a matching face would give, which
needs no geometric search or tolerance. The operator stays conservative, with
a somewhat larger local truncation error on the coupled plane than in the
interior of the mesh. A periodic mesh gives full accuracy on the plane, and
the difference can be measured directly, for instance by comparing the energy
of a domain wall on the plane with that of the same wall in the interior.

The analysis reports which mechanism it used along each direction. For mortar
coupling it names the direction, suggests meshing periodically for full
accuracy, and reports the number of sub-faces together with the check that
the tiling is exact::

    Mortar coupling: 72 boundary faces split into 54 sub-faces
    Mortar coupling: worst relative area mismatch   2.84E-15

The last line compares the summed sub-face area of every parent triangle with
the area of that triangle and reports the worst relative mismatch, a rounding
level number for a mesh that can be coupled.

Mixing the two
~~~~~~~~~~~~~~

The choice is made per direction, so a mesh that is periodic along one direction
and not along another uses node identification where it can and mortar coupling
where it must. A single run can report::

    Periodic boundary conditions were requested along x, but the mesh
    does not match across it: the two planes hold 16 and 34 nodes,
    and only 16 of them could be paired up.
    Falling back to MORTAR COUPLING of the boundary faces.
    ...
    Periodic exchange along y: mesh is periodic, 20 node pairs merged
    Periodic exchange along z: mesh is periodic, 20 node pairs merged
    Mortar coupling: 72 boundary faces split into 54 sub-faces

The node pairing is applied when *every* node on the plane has a partner, so
each direction is handled by one mechanism throughout.

What the geometry has to satisfy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mortar coupling matches up two differently meshed planes; the *geometry*
itself has to be periodic, i.e. the two boundary planes cover the same cross
section. The analysis verifies this through the area check above: where part
of a boundary face has nothing on the opposite plane to couple to, the summed
sub-face area falls short of the parent area, and the analysis stops and says
that the geometry, not just the mesh, has to be periodic. Together with the
thickness requirement listed above, this is what a tetrahedral mesh has to
satisfy to be made periodic.
