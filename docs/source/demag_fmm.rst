Demag field - FMM
=================

Calculating the demagnetization (stray) field is typically the most
computationally demanding part of a micromagnetic simulation. MagTense's
default approach - the fully analytical demagnetization tensor - is exact, with
memory and work that scale as :math:`O(N^2)` in the number of cells.

For such systems MagTense can evaluate the demagnetization field with the
**Fast Multipole Method (FMM)** through the dip-fmm library, which reduces the
cost of a field evaluation to :math:`O(N)`.

.. note::
   The FMM path is **off by default**. It has to be enabled explicitly with
   ``use_cdfmm`` (``use_fmm`` in MATLAB), and it is only available when the
   library was built with ``USE_CDFMM=1``. The earlier FMM3D-based backend
   (``USE_FMM3D``, ``fmm_eps``, ``fmm_cells``, ...) has been removed; in MATLAB
   the names ``use_fmm`` and ``USE_FMM`` now refer to dip-fmm.

Compilation
-----------

dip-fmm is built from a ``dip-fmm`` checkout inside the MagTense repository or
next to it, and installed under ``dip-fmm/local`` by the root Makefile:

.. code-block:: bash

    make python USE_CUDA=1 USE_CVODE=1 USE_MATLAB=0 USE_CDFMM=1

The dip-fmm CUDA backend is enabled only with ``USE_CUDA=1``. A checkout
somewhere else is selected with ``CDFMM_DIR=<path>``.

On Linux ``dip-fmm/local/lib`` must be on ``LD_LIBRARY_PATH`` when running.

On Windows the same command is run from the shell used for the rest of the
build (conda environment activated, then ``vcvars64.bat``). dip-fmm is compiled
with the environment's ``icx`` and, for CUDA, its ``nvcc`` with MSVC as host
compiler, and ``cdfmm_c.dll`` is copied next to the Python extension, so nothing
has to be added to ``PATH``. This needs a dip-fmm that carries its Windows
port (file mapping in the cache layer, the link order of the UCRT and Intel's
``libmmd``). There ``long double`` is the 64-bit ``double``, so the exact
prism and tetrahedron pair tensors are evaluated in double rather than
extended precision. ``python\build_python_fmm.cmd`` sets that shell up
(through ``tools\windows\env.cmd`` of the dip-fmm checkout) and runs the
command; ``matlab\build_mex_fmm.cmd`` does the same for the MEX files
described next.

For MATLAB the libraries are built with ``USE_MATLAB=1 USE_CDFMM=1`` and the
MEX files with ``buildMagTenseMEX('USE_FMM', true, 'fmm_root',
'<dip-fmm>/local', ...)``; ``fmm_root`` defaults to a dip-fmm checkout next to
the MagTense repository. With ``USE_MATLAB=1`` dip-fmm is linked with oneMKL
inside the library (``CDFMM_MKL_LINK=static``) and, on Linux, with its C++
runtime as well: MATLAB brings its own MKL and ``libstdc++`` into the process,
and a dip-fmm that shares the environment's would bind to those instead. Add
``CDFMM_ENABLE_CUDA=ON`` when building the libraries of the NoCUDA solver, so
that one dip-fmm library serves both solver MEX files. On Windows
``buildMagTenseMEX`` copies ``cdfmm_c.dll`` into ``MEX_files``; on Linux the
MEX files get the run path of ``<dip-fmm>/local/lib``.

In the Visual Studio solution every configuration without ``_no_FMM`` in its
name builds ``MagTenseMicroMag`` with dip-fmm (``USE_CDFMM`` is defined and
``DipFmmDemag.f90`` and dip-fmm's ``cdfmm_fortran.f90`` are compiled); the
``_no_FMM`` configurations leave it out. Visual Studio does not build dip-fmm
itself: it expects an installed dip-fmm in ``dip-fmm/local`` of a checkout next
to the MagTense repository, built for MATLAB as described above
(``tools\windows\build.cmd static`` in the dip-fmm repository). The MEX files
are then built with ``buildMagTenseMEX('VS_STUDIO', true, 'USE_FMM', true,
...)``, which links the ``_no_FMM`` libraries when ``USE_FMM`` is false.

How it is used
--------------

* The cells are treated as point dipoles with moment
  :math:`\mathbf{m} = \Delta V \, M_s \, \hat{\mathbf{m}}`, and the field is
  :math:`\mathbf{H} = -\nabla \phi`.
* The FMM plan (octree and translation operators) is built once per solve for
  the fixed cell geometry and reused for every field evaluation. The octree is
  either uniform (fixed depth ``cdfmm_depth``) or adaptive (split by
  ``cdfmm_capacity`` up to ``cdfmm_max_depth``); by default its root is aligned
  so that every leaf holds whole cells.
* The macrogeometry/sample shape correction and the demag-field noise
  (``CV``) are applied exactly as on the dense path.
* CUDA is used when MagTense is built with CUDA and the problem is run with
  ``cuda=True``; otherwise dip-fmm runs on the CPU (oneMKL, or a portable
  fallback). If the dip-fmm plan cannot be created, or an evaluation fails,
  the solve continues with the regular dense demagnetization tensor.

.. warning::
   dip-fmm supports **uniform grids** (``grid_type='uniform'``) and grids of
   axis-aligned prisms (``grid_type='unstructuredPrisms'``, e.g. octree-refined
   meshes) only. The field is computed in single precision (FP32), matching
   MagTense's demagnetization arrays. It is available from Python and from
   MATLAB.

FMM input variables
-------------------

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   * - Python
     - Default
     - Description
   * - ``use_cdfmm``
     - ``False``
     - Use dip-fmm instead of the dense demagnetization tensor.
   * - ``cdfmm_order``
     - 6
     - Expansion order.
   * - ``cdfmm_depth``
     - 4
     - Depth of the uniform tree.
   * - ``cdfmm_basis``
     - ``"spherical"``
     - Expansion basis, ``"spherical"`` or ``"cartesian"``.
   * - ``cdfmm_tree``
     - ``"uniform"``
     - ``"uniform"``, ``"uniform_strict"`` (every occupied box split to
       ``cdfmm_depth``) or ``"adaptive"``.
   * - ``cdfmm_capacity``
     - 32
     - Adaptive tree: maximum number of tiles per leaf before a split.
   * - ``cdfmm_max_depth``
     - 5
     - Adaptive tree: deepest level, 0 to 8.
   * - ``cdfmm_root``
     - ``None``
     - Root cube as (centre x, y, z, half-width) [m]; ``None`` lets MagTense
       choose it.
   * - ``cdfmm_align``
     - ``True``
     - With ``cdfmm_root=None``, choose a root on which every leaf holds whole
       cells, so no cell sticks out of its leaf box.

In MATLAB the same settings are properties of ``DefaultMicroMagProblem``,
named ``fmm_`` instead of ``cdfmm_`` and with integer codes instead of names:
``use_fmm`` (set with ``setUseDipFmm``), ``fmm_order``, ``fmm_depth``,
``fmm_basis`` (0 spherical, 1 Cartesian), ``fmm_tree`` (0 uniform, 1 adaptive,
2 adaptive with ``fmm_cap`` as the occupancy to aim for), ``fmm_cap``
(``cdfmm_capacity``), ``fmm_maxdep`` (``cdfmm_max_depth``), ``fmm_root``
(``[0, 0, 0, -1]`` for none) and ``fmm_align`` (0 or 1). ``"uniform_strict"``
has no MATLAB counterpart.

FMM Python example
------------------

.. code-block:: python

    problem = MicromagProblem(
        res=[64, 64, 64],
        grid_L=[500e-9, 500e-9, 500e-9],
        cuda=True,
        use_cdfmm=True,
        cdfmm_order=6,
        cdfmm_depth=4,
        cdfmm_basis="spherical",
    )

    result = problem.run_simulation(
        t_end=t_end, nt=nt, fct_h_ext=h_ext_fct, nt_h_ext=nt_h_ext
    )

See ``python/examples/micromagnetism/FMM`` for a comparison of dip-fmm against
the dense demagnetization tensor.

FMM MATLAB example
------------------

.. code-block:: matlab

    problem = DefaultMicroMagProblem(64, 64, 64);
    problem.grid_L = [500e-9, 500e-9, 500e-9];
    problem = problem.setUseCuda(true);
    problem = problem.setUseDipFmm(true);
    problem.fmm_order = int32(6);
    problem.fmm_depth = int32(4);

    solution = problem.MagTenseLandauLifshitzSolver_mex(struct(problem), struct());

``matlab/examples/Micromagnetism/FMM/fmm_vs_regular.m`` is the MATLAB
counterpart of the Python comparison.

Tests
-----

``python/examples/micromagnetism/FMM/fmm_smoke_test.py`` is a one-minute
version of the comparison, meant as the check after a build. The demagnetizing
field of dip-fmm at ``t = 0`` is tested against the dense tensor by
``dipfmm_octree_test`` (a graded mesh of prisms, the uniform and the adaptive
tree, the aligned root) and ``dipfmm_tetra_test`` (a tetrahedral mesh), in
``python/examples/micromagnetism/MagTense_tests`` and, with the same problems
and limits, in ``matlab/examples/Micromagnetism/MagTense_tests``. They need a
build with dip-fmm and are therefore not part of ``testMagTenseFunctions``.
