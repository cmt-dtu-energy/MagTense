Demag field - FMM
=================

Calculating the demagnetization (stray) field is typically the most
computationally demanding part of a micromagnetic simulation. MagTense's
default approach - the fully analytical demagnetization tensor - is exact, but
its :math:`O(N^2)` memory and work scaling becomes prohibitive for very large
systems.

For such systems MagTense can evaluate the demagnetization field with the
**Fast Multipole Method (FMM)** through the dip-fmm library, which reduces the
cost of a field evaluation to :math:`O(N)`.

.. note::
   The FMM path is **off by default**. It has to be enabled explicitly with
   ``use_cdfmm``, and it is only available when the library was built with
   ``USE_CDFMM=1``. The earlier FMM3D-based backend (``USE_FMM3D``,
   ``use_fmm``, ``fmm_eps``, ...) has been removed.

Compilation
-----------

dip-fmm is built from the ``dip-fmm`` checkout next to the MagTense sources and
installed under ``dip-fmm/local`` by the root Makefile:

.. code-block:: bash

    make python USE_CUDA=1 USE_CVODE=1 USE_MATLAB=0 USE_CDFMM=1

The automatic dip-fmm build is currently implemented for Linux only, and its
CUDA backend is enabled only with ``USE_CUDA=1``. ``dip-fmm/local/lib`` must be
on ``LD_LIBRARY_PATH`` when running.

How it is used
--------------

* The cells are treated as point dipoles with moment
  :math:`\mathbf{m} = \Delta V \, M_s \, \hat{\mathbf{m}}`, and the field is
  :math:`\mathbf{H} = -\nabla \phi`.
* The FMM plan (octree and translation operators) is built once per solve for
  the fixed uniform cuboid grid and reused for every field evaluation.
* The macrogeometry/sample shape correction and the demag-field noise
  (``CV``) are applied exactly as on the dense path.
* CUDA is used when MagTense is built with CUDA and the problem is run with
  ``cuda=True``; otherwise dip-fmm runs on the CPU (oneMKL, or a portable
  fallback). If the dip-fmm plan cannot be created, or an evaluation fails,
  the solve continues with the regular dense demagnetization tensor.

.. warning::
   dip-fmm currently supports **uniform grids** (``grid_type='uniform'``) only,
   and the field is computed in single precision (FP32), matching MagTense's
   demagnetization arrays. It is available from Python; the MATLAB interface
   does not expose it.

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
     - Fixed octree depth. Aligning the octree with the grid (e.g. a power of
       two number of cells per side) avoids cells straddling box boundaries.
   * - ``cdfmm_basis``
     - ``"spherical"``
     - Expansion basis, ``"spherical"`` or ``"cartesian"``.

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
