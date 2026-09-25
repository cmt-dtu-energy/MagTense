Demag field - FMM
=================

Calculating the demagnetization (stray) field is typically the most
computationally demanding part of a micromagnetic simulation. MagTense's
default approach - the fully analytical demagnetization tensor - is exact, but
its :math:`O(N^2)` memory and work scaling becomes prohibitive for very large
systems.

To address this, MagTense includes an implementation of the **Fast Multipole
Method (FMM)**, which reduces the computational complexity to :math:`O(N)`.

.. note::
   The FMM path is **off by default**. It has to be enabled explicitly with
   ``use_fmm``, and it is ignored altogether unless the library was built with
   ``USE_FMM3D=1``. Opting in explicitly means that a problem does not silently
   change its demagnetization path when the library is rebuilt with FMM
   support.

Implementation and attribution
------------------------------

The FMM acceleration in MagTense is built upon the **FMM3D** library developed
by the **Flatiron Institute**.

* **Core library:** `Flatiron Institute FMM3D <https://github.com/flatironinstitute/fmm3d>`_
* **Linking:** MagTense links against a fork, `Ximtecs/FMM3D
  <https://github.com/Ximtecs/FMM3D>`_, which contains updated build
  configurations and makefiles to allow integration with the MagTense Fortran
  core.

Compilation
-----------

To enable FMM support during the build, include the FMM flag in the make
command. This links the FMM3D libraries and enables the specialised Fortran
modules:

.. code-block:: bash

    make USE_FMM3D=1

On Linux, ``$(MagTense)/external/FMM3D/local`` must be on ``LD_LIBRARY_PATH``
when building. Without ``USE_FMM3D=1`` the FMM source is not compiled at all
and the ``use_fmm`` flag has no effect.

Optimization: persistent tree structure
---------------------------------------

In micromagnetic simulations the spatial distribution of cells is static
throughout the simulation. To exploit this, MagTense uses a *magtense-local*
tree structure that caches the octree setup between consecutive calls to the
solver. By reusing the tree, the cost of re-partitioning space at every time
step is eliminated.

.. note::
   The current implementation requires a fully grown tree, so **ifunif** must
   always be set to **1**.

Near-field evaluation and neighbour tensors
-------------------------------------------

In standard FMM implementations, near-field interactions are handled by a
direct point-to-point (P2P) evaluation. In MagTense this standard P2P
evaluation is **disabled**, and the high-precision analytical demagnetization
tensors are used for all near-field interactions instead:

1. **Neighbour identification:** based on *List 1* from the FMM tree creation,
   MagTense identifies neighbour pairs, i.e. cells within the same or adjacent
   leaf nodes.
2. **Sparse neighbour tensor:** a sparse tensor structure is created to map
   these neighbour interactions.
3. **Analytical calculation:** the exact analytical demagnetization tensor is
   evaluated for every neighbour pair.
4. **Sparse matrix storage:** the values are converted into a sparse matrix. On
   **CPU** this is an Intel MKL sparse matrix; on **GPU** it is stored
   persistently in global memory for **CUDA** evaluation at each time step.

The macrogeometry and sample :ref:`Sample shape correction` are applied on the
FMM path as well, so switching to FMM does not change those terms.

FMM input variables
-------------------

.. list-table::
   :widths: 24 22 8 46
   :header-rows: 1

   * - Python
     - Matlab
     - Type
     - Description
   * - ``use_fmm``
     - ``use_fmm``
     - int
     - **1** to use the FMM demagnetization path, **0** for the dense
       analytical tensor. Default **0**.
   * - ``fmm_eps``
     - ``fmm_eps``
     - float
     - Requested accuracy, which controls the exponential order of the
       plane-wave expansion. Default ``1e-4``.
   * - ``fmm_nterms``
     - ``fmm_nterms``
     - int
     - Multipole expansion order. A negative value (default ``-1``) derives the
       order from ``fmm_eps`` during the tree build.
   * - ``ifunif``
     - ``ifunif``
     - int
     - Tree type. Must be **1** (uniform tree).
   * - ``nlmin``
     - ``nlmin``
     - int
     - Minimum level of the octree hierarchy. Default 1.
   * - ``nlmax``
     - ``nlmax``
     - int
     - Maximum level of the octree hierarchy. For the required uniform tree,
       this parameter controls the tree depth. Default 5.
   * - ``allow_fmm_short_circuit``
     - ``fmm_short``
     - int
     - **1** to fall back to the direct calculation for small problems, **0**
       to force FMM. Default 1.
   * - ``fmm_min_n``
     - ``fmm_min_n``
     - int
     - Cell count below which FMM is disabled when the short circuit is
       allowed. Default 20000.

When the short circuit triggers, the solver reports
``MagTense: problem smaller than fmm_min_n - disabling FMM and using the full
demag tensor`` and continues with the dense tensor. The test is made before the
octree is built, so a small or effectively one-dimensional geometry never
reaches the tree construction.

.. warning::
   The FMM implementation is currently configured for rectangular-prism cells,
   corresponding to the MagTense grid types ``gridTypeUniform`` and
   ``gridTypeUnstructuredPrisms``. Cell dimensions may be supplied either
   through the uniform-grid dimensions ``dx``, ``dy``, and ``dz``, or through
   the per-cell ``grid_abc`` dimensions for unstructured prism grids.

FMM Python example
------------------

.. code-block:: python

    problem.use_fmm = 1
    problem.fmm_eps = 1e-4
    problem.fmm_nterms = 12
    problem.ifunif = 1
    problem.nlmin = 1
    problem.nlmax = 2
    problem.allow_fmm_short_circuit = 1
    problem.fmm_min_n = 20000

    result = problem.run_simulation(
        t_end=t_end, nt=nt, fct_h_ext=h_ext_fct, nt_h_ext=nt_h_ext
    )

FMM Matlab example
------------------

.. code-block:: matlab

    problem.use_fmm    = int32(1);
    problem.fmm_eps    = 1e-4;
    problem.fmm_nterms = int32(12);
    problem.ifunif     = int32(1);
    problem.nlmin      = int32(1);
    problem.nlmax      = int32(2);
    problem.fmm_short  = int32(1);
    problem.fmm_min_n  = int32(20000);
