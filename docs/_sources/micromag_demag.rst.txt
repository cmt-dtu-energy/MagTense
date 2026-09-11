Demagnetization field
========================================

The demagnetization (stray) field is the most computationally expensive part of micromagnetics.
MagTense computes it from the **fully analytical demagnetization tensor** of
the cell shape; a rectangular prism for the uniform and unstructured-prism
grids, and a tetrahedron for the tetrahedral grid. The accuracy of these tensors is documented in
`Bjørk, R. and d'Aquino, M.: Accuracy of the analytical demagnetization tensor
for various geometries, Journal of Magnetism and Magnetic Materials, 587,
171245, 2023 <https://www.sciencedirect.com/science/article/pii/S0304885323008958>`_.

The tensor is symmetric, so it is stored as the six components
:math:`K_{xx}, K_{xy}, K_{xz}, K_{yy}, K_{yz}, K_{zz}`, each a dense
:math:`n_\mathrm{tot} \times n_\mathrm{tot}` matrix in **single precision**.
The memory requirement is therefore approximately
:math:`24\,n_\mathrm{tot}^2` bytes, which is what limits the problem size; the
Matlab interface prints an estimate when the problem struct is created. The
saturation magnetization is folded into the tensor when it is built, so a
spatially varying :math:`M_s` costs nothing extra.

.. note::
   The parameter ``usePres`` (Python ``precision``) is accepted but has no
   effect - the demagnetization tensor is always kept in single precision,
   while everything else is double precision.

----------------------------------------
Switching the demagnetization off
----------------------------------------

Setting ``useDemag`` (Matlab ``setUseDemag``, Python ``usedemag``) to False removes the
demagnetization field from the effective field and skips building the tensor
altogether.

-------------------------------------------
Cell-averaged versus point-evaluated tensor
-------------------------------------------

``useAvgN`` (Python ``useavgn``, default **on**) selects the *volume-averaged*
prism tensor, i.e. the tensor averaged over the receiving cell rather than
evaluated at its centre. This is the more accurate choice for a uniform mesh
and is the default, but it is only implemented for the rectangular prism tile
(tile_type = 8).

Some analytical benchmarks are formulated for the tensor evaluated at
the cell centre, and for those ``useAvgN`` must be switched off. This includes the
macrogeometry and shape-correction tests.

For the ``unstructuredPrisms`` grid the tensor may additionally be averaged
numerically over the receiving cell using ``N_ave = [nx, ny, nz]``, which
evaluates the tensor on an ``nx*ny*nz`` sub-grid inside each receiving prism
and averages. The default ``[1 1 1]`` means no averaging. This is not supported
for tetrahedral meshes.

----------------------------------------
Approximating the tensor
----------------------------------------

The dense tensor may be sparsified, which trades accuracy for memory and speed.
The approximation is selected with ``setMicroMagDemagApproximation`` in Matlab
and with ``demag_approx`` in Python, and the associated cut-off is
``dem_thres``.

.. list-table::
   :widths: 28 8 64
   :header-rows: 1

   * - Name
     - Value
     - Description
   * - ``none`` (Python ``None``)
     - 1
     - Default. The dense tensor is used as is.
   * - ``threshold``
     - 2
     - Every element with :math:`|K| <` ``dem_thres`` is set to zero and the
       six matrices are converted to MKL sparse matrices.
   * - ``threshold_fraction``
     - 4
     - As above, but ``dem_thres`` is interpreted as the *fraction* of tensor
       elements to discard. The cut-off value is found by bisection. A value
       :math:`\geq 1` zeroes the demagnetization field entirely.
   * - ``fft_thres``
     - 3
     - Threshold applied in Fourier space. **Stale, do not use.**
   * - ``fft_threshold_fraction``
     - 5
     - Fractional threshold applied in Fourier space. **Stale, do not use.**

.. warning::
   The two Fourier-space approximations are not maintained. The solver prints a
   loud warning if they are selected: the transformed field is never mapped
   back onto the demagnetization field and the tensor is not corrected for a
   spatially varying :math:`M_s`, so the result is wrong. Use ``none``,
   ``threshold`` or ``threshold_fraction``.

-----------------------------------------
Adding noise to the demagnetization field
-----------------------------------------

``CV`` (Python ``cv``) adds a normally distributed relative error to the
demagnetization field, with ``CV`` being the coefficient of variation, i.e. the
ratio of the standard deviation to the mean. The default of zero disables it.
The draw is controlled by the same random seed as the thermal field, see
:ref:`Reproducible random numbers`.

----------------------------------------
Storing and reusing the tensor
----------------------------------------

Because building the tensor can dominate the run time of a short simulation,
it can be written to disk and read back:

.. list-table::
   :widths: 22 22 56
   :header-rows: 1

   * - Matlab
     - Python
     - Description
   * - ``N_ret``, ``N_file_out``
     - ``N_ret``, ``N_file_out``
     - ``1`` does not return the tensor, ``2`` returns it in memory, and a
       value ``> 2`` writes it to the file named by ``N_file_out``. The value
       is the length of that filename.
   * - ``N_load``, ``N_file_in``
     - ``N_load``, ``N_file_in``
     - The same encoding for loading a previously stored tensor instead of
       computing it.

The Matlab helpers ``setReturnNFilename`` and ``setLoadNFilename`` set the
filename and the matching length in one call. In Python the ``filename``
constructor argument sets all four fields; its default ``"t"`` has length one
and therefore means "do not return and do not load".

----------------------------------------
Parallelism and hardware
----------------------------------------

* ``nThreads`` (Python ``n_threads``) is the number of OpenMP threads used when
  building the demagnetization tensor.
* ``useCuda`` (Matlab ``setUseCuda``, Python ``cuda``) evaluates the
  tensor-vector product on an NVIDIA GPU at every time step. The Python
  interface checks for ``nvidia-smi`` and falls back to the CPU with a warning
  if no GPU is present. In Matlab, ``setUseCuda`` also picks the matching
  MEX-file.
* The **Fast Multipole Method** replaces the dense tensor by an :math:`O(N)`
  evaluation and is described in :ref:`Demag field - FMM`. It is off by
  default and requires a build with ``USE_FMM3D=1``.

.. note::
   The parameter ``demigstp`` (``demag_ignore_steps``), intended to recompute
   the demagnetization tensor only every n'th step of a hysteresis
   calculation, is accepted by both interfaces but is not acted on by the
   current solver.
