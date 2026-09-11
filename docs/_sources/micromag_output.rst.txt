Micromagnetic output
========================================

----------------------------------------
Matlab solution struct
----------------------------------------

The MEX-file returns the solution struct, and optionally the grid information
as a second output:

.. code-block:: matlab

    solution = struct();
    [solution, GridInfo] = problem.MagTenseLandauLifshitzSolver_mex( struct(problem), solution );

The solution struct contains

.. list-table::
   :widths: 16 30 54
   :header-rows: 1

   * - Field
     - Size
     - Description
   * - ``t``
     - ``(nt)``
     - The times at which the solution was stored.
   * - ``M``
     - ``(nt, ntot, nt_Hext, 3)``
     - The reduced magnetization. The third dimension is the applied-field
       index, which is 1 for the dynamic solver and the number of field steps
       for the explicit solver.
   * - ``pts``
     - ``(ntot, 3)``
     - Centre coordinates of the cells, in the same order as the other per-cell
       arrays.
   * - ``H_exc``
     - ``(nt, ntot, nt_Hext, 3)``
     - Exchange field in A/m. Only filled when ``ReturnHall`` is set.
   * - ``H_ext``
     - ``(nt, ntot, nt_Hext, 3)``
     - Applied field in A/m. Only filled when ``ReturnHall`` is set.
   * - ``H_dem``
     - ``(nt, ntot, nt_Hext, 3)``
     - Demagnetization field in A/m. Only filled when ``ReturnHall`` is set.
   * - ``H_ani``
     - ``(nt, ntot, nt_Hext, 3)``
     - Anisotropy field in A/m. Only filled when ``ReturnHall`` is set.
   * - ``n_Hext_acc``
     - scalar
     - Number of accepted applied-field steps. Only meaningful for
       :ref:`Adaptive hysteresis`.

``ReturnHall`` (Python ``usereturnhall``) is off by default. Turning it on
costs four additional ``(nt, ntot, nt_Hext, 3)`` arrays, which for a large
problem is substantial, so leave it off unless the individual field terms are
actually needed.

.. note::
   The four ``H_*`` arrays are **not** the fields the integrator ran on. They
   are recomputed from the stored solution at the ``nt`` output times after the
   integration has finished, because the integrator's last right-hand-side
   evaluation is at an internal step rather than at the requested output state.

The optional ``GridInfo`` output carries the analysed mesh - face normals,
face areas, cell volumes, cell and face centres, the interpolation index
arrays, and the assembled exchange matrix in COO form (``ExchMat_r``,
``ExchMat_c``, ``ExchMat_v``, ``ExchMat_nr``, ``ExchMat_nc``). It is what
``cartesianUnstructuredMeshPlot`` plots.

The Matlab helper ``computeMagneticMomentGeneralMesh`` converts a solution into
the volume-averaged magnetic moment, and ``computeMagneticEnergy`` evaluates
the energy terms.

----------------------------------------
Python result list
----------------------------------------

``run_simulation`` and ``run_hysteresis`` return a list:

.. list-table::
   :widths: 10 22 68
   :header-rows: 1

   * - Index
     - Name
     - Description
   * - 0
     - ``t_out``
     - Output times.
   * - 1
     - ``M_out``
     - Reduced magnetization, ``(nt, ntot, nt_h_ext, 3)``.
   * - 2
     - ``pts``
     - Cell centres, ``(ntot, 3)``.
   * - 3
     - ``H_exc``
     - Exchange field, same shape as ``M_out``.
   * - 4
     - ``H_ext``
     - Applied field.
   * - 5
     - ``H_dem``
     - Demagnetization field.
   * - 6
     - ``H_ani``
     - Anisotropy field.
   * - 7
     - ``Exch_mat_ntot``
     - Number of non-zero entries in the exchange matrix.
   * - 8, 9, 10
     - ``Exch_mat_r``, ``Exch_mat_c``, ``Exch_mat_v``
     - The exchange matrix in COO form, already sliced to
       ``Exch_mat_ntot`` entries.
   * - 11, 12
     - ``Exch_mat_nr``, ``Exch_mat_nc``
     - Number of rows and columns of the exchange matrix.

``run_hysteresis_adaptive`` returns the same twelve entries followed by one
extra element, the number of accepted field steps, and the arrays with an
applied-field dimension are already sliced to that number.

A typical unpacking is

.. code-block:: python

    t_out, M_out = problem.run_simulation(
        t_end=1e-9, nt=200, fct_h_ext=h_ext_fct, nt_h_ext=2000
    )[:2]

    # Volume-averaged reduced magnetization as a function of time
    M_avg = M_out[:, :, 0, :].mean(axis=1)

The helper ``magtense.utils.plot_M_thin_film`` plots the magnetization of a
thin-film problem, and ``magtense.utils.create_plot`` handles the magnetostatic
tile plots.

----------------------------------------
Dry runs
----------------------------------------

Setting ``dummy_run = 1`` builds the problem, the mesh, the exchange operator
and the demagnetization tensor, but skips the time integration. This is useful
for timing the setup phase, for checking the memory footprint of a large
problem, and for extracting the exchange matrix or the demagnetization tensor
without paying for a solve. The returned magnetization arrays are zeroed in
that case.
