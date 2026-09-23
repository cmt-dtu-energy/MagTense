Micromagnetic parameter reference
========================================

This page lists every micromagnetic input parameter, its Matlab name, its
Python name and its default value. The Matlab names are the field names of the
struct that is handed to Fortran, as defined by ``getProblemFieldnames`` in
`MagTenseMicroMagIO.f90 <https://github.com/cmt-dtu-energy/MagTense/blob/master/source/MagTenseMicroMag/MagTenseMicroMagIO.f90>`_.

In Matlab, parameters that map onto an internal integer flag are set through a
``set...`` method rather than by assigning the number; those methods are named
in the tables. Integer-valued fields must be assigned as ``int32``.

In Python, most parameters are constructor arguments of ``MicromagProblem``;
those that are not are marked "attribute" and are assigned after construction.

----------------------------------------
Grid and geometry
----------------------------------------

.. list-table::
   :widths: 18 18 16 48
   :header-rows: 1

   * - Matlab
     - Python
     - Default
     - Description
   * - ``grid_n``
     - ``res``
     - required
     - Number of cells along x, y and z. For unstructured meshes use
       ``[ntot,1,1]``.
   * - ``grid_L``
     - ``grid_L``
     - ``[500e-9,125e-9,3e-9]``
     - Size of the rectangular domain in m. Only used for a uniform grid.
   * - ``grid_type``
     - ``grid_type``
     - ``uniform``
     - ``uniform`` (1), ``tetrahedron`` (2) or ``unstructuredPrisms`` (3).
       Matlab: ``setMicroMagGridType``.
   * - ``grid_pts``
     - ``grid_pts``
     - ``0``
     - ``(ntot,3)`` cell centres for the two unstructured grid types.
   * - ``grid_abc``
     - ``grid_abc``
     - ``0``
     - ``(ntot,3)`` full side lengths of the prisms, ``unstructuredPrisms``
       only.
   * - ``grid_nod``
     - ``grid_nod``
     - ``0``
     - Node coordinates of a tetrahedral mesh.
   * - ``grid_ele``
     - ``grid_ele``
     - ``0``
     - Node indices of each tetrahedral element.
   * - ``grid_nnod``
     - ``grid_nnod``
     - ``0``
     - Number of nodes in the tetrahedral mesh.

See :ref:`Grids and meshes`.

----------------------------------------
Material and anisotropy
----------------------------------------

.. list-table::
   :widths: 18 18 16 48
   :header-rows: 1

   * - Matlab
     - Python
     - Default
     - Description
   * - ``Ms``
     - ``Ms``
     - ``8e5``
     - Saturation magnetization [A/m], per cell.
   * - ``A0``
     - ``A0``
     - ``1.3e-11``
     - Exchange stiffness [J/m], per cell.
   * - ``K0``
     - ``K0``
     - ``0``
     - Uniaxial anisotropy constant [J/m^3], per cell.
   * - ``K1``, ``K2``
     - ``K1``, ``K2``
     - ``0``
     - Cubic anisotropy constants [J/m^3], per cell.
   * - ``K0_arr``
     - ``K0_arr``
     - ``0``
     - ``(ntot,6,3)`` general anisotropy expansion [J/m^3].
   * - ``CrysAxis``
     - ``CrysAxis``
     - identity
     - ``(ntot,3,3)`` local crystal coordinate system.
   * - ``u_ea``
     - ``u_ea``
     - ``0``
     - ``(ntot,3)`` uniaxial easy axis. Python: attribute.
   * - ``alpha``
     - ``alpha``
     - ``4.42e3``
     - Damping constant [m/(A s)]. Zero switches to the ``alphat`` table.
   * - ``gamma``
     - ``gamma``
     - ``2.21e5``
     - Precession constant [m/(A s)]. Zero removes precession.
   * - ``alphat``, ``nt_alpha``
     - ``t_alpha``, ``alpha_fct``
     - ``0``
     - Tabulated time-dependent damping. Matlab: ``setAlpha``.
   * - ``MaxT0``
     - ``max_T0``
     - ``2``
     - Legacy, not used by the current solver.
   * - ``m0``
     - ``m0``
     - random
     - ``(ntot,3)`` initial reduced magnetization. Normalised on input.

See :ref:`Material parameters` and :ref:`Magnetocrystalline anisotropy`.

----------------------------------------
Solver, time and applied field
----------------------------------------

.. list-table::
   :widths: 18 18 16 48
   :header-rows: 1

   * - Matlab
     - Python
     - Default
     - Description
   * - ``solver``
     - ``solver``
     - ``Dynamic``
     - ``Explicit`` (1), ``Dynamic`` (2), ``Implicit`` (3, not implemented).
       Matlab: ``setMicroMagSolver``.
   * - ``ProblemMod``
     - ``prob_mode``
     - ``new``
     - ``new`` (1) or ``old`` (2). Matlab:
       ``setMicroMagProblemMode``.
   * - ``t``, ``nt``
     - ``t_end``, ``nt``
     - ``linspace(0,1,1000)``
     - Output times. Matlab: ``setTime``. In Python these are arguments of
       ``run_simulation``.
   * - ``Hext``, ``nt_Hext``
     - ``fct_h_ext``, ``nt_h_ext``
     - zero field
     - Applied field table of ``[t,Hx,Hy,Hz]`` rows [A/m]. Matlab:
       ``setHext``.
   * - ``tol``
     - ``tol``
     - ``1e-4``
     - Relative tolerance of the ODE solver.
   * - ``thres``
     - ``thres``
     - ``1e-6``
     - Threshold below which a solution component is treated as zero.
   * - ``t_conv``, ``nt_conv``
     - ``nt_conv``
     - ``0``, ``1``
     - Convergence-check times. Matlab:
       ``setConvergenceCheckTime``.
   * - ``conv_tol``
     - ``conv_tol``
     - ``1e-4``
     - Convergence criterion.
   * - ``setTimeDis``
     - ``setTimeDis``
     - ``10``
     - Progress is reported every n'th step.
   * - ``useCVODE``
     - ``cvode``
     - ``0``
     - Use CVODE instead of RKSuite. Matlab: ``setUseCVODE``.
   * - ``useCuda``
     - ``cuda``
     - ``0``
     - Evaluate the demagnetization field on a GPU. Matlab: ``setUseCuda``.

See :ref:`Solvers and time integration`.

----------------------------------------
Adaptive hysteresis parameters
----------------------------------------

.. list-table::
   :widths: 18 22 16 44
   :header-rows: 1

   * - Matlab
     - Python
     - Default
     - Description
   * - ``adaptiveHext``
     - ``hysteresis_solver``
     - ``0`` / ``static``
     - ``1`` / ``adaptive`` enables adaptive field stepping.
   * - ``maxHextSteps``
     - ``max_steps``
     - ``1000``
     - Maximum number of accepted field steps.
   * - ``H_start``, ``H_end``
     - ``H_start``, ``H_end``
     - ``[0 0 0]``
     - Start and end applied field [A/m].
   * - ``dH_initial``
     - ``dH_initial``
     - ``1e3``
     - First field step [A/m].
   * - ``dH_min``, ``dH_max``
     - ``dH_min``, ``dH_max``
     - ``1e1``, ``1e5``
     - Smallest and largest allowed field step [A/m].
   * - ``dH_grow``
     - ``dH_grow``
     - ``1.25`` / ``1.5``
     - Growth factor. The Matlab default is 1.25, the Python method default is
       1.5.
   * - ``dH_shrink``
     - ``dH_shrink``
     - ``0.5`` / ``0.75``
     - Shrink factor. The Matlab default is 0.5, the Python method default is
       0.75.
   * - ``dM_min``
     - ``dM_min``
     - ``1e-3``
     - Grow the step below this magnetization change.
   * - ``dM_target``
     - ``dM_target``
     - ``1e-2``
     - Shrink the step above this magnetization change.
   * - ``dM_reject``
     - ``dM_reject``
     - ``5e-2``
     - Reject and retry the step above this magnetization change.
   * - ``switch_refdH``, ``use_sw_ref``
     - ``switch_refine_dH``
     - ``0``, ``0``
     - Largest step accepted across a magnetization sign change [A/m].

See :ref:`Adaptive hysteresis`.

----------------------------------------
Exchange
----------------------------------------

.. list-table::
   :widths: 18 20 16 46
   :header-rows: 1

   * - Matlab
     - Python
     - Default
     - Description
   * - ``exch_meth``
     - ``exch_meth``
     - ``DirectLaplacianNeumann``
     - ``DirectLaplacianNeumann`` (1) or ``GGNeumann`` (2). Matlab:
       ``setMicroMagExchMethod``.
   * - ``exch_intpn``
     - ``exch_intpn``
     - ``Extended``
     - ``Extended`` (1) or ``Compact`` (2). Matlab:
       ``setMicroMagExchInterpn``.
   * - ``exch_weigh``
     - ``exch_weigh``
     - ``8``
     - Exponent of the inverse-distance face weighting.
   * - ``passExch``
     - ``passexch``
     - ``0``
     - Use an externally supplied exchange matrix. Matlab:
       ``setMicroMagpassExch`` or ``setExchangeMatrixCOO``.
   * - ``exch_nval``, ``exch_nrow``, ``exch_ncol``
     - ``exch_nval``, ``exch_nrow``, ``exch_ncols``
     - ``0`` / ``1``
     - Size of the supplied exchange matrix.
   * - ``exch_val``, ``exch_rows``, ``exch_cols``
     - ``exch_val``, ``exch_rows``, ``exch_cols``
     - ``0``
     - The supplied exchange matrix in COO form, 1-based indices.
   * - \-
     - ``exch_presize``
     - ``12``
     - Python only. Length of the returned exchange-matrix arrays as
       ``exch_presize * ntot``.
   * - ``exchPBC``
     - ``exchPBC``
     - ``[0 0 0]``
     - Periodic exchange boundaries along x, y and z.

See :ref:`Exchange interaction` and :ref:`Periodic exchange boundaries`.

----------------------------------------
Demagnetization
----------------------------------------

.. list-table::
   :widths: 18 20 16 46
   :header-rows: 1

   * - Matlab
     - Python
     - Default
     - Description
   * - ``useDemag``
     - ``usedemag``
     - ``1``
     - Include the demagnetization field. Matlab: ``setUseDemag``.
   * - ``useAvgN``
     - ``useavgn``
     - ``1``
     - Use the cell-averaged prism tensor rather than the tensor at the cell
       centre.
   * - ``N_ave``
     - ``N_ave``
     - ``[1 1 1]``
     - Numerical sub-grid averaging of the tensor over the receiving cell,
       ``unstructuredPrisms`` only.
   * - ``dem_appr``
     - ``demag_approx``
     - ``none``
     - ``none`` (1), ``threshold`` (2), ``fft_thres`` (3),
       ``threshold_fraction`` (4), ``fft_threshold_fraction`` (5). Matlab:
       ``setMicroMagDemagApproximation``.
   * - ``dem_thres``
     - ``dem_thres``
     - ``0``
     - Cut-off used by the threshold approximations.
   * - ``CV``
     - ``cv``
     - ``0``
     - Coefficient of variation of a random error added to the demagnetization
       field.
   * - ``N_ret``, ``N_file_out``
     - ``filename``
     - ``1``, ``'t'``
     - Return or write the demagnetization tensor. Matlab:
       ``setReturnNFilename``.
   * - ``N_load``, ``N_file_in``
     - ``filename``
     - ``1``, ``'t'``
     - Load a stored demagnetization tensor. Matlab: ``setLoadNFilename``.
   * - ``nThreads``
     - ``n_threads``
     - ``1``
     - OpenMP threads used when building the tensor.
   * - ``usePres``
     - ``precision``
     - ``0``
     - Accepted but not used; the tensor is always single precision.
   * - ``demigstp``
     - ``demigstp``
     - ``0``
     - Accepted but not used by the current solver.

See :ref:`Demagnetization field`.

----------------------------------------
Macrogeometry and shape correction
----------------------------------------

.. list-table::
   :widths: 18 18 16 48
   :header-rows: 1

   * - Matlab
     - Python
     - Default
     - Description
   * - ``n_macro``
     - ``n_macro``
     - ``[0 0 0]``
     - Number of domain copies on each side along x, y and z.
   * - ``shiftVec``
     - ``shiftVec``
     - ``[0 0 0]``
     - Spacing between neighbouring copies [m].
   * - ``macroShape``
     - ``macroShape``
     - ``[1 1 1]``
     - Side lengths of the prism approximating the macrogeometry [m].
   * - ``sampleShape``
     - ``sampleShape``
     - ``[1 1 1]``
     - Side lengths of the prism approximating the sample [m].

See :ref:`Periodic boundaries and macrogeometry`.

----------------------------------------
Thermal parameters
----------------------------------------

.. list-table::
   :widths: 18 18 16 48
   :header-rows: 1

   * - Matlab
     - Python
     - Default
     - Description
   * - ``temperature``
     - ``T``
     - ``0``
     - Per-cell temperature [K].
   * - ``rng_seed``
     - ``rng_seed``
     - ``0``
     - Random seed: ``0`` compiler default, ``>0`` deterministic, ``<0`` from
       the clock.

See :ref:`Thermal fluctuations`.

----------------------------------------
FMM demagnetization
----------------------------------------

.. list-table::
   :widths: 18 24 14 44
   :header-rows: 1

   * - Matlab
     - Python
     - Default
     - Description
   * - ``use_fmm``
     - ``use_fmm``
     - ``0``
     - Enable the FMM demagnetization path. Requires ``USE_FMM3D=1``.
   * - ``fmm_cells``
     - ``fmm_cells_per_node``
     - ``100``
     - Maximum cells in a leaf node before splitting.
   * - ``fmm_eps``
     - ``fmm_eps``
     - ``1e-4``
     - Requested FMM accuracy, which sets the plane-wave expansion order.
   * - ``fmm_nterms``
     - ``fmm_nterms``
     - ``-1``
     - Multipole expansion order. Negative means derive it from ``fmm_eps``.
   * - ``ifunif``
     - ``ifunif``
     - ``1``
     - Tree type. Must be 1 (uniform tree).
   * - ``nlmin``, ``nlmax``
     - ``nlmin``, ``nlmax``
     - ``1``, ``5``
     - Minimum and maximum level of the octree.
   * - ``fmm_short``
     - ``allow_fmm_short_circuit``
     - ``1``
     - Allow falling back to the direct calculation for small problems.
   * - ``fmm_min_n``
     - ``fmm_min_n``
     - ``20000``
     - Cell-count threshold above which FMM is used.

See :ref:`Demag field - FMM`.

----------------------------------------
Output, tracing and timing
----------------------------------------

.. list-table::
   :widths: 18 22 16 44
   :header-rows: 1

   * - Matlab
     - Python
     - Default
     - Description
   * - ``ReturnHall``
     - ``usereturnhall``
     - ``0``
     - Return the individual effective-field terms.
   * - ``dummy_run``
     - ``dummy_run``
     - ``0``
     - Build everything but skip the time integration.
   * - ``log_dir``, ``N_log_dir``
     - ``log_dir``
     - ``'logs'``
     - Directory for the trace and timing logs. Matlab:
       ``setLogDirFilename``.
   * - ``timer_log``, ``N_timer_log``
     - ``timer_log_file``
     - ``'timing.log'``
     - Timing log file. Matlab: ``setTimerLogFilename``.
   * - ``trace_log``, ``N_trace_log``
     - ``trace_log_file``
     - ``'trace.log'``
     - Trace log file. Matlab: ``setTraceLogFilename``.
   * - ``window_ena``
     - ``window_enabled``
     - ``1``
     - Flush timing statistics at intervals rather than only at the end.
   * - ``window_int``
     - ``window_interval``
     - ``30.0``
     - Timing window length [s].
   * - ``trace_ena``
     - ``trace_enabled``
     - ``0``
     - Enable the execution trace. Significant performance cost.
   * - ``flush_each``
     - ``flush_each``
     - ``1``
     - Flush the trace file after every entry.
   * - ``trace_verb``
     - ``trace_verbose``
     - ``1``
     - Only log trace events with at least this verbosity.

See :ref:`Micromagnetic output` and
:ref:`Performance, Timing & Trace Logging`.

----------------------------------------
Matlab-only convenience fields
----------------------------------------

The following fields exist on the Matlab problem object but are not passed to
Fortran. They are used by the example scripts to steer plotting and saving:

``SaveTheResult``, ``ShowTheResult``, ``SolverType``, ``DirectoryFilename``,
``SimulationName``, ``FileName``, ``HextFct``, ``FFTdims``, ``ExternalMesh``,
``MeshType``, ``ExternalMeshFileName``, ``DemagTensorFileName`` and
``exch_mat``.
