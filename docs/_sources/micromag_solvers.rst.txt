Solvers and time integration
========================================

----------------------------------------
Solver type
----------------------------------------

The solver is selected with ``setMicroMagSolver`` in Matlab and with the
``solver`` argument in Python.

.. list-table::
   :widths: 20 8 72
   :header-rows: 1

   * - Name
     - Value
     - Description
   * - ``Dynamic``
     - 2
     - Integrates the Landau-Lifshitz equation in time under a single,
       time-varying applied field. The field is linearly interpolated from the
       tabulated values at the requested time.
   * - ``Explicit``
     - 1
     - Treats every row of the applied-field table as a separate, *constant*
       field and relaxes the magnetization to equilibrium for each of them in
       turn, starting from the state reached at the previous field. This is the
       quasi-static mode used for hysteresis loops.
   * - ``Implicit``
     - 3
     - Not implemented. Selecting it prints a message and produces no field
       update.

``ProblemMod`` (Python ``prob_mode``) selects ``new`` (1) or ``old`` (2). Use
``new``, which is the default; ``old`` skips the allocation of the solution
arrays and only makes sense for continuing a solution inside the same process.

----------------------------------------
Time grid and applied field
----------------------------------------

Two independent time arrays are involved:

* the **output times** ``t`` (``nt`` entries), at which the solution is
  returned;
* the **applied-field table** ``Hext`` (``nt_Hext`` rows of
  :math:`[t, H_x, H_y, H_z]`).

In Matlab both are set through helper methods:

.. code-block:: matlab

    problem = problem.setTime( linspace(0,1e-9,200) );          % output times
    HystDir = 1/mu0*[-24.6,4.3,0]/1000;                         % A/m
    problem = problem.setHext( @(t) (t>-1)'.*HystDir, linspace(0,1e-9,2000) );

``setHext(fct, t_Hext)`` evaluates the function handle on the given time array
and stores the result together with the times, so the applied field can be an
arbitrary function of time.

In Python the output times follow from ``t_end`` and ``nt``. The field
function is passed to the run method, which evaluates it on ``nt_h_ext``
uniformly spaced times between 0 and ``t_end``:

.. code-block:: python

    result = problem.run_simulation(
        t_end=1e-9, nt=200, fct_h_ext=h_ext_fct, nt_h_ext=2000
    )

The field function must return an ``(nt_h_ext, 3)`` array.

For the **explicit** solver the same table is read differently: the time column
is ignored and each row is one constant field to relax at. ``nt_h_ext`` is
therefore the number of points on the hysteresis curve, and ``t_end``/``nt``
only control how long each relaxation is integrated and how densely the hysteresis curve is
sampled.

``setTimeDis`` (Python ``setTimeDis``, default 10) controls how often the
solver reports progress: a message is printed every ``setTimeDis``'th
integration step.

----------------------------------------
ODE solver settings
----------------------------------------

Two time integrators are available:

* **RKSuite** - the default explicit Runge-Kutta suite, always available.
* **CVODE** - from `SUNDIALS <https://computing.llnl.gov/projects/sundials/cvode>`_,
  enabled with ``setUseCVODE(true)`` in Matlab or ``cvode=True`` in Python.
  It requires a build with ``USE_CVODE=1`` and is the better choice for stiff
  problems; RKSuite prints a warning suggesting CVODE when it detects
  stiffness.

.. list-table::
   :widths: 18 12 70
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``tol``
     - ``1e-4``
     - Relative tolerance of the ODE solver.
   * - ``thres``
     - ``1e-6``
     - Threshold value: a solution component smaller in magnitude than this is
       treated as zero by RKSuite.
   * - ``t_conv``, ``nt_conv``
     - ``0``, ``1``
     - Times at which the solution is checked for convergence against the
       previous step.
   * - ``conv_tol``
     - ``1e-4``
     - Convergence criterion, i.e. the maximum allowed change in magnetization
       between two steps.

.. note::
   ``t_conv`` and ``conv_tol`` are accepted by both interfaces, but the
   early-exit convergence test inside the RKSuite driver is currently disabled
   in the source, and the CVODE path never receives them. The times in
   ``t_conv`` *are* merged into the integration grid, which matters for the
   thermal field because it is redrawn once per grid time, see
   :ref:`Thermal fluctuations`. Keep ``t_conv`` a subset of ``t`` unless you
   specifically want to add integration times.

========================================
Hysteresis simulations
========================================

A hysteresis loop is computed with the **explicit** solver: a sequence of
constant applied fields is stepped through and the equilibrium magnetization is
found at each of them, starting from the previous equilibrium state. Setting
``gamma = 0`` removes the precession, which makes each relaxation considerably
cheaper when only the equilibrium state matters.

In Matlab this is just the explicit solver with a field table that walks the
loop, as in the standard problem 2 example:

.. code-block:: matlab

    problem = problem.setMicroMagSolver( 'Explicit' );
    HystDir = 1/mu0*[1,1,1]/sqrt(3);
    problem = problem.setHext( @(t) HystDir.*t', linspace(MaxH,-MaxH,40) );
    problem = problem.setTime( linspace(0,40e-9,2) );

In Python the dedicated method ``run_hysteresis`` takes the field table
directly as an ``(n,4)`` array, and requires ``hysteresis_solver='static'``,
which is the default:

.. code-block:: python

    problem = MicromagProblem(res=res, solver='explicit', ...)
    H_ext = np.zeros((40, 4))
    H_ext[:, 1:4] = np.linspace(max_H, -max_H, 40)[:, None] * hyst_dir
    result = problem.run_hysteresis(H_ext)

----------------------------------------
Adaptive hysteresis
----------------------------------------

Choosing the field steps of a loop by hand is wasteful: most of the curve is
smooth and needs few points, while the switching region needs many. The
adaptive hysteresis solver walks from ``H_start`` to ``H_end`` and picks the
step length itself, based on how much the volume-averaged magnetization moved.

The whole accept/reject loop runs inside Fortran in a single call. It is
available from Python through ``run_hysteresis_adaptive`` and requires
``hysteresis_solver='adaptive'`` together with ``solver='explicit'``. In Matlab
the same parameters exist on the problem object (``adaptiveHext``,
``maxHextSteps``, ``H_start``, ``H_end``, ``dH_initial``, ``dH_min``,
``dH_max``, ``dH_grow``, ``dH_shrink``, ``dM_min``, ``dM_target``,
``dM_reject``, ``switch_refdH``, ``use_sw_ref``).

The step metric is the length of the change of the cell-averaged reduced
magnetization vector across the trial field step,

.. math::

    \mathrm{d}M = \left| \langle\mathbf{m}\rangle_\mathrm{trial}
                       - \langle\mathbf{m}\rangle_\mathrm{before} \right| .

The step-control logic is then

* :math:`\mathrm{d}M >` ``dM_reject`` and :math:`\mathrm{d}H >` ``dH_min``:
  the step is **rejected**, ``dH`` is multiplied by ``dH_shrink`` (never below
  ``dH_min``) and the step is retried from the same state.
* :math:`\mathrm{d}M >` ``dM_target``: the step is accepted but ``dH`` is
  reduced for the next step.
* :math:`\mathrm{d}M <` ``dM_min``: the step is accepted and ``dH`` is
  multiplied by ``dH_grow`` for the next step, capped at ``dH_max``.
* Optionally, if ``use_sw_ref`` is set, a step across which the mean
  magnetization *along the field direction* changes sign is rejected whenever
  ``dH`` exceeds ``switch_refdH``. This forces fine sampling right at the
  switching field, which is what a coercivity calculation needs.

At ``dH_min`` a large change is accepted rather than looping forever, and the
solver says so in its progress output. Too many rejected steps in total aborts
the run.

.. list-table::
   :widths: 22 18 60
   :header-rows: 1

   * - Argument
     - Default
     - Description
   * - ``H_start``, ``H_end``
     - \-
     - Start and end field vectors in A/m. They must differ.
   * - ``dH_initial``
     - \-
     - First step length in A/m.
   * - ``dH_min``, ``dH_max``
     - \-
     - Smallest and largest allowed step length in A/m.
   * - ``max_steps``
     - \-
     - Maximum number of accepted field states, excluding the initial one.
   * - ``dM_min``
     - ``1e-3``
     - Below this change the step is grown.
   * - ``dM_target``
     - ``1e-2``
     - Above this change the step is shrunk.
   * - ``dM_reject``
     - ``5e-2``
     - Above this change the step is rejected and retried.
   * - ``dH_grow``
     - ``1.5``
     - Growth factor, must be larger than 1.
   * - ``dH_shrink``
     - ``0.75``
     - Shrink factor, must be between 0 and 1.
   * - ``switch_refine_dH``
     - ``None``
     - When given, the largest step accepted across a sign change of the mean
       magnetization along the field, in A/m.

The output arrays are preallocated to ``max_steps`` in Fortran and sliced on
return, and the number of accepted steps is the last element of the returned
list. A minimal, self-contained example is
`adaptive_hysteresis_minimal.py <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/examples/micromagnetism/adaptive_hysteresis/adaptive_hysteresis_minimal.py>`_.

.. code-block:: python

    problem = MicromagProblem(
        res=[5, 5, 5],
        grid_L=[10e-9, 10e-9, 10e-9],
        solver="explicit",
        hysteresis_solver="adaptive",
        Ms=2.4 / mu0, K0=1e6, A0=7e-12,
        alpha=4000.0, gamma=0.0,
        usereturnhall=True,
    )
    problem.u_ea[:, :] = [0.0, 0.0, 1.0]

    result = problem.run_hysteresis_adaptive(
        H_start=1.0 / mu0 * field_direction,
        H_end=-2.0 / mu0 * field_direction,
        dH_initial=0.5 / mu0,
        dH_min=0.01 / mu0,
        dH_max=0.5 / mu0,
        max_steps=512,
        switch_refine_dH=0.01 / mu0,
    )
    n_accepted = result[-1]
