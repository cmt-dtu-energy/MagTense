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
     - 3
     - Treats every row of the applied-field table as a separate, *constant*
       field and relaxes the magnetization to equilibrium for each of them in
       turn, starting from the state reached at the previous field. This is the
       quasi-static mode used for hysteresis loops. The equilibrium is found
       with the :ref:`Energy minimizer`. The minimizer cannot include the
       thermal field, so a finite temperature in any cell is an error, raised
       just before the Fortran call with a message pointing to ``ExplicitLL``.
   * - ``ExplicitLL``
     - 1
     - Reads the applied-field table as ``Explicit`` does, but finds the
       equilibrium at each field by integrating the Landau-Lifshitz equation in
       time. This is what ``Explicit`` did before the minimizer became its
       default, and the choice for thermal runs. Python
       ``solver="explicit_ll"``.

The names ``Minimizer`` and ``Implicit`` are no longer accepted; use
``Explicit``.

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

For the **explicit** solvers (``Explicit`` and ``ExplicitLL``)
the same table is read differently: the time column is ignored and each row is
one constant field to relax at. ``nt_h_ext`` is therefore the number of points
on the hysteresis curve, and ``t_end``/``nt`` only control how long each
relaxation is integrated (for the minimizer, only when it falls back to the
time integration) and how densely the hysteresis curve is sampled.

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
       state at the previous convergence time.
   * - ``conv_tol``
     - ``1e-4``
     - Convergence criterion: when the largest change of any magnetization
       component between two convergence times is below this, the integration
       stops and the converged state is held at the remaining output times.

The convergence test is what lets a relaxation to equilibrium stop when it is
done instead of always running to the last output time. With the default
``t_conv = 0`` it never fires. To use it, make ``t_conv`` the output times (or
a subset of them) and choose ``conv_tol``; in Python this is

.. code-block:: python

    problem.t_conv = np.linspace(0, t_end, nt)
    problem.nt_conv = nt
    problem.conv_tol = np.repeat(1e-6, nt)

The test is skipped in a thermal run, where the noise never settles. Both the
RKSuite and the CVODE driver apply it; CVODE only steps to the output times, so
a convergence time that is not also an output time is never visited there. The
times in ``t_conv`` *are* merged into the RKSuite integration grid, which
matters for the thermal field because it is redrawn once per grid time, see
:ref:`Thermal fluctuations`. Keep ``t_conv`` a subset of ``t`` unless you
specifically want to add integration times.

========================================
Energy minimizer
========================================

With ``solver = Explicit`` (Python ``solver="explicit"``) the equilibrium at
each constant applied field is found by minimizing the energy directly rather
than by integrating the Landau-Lifshitz equation with a large damping. The
method is steepest descent on the unit sphere with Barzilai-Borwein step
lengths, following Exl et al., *J. Appl. Phys.* **115**, 17D118 (2014). The
descent direction of cell :math:`i` is :math:`\mathbf{m}_i \times
(\mathbf{m}_i \times \mathbf{H}_i)`, the same vector that drives the damping
term of the Landau-Lifshitz equation, and the update

.. math::

    \mathbf{m}^{k+1} = \frac{(1 - \tau^2 |\mathbf{t}|^2/4)\,\mathbf{m}^k
      - \tau\, \mathbf{m}^k \times \mathbf{t}}{1 + \tau^2 |\mathbf{t}|^2/4},
    \qquad \mathbf{t} = \mathbf{m}^k \times \mathbf{H}^k ,

is the exact solution of the midpoint rule, so :math:`|\mathbf{m}| = 1` is
preserved to rounding. Each iteration costs **one** effective-field
evaluation, against about six per accepted step of the RK(4,5) time
integration, and the step length adapts to the local curvature instead of
being bounded by the stability limit of an explicit integrator, which for a
fine grid is set by the exchange stiffness whatever the distance from
equilibrium.

Convergence is declared when the largest torque over the cells,
:math:`\max_i |\mathbf{m}_i \times \mathbf{H}_i| / \max(M_s)`, falls below
``min_tol``. Two safeguards keep the non-monotone Barzilai-Borwein iteration
in check:

* no cell rotates by more than ``min_maxrot`` in one iteration, and
* the energy may not rise above the largest of the last ten iterates by more
  than a small fraction of the energy scale; if it would, the step is halved
  and retried.

If the iteration cap ``min_maxiter`` is hit, or the watchdog stalls, and
``min_fallback`` is set, the Landau-Lifshitz time integration is run once over
the requested time window from the current state and the minimizer is
restarted from where it ends up. This is what carries a hysteresis loop
through a switching event, where the state sits near a saddle point.

A vanishing torque is also what a saddle point looks like, and a symmetric
starting state sits on one: the canonical vortex of standard problem 3, or a
magnetization exactly antiparallel to the applied field. Steepest descent
converges onto such a point, whereas the time integration only leaves it
through rounding noise, slowly. With ``min_saddle_check = 1``, a converged
state is therefore nudged by a random tilt of about half
a degree, common to all cells with a smaller independent part per cell, and
relaxed again for up to a hundred iterations. A minimum keeps the energy above
the unperturbed value throughout and is returned unperturbed; a saddle lets the
energy fall below it, at which point the descent continues to the lower
minimum, which is then checked in the same way, up to three times. The check
costs up to a hundred field evaluations per applied field.

With ``min_saddle_check = 2``, the default, the question is answered rigorously
instead: the
lowest eigenvalue of the energy Hessian in the tangent space of the converged
state is computed by the Lanczos iteration, each step costing one field
evaluation, since the Hessian applied to a tangent displacement is the
finite difference of the torque along it (the energy is quadratic in the
magnetization, up to the anisotropy). The eigenvalue is returned in
``min_eig`` in units of the largest saturation magnetization: positive means a
minimum, negative a saddle, in which case the state is pushed along the
eigenvector, which is the direction of steepest descent out of the saddle,
and relaxed again. Near a switching event the eigenvalue goes to zero, so it
also tells how close a state is to switching. Along a field sweep the
eigenvector of the previous field starts the iteration, which then converges
in a handful of steps; a cold start takes a few tens. On a single-grain loop
the check costs less than the random nudge, and it finds saddles that the
nudge misses, such as the plain flower state of standard problem 3 near the
flower-vortex transition, which relaxes to the lower twisted flower.

.. list-table::
   :widths: 18 12 70
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``min_tol``
     - ``1e-5``
     - Convergence criterion on the largest relative torque.
   * - ``min_maxiter``
     - ``10000``
     - Maximum number of iterations per applied field.
   * - ``min_maxrot``
     - ``0.3``
     - Largest rotation of any cell in one iteration [rad].
   * - ``min_fallback``
     - ``1``
     - Fall back to the time integration when the minimizer stalls (1) or give
       up and report it (0).
   * - ``min_saddle_check``
     - ``2``
     - Nudge a converged state and relax again to make sure it is a minimum
       (1), compute the lowest eigenvalue of the energy Hessian instead (2, see
       below) or accept the state as it is (0). Matlab: ``min_saddle``.
   * - ``min_predictor``
     - ``1``
     - Start the minimizer at each applied field from the secant extrapolation
       of the two previous equilibria (1) instead of from the previous one (0),
       with the first step taken at the step length the previous field ended
       with. Costs no field evaluation; the extrapolation is skipped across a
       switching event. Matlab: ``min_pred``.

The minimizer works with every grid type, with CUDA and with FMM, because it
calls the same field routines as the time integration. It cannot be combined
with a finite temperature, since a stochastic field has no stationary point;
the solve stops with a message. The thermal, ``dynamic`` and time-dependent
``alpha`` settings are ignored by it.

The output has the layout of the time integration: the first output time holds
the starting state and every later one the converged state. In addition every
run, with either method, returns the energies and a few diagnostics per applied
field, see :ref:`Micromagnetic output`: ``n_feval`` is the number of
effective-field evaluations spent relaxing at that field, which is the fair
cost measure between the two methods, and ``min_status`` records whether the
minimizer converged directly, needed the fallback, or failed.

How much cheaper it is depends on the problem. On the shipped examples,
measured in field evaluations per applied field on a CPU:

.. list-table::
   :widths: 40 20 20 20
   :header-rows: 1

   * - Problem
     - LL relaxation
     - Minimizer
     - Ratio
   * - Single 10 nm grain, 5\ :sup:`3` cells, adaptive hysteresis loop,
       1 ns relaxation window
     - 28695
     - 1877
     - 15
   * - Same, LL window extended to 100 ns so that its switching field agrees
       with the minimizer
     - 240624
     - 1877
     - 128
   * - Standard problem 3, 10\ :sup:`3` cells, flower state
     - 988
     - 177
     - 5.6
   * - Standard problem 3, 10\ :sup:`3` cells, vortex state
     - 12660
     - 337
     - 38

Roughly a hundred of the minimizer's evaluations per applied field go into the
saddle check; without it the flower state takes 39 evaluations and the vortex
156. The energies agree to six digits in the vortex case and the minimizer
finds a marginally lower flower energy than the time integration did in its
10 ns window. On the grain the 1 ns time integration is not relaxed near the
switching field and overshoots the Stoner-Wohlfarth value by ten percent; the
minimizer lands within one field step of it. The script
`minimizer_vs_llg.py <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/examples/micromagnetism/minimizer/minimizer_vs_llg.py>`_
reproduces these numbers.

Standard problem 6 shows the same picture on a domain wall pinned at a phase
boundary. The reference values are static depinning fields, the example's time
ramp of 0.02 T per ns gives rate-dependent ones, and the minimizer visiting the
same 201 field values as constant fields lands within the 0.01 T field
resolution of the analytical value in every setting, at a few hundred times
fewer field evaluations:

.. list-table::
   :widths: 16 16 20 20 24
   :header-rows: 1

   * - Setting
     - Analytical [T]
     - LL time ramp [T]
     - Minimizer [T]
     - Field evaluations LL / minimizer
   * - ``akj``
     - 1.568
     - 1.590
     - 1.580
     - 2 872 407 / 18 979
   * - ``ak``
     - 1.089
     - 1.120
     - 1.110
     - 2 850 968 / 17 067
   * - ``aj``
     - 1.206
     - 1.260
     - 1.240
     - 3 049 920 / 17 205
   * - ``a``
     - 0.838
     - 0.870
     - 0.860
     - 2 845 165 / 12 866
   * - ``kj``
     - 1.005
     - 1.020
     - 1.010
     - 11 088 509 / 28 309
   * - ``k``
     - 0.565
     - 0.590
     - 0.570
     - 2 848 863 / 13 982

The time integration remains the method of choice for dynamics, for thermal
runs, and as a robust reference when the minimizer reports a failure.

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

``Explicit`` relaxes at each field with the energy minimizer. Use
``ExplicitLL`` to integrate the Landau-Lifshitz equation instead; see
:ref:`Energy minimizer` for the trade-off.

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
``hysteresis_solver='adaptive'`` together with ``solver='explicit'`` or
``solver='explicit_ll'``. In Matlab
the same parameters exist on the problem object (``adaptiveHext``,
``maxHextSteps``, ``H_start``, ``H_end``, ``dH_initial``, ``dH_min``,
``dH_max``, ``dH_grow``, ``dH_shrink``, ``dM_min``, ``dM_target``,
``dM_reject``, ``switch_refdH``, ``use_sw_ref``).

The sweep first relaxes the starting state ``m0`` at ``H_start`` and stores
the result as the first field state, so every entry of the output, including
the first, is an equilibrium (and the first entry of ``n_feval``,
``min_iter``, ``min_torque``, ``min_status`` and ``min_eig`` describes that relaxation).
Every step is measured against the last accepted state, and measuring the
first one against an ``m0`` far from equilibrium - for example ``m0`` along a
field well away from the easy axis and below the anisotropy field - would make
it look like a switch. There is therefore no need to relax ``m0`` in a
separate solve before the sweep.

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
* Recovery from the floor: once ``dH`` has been driven down to ``dH_min``
  (by rejections or by the ``dM_target`` rule, e.g. around a fast change of
  the magnetization) or to ``switch_refdH`` by the switch refinement, ``dH``
  also grows by ``dH_grow`` after every accepted step with
  :math:`\mathrm{d}M \le` ``dM_target``, not only below ``dM_min``. After a
  switch refinement this waits until the step across the sign change has
  been accepted. The first step above ``dM_target`` after ``dH`` has grown
  ends the recovery and divides ``dH`` by ``dH_grow``, back to the last step
  length that stayed within ``dM_target``, and the rules above take over
  again. Without it, a
  smooth stretch that follows would give a :math:`\mathrm{d}M` between
  ``dM_min`` and ``dM_target`` at the floor and keep the rest of the sweep
  there: in a hard-axis loop that is about a thousand steps of ``dH_min``.
  Until ``dH`` first reaches its floor the step control is exactly the rules
  above.

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
     - ``1.25``
     - Growth factor, must be larger than 1.
   * - ``dH_shrink``
     - ``0.5``
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
