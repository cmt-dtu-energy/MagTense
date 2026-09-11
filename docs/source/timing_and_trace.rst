Performance, Timing & Trace Logging
===================================

MagTense includes a **trace and timing module** (in the ``AuxMT`` sub-project)
for performance profiling of the Fortran core. Note that because parameters are
passed between Matlab/Python and Fortran, boolean flags are represented as
integers: **1 for true, 0 for false**.

Trace and timing variables
--------------------------

.. list-table::
   :widths: 22 20 8 50
   :header-rows: 1

   * - Python
     - Matlab
     - Type
     - Description
   * - ``log_dir``
     - ``log_dir``
     - str
     - Directory that both log files are written to. Default ``"logs"``. In
       Matlab set it with ``setLogDirFilename``, which also stores the string
       length in ``N_log_dir``.
   * - ``timer_log_file``
     - ``timer_log``
     - str
     - Name of the timing log file. Default ``"timing.log"``. Matlab:
       ``setTimerLogFilename``.
   * - ``trace_log_file``
     - ``trace_log``
     - str
     - Name of the trace log file. Default ``"trace.log"``. Matlab:
       ``setTraceLogFilename``.
   * - ``window_enabled``
     - ``window_ena``
     - int
     - **1** to enable windowed timing, **0** to output only at the end.
       Default 1.
   * - ``window_interval``
     - ``window_int``
     - float
     - Timing output frequency in **seconds**. Default 30.
   * - ``trace_enabled``
     - ``trace_ena``
     - int
     - **1** to enable the execution trace. **Warning: significant performance
       cost.** Default 0.
   * - ``flush_each``
     - ``flush_each``
     - int
     - **1** to flush the file after every trace entry: safe but slow. Default
       1.
   * - ``trace_verbose``
     - ``trace_verb``
     - int
     - Only trace events with a verbosity :math:`\geq` this value are logged.
       Default 1.

.. note::
   The log directory is created if it does not exist, on both Windows and
   Unix-like systems. If the creation fails, the module writes
   ``TRACE: failed to create log directory`` to standard error and the run
   continues.

Timing & Windowing
------------------

If ``window_enabled = 1``, the module tracks elapsed time. Once
``window_interval`` is exceeded, the next recorded event triggers a flush of
the accumulated timing data for that window into the timing log. This gives a
running picture of where the time goes in a long simulation, rather than a
single summary at the end.

The instrumented regions correspond to the main phases of the solver, such as
``SolveLandauLifshitzEquation``, ``ComputeDemagfieldTensor``,
``ComputeExchangeTerm3D_Uniform``, ``dmdt_fct``, ``updateDemagfield``,
``updateExchangeTerms`` and ``updateAnisotropy``, so the log shows directly
whether a run is dominated by tensor construction or by the time integration.

Developer integration: trace API
--------------------------------

When adding new Fortran code, use the following paired call structure:

.. code-block:: fortran

    call trace%begin(label, itimer=itime_counter, verbose=int_value)
    ! ... [computational code] ...
    call trace%end(label, itimer=itime_counter, verbose=int_value)

**Requirements:**

* The **label**, **itimer** and **verbose** value **must** be identical in both
  calls.
* The ``itimer`` must be a saved integer within the function scope:

.. code-block:: fortran

    integer, save :: itimer = 0

Verbosity 1 is used for the top-level phases and 2 for the inner routines, so
``trace_verbose = 1`` gives a coarse picture and higher values progressively
more detail.

Python trace example
--------------------

.. code-block:: python

    problem.log_dir = "./logs"
    problem.timer_log_file = "timing.log"
    problem.trace_log_file = "trace.log"

    # Timing window: output stats every 60 seconds
    problem.window_enabled = 1
    problem.window_interval = 60.0

    # Trace configuration
    problem.trace_enabled = 0  # Disabled for speed
    problem.trace_verbose = 2
    problem.flush_each = 0

    result = problem.run_simulation(
        t_end=t_end, nt=nt, fct_h_ext=h_ext_fct, nt_h_ext=nt_h_ext
    )

Matlab trace example
--------------------

.. code-block:: matlab

    problem = problem.setLogDirFilename( 'logs' );
    problem = problem.setTimerLogFilename( 'timing.log' );
    problem = problem.setTraceLogFilename( 'trace.log' );

    problem.window_ena  = int32(1);
    problem.window_int  = 60.0;
    problem.trace_ena   = int32(0);
    problem.flush_each  = int32(1);
    problem.trace_verb  = int32(1);
