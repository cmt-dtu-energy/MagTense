Performance, Timing & Trace Logging
===================================

MagTense includes a **Trace and Timing Module** (located in ``AuxMT``) for performance profiling. Note that because parameters are passed between Python and the Fortran core, boolean flags are represented as integers (**1 for True, 0 for False**).

Input Variables
---------------

.. list-table::
   :widths: 25 10 65
   :header-rows: 1

   * - Variable
     - Type
     - Description
   * - **log_dir**
     - str
     - Directory for **trace logs**.
   * - **timer_log_dir**
     - str
     - Directory for **timing logs**.
   * - **window_enable**
     - int
     - **1** to enable windowed timing; **0** to output only at the end.
   * - **window_interval**
     - float
     - Timing output frequency in **seconds**.
   * - **trace_enable**
     - int
     - **1** to enable execution trace. **Warning: Significant performance hit.**
   * - **flush_each**
     - int
     - **1** to flush file after every trace entry (safe but slow).
   * - **trace_verbose**
     - int
     - Only logs trace events with a verbosity $\ge$ this value.

Timing & Windowing
------------------
If ``window_enable = 1``, the module tracks elapsed time. Once the ``window_interval`` is exceeded, the next recorded event triggers a flush of the accumulated timing data for that window into the ``timer_log_dir``.

Developer Integration: Trace API
--------------------------------
When adding new Fortran code, use the following paired call structure:

.. code-block:: fortran

    call trace%begin(label, itimer=itime_counter, verbose=int_value)
    ! ... [computational code] ...
    call trace%end(label, itimer=itime_counter, verbose=int_value)

**Requirements:**
* The **label**, **itimer**, and **verbose** value **MUST** be identical in both calls.
* The ``itimer`` must be a saved integer within the function scope:
  
.. code-block:: fortran

    integer, save :: itimer = 0

Python Usage Example
--------------------

.. code-block:: python

    problem.log_dir = "./logs/trace"
    problem.timer_log_dir = "./logs/timing"
    
    # Timing window: Output stats every 60 seconds
    problem.window_enable = 1
    problem.window_interval = 60.0 
    
    # Trace configuration
    problem.trace_enable = 0  # Disabled for speed
    problem.trace_verbose = 2
    problem.flush_each = 0
    
    RunMicroMagSimulation(problem)