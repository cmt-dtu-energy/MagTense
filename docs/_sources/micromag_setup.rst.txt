Setting up a micromagnetic problem
========================================

A micromagnetic problem is defined by a single problem object that is filled
with parameters and then handed to the Fortran solver. The Matlab and the
Python interfaces expose the same set of parameters, but they use slightly
different names; the naming used by Fortran is dictated by the subroutine
``getProblemFieldnames`` in ``MagTenseMicroMagIO.f90``.

The authoritative lists of parameters are

* Matlab: `DefaultMicroMagProblem.m <https://github.com/cmt-dtu-energy/MagTense/blob/master/matlab/util/DefaultMicroMagProblem.m>`_
* Python: `micromag.py <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/src/magtense/micromag.py>`_

and a side-by-side listing is given in :ref:`Micromagnetic parameter reference`.

----------------------------------------
Matlab problem setup
----------------------------------------

In Matlab a problem is created with ``DefaultMicroMagProblem(nx,ny,nz)``, which
returns an object with all parameters set to sensible defaults for a uniform
grid with ``nx*ny*nz`` cells. Parameters are then overwritten, and a number of
``set...`` helper methods map human-readable names onto the integer flags that
Fortran expects. Finally the object is converted to a plain struct and passed
to the MEX-file:

.. code-block:: matlab

    addpath('MagTense/matlab/MEX_files');
    addpath('MagTense/matlab/util');

    mu0 = 4*pi*1e-7;
    resolution = [36,9,1];

    problem = DefaultMicroMagProblem( resolution(1), resolution(2), resolution(3) );
    problem.grid_L = [500e-9,125e-9,3e-9];              % m

    % Solver and hardware options
    problem = problem.setMicroMagSolver( 'Dynamic' );
    problem = problem.setUseCuda( true );
    problem = problem.setUseCVODE( false );
    problem.nThreads = int32(8);

    % Material parameters
    problem.A0    = 1.3e-11;                            % J/m
    problem.Ms    = 8e5*ones(prod(resolution),1);       % A/m
    problem.K0    = zeros(prod(resolution),1);          % J/m^3
    problem.alpha = 4.42e3;                             % m/(A s)
    problem.gamma = 2.21e5;                             % m/(A s)

    % Initial state, normalised internally if it is not already
    problem.m0(:,1) = 1/sqrt(3);
    problem.m0(:,2) = 1/sqrt(3);
    problem.m0(:,3) = 1/sqrt(3);

    % Output times and the applied field as a function of time
    problem = problem.setTime( linspace(0,1e-9,200) );
    HystDir = 1/mu0*[-24.6,4.3,0]/1000;
    problem = problem.setHext( @(t) (t>-1)'.*HystDir, linspace(0,1e-9,2000) );

    % Solve
    solution = struct();
    solution = problem.MagTenseLandauLifshitzSolver_mex( struct(problem), solution );

The ``struct(problem)`` call performs a final consistency check before the data
is handed to Fortran: scalar values of ``Ms``, ``A0``, ``K0``, ``K1`` and
``K2`` are expanded to one value per cell, and ``m0`` is normalised if it is
not of unit length (a warning is issued). It also prints an estimate of the
memory needed by the demagnetization tensor.

``problem.MagTenseLandauLifshitzSolver_mex`` is a function handle that
``setUseCuda`` sets, pointing either at ``MagTenseLandauLifshitzSolver_mex`` or
at ``MagTenseLandauLifshitzSolverNoCUDA_mex``.

.. note::
   The constructor sets ``useCuda`` to 0 but does not populate that function
   handle, so ``setUseCuda`` has to be called at least once - with ``false``
   for a CPU run - before the solver can be invoked through
   ``problem.MagTenseLandauLifshitzSolver_mex``. Alternatively, call the
   MEX-file directly by name.

----------------------------------------
Python problem setup
----------------------------------------

In Python the problem is a ``MicromagProblem`` object. Most parameters are
constructor arguments, the remaining ones are plain attributes that can be
assigned after construction. The applied field is passed as a callable to the
run method rather than being stored on the problem:

.. code-block:: python

    import numpy as np
    from magtense.micromag import MicromagProblem

    mu0 = 4 * np.pi * 1e-7
    res = (36, 9, 1)

    problem = MicromagProblem(
        res=res,
        grid_L=[500e-9, 125e-9, 3e-9],
        grid_type="uniform",
        solver="dynamic",
        m0=1 / np.sqrt(3),
        A0=1.3e-11,
        Ms=8e5,
        K0=0.0,
        alpha=4.42e3,
        gamma=2.21e5,
        cuda=False,
        cvode=False,
    )

    h_ext = np.array([-24.6, 4.3, 0]) / 1000 / mu0

    def h_ext_fct(t):
        # Constant in time; must return an (nt_h_ext, 3) array
        return np.outer(np.ones_like(t), h_ext)

    t_out, M_out, pts, H_exc, H_ext, H_dem, H_ani, *_ = problem.run_simulation(
        t_end=1e-9, nt=200, fct_h_ext=h_ext_fct, nt_h_ext=2000
    )

``m0`` accepts a scalar (all components set to that value), an ``(ntot,3)``
array, or ``None``, in which case a random unit vector is drawn for every cell.
As in Matlab, ``m0`` is normalised on assignment and a warning is issued if it
was not already of unit length.

The three run methods are

* ``run_simulation(t_end, nt, fct_h_ext, nt_h_ext)`` - a single solve, used for
  both the dynamic and the explicit solver;
* ``run_hysteresis(H_ext)`` - a predefined sequence of applied fields, see
  :ref:`Hysteresis simulations`;
* ``run_hysteresis_adaptive(...)`` - a field sweep where the solver picks the
  field steps itself, see :ref:`Adaptive hysteresis`.

----------------------------------------
Units and conventions
----------------------------------------

.. list-table::
   :widths: 30 20 50
   :header-rows: 1

   * - Quantity
     - Unit
     - Comment
   * - Lengths, grid size
     - m
     - ``grid_L``, ``grid_pts``, ``grid_abc``
   * - Saturation magnetization :math:`M_s`
     - A/m
     - divide a value in tesla by :math:`\mu_0`
   * - Applied and effective fields
     - A/m
     - ``Hext``, and all returned ``H_*`` arrays
   * - Exchange constant :math:`A_0`
     - J/m
     - often written :math:`A_\mathrm{ex}`
   * - Anisotropy constants :math:`K_0, K_1, K_2`
     - J/m\ :sup:`3`
     -
   * - Damping :math:`\alpha`, precession :math:`\gamma`
     - m/(A s)
     - Landau-Lifshitz form, see :ref:`The equation that is solved`
   * - Time
     - s
     -
   * - Temperature
     - K
     -

The returned magnetization ``M`` is the reduced magnetization
:math:`\mathbf{m}`, i.e. it is normalised to unit length and must be multiplied
by :math:`M_s` to obtain A/m.
