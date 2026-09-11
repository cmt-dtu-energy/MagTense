Thermal fluctuations
========================================

A finite temperature is modelled by including a stochastic field in the effective
field, i.e. by integrating the stochastic Landau-Lifshitz equation. The
temperature is specified per cell, so a temperature gradient or a partially
heated sample can be modelled.

.. list-table::
   :widths: 20 18 62
   :header-rows: 1

   * - Matlab
     - Python
     - Description
   * - ``temperature``
     - ``T``
     - Temperature of every cell in K. A scalar is expanded to all cells. The
       default is zero, i.e. no thermal field.
   * - ``rng_seed``
     - ``rng_seed``
     - Seed for the random number generator, see
       :ref:`Reproducible random numbers`.

The thermal field is switched on automatically whenever the temperature is
non-zero anywhere in the mesh (the test is ``T > 1e-15`` K), and the solver
reports ``Including thermal noise``. Two things then change:

* a Gaussian random field is added to :math:`\mathbf{H}_\mathrm{eff}`, redrawn
  once per step of the integration time grid;
* the magnetization is renormalised to :math:`|\mathbf{m}| = 1` after every
  step, which prevents the thermal kicks from slowly inflating or deflating the
  magnetization.

----------------------------------------
The stochastic field
----------------------------------------

Each Cartesian component of the thermal field in cell :math:`i` is drawn from a
normal distribution with zero mean and standard deviation

.. math::

    \sigma_i = \sqrt{\frac{2 k_B T_i \, \alpha_\mathrm{G}}
                          {\mu_0 \gamma \, V_i \, M_{s,i} \, \Delta t}},

where :math:`V_i` is the volume of the cell, :math:`\Delta t` is the spacing of
the integration time grid, :math:`k_B` is Boltzmann's constant and
:math:`\alpha_\mathrm{G}` is the dimensionless Gilbert damping, recovered from
the Landau-Lifshitz damping constant :math:`\alpha` as the root
:math:`\alpha_\mathrm{G} \leq 1` of
:math:`\alpha_\mathrm{G}^2 - (\gamma/\alpha)\alpha_\mathrm{G} + 1 = 0`.

Two consequences of the :math:`1/\sqrt{\Delta t}` factor are worth keeping in
mind:

* The prefactor assumes a **uniform** time grid. The step used is the total
  simulated time divided by the number of unique times in the integration grid,
  i.e. the union of the output times ``t`` and the convergence-check times
  ``t_conv``. Adding times through ``t_conv`` therefore changes the effective
  temperature unless those times are already in ``t``.
* Because the field scales as :math:`1/\sqrt{V_i}`, a coarser mesh gives
  smaller fluctuations per cell, as it should: each cell then represents a
  larger, more strongly averaged volume.

.. note::
   The thermal field requires the **dynamic** solver. The explicit solver
   computes equilibrium states at a sequence of constant applied fields and has
   no time axis for the fluctuations to live on. The prefactor is still
   allocated for all solver types, so no error is raised, but a thermal run
   with the explicit solver is not meaningful.

----------------------------------------
Reproducible random numbers
----------------------------------------

``rng_seed`` seeds the Fortran random number generator that both the thermal
field and the demagnetization-field noise (``CV``) draw from. Its three regimes
are:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Value
     - Behaviour
   * - ``0`` (default)
     - The compiler default sequence is used. It is identical in every process
       and keeps advancing between solves within one process, so runs are
       neither reproducible nor independent. This preserves the behaviour of
       earlier MagTense versions.
   * - ``> 0``
     - The generator is seeded deterministically from this value. The same
       seed reproduces a run exactly; different seeds give independent runs.
       This is what a reproducible test needs.
   * - ``< 0``
     - The generator is seeded from the system clock, i.e. a fresh realisation
       on every run. This is what independent Monte-Carlo samples need.

.. note::
   The seed is applied inside Fortran. Seeding NumPy in a Python script has no
   effect on the thermal field - only on quantities that Python itself draws,
   such as a randomly initialised ``m0``.

----------------------------------------
Validation
----------------------------------------

The ``temperature_test`` example (available for both Matlab and Python) checks
the implementation against theory. It simulates a collection of
non-interacting cells with no anisotropy, initially all magnetised along z, and
compares the resulting angular distribution with the analytical solution of the
corresponding diffusion problem,

.. math::

    P(t, \theta) = \sin\theta \sum_{n=0}^{\infty} \frac{2n+1}{2}
        P_n(\cos\theta) \, e^{-n(n+1) D_\mathrm{LLG} t},

and it extracts the diffusion constant from
:math:`\langle \cos\theta \rangle = e^{-2 D t}`. Since the cells are
independent, the N cells are N samples of the same random walk, so the
comparison is a Monte-Carlo one and its tolerances are set by sampling noise
rather than by solver accuracy.
