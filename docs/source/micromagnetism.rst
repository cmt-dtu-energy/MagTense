Micromagnetism
========================================

MagTense solves the Landau-Lifshitz (LL) equation on a mesh of tetrahedral or rectangular prism 
cells, using the same fully analytical demagnetization tensor that
is used by the magnetostatic framework. The micromagnetic model is implemented
in Fortran (sub-project ``MagTenseMicroMag``) and is driven either from Matlab
through the ``MagTenseLandauLifshitzSolver_mex`` MEX-file, or from Python
through the :code:`magtense.micromag.MicromagProblem` class.

This part of the manual describes the input parameters and features of the
micromagnetic model. The pages below go through the model one topic at a time;
:ref:`Micromagnetic parameter reference` lists every parameter with its Matlab
name, its Python name and its default value.

.. toctree::
   :maxdepth: 2

   micromag_setup
   micromag_grid
   micromag_material
   micromag_exchange
   micromag_demag
   micromag_solvers
   micromag_pbc
   micromag_thermal
   micromag_output
   micromag_parameters
   micromag_examples

========================================
The equation that is solved
========================================

MagTense integrates the Landau-Lifshitz equation in the form

.. math::

    \frac{\partial \mathbf{m}}{\partial t} = -\gamma \, \mathbf{m} \times \mathbf{H}_\mathrm{eff}
    - \alpha \, \mathbf{m} \times \left( \mathbf{m} \times \mathbf{H}_\mathrm{eff} \right),

where :math:`\mathbf{m} = \mathbf{M}/M_s` is the reduced magnetization of a
cell, i.e. a unit vector. Both the precession constant :math:`\gamma` and the
damping constant :math:`\alpha` are given in units of :math:`\mathrm{m/(A\,s)}`,
i.e. this is the Landau-Lifshitz form and *not* the Gilbert form. They are
related to the dimensionless Gilbert damping :math:`\alpha_\mathrm{G}` by

.. math::

    \alpha = \frac{\gamma \, \alpha_\mathrm{G}}{1 + \alpha_\mathrm{G}^2}.

.. note::
   Setting :math:`\gamma = 0` removes the precession term. Several of the
   quasi-static examples do this - among them the relaxation stage of standard
   problem 4 and the adaptive hysteresis example - since only the relaxed state
   is of interest there and dropping the precession makes the relaxation much
   faster.

The effective field is the sum of five contributions, all in :math:`\mathrm{A/m}`,

.. math::

    \mathbf{H}_\mathrm{eff} = \mathbf{H}_\mathrm{ext} + \mathbf{H}_\mathrm{exc}
    + \mathbf{H}_\mathrm{dem} + \mathbf{H}_\mathrm{ani} + \mathbf{H}_\mathrm{th},

being the applied (external) field, the exchange field, the demagnetization
field, the magnetocrystalline anisotropy field and the stochastic thermal
field. Each term is described on its own page below. The thermal term is only
included when a non-zero temperature is specified, see
:ref:`Thermal fluctuations`.

All quantities are in SI units. Note in particular that fields and
magnetizations are in :math:`\mathrm{A/m}` and not in tesla; multiply by
:math:`\mu_0` to convert.