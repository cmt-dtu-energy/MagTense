Micromagnetic examples and validation
========================================

The examples below ship with MagTense and are the recommended starting point.
The Matlab and the Python versions of each test are written to mirror each
other, so the same physics is checked from both interfaces.

----------------------------------------
muMag standard problems
----------------------------------------

.. list-table::
   :widths: 16 42 42
   :header-rows: 1

   * - Problem
     - Matlab
     - Python
   * - Standard problem 2
     - `Standard_problem_2.m <https://github.com/cmt-dtu-energy/MagTense/blob/master/matlab/examples/Micromagnetism/mumag_micromag_Std_problem_2/Standard_problem_2.m>`_
     - \-
   * - Standard problem 3
     - `Standard_problem_3.m <https://github.com/cmt-dtu-energy/MagTense/blob/master/matlab/examples/Micromagnetism/mumag_micromag_Std_problem_3/Standard_problem_3.m>`_,
       ``Standard_problem_3_unstructured_cart.m``,
       ``Standard_problem_3_tetra.m``
     - `std_problem_3.py <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/examples/micromagnetism/mumag_micromag_Std_problem_3/std_problem_3.py>`_
   * - Standard problem 4
     - `Standard_problem_4.m <https://github.com/cmt-dtu-energy/MagTense/blob/master/matlab/examples/Micromagnetism/mumag_micromag_Std_problem_4/Standard_problem_4.m>`_
     - `std_problem_4.py <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/examples/micromagnetism/mumag_micromag_Std_problem_4/std_problem_4.py>`_
   * - Standard problem 6
     - `Standard_problem_6.m <https://github.com/cmt-dtu-energy/MagTense/blob/master/matlab/examples/Micromagnetism/mumag_micromag_Std_problem_6/Standard_problem_6.m>`_
     - `std_problem_6.py <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/examples/micromagnetism/mumag_micromag_Std_problem_6/std_problem_6.py>`_

What each of these exercises do:

* **Standard problem 2** - the quasi-static hysteresis loop of a rectangular
  bar, computed with the **explicit** solver and swept over the width of the
  bar in units of the exchange length. This is the example to copy for a
  hysteresis calculation.
* **Standard problem 3** - the flower/vortex energy crossover of a cube as a
  function of its size in units of the exchange length. It is run on a uniform
  grid, on an unstructured Cartesian mesh and on a tetrahedral mesh, so it is
  also the example to copy for either kind of unstructured mesh.
  ``Standard_problem_3_tetra.m`` is the shortest illustration of a tetrahedral
  problem: it meshes the cube with ``CreateTetraMesh``, hands the mesh over with
  a single call to ``setMicroMagGridTetrahedron``, and lets MagTense do the rest.
* **Standard problem 4** - the switching dynamics of a thin film under a
  reversed field, compared against the published mean solutions. It is a
  two-stage calculation: first an s-state is relaxed, then that state is used
  as the initial condition for the **dynamic** run.
* **Standard problem 6** - domain-wall pinning at a phase boundary, which is
  the test of the spatially varying :math:`A_0`, :math:`K_0` and :math:`M_s`
  and of the modified exchange stencil.

----------------------------------------
Feature tests
----------------------------------------

These live in ``matlab/examples/Micromagnetism/MagTense_tests`` and in
``python/examples/micromagnetism/MagTense_tests``, and each of them compares
against an analytical result rather than against a reference simulation.

.. list-table::
   :widths: 26 74
   :header-rows: 1

   * - Test
     - What it checks
   * - ``macrogeometry_PBC_test``
     - Periodic boundaries through the macrogeometry method. Two cubes are
       repeated along one axis and the critical spacing at which the effective
       anisotropy vanishes is compared with the analytical prediction. The
       whole setup is rotated so that the axis of periodicity becomes x, y and
       z in turn, which validates all three directions. Two further spacings
       act as controls. See :ref:`Macrogeometry`.
   * - ``periodic_exchange_test``
     - Periodic exchange coupling, on a uniform grid, on an unstructured mesh
       and on a grain mesh. Exchange-coupled moments across the periodic
       boundary must end up identical. See
       :ref:`Periodic exchange boundaries`.
   * - ``shape_correction_test``
     - Both the shape anisotropy (shape-dependent demagnetisation) and the magnetocrystalline anisotropy give uniaxial anisotropy energies. The shape correction field rewrites the shape anisotropy so it precisely cancels the magnetocrystalline contribution. Consequently, in a successful test the magnetization stays put. See
       :ref:`Sample shape correction`.
   * - ``temperature_test``
     - Thermal fluctuations against the analytical angular diffusion of
       non-interacting moments. See :ref:`Thermal fluctuations`.

----------------------------------------
Other examples
----------------------------------------

* `CoFe_nanopillar.m <https://github.com/cmt-dtu-energy/MagTense/tree/master/matlab/examples/Micromagnetism/CoFe%20nanopillar>`_
  - a nanopillar on an externally generated unstructured prism mesh.
* ``python/examples/micromagnetism/adaptive_hysteresis/adaptive_hysteresis_minimal.py``
  - the smallest possible adaptive hysteresis run, with no plotting or
  post-processing. ``adaptive_hysteresis_demo.py`` in the same folder adds a
  comparison against a fixed-step loop.

----------------------------------------
Running the test suites
----------------------------------------

Both interfaces ship a runner that executes the examples, turns each of them
into a set of pass/fail checks and prints a summary table.

In Matlab, from ``matlab/util``:

.. code-block:: matlab

    testMagTenseFunctions

The tests to run can be restricted with the environment variables
``MAGTENSE_TESTS`` and ``MAGTENSE_SKIP``, and the CUDA and CVODE variants are
selected with ``MAGTENSE_TEST_CUDA`` and ``MAGTENSE_TEST_CVODE``. Standard
problem 6 alone accounts for most of the running time.

In Python, from ``python/examples/micromagnetism``:

.. code-block:: bash

    python testMagTenseFunctions.py

.. code-block:: bash

    python testMagTenseFunctions.py --list

.. code-block:: bash

    python testMagTenseFunctions.py --tests temperature_test,std_problem_4

.. code-block:: bash

    python testMagTenseFunctions.py --skip std_problem_6

``--include-slow`` adds standard problem 3, which is very slow. Each example
lives in its own directory and is run from there, so the figures and timer logs
it produces land beside it, and the overview figure of the suite is written next
to ``testMagTenseFunctions.py``. There are in addition ``pytest`` suites in
``python/examples/magnetostatics`` and
``python/examples/micromagnetism/MagTense_tests``.
