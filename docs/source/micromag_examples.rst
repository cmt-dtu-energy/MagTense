Micromagnetic examples and validation
========================================

The examples below ship with MagTense and are the recommended starting point.
The Matlab and the Python versions of each test are written to mirror each
other, so the same physics is checked from both interfaces. The standard
problems take the same options in both languages, with the same defaults:

.. list-table::
   :widths: 30 30 40
   :header-rows: 1

   * - Matlab option
     - Python argument
     - Meaning
   * - ``use_CUDA``
     - ``cuda``
     - Use CUDA for the calculations (default true)
   * - ``use_CVODE``
     - ``cvode``
     - Use CVODE for the time integration (default false)
   * - ``ShowTheResult``
     - ``plotting``, ``figpath``
     - Show the result, or in Python save it to ``figpath``
   * - ``use_minimizer``
     - ``use_minimizer``
     - Relax with the :ref:`Energy minimizer` instead of the time integration
   * - ``use_adaptive``
     - ``use_adaptive``
     - Problem 2: adaptive field stepping instead of the fixed table
   * - ``mesh_type``
     - ``mesh_type``
     - ``'uniform'``, ``'unstructuredPrisms'`` or, in problem 3, ``'tetrahedron'``
   * - ``mesh_file``
     - ``mesh_file``
     - The text file with the unstructured Cartesian mesh
   * - ``mesh_res_param``
     - ``mesh_res_param``
     - Problem 3: tetrahedra per edge length
   * - ``use_AvgN``
     - ``use_avgn``
     - Problem 4: the averaged prism tensor for the demag field
   * - ``cart_dir``, ``TwoDsim``, ``TwoDsize``
     - ``cart_dir``, ``two_d_sim``, ``two_d_size``
     - Problem 6: the axis of the sample and the 2D strip

The unstructured Cartesian meshes are text files, one prism per row with its
centre and its side lengths, and both languages read the same files from
``documentation/examples_mumag_validation``.

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
     - `std_problem_2.py <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/examples/micromagnetism/mumag_micromag_Std_problem_2/std_problem_2.py>`_
   * - Standard problem 3
     - `Standard_problem_3.m <https://github.com/cmt-dtu-energy/MagTense/blob/master/matlab/examples/Micromagnetism/mumag_micromag_Std_problem_3/Standard_problem_3.m>`_
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
  function of its size in units of the exchange length. The ``mesh_type``
  option runs it on a uniform grid (``'uniform'``, the default), on an
  unstructured Cartesian mesh read from a text file
  (``'unstructuredPrisms'``) or on a tetrahedral mesh (``'tetrahedron'``), so
  it is also the example to copy for either kind of unstructured mesh. The
  tetrahedral branch is the shortest illustration of a tetrahedral problem: it
  meshes the cube with ``CreateTetraMesh`` (Matlab, PDE Toolbox) or
  ``magtense.utils.create_tetra_mesh`` (Python), hands the mesh over with a
  single call to ``setMicroMagGridTetrahedron`` or the ``grid_nod`` and
  ``grid_ele`` arguments of ``MicromagProblem``, and lets MagTense do the rest.
* **Standard problem 4** - the switching dynamics of a thin film under a
  reversed field, compared against the published mean solutions. It is a
  two-stage calculation: first an s-state is relaxed, then that state is used
  as the initial condition for the **dynamic** run. With
  ``mesh_type = 'unstructuredPrisms'`` it runs on an unstructured Cartesian
  mesh, and the exchange operator MagTense builds in the first stage is handed
  to the second, so the mesh is analysed only once.
* **Standard problem 6** - domain-wall pinning at a phase boundary, which is
  the test of the spatially varying :math:`A_0`, :math:`K_0` and :math:`M_s`
  and of the modified exchange stencil. ``cart_dir`` orients the sample along
  x, y or z, which must not change the result, and ``mesh_type`` runs the same
  chain of cells as unstructured prisms.

Standard problems 2, 3 and 6 take a ``use_minimizer`` switch in both languages
(``options.use_minimizer`` in Matlab, the ``use_minimizer`` argument in Python).
With it set, the equilibrium at each field is found by the
:ref:`Energy minimizer` instead of the Landau-Lifshitz time integration, and
the scripts print the number of effective-field evaluations spent either way,
which is the cost to compare. Problem 3 reads its energies from the ``E``
output of the solver in both cases. In problem 6 the switch replaces the time
ramp of the applied field by the same field values visited as constant
fields, so the depinning field it reports is the static one, which is what the
analytical values of the reference paper are; the ramp gives a rate-dependent
one. Problem 2 additionally takes ``use_adaptive``, which replaces the fixed
table of 40 fields by the :ref:`Adaptive hysteresis` stepping, refining the
field step to 0.0005 T across the coercive field.

The Python script
`minimizer_vs_llg.py <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/examples/micromagnetism/minimizer/minimizer_vs_llg.py>`_
runs the single-grain hysteresis loop and standard problem 3 with both the
Landau-Lifshitz relaxation and the :ref:`Energy minimizer`, and prints the
switching fields, the energies and the number of field evaluations each method
needed.

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
into a set of pass/fail checks and prints a summary table. The two run the
same tests - the magnetostatic field of every tile type against FEM, and the
micromagnetic tests and standard problems - under the same names, with the
same checks and limits.

In Matlab, from ``matlab/util``:

.. code-block:: matlab

    testMagTenseFunctions

The tests to run can be restricted with the environment variables
``MAGTENSE_TESTS`` and ``MAGTENSE_SKIP``, and the CUDA and CVODE variants are
selected with ``MAGTENSE_TEST_CUDA`` and ``MAGTENSE_TEST_CVODE``.
``MAGTENSE_INCLUDE_SLOW=1`` adds standard problem 3 by time integration,
which is very slow. Standard problem 6 accounts for most of the running time
of the rest.

In Python, from ``python/util``:

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
it produces land beside it, and the overview figure of the suite is written to
``python/util/results``, as the Matlab suite writes to ``matlab/util/results``.
There are in addition ``pytest`` suites in
``python/examples/magnetostatics`` and
``python/examples/micromagnetism/MagTense_tests``.
