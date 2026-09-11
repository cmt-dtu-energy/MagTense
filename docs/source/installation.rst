Installation
==============================================

MagTense consists of a Fortran core with a Matlab MEX interface and a Python
interface. Neither interface requires building anything from source: prebuilt
MEX-files are attached to every release, and the Python interface is available
as a wheel on PyPI.

MagTense is tested on Linux and on Windows 11. macOS is not supported at the
moment.

==============================================
Python
==============================================

----------------------------------------------
Installing the Python package
----------------------------------------------

Requires **Python** :math:`\geq` **3.12** - wheels are published for Python 3.12, 3.13 and 3.14:

.. code-block:: bash

    pip install magtense

The wheel ships the compiled Fortran core together with the Intel MKL and
Intel Fortran runtime it needs, so no compiler is required.

For GPU support, the CUDA runtime libraries have to be present as well:

.. code-block:: bash

    pip install nvidia-cuda-runtime nvidia-cublas nvidia-cusparse nvidia-nvjitlink

----------------------------------------------
Usage
----------------------------------------------

The magnetostatic framework lives in ``magtense.magstatics`` and the
micromagnetic framework in ``magtense.micromag``:

.. code-block:: python

    from magtense.magstatics import Tiles, run_simulation
    from magtense.micromag import MicromagProblem

Examples for both are in
`python/examples <https://github.com/cmt-dtu-energy/MagTense/tree/master/python/examples>`_.

==============================================
Matlab
==============================================

----------------------------------------------
Prerequisites
----------------------------------------------

* Matlab :math:`\geq` 2023a
* An installation of the CUDA toolkit, if the CUDA-enabled MEX-files are to be
  used

MagTense is directly usable in Matlab by downloading the already compiled
`MEX-files <https://github.com/cmt-dtu-energy/MagTense/releases>`_, which are
provided for both Windows and Linux. Add the folder holding the MEX-files and
the ``matlab/util`` folder to the Matlab path:

.. code-block:: matlab

    addpath('MagTense/matlab/MEX_files');
    addpath('MagTense/matlab/util');

The MEX-files are

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - MEX-file
     - Purpose
   * - ``MagTenseLandauLifshitzSolver_mex``
     - Micromagnetic solver, CUDA-enabled
   * - ``MagTenseLandauLifshitzSolverNoCUDA_mex``
     - Micromagnetic solver without CUDA
   * - ``IterateMagnetization_mex``
     - Iterate the magnetization of soft magnetic tiles
   * - ``getHFromTiles_mex``
     - H-field from a set of tiles
   * - ``getNFromTile_mex``
     - Demagnetization tensor of a single tile
   * - ``getMagForce_mex``
     - Magnetic force on a tile

Examples are in
`matlab/examples <https://github.com/cmt-dtu-energy/MagTense/tree/master/matlab/examples>`_.

==============================================
Building from source
==============================================

Building is only necessary in order to modify the Fortran core, or to enable an
optional component that the prebuilt binaries do not carry. The optional
components are selected with make flags:

.. list-table::
   :widths: 24 76
   :header-rows: 1

   * - Flag
     - Effect
   * - ``USE_CUDA=1``
     - GPU acceleration of the demagnetization field through CUDA.
   * - ``USE_CVODE=1``
     - The CVODE time integrator from SUNDIALS, see
       :ref:`ODE solver settings`.
   * - ``USE_FMM3D=1``
     - The FMM demagnetization path, see :ref:`Demag field - FMM`.
   * - ``USE_MATLAB=1``
     - Compile the Fortran objects needed by the Matlab MEX-files.

A typical Linux build of the Python module with everything enabled is

.. code-block:: bash

    conda env create -n magtense-env -f python/.build/env-313-linux.yml
    conda activate magtense-env
    make python USE_CUDA=1 USE_CVODE=1 USE_MATLAB=0 USE_FMM3D=0
    python -m pip install -e ./python

The environment files come in one variant per Python version and platform,
``env-312``, ``env-313`` and ``env-314`` times ``-linux`` and ``-win``.

The full, up-to-date build instructions - including the Conda environment
files, the CVODE and CUDA prerequisites, the Windows toolchain and the
Visual Studio solution ``MagTense.sln`` - are kept next to the code, since they
change with the compiler versions:

* Python: `python/README.md <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/README.md>`_
* Matlab: `matlab/README.md <https://github.com/cmt-dtu-energy/MagTense/blob/master/matlab/README.md>`_
* The `GitHub workflow files <https://github.com/cmt-dtu-energy/MagTense/tree/master/.github/workflows>`_
  are always a working reference, since they are what builds and tests every
  commit.

For Matlab, the Fortran objects are built first and the MEX-files are then
produced from Matlab itself with
`buildMagTenseMEX.m <https://github.com/cmt-dtu-energy/MagTense/blob/master/matlab/buildMagTenseMEX.m>`_:

.. code-block:: matlab

    mex -setup FORTRAN
    buildMagTenseMEX('USE_RELEASE', true, 'USE_CUDA', true, 'USE_CVODE', false)
