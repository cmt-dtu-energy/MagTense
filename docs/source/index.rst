Welcome to MagTense!
==============================================

A magnetostatic and micromagnetic calculation framework.

The magnetostatic framework is fully implemented in Fortran and has both a
Matlab MEX interface and a Python interface. It calculates the demagnetization
tensor fully analytically for cylinders, pieces of cylinders, prisms, circular
pieces, tetrahedra, spheres and spheroids.

The micromagnetic framework solves the Landau-Lifshitz equation using that same
analytical demagnetization tensor. It supports uniform, unstructured-prism and
tetrahedral meshes, spatially varying material parameters, general
magnetocrystalline anisotropy, thermal fluctuations, periodic boundary
conditions and adaptive hysteresis calculations.

MagTense uses Intel MKL, and can optionally be built with CUDA for GPU
acceleration, with CVODE as the time integrator and with FMM3D for an
:math:`O(N)` demagnetization field.

| Github: `http://github.com/cmt-dtu-energy/MagTense <http://github.com/cmt-dtu-energy/MagTense>`_
| Webpage: `http://www.magtense.org <http://www.magtense.org>`_

Getting started
----------------------------------------------

.. code-block:: bash

    pip install magtense

For Matlab, download the compiled MEX-files from the
`releases <https://github.com/cmt-dtu-energy/MagTense/releases>`_. See
:ref:`Installation` for details, and
:ref:`Setting up a micromagnetic problem` for a first micromagnetic
calculation.

==============================================
Content
==============================================

.. toctree::
   :maxdepth: 2

   about
   installation
   calculations
   micromagnetism
   demag_fmm
   timing_and_trace
   documentation
   theory
   publications
   gallery
