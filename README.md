# MagTense version 2021

MagTense consists of both a magnetostatic and a micromagnetism calculation framework.

The magnetostatic framework can calculate the magnetic field from objects shaped as cylinders, pieces of cylinders, prisms, circular pieces and tetrahedrons. This is done using a fully analytical calculation of the demagnetization tensor. The frameworkis fully implemented in Fortran and has both a Matlab MEX interface and a Python interface.

The micromagnetism framework solves the Landau–Lifshitz equation. The frameworkis fully implemented in Fortran and has a Matlab MEX interface, as well as an older Matlab implementation. The micromagnetism framework utilizes the magnetostatic framework for calculating the demagnetizating field.

MagTense is directly useable in Matlab on Windows by downloading the already compiled MEX-files in matlab\MEX with no compilation required. The files are directly useable with no compilation required, although Matlab 2018a or greater is required. Examples of how to calculate magnetostatic and micromagnetic problems using the Matlab interface can be found in matlab\examples. Documentation on how to use the Python interface can be found in the Techmanual below.

If you want to compile MagTense a Visual Studio project file for Windows, MagTense.sln, is available, as well as a Matlab function to build the MEX-files, buildMagTenseMEX.m. MagTense utilizies Intel MKL for the micromagnetic simlations and can also utilize CUDA and CVODE.

The webpage of the code is available at https://www.magtense.org
The TechManual on the code is available at https://cmt-dtu-energy.github.io/MagTense/index.html

Currently the code is in development. The main features being worked on at the moment are:
- Proper code documentation
- Non-uniform grids
