# MagTense initial version 2019

MagTense consists of both a magnetostatic and a micromagnetism calculation framework

The magnetostatic framework is fully implemented in Fortran and has a Matlab MEX interface, as well as python interface (in progress)

The micromagnetism model currently only exists in Matlab, but utilizes the magnetostatic framework for calculating the demagnetization tensor.

Examples of how to calculate magnetostatic problems can be found in matlab\examples .

The micromagnetism code is currently in full development. An initial example is file "Script_3D_TestDynamicsStdProbl4.m" in the micromagnetism folder.

Currently the code is in massive development. The main features being worked on at the moment are:
- Proper code documentation
- Example code for both Matlab and python
- A transfer of the micromagnetism model to Fortran
