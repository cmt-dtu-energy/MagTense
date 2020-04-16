# MagTense version 2020

MagTense consists of both a magnetostatic and a micromagnetism calculation framework.

The magnetostatic framework is fully implemented in Fortran and has a Matlab MEX interface, as well as python interface.

The micromagnetism model utilizes the magnetostatic framework for calculating the demagnetization tensor.

Examples of how to calculate magnetostatic and micromagnetic problems can be found in matlab\examples .

Currently the code is in massive development. The main features being worked on at the moment are:
- Proper code documentation
- Solidifying the input/output of the micromagnetic model
- Benchmarking the code
