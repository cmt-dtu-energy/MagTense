# MagTense version 2020

MagTense consists of both a magnetostatic and a micromagnetism calculation framework.

The magnetostatic framework is fully implemented in Fortran and has a Matlab MEX interface, as well as python interface.

The micromagnetism model utilizes the magnetostatic framework for calculating the demagnetization tensor.

MagTense is directly useable in Matlab by downloading the already compiled MEX-files in matlab\MEX. The files are directly useable with no compilation required, although Matlab 2018a or greater is required. Examples of how to calculate magnetostatic and micromagnetic problems using the Matlab interface can be found in matlab\examples.  

MagTense supports both Intel MKL, CUDA and CVODE.

The webpage of the code is available at https://www.magtense.org
The TechManual on the code is available at https://cmt-dtu-energy.github.io/MagTense/index.html

Currently the code is in massive development. The main features being worked on at the moment are:
- Proper code documentation
- Solidifying the input/output of the micromagnetic model
- Benchmarking the code
