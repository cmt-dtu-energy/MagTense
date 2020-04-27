Micromagnetism
========================================

========================================
Basic input
========================================
The matlab file *DefaultMicroMagProblem.m* contains an updated list of all parameters that can be specified in a micromagnetic problem.

For now, the micromagnetic model automatically generates a rectangular grid.

The common parameters are:
* Grid resolution
* Gird size
* Exchange term constant
* Saturation magnetization
* Anisotropy constant
* Dampning constant, :math:`\alpha`
* Precession constant, :math:`\gamma`
  
========================================
Dynamic simulations
========================================
An example of how to run a dynamic simulation is given in matlab/examples/NIST_micromag_Std_problem_2

========================================
Hysteresis loop simulation
========================================
An example of how to run a hysteresis loop simulation is given in matlab/examples/NIST_micromag_Std_problem_4


