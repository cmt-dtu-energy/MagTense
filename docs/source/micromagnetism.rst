Micromagnetism
========================================

========================================
Basic input
========================================
The matlab file `DefaultMicroMagProblem.m <https://github.com/cmt-dtu-energy/MagTense/blob/master/matlab/util/DefaultMicroMagProblem.m>`_ and similarly the python file `micromag.py <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/src/magtense/micromag.py>`_ contains an updated list of all parameters that can be specified in a micromagnetic problem.

The micromagnetic model automatically generates a rectangular grid, if not provided with a grid. See the example files for how to do this.

The common parameters are:

* Grid resolution
* Gird size
* Exchange term constant
* Saturation magnetization
* Anisotropy constant
* Damping constant, :math:`\alpha`
* Precession constant, :math:`\gamma`


========================================
Dynamic simulations
========================================
An example of how to run a dynamic simulation in Matlab is given `here <https://github.com/cmt-dtu-energy/MagTense/blob/master/matlab/examples/Micromagnetism/mumag_micromag_Std_problem_2/Standard_problem_2.m>`_.

========================================
Hysteresis loop simulation
========================================
An example of how to run a hysteresis loop simulation in Matlab is given `here <https://github.com/cmt-dtu-energy/MagTense/blob/master/matlab/examples/Micromagnetism/mumag_micromag_Std_problem_4/Standard_problem_4.m>`_ and in python `here <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/examples/micromagnetism/std_problem_4.py>`_.
