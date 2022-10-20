Micromagnetism
========================================

========================================
Basic input
========================================
The matlab file `DefaultMicroMagProblem.m <https://github.com/cmt-dtu-energy/MagTense/blob/micromag_python_rebase/matlab/util/DefaultMicroMagProblem.m>`_ contains an updated list of all parameters that can be specified in a micromagnetic problem.

For now, the micromagnetic model automatically generates a rectangular grid.

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
An example of how to run a dynamic simulation is given `here <https://github.com/cmt-dtu-energy/MagTense/blob/micromag_python_rebase/matlab/examples/Micromagnetism/mumag_micromag_Std_problem_2/Standard_problem_2.m>`_.

========================================
Hysteresis loop simulation
========================================
An example of how to run a hysteresis loop simulation is given `here <https://github.com/cmt-dtu-energy/MagTense/blob/micromag_python_rebase/matlab/examples/Micromagnetism/mumag_micromag_Std_problem_4/Standard_problem_4.m>`_.
