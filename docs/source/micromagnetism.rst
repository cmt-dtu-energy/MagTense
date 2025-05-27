Micromagnetism
========================================

========================================
Basic input
========================================
The matlab file `DefaultMicroMagProblem.m <https://github.com/cmt-dtu-energy/MagTense/blob/master/matlab/util/DefaultMicroMagProblem.m>`_ and similarly the python file `micromag.py <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/src/magtense/micromag.py>`_ contains an updated list of all parameters that can be specified in a micromagnetic problem.

The micromagnetic model automatically generates a rectangular grid, if not provided with a grid. See the example files for how to do this.

The common parameters are:

* Grid resolution
* Grid size
* Exchange constant
* Saturation magnetization
* Anisotropy constant
* Damping constant, :math:`\alpha`
* Precession constant, :math:`\gamma`


========================================
Anisotropy constants
========================================
The anisotropy can be specified as a local uniaxial anisotropy, :math:`problem.K_0(i)` at every grid point along with an anisotropy vector, :math:`problem.u_{ea}(i,:)=[easyX(i),easyY(i),easyZ(i)]` in matlab notation. The anisotropy can also be specified as a cubic anisotropy with :math:`problem.K_1(i)` and :math:`problem.K_2(i)` at every grid point along with a local coordinate system:

.. math::
    problem.CrystalAxis(i,:,:)
    = 
    \begin{pmatrix}
    x_1 &  x_2 &  x_3 \\
	y_1 &  y_2 &  y_3 \\
	z_1 &  z_2 &  z_3 
    \end{pmatrix}
	
However, if the anisotropy is not uniaxial, a general matrix formulation with the anisotropy energy density given by [following Eq. (2) in https://doi.org/10.1088/1361-665X/aafff8]:

.. math:: 
    f_{an}(\mathbf{m}) = \alpha_{1,x}m_x^2 +  \alpha_{1,y}m_y^2 +  \alpha_{1,z}m_z^2 + \alpha_{11,x}m_x^4 + \alpha_{11,y}m_y^4 + \alpha_{11,z}m_z^4 \\
	+ \alpha_{12,x}m_y^2m_z^2 + \alpha_{12,y}m_z^2m_x^2 + \alpha_{12,z}m_x^2m_y^2 + \alpha_{111,x}m_x^6 + \alpha_{111,y}m_y^6 + \alpha_{111,z}m_z^6 \\
	+ \alpha_{112,x}m_x^4(m_y^2 + m_z^2) + \alpha_{112,y}m_y^4(m_z^2 + m_x^2) + \alpha_{112,z}m_z^4(m_x^2 + m_y^2) \\
	+ \alpha_{123}(m_x^2m_y^2m_z^2)

In this case the anisotropy constants are specified locally at each grid point in a 6x3 matrix as 

.. math::
    K_{0,arr}
    = 
    \begin{pmatrix}
    \alpha_{1,x} &   \alpha_{1,y} &   \alpha_{1,z} \\
	\alpha_{11,x} &  \alpha_{11,y} &  \alpha_{11,z} \\
	\alpha_{12,x} &  \alpha_{12,y} &  \alpha_{12,z} \\
	\alpha_{111,x} &  \alpha_{111,y} &  \alpha_{111,z} \\
	\alpha_{112,x} &  \alpha_{112,y} &  \alpha_{112,z} \\
	\alpha_{123,x} &  0 &  0
    \end{pmatrix}
	
In the above matrix notation, a uniaxial anisotropy in the z-direction is given by:

.. math::
   K_{0,arr}
    =  
	\begin{pmatrix}
	0 & 0 & K_0 \\
    0 & 0 & 0  \\
    0 & 0 & 0  \\
    0 & 0 & 0  \\
    0 & 0 & 0  \\
    0 & 0 & 0  
    \end{pmatrix}

and a cubic anisotropy is given by:

.. math::
   K_{0,arr}
    =  
	\begin{pmatrix}
    0 & 0 & 0 \\
    0 & 0 & 0  \\
    K_1 & K_1 & K_1 \\ 
    0 & 0 & 0 \\
    0 & 0 & 0 \\
    K_2 & 0 & 0 
    \end{pmatrix}

For the general anisotropy the crystal axis must also be specified in a format similar to the cubic anisotropy.

========================================
Dynamic simulations
========================================
An example of how to run a dynamic simulation in Matlab is given `here <https://github.com/cmt-dtu-energy/MagTense/blob/master/matlab/examples/Micromagnetism/mumag_micromag_Std_problem_2/Standard_problem_2.m>`_.

========================================
Hysteresis loop simulation
========================================
An example of how to run a hysteresis loop simulation in Matlab is given `here <https://github.com/cmt-dtu-energy/MagTense/blob/master/matlab/examples/Micromagnetism/mumag_micromag_Std_problem_4/Standard_problem_4.m>`_ and in python `here <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/examples/micromagnetism/std_problem_4.py>`_.
