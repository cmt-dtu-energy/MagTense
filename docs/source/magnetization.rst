Magnetization
========================================

::

    integer :: magnetType
    real,dimension(3) :: u_ea, u_oa1, u_oa2
    real :: mu_r_ea, mu_r_oa
    real :: Mrem
    real,dimension(3) :: M
    integer :: stateFunctionIndex

==============================================
Magnet Type
==============================================
* 1 = hard magnet
* 2 = soft magnet with state function
* 3 = soft magnet with constant permeability

The tile can have different magnetic properties.
It can be a hard magnet, which has a constant magnetization during the simulation.
Soft magnets can be defined either with a state function (hysteresis loop), 
which relates the magnetization of the tile to the currently H-field.
The third option describes a softmagnet with a constant permeability 
:math:`\mu \equiv \mu_0 ( 1 + \chi_m)` and therefore a linear relationship
for the magnetization with the H-field: :math:`M = \chi_m H`

==============================================
Easy axis and other axes
==============================================
The easy axis defines the direction of preferred magnetization in the magnetic tile.
With the other axes a orthonormal set is formed.
In the Matlab and python interfaces, the easy axis of the magnetization of 
is defined according to the global three-dimensional global coordinate 
system with two angles: polar angle :math:`\theta`, and azimuth :math:`\phi`, where 
:math:`\theta \in [0, \pi]` and :math:`\phi \in [0, 2 \pi]` .

.. math::
    \begin{pmatrix}
    u\_ea_x \\ 
    u\_ea_y \\ 
    u\_ea_z
    \end{pmatrix} 
    = 
    \begin{pmatrix}
    sin(\theta) cos(\phi)\\ 
    sin(\theta) sin(\phi) \\ 
    cos(\theta)
    \end{pmatrix}

.. math::
    \begin{pmatrix}
    u\_oa1_x \\ 
    u\_oa1_y \\ 
    u\_oa1_z
    \end{pmatrix} 
    = 
    \begin{pmatrix}
    sin(\theta) sin(\phi)\\ 
    -sin(\theta) cos(\phi) \\ 
    0
    \end{pmatrix}

.. math::
    \begin{pmatrix}
    u\_oa2_x \\ 
    u\_oa2_y \\ 
    u\_oa2_z
    \end{pmatrix} 
    = 
    \begin{pmatrix}
    \frac{1}{2}sin(2\theta) cos(\phi)\\ 
    \frac{1}{2}sin(2\theta) sin(\phi) \\ 
    -\sin^{2}(\theta)
    \end{pmatrix}

The permeability of the easy axis :math:`\mu_{r\_ea}`
and the permeability of the other axes :math:`\mu_{r\_oa}`
can be specified individually.
.. note:: Its default value is 1.

==============================================
Remanent Magnetization :math:`M_{rem}`
==============================================
The value of :math:`M_{rem}` specifies magnetization, which is 
remanent in the direction of the easy axis of the magnetic tile.
.. note:: Its unit is measured in :math:`A/m` in MagTense.

==============================================
Magnetization :math:`M`
==============================================
The magnetization of the tile in the global three-dimensional
global coordinate system.
Its quantity is calculated in MagTense and the effect on the
other existing tiles in the model is done through iterations.
For a hard (permanent) magnet with :math:`\mu_{r\_ea} = \mu_{r\_oa} = 1`,
the magnetization of the considered tile :math:`M = M_{rem} * u\_ea`.

==============================================
State function index
==============================================
For a soft magnet, whose magnetization is dependent on the H-field
in the form of a hysteresis loop, a 'state function'_ has to be delivered 
to MagTense. The state function index determines, which function of 
the possibly multiple given state functions belongs to that specific tile.
