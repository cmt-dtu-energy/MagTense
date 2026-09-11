Material parameters
========================================

All material parameters may be specified per cell, i.e. as an array with one
entry per grid cell, which is what makes it possible to model composite
materials, grain structures and phase boundaries. In Matlab a scalar value is
expanded to a full array by ``struct(problem)``; in Python a scalar passed to
the constructor is expanded on assignment.

.. list-table::
   :widths: 14 14 12 60
   :header-rows: 1

   * - Matlab
     - Python
     - Unit
     - Description
   * - ``Ms``
     - ``Ms``
     - A/m
     - Saturation magnetization. Enters the demagnetization field directly and
       scales the exchange and anisotropy fields.
   * - ``A0``
     - ``A0``
     - J/m
     - Exchange stiffness constant. Where two cells with different values
       share a face, the coupling between them is the harmonic mean of the
       two, which can be overridden as seen in
       :ref:`Exchange across a material interface`.
   * - ``K0``
     - ``K0``
     - J/m\ :sup:`3`
     - Uniaxial anisotropy constant, used together with the easy axis ``u_ea``.
   * - ``K1``, ``K2``
     - ``K1``, ``K2``
     - J/m\ :sup:`3`
     - Cubic anisotropy constants, used together with ``CrysAxis``.
   * - ``K0_arr``
     - ``K0_arr``
     - J/m\ :sup:`3`
     - General anisotropy expansion, an ``(ntot,6,3)`` array.
   * - ``CrysAxis``
     - ``CrysAxis``
     - \-
     - Local crystal coordinate system, an ``(ntot,3,3)`` array.
   * - ``u_ea``
     - ``u_ea``
     - \-
     - Uniaxial easy-axis direction of every cell, an ``(ntot,3)`` array.
   * - ``alpha``
     - ``alpha``
     - m/(A s)
     - Damping constant of the Landau-Lifshitz equation.
   * - ``gamma``
     - ``gamma``
     - m/(A s)
     - Precession constant. Set to zero to remove precession.
   * - ``temperature``
     - ``T``
     - K
     - Per-cell temperature, see :ref:`Thermal fluctuations`.

Internally the solver forms the scaled coefficients
:math:`J = A_0/(\mu_0 M_s)`, :math:`K = K_0/(\mu_0 M_s)` and
:math:`M = M_s`, so that all field terms come out in A/m.

.. note::
   For the ``tetrahedron`` and ``unstructuredPrisms`` grids the exchange
   prefactor is normalised by the largest value of :math:`A_0` in the mesh, and
   the spatial variation of :math:`A_0` is instead carried inside the exchange
   operator itself. On a uniform grid the local :math:`A_0` enters the
   prefactor directly.

The exchange length that sets the required cell size follows from :math:`A_0`
and :math:`M_s`,

.. math::

    l_\mathrm{ex} = \sqrt{\frac{A_0}{\tfrac{1}{2}\mu_0 M_s^2}} .

========================================
Magnetocrystalline anisotropy
========================================

MagTense offers three ways of specifying the anisotropy. They are **mutually
exclusive in the following sense**: the uniaxial constant ``K0`` cannot be
combined with ``K1``, ``K2`` or ``K0_arr``. Doing so aborts the solve with

::

    MagTense: the uniaxial anisotropy K0 cannot be combined with K1, K2 or K0_arr
    MagTense: express the uniaxial term through K0_arr(:,1,:) instead, or set K0 to zero

``K1``/``K2`` and ``K0_arr`` *are* compatible - they are expressed in the same
rotated frame and are summed.

----------------------------------------
Uniaxial anisotropy
----------------------------------------

Specify a local anisotropy constant ``K0(i)`` at every grid point together with
an easy axis ``u_ea(i,:) = [easyX(i), easyY(i), easyZ(i)]``. The anisotropy
field is then

.. math::

    \mathbf{H}_\mathrm{ani} = \frac{2 K_0}{\mu_0 M_s}
        \left( \mathbf{u}_\mathrm{ea}\cdot\mathbf{m} \right) \mathbf{u}_\mathrm{ea},

i.e. a positive :math:`K_0` makes :math:`\mathbf{u}_\mathrm{ea}` an easy axis.
This path does not use ``CrysAxis``.

----------------------------------------
Cubic anisotropy
----------------------------------------

Specify ``K1(i)`` and ``K2(i)`` at every grid point together with a local
coordinate system

.. math::
    \mathrm{CrysAxis}(i,:,:)
    =
    \begin{pmatrix}
    x_1 &  x_2 &  x_3 \
    y_1 &  y_2 &  y_3 \
    z_1 &  z_2 &  z_3
    \end{pmatrix},

whose rows are the crystal axes of cell :math:`i` expressed in the global
coordinate system. The default is the three Cartesian axes. The energy density
follows the standard convention

.. math::

    E_\mathrm{cubic} = K_1 \left( m_x^2m_y^2 + m_y^2m_z^2 + m_z^2m_x^2 \right)
        + K_2\, m_x^2m_y^2m_z^2 ,

with :math:`m` expressed in the local crystal frame, so that :math:`K_1 > 0`
gives easy axes along :math:`\langle 100 \rangle`.

----------------------------------------
General anisotropy
----------------------------------------

If the anisotropy is neither uniaxial nor cubic, a general matrix formulation
is available. The anisotropy energy density is written as a polynomial
expansion [following Eq. (2) in https://doi.org/10.1088/1361-665X/aafff8,
generalized to independent coordinates]

.. math::
    E_\mathrm{an}(\mathbf{m}) = - \big[ \alpha_{1,x}m_x^2 +  \alpha_{1,y}m_y^2 +  \alpha_{1,z}m_z^2 + \alpha_{11,x}m_x^4 + \alpha_{11,y}m_y^4 + \alpha_{11,z}m_z^4 \
    + \alpha_{12,x}m_y^2m_z^2 + \alpha_{12,y}m_z^2m_x^2 + \alpha_{12,z}m_x^2m_y^2 + \alpha_{111,x}m_x^6 + \alpha_{111,y}m_y^6 + \alpha_{111,z}m_z^6 \
    + \alpha_{112,x}m_x^4(m_y^2 + m_z^2) + \alpha_{112,y}m_y^4(m_z^2 + m_x^2) + \alpha_{112,z}m_z^4(m_x^2 + m_y^2) \
    + \alpha_{123}(m_x^2m_y^2m_z^2) \big],

again with :math:`\mathbf{m}` in the local crystal frame given by
``CrysAxis``. The anisotropy field is
:math:`\mathbf{H}_\mathrm{ani} = -\frac{1}{\mu_0 M_s}\,\partial E_\mathrm{an}/\partial \mathbf{m}`.
The overall minus sign in front of the bracket is what makes a *positive*
coefficient an easy direction, consistent with the uniaxial convention above.

The coefficients are specified per grid point as a 6x3 matrix

.. math::
    K_{0,arr}(i,:,:)
    =
    \begin{pmatrix}
    \alpha_{1,x} &   \alpha_{1,y} &   \alpha_{1,z} \
    \alpha_{11,x} &  \alpha_{11,y} &  \alpha_{11,z} \
    \alpha_{12,x} &  \alpha_{12,y} &  \alpha_{12,z} \
    \alpha_{111,x} &  \alpha_{111,y} &  \alpha_{111,z} \
    \alpha_{112,x} &  \alpha_{112,y} &  \alpha_{112,z} \
    \alpha_{123} &  0 &  0
    \end{pmatrix}.

In this notation a uniaxial anisotropy along the local z-direction is

.. math::
   K_{0,arr}
    =
    \begin{pmatrix}
    0 & 0 & K_0 \
    0 & 0 & 0  \
    0 & 0 & 0  \
    0 & 0 & 0  \
    0 & 0 & 0  \
    0 & 0 & 0
    \end{pmatrix},

and a cubic anisotropy with the standard sign convention is

.. math::
   K_{0,arr}
    =
    \begin{pmatrix}
    0 & 0 & 0 \
    0 & 0 & 0  \
    -K_1 & -K_1 & -K_1 \
    0 & 0 & 0 \
    0 & 0 & 0 \
    -K_2 & 0 & 0
    \end{pmatrix},

which is exactly what the solver inserts when ``K1`` and ``K2`` are given
directly. Because the two are added, an arbitrary higher-order term may be
placed in ``K0_arr`` while the cubic part is still given through ``K1`` and
``K2``.

========================================
Time-dependent damping
========================================

The damping constant may be made a function of time. This is used to speed up
relaxation towards an equilibrium state: a large damping early on kills the
precession quickly, and it is then lowered towards the physical value.

If ``alpha`` is set to zero, the solver instead interpolates the tabulated
values in ``alphat``. In Matlab this table is filled by

.. code-block:: matlab

    problem = problem.setAlpha( @(t) alpha_of_t(t), t_alpha );

and in Python through the ``t_alpha`` and ``alpha_fct`` constructor arguments.
If ``alpha`` is non-zero it is used as a constant and the table is ignored.

.. note::
   The parameter ``MaxT0`` is accepted by both interfaces for backward
   compatibility but is not used by the current solver. Use the tabulated
   ``alphat`` above instead.
