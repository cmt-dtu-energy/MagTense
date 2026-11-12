module TileConeTensor
    !::Demagnetization tensor of a homogeneously magnetized right circular cone of base radius Rb and height h,
    !::from R. Bjoerk, "The magnetic field of a homogeneously magnetized cone". The cone stands on its base: the
    !::centre of the base is at the origin, the axis is the z-axis and the apex is at (0,0,h).
    !::The tensor follows the MagTense convention H = N.M, i.e. it is minus the tensor of the paper.
    !::
    !::The paper has the apex at the origin and the base in the plane z = h. A point (x,y,z) of the frame used
    !::here is the point (x,y,h-z) of the paper, and under this reflection only the xz and yz elements of the
    !::tensor change sign (getN_cone). All other routines of the module work in the frame of the paper, with the
    !::half-opening angle theta = atan(Rb/h), st = sin(theta) = Rb/ell, ct = cos(theta) = h/ell and the slant
    !::length ell = sqrt(Rb**2 + h**2).
    !::
    !::All elements are elementary functions and complete elliptic integrals of the three kinds, except for
    !::one function, I0, the integral over the azimuth of the potential of a generator of the cone, which is
    !::evaluated with QUADPACK (qags). The elliptic integrals are evaluated through the Carlson symmetric forms
    !::of SpecialFunctions, with kc2 = 1 - m:
    !::  K(m) = rf(0,kc2,1),   K(m) - E(m) = m/3 rd(0,kc2,1),   Pi(n|m) = K(m) + n/3 rj(0,kc2,1,1-n).
    !::The characteristic n- of the paper lies in (0,1], so the last form is a sum of positive terms. The
    !::characteristic n+ is negative and grows without bound towards the mirror cone, where that form would be a
    !::difference of nearly equal terms. It is therefore transformed with (Abramowitz and Stegun 17.7.17)
    !::  Pi(n|m) = -n (1-m) / ((1-n)(m-n)) Pi(n2|m) + m/(m-n) K(m),   n2 = (m-n)/(1-n),
    !::which with n2 in (0,1) is again a sum of positive terms (cone_closed).
    !::
    !::Three sets of points need separate treatment:
    !::  1. Near the axis the closed form contains inverse powers of the distance to the axis and loses about two
    !::     digits per decade, so there the tensor is evaluated from the expansion of the potential about the axis
    !::     (cone_axis). The iteration evaluates the field of a tile at its own centre, which lies on the axis.
    !::  2. On the cone through the lateral surface, i.e. on the lateral surface itself and on its continuation
    !::     beyond the base, and on its mirror cone below the apex, single terms of the closed form diverge
    !::     although their sum is finite. Within a small angle (seen from the apex) of these cones the tensor is
    !::     extrapolated quadratically from three points at one, two and three times that angle, on the side of
    !::     the field point. A point on the lateral surface takes the limit from outside.
    !::  3. On the rim of the base and at the apex the field diverges logarithmically. Floors give a large but
    !::     finite tensor there.
    !::
    !::Accuracy, measured against the closed form evaluated in 80-digit arithmetic for cones with h/Rb = 0.25,
    !::1.5 and 8, at distances to the apex from 1e-9 to 60 slant lengths (absolute error of the worst element):
    !::  generic points                              1e-13
    !::  near the axis                               4e-12, 8e-11 for the slender cone
    !::  within 1e-4 rad of the lateral surface      3e-10
    !::  within 1e-4 rad of the mirror cone          6e-12
    !::  at 1e-3 slant lengths from the rim          2e-8
    !::The error near the axis and near the lateral surface grows as the cone becomes more slender.
    use SpecialFunctions, only: rf, rd, rj
    use integrationDataTypes
    use QuadPack, only: qags
    use TileTensorHelperFunctions, only: pi
    implicit none

    private
    public :: getN_cone

    !::The expansion about the axis is used instead of the closed form for rho < r_sw * a_conv, where a_conv is
    !::the radius of convergence of the expansion and r_sw = cone_axis_tol * (c/(st*a_conv))**(2/7), at most
    !::cone_axis_max, see getN_cone_apex.
    real, parameter :: cone_axis_tol = 7e-3, cone_axis_max = 0.1

    !::Angle, seen from the apex, to the cone through the lateral surface and to its mirror cone, below which
    !::the tensor is extrapolated from three points further away, see point 2 above. The angle is
    !::cone_surf_tol * s**(3/4), with s the angle over which the field varies (see getN_cone_apex), which
    !::goes to zero at the rim of the base, but not below cone_surf_min.
    real, parameter :: cone_surf_tol = 1.5e-4, cone_surf_min = 1e-9

    !::Floor on kc2 = 1 - m. On the rim of the base K(m) diverges logarithmically; rd requires its arguments to
    !::sum to at least 1e-25.
    real, parameter :: cone_ksq_min = 1e-24

    !::Distance to the apex, relative to the slant length, below which the tensor is evaluated on the axis at
    !::that distance below the apex.
    real, parameter :: cone_apex_tol = 1e-12

    !::Tolerances of the quadrature of I0, which stops at an error of max( epsabs, epsrel * |I0| )
    real, parameter :: cone_I0_epsabs = 1e-13, cone_I0_epsrel = 1e-13

    !::Below this parameter m the function P1 = K - 2/3 rd is summed from its power series
    real, parameter :: cone_m_series = 0.05

    contains

    !::The demagnetization tensor N (H = N.M) at the point (x,y,z) of a cone of base radius Rb and height h that
    !::stands on its base, with the centre of the base at the origin and the apex at (0,0,h).
    !::Points on the surface count as outside.
    subroutine getN_cone( Rb, h, x, y, z, N )
    real,intent(in) :: Rb, h, x, y, z
    real,dimension(3,3),intent(inout) :: N

        if ( Rb .le. 0. .OR. h .le. 0. ) then
            N(:,:) = 0.
            return
        endif

        call getN_cone_apex( Rb, h, x, y, h - z, N )

        !::The reflection z -> h - z between the two frames changes the sign of the xz and yz elements
        N(1,3) = -N(1,3)
        N(3,1) = -N(3,1)
        N(2,3) = -N(2,3)
        N(3,2) = -N(3,2)

    end subroutine getN_cone

    !::The tensor (H = N.M) in the frame of the paper: apex at the origin, base in the plane z = h.
    !::Selects between the expansion about the axis, the closed form, and the extrapolation near the two cones
    !::through the apex.
    subroutine getN_cone_apex( Rb, h, x, y, z, N )
    real,intent(in) :: Rb, h, x, y, z
    real,dimension(3,3),intent(inout) :: N
    real :: ell, st, ct, theta, rho, c, rho_rim, a_conv, r_sw, cphi, sphi, sin_m, sin_p, side, vth0, u
    real :: tol, tol_c, tol_m, s_c
    real,dimension(3,3) :: N1, N2, N3
    logical :: near_cone, near_mirror

        ell = sqrt( Rb**2 + h**2 )
        st = Rb / ell
        ct = h / ell

        rho = sqrt( x**2 + y**2 )
        c = sqrt( rho**2 + z**2 )           !< distance to the apex

        !::At the apex the field diverges logarithmically
        if ( c .lt. cone_apex_tol * ell ) then
            call cone_axis( Rb, h, ell, st, ct, 0., 0., -cone_apex_tol * ell, N )
            return
        endif

        !::The expansion about the axis converges out to the nearest singularity of the field of the axis
        !::point: the apex, the rim of the base and, for an axis point inside the cone, the lateral surface.
        !::With a_conv the distance to that singularity, the truncation error of the expansion grows as
        !::(rho/a_conv)**5, while the rounding error of the closed form grows as (c/(rho*st))**2. The two are
        !::equal at rho = r_sw * a_conv.
        rho_rim = sqrt( ( rho - Rb )**2 + ( z - h )**2 )
        a_conv = min( c, rho_rim )
        if ( z .gt. 0. .AND. z .lt. h ) a_conv = min( a_conv, z * st / ct )
        if ( a_conv .gt. 0. ) then
            r_sw = min( cone_axis_tol * ( c / ( st * a_conv ) )**( 2. / 7. ), cone_axis_max )
            if ( rho .lt. r_sw * a_conv ) then
                call cone_axis( Rb, h, ell, st, ct, x, y, z, N )
                return
            endif
        endif

        cphi = x / rho
        sphi = y / rho

        !::With vartheta the polar angle of the field point, sin_m = sin(vartheta - theta) vanishes on the cone
        !::through the lateral surface and sin_p = sin(vartheta + theta) on its mirror cone below the apex.
        !::The half-width of the band around each cone follows the angle s over which the field varies there:
        !::the truncation error of the extrapolation grows as (tol/s)**3 and the rounding error of the closed
        !::form as 1/tol, which balance for tol proportional to s**(3/4). For the cone through the lateral
        !::surface s is the smallest of the angle to the axis, to the rim of the base and to the plane of the
        !::base, the last of which also keeps the points used for a field point inside the cone below the base.
        !::For the mirror cone s is the angle to the cone itself and to the axis.
        sin_m = ( rho * ct - z * st ) / c
        sin_p = ( rho * ct + z * st ) / c
        s_c = min( st, rho_rim / c, abs( h - z ) / c )
        tol_c = max( cone_surf_tol * s_c**0.75, cone_surf_min )
        tol_m = cone_surf_tol * min( 1., 2. * st * ct )**0.75
        near_cone = abs( sin_m ) .lt. tol_c .AND. ( z * ct + rho * st ) .gt. 0.
        near_mirror = abs( sin_p ) .lt. tol_m .AND. ( z * ct - rho * st ) .lt. 0.

        if ( near_cone .OR. near_mirror ) then
            theta = atan2( Rb, h )
            if ( near_cone ) then
                !::side = +1 outside the cone, also for a point on it
                side = sign( 1., sin_m )
                vth0 = theta
                tol = tol_c
                u = abs( asin( sin_m ) ) / tol
            else
                !::vartheta - (pi - theta) = -asin( sin_p )
                side = sign( 1., -sin_p )
                vth0 = pi - theta
                tol = tol_m
                u = abs( asin( sin_p ) ) / tol
            endif
            !::Quadratic extrapolation from the polar angles vth0 + k*side*tol, k = 1,2,3, to the angle
            !::vth0 + u*side*tol of the field point, with u between 0 (on the cone) and 1
            call cone_closed_polar( Rb, h, ell, st, ct, c, vth0 + side * tol, cphi, sphi, N1 )
            call cone_closed_polar( Rb, h, ell, st, ct, c, vth0 + 2 * side * tol, cphi, sphi, N2 )
            call cone_closed_polar( Rb, h, ell, st, ct, c, vth0 + 3 * side * tol, cphi, sphi, N3 )
            N = 0.5 * ( u - 2. ) * ( u - 3. ) * N1 - ( u - 1. ) * ( u - 3. ) * N2 + 0.5 * ( u - 1. ) * ( u - 2. ) * N3
        else
            call cone_closed( Rb, h, ell, st, ct, rho, z, cphi, sphi, N )
        endif

    end subroutine getN_cone_apex

    !::The closed form at the point with the distance c to the apex, the polar angle vth and the azimuth given
    !::by its cosine and sine
    subroutine cone_closed_polar( Rb, h, ell, st, ct, c, vth, cphi, sphi, N )
    real,intent(in) :: Rb, h, ell, st, ct, c, vth, cphi, sphi
    real,dimension(3,3),intent(inout) :: N

        call cone_closed( Rb, h, ell, st, ct, c * sin( vth ), c * cos( vth ), cphi, sphi, N )

    end subroutine cone_closed_polar

    !::The closed form of the paper at the point with the cylindrical coordinates rho, z and the azimuth given by
    !::cphi = cos(phi), sphi = sin(phi), in the frame of the paper. Not valid on the axis, on the cone through
    !::the lateral surface or on its mirror cone.
    !::
    !::Notation of the paper: Aq = A = rho sin(theta), Bq = B = z cos(theta), g_m = g-, g_p = g+, e_m = e-,
    !::e_p = e+, h_m = h-, h_p = h+, k_m = k-, k_p = k+, r_far = varrho+, and calA, calB, calC, calD the four
    !::azimuthal functions. The differences that vanish on the two cones are formed without cancellation:
    !::  t_m = c - B - A = (rho ct - z st)**2 / (c + B + A),   g_p = c + B - A = (rho ct + z st)**2 / g_m,
    !::  e_m = sqrt( t_m g_m ),   e_p = sqrt( g_p (c + B + A) ).
    !::The two elliptic integrals of the third kind enter only as Pi(n-|m)/g- and Pi(n+|m)/g+, which are
    !::  Pi(n-|m)/g- = K/g- + 2 A rj(0,kc2,1,t_m/g-) / (3 g-**2)
    !::  Pi(n+|m)/g+ = K (2 A kc2 + m c_p) / (c_p (m g+ + 2 A)) + 2 A kc2 rj(0,kc2,1,kc2 g+/c_p) / (3 c_p**2)
    !::with c_p = c + B + A, the second from the transformation given at the top of the module. The complete
    !::integral of the second kind is eliminated with E = K - m/3 rd, which removes the terms in 1/m**2 of the
    !::paper and leaves P1 = K - 2/3 rd, a function of order m.
    subroutine cone_closed( Rb, h, ell, st, ct, rho, z, cphi, sphi, N )
    real,intent(in) :: Rb, h, ell, st, ct, rho, z, cphi, sphi
    real,dimension(3,3),intent(inout) :: N
    real :: c, Aq, Bq, c_p, g_m, g_p, t_m, e_m, e_p, r_far, kc2, m, elK, elD, P1, Qm, Qp
    real :: h_m, h_p, k_m, k_p, calA, calB, calC, calD, I0, Pq, Qq, Rq, fac, c2p, s2p
    real :: abserr
    integer(8) :: neval, ier
    type(dataCollectionBase), target :: dat
    procedure (f_int_dat), pointer :: f_ptr => null ()

        c = sqrt( rho**2 + z**2 )
        Aq = rho * st
        Bq = z * ct

        c_p = c + Bq + Aq
        g_m = c - Bq + Aq
        t_m = ( rho * ct - z * st )**2 / c_p
        g_p = ( rho * ct + z * st )**2 / g_m
        e_m = sqrt( t_m * g_m )
        e_p = sqrt( g_p * c_p )

        !::Parameter of the ring formed by the rim of the base; r_far is the distance to its far point
        r_far = sqrt( ( rho + Rb )**2 + ( z - h )**2 )
        kc2 = ( ( rho - Rb )**2 + ( z - h )**2 ) / r_far**2
        if ( kc2 .lt. cone_ksq_min ) kc2 = cone_ksq_min
        m = 1 - kc2

        elK = rf( 0., kc2, 1. )
        elD = rd( 0., kc2, 1. )
        P1 = cone_P1( m, elK, elD )

        !::Qm = (c - ell) Pi(n-|m)/g-  and  Qp = (c + ell) Pi(n+|m)/g+
        Qm = ( c - ell ) * ( elK / g_m + 2 * Aq * rj( 0., kc2, 1., max( t_m / g_m, 1e-90 ) ) / ( 3 * g_m**2 ) )
        Qp = ( c + ell ) * ( elK * ( 2 * Aq * kc2 + m * c_p ) / ( c_p * ( m * g_p + 2 * Aq ) ) &
                             + 2 * Aq * kc2 * rj( 0., kc2, 1., max( kc2 * g_p / c_p, 1e-90 ) ) / ( 3 * c_p**2 ) )

        h_m = Aq**2 - Bq * ( Bq - c )
        h_p = Aq**2 - Bq * ( Bq + c )
        k_m = ( c - Bq )**2 * ( c + Bq ) - Aq**2 * Bq
        k_p = ( c + Bq )**2 * ( c - Bq ) + Aq**2 * Bq

        calA = pi / Aq * ( ( c - Bq ) / e_m - ( c + Bq ) / e_p ) &
             - 2 / ( Aq * r_far ) * ( ( c - Bq ) * Qm - ( c + Bq ) * Qp + 2 * ell * elK )

        calD = pi * c * ( 1 / e_m - 1 / e_p ) - 2 * c / r_far * ( Qm - Qp ) - 8 * ell * elK / r_far

        calB = pi / Aq * ( h_m / e_m + h_p / e_p ) - 2 / ( Aq * r_far ) * ( h_m * Qm + h_p * Qp ) &
             + 4 * ell / r_far * P1 - 4 * Bq * ell / ( Aq * r_far ) * elK

        calC = pi / Aq**2 * ( k_m / e_m - k_p / e_p ) - 2 / ( Aq**2 * r_far ) * ( k_m * Qm - k_p * Qp ) &
             + 4 / ( Aq**2 * r_far ) * ( -Aq * c**2 - 2 * Aq**2 * ell + Aq * Bq * ell + Bq**2 * ell - c**2 * ell ) * elK &
             + 16 * ell / ( m * r_far ) * P1 + 8 / ( 3 * Aq * r_far ) * ( c**2 + 2 * Aq * ell - Bq * ell ) * elD

        !::I0 = int_0^{2 pi} ln( (c + Delta_ell + ell)/(c + Delta_ell - ell) ) dpsi, an even function of psi
        dat%x1 = ell
        dat%x2 = Aq
        dat%y1 = Bq
        dat%y2 = t_m
        dat%z1 = g_p
        f_ptr => cone_I0_integrand
        call qags( f_ptr, dat, 0., pi, cone_I0_epsabs, cone_I0_epsrel, I0, abserr, neval, ier )
        I0 = 2 * I0

        Pq = rho * calA - st * calC
        Qq = z * calA - ct * calB
        Rq = rho * calA - st * ( calD + I0 )

        c2p = cphi**2 - sphi**2
        s2p = 2 * sphi * cphi
        fac = st * ct / ( 8 * pi )

        N(1,1) = fac * ( Rq + Pq * c2p )
        N(2,2) = fac * ( Rq - Pq * c2p )
        N(1,2) = fac * Pq * s2p
        N(2,1) = N(1,2)
        N(1,3) = 2 * fac * Qq * cphi
        N(3,1) = N(1,3)
        N(2,3) = 2 * fac * Qq * sphi
        N(3,2) = N(2,3)
        N(3,3) = -2 * fac * Rq
        !::Inside the cone H = B/mu0 - M
        if ( rho * ct - z * st .lt. 0. .AND. z .lt. h ) N(3,3) = N(3,3) - 1

    end subroutine cone_closed

    !::P1 = K(m) - 2/3 rd(0,1-m,1) = ( 2 E(m) - (2 - m) K(m) ) / m. Both terms tend to pi/2 for m -> 0, so for
    !::small m the function is summed from its power series,
    !::  P1 = -pi/2 sum_{j>=1} a_j j/(j+1) m**j,   a_j = ( (2j-1)!!/(2j)!! )**2.
    function cone_P1( m, elK, elD )
    real :: cone_P1
    real,intent(in) :: m, elK, elD
    real :: a, term, mj
    integer :: j

        if ( m .ge. cone_m_series ) then
            cone_P1 = elK - 2 * elD / 3
        else
            a = 1.
            mj = 1.
            cone_P1 = 0.
            do j=1,30
                a = a * ( real( 2 * j - 1 ) / real( 2 * j ) )**2
                mj = mj * m
                term = a * real( j ) / real( j + 1 ) * mj
                cone_P1 = cone_P1 + term
                if ( term .lt. 1e-17 * cone_P1 ) exit
            enddo
            cone_P1 = -pi / 2 * cone_P1
        endif

    end function cone_P1

    !::The integrand of I0, the potential of the generator at the azimuth psi relative to the field point,
    !::ln( (ell - b + Delta)/(c - b) ) with b = B + A cos(psi) the projection of the field point on the
    !::generator and Delta its distance to the rim point of the generator. For b > ell the equal form
    !::ln( (c + b)/(Delta + b - ell) ) is used. c - b and c + b are formed from t_m = c - B - A and
    !::g_p = c + B - A, so that they keep their accuracy where they vanish.
    function cone_I0_integrand( psi, dat )
    real :: cone_I0_integrand
    real,intent(in) :: psi
    class(dataCollectionBase), intent(inout), target :: dat
    real :: ell, Aq, Bq, c_minus_b, c_plus_b, b, Delta

        ell = dat%x1
        Aq = dat%x2
        Bq = dat%y1
        c_minus_b = dat%y2 + 2 * Aq * sin( psi / 2 )**2
        c_plus_b = dat%z1 + 2 * Aq * cos( psi / 2 )**2
        b = Bq + Aq * cos( psi )
        Delta = sqrt( ( ell - b )**2 + c_minus_b * c_plus_b )

        if ( ell - b .ge. 0. ) then
            cone_I0_integrand = log( ( ell - b + Delta ) / max( c_minus_b, tiny( 1. ) ) )
        else
            cone_I0_integrand = log( c_plus_b / max( Delta + b - ell, tiny( 1. ) ) )
        endif

    end function cone_I0_integrand

    !::The tensor close to the axis, in the frame of the paper, from the expansion of the potential about the
    !::axis. With S(z) twice the axial flux density on the axis for a unit axial magnetization, the on-axis
    !::field of a unit magnetization is H_z = S/2, minus one inside the cone, for the axial and H_rho = -S/4 for
    !::the transverse direction, by the trace of the tensor. Off the axis the potentials are
    !::Phi = Phi0 - rho**2 Phi0''/4 + rho**4 Phi0''''/64 (axial) and
    !::Phi = cos(phi) (rho g - rho**3 g''/8 + rho**5 g''''/192) (transverse, g = S/4), as for the cylinder
    !::(getN_cylinder_axis). The error is of the sixth order in the distance to the axis relative to the
    !::radius of convergence. Valid on the axis itself, except at the apex.
    !::
    !::From the on-axis expressions of the paper, with u = h - z, q = sqrt(Rb**2 + u**2) the distance to the
    !::rim, w = ell - z ct + q and sg = sign(z),
    !::  S = sg ct**2 + ( ct**2 u - h st**2 )/q + st**2 ct ( ln(w) - ln|z| - ln(1 - sg ct) ).
    !::T1 is the second and T2 = ln(w) - ln|z| the third term; T1a..T1d and T2a..T2d are their first to fourth
    !::derivatives with respect to u, and the derivatives with respect to z are (-1)**n times these.
    subroutine cone_axis( Rb, h, ell, st, ct, x, y, z, N )
    real,intent(in) :: Rb, h, ell, st, ct, x, y, z
    real,dimension(3,3),intent(inout) :: N
    real :: u, q, w, sg, a, b, kap, R2, rho2
    real :: T1, T1a, T1b, T1c, T1d, w1, w2, w3, w4, L1, L2, L3, L4, T2, T2a, T2b, T2c, T2d
    real :: S, S1, S2, S3, S4

        u = h - z
        R2 = Rb**2
        q = sqrt( R2 + u**2 )
        w = ell * st**2 + u * ct + q
        sg = sign( 1., z )
        a = ct**2
        b = -h * st**2
        kap = st**2 * ct
        rho2 = x**2 + y**2

        T1  = ( a * u + b ) / q
        T1a = ( a * R2 - b * u ) / q**3
        T1b = ( 2 * b * u**2 - 3 * a * R2 * u - b * R2 ) / q**5
        T1c = ( -6 * b * u**3 + 12 * a * R2 * u**2 + 9 * b * R2 * u - 3 * a * R2**2 ) / q**7
        T1d = ( 24 * b * u**4 - 60 * a * R2 * u**3 - 72 * b * R2 * u**2 + 45 * a * R2**2 * u + 9 * b * R2**2 ) / q**9

        !::Derivatives of w and of ln(w) with respect to u
        w1 = ct + u / q
        w2 = R2 / q**3
        w3 = -3 * R2 * u / q**5
        w4 = 3 * R2 * ( 4 * u**2 - R2 ) / q**7
        L1 = w1 / w
        L2 = w2 / w - L1**2
        L3 = w3 / w - 3 * L1 * w2 / w + 2 * L1**3
        L4 = w4 / w - 4 * L1 * w3 / w - 3 * ( w2 / w )**2 + 12 * L1**2 * w2 / w - 6 * L1**4

        T2  = log( w ) - log( abs( z ) )
        T2a = L1 + 1 / z
        T2b = L2 + 1 / z**2
        T2c = L3 + 2 / z**3
        T2d = L4 + 6 / z**4

        !::S and its first four derivatives with respect to z
        S  = sg * a + T1 + kap * ( T2 - log( 1 - sg * ct ) )
        S1 = -( T1a + kap * T2a )
        S2 = T1b + kap * T2b
        S3 = -( T1c + kap * T2c )
        S4 = T1d + kap * T2d

        N(1,1) = -S / 4 + S2 * ( 3 * x**2 + y**2 ) / 32 - S4 * rho2 * ( 5 * x**2 + y**2 ) / 768
        N(2,2) = -S / 4 + S2 * ( 3 * y**2 + x**2 ) / 32 - S4 * rho2 * ( 5 * y**2 + x**2 ) / 768
        N(1,2) = x * y * ( S2 / 16 - rho2 * S4 / 192 )
        N(2,1) = N(1,2)
        N(1,3) = -x * ( S1 / 4 - rho2 * S3 / 32 )
        N(3,1) = N(1,3)
        N(2,3) = -y * ( S1 / 4 - rho2 * S3 / 32 )
        N(3,2) = N(2,3)
        N(3,3) = S / 2 - rho2 * S2 / 8 + rho2**2 * S4 / 128
        !::Inside the cone H = B/mu0 - M. The points of this branch lie well inside or well outside the cone.
        if ( z .gt. 0. .AND. z .lt. h .AND. sqrt( rho2 ) * ct - z * st .lt. 0. ) N(3,3) = N(3,3) - 1

    end subroutine cone_axis

end module TileConeTensor
