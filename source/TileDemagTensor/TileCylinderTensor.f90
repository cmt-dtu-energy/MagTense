module TileCylinderTensor
    !::Demagnetization tensor of a full cylinder of radius R and half-length Lh, centred at the origin with its
    !::axis along z, from the closed-form field of Caciagli, Baars, Philipse and Kuipers, "Exact expression
    !::for the magnetic field of a finite cylinder with arbitrary uniform magnetization", J. Magn. Magn.
    !::Mater. 456 (2018) 423, Eqs. (3)-(5) for the axial and Eqs. (20)-(22) for the transverse magnetization.
    !::The tensor follows the MagTense convention H = N.M: its columns are the field of a unit magnetization
    !::along x, y and z, and inside the cylinder the axial column is H = B/mu0 - M.
    !::
    !::The complete elliptic integrals of the paper are evaluated through the Carlson symmetric forms of
    !::SpecialFunctions, which also removes the cancellations in the paper's auxiliary functions:
    !::  K(m) = rf(0,k2,1),  K(m) - E(m) = m/3 rd(0,k2,1),  Pi(n|m) - K(m) = n/3 rj(0,k2,1,gamma2)
    !::with m = 1 - k2 and n = 1 - gamma2 (the paper's K, E, P are K(1-k2), E(1-k2), Pi(1-gamma2|1-k2)), so
    !::  P1 = K - 2/3 rd
    !::  P2 = K - gamma (1 + gamma) rj / 3
    !::  P3 = ( rd - gamma2 rj ) / 3
    !::  P4 = gamma (1 + gamma2) rj / 3 - (1 + gamma) K + 2/3 rd
    use SpecialFunctions
    use TileTensorHelperFunctions
    implicit none

    !::Relative distance from the axis, rho/R, below which the tensor is evaluated from its expansion about
    !::the axis (getN_cylinder_axis) instead of the closed form. The closed form contains 1/rho and loses
    !::about two digits per decade below rho ~ 1e-2 R (error ~1e-12 at 1e-2 R, ~1e-9 at 1e-4 R, up to a
    !::factor 50 worse for a cylinder much longer than its radius), while the expansion, fourth order in rho,
    !::has an error of about (rho/R)**6 (~5e-12 at 1e-2 R, ~1e-13 at 5e-3 R).
    real, parameter :: cyl_axis_tol = 1e-2

    !::Floor on |gamma| = |rho - R| / (rho + R). On the lateral surface the elliptic integral of the third
    !::kind diverges and the transverse field jumps by (M.n) n, so points this close to the surface take the
    !::limit from outside, which is also the side the inside test (strict inequalities) assigns them to.
    real, parameter :: cyl_gamma_min = 1e-12

    !::Floor on k2. On the edge circle rho = R, |z| = Lh the integral K(1-k2) diverges logarithmically; the
    !::floor gives a large but finite tensor there. rd requires its arguments to sum to at least 1e-25.
    real, parameter :: cyl_ksq_min = 1e-24

    contains

    !::The demagnetization tensor N (H = N.M) at the point (x,y,z) of a cylinder of radius R and half-length
    !::Lh, centred at the origin with its axis along z. Points on the surface count as outside.
    subroutine getN_cylinder( R, Lh, x, y, z, N )
    real,intent(in) :: R, Lh, x, y, z
    real,dimension(3,3),intent(inout) :: N
    real :: rho, c, s, gam, A1, A2, A3, A4, fac
    real :: alp, bep, P1p, P2p, P3p, P4p, alm, bem, P1m, P2m, P3m, P4m

        rho = sqrt( x**2 + y**2 )

        if ( rho .lt. cyl_axis_tol * R ) then
            call getN_cylinder_axis( R, Lh, x, y, z, N )
            return
        endif

        c = x / rho
        s = y / rho
        gam = ( rho - R ) / ( rho + R )
        if ( abs(gam) .lt. cyl_gamma_min ) gam = cyl_gamma_min  !< outside limit, see cyl_gamma_min

        call cylinder_P_functions( z + Lh, rho, R, gam, alp, bep, P1p, P2p, P3p, P4p )
        call cylinder_P_functions( z - Lh, rho, R, gam, alm, bem, P1m, P2m, P3m, P4m )

        !::The differences between the two ends that the field expressions share, Eqs. (3) and (21)
        A1 = alp * P1p - alm * P1m
        A2 = bep * P2p - bem * P2m
        A3 = bep * P3p - bem * P3m
        A4 = bep * P4p - bem * P4m

        !::Transverse magnetization along x, Eq. (21): H_rho = R cos(phi)/(2 pi rho) A4, H_phi = R sin(phi)/(pi rho) A3,
        !::H_z = R cos(phi)/pi A1, converted to Cartesian components; the y column is the same with phi - pi/2
        fac = R / ( 2 * pi * rho )
        N(1,1) = fac * ( A4 * c**2 - 2 * A3 * s**2 )
        N(2,2) = fac * ( A4 * s**2 - 2 * A3 * c**2 )
        N(1,2) = fac * ( A4 + 2 * A3 ) * s * c
        N(2,1) = N(1,2)
        N(3,1) = R / pi * A1 * c
        N(3,2) = R / pi * A1 * s
        !::Axial magnetization, Eq. (3): B_rho/mu0 = R/pi A1, which equals the z-row above, and B_z/mu0 = R/(pi (rho+R)) A2
        N(1,3) = N(3,1)
        N(2,3) = N(3,2)
        N(3,3) = R / ( pi * ( rho + R ) ) * A2
        !::Inside the cylinder H = B/mu0 - M. The test agrees with the floor on gamma above: a point within the
        !::floor of the lateral surface takes the outside limit of B_z, which jumps by M_z across the surface
        if ( gam .lt. -cyl_gamma_min .AND. abs(z) .lt. Lh ) N(3,3) = N(3,3) - 1

    end subroutine getN_cylinder

    !::The auxiliary functions of the paper for one end of the cylinder, xi = z + Lh or z - Lh: alpha and beta
    !::of Eq. (5) and P1 to P4 of Eqs. (4), (20) and (22), in the Carlson forms given at the top of the module.
    subroutine cylinder_P_functions( xi, rho, R, gam, alpha, beta, P1, P2, P3, P4 )
    real,intent(in) :: xi, rho, R, gam
    real,intent(out) :: alpha, beta, P1, P2, P3, P4
    real :: ksq, gamsq, elK, elD, elJ

        ksq = ( xi**2 + ( rho - R )**2 ) / ( xi**2 + ( rho + R )**2 )
        if ( ksq .lt. cyl_ksq_min ) ksq = cyl_ksq_min
        gamsq = gam**2

        elK = rf( 0., ksq, 1. )
        elD = rd( 0., ksq, 1. )
        elJ = rj( 0., ksq, 1., gamsq )

        P1 = elK - 2 * elD / 3
        P2 = elK - gam * ( 1 + gam ) * elJ / 3
        P3 = ( elD - gamsq * elJ ) / 3
        P4 = gam * ( 1 + gamsq ) * elJ / 3 - ( 1 + gam ) * elK + 2 * elD / 3

        alpha = 1 / sqrt( xi**2 + ( rho + R )**2 )
        beta = xi * alpha

    end subroutine cylinder_P_functions

    !::The tensor close to the axis from the expansion of the potential about the axis. With
    !::S(z) = xi+/sqrt(xi+**2 + R**2) - xi-/sqrt(xi-**2 + R**2) the on-axis field of a unit magnetization is
    !::H_z = S/2 for the axial (the solenoid, Eq. (3) at rho = 0) and H_rho = -S/4 for the transverse
    !::direction (Eq. (26)). Off the axis the potentials are Phi = Phi0 - rho**2 Phi0''/4 + rho**4 Phi0''''/64
    !::(axial) and Phi = cos(phi) (rho g - rho**3 g''/8 + rho**5 g''''/192) (transverse, g = S/4), which are
    !::the solutions of Laplace's equation with these axis values. The error is of order (rho/R)**6 in the
    !::diagonal and (rho/R)**5 in the off-diagonal entries. Valid on the axis itself.
    subroutine getN_cylinder_axis( R, Lh, x, y, z, N )
    real,intent(in) :: R, Lh, x, y, z
    real,dimension(3,3),intent(inout) :: N
    real :: xp, xm, dp, dm, R2, S, S1, S2, S3, S4, rho2

        xp = z + Lh
        xm = z - Lh
        dp = sqrt( xp**2 + R**2 )
        dm = sqrt( xm**2 + R**2 )
        R2 = R**2
        rho2 = x**2 + y**2

        !::S and its first four derivatives with respect to z
        S  = xp / dp - xm / dm
        S1 = R2 * ( 1 / dp**3 - 1 / dm**3 )
        S2 = -3 * R2 * ( xp / dp**5 - xm / dm**5 )
        S3 = -3 * R2 * ( ( R2 - 4 * xp**2 ) / dp**7 - ( R2 - 4 * xm**2 ) / dm**7 )
        S4 = 15 * R2 * ( xp * ( 3 * R2 - 4 * xp**2 ) / dp**9 - xm * ( 3 * R2 - 4 * xm**2 ) / dm**9 )

        N(1,1) = -S / 4 + S2 * ( 3 * x**2 + y**2 ) / 32 - S4 * rho2 * ( 5 * x**2 + y**2 ) / 768
        N(2,2) = -S / 4 + S2 * ( 3 * y**2 + x**2 ) / 32 - S4 * rho2 * ( 5 * y**2 + x**2 ) / 768
        N(1,2) = x * y * ( S2 / 16 - rho2 * S4 / 192 )
        N(2,1) = N(1,2)
        N(1,3) = -x * ( S1 / 4 - rho2 * S3 / 32 )
        N(3,1) = N(1,3)
        N(2,3) = -y * ( S1 / 4 - rho2 * S3 / 32 )
        N(3,2) = N(2,3)
        N(3,3) = S / 2 - rho2 * S2 / 8 + rho2**2 * S4 / 128
        if ( sqrt(rho2) .lt. R .AND. abs(z) .lt. Lh ) N(3,3) = N(3,3) - 1

    end subroutine getN_cylinder_axis

end module TileCylinderTensor
