!---------------------------------------------------------------------------
!> @brief
!> Constant-density Galerkin Laplace integral over a pair of flat triangles,
!>
!>     I = \int_{T_a} \int_{T_b} 1 / |r_a - r_b| dS_b dS_a .
!>
!> The integral is evaluated in closed form by a Fortran port of the
!> constant-density branch of the dimensional-reduction hierarchy
!> I3 -> I2s/I2t -> I1 -> I0 of
!>
!>   N. A. Gumerov, S. Kaneko and R. Duraiswami,
!>   "Recursive computation of the multipole expansions of layer potential
!>    integrals over simplices for efficient fast multipole accelerated
!>    boundary elements",
!>   SIAM J. Sci. Comput. 46 (2024), DOI 10.1137/23M1547688,
!>
!> following the same reduction used by the (MIT licensed) reference
!> implementation GalerkinLaplaceTriGS.m and by its C++ port in the dip-fmm
!> project (src/geometry/primitives/tetrahedron.cpp).
!>
!> ACCURACY. The reduction is exact in exact arithmetic, but its terms carry
!> coefficients that grow without bound when the geometry is nearly degenerate
!> (nearly parallel planes, nearly collinear or coplanar edge vectors), and it
!> drops or snaps the corresponding directions with fixed tolerances. In double
!> precision this loses all digits for parallel faces with non-aligned edges,
!> e.g. the top and bottom faces of a tetrahedral film mesh (validation of
!> 19 Oct 2026 against an independent Sauter-Schwab quadrature, FMM-article
!> repository, tet-tet-article/data/REPORT.md "A1"). Therefore:
!>   - a geometric screen leaves well-conditioned face pairs to the former
!>     double-precision reduction (module TileTriangleTriangleDouble); all
!>     others are reduced in quadruple precision (TileTriangleTriangleQuad)
!>     with tolerances chosen for it. Both modules include the same reduction
!>     source, TileTriangleTriangleReduction.inc;
!>   - a Gram-Schmidt rank decision close to the rank tolerance, where neither
!>     dropping nor keeping the direction is reliable, flags the face pair,
!>     which is then integrated numerically instead (Sauter-Schwab quadrature
!>     for touching triangles, Gauss quadrature otherwise);
!>   - triangles without a common vertex whose planes are nearly, but not
!>     exactly, parallel are integrated by Gauss quadrature directly: there the
!>     reduction expands the offset in a nearly degenerate basis whatever the
!>     precision.
!> The quadratures are in double precision; their integrand is positive, so
!> they do not cancel.
!>
!> The interface is double precision (kind ttReal) so that the result does not
!> depend on whether the build promotes the default real kind (ifx -real-size
!> 64, gfortran -fdefault-real-8).
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!> @brief The triangle-pair reduction in quadruple precision, for face pairs
!> the geometric screen of triangle_triangle_laplace_integral does not clear.
!> Tolerances chosen for this precision (FMM-article tet-tet-tests, REPORT
!> "A1"): coefficient drop 1e-22, height snap 1e-12, I0 zero test 1e-24, rank
!> 1e-23 (squared, relative); a rank decision with a squared relative norm in
!> [1e-26, 1e-20) sets the guard; face pairs within 1e-8 of parallel are
!> projected onto parallel planes.
!---------------------------------------------------------------------------
module TileTriangleTriangleQuad
    use, intrinsic :: ieee_arithmetic
    implicit none
    private
    integer,parameter :: dp = selected_real_kind(15, 307)
    integer,parameter :: wk = selected_real_kind(33, 4931)
    integer,parameter :: wi = wk
    real(wk),parameter :: zeroTol = 1.0e-22_wk
    real(wi),parameter :: i0Z2 = 1.0e-24_wi
    real(wk),parameter :: relResidual = 1.0e-12_wk
    real(wk),parameter :: rankRel2 = 1.0e-23_wk
    real(wk),parameter :: guardLow = 1.0e-26_wk
    real(wk),parameter :: guardHigh = 1.0e-20_wk
    real(wk),parameter :: nearParallel = 1.0e-8_wk
    public :: reduce_triangle_pair
    contains
    include 'TileTriangleTriangleReduction.inc'
end module TileTriangleTriangleQuad


!---------------------------------------------------------------------------
!> @brief The triangle-pair reduction in double precision with the former
!> double-precision tolerances (I0 in quadruple precision), for the
!> well-conditioned face pairs the geometric screen clears; the guard is off.
!---------------------------------------------------------------------------
module TileTriangleTriangleDouble
    use, intrinsic :: ieee_arithmetic
    implicit none
    private
    integer,parameter :: dp = selected_real_kind(15, 307)
    integer,parameter :: wk = dp
    integer,parameter :: wi = selected_real_kind(33, 4931)
    real(wk),parameter :: zeroTol = 1.0e-14_wk
    real(wi),parameter :: i0Z2 = 1.0e-28_wi
    real(wk),parameter :: relResidual = 1.0e-10_wk
    real(wk),parameter :: rankRel2 = 256._wk * epsilon( 1._wk )
    real(wk),parameter :: guardLow = 0._wk
    real(wk),parameter :: guardHigh = 0._wk
    real(wk),parameter :: nearParallel = 1.0e-3_wk
    public :: reduce_triangle_pair
    contains
    include 'TileTriangleTriangleReduction.inc'
end module TileTriangleTriangleDouble


module TileTriangleTriangle
    use, intrinsic :: ieee_arithmetic
    use TileTriangleTriangleQuad, only : reduce_quad => reduce_triangle_pair
    use TileTriangleTriangleDouble, only : reduce_double => reduce_triangle_pair
    implicit none

    private

    !> Precision of the interface and of the quadratures.
    integer,parameter,public :: ttReal = selected_real_kind(15, 307)

    !> The geometric screen: a face pair whose smallest non-exact measure of
    !> near-degeneracy (angle between the normals, between an edge and the other
    !> plane, between two edges, height of a vertex over the other plane;
    !> measures below screenExact are exact degeneracies) is at least screenBand
    !> is reduced in double precision (TileTriangleTriangleDouble), all others in
    !> quadruple precision with the guard (TileTriangleTriangleQuad). On the
    !> 440-pair validation set the double-precision reduction is accurate to
    !> 3.6e-11 on every face pair it clears (26 % are not cleared); a band of
    !> 1e-3 would have let 2.8e-9 through.
    real(ttReal),parameter :: screenBand = 1.0e-2_ttReal
    real(ttReal),parameter :: screenExact = 1.0e-13_ttReal

    !> Separated triangles (see separatedVertexTol) whose normals differ by
    !> more than round-off but less than this are integrated by Gauss
    !> quadrature.
    real(ttReal),parameter :: parallelQuadrature = 1.0e-2_ttReal

    !> Vertices closer than this fraction of the larger edge are common. It is
    !> loose on purpose: rounded mesh coordinates (dip-fmm's canonical grid
    !> moves shared vertices by up to ~1e-9 of the root) leave common vertices
    !> slightly apart, and the Sauter-Schwab maps tolerate such a mismatch.
    real(ttReal),parameter :: sharedVertexTol = 1.0e-6_ttReal

    !> Triangles count as separated, and may be integrated by Gauss
    !> quadrature, only when every vertex pair is further apart than this
    !> fraction of the larger edge: a nearly touching pair needs more
    !> splitting than maxSplitDepth and stays with the closed form, whose
    !> near-parallel projection and height snapping exist for that case.
    real(ttReal),parameter :: separatedVertexTol = 1.0e-3_ttReal

    !> Gauss quadrature of separated triangles: q x q collapsed Gauss points per
    !> triangle, the larger triangle split into four until the gap between the
    !> bounding spheres is at least separationRatio times the larger radius.
    integer,parameter :: regularOrder = 10
    real(ttReal),parameter :: separationRatio = 1.0_ttReal
    integer,parameter :: maxSplitDepth = 12

    !> Sauter-Schwab quadrature of touching triangles: q**4 Gauss points per
    !> subdomain (relative error 1e-9 to 1e-12 on mesh face pairs).
    integer,parameter :: ssOrder = 24

    !> Evaluation path, returned through the optional argument of
    !> triangle_triangle_laplace_integral.
    integer,parameter,public :: ttClosedForm = 0
    integer,parameter,public :: ttGaussParallel = 1
    integer,parameter,public :: ttGuardQuadrature = 2

    public :: triangle_triangle_laplace_integral
    public :: gauss_legendre01

    !> Measurement only (branch tet-pair-stats): face pairs per evaluation path,
    !> 1 separated near-parallel Gauss, 2 FP64 (screened) closed form, 3 quadruple
    !> closed form, 4 guard -> Sauter-Schwab, 5 guard -> split Gauss, 6 guard kept.
    integer(8),public :: ttPathCount(6) = 0
    !> Counting is off unless tet_pair_stats_reset finds MAGTENSE_TETPAIR_STATS set,
    !> so production runs pay neither the atomics nor their contention.
    logical,public :: ttStatsOn = .false.

    contains

    !---------------------------------------------------------------------------
    !> @brief Galerkin Laplace integral over a pair of triangles.
    !> @param[in] triA first triangle, triA(:,k) is the k'th vertex
    !> @param[in] triB second triangle, triB(:,k) is the k'th vertex
    !> @param[out] val the value of the double surface integral of 1/R
    !> @param[out] ok .false. if either triangle is degenerate or the reduction
    !>             produced a non-finite value; val is then not usable
    !> @param[out] method optional, the evaluation path (ttClosedForm, ...)
    !>
    !> The result is symmetric in its two arguments and invariant under a
    !> cyclic or reflecting permutation of the vertices of either triangle.
    !---------------------------------------------------------------------------
    subroutine triangle_triangle_laplace_integral( triA, triB, val, ok, method )
    real(ttReal),dimension(3,3),intent(in) :: triA, triB
    real(ttReal),intent(out) :: val
    logical,intent(out) :: ok
    integer,intent(out),optional :: method

    real(ttReal),dimension(3,3) :: first, second
    real(ttReal) :: scaleA, scaleB, valR
    real(ttReal),dimension(3) :: nA, nB
    real(ttReal) :: pm, edge, closest
    integer :: i, j, nShared
    logical :: guard, separated

        val = 0._ttReal
        ok = .true.
        if ( present(method) ) method = ttClosedForm

        if ( .not. all( ieee_is_finite( triA ) ) .or. .not. all( ieee_is_finite( triB ) ) ) then
            ok = .false.
            return
        endif

        nA = cross_dp( triA(:,2) - triA(:,1), triA(:,3) - triA(:,1) )
        nB = cross_dp( triB(:,2) - triB(:,1), triB(:,3) - triB(:,1) )
        if ( .not. ( norm2( nA ) .gt. 0._ttReal ) .or. .not. ( norm2( nB ) .gt. 0._ttReal ) ) then
            ok = .false.
            return
        endif
        nA = nA / norm2( nA )
        nB = nB / norm2( nB )
        pm = min( norm2( nA + nB ), norm2( nA - nB ) )

        edge = max( triangle_scale_dp( triA ), triangle_scale_dp( triB ) )
        nShared = 0
        closest = huge( 1._ttReal )
        do i = 1, 3
            do j = 1, 3
                closest = min( closest, norm2( triA(:,i) - triB(:,j) ) )
                if ( norm2( triA(:,i) - triB(:,j) ) .le. sharedVertexTol * edge ) nShared = nShared + 1
            enddo
        enddo
        separated = closest .gt. separatedVertexTol * edge

        !Separated, nearly parallel triangles: Gauss quadrature (see the header).
        if ( separated .and. pm .gt. 1.0e-14_ttReal .and. pm .lt. parallelQuadrature ) then
            val = triangle_triangle_regular( triA, triB, 0 )
            if ( present(method) ) method = ttGaussParallel
            if ( ttStatsOn ) then
                !$omp atomic
                ttPathCount(1) = ttPathCount(1) + 1
            endif
            if ( .not. ieee_is_finite( val ) ) ok = .false.
            return
        endif

        !The reduction is symmetric in its two triangles. Put the larger one
        !first, which is the ordering the reference reduction is written for.
        scaleA = triangle_scale_dp( triA )
        scaleB = triangle_scale_dp( triB )
        if ( scaleA .lt. scaleB ) then
            first = triB
            second = triA
        else
            first = triA
            second = triB
        endif

        if ( smallest_degeneracy( triA, triB, nA, nB, pm, edge ) .ge. screenBand ) then
            !Well-conditioned: the double-precision reduction.
            call reduce_double( first, second, edge, valR, ok, guard )
            val = valR
            if ( ttStatsOn ) then
                !$omp atomic
                ttPathCount(2) = ttPathCount(2) + 1
            endif
            if ( .not. ieee_is_finite( val ) ) ok = .false.
            return
        endif

        call reduce_quad( first, second, edge, valR, ok, guard )

        if ( guard .and. nShared .gt. 0 ) then
            !An unreliable rank decision: integrate numerically instead.
            ok = .true.
            val = triangle_triangle_sauter_schwab( triA, triB, edge )
            if ( present(method) ) method = ttGuardQuadrature
            if ( ttStatsOn ) then
                !$omp atomic
                ttPathCount(4) = ttPathCount(4) + 1
            endif
        elseif ( guard .and. separated ) then
            ok = .true.
            val = triangle_triangle_regular( triA, triB, 0 )
            if ( present(method) ) method = ttGuardQuadrature
            if ( ttStatsOn ) then
                !$omp atomic
                ttPathCount(5) = ttPathCount(5) + 1
            endif
        else
            !A guarded pair that is neither touching nor separated (nearly
            !touching without a common vertex) keeps the closed form: no
            !fallback here resolves its gap.
            val = valR
            if ( guard ) then
                if ( ttStatsOn ) then
                    !$omp atomic
                    ttPathCount(6) = ttPathCount(6) + 1
                endif
            else
                if ( ttStatsOn ) then
                    !$omp atomic
                    ttPathCount(3) = ttPathCount(3) + 1
                endif
            endif
        endif

        if ( .not. ieee_is_finite( val ) ) ok = .false.

    end subroutine triangle_triangle_laplace_integral


    !---------------------------------------------------------------------------
    !> @brief The smallest measure of near-degeneracy of a triangle pair that is
    !> not an exact degeneracy (see screenBand); huge() if all are exact.
    !---------------------------------------------------------------------------
    function smallest_degeneracy( triA, triB, nA, nB, pm, edge ) result( smallest )
    real(ttReal),dimension(3,3),intent(in) :: triA, triB
    real(ttReal),dimension(3),intent(in) :: nA, nB
    real(ttReal),intent(in) :: pm, edge
    real(ttReal) :: smallest
    real(ttReal),dimension(3,3) :: eA, eB
    integer :: i, j

        smallest = huge( 1._ttReal )
        call consider( pm )
        do i = 1, 3
            eA(:,i) = triA(:,modulo(i,3)+1) - triA(:,i)
            eA(:,i) = eA(:,i) / norm2( eA(:,i) )
            eB(:,i) = triB(:,modulo(i,3)+1) - triB(:,i)
            eB(:,i) = eB(:,i) / norm2( eB(:,i) )
        enddo
        do i = 1, 3
            !An edge nearly parallel to the other triangle's plane.
            call consider( abs( dot_product( nA, eB(:,i) ) ) )
            call consider( abs( dot_product( nB, eA(:,i) ) ) )
            !A vertex nearly in the other triangle's plane.
            call consider( abs( dot_product( nA, triB(:,i) - triA(:,1) ) ) / edge )
            call consider( abs( dot_product( nB, triA(:,i) - triB(:,1) ) ) / edge )
            do j = 1, 3
                !Two nearly parallel edges.
                call consider( norm2( cross_dp( eA(:,i), eB(:,j) ) ) )
            enddo
        enddo

    contains

        subroutine consider( measure )
        real(ttReal),intent(in) :: measure
            if ( measure .ge. screenExact ) smallest = min( smallest, measure )
        end subroutine consider

    end function smallest_degeneracy


    function triangle_scale_dp( tri ) result( s )
    real(ttReal),dimension(3,3),intent(in) :: tri
    real(ttReal) :: s
        s = max( norm2( tri(:,2) - tri(:,1) ), &
                 norm2( tri(:,3) - tri(:,2) ), &
                 norm2( tri(:,1) - tri(:,3) ) )
    end function triangle_scale_dp




    function cross_dp( a, b ) result( c )
    real(ttReal),dimension(3),intent(in) :: a, b
    real(ttReal),dimension(3) :: c
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_dp


















    !> Gauss-Legendre nodes and weights on [0,1] (Newton iteration on the
    !> Legendre recurrence; the weights sum to one).
    subroutine gauss_legendre01( q, x, w )
    integer,intent(in) :: q
    real(ttReal),dimension(q),intent(out) :: x, w
    real(ttReal),parameter :: pi = 3.14159265358979323846264338327950288_ttReal
    real(ttReal) :: z, z1, p1, p2, p3, pp
    integer :: i, j, it
        do i = 1, q
            z = cos( pi * ( real( i, ttReal ) - 0.25_ttReal ) / ( real( q, ttReal ) + 0.5_ttReal ) )
            do it = 1, 100
                p1 = 1._ttReal
                p2 = 0._ttReal
                do j = 1, q
                    p3 = p2
                    p2 = p1
                    p1 = ( real( 2*j-1, ttReal ) * z * p2 - real( j-1, ttReal ) * p3 ) / real( j, ttReal )
                enddo
                pp = real( q, ttReal ) * ( z * p1 - p2 ) / ( z*z - 1._ttReal )
                z1 = z
                z = z1 - p1 / pp
                if ( abs( z - z1 ) .le. 4._ttReal * epsilon( z ) ) exit
            enddo
            x(i) = 0.5_ttReal * ( 1._ttReal - z )
            w(i) = 1._ttReal / ( ( 1._ttReal - z*z ) * pp*pp )
        enddo
    end subroutine gauss_legendre01


    !> Twice the area of a triangle.
    function twice_area( tri ) result( a2 )
    real(ttReal),dimension(3,3),intent(in) :: tri
    real(ttReal) :: a2
        a2 = norm2( cross_dp( tri(:,2) - tri(:,1), tri(:,3) - tri(:,2) ) )
    end function twice_area


    !---------------------------------------------------------------------------
    !> @brief Gauss quadrature of the pair integral for triangles without a
    !> common point: collapsed Gauss rules of regularOrder**2 points on both
    !> triangles, the larger triangle split into four (recursively) until the
    !> two are well separated (see separationRatio).
    !---------------------------------------------------------------------------
    recursive function triangle_triangle_regular( triA, triB, depth ) result( val )
    real(ttReal),dimension(3,3),intent(in) :: triA, triB
    integer,intent(in) :: depth
    real(ttReal) :: val

    real(ttReal),dimension(3) :: ca, cb
    real(ttReal),dimension(3,3) :: sub
    real(ttReal),dimension(3,3,4) :: children
    real(ttReal),dimension(regularOrder) :: x, w
    real(ttReal),dimension(3,regularOrder*regularOrder) :: pa, pb
    real(ttReal),dimension(regularOrder*regularOrder) :: wt
    real(ttReal) :: ra, rb, gap, acc
    integer :: i, j, k, n

        ca = sum( triA, dim=2 ) / 3._ttReal
        cb = sum( triB, dim=2 ) / 3._ttReal
        ra = max( norm2( triA(:,1) - ca ), norm2( triA(:,2) - ca ), norm2( triA(:,3) - ca ) )
        rb = max( norm2( triB(:,1) - cb ), norm2( triB(:,2) - cb ), norm2( triB(:,3) - cb ) )
        gap = norm2( ca - cb ) - ra - rb

        if ( gap .lt. separationRatio * max( ra, rb ) .and. depth .lt. maxSplitDepth ) then
            val = 0._ttReal
            if ( ra .ge. rb ) then
                call split_triangle( triA, children )
                do k = 1, 4
                    sub = children(:,:,k)
                    val = val + triangle_triangle_regular( sub, triB, depth + 1 )
                enddo
            else
                call split_triangle( triB, children )
                do k = 1, 4
                    sub = children(:,:,k)
                    val = val + triangle_triangle_regular( triA, sub, depth + 1 )
                enddo
            endif
            return
        endif

        call gauss_legendre01( regularOrder, x, w )
        n = 0
        do i = 1, regularOrder
            do j = 1, regularOrder
                n = n + 1
                !chi(s, s t) = P1 + s (P2 - P1) + s t (P3 - P2), Jacobian s * 2 area
                pa(:,n) = triA(:,1) + x(i) * ( triA(:,2) - triA(:,1) ) + x(i) * x(j) * ( triA(:,3) - triA(:,2) )
                pb(:,n) = triB(:,1) + x(i) * ( triB(:,2) - triB(:,1) ) + x(i) * x(j) * ( triB(:,3) - triB(:,2) )
                wt(n) = w(i) * w(j) * x(i)
            enddo
        enddo
        acc = 0._ttReal
        do i = 1, n
            do j = 1, n
                acc = acc + wt(i) * wt(j) / norm2( pa(:,i) - pb(:,j) )
            enddo
        enddo
        val = twice_area( triA ) * twice_area( triB ) * acc

    end function triangle_triangle_regular


    !> The four midpoint children of a triangle.
    subroutine split_triangle( tri, children )
    real(ttReal),dimension(3,3),intent(in) :: tri
    real(ttReal),dimension(3,3,4),intent(out) :: children
    real(ttReal),dimension(3) :: m12, m23, m31
        m12 = 0.5_ttReal * ( tri(:,1) + tri(:,2) )
        m23 = 0.5_ttReal * ( tri(:,2) + tri(:,3) )
        m31 = 0.5_ttReal * ( tri(:,3) + tri(:,1) )
        children(:,1,1) = tri(:,1)
        children(:,2,1) = m12
        children(:,3,1) = m31
        children(:,1,2) = m12
        children(:,2,2) = tri(:,2)
        children(:,3,2) = m23
        children(:,1,3) = m31
        children(:,2,3) = m23
        children(:,3,3) = tri(:,3)
        children(:,1,4) = m12
        children(:,2,4) = m23
        children(:,3,4) = m31
    end subroutine split_triangle


    !---------------------------------------------------------------------------
    !> @brief Sauter-Schwab quadrature of the pair integral for triangles with a
    !> common vertex, edge or face (S. Sauter and C. Schwab, Boundary Element
    !> Methods, Springer 2011, Sec. 5.2). The singularity is removed by
    !> transformations of the four-dimensional integration domain, after which
    !> tensor Gauss rules of ssOrder points per direction converge
    !> exponentially. Reference triangle {0 <= x2 <= x1 <= 1} with
    !> chi(x) = P1 + x1 (P2 - P1) + x2 (P3 - P2) for both triangles.
    !---------------------------------------------------------------------------
    function triangle_triangle_sauter_schwab( triA, triB, edge ) result( val )
    real(ttReal),dimension(3,3),intent(in) :: triA, triB
    real(ttReal),intent(in) :: edge
    real(ttReal) :: val

    real(ttReal),dimension(3,3) :: A, B
    real(ttReal),dimension(ssOrder) :: x, w
    real(ttReal),dimension(2,6) :: xa, xb
    real(ttReal),dimension(6) :: jac
    integer,dimension(3) :: ia, ib
    integer,dimension(2,3) :: pairs
    integer :: i, j, nShared, nMaps, i1, i2, i3, i4, m
    real(ttReal) :: xi, e1, e2, e3, wq4, acc

        nShared = 0
        do i = 1, 3
            do j = 1, 3
                if ( norm2( triA(:,i) - triB(:,j) ) .le. sharedVertexTol * edge .and. nShared .lt. 3 ) then
                    nShared = nShared + 1
                    pairs(1,nShared) = i
                    pairs(2,nShared) = j
                endif
            enddo
        enddo

        !Common vertices first, in matching order.
        select case ( nShared )
        case ( 3 )
            ia = [ 1, 2, 3 ]
            ib = pairs(2,1:3)
        case ( 2 )
            ia = [ pairs(1,1), pairs(1,2), 6 - pairs(1,1) - pairs(1,2) ]
            ib = [ pairs(2,1), pairs(2,2), 6 - pairs(2,1) - pairs(2,2) ]
        case default
            ia = [ pairs(1,1), modulo( pairs(1,1), 3 ) + 1, modulo( pairs(1,1) + 1, 3 ) + 1 ]
            ib = [ pairs(2,1), modulo( pairs(2,1), 3 ) + 1, modulo( pairs(2,1) + 1, 3 ) + 1 ]
        end select
        A = triA(:,ia)
        B = triB(:,ib)

        call gauss_legendre01( ssOrder, x, w )
        acc = 0._ttReal
        do i1 = 1, ssOrder
            xi = x(i1)
            do i2 = 1, ssOrder
                e1 = x(i2)
                do i3 = 1, ssOrder
                    e2 = x(i3)
                    do i4 = 1, ssOrder
                        e3 = x(i4)
                        wq4 = w(i1) * w(i2) * w(i3) * w(i4)
                        select case ( nShared )
                        case ( 3 )
                            nMaps = 6
                            xa(:,1) = [ xi, xi * ( 1._ttReal - e1 + e1*e2 ) ]
                            xb(:,1) = [ xi * ( 1._ttReal - e1*e2*e3 ), xi * ( 1._ttReal - e1 ) ]
                            xa(:,3) = [ xi, xi * e1 * ( 1._ttReal - e2 + e2*e3 ) ]
                            xb(:,3) = [ xi * ( 1._ttReal - e1*e2 ), xi * e1 * ( 1._ttReal - e2 ) ]
                            xa(:,5) = [ xi * ( 1._ttReal - e1*e2*e3 ), xi * e1 * ( 1._ttReal - e2*e3 ) ]
                            xb(:,5) = [ xi, xi * e1 * ( 1._ttReal - e2 ) ]
                            xa(:,2) = xb(:,1)
                            xb(:,2) = xa(:,1)
                            xa(:,4) = xb(:,3)
                            xb(:,4) = xa(:,3)
                            xa(:,6) = xb(:,5)
                            xb(:,6) = xa(:,5)
                            jac = xi**3 * e1**2 * e2
                        case ( 2 )
                            nMaps = 5
                            xa(:,1) = [ xi, xi * e1 * e3 ]
                            xb(:,1) = [ xi * ( 1._ttReal - e1*e2 ), xi * e1 * ( 1._ttReal - e2 ) ]
                            xa(:,2) = [ xi, xi * e1 ]
                            xb(:,2) = [ xi * ( 1._ttReal - e1*e2*e3 ), xi * e1 * e2 * ( 1._ttReal - e3 ) ]
                            xa(:,3) = [ xi * ( 1._ttReal - e1*e2 ), xi * e1 * ( 1._ttReal - e2 ) ]
                            xb(:,3) = [ xi, xi * e1 * e2 * e3 ]
                            xa(:,4) = [ xi * ( 1._ttReal - e1*e2*e3 ), xi * e1 * e2 * ( 1._ttReal - e3 ) ]
                            xb(:,4) = [ xi, xi * e1 ]
                            xa(:,5) = [ xi * ( 1._ttReal - e1*e2*e3 ), xi * e1 * ( 1._ttReal - e2*e3 ) ]
                            xb(:,5) = [ xi, xi * e1 * e2 ]
                            jac(1) = xi**3 * e1**2
                            jac(2:5) = xi**3 * e1**2 * e2
                        case default
                            nMaps = 2
                            xa(:,1) = [ xi, xi * e1 ]
                            xb(:,1) = [ xi * e2, xi * e2 * e3 ]
                            xa(:,2) = xb(:,1)
                            xb(:,2) = xa(:,1)
                            jac(1:2) = xi**3 * e2
                        end select
                        do m = 1, nMaps
                            acc = acc + wq4 * jac(m) / norm2( &
                                ( A(:,1) + xa(1,m) * ( A(:,2) - A(:,1) ) + xa(2,m) * ( A(:,3) - A(:,2) ) ) &
                              - ( B(:,1) + xb(1,m) * ( B(:,2) - B(:,1) ) + xb(2,m) * ( B(:,3) - B(:,2) ) ) )
                        enddo
                    enddo
                enddo
            enddo
        enddo
        val = twice_area( A ) * twice_area( B ) * acc

    end function triangle_triangle_sauter_schwab

end module TileTriangleTriangle
