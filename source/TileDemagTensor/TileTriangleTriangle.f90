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
!>   - the reduction is evaluated in quadruple precision (kind wq), with
!>     tolerances chosen for it (see the parameters below);
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
module TileTriangleTriangle
    use, intrinsic :: ieee_arithmetic
    implicit none

    private

    !> Precision of the interface and of the quadratures.
    integer,parameter,public :: ttReal = selected_real_kind(15, 307)

    !> Working precision of the analytical reduction (quadruple). Its cost is
    !> about twice that of the double-precision reduction (ifx, 19 Oct 2026).
    integer,parameter :: wq = selected_real_kind(33, 4931)

    !> Relative threshold below which an expansion coefficient is treated as
    !> exactly zero and its term is dropped from the reduction. Dropping the
    !> term is what keeps the reduction finite for touching triangles. The
    !> double-precision value 1e-14 dropped real terms of nearly degenerate pairs
    !> (errors up to 7e-8); 1e-22 is still far above quadruple round-off.
    real(wq),parameter :: zeroTol = 1.0e-22_wq

    !> Squared relative threshold used by the I0 case selection to decide which
    !> of the four heights are identically zero.
    real(wq),parameter :: i0Z2 = 1.0e-24_wq

    !> Heights below this fraction of the vector scale are zero. A nearly shared
    !> vertex of double-precision input leaves heights of a few double ulps,
    !> which must be snapped; the I0 forms lose about eps/h**2, which this keeps
    !> below 1e-10 in quadruple precision.
    real(wq),parameter :: relResidual = 1.0e-12_wq

    !> Gram-Schmidt rank tolerance, squared and relative to the largest basis
    !> vector. Mesh geometry produces nearly dependent directions at every
    !> scale; directions down to 1e-23 must be kept and those below it dropped
    !> on the validation set (a decade either way fails on some pairs).
    real(wq),parameter :: rankRel2 = 1.0e-23_wq

    !> A rank decision with a squared relative norm in this band is unreliable:
    !> the face pair is then integrated numerically (about 2 % of the touching
    !> face pairs of an unstructured film mesh).
    real(wq),parameter :: guardLow = 1.0e-26_wq
    real(wq),parameter :: guardHigh = 1.0e-20_wq

    !> Face pairs whose normals differ by less than this are projected onto a
    !> common plane and evaluated by the parallel branch. The projection costs
    !> about 0.1 x angle, the general branch's round-off grows like
    !> eps / angle**3; both are below 1e-9 here in quadruple precision.
    real(wq),parameter :: nearParallel = 1.0e-8_wq

    !> Triangles without a common vertex whose normals differ by more than
    !> round-off but less than this are integrated by Gauss quadrature.
    real(ttReal),parameter :: parallelQuadrature = 1.0e-2_ttReal

    !> Vertices closer than this fraction of the larger edge are common.
    real(ttReal),parameter :: sharedVertexTol = 1.0e-10_ttReal

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

    real(wq),dimension(3,3) :: first, second, normA, normB
    real(wq) :: scaleA, scaleB, scaleVal, invScale, valQ
    real(ttReal),dimension(3) :: nA, nB
    real(ttReal) :: pm, edge
    integer :: i, j, nShared
    logical :: guard

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
        do i = 1, 3
            do j = 1, 3
                if ( norm2( triA(:,i) - triB(:,j) ) .le. sharedVertexTol * edge ) nShared = nShared + 1
            enddo
        enddo

        !Separated, nearly parallel triangles: Gauss quadrature (see the header).
        if ( nShared .eq. 0 .and. pm .gt. 1.0e-14_ttReal .and. pm .lt. parallelQuadrature ) then
            val = triangle_triangle_regular( triA, triB, 0 )
            if ( present(method) ) method = ttGaussParallel
            if ( .not. ieee_is_finite( val ) ) ok = .false.
            return
        endif

        first = real( triA, wq )
        second = real( triB, wq )

        scaleA = triangle_scale( first )
        scaleB = triangle_scale( second )

        !The core is symmetric in its two triangles. Put the larger one first,
        !which is the ordering the reference reduction is written for, and do so
        !before the common origin is chosen so the normalisation stays one-pass.
        if ( scaleA .lt. scaleB ) then
            normA = first
            first = second
            second = normA
            scaleVal = scaleA
            scaleA = scaleB
            scaleB = scaleVal
        endif

        scaleVal = max( scaleA, scaleB )
        invScale = 1._wq / scaleVal

        !Translate to the first vertex of the first triangle and rescale to unit
        !size. Both operations are exact up to rounding and the integral scales
        !as length**3, so this only improves the conditioning of the reduction.
        do i = 1, 3
            normA(:,i) = ( first(:,i)  - first(:,1) ) * invScale
            normB(:,i) = ( second(:,i) - first(:,1) ) * invScale
        enddo

        guard = .false.
        call triangle_triangle_core( normA, normB, valQ, ok, guard )

        if ( guard ) then
            !An unreliable rank decision: integrate numerically instead.
            ok = .true.
            if ( nShared .gt. 0 ) then
                val = triangle_triangle_sauter_schwab( triA, triB, edge )
            else
                val = triangle_triangle_regular( triA, triB, 0 )
            endif
            if ( present(method) ) method = ttGuardQuadrature
        else
            val = real( scaleVal**3 * valQ, ttReal )
        endif

        if ( .not. ieee_is_finite( val ) ) ok = .false.

    end subroutine triangle_triangle_laplace_integral


    !> Largest edge length of a triangle, used as its characteristic size.
    function triangle_scale( tri ) result( s )
    real(wq),dimension(3,3),intent(in) :: tri
    real(wq) :: s
        s = max( norm2( tri(:,2) - tri(:,1) ), &
                 norm2( tri(:,3) - tri(:,2) ), &
                 norm2( tri(:,1) - tri(:,3) ) )
    end function triangle_scale


    function triangle_scale_dp( tri ) result( s )
    real(ttReal),dimension(3,3),intent(in) :: tri
    real(ttReal) :: s
        s = max( norm2( tri(:,2) - tri(:,1) ), &
                 norm2( tri(:,3) - tri(:,2) ), &
                 norm2( tri(:,1) - tri(:,3) ) )
    end function triangle_scale_dp


    function cross3( a, b ) result( c )
    real(wq),dimension(3),intent(in) :: a, b
    real(wq),dimension(3) :: c
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross3


    function cross_dp( a, b ) result( c )
    real(ttReal),dimension(3),intent(in) :: a, b
    real(ttReal),dimension(3) :: c
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_dp


    !---------------------------------------------------------------------------
    !> @brief The reduction itself, for two triangles of order-unity size.
    !>
    !> a1, a2 span the first triangle and a3, a4 the second one; e4 is the
    !> offset between their first vertices. The non-parallel branch reduces the
    !> four-dimensional integral with a vanishing out-of-plane height, while the
    !> parallel branch carries the (constant) separation of the two planes as
    !> the height h4 of the reduction. guard is set when a rank decision of the
    !> reduction is unreliable (see guardLow).
    !---------------------------------------------------------------------------
    recursive subroutine triangle_triangle_core( first, second, val, ok, guard )
    real(wq),dimension(3,3),intent(in) :: first, second
    real(wq),intent(out) :: val
    logical,intent(inout) :: ok, guard

    real(wq),dimension(3) :: a1, a2, a3, a4, e4, e4x, crossA, crossB, nA, nB, proj, nMean
    real(wq),dimension(3,4) :: basis
    real(wq),dimension(4) :: coeff
    real(wq),dimension(3,3) :: projA, projB
    real(wq) :: twiceAreaA, twiceAreaB, residual, parallelMeasure, h4, reduced
    real(wq) :: s0, s1, s2, s3, i31, i32, i33, i34, i35, i36
    integer :: k

        val = 0._wq

        a1 = first(:,2)  - first(:,1)
        a2 = first(:,3)  - first(:,1)
        a3 = second(:,1) - second(:,2)
        a4 = second(:,1) - second(:,3)
        e4 = first(:,1)  - second(:,1)

        crossA = cross3( a1, a2 )
        crossB = cross3( a3, a4 )
        twiceAreaA = norm2( crossA )
        twiceAreaB = norm2( crossB )

        if ( .not. ( twiceAreaA .gt. 0._wq ) .or. .not. ( twiceAreaB .gt. 0._wq ) ) then
            ok = .false.
            return
        endif

        basis(:,1) = a1
        basis(:,2) = a2
        basis(:,3) = a3
        basis(:,4) = a4
        call expand_gram_schmidt( basis, 4, e4, coeff, proj, residual, guard )
        s0 = coeff(1)
        s1 = coeff(2)
        s2 = coeff(3)
        s3 = coeff(4)
        e4x = proj

        nA = crossA / twiceAreaA
        nB = crossB / twiceAreaB
        parallelMeasure = min( norm2( nA + nB ), norm2( nA - nB ) )

        !Nearly parallel faces (coplanar faces of a mesh after rounding, or
        !genuinely almost parallel ones): project both triangles onto planes
        !normal to their mean normal, each through its own centroid, and let the
        !parallel branch evaluate them. That costs O(angle) and is symmetric in
        !the two triangles; the general branch cancels like eps / angle**3 here.
        if ( parallelMeasure .gt. zeroTol .and. parallelMeasure .le. nearParallel ) then
            if ( dot_product( nA, nB ) .ge. 0._wq ) then
                nMean = nA + nB
            else
                nMean = nA - nB
            endif
            nMean = nMean / norm2( nMean )
            do k = 1, 3
                projA(:,k) = first(:,k) - dot_product( first(:,k) - sum( first, dim=2 ) / 3._wq, nMean ) * nMean
                projB(:,k) = second(:,k) - dot_product( second(:,k) - sum( second, dim=2 ) / 3._wq, nMean ) * nMean
            enddo
            call triangle_triangle_core( projA, projB, val, ok, guard )
            return
        endif

        if ( parallelMeasure .gt. zeroTol ) then
            !The two triangle planes intersect: the reduction closes in the
            !plane of each triangle and there is no out-of-plane height.
            h4 = 0._wq
            call triangle_i3( a1, a2, a3,       e4,      h4, i31, ok, guard )
            call triangle_i3( a1, a2, a4,       e4,      h4, i32, ok, guard )
            call triangle_i3( a1, a2, a4 - a3,  e4 + a3, h4, i33, ok, guard )
            call triangle_i3( a3, a4, a1,       e4,      h4, i34, ok, guard )
            call triangle_i3( a3, a4, a2,       e4,      h4, i35, ok, guard )
            call triangle_i3( a3, a4, a2 - a1,  e4 + a1, h4, i36, ok, guard )
            reduced = ( 1._wq + s2 + s3 ) * i33 - s3 * i31 - s2 * i32 &
                    + ( 1._wq + s0 + s1 ) * i36 - s1 * i34 - s0 * i35
        else
            !Parallel planes: the component of e4 normal to the common plane is
            !the residual of the expansion and enters the reduction as a
            !constant height. The projection of e4 into the plane is the only
            !in-plane offset the reduction needs.
            h4 = residual
            call triangle_i3( a1, a2, a4 - a3,  e4x + a3, h4, i33, ok, guard )
            call triangle_i3( a3, a4, a1,       e4x,      h4, i34, ok, guard )
            call triangle_i3( a3, a4, a2,       e4x,      h4, i35, ok, guard )
            call triangle_i3( a3, a4, a2 - a1,  e4x + a1, h4, i36, ok, guard )
            reduced = i33 + ( 1._wq + s0 + s1 ) * i36 - s1 * i34 - s0 * i35
        endif

        val = 4._wq * ( 0.5_wq * twiceAreaA ) * ( 0.5_wq * twiceAreaB ) * reduced

    end subroutine triangle_triangle_core


    !---------------------------------------------------------------------------
    !> @brief Small Gram-Schmidt expansion with explicit triangular back
    !> substitution.
    !>
    !> This mirrors the reference expan_GrammSchmidt routine. Keeping the back
    !> substitution explicit matters: I1/I2/I3 need the coefficients of val in
    !> the ORIGINAL (generally non-orthogonal) vectors, not the coefficients in
    !> the orthogonalised basis. A rank-deficient direction is detected with a
    !> scale-aware tolerance, zeroed, and skipped by the back substitution, so
    !> a degenerate configuration reduces rather than divides by zero.
    !>
    !> @param[out] coeff expansion coefficients of val in vectors(:,1:n)
    !> @param[out] projected the part of val spanned by vectors(:,1:n)
    !> @param[out] residual the norm of the unspanned remainder of val
    !> @param[inout] guard set if a rank decision falls in the guard band
    !---------------------------------------------------------------------------
    subroutine expand_gram_schmidt( vectors, n, val, coeff, projected, residual, guard )
    real(wq),dimension(3,4),intent(in) :: vectors
    integer,intent(in) :: n
    real(wq),dimension(3),intent(in) :: val
    real(wq),dimension(4),intent(out) :: coeff
    real(wq),dimension(3),intent(out) :: projected
    real(wq),intent(out) :: residual
    logical,intent(inout) :: guard

    real(wq),dimension(3,4) :: ortho
    real(wq),dimension(4) :: onorm2, residProj
    real(wq),dimension(4,4) :: cmat
    real(wq) :: vecScale, basisScale2, rankTol, residTol, nrm2, len, ratio
    integer :: i, col, row

        ortho = 0._wq
        onorm2 = 0._wq
        cmat = 0._wq
        coeff = 0._wq
        residProj = 0._wq

        vecScale = norm2( val )
        basisScale2 = 0._wq
        do i = 1, n
            len = norm2( vectors(:,i) )
            vecScale = max( vecScale, len )
            basisScale2 = max( basisScale2, len*len )
        enddo
        basisScale2 = max( tiny(1._wq), basisScale2 )
        rankTol = rankRel2 * basisScale2

        do col = 1, n
            ortho(:,col) = vectors(:,col)
            do row = 1, col-1
                cmat(row,col) = dot_product( ortho(:,row), vectors(:,col) ) / onorm2(row)
                ortho(:,col) = ortho(:,col) - cmat(row,col) * ortho(:,row)
            enddo
            nrm2 = dot_product( ortho(:,col), ortho(:,col) )
            ratio = nrm2 / basisScale2
            if ( ratio .ge. guardLow .and. ratio .lt. guardHigh ) guard = .true.
            if ( .not. ( nrm2 .gt. rankTol ) ) then
                ortho(:,col) = 0._wq
                onorm2(col) = 1._wq
                cmat(col,col) = 0._wq
            else
                onorm2(col) = nrm2
                cmat(col,col) = dot_product( ortho(:,col), vectors(:,col) ) / nrm2
            endif
        enddo

        do col = n, 1, -1
            if ( cmat(col,col) .eq. 0._wq ) cycle
            coeff(col) = ( dot_product( val, ortho(:,col) ) - residProj(col) ) &
                       / ( cmat(col,col) * onorm2(col) )
            do row = 1, col-1
                residProj(row) = residProj(row) + coeff(col) * cmat(row,col) * onorm2(row)
            enddo
        enddo

        projected = 0._wq
        do col = 1, n
            projected = projected + coeff(col) * vectors(:,col)
        enddo

        residual = norm2( val - projected )
        residTol = relResidual * max( 1._wq, vecScale )
        if ( residual .le. residTol ) residual = 0._wq

    end subroutine expand_gram_schmidt


    !---------------------------------------------------------------------------
    !> @brief Third-level reduction. Reduces the pair integral to five
    !> second-level integrals over the edges spanned by first, second, third.
    !---------------------------------------------------------------------------
    subroutine triangle_i3( first, second, third, offset, h4, res, ok, guard )
    real(wq),dimension(3),intent(in) :: first, second, third, offset
    real(wq),intent(in) :: h4
    real(wq),intent(out) :: res
    logical,intent(inout) :: ok, guard

    real(wq),dimension(3,4) :: basis
    real(wq),dimension(4) :: coeff
    real(wq),dimension(3) :: proj
    real(wq) :: residual, s0, s1, s2, term

        res = 0._wq
        basis = 0._wq
        basis(:,1) = first
        basis(:,2) = second
        basis(:,3) = third
        call expand_gram_schmidt( basis, 3, offset, coeff, proj, residual, guard )
        s0 = coeff(1)
        s1 = coeff(2)
        s2 = coeff(3)

        if ( abs( 1._wq + s0 + s1 ) .gt. zeroTol ) then
            call triangle_i2s( first - second, third, proj + second, residual, h4, term, ok, guard )
            res = res + ( 1._wq + s0 + s1 ) * term
        endif
        if ( abs( s0 ) .gt. zeroTol ) then
            call triangle_i2s( second, third, proj, residual, h4, term, ok, guard )
            res = res - s0 * term
        endif
        if ( abs( s1 ) .gt. zeroTol ) then
            call triangle_i2s( first, third, proj, residual, h4, term, ok, guard )
            res = res - s1 * term
        endif
        if ( abs( 1._wq + s2 ) .gt. zeroTol ) then
            call triangle_i2t( first, second, proj + third, residual, h4, term, ok, guard )
            res = res + ( 1._wq + s2 ) * term
        endif
        if ( abs( s2 ) .gt. zeroTol ) then
            call triangle_i2t( first, second, proj, residual, h4, term, ok, guard )
            res = res - s2 * term
        endif

    end subroutine triangle_i3


    !> Second-level reduction over a parallelogram-type pair of directions.
    subroutine triangle_i2s( first, second, offset, h3, h4, res, ok, guard )
    real(wq),dimension(3),intent(in) :: first, second, offset
    real(wq),intent(in) :: h3, h4
    real(wq),intent(out) :: res
    logical,intent(inout) :: ok, guard

    real(wq),dimension(3,4) :: basis
    real(wq),dimension(4) :: coeff
    real(wq),dimension(3) :: proj
    real(wq) :: residual, s0, s1, term

        res = 0._wq
        basis = 0._wq
        basis(:,1) = first
        basis(:,2) = second
        call expand_gram_schmidt( basis, 2, offset, coeff, proj, residual, guard )
        s0 = coeff(1)
        s1 = coeff(2)

        if ( abs( 1._wq + s0 ) .gt. zeroTol ) then
            call triangle_i1( second, proj + first, residual, h3, h4, term, ok, guard )
            res = res + ( 1._wq + s0 ) * term
        endif
        if ( abs( s0 ) .gt. zeroTol ) then
            call triangle_i1( second, proj, residual, h3, h4, term, ok, guard )
            res = res - s0 * term
        endif
        if ( abs( 1._wq + s1 ) .gt. zeroTol ) then
            call triangle_i1( first, proj + second, residual, h3, h4, term, ok, guard )
            res = res + ( 1._wq + s1 ) * term
        endif
        if ( abs( s1 ) .gt. zeroTol ) then
            call triangle_i1( first, proj, residual, h3, h4, term, ok, guard )
            res = res - s1 * term
        endif

    end subroutine triangle_i2s


    !> Second-level reduction over a triangular-type pair of directions.
    subroutine triangle_i2t( first, second, offset, h3, h4, res, ok, guard )
    real(wq),dimension(3),intent(in) :: first, second, offset
    real(wq),intent(in) :: h3, h4
    real(wq),intent(out) :: res
    logical,intent(inout) :: ok, guard

    real(wq),dimension(3,4) :: basis
    real(wq),dimension(4) :: coeff
    real(wq),dimension(3) :: proj
    real(wq) :: residual, s0, s1, term

        res = 0._wq
        basis = 0._wq
        basis(:,1) = first
        basis(:,2) = second
        call expand_gram_schmidt( basis, 2, offset, coeff, proj, residual, guard )
        s0 = coeff(1)
        s1 = coeff(2)

        if ( abs( s0 ) .gt. zeroTol ) then
            call triangle_i1( second, proj, residual, h3, h4, term, ok, guard )
            res = res - s0 * term
        endif
        if ( abs( s1 ) .gt. zeroTol ) then
            call triangle_i1( first, proj, residual, h3, h4, term, ok, guard )
            res = res - s1 * term
        endif
        if ( abs( 1._wq + s0 + s1 ) .gt. zeroTol ) then
            call triangle_i1( first - second, proj + second, residual, h3, h4, term, ok, guard )
            res = res + ( 1._wq + s0 + s1 ) * term
        endif

    end subroutine triangle_i2t


    !> First-level reduction along a single edge direction.
    subroutine triangle_i1( vec, offset, h2, h3, h4, res, ok, guard )
    real(wq),dimension(3),intent(in) :: vec, offset
    real(wq),intent(in) :: h2, h3, h4
    real(wq),intent(out) :: res
    logical,intent(inout) :: ok, guard

    real(wq),dimension(3,4) :: basis
    real(wq),dimension(4) :: coeff
    real(wq),dimension(3) :: proj
    real(wq) :: residual, len, c, term

        res = 0._wq
        len = norm2( vec )
        !A vanishing edge contributes nothing; returning here is what keeps the
        !degenerate reductions out of the I0 primitive.
        if ( .not. ( len .gt. 0._wq ) ) return

        basis = 0._wq
        basis(:,1) = vec
        call expand_gram_schmidt( basis, 1, offset, coeff, proj, residual, guard )
        c = coeff(1)

        if ( abs( 1._wq + c ) .gt. zeroTol ) then
            call triangle_i0( abs( 1._wq + c ) * len, [ residual, h2, h3, h4 ], term, ok )
            res = res + ( 1._wq + c ) * term
        endif
        if ( abs( c ) .gt. zeroTol ) then
            call triangle_i0( abs( c ) * len, [ residual, h2, h3, h4 ], term, ok )
            res = res - c * term
        endif

    end subroutine triangle_i1


    !> Helper for the I0 case analysis, the reference Phi2 function.
    function i0_phi2( h, h1, p, hh, rad ) result( phi2 )
    real(wq),intent(in) :: h, h1, p, hh, rad
    real(wq) :: phi2
        if ( h1 * p .eq. 0._wq ) then
            phi2 = 1._wq / ( hh + h * rad )
        else
            phi2 = atan( ( h1 * p ) / ( hh + h * rad ) ) / ( h1 * p )
        endif
    end function i0_phi2


    !---------------------------------------------------------------------------
    !> @brief The closed-form primitive at the bottom of the reduction.
    !>
    !> @param[in] p the (signed) reduced edge length; only its modulus is used
    !> @param[in] hIn the four reduction heights. At most two of them are
    !>            non-zero for a valid reduction, so the two smallest are set
    !>            to zero before the case analysis, exactly as in the reference.
    !>
    !> The branches are the reference I0 cases, plus the analytical h1 -> 0
    !> limits of the two cases whose reference expressions are 0/0 there. The
    !> final branch is the reference's prohibited case: reaching it means the
    !> preceding reduction lost its expected rank, so it is reported through ok
    !> rather than being given an ad-hoc value that would hide the error.
    !---------------------------------------------------------------------------
    subroutine triangle_i0( p, hIn, ans, ok )
    real(wq),intent(in) :: p
    real(wq),dimension(4),intent(in) :: hIn
    real(wq),intent(out) :: ans
    logical,intent(inout) :: ok

    real(wq),dimension(4) :: h
    integer,dimension(4) :: ord
    integer :: i, j, key
    real(wq) :: h1, h2, h3, h4, pAbs, hh, hh1, rad, rad1, rad2, hNorm
    real(wq) :: phi1, phi2, phi3, phi4, q, lg

        ans = 0._wq
        h = hIn

        !Stable ascending index sort of the four heights, then zero the two
        !smallest ones. Ties are broken by index so the result is deterministic.
        ord = [ 1, 2, 3, 4 ]
        do i = 2, 4
            key = ord(i)
            j = i - 1
            do while ( j .ge. 1 )
                if ( h(ord(j)) .gt. h(key) ) then
                    ord(j+1) = ord(j)
                    j = j - 1
                else
                    exit
                endif
            enddo
            ord(j+1) = key
        enddo
        h(ord(1)) = 0._wq
        h(ord(2)) = 0._wq

        h1 = h(1)
        h2 = h(2)
        h3 = h(3)
        h4 = h(4)

        pAbs = abs( p )
        hh = h1*h1 + h2*h2 + h3*h3 + h4*h4

        !Only reachable through a vanishing edge in a degenerate triangle; the
        !valid reductions drop such a term before they reach I0.
        if ( pAbs .eq. 0._wq .and. hh .eq. 0._wq ) return

        rad = sqrt( pAbs*pAbs + hh )

        !The complete zero-height branch of the reference reduction, not merely
        !an alternative way of evaluating Phi1.
        if ( hh .lt. i0Z2 * pAbs * pAbs ) then
            ans = log( pAbs ) / ( 6._wq * pAbs )
            return
        endif

        if ( pAbs .eq. 0._wq ) then
            phi1 = 1._wq / sqrt( hh )
        else
            !asinh(p/sqrt(hh))/p is the logarithm of the reference I0 written in
            !a form that stays well behaved as p -> 0.
            phi1 = asinh( pAbs / sqrt( hh ) ) / pAbs
        endif

        hh1 = hh - h1*h1

        if ( h1*h1 + h2*h2 + h4*h4 .lt. i0Z2 * hh ) then
            !Rank-deficient reductions can leave h3 as the only non-zero height.
            !The reference case-4 expression is 0/0 at h1 = 0; this is its
            !analytical h1 -> 0 limit.
            q = hypot( pAbs, h3 )
            lg = asinh( pAbs / h3 )
            ans = ( lg * ( 2._wq*q*q - 3._wq*h3*h3 ) / ( 2._wq * pAbs**3 ) &
                  + ( 4._wq*h3 - 3._wq*q ) / ( 2._wq * pAbs*pAbs ) ) / 6._wq

        elseif ( h1*h1 + h2*h2 + h3*h3 .lt. i0Z2 * hh ) then
            !Likewise the h1 -> 0 limit of the reference case-6 expression, for
            !h4 being the only non-zero height.
            q = hypot( pAbs, h4 )
            lg = asinh( pAbs / h4 )
            ans = ( lg * ( 2._wq*q*q - 5._wq*h4*h4 ) / ( 2._wq * pAbs**3 ) &
                  + 3._wq * ( 2._wq*h4 - q ) / ( 2._wq * pAbs*pAbs ) &
                  - ( q + 2._wq*h4 ) / ( 3._wq * ( q + h4 )**2 ) ) / 6._wq

        elseif ( hh1 .lt. i0Z2 * hh ) then
            ans = phi1 / 6._wq

        elseif ( h3*h3 + h4*h4 .lt. i0Z2 * hh ) then
            ans = ( phi1 - h2 * i0_phi2( h2, h1, pAbs, hh, rad ) ) / 6._wq

        else
            rad1 = sqrt( pAbs*pAbs + h1*h1 )

            if ( h2*h2 + h4*h4 .lt. i0Z2 * hh ) then
                phi2 = i0_phi2( h3, h1, pAbs, hh, rad )
                phi3 = 0.5_wq * hh1 / rad1 * log( ( rad1 + rad )**2 / hh1 )
                ans = ( ( h1*h1 - h3*h3 ) * phi1 - 2._wq*h1*h1*h3 * phi2 + phi3 ) &
                    / ( 6._wq * h1*h1 )

            elseif ( h2*h2 + h3*h3 .lt. i0Z2 * hh ) then
                phi2 = i0_phi2( h4, h1, pAbs, hh, rad )
                phi3 = 0.5_wq / rad1 * log( ( rad1 + rad )**2 / hh1 )
                ans = ( ( h1*h1 - 3._wq*h4*h4 ) * phi1 &
                      - h4 * ( 3._wq*h1*h1 - h4*h4 ) * phi2 &
                      + 3._wq*h4*h4 * phi3 - h4*h4 / ( rad + h4 ) ) &
                    / ( 6._wq * h1*h1 )

            elseif ( h1*h1 + h4*h4 .lt. i0Z2 * hh ) then
                hNorm = sqrt( hh )
                rad2 = sqrt( pAbs*pAbs + h2*h2 )
                phi4 = h3*h3 / ( h2 * pAbs*pAbs ) &
                     * ( ( rad2 / h2 ) * log( ( rad2 + rad ) / h3 ) - log( ( h2 + hNorm ) / h3 ) )
                ans = ( hh / ( h2*h2 ) * phi1 - 1._wq / ( rad + hNorm ) - phi4 ) / 6._wq

            elseif ( h1*h1 + h3*h3 .lt. i0Z2 * hh ) then
                hNorm = sqrt( hh )
                rad2 = sqrt( pAbs*pAbs + h2*h2 )
                phi2 = atan( h2 * pAbs / ( hh + h4 * rad ) ) / pAbs
                phi4 = h4*h4 / ( h2 * pAbs*pAbs ) &
                     * ( ( rad2 / h2 ) * log( ( rad2 + rad ) / h4 ) - log( ( h2 + hNorm ) / h4 ) )
                ans = ( ( 1._wq + 3._wq*h4*h4 / ( h2*h2 ) ) * phi1 &
                      - 2._wq * ( h4 / h2 )**3 * phi2 - 3._wq * phi4 &
                      + ( 2._wq*h4*h4 - h2*h2 ) / ( h2*h2 * ( rad + hNorm ) ) ) / 6._wq

            else
                ok = .false.
                ans = ieee_value( 1._wq, ieee_quiet_nan )
            endif
        endif

    end subroutine triangle_i0


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
