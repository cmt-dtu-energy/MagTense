!---------------------------------------------------------------------------
!> @brief
!> Analytical, constant-density Galerkin Laplace integral over a pair of flat
!> triangles,
!>
!>     I = \int_{T_a} \int_{T_b} 1 / |r_a - r_b| dS_b dS_a ,
!>
!> evaluated in closed form. No numerical quadrature is used anywhere in this
!> module.
!>
!> The implementation is a Fortran port of the constant-density branch of the
!> dimensional-reduction hierarchy I3 -> I2s/I2t -> I1 -> I0 of
!>
!>   N. A. Gumerov, S. Kaneko and R. Duraiswami,
!>   "Recursive computation of the multipole expansions of layer potential
!>    integrals over simplices for efficient fast multipole accelerated
!>    boundary elements",
!>   SIAM J. Sci. Comput. 46 (2024), DOI 10.1137/23M1547688,
!>
!> following the same reduction used by the (MIT licensed) reference
!> implementation GalerkinLaplaceTriGS.m and by the validated C++ port of it in
!> the dip-fmm project (src/geometry/primitives/tetrahedron.cpp).
!>
!> The integral is finite for every pair of non-degenerate triangles, including
!> coincident, edge-sharing and vertex-sharing pairs, and the parallel-plane
!> configuration is handled by its own reduction branch.
!>
!> All arithmetic is done in an explicitly declared double-precision kind so
!> that the result does not depend on whether the build promotes the default
!> real kind (ifx -real-size 64, gfortran -fdefault-real-8).
!---------------------------------------------------------------------------
module TileTriangleTriangle
    use, intrinsic :: ieee_arithmetic
    implicit none

    private

    !> Working precision of the analytical triangle-triangle machinery.
    integer,parameter,public :: ttReal = selected_real_kind(15, 307)

    !> Relative threshold below which an expansion coefficient is treated as
    !> exactly zero and its term is dropped from the reduction. Dropping the
    !> term is what keeps the reduction finite for touching triangles, so this
    !> is part of the algorithm and not a tuning knob.
    real(ttReal),parameter :: zeroTol = 1.0e-14_ttReal

    !> Squared relative threshold used by the I0 case selection to decide which
    !> of the four heights are identically zero.
    real(ttReal),parameter :: i0Zero2 = 1.0e-28_ttReal

    public :: triangle_triangle_laplace_integral

    contains

    !---------------------------------------------------------------------------
    !> @brief Analytical Galerkin Laplace integral over a pair of triangles.
    !> @param[in] triA first triangle, triA(:,k) is the k'th vertex
    !> @param[in] triB second triangle, triB(:,k) is the k'th vertex
    !> @param[out] val the value of the double surface integral of 1/R
    !> @param[out] ok .false. if either triangle is degenerate or the reduction
    !>             produced a non-finite value; val is then not usable
    !>
    !> The result is symmetric in its two arguments and invariant under a
    !> cyclic or reflecting permutation of the vertices of either triangle.
    !---------------------------------------------------------------------------
    subroutine triangle_triangle_laplace_integral( triA, triB, val, ok )
    real(ttReal),dimension(3,3),intent(in) :: triA, triB
    real(ttReal),intent(out) :: val
    logical,intent(out) :: ok

    real(ttReal),dimension(3,3) :: first, second, normA, normB
    real(ttReal) :: scaleA, scaleB, scaleVal, invScale
    integer :: i

        val = 0._ttReal
        ok = .true.

        if ( .not. all( ieee_is_finite( triA ) ) .or. .not. all( ieee_is_finite( triB ) ) ) then
            ok = .false.
            return
        endif

        first = triA
        second = triB

        scaleA = triangle_scale( first )
        scaleB = triangle_scale( second )

        if ( .not. ( scaleA .gt. 0._ttReal ) .or. .not. ( scaleB .gt. 0._ttReal ) ) then
            ok = .false.
            return
        endif

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
        invScale = 1._ttReal / scaleVal

        !Translate to the first vertex of the first triangle and rescale to unit
        !size. Both operations are exact up to rounding and the integral scales
        !as length**3, so this only improves the conditioning of the reduction.
        do i = 1, 3
            normA(:,i) = ( first(:,i)  - first(:,1) ) * invScale
            normB(:,i) = ( second(:,i) - first(:,1) ) * invScale
        enddo

        call triangle_triangle_core( normA, normB, val, ok )

        val = scaleVal**3 * val

        if ( .not. ieee_is_finite( val ) ) ok = .false.

    end subroutine triangle_triangle_laplace_integral


    !> Largest edge length of a triangle, used as its characteristic size.
    function triangle_scale( tri ) result( s )
    real(ttReal),dimension(3,3),intent(in) :: tri
    real(ttReal) :: s
        s = max( norm2( tri(:,2) - tri(:,1) ), &
                 norm2( tri(:,3) - tri(:,2) ), &
                 norm2( tri(:,1) - tri(:,3) ) )
    end function triangle_scale


    function cross3( a, b ) result( c )
    real(ttReal),dimension(3),intent(in) :: a, b
    real(ttReal),dimension(3) :: c
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross3


    !---------------------------------------------------------------------------
    !> @brief The reduction itself, for two triangles of order-unity size.
    !>
    !> a1, a2 span the first triangle and a3, a4 the second one; e4 is the
    !> offset between their first vertices. The non-parallel branch reduces the
    !> four-dimensional integral with a vanishing out-of-plane height, while the
    !> parallel branch carries the (constant) separation of the two planes as
    !> the height h4 of the reduction.
    !---------------------------------------------------------------------------
    subroutine triangle_triangle_core( first, second, val, ok )
    real(ttReal),dimension(3,3),intent(in) :: first, second
    real(ttReal),intent(out) :: val
    logical,intent(inout) :: ok

    real(ttReal),dimension(3) :: a1, a2, a3, a4, e4, e4x, crossA, crossB, nA, nB, proj
    real(ttReal),dimension(3,4) :: basis
    real(ttReal),dimension(4) :: coeff
    real(ttReal) :: twiceAreaA, twiceAreaB, residual, parallelMeasure, h4, reduced
    real(ttReal) :: s0, s1, s2, s3, i31, i32, i33, i34, i35, i36

        val = 0._ttReal

        a1 = first(:,2)  - first(:,1)
        a2 = first(:,3)  - first(:,1)
        a3 = second(:,1) - second(:,2)
        a4 = second(:,1) - second(:,3)
        e4 = first(:,1)  - second(:,1)

        crossA = cross3( a1, a2 )
        crossB = cross3( a3, a4 )
        twiceAreaA = norm2( crossA )
        twiceAreaB = norm2( crossB )

        if ( .not. ( twiceAreaA .gt. 0._ttReal ) .or. .not. ( twiceAreaB .gt. 0._ttReal ) ) then
            ok = .false.
            return
        endif

        basis(:,1) = a1
        basis(:,2) = a2
        basis(:,3) = a3
        basis(:,4) = a4
        call expand_gram_schmidt( basis, 4, e4, coeff, proj, residual )
        s0 = coeff(1)
        s1 = coeff(2)
        s2 = coeff(3)
        s3 = coeff(4)
        e4x = proj

        nA = crossA / twiceAreaA
        nB = crossB / twiceAreaB
        parallelMeasure = min( norm2( nA + nB ), norm2( nA - nB ) )

        if ( parallelMeasure .gt. zeroTol ) then
            !The two triangle planes intersect: the reduction closes in the
            !plane of each triangle and there is no out-of-plane height.
            h4 = 0._ttReal
            call triangle_i3( a1, a2, a3,       e4,      h4, i31, ok )
            call triangle_i3( a1, a2, a4,       e4,      h4, i32, ok )
            call triangle_i3( a1, a2, a4 - a3,  e4 + a3, h4, i33, ok )
            call triangle_i3( a3, a4, a1,       e4,      h4, i34, ok )
            call triangle_i3( a3, a4, a2,       e4,      h4, i35, ok )
            call triangle_i3( a3, a4, a2 - a1,  e4 + a1, h4, i36, ok )
            reduced = ( 1._ttReal + s2 + s3 ) * i33 - s3 * i31 - s2 * i32 &
                    + ( 1._ttReal + s0 + s1 ) * i36 - s1 * i34 - s0 * i35
        else
            !Parallel planes: the component of e4 normal to the common plane is
            !the residual of the expansion and enters the reduction as a
            !constant height. The projection of e4 into the plane is the only
            !in-plane offset the reduction needs.
            h4 = residual
            call triangle_i3( a1, a2, a4 - a3,  e4x + a3, h4, i33, ok )
            call triangle_i3( a3, a4, a1,       e4x,      h4, i34, ok )
            call triangle_i3( a3, a4, a2,       e4x,      h4, i35, ok )
            call triangle_i3( a3, a4, a2 - a1,  e4x + a1, h4, i36, ok )
            reduced = i33 + ( 1._ttReal + s0 + s1 ) * i36 - s1 * i34 - s0 * i35
        endif

        val = 4._ttReal * ( 0.5_ttReal * twiceAreaA ) * ( 0.5_ttReal * twiceAreaB ) * reduced

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
    !---------------------------------------------------------------------------
    subroutine expand_gram_schmidt( vectors, n, val, coeff, projected, residual )
    real(ttReal),dimension(3,4),intent(in) :: vectors
    integer,intent(in) :: n
    real(ttReal),dimension(3),intent(in) :: val
    real(ttReal),dimension(4),intent(out) :: coeff
    real(ttReal),dimension(3),intent(out) :: projected
    real(ttReal),intent(out) :: residual

    real(ttReal),dimension(3,4) :: ortho
    real(ttReal),dimension(4) :: onorm2, residProj
    real(ttReal),dimension(4,4) :: cmat
    real(ttReal) :: vecScale, basisScale, rankTol, residTol, nrm2, len
    integer :: i, col, row

        ortho = 0._ttReal
        onorm2 = 0._ttReal
        cmat = 0._ttReal
        coeff = 0._ttReal
        residProj = 0._ttReal

        vecScale = norm2( val )
        basisScale = 0._ttReal
        do i = 1, n
            len = norm2( vectors(:,i) )
            vecScale = max( vecScale, len )
            basisScale = max( basisScale, len )
        enddo
        rankTol = 256._ttReal * epsilon(1._ttReal) * max( tiny(1._ttReal), basisScale**2 )

        do col = 1, n
            ortho(:,col) = vectors(:,col)
            do row = 1, col-1
                cmat(row,col) = dot_product( ortho(:,row), vectors(:,col) ) / onorm2(row)
                ortho(:,col) = ortho(:,col) - cmat(row,col) * ortho(:,row)
            enddo
            nrm2 = dot_product( ortho(:,col), ortho(:,col) )
            if ( .not. ( nrm2 .gt. rankTol ) ) then
                ortho(:,col) = 0._ttReal
                onorm2(col) = 1._ttReal
                cmat(col,col) = 0._ttReal
            else
                onorm2(col) = nrm2
                cmat(col,col) = dot_product( ortho(:,col), vectors(:,col) ) / nrm2
            endif
        enddo

        do col = n, 1, -1
            if ( cmat(col,col) .eq. 0._ttReal ) cycle
            coeff(col) = ( dot_product( val, ortho(:,col) ) - residProj(col) ) &
                       / ( cmat(col,col) * onorm2(col) )
            do row = 1, col-1
                residProj(row) = residProj(row) + coeff(col) * cmat(row,col) * onorm2(row)
            enddo
        enddo

        projected = 0._ttReal
        do col = 1, n
            projected = projected + coeff(col) * vectors(:,col)
        enddo

        residual = norm2( val - projected )
        residTol = 256._ttReal * epsilon(1._ttReal) * max( 1._ttReal, vecScale )
        if ( residual .le. residTol ) residual = 0._ttReal

    end subroutine expand_gram_schmidt


    !---------------------------------------------------------------------------
    !> @brief Third-level reduction. Reduces the pair integral to five
    !> second-level integrals over the edges spanned by first, second, third.
    !---------------------------------------------------------------------------
    subroutine triangle_i3( first, second, third, offset, h4, res, ok )
    real(ttReal),dimension(3),intent(in) :: first, second, third, offset
    real(ttReal),intent(in) :: h4
    real(ttReal),intent(out) :: res
    logical,intent(inout) :: ok

    real(ttReal),dimension(3,4) :: basis
    real(ttReal),dimension(4) :: coeff
    real(ttReal),dimension(3) :: proj
    real(ttReal) :: residual, s0, s1, s2, term

        res = 0._ttReal
        basis = 0._ttReal
        basis(:,1) = first
        basis(:,2) = second
        basis(:,3) = third
        call expand_gram_schmidt( basis, 3, offset, coeff, proj, residual )
        s0 = coeff(1)
        s1 = coeff(2)
        s2 = coeff(3)

        if ( abs( 1._ttReal + s0 + s1 ) .gt. zeroTol ) then
            call triangle_i2s( first - second, third, proj + second, residual, h4, term, ok )
            res = res + ( 1._ttReal + s0 + s1 ) * term
        endif
        if ( abs( s0 ) .gt. zeroTol ) then
            call triangle_i2s( second, third, proj, residual, h4, term, ok )
            res = res - s0 * term
        endif
        if ( abs( s1 ) .gt. zeroTol ) then
            call triangle_i2s( first, third, proj, residual, h4, term, ok )
            res = res - s1 * term
        endif
        if ( abs( 1._ttReal + s2 ) .gt. zeroTol ) then
            call triangle_i2t( first, second, proj + third, residual, h4, term, ok )
            res = res + ( 1._ttReal + s2 ) * term
        endif
        if ( abs( s2 ) .gt. zeroTol ) then
            call triangle_i2t( first, second, proj, residual, h4, term, ok )
            res = res - s2 * term
        endif

    end subroutine triangle_i3


    !> Second-level reduction over a parallelogram-type pair of directions.
    subroutine triangle_i2s( first, second, offset, h3, h4, res, ok )
    real(ttReal),dimension(3),intent(in) :: first, second, offset
    real(ttReal),intent(in) :: h3, h4
    real(ttReal),intent(out) :: res
    logical,intent(inout) :: ok

    real(ttReal),dimension(3,4) :: basis
    real(ttReal),dimension(4) :: coeff
    real(ttReal),dimension(3) :: proj
    real(ttReal) :: residual, s0, s1, term

        res = 0._ttReal
        basis = 0._ttReal
        basis(:,1) = first
        basis(:,2) = second
        call expand_gram_schmidt( basis, 2, offset, coeff, proj, residual )
        s0 = coeff(1)
        s1 = coeff(2)

        if ( abs( 1._ttReal + s0 ) .gt. zeroTol ) then
            call triangle_i1( second, proj + first, residual, h3, h4, term, ok )
            res = res + ( 1._ttReal + s0 ) * term
        endif
        if ( abs( s0 ) .gt. zeroTol ) then
            call triangle_i1( second, proj, residual, h3, h4, term, ok )
            res = res - s0 * term
        endif
        if ( abs( 1._ttReal + s1 ) .gt. zeroTol ) then
            call triangle_i1( first, proj + second, residual, h3, h4, term, ok )
            res = res + ( 1._ttReal + s1 ) * term
        endif
        if ( abs( s1 ) .gt. zeroTol ) then
            call triangle_i1( first, proj, residual, h3, h4, term, ok )
            res = res - s1 * term
        endif

    end subroutine triangle_i2s


    !> Second-level reduction over a triangular-type pair of directions.
    subroutine triangle_i2t( first, second, offset, h3, h4, res, ok )
    real(ttReal),dimension(3),intent(in) :: first, second, offset
    real(ttReal),intent(in) :: h3, h4
    real(ttReal),intent(out) :: res
    logical,intent(inout) :: ok

    real(ttReal),dimension(3,4) :: basis
    real(ttReal),dimension(4) :: coeff
    real(ttReal),dimension(3) :: proj
    real(ttReal) :: residual, s0, s1, term

        res = 0._ttReal
        basis = 0._ttReal
        basis(:,1) = first
        basis(:,2) = second
        call expand_gram_schmidt( basis, 2, offset, coeff, proj, residual )
        s0 = coeff(1)
        s1 = coeff(2)

        if ( abs( s0 ) .gt. zeroTol ) then
            call triangle_i1( second, proj, residual, h3, h4, term, ok )
            res = res - s0 * term
        endif
        if ( abs( s1 ) .gt. zeroTol ) then
            call triangle_i1( first, proj, residual, h3, h4, term, ok )
            res = res - s1 * term
        endif
        if ( abs( 1._ttReal + s0 + s1 ) .gt. zeroTol ) then
            call triangle_i1( first - second, proj + second, residual, h3, h4, term, ok )
            res = res + ( 1._ttReal + s0 + s1 ) * term
        endif

    end subroutine triangle_i2t


    !> First-level reduction along a single edge direction.
    subroutine triangle_i1( vec, offset, h2, h3, h4, res, ok )
    real(ttReal),dimension(3),intent(in) :: vec, offset
    real(ttReal),intent(in) :: h2, h3, h4
    real(ttReal),intent(out) :: res
    logical,intent(inout) :: ok

    real(ttReal),dimension(3,4) :: basis
    real(ttReal),dimension(4) :: coeff
    real(ttReal),dimension(3) :: proj
    real(ttReal) :: residual, len, c, term

        res = 0._ttReal
        len = norm2( vec )
        !A vanishing edge contributes nothing; returning here is what keeps the
        !degenerate reductions out of the I0 primitive.
        if ( .not. ( len .gt. 0._ttReal ) ) return

        basis = 0._ttReal
        basis(:,1) = vec
        call expand_gram_schmidt( basis, 1, offset, coeff, proj, residual )
        c = coeff(1)

        if ( abs( 1._ttReal + c ) .gt. zeroTol ) then
            call triangle_i0( abs( 1._ttReal + c ) * len, [ residual, h2, h3, h4 ], term, ok )
            res = res + ( 1._ttReal + c ) * term
        endif
        if ( abs( c ) .gt. zeroTol ) then
            call triangle_i0( abs( c ) * len, [ residual, h2, h3, h4 ], term, ok )
            res = res - c * term
        endif

    end subroutine triangle_i1


    !> Helper for the I0 case analysis, the reference Phi2 function.
    function i0_phi2( h, h1, p, hh, rad ) result( phi2 )
    real(ttReal),intent(in) :: h, h1, p, hh, rad
    real(ttReal) :: phi2
        if ( h1 * p .eq. 0._ttReal ) then
            phi2 = 1._ttReal / ( hh + h * rad )
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
    real(ttReal),intent(in) :: p
    real(ttReal),dimension(4),intent(in) :: hIn
    real(ttReal),intent(out) :: ans
    logical,intent(inout) :: ok

    real(ttReal),dimension(4) :: h
    integer,dimension(4) :: ord
    integer :: i, j, key
    real(ttReal) :: h1, h2, h3, h4, pAbs, hh, hh1, rad, rad1, rad2, hNorm
    real(ttReal) :: phi1, phi2, phi3, phi4, q, lg

        ans = 0._ttReal
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
        h(ord(1)) = 0._ttReal
        h(ord(2)) = 0._ttReal

        h1 = h(1)
        h2 = h(2)
        h3 = h(3)
        h4 = h(4)

        pAbs = abs( p )
        hh = h1*h1 + h2*h2 + h3*h3 + h4*h4

        !Only reachable through a vanishing edge in a degenerate triangle; the
        !valid reductions drop such a term before they reach I0.
        if ( pAbs .eq. 0._ttReal .and. hh .eq. 0._ttReal ) return

        rad = sqrt( pAbs*pAbs + hh )

        !The complete zero-height branch of the reference reduction, not merely
        !an alternative way of evaluating Phi1.
        if ( hh .lt. i0Zero2 * pAbs * pAbs ) then
            ans = log( pAbs ) / ( 6._ttReal * pAbs )
            return
        endif

        if ( pAbs .eq. 0._ttReal ) then
            phi1 = 1._ttReal / sqrt( hh )
        else
            !asinh(p/sqrt(hh))/p is the logarithm of the reference I0 written in
            !a form that stays well behaved as p -> 0.
            phi1 = asinh( pAbs / sqrt( hh ) ) / pAbs
        endif

        hh1 = hh - h1*h1

        if ( h1*h1 + h2*h2 + h4*h4 .lt. i0Zero2 * hh ) then
            !Rank-deficient reductions can leave h3 as the only non-zero height.
            !The reference case-4 expression is 0/0 at h1 = 0; this is its
            !analytical h1 -> 0 limit.
            q = hypot( pAbs, h3 )
            lg = asinh( pAbs / h3 )
            ans = ( lg * ( 2._ttReal*q*q - 3._ttReal*h3*h3 ) / ( 2._ttReal * pAbs**3 ) &
                  + ( 4._ttReal*h3 - 3._ttReal*q ) / ( 2._ttReal * pAbs*pAbs ) ) / 6._ttReal

        elseif ( h1*h1 + h2*h2 + h3*h3 .lt. i0Zero2 * hh ) then
            !Likewise the h1 -> 0 limit of the reference case-6 expression, for
            !h4 being the only non-zero height.
            q = hypot( pAbs, h4 )
            lg = asinh( pAbs / h4 )
            ans = ( lg * ( 2._ttReal*q*q - 5._ttReal*h4*h4 ) / ( 2._ttReal * pAbs**3 ) &
                  + 3._ttReal * ( 2._ttReal*h4 - q ) / ( 2._ttReal * pAbs*pAbs ) &
                  - ( q + 2._ttReal*h4 ) / ( 3._ttReal * ( q + h4 )**2 ) ) / 6._ttReal

        elseif ( hh1 .lt. i0Zero2 * hh ) then
            ans = phi1 / 6._ttReal

        elseif ( h3*h3 + h4*h4 .lt. i0Zero2 * hh ) then
            ans = ( phi1 - h2 * i0_phi2( h2, h1, pAbs, hh, rad ) ) / 6._ttReal

        else
            rad1 = sqrt( pAbs*pAbs + h1*h1 )

            if ( h2*h2 + h4*h4 .lt. i0Zero2 * hh ) then
                phi2 = i0_phi2( h3, h1, pAbs, hh, rad )
                phi3 = 0.5_ttReal * hh1 / rad1 * log( ( rad1 + rad )**2 / hh1 )
                ans = ( ( h1*h1 - h3*h3 ) * phi1 - 2._ttReal*h1*h1*h3 * phi2 + phi3 ) &
                    / ( 6._ttReal * h1*h1 )

            elseif ( h2*h2 + h3*h3 .lt. i0Zero2 * hh ) then
                phi2 = i0_phi2( h4, h1, pAbs, hh, rad )
                phi3 = 0.5_ttReal / rad1 * log( ( rad1 + rad )**2 / hh1 )
                ans = ( ( h1*h1 - 3._ttReal*h4*h4 ) * phi1 &
                      - h4 * ( 3._ttReal*h1*h1 - h4*h4 ) * phi2 &
                      + 3._ttReal*h4*h4 * phi3 - h4*h4 / ( rad + h4 ) ) &
                    / ( 6._ttReal * h1*h1 )

            elseif ( h1*h1 + h4*h4 .lt. i0Zero2 * hh ) then
                hNorm = sqrt( hh )
                rad2 = sqrt( pAbs*pAbs + h2*h2 )
                phi4 = h3*h3 / ( h2 * pAbs*pAbs ) &
                     * ( ( rad2 / h2 ) * log( ( rad2 + rad ) / h3 ) - log( ( h2 + hNorm ) / h3 ) )
                ans = ( hh / ( h2*h2 ) * phi1 - 1._ttReal / ( rad + hNorm ) - phi4 ) / 6._ttReal

            elseif ( h1*h1 + h3*h3 .lt. i0Zero2 * hh ) then
                hNorm = sqrt( hh )
                rad2 = sqrt( pAbs*pAbs + h2*h2 )
                phi2 = atan( h2 * pAbs / ( hh + h4 * rad ) ) / pAbs
                phi4 = h4*h4 / ( h2 * pAbs*pAbs ) &
                     * ( ( rad2 / h2 ) * log( ( rad2 + rad ) / h4 ) - log( ( h2 + hNorm ) / h4 ) )
                ans = ( ( 1._ttReal + 3._ttReal*h4*h4 / ( h2*h2 ) ) * phi1 &
                      - 2._ttReal * ( h4 / h2 )**3 * phi2 - 3._ttReal * phi4 &
                      + ( 2._ttReal*h4*h4 - h2*h2 ) / ( h2*h2 * ( rad + hNorm ) ) ) / 6._ttReal

            else
                ok = .false.
                ans = ieee_value( 1._ttReal, ieee_quiet_nan )
            endif
        endif

    end subroutine triangle_i0

end module TileTriangleTriangle
