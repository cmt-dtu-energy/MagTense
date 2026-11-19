! Triangle-pair reduction of TileTriangleTriangle, included by the quadruple- and
! double-precision modules at the top of TileTriangleTriangle.f90. It works in the
! including module's kind wk (I0 in kind wi) with its tolerances zeroTol, i0Z2,
! relResidual, rankRel2, guardLow, guardHigh and nearParallel.

    !---------------------------------------------------------------------------
    !> @brief Reduce an ordered triangle pair (larger triangle first) in the
    !> working precision of this module: translate to the first vertex, scale to
    !> unit size, reduce, and rescale by scale**3 (the integral is homogeneous of
    !> degree 3 in length).
    !---------------------------------------------------------------------------
    subroutine reduce_triangle_pair( first, second, scaleVal, val, ok, guard )
    real(dp),dimension(3,3),intent(in) :: first, second
    real(dp),intent(in) :: scaleVal
    real(dp),intent(out) :: val
    logical,intent(out) :: ok, guard

    real(wk),dimension(3,3) :: normA, normB
    real(wk) :: invScale, valW
    integer :: i

        ok = .true.
        guard = .false.
        invScale = 1._wk / real( scaleVal, wk )
        do i = 1, 3
            normA(:,i) = ( real( first(:,i), wk ) - real( first(:,1), wk ) ) * invScale
            normB(:,i) = ( real( second(:,i), wk ) - real( first(:,1), wk ) ) * invScale
        enddo
        call triangle_triangle_core( normA, normB, valW, ok, guard )
        val = real( real( scaleVal, wk )**3 * valW, dp )

    end subroutine reduce_triangle_pair


    function cross3( a, b ) result( c )
    real(wk),dimension(3),intent(in) :: a, b
    real(wk),dimension(3) :: c
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
    !> the height h4 of the reduction. guard is set when a rank decision of the
    !> reduction is unreliable (see guardLow).
    !---------------------------------------------------------------------------
    recursive subroutine triangle_triangle_core( first, second, val, ok, guard )
    real(wk),dimension(3,3),intent(in) :: first, second
    real(wk),intent(out) :: val
    logical,intent(inout) :: ok, guard

    real(wk),dimension(3) :: a1, a2, a3, a4, e4, e4x, crossA, crossB, nA, nB, proj, nMean
    real(wk),dimension(3,4) :: basis
    real(wk),dimension(4) :: coeff
    real(wk),dimension(3,3) :: projA, projB
    real(wk) :: twiceAreaA, twiceAreaB, residual, parallelMeasure, h4, reduced
    real(wk) :: s0, s1, s2, s3, i31, i32, i33, i34, i35, i36
    integer :: k

        val = 0._wk

        a1 = first(:,2)  - first(:,1)
        a2 = first(:,3)  - first(:,1)
        a3 = second(:,1) - second(:,2)
        a4 = second(:,1) - second(:,3)
        e4 = first(:,1)  - second(:,1)

        crossA = cross3( a1, a2 )
        crossB = cross3( a3, a4 )
        twiceAreaA = norm2( crossA )
        twiceAreaB = norm2( crossB )

        if ( .not. ( twiceAreaA .gt. 0._wk ) .or. .not. ( twiceAreaB .gt. 0._wk ) ) then
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
            if ( dot_product( nA, nB ) .ge. 0._wk ) then
                nMean = nA + nB
            else
                nMean = nA - nB
            endif
            nMean = nMean / norm2( nMean )
            do k = 1, 3
                projA(:,k) = first(:,k) - dot_product( first(:,k) - sum( first, dim=2 ) / 3._wk, nMean ) * nMean
                projB(:,k) = second(:,k) - dot_product( second(:,k) - sum( second, dim=2 ) / 3._wk, nMean ) * nMean
            enddo
            call triangle_triangle_core( projA, projB, val, ok, guard )
            return
        endif

        if ( parallelMeasure .gt. zeroTol ) then
            !The two triangle planes intersect: the reduction closes in the
            !plane of each triangle and there is no out-of-plane height.
            h4 = 0._wk
            call triangle_i3( a1, a2, a3,       e4,      h4, i31, ok, guard )
            call triangle_i3( a1, a2, a4,       e4,      h4, i32, ok, guard )
            call triangle_i3( a1, a2, a4 - a3,  e4 + a3, h4, i33, ok, guard )
            call triangle_i3( a3, a4, a1,       e4,      h4, i34, ok, guard )
            call triangle_i3( a3, a4, a2,       e4,      h4, i35, ok, guard )
            call triangle_i3( a3, a4, a2 - a1,  e4 + a1, h4, i36, ok, guard )
            reduced = ( 1._wk + s2 + s3 ) * i33 - s3 * i31 - s2 * i32 &
                    + ( 1._wk + s0 + s1 ) * i36 - s1 * i34 - s0 * i35
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
            reduced = i33 + ( 1._wk + s0 + s1 ) * i36 - s1 * i34 - s0 * i35
        endif

        val = 4._wk * ( 0.5_wk * twiceAreaA ) * ( 0.5_wk * twiceAreaB ) * reduced

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
    real(wk),dimension(3,4),intent(in) :: vectors
    integer,intent(in) :: n
    real(wk),dimension(3),intent(in) :: val
    real(wk),dimension(4),intent(out) :: coeff
    real(wk),dimension(3),intent(out) :: projected
    real(wk),intent(out) :: residual
    logical,intent(inout) :: guard

    real(wk),dimension(3,4) :: ortho
    real(wk),dimension(4) :: onorm2, residProj
    real(wk),dimension(4,4) :: cmat
    real(wk) :: vecScale, basisScale2, rankTol, residTol, nrm2, len, ratio
    integer :: i, col, row

        ortho = 0._wk
        onorm2 = 0._wk
        cmat = 0._wk
        coeff = 0._wk
        residProj = 0._wk

        vecScale = norm2( val )
        basisScale2 = 0._wk
        do i = 1, n
            len = norm2( vectors(:,i) )
            vecScale = max( vecScale, len )
            basisScale2 = max( basisScale2, len*len )
        enddo
        basisScale2 = max( tiny(1._wk), basisScale2 )
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
                ortho(:,col) = 0._wk
                onorm2(col) = 1._wk
                cmat(col,col) = 0._wk
            else
                onorm2(col) = nrm2
                cmat(col,col) = dot_product( ortho(:,col), vectors(:,col) ) / nrm2
            endif
        enddo

        do col = n, 1, -1
            if ( cmat(col,col) .eq. 0._wk ) cycle
            coeff(col) = ( dot_product( val, ortho(:,col) ) - residProj(col) ) &
                       / ( cmat(col,col) * onorm2(col) )
            do row = 1, col-1
                residProj(row) = residProj(row) + coeff(col) * cmat(row,col) * onorm2(row)
            enddo
        enddo

        projected = 0._wk
        do col = 1, n
            projected = projected + coeff(col) * vectors(:,col)
        enddo

        residual = norm2( val - projected )
        residTol = relResidual * max( 1._wk, vecScale )
        if ( residual .le. residTol ) residual = 0._wk

    end subroutine expand_gram_schmidt


    !---------------------------------------------------------------------------
    !> @brief Third-level reduction. Reduces the pair integral to five
    !> second-level integrals over the edges spanned by first, second, third.
    !---------------------------------------------------------------------------
    subroutine triangle_i3( first, second, third, offset, h4, res, ok, guard )
    real(wk),dimension(3),intent(in) :: first, second, third, offset
    real(wk),intent(in) :: h4
    real(wk),intent(out) :: res
    logical,intent(inout) :: ok, guard

    real(wk),dimension(3,4) :: basis
    real(wk),dimension(4) :: coeff
    real(wk),dimension(3) :: proj
    real(wk) :: residual, s0, s1, s2, term

        res = 0._wk
        basis = 0._wk
        basis(:,1) = first
        basis(:,2) = second
        basis(:,3) = third
        call expand_gram_schmidt( basis, 3, offset, coeff, proj, residual, guard )
        s0 = coeff(1)
        s1 = coeff(2)
        s2 = coeff(3)

        if ( abs( 1._wk + s0 + s1 ) .gt. zeroTol ) then
            call triangle_i2s( first - second, third, proj + second, residual, h4, term, ok, guard )
            res = res + ( 1._wk + s0 + s1 ) * term
        endif
        if ( abs( s0 ) .gt. zeroTol ) then
            call triangle_i2s( second, third, proj, residual, h4, term, ok, guard )
            res = res - s0 * term
        endif
        if ( abs( s1 ) .gt. zeroTol ) then
            call triangle_i2s( first, third, proj, residual, h4, term, ok, guard )
            res = res - s1 * term
        endif
        if ( abs( 1._wk + s2 ) .gt. zeroTol ) then
            call triangle_i2t( first, second, proj + third, residual, h4, term, ok, guard )
            res = res + ( 1._wk + s2 ) * term
        endif
        if ( abs( s2 ) .gt. zeroTol ) then
            call triangle_i2t( first, second, proj, residual, h4, term, ok, guard )
            res = res - s2 * term
        endif

    end subroutine triangle_i3


    !> Second-level reduction over a parallelogram-type pair of directions.
    subroutine triangle_i2s( first, second, offset, h3, h4, res, ok, guard )
    real(wk),dimension(3),intent(in) :: first, second, offset
    real(wk),intent(in) :: h3, h4
    real(wk),intent(out) :: res
    logical,intent(inout) :: ok, guard

    real(wk),dimension(3,4) :: basis
    real(wk),dimension(4) :: coeff
    real(wk),dimension(3) :: proj
    real(wk) :: residual, s0, s1, term

        res = 0._wk
        basis = 0._wk
        basis(:,1) = first
        basis(:,2) = second
        call expand_gram_schmidt( basis, 2, offset, coeff, proj, residual, guard )
        s0 = coeff(1)
        s1 = coeff(2)

        if ( abs( 1._wk + s0 ) .gt. zeroTol ) then
            call triangle_i1( second, proj + first, residual, h3, h4, term, ok, guard )
            res = res + ( 1._wk + s0 ) * term
        endif
        if ( abs( s0 ) .gt. zeroTol ) then
            call triangle_i1( second, proj, residual, h3, h4, term, ok, guard )
            res = res - s0 * term
        endif
        if ( abs( 1._wk + s1 ) .gt. zeroTol ) then
            call triangle_i1( first, proj + second, residual, h3, h4, term, ok, guard )
            res = res + ( 1._wk + s1 ) * term
        endif
        if ( abs( s1 ) .gt. zeroTol ) then
            call triangle_i1( first, proj, residual, h3, h4, term, ok, guard )
            res = res - s1 * term
        endif

    end subroutine triangle_i2s


    !> Second-level reduction over a triangular-type pair of directions.
    subroutine triangle_i2t( first, second, offset, h3, h4, res, ok, guard )
    real(wk),dimension(3),intent(in) :: first, second, offset
    real(wk),intent(in) :: h3, h4
    real(wk),intent(out) :: res
    logical,intent(inout) :: ok, guard

    real(wk),dimension(3,4) :: basis
    real(wk),dimension(4) :: coeff
    real(wk),dimension(3) :: proj
    real(wk) :: residual, s0, s1, term

        res = 0._wk
        basis = 0._wk
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
        if ( abs( 1._wk + s0 + s1 ) .gt. zeroTol ) then
            call triangle_i1( first - second, proj + second, residual, h3, h4, term, ok, guard )
            res = res + ( 1._wk + s0 + s1 ) * term
        endif

    end subroutine triangle_i2t


    !> First-level reduction along a single edge direction.
    subroutine triangle_i1( vec, offset, h2, h3, h4, res, ok, guard )
    real(wk),dimension(3),intent(in) :: vec, offset
    real(wk),intent(in) :: h2, h3, h4
    real(wk),intent(out) :: res
    logical,intent(inout) :: ok, guard

    real(wk),dimension(3,4) :: basis
    real(wk),dimension(4) :: coeff
    real(wk),dimension(3) :: proj
    real(wk) :: residual, len, c, term

        res = 0._wk
        len = norm2( vec )
        !A vanishing edge contributes nothing; returning here is what keeps the
        !degenerate reductions out of the I0 primitive.
        if ( .not. ( len .gt. 0._wk ) ) return

        basis = 0._wk
        basis(:,1) = vec
        call expand_gram_schmidt( basis, 1, offset, coeff, proj, residual, guard )
        c = coeff(1)

        if ( abs( 1._wk + c ) .gt. zeroTol ) then
            call triangle_i0( abs( 1._wk + c ) * len, [ residual, h2, h3, h4 ], term, ok )
            res = res + ( 1._wk + c ) * term
        endif
        if ( abs( c ) .gt. zeroTol ) then
            call triangle_i0( abs( c ) * len, [ residual, h2, h3, h4 ], term, ok )
            res = res - c * term
        endif

    end subroutine triangle_i1


    !> Helper for the I0 case analysis, the reference Phi2 function.
    function i0_phi2( h, h1, p, hh, rad ) result( phi2 )
    real(wi),intent(in) :: h, h1, p, hh, rad
    real(wi) :: phi2
        if ( h1 * p .eq. 0._wi ) then
            phi2 = 1._wi / ( hh + h * rad )
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
    real(wk),intent(in) :: p
    real(wk),dimension(4),intent(in) :: hIn
    real(wk),intent(out) :: ans
    logical,intent(inout) :: ok

    real(wi),dimension(4) :: h
    integer,dimension(4) :: ord
    integer :: i, j, key
    real(wi) :: h1, h2, h3, h4, pAbs, hh, hh1, rad, rad1, rad2, hNorm
    real(wi) :: phi1, phi2, phi3, phi4, q, lg, res

        res = 0._wi
        ans = 0._wk
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
        h(ord(1)) = 0._wi
        h(ord(2)) = 0._wi

        h1 = h(1)
        h2 = h(2)
        h3 = h(3)
        h4 = h(4)

        pAbs = abs( p )
        hh = h1*h1 + h2*h2 + h3*h3 + h4*h4

        !Only reachable through a vanishing edge in a degenerate triangle; the
        !valid reductions drop such a term before they reach I0.
        if ( pAbs .eq. 0._wi .and. hh .eq. 0._wi ) return

        rad = sqrt( pAbs*pAbs + hh )

        !The complete zero-height branch of the reference reduction, not merely
        !an alternative way of evaluating Phi1.
        if ( hh .lt. i0Z2 * pAbs * pAbs ) then
            ans = real( log( pAbs ) / ( 6._wi * pAbs ), wk )
            return
        endif

        if ( pAbs .eq. 0._wi ) then
            phi1 = 1._wi / sqrt( hh )
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
            res = ( lg * ( 2._wi*q*q - 3._wi*h3*h3 ) / ( 2._wi * pAbs**3 ) &
                  + ( 4._wi*h3 - 3._wi*q ) / ( 2._wi * pAbs*pAbs ) ) / 6._wi

        elseif ( h1*h1 + h2*h2 + h3*h3 .lt. i0Z2 * hh ) then
            !Likewise the h1 -> 0 limit of the reference case-6 expression, for
            !h4 being the only non-zero height.
            q = hypot( pAbs, h4 )
            lg = asinh( pAbs / h4 )
            res = ( lg * ( 2._wi*q*q - 5._wi*h4*h4 ) / ( 2._wi * pAbs**3 ) &
                  + 3._wi * ( 2._wi*h4 - q ) / ( 2._wi * pAbs*pAbs ) &
                  - ( q + 2._wi*h4 ) / ( 3._wi * ( q + h4 )**2 ) ) / 6._wi

        elseif ( hh1 .lt. i0Z2 * hh ) then
            res = phi1 / 6._wi

        elseif ( h3*h3 + h4*h4 .lt. i0Z2 * hh ) then
            res = ( phi1 - h2 * i0_phi2( h2, h1, pAbs, hh, rad ) ) / 6._wi

        else
            rad1 = sqrt( pAbs*pAbs + h1*h1 )

            if ( h2*h2 + h4*h4 .lt. i0Z2 * hh ) then
                phi2 = i0_phi2( h3, h1, pAbs, hh, rad )
                phi3 = 0.5_wi * hh1 / rad1 * log( ( rad1 + rad )**2 / hh1 )
                res = ( ( h1*h1 - h3*h3 ) * phi1 - 2._wi*h1*h1*h3 * phi2 + phi3 ) &
                    / ( 6._wi * h1*h1 )

            elseif ( h2*h2 + h3*h3 .lt. i0Z2 * hh ) then
                phi2 = i0_phi2( h4, h1, pAbs, hh, rad )
                phi3 = 0.5_wi / rad1 * log( ( rad1 + rad )**2 / hh1 )
                res = ( ( h1*h1 - 3._wi*h4*h4 ) * phi1 &
                      - h4 * ( 3._wi*h1*h1 - h4*h4 ) * phi2 &
                      + 3._wi*h4*h4 * phi3 - h4*h4 / ( rad + h4 ) ) &
                    / ( 6._wi * h1*h1 )

            elseif ( h1*h1 + h4*h4 .lt. i0Z2 * hh ) then
                hNorm = sqrt( hh )
                rad2 = sqrt( pAbs*pAbs + h2*h2 )
                phi4 = h3*h3 / ( h2 * pAbs*pAbs ) &
                     * ( ( rad2 / h2 ) * log( ( rad2 + rad ) / h3 ) - log( ( h2 + hNorm ) / h3 ) )
                res = ( hh / ( h2*h2 ) * phi1 - 1._wi / ( rad + hNorm ) - phi4 ) / 6._wi

            elseif ( h1*h1 + h3*h3 .lt. i0Z2 * hh ) then
                hNorm = sqrt( hh )
                rad2 = sqrt( pAbs*pAbs + h2*h2 )
                phi2 = atan( h2 * pAbs / ( hh + h4 * rad ) ) / pAbs
                phi4 = h4*h4 / ( h2 * pAbs*pAbs ) &
                     * ( ( rad2 / h2 ) * log( ( rad2 + rad ) / h4 ) - log( ( h2 + hNorm ) / h4 ) )
                res = ( ( 1._wi + 3._wi*h4*h4 / ( h2*h2 ) ) * phi1 &
                      - 2._wi * ( h4 / h2 )**3 * phi2 - 3._wi * phi4 &
                      + ( 2._wi*h4*h4 - h2*h2 ) / ( h2*h2 * ( rad + hNorm ) ) ) / 6._wi

            else
                ok = .false.
                res = ieee_value( 1._wi, ieee_quiet_nan )
            endif
        endif

        ans = real( res, wk )

    end subroutine triangle_i0
