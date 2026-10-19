!---------------------------------------------------------------------------
!> @brief
!> Analytically exact, target-volume-averaged demagnetization tensor between a
!> uniformly magnetized SOURCE tetrahedron and a TARGET tetrahedron.
!>
!> This is the finite-target counterpart of getN_tensor_tetrahedron in
!> TileNComponents, which evaluates the field of a tetrahedral source at a
!> single POINT. Both routines are kept: the point routine is unchanged and
!> keeps its meaning, this one adds the tetrahedral-target average.
!>
!> The returned tensor satisfies the usual MagTense convention
!>
!>     < H_i >_{V_t} = N_ij M_j
!>
!> i.e. it maps the MAGNETIZATION of the source to the field averaged over the
!> volume of the target. Writing the field of the uniformly magnetized source
!> through its surface magnetic charge and applying the divergence theorem over
!> the target volume turns the volume average into a double surface integral,
!>
!>     N_ij = -1/(4 pi V_t) sum_{a=1..4} sum_{b=1..4} n^t_{a,i} n^s_{b,j} I_ab
!>
!>     I_ab = \int_{S^t_a} \int_{S^s_b} 1/|r_t - r_s| dS_s dS_t
!>
!> with n^t and n^s the OUTWARD face normals of target and source. The sixteen
!> scalar face-pair integrals are evaluated by TileTriangleTriangle (closed
!> form in quadruple precision, with numerical fallbacks for nearly degenerate
!> face pairs) for touching, self and near pairs. Pairs separated by at least
!> 1.5 summed circumradii are instead averaged by a collapsed Gauss rule of the
!> analytical point field over the target (see gaussFactor), where the closed
!> form loses its accuracy through cancellation.
!>
!> NOTE ON NORMALISATION. The equivalent routine in dip-fmm maps the TOTAL
!> SOURCE MOMENT m = V_s M to the field and therefore carries the prefactor
!> -1/(4 pi V_s V_t). MagTense maps the magnetization M directly, so the
!> prefactor here is -1/(4 pi V_t) and the division by the SOURCE volume must
!> NOT be applied.
!---------------------------------------------------------------------------
module TileTetrahedronPairTensor
    use TileTriangleTriangle
    use, intrinsic :: ieee_arithmetic
    implicit none

    private

    real(ttReal),parameter :: fourPi = 12.566370614359172953850573533118_ttReal

    !> Separation, in units of the summed circumradii, from which the pair is
    !> averaged by a collapsed Gauss rule of the exact point field of the
    !> source over the target instead of formed from the analytical face-pair
    !> integrals. The analytical reduction cancels catastrophically for
    !> separated pairs: on random irregular tetrahedra its error grows from
    !> 3e-11 at 1.5 to 2e-8 at eight, and for parallel faces with non-aligned
    !> edges (the top and bottom faces of a tetrahedral film mesh) it lost all
    !> digits by ten. Below 1.5 the source's singularity is too close for a
    !> fixed rule. Values and rule orders are dip-fmm's (tetrahedron_detail.hpp):
    !> 7 points per direction from 1.5, 6 from 2 and 5 from 4, each with a worst
    !> relative error of about 2e-10 in its range.
    real(ttReal),parameter :: gaussFactor = 1.5_ttReal
    real(ttReal),parameter :: gaussSixFactor = 2.0_ttReal
    real(ttReal),parameter :: gaussFiveFactor = 4.0_ttReal

    !> The four triangular faces of a tetrahedron, as vertex triples. Face k is
    !> the face opposite vertex k, so the vertex used to orient the outward
    !> normal of face k is vertex k itself.
    integer,dimension(3,4),parameter :: tetFaceVertex = reshape( &
        [ 2, 3, 4,   1, 4, 3,   1, 2, 4,   1, 3, 2 ], [ 3, 4 ] )

    !> Error codes returned through the optional ierr argument.
    integer,parameter,public :: tetPairOk = 0
    integer,parameter,public :: tetPairBadSource = 1
    integer,parameter,public :: tetPairBadTarget = 2
    integer,parameter,public :: tetPairNotFinite = 3

    public :: getN_tensor_tetrahedron_tetrahedron
    public :: tetrahedron_pair_volume

    contains

    !---------------------------------------------------------------------------
    !> @brief Exact tetrahedron-source to tetrahedron-target averaged tensor.
    !>
    !> @param[in] vert_source (3,4) vertices of the source tetrahedron, in
    !>            absolute coordinates, vert_source(:,k) being the k'th vertex
    !> @param[in] vert_target (3,4) vertices of the target tetrahedron, same
    !>            layout and the same coordinate system as vert_source
    !> @param[out] N the (3,3) tensor with < H >_{V_t} = N M
    !> @param[out] ierr optional status, tetPairOk on success
    !>
    !> The result is independent of the orientation and of any permutation of
    !> the vertices of either tetrahedron, because the face normals are made
    !> outward explicitly rather than inherited from the vertex ordering, and it
    !> is invariant under a common translation of the two tetrahedra.
    !---------------------------------------------------------------------------
    subroutine getN_tensor_tetrahedron_tetrahedron( vert_source, vert_target, N, ierr )
    real(ttReal),dimension(3,4),intent(in) :: vert_source, vert_target
    real(ttReal),dimension(3,3),intent(out) :: N
    integer,intent(out),optional :: ierr

    real(ttReal),dimension(3,4) :: vs, vt
    real(ttReal),dimension(3,3,4) :: faceS, faceT
    real(ttReal),dimension(3,4) :: normalS, normalT
    real(ttReal),dimension(3,3) :: acc
    real(ttReal),dimension(3) :: shift
    real(ttReal) :: volS, volT, integral, prefactor, avg, dr
    integer :: a, b, i, j, status, q
    logical :: ok

        N = 0._ttReal
        status = tetPairOk

        !Work relative to the source centroid. This is an exact translation of
        !the geometry and keeps all coordinates of the order of the element size
        !even for a mesh sitting far from the origin, where the difference of
        !two large absolute coordinates would otherwise lose precision.
        shift = sum( vert_source, dim=2 ) / 4._ttReal
        do i = 1, 4
            vs(:,i) = vert_source(:,i) - shift
            vt(:,i) = vert_target(:,i) - shift
        enddo

        call tetrahedron_pair_volume( vs, volS, ok )
        if ( .not. ok ) then
            if ( present(ierr) ) ierr = tetPairBadSource
            return
        endif

        call tetrahedron_pair_volume( vt, volT, ok )
        if ( .not. ok ) then
            if ( present(ierr) ) ierr = tetPairBadTarget
            return
        endif

        !Separated pairs: Gauss rule of the analytical point field (see
        !gaussFactor). The source centroid is the origin here.
        dr = norm2( sum( vt, dim=2 ) / 4._ttReal ) &
           / ( tetrahedron_circumradius( vs ) + tetrahedron_circumradius( vt ) )
        if ( dr .ge. gaussFactor .and. .not. tetrahedra_share_vertex( vs, vt ) ) then
            q = 7
            if ( dr .ge. gaussSixFactor ) q = 6
            if ( dr .ge. gaussFiveFactor ) q = 5
            call gauss_pair_tensor( vs, vt, q, N )
            if ( .not. all( ieee_is_finite( N ) ) ) then
                N = 0._ttReal
                if ( present(ierr) ) ierr = tetPairNotFinite
                return
            endif
            do i = 1, 3
                do j = i+1, 3
                    avg = 0.5_ttReal * ( N(i,j) + N(j,i) )
                    N(i,j) = avg
                    N(j,i) = avg
                enddo
            enddo
            if ( present(ierr) ) ierr = tetPairOk
            return
        endif

        call tetrahedron_faces( vs, faceS, normalS, ok )
        if ( .not. ok ) then
            if ( present(ierr) ) ierr = tetPairBadSource
            return
        endif

        call tetrahedron_faces( vt, faceT, normalT, ok )
        if ( .not. ok ) then
            if ( present(ierr) ) ierr = tetPairBadTarget
            return
        endif

        !The sixteen analytical face-pair integrals.
        acc = 0._ttReal
        do a = 1, 4
            do b = 1, 4
                call triangle_triangle_laplace_integral( faceT(:,:,a), faceS(:,:,b), integral, ok )
                if ( .not. ok ) then
                    if ( present(ierr) ) ierr = tetPairNotFinite
                    N = 0._ttReal
                    return
                endif
                do i = 1, 3
                    do j = 1, 3
                        acc(i,j) = acc(i,j) + integral * normalT(i,a) * normalS(j,b)
                    enddo
                enddo
            enddo
        enddo

        !MagTense normalisation: magnetization to field, averaged over the
        !TARGET volume only. No division by the source volume.
        prefactor = -1._ttReal / ( fourPi * volT )
        N = prefactor * acc

        if ( .not. all( ieee_is_finite( N ) ) ) then
            N = 0._ttReal
            if ( present(ierr) ) ierr = tetPairNotFinite
            return
        endif

        !The tensor is symmetric analytically. The sixteen face integrals are
        !evaluated independently, so the off-diagonal pairs can differ in the
        !last bits; averaging them removes that round-off and nothing else.
        do i = 1, 3
            do j = i+1, 3
                avg = 0.5_ttReal * ( N(i,j) + N(j,i) )
                N(i,j) = avg
                N(j,i) = avg
            enddo
        enddo

        if ( present(ierr) ) ierr = status

    end subroutine getN_tensor_tetrahedron_tetrahedron


    !---------------------------------------------------------------------------
    !> @brief Positive volume of a tetrahedron, with a validity check.
    !>
    !> @param[out] vol the positive volume
    !> @param[out] ok .false. for a non-finite or degenerate (rank deficient)
    !>             tetrahedron. The tolerance is scale aware so that valid
    !>             geometry is accepted over a wide range of physical scales
    !>             while a round-off sized volume is rejected.
    !---------------------------------------------------------------------------
    subroutine tetrahedron_pair_volume( vert, vol, ok )
    real(ttReal),dimension(3,4),intent(in) :: vert
    real(ttReal),intent(out) :: vol
    logical,intent(out) :: ok

    real(ttReal),dimension(3) :: a, b, c
    real(ttReal) :: scaleVal, tol
    integer :: i, j

        vol = 0._ttReal
        ok = .false.

        if ( .not. all( ieee_is_finite( vert ) ) ) return

        scaleVal = 0._ttReal
        do i = 1, 4
            do j = i+1, 4
                scaleVal = max( scaleVal, norm2( vert(:,i) - vert(:,j) ) )
            enddo
        enddo
        if ( .not. ( scaleVal .gt. 0._ttReal ) ) return

        a = vert(:,2) - vert(:,1)
        b = vert(:,3) - vert(:,1)
        c = vert(:,4) - vert(:,1)

        vol = abs( a(1) * ( b(2)*c(3) - b(3)*c(2) ) &
                 + a(2) * ( b(3)*c(1) - b(1)*c(3) ) &
                 + a(3) * ( b(1)*c(2) - b(2)*c(1) ) ) / 6._ttReal

        tol = 128._ttReal * epsilon(1._ttReal) * scaleVal**3
        if ( .not. ( vol .gt. tol ) .or. .not. ieee_is_finite( vol ) ) then
            vol = 0._ttReal
            return
        endif

        ok = .true.

    end subroutine tetrahedron_pair_volume


    !---------------------------------------------------------------------------
    !> @brief The four faces of a tetrahedron with their OUTWARD unit normals.
    !>
    !> The normal of a face is flipped, if needed, so that it points away from
    !> the remaining vertex. That makes the result independent of the handedness
    !> of the supplied vertex ordering.
    !---------------------------------------------------------------------------
    subroutine tetrahedron_faces( vert, faces, normals, ok )
    real(ttReal),dimension(3,4),intent(in) :: vert
    real(ttReal),dimension(3,3,4),intent(out) :: faces
    real(ttReal),dimension(3,4),intent(out) :: normals
    logical,intent(out) :: ok

    real(ttReal),dimension(3) :: p1, p2, p3, nrm
    real(ttReal) :: len
    integer :: f, k

        faces = 0._ttReal
        normals = 0._ttReal
        ok = .false.

        do f = 1, 4
            do k = 1, 3
                faces(:,k,f) = vert(:,tetFaceVertex(k,f))
            enddo

            p1 = faces(:,1,f)
            p2 = faces(:,2,f)
            p3 = faces(:,3,f)

            nrm(1) = ( p2(2)-p1(2) ) * ( p3(3)-p1(3) ) - ( p2(3)-p1(3) ) * ( p3(2)-p1(2) )
            nrm(2) = ( p2(3)-p1(3) ) * ( p3(1)-p1(1) ) - ( p2(1)-p1(1) ) * ( p3(3)-p1(3) )
            nrm(3) = ( p2(1)-p1(1) ) * ( p3(2)-p1(2) ) - ( p2(2)-p1(2) ) * ( p3(1)-p1(1) )

            !Face f is opposite vertex f, so orient away from that vertex.
            if ( dot_product( nrm, vert(:,f) - p1 ) .gt. 0._ttReal ) nrm = -nrm

            len = norm2( nrm )
            if ( .not. ( len .gt. 0._ttReal ) .or. .not. ieee_is_finite( len ) ) return

            normals(:,f) = nrm / len
        enddo

        ok = .true.

    end subroutine tetrahedron_faces


    !> Circumradius of a non-degenerate tetrahedron.
    function tetrahedron_circumradius( vert ) result( radius )
    real(ttReal),dimension(3,4),intent(in) :: vert
    real(ttReal) :: radius
    real(ttReal),dimension(3) :: u, v, w, centre
        u = vert(:,2) - vert(:,1)
        v = vert(:,3) - vert(:,1)
        w = vert(:,4) - vert(:,1)
        centre = ( dot_product( u, u ) * cross_product( v, w ) + dot_product( v, v ) * cross_product( w, u ) &
                 + dot_product( w, w ) * cross_product( u, v ) ) / ( 2._ttReal * dot_product( u, cross_product( v, w ) ) )
        radius = norm2( centre )
    end function tetrahedron_circumradius


    function cross_product( a, b ) result( c )
    real(ttReal),dimension(3),intent(in) :: a, b
    real(ttReal),dimension(3) :: c
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product


    !> .true. if the two tetrahedra have a vertex in common (to round-off of the
    !> element size): their point field is singular on the target, so such a
    !> pair is never averaged by quadrature.
    function tetrahedra_share_vertex( vs, vt ) result( shared )
    real(ttReal),dimension(3,4),intent(in) :: vs, vt
    logical :: shared
    real(ttReal) :: scaleVal
    integer :: i, j
        scaleVal = 0._ttReal
        do i = 1, 4
            do j = i+1, 4
                scaleVal = max( scaleVal, norm2( vs(:,i) - vs(:,j) ), norm2( vt(:,i) - vt(:,j) ) )
            enddo
        enddo
        shared = .false.
        do i = 1, 4
            do j = 1, 4
                if ( norm2( vs(:,i) - vt(:,j) ) .le. 1.0e-10_ttReal * scaleVal ) shared = .true.
            enddo
        enddo
    end function tetrahedra_share_vertex


    !---------------------------------------------------------------------------
    !> @brief Target average of the analytical point tensor of the source by a
    !> collapsed (Duffy) Gauss rule of q points per direction, q**3 in all.
    !>
    !> The point tensor is tet_point_tensor, with the same sign convention as
    !> the analytical pair tensor. MagTense's getN_Triangle is not used here: its
    !> differences of atanh and atan terms lose up to 1e-3 relative accuracy at
    !> twenty element sizes (validation of 19 Oct 2026, FMM-article repository,
    !> tet-tet-article/data/REPORT.md "A1").
    !---------------------------------------------------------------------------
    subroutine gauss_pair_tensor( vs, vt, q, N )
    real(ttReal),dimension(3,4),intent(in) :: vs, vt
    integer,intent(in) :: q
    real(ttReal),dimension(3,3),intent(out) :: N

    real(ttReal),dimension(q) :: x, w
    real(ttReal),dimension(4) :: lam
    real(ttReal),dimension(3,3,4) :: faces
    real(ttReal),dimension(3,4) :: normals
    real(ttReal),dimension(3,3) :: nPt
    real(ttReal),dimension(3) :: pt
    real(ttReal) :: r, s, t, wt, wsum
    integer :: a, b, c
    logical :: ok

        N = 0._ttReal
        call tetrahedron_faces( vs, faces, normals, ok )
        if ( .not. ok ) then
            N = ieee_value( 1._ttReal, ieee_quiet_nan )
            return
        endif
        call gauss_legendre01( q, x, w )
        wsum = 0._ttReal
        do a = 1, q
            r = x(a)
            do b = 1, q
                s = x(b)
                do c = 1, q
                    t = x(c)
                    wt = 6._ttReal * r*r * s * w(a) * w(b) * w(c)
                    lam = [ 1._ttReal - r, r * ( 1._ttReal - s ), r * s * ( 1._ttReal - t ), r * s * t ]
                    pt = matmul( vt, lam )
                    call tet_point_tensor( faces, normals, pt, nPt )
                    N = N + wt * nPt
                    wsum = wsum + wt
                enddo
            enddo
        enddo
        N = N / wsum

    end subroutine gauss_pair_tensor


    !---------------------------------------------------------------------------
    !> @brief Field tensor of a uniformly magnetised tetrahedron at a point
    !> outside it, H = N M, in the face-charge form
    !>
    !>     N = 1/(4 pi) sum_f ( Omega_f n_f + sum_e L_e m_e ) n_f^T ,
    !>
    !> Omega_f the solid angle of face f signed with its outward normal n_f
    !> (Van Oosterom-Strackee), m_e the outward in-plane normal of edge e and
    !> L_e = ln( (R1+R2+l)/(R1+R2-l) ). Both are evaluated without cancellation:
    !> L_e as 2 atanh(l/(R1+R2)) away from the edge, and with R1+R2-l formed as
    !> ((R1+R2)**2 - l**2)/(R1+R2+l) near it.
    !---------------------------------------------------------------------------
    subroutine tet_point_tensor( faces, normals, x, N )
    real(ttReal),dimension(3,3,4),intent(in) :: faces
    real(ttReal),dimension(3,4),intent(in) :: normals
    real(ttReal),dimension(3),intent(in) :: x
    real(ttReal),dimension(3,3),intent(out) :: N

    real(ttReal),dimension(3,3) :: p
    real(ttReal),dimension(3) :: nrm, a, b, c, vec, e, m, tmp, rr
    real(ttReal) :: det, den, le, ssum, lg
    integer :: f, i, j, k

        N = 0._ttReal
        do f = 1, 4
            p = faces(:,:,f)
            nrm = normals(:,f)
            !Vertices counter-clockwise about the outward normal.
            if ( dot_product( cross_product( p(:,2) - p(:,1), p(:,3) - p(:,1) ), nrm ) .lt. 0._ttReal ) then
                tmp = p(:,2)
                p(:,2) = p(:,3)
                p(:,3) = tmp
            endif
            a = p(:,1) - x
            b = p(:,2) - x
            c = p(:,3) - x
            rr = [ norm2( a ), norm2( b ), norm2( c ) ]
            det = dot_product( a, cross_product( p(:,2) - p(:,1), p(:,3) - p(:,1) ) )
            den = rr(1) * rr(2) * rr(3) + dot_product( a, b ) * rr(3) &
                + dot_product( a, c ) * rr(2) + dot_product( b, c ) * rr(1)
            vec = -2._ttReal * atan2( det, den ) * nrm
            do i = 1, 3
                j = modulo( i, 3 ) + 1
                e = p(:,j) - p(:,i)
                le = norm2( e )
                m = cross_product( e, nrm ) / le
                ssum = rr(i) + rr(j)
                if ( le .lt. 0.5_ttReal * ssum ) then
                    lg = 2._ttReal * atanh( le / ssum )
                else
                    lg = log( ( ssum + le ) / ( ( ssum*ssum - le*le ) / ( ssum + le ) ) )
                endif
                vec = vec + lg * m
            enddo
            do k = 1, 3
                N(:,k) = N(:,k) + vec * nrm(k)
            enddo
        enddo
        N = N / fourPi

    end subroutine tet_point_tensor

end module TileTetrahedronPairTensor
