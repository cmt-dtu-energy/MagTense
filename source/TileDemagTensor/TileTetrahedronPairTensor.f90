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
!> scalar face-pair integrals are evaluated in closed form by
!> TileTriangleTriangle; there is no numerical quadrature anywhere in this path.
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
    real(ttReal) :: volS, volT, integral, prefactor, avg
    integer :: a, b, i, j, status
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

end module TileTetrahedronPairTensor
