!>-----------------------------------------
!> Cell-centred exchange operator for unstructured grids (tetrahedra and prisms).
!>
!> The exchange operator div( A grad(m) ) at the centroid of cell i is approximated as
!>
!>     L m_i = sum_j c_ij A_ij ( m_j - m_i ),
!>
!> where the c_ij come from a weighted least squares fit of
!>
!>     phi_j - phi_i = g.d_ij + 1/2 d_ij^T H d_ij        (9 unknowns: g and the symmetric H)
!>
!> over the stencil of i, and sum_j c_ij phi_j is the trace of the fitted Hessian H. The fit is
!> exact for quadratic fields on any stencil that determines it, so the operator is consistent on
!> irregular meshes, unlike a face-based Direct Laplacian whose face gradients are only exact on
!> symmetric (axis aligned, uniform) stencils.
!>
!> Neighbours.
!>   Tetrahedra: every cell that shares a node with i. Nodes are identified by their coordinates
!>   (tolerance 1e-6 of the shortest edge; the neighbouring keys are searched too, so rounding of
!>   the key cannot split two coinciding nodes), so duplicated nodes and the nodes of the two
!>   planes of a periodic direction are merged.
!>   Prisms: every cell whose closed box touches the box of i (faces, edges, corners, also across
!>   hanging faces of a refined grid and for arbitrary cell widths), found with a bucket grid.
!>   A stencil with fewer than minStencil cells is grown by the neighbours of its neighbours.
!>   The offset d_ij of every neighbour is the actual one through the shared node (tetrahedra) or
!>   the touching periodic image (prisms), never a minimum image, so it is right for any period.
!>
!> Boundaries.
!>   Free surfaces (homogeneous Neumann, n.grad(m) = 0) are imposed in the fit itself: every
!>   boundary face of i or of a stencil cell adds the row n.( g + H d_f ) = 0, with d_f the offset
!>   of the face centre. This holds for any boundary shape and at edges and corners, and needs no
!>   ghost cells or mirror planes. A boundary face is a tetrahedron face owned by one live cell,
!>   or the part of a prism face not covered by live neighbours. Faces towards cells with A0 = 0
!>   (non-magnetic) are free surfaces too.
!>   Periodic coupling is used only along the directions with exchPBC set (period = extent of the
!>   mesh): there the two planes are matched, so they are not boundaries. Without exchPBC every
!>   outer face is a free surface.
!>
!> Exchange stiffness: A_ij is the interface value A_int(phase_i,phase_j) when given for a pair of
!> different phases, else the harmonic mean of A_i and A_j; cells with A0 = 0 carry no exchange and
!> are left out of all stencils.
!>
!> Symmetric operator SUSHI-g (computeExchangeSushi, exch_meth 'sushi').
!> The cell Hessian operator is consistent but not symmetric in the volume weighted inner product,
!> so its field is not the gradient of an energy: undamped dynamics do not conserve the exchange
!> energy and damped dynamics need not decrease it. SUSHI-g derives the operator from a discrete
!> energy E(u) = u^T K u, K symmetric positive semidefinite with constants in its kernel, and
!> L = -V^-1 K, so V L is symmetric, the field is exactly -dE/dm and the energy is conserved
!> (alpha = 0) or decreases (alpha > 0). The energy is the SUSHI scheme (Eymard, Gallouet, Herbin,
!> IMA J. Numer. Anal. 30 (2010) 1009) with cell and face values:
!>   g_K  = 1/|K| sum_s |s| (u_s - u_K) n_Ks,  R_Ks = u_s - u_K - g_K.(x_s - x_K),
!>   E_K  = A_K sum_s |s| d_Ks/3 | g_K + n_Ks R_Ks stab / max(d_Ks, 0.3 h_K) |^2,
!> exact (R = 0, g = grad u) for affine u. Interior face values are interpolated from cell values,
!> exact for affine fields and well conditioned (two-point when the face point lies between the
!> centroids, else the face neighbours of both cells, else their vertex neighbours; least squares
!> weights V_j^(1/2) |d_j|^-4, so that small cells carry small weights); free-surface
!> face values are eliminated by minimising E_K, which is the natural homogeneous Neumann
!> condition. stab = sqrt(3) for prisms, where the operator is then identical to the 7-point
!> Laplacian with mirror Neumann boundaries on uniform grids, and 1 for tetrahedra. Every cell's
!> energy carries its own A_K; interface exchange values are not used by this operator.
!---------------------------------------------------------------------------
module CellHessianExchange
    use MKL_SPBLAS
    use MicroMagParameters
    use IO_GENERAL
    use UTIL_MICROMAG
    implicit none
    private

    integer, parameter :: i8 = selected_int_kind(18)
    integer, parameter :: minStencil = 12
    integer, parameter :: nUnk = 9
    real(dp), parameter :: weightExp = 2._dp
    real(dp), parameter :: rcondFit = 1.e-10_dp

    !> SUSHI-g: stabilisation (boxes: identical to the 7-point Laplacian on uniform grids;
    !> tetrahedra: from the parameter study), its distance floor (fraction of the cell size), the
    !> largest face-neighbour stencil used whole and the conditioning bound of the interpolation
    real(dp), parameter :: stabTet = 1._dp
    real(dp), parameter :: stabBox = 1.7320508075688772_dp
    real(dp), parameter :: capStab = 0.3_dp
    real(dp), parameter :: bMax = 2.5_dp
    integer, parameter :: maxFull = 8

    !> Open addressing hash table for triples of 64 bit integer keys
    type keyTable
        integer :: nslot = 0
        integer(i8), allocatable :: key(:,:)
        integer, allocatable :: val(:)
        logical, allocatable :: used(:)
    end type keyTable

    !> Interpolation of a face value from cell values
    type faceStencil
        integer, allocatable :: idx(:)
        real(dp), allocatable :: w(:)
    end type faceStencil

    public :: computeExchangeCellHessian, computeExchangeSushi

    contains

    !>-----------------------------------------
    !> @brief Builds the exchange matrix of an unstructured grid with the cell Hessian operator.
    !> @param[in] grid the grid (tetrahedra: nodes, elements; prisms: pts, abc)
    !> @param[in] pbc periodic exchange along x, y, z
    !> @param[in] A0n exchange stiffness per cell, normalised as for the Direct Laplacian
    !> @param[in] phase material index per cell
    !> @param[in] Aint interface exchange per phase pair, same normalisation, negative = unset
    !> @param[inout] GridInfo receives the COO copy of the matrix
    !> @param[out] Aexch_matrix the CSR exchange matrix
    !---------------------------------------------------------------------------
    subroutine computeExchangeCellHessian( grid, pbc, A0n, phase, Aint, GridInfo, Aexch_matrix )
    type(MicroMagGrid), intent(in) :: grid
    logical, intent(in) :: pbc(3)
    real(dp), intent(in) :: A0n(:)
    integer, intent(in) :: phase(:)
    real(dp), intent(in) :: Aint(:,:)
    type(MicroMagGridInfo), intent(inout) :: GridInfo
    type(sparse_matrix_t), intent(out) :: Aexch_matrix

    real(dp), allocatable :: cen(:,:), bfC(:,:), bfN(:,:), nbOff(:,:), stOff(:,:), cooV(:)
    integer, allocatable :: nbPtr(:), nbIdx(:), bfOwner(:), bfPtr(:), bfIdx(:)
    integer, allocatable :: stPtr(:), stIdx(:), cooR(:), cooC(:)
    logical, allocatable :: live(:)
    real(dp) :: per(3), lo(3), hi(3)
    integer :: N, i, nGrown, nRankDef, nBad, stMin, stMax, nbf
    character(len=256) :: msg

        if ( grid%gridType .eq. gridTypeTetrahedron ) then
            N = size(grid%elements, 2)
        else
            N = size(grid%pts, 1)
        endif
        allocate( live(N) )
        live = A0n(1:N) .gt. 0._dp

        !--- centroids, extent and periods
        allocate( cen(3,N) )
        if ( grid%gridType .eq. gridTypeTetrahedron ) then
            do i = 1, N
                cen(:,i) = sum( grid%nodes(:, grid%elements(1:4,i)), dim=2 ) / 4._dp
            enddo
            lo = minval( grid%nodes, dim=2 )
            hi = maxval( grid%nodes, dim=2 )
        else
            cen = transpose( grid%pts(1:N,:) )
            lo = minval( grid%pts(1:N,:) - 0.5_dp * grid%abc(1:N,:), dim=1 )
            hi = maxval( grid%pts(1:N,:) + 0.5_dp * grid%abc(1:N,:), dim=1 )
        endif
        per = 0._dp
        where ( pbc ) per = hi - lo

        !--- neighbours (with their actual offsets) and boundary faces
        if ( grid%gridType .eq. gridTypeTetrahedron ) then
            call tetNeighboursAndFaces( grid%nodes, grid%elements(1:4,1:N), cen, live, pbc, per, lo, &
                                        nbPtr, nbIdx, nbOff, bfOwner, bfC, bfN, nbf, nBad )
            if ( nBad .gt. 0 ) then
                write(msg,'(A,I0,A)') 'Cell Hessian exchange: ', nBad, &
                    ' faces are shared by more than two tetrahedra - the mesh is not conforming'
                call displayGUIMessage( trim(msg) )
            endif
        else
            call prismNeighboursAndFaces( grid%pts(1:N,:), grid%abc(1:N,:), live, per, lo, hi, &
                                          nbPtr, nbIdx, nbOff, bfOwner, bfC, bfN, nbf )
        endif

        !--- stencils: live neighbours, grown to the second ring when too small
        call buildStencils( N, live, nbPtr, nbIdx, nbOff, stPtr, stIdx, stOff, nGrown, stMin, stMax )

        !--- boundary faces grouped by owner
        call groupByOwner( N, bfOwner(1:nbf), bfPtr, bfIdx )

        !--- the fits
        call assemble( N, cen, live, A0n, phase, Aint, stPtr, stIdx, stOff, bfPtr, bfIdx, bfC, bfN, &
                       cooR, cooC, cooV, nRankDef )

        call create_CSR_matrix( cooR, cooC, cooV, N, N, 0._dp, Aexch_matrix )
        call create_COO_values_from_CSR( Aexch_matrix, GridInfo )

        write(msg,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') 'Cell Hessian exchange: ', N, ' cells, stencil ', &
            stMin, '-', stMax, ' cells, ', nbf, ' boundary faces, ', nGrown, &
            ' stencils grown, rank deficient fits: ', nRankDef
        call displayGUIMessage( trim(msg) )
        if ( any(pbc) ) call displayGUIMessage( 'Cell Hessian exchange: periodic along the exchPBC directions only' )

    end subroutine computeExchangeCellHessian


    !------------------------------------------------------------------ growable buffers
    subroutine pushInt( a, n, v )
    integer, allocatable, intent(inout) :: a(:)
    integer, intent(inout) :: n
    integer, intent(in) :: v
    integer, allocatable :: t(:)
        if ( .not. allocated(a) ) allocate( a(1024) )
        if ( n .ge. size(a) ) then
            allocate( t(2*size(a)) )
            t(1:n) = a(1:n)
            call move_alloc( t, a )
        endif
        n = n + 1
        a(n) = v
    end subroutine pushInt

    subroutine pushVec( a, n, v )
    real(dp), allocatable, intent(inout) :: a(:,:)
    integer, intent(in) :: n          !> the slot to write (the matching pushInt has advanced to it)
    real(dp), intent(in) :: v(3)
    real(dp), allocatable :: t(:,:)
        if ( .not. allocated(a) ) allocate( a(3,1024) )
        if ( n .gt. size(a,2) ) then
            allocate( t(3, 2*size(a,2)) )
            t(:,1:size(a,2)) = a
            call move_alloc( t, a )
        endif
        a(:,n) = v
    end subroutine pushVec


    !>-----------------------------------------
    !> Tetrahedra: canonical nodes, node-sharing neighbours with offsets, boundary faces
    !---------------------------------------------------------------------------
    subroutine tetNeighboursAndFaces( nodes, ele, cen, live, pbc, per, lo, nbPtr, nbIdx, nbOff, &
                                      bfOwner, bfC, bfN, nbf, nBad )
    real(dp), intent(in) :: nodes(:,:), cen(:,:), per(3), lo(3)
    integer, intent(in) :: ele(:,:)
    logical, intent(in) :: live(:), pbc(3)
    integer, allocatable, intent(out) :: nbPtr(:), nbIdx(:), bfOwner(:)
    real(dp), allocatable, intent(out) :: nbOff(:,:), bfC(:,:), bfN(:,:)
    integer, intent(out) :: nbf, nBad

    integer :: N, i, k, a, b, j, v, c, slot, nCanon, f, nf, owner, nnz, kk
    integer, allocatable :: ndCnt(:), ndPtr(:), ndItem(:), stamp(:)
    integer, allocatable :: fCount(:), fOwner(:), fLoc(:), fOther(:), Ec(:,:)
    real(dp) :: p(3,3), nrm(3), fc(3), d(3)
    integer(i8) :: q(3)
    type(keyTable) :: ftab
    integer :: tri(3)

        N = size(ele, 2)
        call tetCanonicalNodes( nodes, ele, pbc, per, lo, Ec, nCanon )

        !Cells of each canonical node
        allocate( ndCnt(nCanon), ndPtr(nCanon+1) )
        ndCnt = 0
        do i = 1, N
            do a = 1, 4
                ndCnt(Ec(a,i)) = ndCnt(Ec(a,i)) + 1
            enddo
        enddo
        ndPtr(1) = 1
        do v = 1, nCanon
            ndPtr(v+1) = ndPtr(v) + ndCnt(v)
        enddo
        allocate( ndItem(ndPtr(nCanon+1)-1) )
        ndCnt = 0
        do i = 1, N
            do a = 1, 4
                v = Ec(a,i)
                ndItem(ndPtr(v) + ndCnt(v)) = i
                ndCnt(v) = ndCnt(v) + 1
            enddo
        enddo

        !Node-sharing neighbours, deduplicated with a stamp. The offset goes through the shared
        !node, using its coordinates as each of the two cells has it: right for any period.
        allocate( nbPtr(N+1), stamp(N) )
        stamp = 0
        nbPtr(1) = 1
        nnz = 0
        do i = 1, N
            stamp(i) = i
            do a = 1, 4
                v = Ec(a,i)
                do k = ndPtr(v), ndPtr(v+1)-1
                    j = ndItem(k)
                    if ( stamp(j) .eq. i ) cycle
                    stamp(j) = i
                    do b = 1, 4
                        if ( Ec(b,j) .eq. v ) exit
                    enddo
                    d = ( cen(:,j) - nodes(:,ele(b,j)) ) - ( cen(:,i) - nodes(:,ele(a,i)) )
                    call pushInt( nbIdx, nnz, j )
                    call pushVec( nbOff, nnz, d )
                enddo
            enddo
            nbPtr(i+1) = nnz + 1
        enddo
        if ( .not. allocated(nbIdx) ) allocate( nbIdx(1), nbOff(3,1) )

        !Faces by sorted canonical node triple; owned by one live cell, or between a live and a
        !dead cell, means a free surface of the live one
        nf = 4 * N
        allocate( fCount(nf), fOwner(nf), fLoc(nf), fOther(nf) )
        fCount = 0
        call tableInit( ftab, nf )
        f = 0
        do i = 1, N
            do k = 1, 4
                tri = pack( Ec(:,i), [ (a .ne. k, a = 1, 4) ] )
                call sort3( tri )
                q = int( tri, i8 )
                slot = tableFind( ftab, q )
                if ( slot .gt. 0 ) then
                    j = ftab%val(slot)
                    fCount(j) = fCount(j) + 1
                    fOther(j) = i
                else
                    f = f + 1
                    call tableInsert( ftab, q, f )
                    fCount(f) = 1
                    fOwner(f) = i
                    fLoc(f) = k
                    fOther(f) = 0
                endif
            enddo
        enddo
        nBad = count( fCount(1:f) .gt. 2 )

        nbf = 0
        allocate( bfOwner(1024), bfC(3,1024), bfN(3,1024) )
        do j = 1, f
            owner = 0
            if ( fCount(j) .eq. 1 ) then
                if ( live(fOwner(j)) ) owner = fOwner(j)
            elseif ( fCount(j) .eq. 2 ) then
                if ( live(fOwner(j)) .and. .not. live(fOther(j)) ) owner = fOwner(j)
                if ( live(fOther(j)) .and. .not. live(fOwner(j)) ) owner = -fOther(j)
            endif
            if ( owner .eq. 0 ) cycle
            if ( owner .gt. 0 ) then
                i = owner
                k = fLoc(j)
            else
                !the second cell: its local face with the same canonical triple
                i = -owner
                do k = 1, 4
                    tri = pack( Ec(:,i), [ (a .ne. k, a = 1, 4) ] )
                    call sort3( tri )
                    slot = tableFind( ftab, int(tri,i8) )
                    if ( slot .gt. 0 ) then
                        if ( ftab%val(slot) .eq. j ) exit
                    endif
                enddo
            endif
            c = 0
            do a = 1, 4
                if ( a .eq. k ) cycle
                c = c + 1
                p(:,c) = nodes(:,ele(a,i))
            enddo
            fc = sum( p, dim=2 ) / 3._dp
            nrm = cross( p(:,2) - p(:,1), p(:,3) - p(:,1) )
            nrm = nrm / norm2( nrm )
            if ( dot_product( nrm, nodes(:,ele(k,i)) - fc ) .gt. 0._dp ) nrm = -nrm
            kk = nbf
            call pushInt( bfOwner, nbf, i )
            call pushVec( bfC, nbf, fc )
            call pushVec( bfN, nbf, nrm )
        enddo

    end subroutine tetNeighboursAndFaces


    !>-----------------------------------------
    !> Tetrahedra: canonical node ids per element corner. Nodes closer than 1e-6 of the shortest
    !> edge are one node; along periodic directions the two boundary planes are identified.
    !---------------------------------------------------------------------------
    subroutine tetCanonicalNodes( nodes, ele, pbc, per, lo, Ec, nCanon )
    real(dp), intent(in) :: nodes(:,:), per(3), lo(3)
    integer, intent(in) :: ele(:,:)
    logical, intent(in) :: pbc(3)
    integer, allocatable, intent(out) :: Ec(:,:)
    integer, intent(out) :: nCanon

    integer :: N, M, i, k, a, b, v, c, slot, da, db, dc, s1, found
    integer, allocatable :: canon(:)
    real(dp) :: tol, emin, X(3), d(3)
    integer(i8) :: q(3), kq(3), Kp(3)
    type(keyTable) :: tab

        N = size(ele, 2)
        M = size(nodes, 2)

        !Tolerance for identifying nodes: a millionth of the shortest edge
        emin = huge(1._dp)
        do i = 1, N
            do a = 1, 3
                do b = a+1, 4
                    emin = min( emin, norm2( nodes(:,ele(a,i)) - nodes(:,ele(b,i)) ) )
                enddo
            enddo
        enddo
        tol = 1.e-6_dp * emin
        Kp = 0
        do k = 1, 3
            if ( pbc(k) ) Kp(k) = nint( per(k) / tol, i8 )
        enddo

        !Canonical node ids. The coordinate (wrapped over the period along periodic directions and
        !shifted by tol/4 so that both periodic planes fall inside key 0) is rounded to a key; the
        !27 neighbouring keys, wrapped along periodic directions, are searched as well.
        allocate( canon(M) )
        call tableInit( tab, M )
        nCanon = 0
        do v = 1, M
            X = nodes(:,v) - lo
            do k = 1, 3
                if ( pbc(k) ) X(k) = modulo( X(k) + 0.25_dp * tol, per(k) )
            enddo
            q = nint( X / tol, i8 )
            found = 0
            do da = -1, 1
                do db = -1, 1
                    do dc = -1, 1
                        if ( found .gt. 0 ) cycle
                        kq = q + [ da, db, dc ]
                        do s1 = 1, 3
                            if ( pbc(s1) ) then
                                if ( kq(s1) .lt. 0 ) kq(s1) = kq(s1) + Kp(s1)
                                if ( kq(s1) .ge. Kp(s1) ) kq(s1) = kq(s1) - Kp(s1)
                            endif
                        enddo
                        slot = tableFind( tab, kq )
                        if ( slot .gt. 0 ) then
                            c = tab%val(slot)
                            d = nodes(:,c) - nodes(:,v)
                            do s1 = 1, 3
                                if ( pbc(s1) ) d(s1) = d(s1) - per(s1) * anint( d(s1) / per(s1) )
                            enddo
                            if ( maxval( abs(d) ) .le. tol ) found = c
                        endif
                    enddo
                enddo
            enddo
            if ( found .gt. 0 ) then
                canon(v) = canon(found)
            else
                nCanon = nCanon + 1
                canon(v) = nCanon
                call tableInsert( tab, q, v )
            endif
        enddo
        allocate( Ec(4,N) )
        do i = 1, N
            Ec(:,i) = canon(ele(:,i))
        enddo
    end subroutine tetCanonicalNodes


    !>-----------------------------------------
    !> Prisms: touching boxes (bucket grid, explicit periodic images) and uncovered face parts
    !---------------------------------------------------------------------------
    subroutine prismNeighboursAndFaces( pts, abc, live, per, lo, hi, nbPtr, nbIdx, nbOff, &
                                        bfOwner, bfC, bfN, nbf )
    real(dp), intent(in) :: pts(:,:), abc(:,:), per(3), lo(3), hi(3)
    logical, intent(in) :: live(:)
    integer, allocatable, intent(out) :: nbPtr(:), nbIdx(:), bfOwner(:)
    real(dp), allocatable, intent(out) :: nbOff(:,:), bfC(:,:), bfN(:,:)
    integer, intent(out) :: nbf

    integer :: N, i, j, k, ax, s, t1, t2, b, bi(3), nbk(3), nnz, first, ii, jj, kk, nbuck
    integer, allocatable :: bucketOf(:), bPtr(:), bItem(:), bCnt(:), stamp(:)
    real(dp) :: bs, tol, d(3), cover, area, ov, lo_, hi_, fc(3), nrm(3)

        N = size(pts, 1)
        tol = 1.e-6_dp * minval( abc )
        bs = maxval( abc )
        nbk = max( 1, int( (hi - lo) / bs ) )
        nbuck = product( nbk )
        allocate( bucketOf(N), bCnt(nbuck), bPtr(nbuck+1), bItem(N), stamp(N) )
        bCnt = 0
        do i = 1, N
            bi = min( nbk - 1, max( 0, int( (pts(i,:) - lo) / bs ) ) )
            bucketOf(i) = 1 + bi(1) + nbk(1) * ( bi(2) + nbk(2) * bi(3) )
            bCnt(bucketOf(i)) = bCnt(bucketOf(i)) + 1
        enddo
        bPtr(1) = 1
        do b = 1, nbuck
            bPtr(b+1) = bPtr(b) + bCnt(b)
        enddo
        bCnt = 0
        do i = 1, N
            b = bucketOf(i)
            bItem(bPtr(b) + bCnt(b)) = i
            bCnt(b) = bCnt(b) + 1
        enddo

        allocate( nbPtr(N+1), bfOwner(1024), bfC(3,1024), bfN(3,1024) )
        nbPtr(1) = 1
        nnz = 0
        nbf = 0
        stamp = 0
        do i = 1, N
            bi = min( nbk - 1, max( 0, int( (pts(i,:) - lo) / bs ) ) )
            first = nnz + 1
            do ii = -1, 1
              do jj = -1, 1
                do kk = -1, 1
                    call visitBucket( bi + [ ii, jj, kk ], [ ii, jj, kk ] )
                enddo
              enddo
            enddo
            nbPtr(i+1) = nnz + 1
            if ( .not. allocated(nbIdx) ) allocate( nbIdx(1024), nbOff(3,1024) )

            !Faces: the part of each of the six faces not covered by live face neighbours
            if ( .not. live(i) ) cycle
            do ax = 1, 3
                select case ( ax )
                case ( 1 ); t1 = 2; t2 = 3
                case ( 2 ); t1 = 1; t2 = 3
                case default; t1 = 1; t2 = 2
                end select
                do s = -1, 1, 2
                    cover = 0._dp
                    do k = first, nnz
                        j = nbIdx(k)
                        if ( .not. live(j) ) cycle
                        d = nbOff(:,k)
                        if ( abs( s * d(ax) - 0.5_dp * ( abc(i,ax) + abc(j,ax) ) ) .gt. tol ) cycle
                        lo_ = max( -0.5_dp * abc(i,t1), d(t1) - 0.5_dp * abc(j,t1) )
                        hi_ = min(  0.5_dp * abc(i,t1), d(t1) + 0.5_dp * abc(j,t1) )
                        ov = max( 0._dp, hi_ - lo_ )
                        lo_ = max( -0.5_dp * abc(i,t2), d(t2) - 0.5_dp * abc(j,t2) )
                        hi_ = min(  0.5_dp * abc(i,t2), d(t2) + 0.5_dp * abc(j,t2) )
                        cover = cover + ov * max( 0._dp, hi_ - lo_ )
                    enddo
                    area = abc(i,t1) * abc(i,t2)
                    if ( cover .lt. area * ( 1._dp - 1.e-9_dp ) ) then
                        fc = pts(i,:)
                        fc(ax) = fc(ax) + s * 0.5_dp * abc(i,ax)
                        nrm = 0._dp
                        nrm(ax) = real( s, dp )
                        call pushInt( bfOwner, nbf, i )
                        call pushVec( bfC, nbf, fc )
                        call pushVec( bfN, nbf, nrm )
                    endif
                enddo
            enddo
        enddo

    contains

        !> Cells of one bucket (wrapped along periodic directions; all buckets of a direction
        !> with fewer than three, visited once) that touch cell i, through any periodic image
        subroutine visitBucket( bjIn, offIn )
        integer, intent(in) :: bjIn(3), offIn(3)
        integer :: lo3(3), hi3(3), x, y, z, bb, m, sx, sy, sz, nsh(3)
        real(dp) :: shift(3), dd(3)
        logical :: hit
            lo3 = bjIn
            hi3 = bjIn
            do m = 1, 3
                if ( nbk(m) .lt. 3 ) then
                    if ( offIn(m) .ne. 0 ) return
                    lo3(m) = 0
                    hi3(m) = nbk(m) - 1
                elseif ( bjIn(m) .lt. 0 .or. bjIn(m) .ge. nbk(m) ) then
                    if ( per(m) .le. 0._dp ) return
                    lo3(m) = modulo( bjIn(m), nbk(m) )
                    hi3(m) = lo3(m)
                endif
            enddo
            nsh = 0
            do m = 1, 3
                if ( per(m) .gt. 0._dp ) nsh(m) = 1
            enddo
            do z = lo3(3), hi3(3)
              do y = lo3(2), hi3(2)
                do x = lo3(1), hi3(1)
                    bb = 1 + x + nbk(1) * ( y + nbk(2) * z )
                    do m = bPtr(bb), bPtr(bb+1) - 1
                        j = bItem(m)
                        if ( j .eq. i .or. stamp(j) .eq. i ) cycle
                        hit = .false.
                        do sz = -nsh(3), nsh(3)
                          do sy = -nsh(2), nsh(2)
                            do sx = -nsh(1), nsh(1)
                                if ( hit ) cycle
                                shift = [ sx * per(1), sy * per(2), sz * per(3) ]
                                dd = pts(j,:) + shift - pts(i,:)
                                if ( all( abs(dd) - 0.5_dp * ( abc(j,:) + abc(i,:) ) .le. tol ) ) then
                                    hit = .true.
                                    stamp(j) = i
                                    call pushInt( nbIdx, nnz, j )
                                    call pushVec( nbOff, nnz, dd )
                                endif
                            enddo
                          enddo
                        enddo
                    enddo
                enddo
              enddo
            enddo
        end subroutine visitBucket

    end subroutine prismNeighboursAndFaces


    !>-----------------------------------------
    !> Stencils of live cells: live neighbours, plus their live neighbours when too few
    !---------------------------------------------------------------------------
    subroutine buildStencils( N, live, nbPtr, nbIdx, nbOff, stPtr, stIdx, stOff, nGrown, stMin, stMax )
    integer, intent(in) :: N, nbPtr(:), nbIdx(:)
    real(dp), intent(in) :: nbOff(:,:)
    logical, intent(in) :: live(:)
    integer, allocatable, intent(out) :: stPtr(:), stIdx(:)
    real(dp), allocatable, intent(out) :: stOff(:,:)
    integer, intent(out) :: nGrown, stMin, stMax
    integer, allocatable :: stamp(:), tmp(:)
    real(dp), allocatable :: toff(:,:)
    integer :: i, k, j, k2, j2, n1, n0, nnz, m
        allocate( stamp(N), tmp(N), toff(3,N), stPtr(N+1) )
        stamp = 0
        nGrown = 0
        stMin = huge(1)
        stMax = 0
        nnz = 0
        stPtr(1) = 1
        do i = 1, N
            n1 = 0
            if ( live(i) ) then
                stamp(i) = i
                do k = nbPtr(i), nbPtr(i+1) - 1
                    j = nbIdx(k)
                    if ( live(j) .and. stamp(j) .ne. i ) then
                        stamp(j) = i
                        n1 = n1 + 1
                        tmp(n1) = j
                        toff(:,n1) = nbOff(:,k)
                    endif
                enddo
                if ( n1 .lt. minStencil ) then
                    nGrown = nGrown + 1
                    n0 = n1
                    do m = 1, n0
                        j = tmp(m)
                        do k2 = nbPtr(j), nbPtr(j+1) - 1
                            j2 = nbIdx(k2)
                            if ( live(j2) .and. stamp(j2) .ne. i ) then
                                stamp(j2) = i
                                n1 = n1 + 1
                                tmp(n1) = j2
                                toff(:,n1) = toff(:,m) + nbOff(:,k2)
                            endif
                        enddo
                    enddo
                endif
                stMin = min( stMin, n1 )
                stMax = max( stMax, n1 )
            endif
            do m = 1, n1
                call pushInt( stIdx, nnz, tmp(m) )
                call pushVec( stOff, nnz, toff(:,m) )
            enddo
            stPtr(i+1) = nnz + 1
        enddo
        if ( .not. allocated(stIdx) ) then
            allocate( stIdx(1), stOff(3,1) )
        endif
    end subroutine buildStencils


    subroutine groupByOwner( N, owner, ptr, idx )
    integer, intent(in) :: N, owner(:)
    integer, allocatable, intent(out) :: ptr(:), idx(:)
    integer, allocatable :: cnt(:)
    integer :: f
        allocate( cnt(N), ptr(N+1), idx(max(1,size(owner))) )
        cnt = 0
        do f = 1, size(owner)
            cnt(owner(f)) = cnt(owner(f)) + 1
        enddo
        ptr(1) = 1
        do f = 1, N
            ptr(f+1) = ptr(f) + cnt(f)
        enddo
        cnt = 0
        do f = 1, size(owner)
            idx(ptr(owner(f)) + cnt(owner(f))) = f
            cnt(owner(f)) = cnt(owner(f)) + 1
        enddo
    end subroutine groupByOwner


    !>-----------------------------------------
    !> The constrained least squares fit of every live cell and the matrix entries
    !---------------------------------------------------------------------------
    subroutine assemble( N, cen, live, A0n, phase, Aint, stPtr, stIdx, stOff, bfPtr, bfIdx, bfC, bfN, &
                         cooR, cooC, cooV, nRankDef )
    integer, intent(in) :: N, phase(:), stPtr(:), stIdx(:), bfPtr(:), bfIdx(:)
    real(dp), intent(in) :: cen(:,:), A0n(:), Aint(:,:), bfC(:,:), bfN(:,:), stOff(:,:)
    logical, intent(in) :: live(:)
    integer, allocatable, intent(out) :: cooR(:), cooC(:)
    real(dp), allocatable, intent(out) :: cooV(:)
    integer, intent(out) :: nRankDef

    integer, allocatable :: rowStart(:)
    integer :: i, k, j, nd, nc, nr, f, o, kk, info, rank, lwork, ldb, e
    real(dp), allocatable :: A(:,:), B(:,:), S(:), work(:), w(:), lap(:)
    integer, allocatable :: iwork(:)
    real(dp) :: d(3), df(3), nrm(3), dOwner(3), scl, ell, wq(1), Aij, diag
    integer :: iwq(1)

        allocate( rowStart(N+1) )
        rowStart(1) = 1
        do i = 1, N
            rowStart(i+1) = rowStart(i) + merge( stPtr(i+1) - stPtr(i) + 1, 0, live(i) )
        enddo
        allocate( cooR(rowStart(N+1)-1), cooC(rowStart(N+1)-1), cooV(rowStart(N+1)-1) )
        nRankDef = 0

        !$OMP PARALLEL DO SCHEDULE(dynamic,64) DEFAULT(shared) &
        !$OMP PRIVATE(i,k,j,nd,nc,nr,f,o,kk,info,rank,lwork,ldb,e,A,B,S,work,w,lap,iwork,d,df,nrm,dOwner,scl,ell,wq,Aij,diag,iwq) &
        !$OMP REDUCTION(+:nRankDef)
        do i = 1, N
            if ( .not. live(i) ) cycle
            nd = stPtr(i+1) - stPtr(i)
            !boundary faces of i and of its stencil cells
            nc = bfPtr(i+1) - bfPtr(i)
            do k = stPtr(i), stPtr(i+1) - 1
                nc = nc + bfPtr(stIdx(k)+1) - bfPtr(stIdx(k))
            enddo
            nr = nd + nc
            ldb = max( nr, nUnk )
            allocate( A(nr,nUnk), B(ldb,max(nd,1)), S(nUnk), w(nr), lap(max(nd,1)) )

            !length scale of the stencil, for conditioning
            scl = 0._dp
            do k = 1, nd
                scl = scl + sum( abs( stOff(:,stPtr(i)+k-1) ) )
            enddo
            scl = max( scl / max( 3 * nd, 1 ), tiny(1._dp) )

            do k = 1, nd
                d = stOff(:,stPtr(i)+k-1) / scl
                A(k,:) = [ d(1), d(2), d(3), 0.5_dp*d(1)**2, 0.5_dp*d(2)**2, 0.5_dp*d(3)**2, &
                           d(1)*d(2), d(1)*d(3), d(2)*d(3) ]
                w(k) = norm2(d) ** (-weightExp)
            enddo
            !Neumann rows n.(g + H d_f) = 0, scaled by |d_f| so that they have the dimension of
            !the data rows. The face offset goes through its owner's offset.
            kk = nd
            do e = 0, nd
                if ( e .eq. 0 ) then
                    o = i
                    dOwner = 0._dp
                else
                    o = stIdx(stPtr(i)+e-1)
                    dOwner = stOff(:,stPtr(i)+e-1)
                endif
                do f = bfPtr(o), bfPtr(o+1) - 1
                    kk = kk + 1
                    df = ( dOwner + bfC(:,bfIdx(f)) - cen(:,o) ) / scl
                    nrm = bfN(:,bfIdx(f))
                    ell = max( norm2(df), 1.e-3_dp )
                    A(kk,:) = ell * [ nrm(1), nrm(2), nrm(3), nrm(1)*df(1), nrm(2)*df(2), nrm(3)*df(3), &
                                      nrm(1)*df(2) + nrm(2)*df(1), nrm(1)*df(3) + nrm(3)*df(1), &
                                      nrm(2)*df(3) + nrm(3)*df(2) ]
                    w(kk) = ell ** (-weightExp)
                enddo
            enddo

            w = sqrt( w )
            do k = 1, nr
                A(k,:) = w(k) * A(k,:)
            enddo
            B = 0._dp
            do k = 1, nd
                B(k,k) = w(k)
            enddo

            call dgelsd( nr, nUnk, nd, A, nr, B, ldb, S, rcondFit, rank, wq, -1, iwq, info )
            lwork = int( wq(1) )
            allocate( work(lwork), iwork(max(1,iwq(1))) )
            call dgelsd( nr, nUnk, nd, A, nr, B, ldb, S, rcondFit, rank, work, lwork, iwork, info )
            if ( info .ne. 0 .or. rank .lt. nUnk ) nRankDef = nRankDef + 1

            !trace of the fitted Hessian: rows 4-6 of the solution, back in physical units
            lap(1:nd) = ( B(4,1:nd) + B(5,1:nd) + B(6,1:nd) ) / scl**2

            diag = 0._dp
            do k = 1, nd
                j = stIdx(stPtr(i)+k-1)
                Aij = pairExchange( i, j, A0n, phase, Aint )
                cooR(rowStart(i)+k-1) = i
                cooC(rowStart(i)+k-1) = j
                cooV(rowStart(i)+k-1) = lap(k) * Aij
                diag = diag - lap(k) * Aij
            enddo
            cooR(rowStart(i+1)-1) = i
            cooC(rowStart(i+1)-1) = i
            cooV(rowStart(i+1)-1) = diag

            deallocate( A, B, S, w, lap, work, iwork )
        enddo
        !$OMP END PARALLEL DO

    end subroutine assemble


    !=================================================================== symmetric operator (SUSHI-g)

    !>-----------------------------------------
    !> @brief Builds the exchange matrix of an unstructured grid with the symmetric, energy based
    !> operator SUSHI-g (see the module header).
    !> @param[in] grid the grid (tetrahedra: nodes, elements; prisms: pts, abc)
    !> @param[in] pbc periodic exchange along x, y, z
    !> @param[in] A0n exchange stiffness per cell, normalised as for the Direct Laplacian
    !> @param[in] phase material index per cell
    !> @param[in] Aint interface exchange per phase pair (not used by this operator, see header)
    !> @param[inout] GridInfo receives the COO copy of the matrix
    !> @param[out] Aexch_matrix the CSR exchange matrix
    !---------------------------------------------------------------------------
    subroutine computeExchangeSushi( grid, pbc, A0n, phase, Aint, GridInfo, Aexch_matrix )
    type(MicroMagGrid), intent(in) :: grid
    logical, intent(in) :: pbc(3)
    real(dp), intent(in) :: A0n(:)
    integer, intent(in) :: phase(:)
    real(dp), intent(in) :: Aint(:,:)
    type(MicroMagGridInfo), intent(inout) :: GridInfo
    type(sparse_matrix_t), intent(out) :: Aexch_matrix

    real(dp), allocatable :: cen(:,:), vol(:), nbOff(:,:), bfC(:,:), bfN(:,:), fnOff(:,:)
    real(dp), allocatable :: fA(:), fX(:,:), fN(:,:), fD(:,:), vals(:), cooV(:)
    integer, allocatable :: nbPtr(:), nbIdx(:), bfOwner(:), fK(:), fL(:), fnPtr(:), fnIdx(:)
    integer, allocatable :: cfPtr(:), cfIdx(:), uPtr(:), uIdx(:), rowPtr(:), colIdx(:)
    integer, allocatable :: cooR(:), cooC(:), fKind(:), fOk(:), nU(:), tmpU(:)
    type(faceStencil), allocatable :: fs(:)
    logical, allocatable :: live(:)
    real(dp) :: per(3), lo(3), hi(3), stab, p(3,4)
    integer :: N, i, f, nF, nbf, nBad, nNeg, k, nnz, nInt
    character(len=256) :: msg

        if ( grid%gridType .eq. gridTypeTetrahedron ) then
            N = size(grid%elements, 2)
            stab = stabTet
        else
            N = size(grid%pts, 1)
            stab = stabBox
        endif
        allocate( live(N) )
        live = A0n(1:N) .gt. 0._dp

        !--- centroids, volumes, extent and periods
        allocate( cen(3,N), vol(N) )
        if ( grid%gridType .eq. gridTypeTetrahedron ) then
            do i = 1, N
                p = grid%nodes(:, grid%elements(1:4,i))
                cen(:,i) = sum( p, dim=2 ) / 4._dp
                vol(i) = abs( dot_product( p(:,2) - p(:,1), cross( p(:,3) - p(:,1), p(:,4) - p(:,1) ) ) ) / 6._dp
            enddo
            lo = minval( grid%nodes, dim=2 )
            hi = maxval( grid%nodes, dim=2 )
        else
            cen = transpose( grid%pts(1:N,:) )
            vol = product( grid%abc(1:N,:), dim=2 )
            lo = minval( grid%pts(1:N,:) - 0.5_dp * grid%abc(1:N,:), dim=1 )
            hi = maxval( grid%pts(1:N,:) + 0.5_dp * grid%abc(1:N,:), dim=1 )
        endif
        per = 0._dp
        where ( pbc ) per = hi - lo

        !--- vertex neighbours (node sharing / touching, with their actual offsets) and the faces
        nBad = 0
        if ( grid%gridType .eq. gridTypeTetrahedron ) then
            call tetNeighboursAndFaces( grid%nodes, grid%elements(1:4,1:N), cen, live, pbc, per, lo, &
                                        nbPtr, nbIdx, nbOff, bfOwner, bfC, bfN, nbf, nBad )
            call tetFaceList( grid%nodes, grid%elements(1:4,1:N), cen, live, pbc, per, lo, &
                              nF, fK, fL, fA, fX, fN, fD )
        else
            call prismNeighboursAndFaces( grid%pts(1:N,:), grid%abc(1:N,:), live, per, lo, hi, &
                                          nbPtr, nbIdx, nbOff, bfOwner, bfC, bfN, nbf )
            call prismFaceList( grid%pts(1:N,:), grid%abc(1:N,:), live, nbPtr, nbIdx, nbOff, &
                                nF, fK, fL, fA, fX, fN, fD )
        endif
        if ( nBad .gt. 0 ) then
            write(msg,'(A,I0,A)') 'Symmetric exchange: ', nBad, &
                ' faces are shared by more than two tetrahedra - the mesh is not conforming'
            call displayGUIMessage( trim(msg) )
        endif

        !--- faces of every cell (signed: + as K, - as L) and the face neighbours
        call cellFaces( N, nF, fK, fL, fD, cfPtr, cfIdx, fnPtr, fnIdx, fnOff )

        !--- interpolation of the interior face values
        allocate( fs(nF), fKind(nF), fOk(nF) )
        fKind = 0
        fOk = 1
        !$OMP PARALLEL DO SCHEDULE(dynamic,256) DEFAULT(shared) PRIVATE(f)
        do f = 1, nF
            if ( fL(f) .eq. 0 ) cycle
            call faceInterpolation( f, cen, vol, live, fK, fL, fX, fD, fnPtr, fnIdx, fnOff, nbPtr, nbIdx, nbOff, &
                                    fs(f), fKind(f), fOk(f) )
        enddo
        !$OMP END PARALLEL DO

        !--- the cells each cell's energy couples (its own and those of its face interpolations)
        allocate( nU(N), uPtr(N+1) )
        nU = 0
        !$OMP PARALLEL DO SCHEDULE(dynamic,256) DEFAULT(shared) PRIVATE(i,tmpU)
        do i = 1, N
            if ( .not. live(i) ) cycle
            call cellUnion( i, cfPtr, cfIdx, fL, fs, tmpU, nU(i) )
        enddo
        !$OMP END PARALLEL DO
        uPtr(1) = 1
        do i = 1, N
            uPtr(i+1) = uPtr(i) + nU(i)
        enddo
        allocate( uIdx(max(1,uPtr(N+1)-1)) )
        !$OMP PARALLEL DO SCHEDULE(dynamic,256) DEFAULT(shared) PRIVATE(i,tmpU,k)
        do i = 1, N
            if ( .not. live(i) ) cycle
            call cellUnion( i, cfPtr, cfIdx, fL, fs, tmpU, k )
            uIdx(uPtr(i):uPtr(i+1)-1) = tmpU(1:k)
        enddo
        !$OMP END PARALLEL DO

        !--- sparsity pattern, then the element matrices summed into it
        call sushiPattern( N, uPtr, uIdx, rowPtr, colIdx )
        allocate( vals(rowPtr(N+1)-1) )
        vals = 0._dp
        nNeg = 0
        !$OMP PARALLEL DO SCHEDULE(dynamic,64) DEFAULT(shared) PRIVATE(i) REDUCTION(+:nNeg)
        do i = 1, N
            if ( .not. live(i) ) cycle
            call cellMatrix( i, cen, vol, A0n, stab, cfPtr, cfIdx, fK, fL, fA, fX, fN, fD, fs, &
                             uIdx(uPtr(i):uPtr(i+1)-1), rowPtr, colIdx, vals, nNeg )
        enddo
        !$OMP END PARALLEL DO

        !--- exact symmetry, constants in the kernel, and L = -V^-1 K
        call symmetriseAndScale( N, vol, rowPtr, colIdx, vals )

        nnz = rowPtr(N+1) - 1
        allocate( cooR(nnz), cooC(nnz), cooV(nnz) )
        do i = 1, N
            cooR(rowPtr(i):rowPtr(i+1)-1) = i
        enddo
        cooC = colIdx(1:nnz)
        cooV = vals(1:nnz)
        call create_CSR_matrix( cooR, cooC, cooV, N, N, 0._dp, Aexch_matrix )
        call create_COO_values_from_CSR( Aexch_matrix, GridInfo )

        nInt = count( fL(1:nF) .gt. 0 )
        write(msg,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,F6.1,A)') 'Symmetric exchange (SUSHI-g): ', N, &
            ' cells, ', nInt, ' interior faces (two-point ', count( fKind .eq. 1 ), ', face neighbours ', &
            count( fKind .eq. 2 ), ', vertex neighbours ', count( fKind .eq. 3 ), '), ', nF - nInt, &
            ' boundary faces, ', real(nnz,dp) / max(1, count(live)), ' entries per row'
        call displayGUIMessage( trim(msg) )
        if ( count( fOk .eq. 0 ) .gt. 0 .or. nNeg .gt. 0 ) then
            write(msg,'(A,I0,A,I0,A)') 'Symmetric exchange: ', count( fOk .eq. 0 ), &
                ' face interpolations not affine exact, ', nNeg, &
                ' cell faces with the centroid on or behind the face (inverted cells?)'
            call displayGUIMessage( trim(msg) )
        endif
        if ( size(Aint,1) .gt. 1 ) then
            if ( any( Aint .ge. 0._dp ) ) call displayGUIMessage( 'Symmetric exchange: interface exchange ' // &
                'values are not used, every cell carries its own A' )
        endif
        if ( any(pbc) ) call displayGUIMessage( 'Symmetric exchange: periodic along the exchPBC directions only' )

    end subroutine computeExchangeSushi


    !>-----------------------------------------
    !> Tetrahedra: all faces. Interior faces (two live cells) carry the offset from the centroid of
    !> K to the image of the centroid of L adjacent through the face; a face of one live cell (outer
    !> surface or towards a cell with A0 = 0) is a free surface.
    !---------------------------------------------------------------------------
    subroutine tetFaceList( nodes, ele, cen, live, pbc, per, lo, nF, fK, fL, fA, fX, fN, fD )
    real(dp), intent(in) :: nodes(:,:), cen(:,:), per(3), lo(3)
    integer, intent(in) :: ele(:,:)
    logical, intent(in) :: live(:), pbc(3)
    integer, intent(out) :: nF
    integer, allocatable, intent(out) :: fK(:), fL(:)
    real(dp), allocatable, intent(out) :: fA(:), fX(:,:), fN(:,:), fD(:,:)

    integer, allocatable :: Ec(:,:), cnt(:), own(:), loc(:), oth(:)
    integer :: N, nCanon, i, k, a, b, c, j, slot, nt, K1, L1, kf, v, tri(3)
    real(dp) :: p(3,3), nrm(3)
    type(keyTable) :: ftab

        N = size(ele, 2)
        call tetCanonicalNodes( nodes, ele, pbc, per, lo, Ec, nCanon )
        allocate( cnt(4*N), own(4*N), loc(4*N), oth(4*N) )
        call tableInit( ftab, 4*N )
        nt = 0
        do i = 1, N
            do k = 1, 4
                tri = pack( Ec(:,i), [ (a .ne. k, a = 1, 4) ] )
                call sort3( tri )
                slot = tableFind( ftab, int(tri,i8) )
                if ( slot .gt. 0 ) then
                    j = ftab%val(slot)
                    cnt(j) = cnt(j) + 1
                    oth(j) = i
                else
                    nt = nt + 1
                    call tableInsert( ftab, int(tri,i8), nt )
                    cnt(nt) = 1
                    own(nt) = i
                    loc(nt) = k
                    oth(nt) = 0
                endif
            enddo
        enddo

        allocate( fK(nt), fL(nt), fA(nt), fX(3,nt), fN(3,nt), fD(3,nt) )
        nF = 0
        do j = 1, nt
            K1 = 0
            L1 = 0
            if ( cnt(j) .eq. 1 ) then
                if ( live(own(j)) ) K1 = own(j)
            elseif ( cnt(j) .eq. 2 ) then
                if ( live(own(j)) ) then
                    K1 = own(j)
                    if ( live(oth(j)) ) L1 = oth(j)
                elseif ( live(oth(j)) ) then
                    K1 = oth(j)
                endif
            endif
            if ( K1 .eq. 0 ) cycle
            !local face of K1 with this canonical triple
            if ( K1 .eq. own(j) ) then
                kf = loc(j)
            else
                do kf = 1, 4
                    tri = pack( Ec(:,K1), [ (a .ne. kf, a = 1, 4) ] )
                    call sort3( tri )
                    slot = tableFind( ftab, int(tri,i8) )
                    if ( ftab%val(slot) .eq. j ) exit
                enddo
            endif
            c = 0
            do a = 1, 4
                if ( a .eq. kf ) cycle
                c = c + 1
                p(:,c) = nodes(:,ele(a,K1))
            enddo
            nF = nF + 1
            fK(nF) = K1
            fL(nF) = L1
            fX(:,nF) = sum( p, dim=2 ) / 3._dp
            nrm = cross( p(:,2) - p(:,1), p(:,3) - p(:,1) )
            fA(nF) = 0.5_dp * norm2( nrm )
            nrm = nrm / norm2( nrm )
            if ( dot_product( nrm, nodes(:,ele(kf,K1)) - fX(:,nF) ) .gt. 0._dp ) nrm = -nrm
            fN(:,nF) = nrm
            fD(:,nF) = 0._dp
            if ( L1 .gt. 0 ) then
                !offset through a shared node, with its coordinates as each of the two cells has it
                a = merge( 2, 1, kf .eq. 1 )
                v = Ec(a,K1)
                do b = 1, 4
                    if ( Ec(b,L1) .eq. v ) exit
                enddo
                fD(:,nF) = ( cen(:,L1) - nodes(:,ele(b,L1)) ) - ( cen(:,K1) - nodes(:,ele(a,K1)) )
            endif
        enddo
    end subroutine tetFaceList


    !>-----------------------------------------
    !> Prisms: the overlap rectangle of every pair of live face neighbours is an interior face
    !> (also across hanging faces); the uncovered remainder of a cell face is a free surface with
    !> the centroid of the remainder. Edge and corner contacts are not faces.
    !---------------------------------------------------------------------------
    subroutine prismFaceList( pts, abc, live, nbPtr, nbIdx, nbOff, nF, fK, fL, fA, fX, fN, fD )
    real(dp), intent(in) :: pts(:,:), abc(:,:), nbOff(:,:)
    integer, intent(in) :: nbPtr(:), nbIdx(:)
    logical, intent(in) :: live(:)
    integer, intent(out) :: nF
    integer, allocatable, intent(out) :: fK(:), fL(:)
    real(dp), allocatable, intent(out) :: fA(:), fX(:,:), fN(:,:), fD(:,:)

    integer :: N, i, j, k, ax, s, t1, t2, cap
    real(dp) :: tol, d(3), lo1, hi1, lo2, hi2, ov, cover, mom(3), area, xs(3), fc(3), nrm(3)

        N = size(pts, 1)
        tol = 1.e-6_dp * minval( abc )
        cap = 6 * N + 16
        allocate( fK(cap), fL(cap), fA(cap), fX(3,cap), fN(3,cap), fD(3,cap) )
        nF = 0
        do i = 1, N
            if ( .not. live(i) ) cycle
            do ax = 1, 3
                select case ( ax )
                case ( 1 ); t1 = 2; t2 = 3
                case ( 2 ); t1 = 1; t2 = 3
                case default; t1 = 1; t2 = 2
                end select
                do s = -1, 1, 2
                    nrm = 0._dp
                    nrm(ax) = real( s, dp )
                    cover = 0._dp
                    mom = 0._dp
                    do k = nbPtr(i), nbPtr(i+1) - 1
                        j = nbIdx(k)
                        if ( .not. live(j) ) cycle
                        d = nbOff(:,k)
                        if ( abs( s * d(ax) - 0.5_dp * ( abc(i,ax) + abc(j,ax) ) ) .gt. tol ) cycle
                        lo1 = max( -0.5_dp * abc(i,t1), d(t1) - 0.5_dp * abc(j,t1) )
                        hi1 = min(  0.5_dp * abc(i,t1), d(t1) + 0.5_dp * abc(j,t1) )
                        lo2 = max( -0.5_dp * abc(i,t2), d(t2) - 0.5_dp * abc(j,t2) )
                        hi2 = min(  0.5_dp * abc(i,t2), d(t2) + 0.5_dp * abc(j,t2) )
                        if ( hi1 - lo1 .le. tol .or. hi2 - lo2 .le. tol ) cycle
                        ov = ( hi1 - lo1 ) * ( hi2 - lo2 )
                        xs = pts(i,:)
                        xs(ax) = xs(ax) + s * 0.5_dp * abc(i,ax)
                        xs(t1) = xs(t1) + 0.5_dp * ( lo1 + hi1 )
                        xs(t2) = xs(t2) + 0.5_dp * ( lo2 + hi2 )
                        cover = cover + ov
                        mom = mom + ov * xs
                        if ( i .lt. j ) then
                            call addFace( i, j, ov, xs, nrm, d )
                        endif
                    enddo
                    area = abc(i,t1) * abc(i,t2)
                    if ( cover .lt. area * ( 1._dp - 1.e-9_dp ) ) then
                        fc = pts(i,:)
                        fc(ax) = fc(ax) + s * 0.5_dp * abc(i,ax)
                        xs = ( area * fc - mom ) / ( area - cover )
                        d = 0._dp
                        call addFace( i, 0, area - cover, xs, nrm, d )
                    endif
                enddo
            enddo
        enddo

    contains

        subroutine addFace( K1, L1, ar, x, n, dd )
        integer, intent(in) :: K1, L1
        real(dp), intent(in) :: ar, x(3), n(3), dd(3)
        integer, allocatable :: ti(:)
        real(dp), allocatable :: tr(:), tv(:,:)
            if ( nF .ge. size(fK) ) then
                allocate( ti(2*size(fK)) ); ti(1:nF) = fK(1:nF); call move_alloc( ti, fK )
                allocate( ti(2*size(fL)) ); ti(1:nF) = fL(1:nF); call move_alloc( ti, fL )
                allocate( tr(2*size(fA)) ); tr(1:nF) = fA(1:nF); call move_alloc( tr, fA )
                allocate( tv(3,2*size(fX,2)) ); tv(:,1:nF) = fX(:,1:nF); call move_alloc( tv, fX )
                allocate( tv(3,2*size(fN,2)) ); tv(:,1:nF) = fN(:,1:nF); call move_alloc( tv, fN )
                allocate( tv(3,2*size(fD,2)) ); tv(:,1:nF) = fD(:,1:nF); call move_alloc( tv, fD )
            endif
            nF = nF + 1
            fK(nF) = K1
            fL(nF) = L1
            fA(nF) = ar
            fX(:,nF) = x
            fN(:,nF) = n
            fD(:,nF) = dd
        end subroutine addFace

    end subroutine prismFaceList


    !>-----------------------------------------
    !> Faces of every cell (signed face ids: + the cell is K, - it is L) and the face neighbours of
    !> every cell with their offsets
    !---------------------------------------------------------------------------
    subroutine cellFaces( N, nF, fK, fL, fD, cfPtr, cfIdx, fnPtr, fnIdx, fnOff )
    integer, intent(in) :: N, nF, fK(:), fL(:)
    real(dp), intent(in) :: fD(:,:)
    integer, allocatable, intent(out) :: cfPtr(:), cfIdx(:), fnPtr(:), fnIdx(:)
    real(dp), allocatable, intent(out) :: fnOff(:,:)
    integer, allocatable :: cnt(:), ncnt(:)
    integer :: f
        allocate( cnt(N), ncnt(N), cfPtr(N+1), fnPtr(N+1) )
        cnt = 0
        ncnt = 0
        do f = 1, nF
            cnt(fK(f)) = cnt(fK(f)) + 1
            if ( fL(f) .gt. 0 ) then
                cnt(fL(f)) = cnt(fL(f)) + 1
                ncnt(fK(f)) = ncnt(fK(f)) + 1
                ncnt(fL(f)) = ncnt(fL(f)) + 1
            endif
        enddo
        cfPtr(1) = 1
        fnPtr(1) = 1
        do f = 1, N
            cfPtr(f+1) = cfPtr(f) + cnt(f)
            fnPtr(f+1) = fnPtr(f) + ncnt(f)
        enddo
        allocate( cfIdx(max(1,cfPtr(N+1)-1)), fnIdx(max(1,fnPtr(N+1)-1)), fnOff(3,max(1,fnPtr(N+1)-1)) )
        cnt = 0
        ncnt = 0
        do f = 1, nF
            cfIdx(cfPtr(fK(f)) + cnt(fK(f))) = f
            cnt(fK(f)) = cnt(fK(f)) + 1
            if ( fL(f) .gt. 0 ) then
                cfIdx(cfPtr(fL(f)) + cnt(fL(f))) = -f
                cnt(fL(f)) = cnt(fL(f)) + 1
                fnIdx(fnPtr(fK(f)) + ncnt(fK(f))) = fL(f)
                fnOff(:,fnPtr(fK(f)) + ncnt(fK(f))) = fD(:,f)
                ncnt(fK(f)) = ncnt(fK(f)) + 1
                fnIdx(fnPtr(fL(f)) + ncnt(fL(f))) = fK(f)
                fnOff(:,fnPtr(fL(f)) + ncnt(fL(f))) = -fD(:,f)
                ncnt(fL(f)) = ncnt(fL(f)) + 1
            endif
        enddo
    end subroutine cellFaces


    !>-----------------------------------------
    !> Interpolation weights of the value at the centroid of interior face f, exact for affine
    !> fields: two-point when the face point lies on the segment between the two centroids; else
    !> the two cells and their face neighbours (|d|^-4 weighted least squares, minimum norm), all of
    !> them when there are at most maxFull cells, otherwise the nearest-first prefixes until the
    !> weights are well conditioned (sum |beta| <= bMax); else the vertex neighbours of both cells.
    !> kind = 1, 2, 3 for the three cases; ok = 0 if not even the last one is affine exact.
    !---------------------------------------------------------------------------
    subroutine faceInterpolation( f, cen, vol, live, fK, fL, fX, fD, fnPtr, fnIdx, fnOff, nbPtr, nbIdx, nbOff, &
                                  st, kind, ok )
    integer, intent(in) :: f, fK(:), fL(:), fnPtr(:), fnIdx(:), nbPtr(:), nbIdx(:)
    real(dp), intent(in) :: cen(:,:), vol(:), fX(:,:), fD(:,:), fnOff(:,:), nbOff(:,:)
    logical, intent(in) :: live(:)
    type(faceStencil), intent(out) :: st
    integer, intent(out) :: kind, ok

    integer :: K1, L1, mx, nc, n, e, c, kk, j, a, nlo, tj
    integer, allocatable :: id(:), cid(:)
    real(dp), allocatable :: Y(:,:), cy(:,:), beta(:), dist(:)
    real(dp) :: xK(3), xL(3), base(3), ty(3), td
    logical :: good

        K1 = fK(f)
        L1 = fL(f)
        xK = cen(:,K1) - fX(:,f)
        xL = xK + fD(:,f)
        mx = 2 + ( fnPtr(K1+1) - fnPtr(K1) ) + ( fnPtr(L1+1) - fnPtr(L1) ) &
               + ( nbPtr(K1+1) - nbPtr(K1) ) + ( nbPtr(L1+1) - nbPtr(L1) )
        allocate( id(mx), cid(mx), Y(3,mx), cy(3,mx), beta(mx), dist(mx) )
        id(1) = K1
        id(2) = L1
        Y(:,1) = xK
        Y(:,2) = xL

        !two-point
        call wlsWeights( Y, vol(id(1:2)), 2, beta, good )
        if ( good .and. sum( abs( beta(1:2) ) ) .le. 1._dp + 1.e-9_dp ) then
            call store( 2, 1, 1 )
            return
        endif

        !face neighbours of both cells, nearest first
        nc = 0
        do e = 1, 2
            if ( e .eq. 1 ) then
                c = K1
                base = xK
            else
                c = L1
                base = xL
            endif
            do kk = fnPtr(c), fnPtr(c+1) - 1
                j = fnIdx(kk)
                if ( j .eq. K1 .or. j .eq. L1 ) cycle
                if ( any( cid(1:nc) .eq. j ) ) cycle
                nc = nc + 1
                cid(nc) = j
                cy(:,nc) = base + fnOff(:,kk)
                dist(nc) = norm2( cy(:,nc) )
            enddo
        enddo
        do a = 2, nc
            tj = cid(a)
            ty = cy(:,a)
            td = dist(a)
            kk = a - 1
            do while ( kk .ge. 1 )
                !nearest first, equal distances (relative 1e-9) by cell index
                if ( abs( dist(kk) - td ) .le. 1.e-9_dp * max( dist(kk), td ) ) then
                    if ( cid(kk) .lt. tj ) exit
                elseif ( dist(kk) .lt. td ) then
                    exit
                endif
                cid(kk+1) = cid(kk)
                cy(:,kk+1) = cy(:,kk)
                dist(kk+1) = dist(kk)
                kk = kk - 1
            enddo
            cid(kk+1) = tj
            cy(:,kk+1) = ty
            dist(kk+1) = td
        enddo
        if ( nc + 2 .le. maxFull ) then
            nlo = nc
        else
            nlo = min( 2, nc )
        endif
        do n = nlo, nc
            if ( n .lt. 1 ) cycle
            id(3:2+n) = cid(1:n)
            Y(:,3:2+n) = cy(:,1:n)
            call wlsWeights( Y, vol(id(1:2+n)), 2+n, beta, good )
            if ( good .and. sum( abs( beta(1:2+n) ) ) .le. bMax ) then
                call store( 2+n, 2, 1 )
                return
            endif
        enddo

        !vertex neighbours of both cells
        n = 2
        do e = 1, 2
            if ( e .eq. 1 ) then
                c = K1
                base = xK
            else
                c = L1
                base = xL
            endif
            do kk = nbPtr(c), nbPtr(c+1) - 1
                j = nbIdx(kk)
                if ( .not. live(j) ) cycle
                if ( any( id(1:n) .eq. j ) ) cycle
                n = n + 1
                id(n) = j
                Y(:,n) = base + nbOff(:,kk)
            enddo
        enddo
        call wlsWeights( Y, vol(id(1:n)), n, beta, good )
        call store( n, 3, merge( 1, 0, good ) )

    contains

        subroutine store( m, k_, ok_ )
        integer, intent(in) :: m, k_, ok_
            allocate( st%idx(m), st%w(m) )
            st%idx = id(1:m)
            st%w = beta(1:m)
            kind = k_
            ok = ok_
        end subroutine store

    end subroutine faceInterpolation


    !>-----------------------------------------
    !> Minimum-norm weights beta (sum beta_j phi_j = phi(0) for affine phi) from the points Y
    !> (relative to the target point), weights (v_j / mean v)^(1/2) |y|^-4 after scaling by the rms
    !> distance. The volume factor keeps small cells (slivers) from entering the face values of
    !> their neighbours with O(1) weights, which would make L = -V^-1 K very stiff.
    !---------------------------------------------------------------------------
    subroutine wlsWeights( Y, v, n, beta, ok )
    real(dp), intent(in) :: Y(:,:), v(:)
    integer, intent(in) :: n
    real(dp), intent(out) :: beta(:)
    logical, intent(out) :: ok
    real(dp) :: sc, y3(3), g(4), M4(4,4), c(4), w, s(3), vm
    integer :: j, a
        sc = sqrt( sum( Y(:,1:n)**2 ) / n )
        vm = sum( v(1:n) ) / n
        M4 = 0._dp
        do j = 1, n
            y3 = Y(:,j) / sc
            w = sqrt( v(j) / vm ) / ( sum( y3**2 ) + 1.e-6_dp )**2
            g = [ 1._dp, y3 ]
            do a = 1, 4
                M4(:,a) = M4(:,a) + w * g * g(a)
            enddo
        enddo
        call symPinv( M4, 4, 1.e-10_dp )
        c = M4(:,1)
        s = 0._dp
        do j = 1, n
            y3 = Y(:,j) / sc
            w = sqrt( v(j) / vm ) / ( sum( y3**2 ) + 1.e-6_dp )**2
            beta(j) = w * ( c(1) + dot_product( y3, c(2:4) ) )
            s = s + beta(j) * y3
        enddo
        ok = abs( sum( beta(1:n) ) - 1._dp ) .lt. 1.e-9_dp .and. norm2( s ) .lt. 1.e-9_dp
    end subroutine wlsWeights


    !> Pseudo-inverse of a symmetric matrix in place (eigenvalues below rcond * max |eigenvalue| dropped)
    subroutine symPinv( A, n, rcond )
    integer, intent(in) :: n
    real(dp), intent(inout) :: A(n,n)
    real(dp), intent(in) :: rcond
    real(dp) :: lam(n), V(n,n), work(max(1,66*n)), lmax
    integer :: info, k
        V = A
        call dsyev( 'V', 'U', n, V, n, lam, work, size(work), info )
        lmax = maxval( abs(lam) )
        A = 0._dp
        if ( info .ne. 0 .or. lmax .le. 0._dp ) return
        do k = 1, n
            if ( abs(lam(k)) .gt. rcond * lmax ) then
                A = A + spread( V(:,k), 2, n ) * spread( V(:,k), 1, n ) / lam(k)
            endif
        enddo
    end subroutine symPinv


    !> The cells coupled by the energy of cell i: i and the interpolation stencils of its interior faces
    subroutine cellUnion( i, cfPtr, cfIdx, fL, fs, U, nU )
    integer, intent(in) :: i, cfPtr(:), cfIdx(:), fL(:)
    type(faceStencil), intent(in) :: fs(:)
    integer, allocatable, intent(inout) :: U(:)
    integer, intent(out) :: nU
    integer :: k, f, m, j, cap
        cap = 1
        do k = cfPtr(i), cfPtr(i+1) - 1
            f = abs( cfIdx(k) )
            if ( fL(f) .gt. 0 ) cap = cap + size( fs(f)%idx )
        enddo
        if ( allocated(U) ) deallocate( U )
        allocate( U(cap) )
        nU = 1
        U(1) = i
        do k = cfPtr(i), cfPtr(i+1) - 1
            f = abs( cfIdx(k) )
            if ( fL(f) .eq. 0 ) cycle
            do m = 1, size( fs(f)%idx )
                j = fs(f)%idx(m)
                if ( any( U(1:nU) .eq. j ) ) cycle
                nU = nU + 1
                U(nU) = j
            enddo
        enddo
    end subroutine cellUnion


    !>-----------------------------------------
    !> CSR pattern of K = sum_i P_i^T S_i P_i: row a holds the union of U_i over the cells i whose
    !> U_i contains a (symmetric by construction), columns sorted
    !---------------------------------------------------------------------------
    subroutine sushiPattern( N, uPtr, uIdx, rowPtr, colIdx )
    integer, intent(in) :: N, uPtr(:), uIdx(:)
    integer, allocatable, intent(out) :: rowPtr(:), colIdx(:)
    integer, allocatable :: rcPtr(:), rcIdx(:), cnt(:), stamp(:), tmp(:)
    integer :: i, k, a, m, c, nnz, nr, b, t
        allocate( cnt(N), rcPtr(N+1), stamp(N), tmp(N), rowPtr(N+1) )
        cnt = 0
        do i = 1, N
            do k = uPtr(i), uPtr(i+1) - 1
                cnt(uIdx(k)) = cnt(uIdx(k)) + 1
            enddo
        enddo
        rcPtr(1) = 1
        do a = 1, N
            rcPtr(a+1) = rcPtr(a) + cnt(a)
        enddo
        allocate( rcIdx(max(1,rcPtr(N+1)-1)) )
        cnt = 0
        do i = 1, N
            do k = uPtr(i), uPtr(i+1) - 1
                a = uIdx(k)
                rcIdx(rcPtr(a) + cnt(a)) = i
                cnt(a) = cnt(a) + 1
            enddo
        enddo
        !two passes over the rows: count, then fill and sort
        stamp = 0
        rowPtr(1) = 1
        do a = 1, N
            nr = 0
            do m = rcPtr(a), rcPtr(a+1) - 1
                i = rcIdx(m)
                do k = uPtr(i), uPtr(i+1) - 1
                    c = uIdx(k)
                    if ( stamp(c) .eq. a ) cycle
                    stamp(c) = a
                    nr = nr + 1
                enddo
            enddo
            rowPtr(a+1) = rowPtr(a) + nr
        enddo
        nnz = rowPtr(N+1) - 1
        allocate( colIdx(max(1,nnz)) )
        stamp = 0
        do a = 1, N
            nr = 0
            do m = rcPtr(a), rcPtr(a+1) - 1
                i = rcIdx(m)
                do k = uPtr(i), uPtr(i+1) - 1
                    c = uIdx(k)
                    if ( stamp(c) .eq. a ) cycle
                    stamp(c) = a
                    nr = nr + 1
                    tmp(nr) = c
                enddo
            enddo
            do b = 2, nr
                t = tmp(b)
                k = b - 1
                do while ( k .ge. 1 )
                    if ( tmp(k) .le. t ) exit
                    tmp(k+1) = tmp(k)
                    k = k - 1
                enddo
                tmp(k+1) = t
            enddo
            colIdx(rowPtr(a):rowPtr(a+1)-1) = tmp(1:nr)
        enddo
    end subroutine sushiPattern


    !>-----------------------------------------
    !> Energy matrix of cell i, A_i P^T S_c P, added into the CSR values.
    !>   local unknowns: u_i and the face values u_s of the faces of i
    !>   g      = 1/|K| sum_s |s| (u_s - u_i) n_s                    (exact for affine u)
    !>   R_s    = u_s - u_i - g.(x_s - x_i)
    !>   S      = sum_s |s| d_s / 3 | g + n_s R_s stab / max(d_s, cap h_i) |^2
    !>   boundary face values (free surfaces) are eliminated by minimising the energy (Schur
    !>   complement), interior face values are their interpolations u_s = sum_j beta_j u_j
    !---------------------------------------------------------------------------
    subroutine cellMatrix( i, cen, vol, A0n, stab, cfPtr, cfIdx, fK, fL, fA, fX, fN, fD, fs, U, &
                           rowPtr, colIdx, vals, nNeg )
    integer, intent(in) :: i, cfPtr(:), cfIdx(:), fK(:), fL(:), U(:), rowPtr(:), colIdx(:)
    real(dp), intent(in) :: cen(:,:), vol(:), A0n(:), stab, fA(:), fX(:,:), fN(:,:), fD(:,:)
    type(faceStencil), intent(in) :: fs(:)
    real(dp), intent(inout) :: vals(:)
    integer, intent(inout) :: nNeg

    integer :: nf, a, b, k, f, sf, nk, nbd, im, nU, ir, row, col, pos
    integer, allocatable :: keep(:), bnd(:)
    real(dp), allocatable :: G(:,:), S(:,:), rr(:,:), nn(:,:), ar(:), dd(:), Rv(:), Bm(:,:)
    real(dp), allocatable :: Sc(:,:), Sbb(:,:), Skb(:,:), P(:,:), Mm(:,:)
    real(dp) :: hK, deff, v

        nf = cfPtr(i+1) - cfPtr(i)
        nU = size(U)
        allocate( G(3,0:nf), S(0:nf,0:nf), rr(3,nf), nn(3,nf), ar(nf), dd(nf), Rv(0:nf), Bm(3,0:nf) )
        allocate( keep(nf+1), bnd(nf) )
        nk = 1
        keep(1) = 0
        nbd = 0
        do a = 1, nf
            sf = cfIdx(cfPtr(i)+a-1)
            f = abs(sf)
            if ( sf .gt. 0 ) then
                nn(:,a) = fN(:,f)
                rr(:,a) = fX(:,f) - cen(:,i)
            else
                nn(:,a) = -fN(:,f)
                rr(:,a) = fX(:,f) - cen(:,fK(f)) - fD(:,f)
            endif
            ar(a) = fA(f)
            dd(a) = dot_product( nn(:,a), rr(:,a) )
            if ( dd(a) .le. 0._dp ) nNeg = nNeg + 1
            if ( fL(f) .eq. 0 ) then
                nbd = nbd + 1
                bnd(nbd) = a
            else
                nk = nk + 1
                keep(nk) = a
            endif
        enddo
        hK = maxval( norm2( rr, dim=1 ) )

        G = 0._dp
        do a = 1, nf
            G(:,a) = ar(a) * nn(:,a) / vol(i)
            G(:,0) = G(:,0) - G(:,a)
        enddo
        S = 0._dp
        do a = 1, nf
            Rv = -matmul( rr(:,a), G )
            Rv(a) = Rv(a) + 1._dp
            Rv(0) = Rv(0) - 1._dp
            deff = max( dd(a), capStab * hK )
            do b = 0, nf
                Bm(:,b) = G(:,b) + nn(:,a) * Rv(b) * ( stab / deff )
            enddo
            S = S + ( ar(a) * max( dd(a), 0._dp ) / 3._dp ) * matmul( transpose(Bm), Bm )
        enddo

        !static condensation of the free-surface face values
        allocate( Sc(nk,nk) )
        Sc = S(keep(1:nk), keep(1:nk))
        if ( nbd .gt. 0 ) then
            allocate( Sbb(nbd,nbd), Skb(nk,nbd) )
            Sbb = S(bnd(1:nbd), bnd(1:nbd))
            Skb = S(keep(1:nk), bnd(1:nbd))
            call symPinv( Sbb, nbd, 1.e-12_dp )
            Sc = Sc - matmul( Skb, matmul( Sbb, transpose(Skb) ) )
        endif

        !kept unknowns in terms of cell values
        allocate( P(nk,nU), Mm(nU,nU) )
        P = 0._dp
        P(1,1) = 1._dp
        do ir = 2, nk
            a = keep(ir)
            f = abs( cfIdx(cfPtr(i)+a-1) )
            do im = 1, size( fs(f)%idx )
                do k = 1, nU
                    if ( U(k) .eq. fs(f)%idx(im) ) exit
                enddo
                P(ir,k) = P(ir,k) + fs(f)%w(im)
            enddo
        enddo
        Mm = A0n(i) * matmul( transpose(P), matmul( Sc, P ) )
        Mm = 0.5_dp * ( Mm + transpose(Mm) )

        do a = 1, nU
            row = U(a)
            do b = 1, nU
                col = U(b)
                pos = findCol( rowPtr(row), rowPtr(row+1) - 1, col )
                v = Mm(a,b)
                !$OMP ATOMIC
                vals(pos) = vals(pos) + v
            enddo
        enddo

    contains

        integer function findCol( lo_, hi_, c ) result( p )
        integer, intent(in) :: lo_, hi_, c
        integer :: l, h, mid
            l = lo_
            h = hi_
            p = -1
            do while ( l .le. h )
                mid = ( l + h ) / 2
                if ( colIdx(mid) .eq. c ) then
                    p = mid
                    return
                elseif ( colIdx(mid) .lt. c ) then
                    l = mid + 1
                else
                    h = mid - 1
                endif
            enddo
        end function findCol

    end subroutine cellMatrix


    !>-----------------------------------------
    !> K = (K + K^T)/2 entry by entry (removes the round-off of the parallel summation), the
    !> diagonal set to minus the sum of the off-diagonal entries of its row (K 1 = 0 exactly),
    !> then the exchange matrix L = -V^-1 K
    !---------------------------------------------------------------------------
    subroutine symmetriseAndScale( N, vol, rowPtr, colIdx, vals )
    integer, intent(in) :: N, rowPtr(:), colIdx(:)
    real(dp), intent(in) :: vol(:)
    real(dp), intent(inout) :: vals(:)
    integer :: a, k, b, l, h, mid, pos
    real(dp) :: avg, offs
        do a = 1, N
            do k = rowPtr(a), rowPtr(a+1) - 1
                b = colIdx(k)
                if ( b .le. a ) cycle
                l = rowPtr(b)
                h = rowPtr(b+1) - 1
                pos = -1
                do while ( l .le. h )
                    mid = ( l + h ) / 2
                    if ( colIdx(mid) .eq. a ) then
                        pos = mid
                        exit
                    elseif ( colIdx(mid) .lt. a ) then
                        l = mid + 1
                    else
                        h = mid - 1
                    endif
                enddo
                avg = 0.5_dp * ( vals(k) + vals(pos) )
                vals(k) = avg
                vals(pos) = avg
            enddo
        enddo
        do a = 1, N
            offs = 0._dp
            pos = -1
            do k = rowPtr(a), rowPtr(a+1) - 1
                if ( colIdx(k) .eq. a ) then
                    pos = k
                else
                    offs = offs + vals(k)
                endif
            enddo
            if ( pos .gt. 0 ) vals(pos) = -offs
            do k = rowPtr(a), rowPtr(a+1) - 1
                vals(k) = -vals(k) / vol(a)
            enddo
        enddo
    end subroutine symmetriseAndScale


    !> Exchange stiffness of a pair: interface value for two phases when given, else harmonic mean
    pure function pairExchange( i, j, A0n, phase, Aint ) result( Aij )
    integer, intent(in) :: i, j, phase(:)
    real(dp), intent(in) :: A0n(:), Aint(:,:)
    real(dp) :: Aij
        Aij = -1._dp
        if ( size(Aint,1) .gt. 1 ) then
            if ( phase(i) .ne. phase(j) ) Aij = Aint( phase(i), phase(j) )
        endif
        if ( Aij .lt. 0._dp ) then
            if ( A0n(i) + A0n(j) .gt. 0._dp ) then
                Aij = 2._dp * A0n(i) * A0n(j) / ( A0n(i) + A0n(j) )
            else
                Aij = 0._dp
            endif
        endif
    end function pairExchange


    pure function cross( a, b ) result( c )
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
        c = [ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) ]
    end function cross


    pure subroutine sort3( t )
    integer, intent(inout) :: t(3)
    integer :: s
        if ( t(1) .gt. t(2) ) then; s = t(1); t(1) = t(2); t(2) = s; endif
        if ( t(2) .gt. t(3) ) then; s = t(2); t(2) = t(3); t(3) = s; endif
        if ( t(1) .gt. t(2) ) then; s = t(1); t(1) = t(2); t(2) = s; endif
    end subroutine sort3


    !------------------------------------------------------------------ hash table
    subroutine tableInit( t, nItems )
    type(keyTable), intent(out) :: t
    integer, intent(in) :: nItems
        t%nslot = 16
        do while ( t%nslot .lt. 2 * max( nItems, 1 ) )
            t%nslot = 2 * t%nslot
        enddo
        allocate( t%key(3,t%nslot), t%val(t%nslot), t%used(t%nslot) )
        t%used = .false.
    end subroutine tableInit

    pure function tableHash( t, k ) result( h )
    type(keyTable), intent(in) :: t
    integer(i8), intent(in) :: k(3)
    integer :: h
    integer(i8), parameter :: P = 1048573_i8
        h = int( modulo( modulo(k(1),P) * 73856093_i8 + modulo(k(2),P) * 19349663_i8 &
                         + modulo(k(3),P) * 83492791_i8, int(t%nslot,i8) ) ) + 1
    end function tableHash

    function tableFind( t, k ) result( slot )
    type(keyTable), intent(in) :: t
    integer(i8), intent(in) :: k(3)
    integer :: slot, h
        h = tableHash( t, k )
        slot = 0
        do while ( t%used(h) )
            if ( all( t%key(:,h) .eq. k ) ) then
                slot = h
                return
            endif
            h = h + 1
            if ( h .gt. t%nslot ) h = 1
        enddo
    end function tableFind

    subroutine tableInsert( t, k, v )
    type(keyTable), intent(inout) :: t
    integer(i8), intent(in) :: k(3)
    integer, intent(in) :: v
    integer :: h
        h = tableHash( t, k )
        do while ( t%used(h) )
            h = h + 1
            if ( h .gt. t%nslot ) h = 1
        enddo
        t%used(h) = .true.
        t%key(:,h) = k
        t%val(h) = v
    end subroutine tableInsert

end module CellHessianExchange
