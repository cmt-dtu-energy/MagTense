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

    !> Open addressing hash table for triples of 64 bit integer keys
    type keyTable
        integer :: nslot = 0
        integer(i8), allocatable :: key(:,:)
        integer, allocatable :: val(:)
        logical, allocatable :: used(:)
    end type keyTable

    public :: computeExchangeCellHessian

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

    integer :: N, M, i, k, a, b, j, v, c, slot, nCanon, f, nf, owner, nnz, kk
    integer, allocatable :: canon(:), ndCnt(:), ndPtr(:), ndItem(:), stamp(:)
    integer, allocatable :: fCount(:), fOwner(:), fLoc(:), fOther(:), Ec(:,:)
    real(dp) :: tol, emin, X(3), p(3,3), nrm(3), fc(3), d(3)
    integer(i8) :: q(3), kq(3), Kp(3)
    type(keyTable) :: tab, ftab
    integer :: tri(3), da, db, dc, s1, found

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
