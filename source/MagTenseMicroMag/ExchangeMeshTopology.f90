!>-----------------------------------------
!> Mesh topology of unstructured grids (tetrahedra and prisms) for the cell-centred exchange
!> operators (ExchangeCellHessian, ExchangeSushi): cell geometry, neighbours with their periodic
!> offsets, boundary faces, the face list with interior/free-surface faces, and the faces and
!> face neighbours of every cell.
!>
!> Neighbours.
!>   Tetrahedra: every cell that shares a node with i. Nodes are identified by their coordinates
!>   (tolerance 1e-6 of the shortest edge; the neighbouring keys are searched too, so rounding of
!>   the key cannot split two coinciding nodes), so duplicated nodes and the nodes of the two
!>   planes of a periodic direction are merged.
!>   Prisms: every cell whose closed box touches the box of i (faces, edges, corners, also across
!>   hanging faces of a refined grid and for arbitrary cell widths), found with a bucket grid.
!>   The offset d_ij of every neighbour is the actual one through the shared node (tetrahedra) or
!>   the touching periodic image (prisms), never a minimum image, so it is right for any period.
!>
!> Boundary faces: a tetrahedron face owned by one live cell, or the part of a prism face not
!> covered by live neighbours. Faces towards cells with A0 = 0 (non-magnetic) are free surfaces.
!> Periodic coupling only along the directions with exchPBC set (period = extent of the mesh).
!---------------------------------------------------------------------------
module ExchangeMeshTopology
    use MKL_SPBLAS
    use MicroMagParameters
    use IO_GENERAL
    use UTIL_MICROMAG
    implicit none
    private

    integer, parameter :: i8 = selected_int_kind(18)

    !> Open addressing hash table for triples of 64 bit integer keys
    type keyTable
        integer :: nslot = 0
        integer(i8), allocatable :: key(:,:)
        integer, allocatable :: val(:)
        logical, allocatable :: used(:)
    end type keyTable

    public :: cellGeometry, tetNeighboursAndFaces, prismNeighboursAndFaces, tetFaceList, prismFaceList, &
              cellFaces, groupByOwner, pushInt, pushVec, cross, sort3

    contains

    !>-----------------------------------------
    !> Cells of an unstructured grid: their number, the live cells (A0 > 0), centroids, volumes,
    !> the extent of the mesh and the period along the exchPBC directions (zero elsewhere)
    !---------------------------------------------------------------------------
    subroutine cellGeometry( grid, pbc, A0n, N, live, cen, vol, lo, hi, per )
    type(MicroMagGrid), intent(in) :: grid
    logical, intent(in) :: pbc(3)
    real(dp), intent(in) :: A0n(:)
    integer, intent(out) :: N
    logical, allocatable, intent(out) :: live(:)
    real(dp), allocatable, intent(out) :: cen(:,:), vol(:)
    real(dp), intent(out) :: lo(3), hi(3), per(3)
    real(dp) :: p(3,4)
    integer :: i

        if ( grid%gridType .eq. gridTypeTetrahedron ) then
            N = size(grid%elements, 2)
        else
            N = size(grid%pts, 1)
        endif
        allocate( live(N) )
        live = A0n(1:N) .gt. 0._dp

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

    end subroutine cellGeometry


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

end module ExchangeMeshTopology
