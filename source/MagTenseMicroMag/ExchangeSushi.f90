!>-----------------------------------------
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

!>
!> Faces, face neighbours, vertex neighbours and cell geometry come from ExchangeMeshTopology.
!---------------------------------------------------------------------------
module ExchangeSushi
    use MKL_SPBLAS
    use MicroMagParameters
    use IO_GENERAL
    use UTIL_MICROMAG
    use ExchangeMeshTopology
    implicit none
    private

    !> SUSHI-g: stabilisation (boxes: identical to the 7-point Laplacian on uniform grids;
    !> tetrahedra: from the parameter study), its distance floor (fraction of the cell size), the
    !> largest face-neighbour stencil used whole and the conditioning bound of the interpolation
    real(dp), parameter :: stabTet = 1._dp
    real(dp), parameter :: stabBox = 1.7320508075688772_dp
    real(dp), parameter :: capStab = 0.3_dp
    real(dp), parameter :: bMax = 2.5_dp
    integer, parameter :: maxFull = 8

    !> Interpolation of a face value from cell values
    type faceStencil
        integer, allocatable :: idx(:)
        real(dp), allocatable :: w(:)
    end type faceStencil

    public :: computeExchangeSushi

    contains

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
    real(dp) :: per(3), lo(3), hi(3), stab
    integer :: N, i, f, nF, nbf, nBad, nNeg, k, nnz, nInt
    character(len=256) :: msg

        !--- live cells, centroids, volumes, extent and periods
        call cellGeometry( grid, pbc, A0n, N, live, cen, vol, lo, hi, per )
        if ( grid%gridType .eq. gridTypeTetrahedron ) then
            stab = stabTet
        else
            stab = stabBox
        endif

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

end module ExchangeSushi
