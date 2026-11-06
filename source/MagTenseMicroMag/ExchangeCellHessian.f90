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
!> Stencil: the neighbours of i (see ExchangeMeshTopology); a stencil with fewer than minStencil
!> cells is grown by the neighbours of its neighbours.
!>
!> Neighbours, boundary faces and cell geometry come from ExchangeMeshTopology.
!---------------------------------------------------------------------------
module ExchangeCellHessian
    use MKL_SPBLAS
    use MicroMagParameters
    use IO_GENERAL
    use UTIL_MICROMAG
    use ExchangeMeshTopology
    implicit none
    private

    integer, parameter :: minStencil = 12
    integer, parameter :: nUnk = 9
    real(dp), parameter :: weightExp = 2._dp
    real(dp), parameter :: rcondFit = 1.e-10_dp

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

    real(dp), allocatable :: cen(:,:), vol(:), bfC(:,:), bfN(:,:), nbOff(:,:), stOff(:,:), cooV(:)
    integer, allocatable :: nbPtr(:), nbIdx(:), bfOwner(:), bfPtr(:), bfIdx(:)
    integer, allocatable :: stPtr(:), stIdx(:), cooR(:), cooC(:)
    logical, allocatable :: live(:)
    real(dp) :: per(3), lo(3), hi(3)
    integer :: N, nGrown, nRankDef, nBad, stMin, stMax, nbf
    character(len=256) :: msg

        !--- live cells, centroids, extent and periods
        call cellGeometry( grid, pbc, A0n, N, live, cen, vol, lo, hi, per )

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

end module ExchangeCellHessian
