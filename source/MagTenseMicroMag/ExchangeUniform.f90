!>-----------------------------------------
!> Exchange operator of a uniform grid: the 7-point finite-difference Laplacian with the
!> position dependent exchange stiffness of Heistracher et al. (ComputeExchangeTerm3D_Uniform),
!> and the zero operator used when there is no exchange at all (ComputeZeroExchangeOperator).
!---------------------------------------------------------------------------
module ExchangeUniform
    use MKL_SPBLAS
    use MicroMagParameters
    use IO_GENERAL
    use UTIL_MICROMAG
    use trace_mod
    implicit none
    private

    public :: ComputeExchangeTerm3D_Uniform, ComputeZeroExchangeOperator

    contains

    !>-----------------------------------------
    !> @author Rasmus Bjørk, rabj@dtu.dk, DTU, 2026
    !> @brief
    !> Calculates the exhange term matrix on a uniform grid
    !> which means it produces the differential operator d^2/dx^2 + d^2/dy^2 + d^2/dz^2 and returns this in the sparse matrix A
    !> Includes position dependent exchange stiffness A0 through the modified stencil described in 
    !> Heistracher et al. "Proposal for a micromagnetic standard problem: Domain wall pinning at phase boundaries", DOI: 10.1016/j.jmmm.2021.168875
    !---------------------------------------------------------------------------   
    subroutine ComputeExchangeTerm3D_Uniform( grid, A, problem, solution )
    type(MicroMagGrid),intent(in) :: grid              !> Struct containing the grid information    
    type(sparse_matrix_t),intent(inout) :: A           !> The returned matrix from the sparse matrix creator
    type(MicroMagProblem),intent(inout) :: problem     !> Problem data structure
    type(MicroMagSolution),intent(inout) :: solution   !> Solution data structure
        
    integer :: stat                                    !> Status value for the various sparse matrix operations        
    type(MagTenseSparse_d) :: d2dx2, d2dy2, d2dz2      !> Sparse matrices for the double derivatives with respect to x, y and z, respectively.
    type(sparse_matrix_t) :: tmp                       !> Temporary sparse matrices used for internal calculations
    integer :: ind,ntot,nNonZero,colInd,rowInd         !> Internal counter for indexing, the total no. of elements in the current sparse matrix being manipulated
    integer :: i,j,k,nx,ny,nz                          !> For-loop counters
    type(matrix_descr) :: descr                        !> Describes a sparse matrix operation
    real(DP) :: const
    character*(100) :: prog_str 
    integer, save :: itimer = 0
    logical :: PBCx, PBCy, PBCz                        !> Whether to have periodic boundary conditions along x, y and z respectively
    
    call trace%begin( "ComputeExchangeTerm3D_Uniform", itimer=itimer, verbose=1 )
    
    !Find the three sparse matrices for the individual directions. Then add them to get the total matrix
    !It is assumed that the magnetization vector to operate on is in fact a single column of Mx, My or Mz respectively.

    !The sparse matrices are encoded using an older version of Intels Compressed Sparse Row (CSR) format
    !cols(i) and values(i) give the column index and value of the i'th non-zero element
    !rows_start(i) and rows_end(i)-1 give the first and last index corresponding to the i'th row
    !If rows_start(i) = rows_end(i), the row is empty
    !For instance the matrix [[1, 0, 2], [0, -1, 0], [0, 0, 0], [-2, 5, 4]] has cols = (1, 3, 2, 1, 2, 3),
    !values = (1, 2, -1, -2, 5, 4), rows_start = (1, 3, 4, 4) and rows_end = (3, 4, 4, 7)
    !Note how size(rows_start) and size(rows_end) equal the number of rows, while size(cols) and size(values) equal the number of non-zero elements
    
    nx = grid%nx
    ny = grid%ny
    nz = grid%nz
    ntot = nx*ny*nz     ! Total number of micromagnetic tiles
    PBCx = problem%macrogrid%exchPBC(1)
    PBCy = problem%macrogrid%exchPBC(2)
    PBCz = problem%macrogrid%exchPBC(3)

    !A periodic direction needs at least three cells. With exactly two, the wrap-around neighbour
    !and the ordinary neighbour of a boundary cell are the same cell, so the row would carry two
    !entries with the same column index, which is not a valid CSR row. The unstructured mesh
    !analysis already refuses this case; the uniform grid used to build the degenerate matrix
    !silently. A single cell along a direction is fine - the direction is simply dropped below,
    !which is the correct periodic answer for one cell.
    if ( ( PBCx .and. nx .eq. 2 ) .or. ( PBCy .and. ny .eq. 2 ) .or. ( PBCz .and. nz .eq. 2 ) ) then
        call displayGUIMessage( 'MagTense: periodic exchange requires at least three cells along a periodic direction' )
        call displayGUIMessage( 'MagTense: with two cells a tile becomes its own neighbour through the periodic image' )
        call abortMagTense( 'ComputeExchangeTerm3D_Uniform: too few cells along a periodic direction' )
    endif


    !----------------------------------d^2dx^2 begins -----------------------------!
    if ( nx .gt. 1 ) then
        ! Make the d^2/dx^2 matrix
        if ( PBCx ) then
            nNonZero = 3*ntot             ! For each tile, there are 3 non-zero elements corresponding to itself and the two neighbours along x
        else
            nNonZero = 3*ntot - 2*ny*nz   ! With a vacuum boundary, the 2*ny*nz tiles at the two end surfaces only have one neighbour each
        end if
        allocate(d2dx2%values(nNonZero),d2dx2%cols(nNonZero),d2dx2%rows_start(ntot),d2dx2%rows_end(ntot))
    
        ind = 1
        rowInd = 1
        colInd = 1
    
        do k=1,nz
            do j=1,ny
                ! The left boundary
                if ( PBCx ) then
                    ! Entries are emitted in ascending column order: the wrap-around neighbour of
                    ! the left boundary is the rightmost tile and therefore comes last, not first.
                    ! MKL's sparse routines take sorted column indices within a row as given.

                    ! Current tile
                    d2dx2%values(ind) = -( exchCoeff( problem, 1,j,k, nx,j,k ) + exchCoeff( problem, 1,j,k, 2,j,k ) )
                    d2dx2%cols(ind) = colInd
                    d2dx2%rows_start(rowInd) = ind
                    ind = ind + 1

                    ! Coupling to the right
                    d2dx2%values(ind) = exchCoeff( problem, 1,j,k, 2,j,k )
                    d2dx2%cols(ind) = colInd + 1
                    ind = ind + 1

                    ! Coupling to the left (loops around and connects with rightmost tile)
                    d2dx2%values(ind) = exchCoeff( problem, 1,j,k, nx,j,k )
                    d2dx2%cols(ind) = colInd + nx - 1
                    d2dx2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                    ind = ind + 1
                else
                    !The left boundary
                    const = exchCoeff( problem, 1,j,k, 2,j,k )

                    d2dx2%values(ind) = -const
                    d2dx2%cols(ind) = colInd
                    d2dx2%rows_start(rowInd) = ind            
                    ind = ind + 1
                
                    d2dx2%values(ind) = const
                    d2dx2%cols(ind) = colInd + 1
                    d2dx2%rows_end(rowInd) = ind+1
                    rowInd = rowInd + 1
                    ind = ind + 1
                end if
                        
            
                !Go through one row at a time
                do i=2,nx-1
                            
                    !Left-most point
                    d2dx2%values( ind ) = exchCoeff( problem, i,j,k, i-1,j,k )
                    
                    d2dx2%cols(ind) = colInd
                    !update where the row starts
                    d2dx2%rows_start(rowInd) = ind                           
                    ind = ind + 1
                
                    !Center point
                    d2dx2%values( ind ) = -( exchCoeff( problem, i,j,k, i-1,j,k ) + exchCoeff( problem, i,j,k, i+1,j,k ) )
                    d2dx2%cols(ind) = colInd + 1
                    ind = ind + 1
                
                    !Right-most point
                    d2dx2%values( ind ) = exchCoeff( problem, i,j,k, i+1,j,k )
                    d2dx2%cols(ind) = colInd + 2
                    d2dx2%rows_end(rowInd) = ind+1
                    rowInd = rowInd + 1
                    ind = ind + 1

                    colInd = colInd + 1
                enddo            

                ! The right boundary
                if ( PBCx ) then
                    ! Ascending column order: here the wrap-around neighbour is the leftmost tile,
                    ! so it is the smallest column index and comes first.

                    ! Coupling to the right (loops around and connects with leftmost tile)
                    d2dx2%values(ind) = exchCoeff( problem, nx,j,k, 1,j,k )
                    d2dx2%cols(ind) = colInd - nx + 2
                    ! Update where the row starts
                    d2dx2%rows_start(rowInd) = ind
                    ind = ind + 1

                    ! Coupling to the left
                    d2dx2%values(ind) = exchCoeff( problem, nx,j,k, nx-1,j,k )
                    d2dx2%cols(ind) = colInd
                    ind = ind + 1

                    ! Current tile
                    d2dx2%values(ind) = -( exchCoeff( problem, nx,j,k, nx-1,j,k ) + exchCoeff( problem, nx,j,k, 1,j,k ) )
                    d2dx2%cols(ind) = colInd + 1
                    d2dx2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                    ind = ind + 1

                    colInd = colInd + 2
                else
                    const = exchCoeff( problem, nx,j,k, nx-1,j,k )
                    !The right boundary
                    d2dx2%values(ind) = const 
                    d2dx2%cols(ind) = colInd
                    !update where the row starts
                    d2dx2%rows_start(rowInd) = ind            
                    ind = ind + 1
            
                    d2dx2%values(ind) = -const
                    d2dx2%cols(ind) = colInd+1
                    d2dx2%rows_end(rowInd) = ind+1
                    rowInd = rowInd + 1
                    ind = ind + 1     
            
                    colInd = colInd + 2
                endif
            enddo
        enddo
    
        !Multiply by the discretization
        d2dx2%values = d2dx2%values * 1.0/grid%dx**2
        
        ! Create the sparse matrix for the d^2dx^2
        stat = mkl_sparse_d_create_csr ( d2dx2%A, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, d2dx2%rows_start, d2dx2%rows_end, d2dx2%cols, d2dx2%values)
    
    endif
    
    !----------------------------------d^2dx^2 ends ---------------------------- -!
    
    
    !----------------------------------d^2dy^2 begins ----------------------------!
    if ( ny .gt. 1 ) then
        ! Make the d^2/dy^2 matrix
        if ( PBCy ) then
            nNonZero = 3*ntot             ! For each tile, there are 3 non-zero elements corresponding to itself and the two neighbours along y
        else
            nNonZero = 3*ntot - 2*nx*nz   ! With a vacuum boundary, the 2*nx*nz tiles at the two end surfaces only have one neighbour each
        end if
        allocate(d2dy2%values(nNonZero),d2dy2%cols(nNonZero),d2dy2%rows_start(ntot),d2dy2%rows_end(ntot))
    
        ind = 1
        rowInd = 1
        colInd = 1
        do k=1,nz
            !The bottom boundary
             ! The backwards boundary (Imagine standing on the origin, eyes pointed along the positive y-axis)
            if ( PBCy ) then
                do i=1,nx
                    ! Ascending column order: the wrap-around neighbour is the forwardmost tile and
                    ! therefore has the largest column index, so it is emitted last.

                    ! Current tile
                    d2dy2%values(ind) = -( exchCoeff( problem, i,1,k, i,ny,k ) + exchCoeff( problem, i,1,k, i,2,k ) )
                    d2dy2%cols(ind) = colInd
                    d2dy2%rows_start(rowInd) = ind
                    ind = ind + 1  ! Increment to next element

                    ! Forwards coupling
                    d2dy2%values(ind) = exchCoeff( problem, i,1,k, i,2,k )
                    d2dy2%cols(ind) = colInd + nx
                    ind = ind + 1  ! Increment to next element

                    ! Backwards coupling (loops around and couples to forwardmost tile)
                    d2dy2%values(ind) = exchCoeff( problem, i,1,k, i,ny,k )
                    d2dy2%cols(ind) = colInd + nx*(ny-1)
                    d2dy2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                    ind = ind + 1  ! Increment to next element
                    colInd = colInd + 1
                end do
            else
            
                do i=1,nx
                    const = exchCoeff( problem, i,1,k, i,2,k )

                    d2dy2%values(ind) = -const
                    d2dy2%cols(ind) = colInd
                    d2dy2%rows_start(rowInd) = ind
            
                    !increment to next element
                    ind = ind + 1
            
                    d2dy2%values(ind) = const
                    d2dy2%cols(ind) = colInd + nx
                    d2dy2%rows_end(rowInd) = ind+1
                    rowInd = rowInd + 1
            
                    !increment to next element
                    ind = ind + 1
            
                    colInd = colInd + 1
                enddo
            endif
            
            ! Everything in between
            do j=2,ny-1
                do i=1,nx
                
                    !lower value
                    d2dy2%values(ind) = exchCoeff( problem, i,j,k, i,j-1,k )
                    d2dy2%cols(ind) = colInd-nx
                    d2dy2%rows_start(rowInd) = ind
                    ind = ind + 1  ! Increment to next element
                
                    !central value
                    d2dy2%values(ind) = -( exchCoeff( problem, i,j,k, i,j-1,k ) + exchCoeff( problem, i,j,k, i,j+1,k ) )
                    d2dy2%cols(ind) = colInd
                    ind = ind + 1  ! Increment to next element
            
            
                    !upper value
                    d2dy2%values(ind) = exchCoeff( problem, i,j,k, i,j+1,k )
                    d2dy2%cols(ind) = colInd + nx
                    d2dy2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                    ind = ind + 1  ! Increment to next element
                    colInd = colInd + 1
                end do
            end do
        
            !The forwards boundary
            if ( PBCy ) then
                do i=1,nx
                    ! Ascending column order: the wrap-around neighbour is the backmost tile and
                    ! therefore has the smallest column index, so it is emitted first.

                    ! Forwards coupling (loops around and couples to backmost tile)
                    d2dy2%values(ind) = exchCoeff( problem, i,ny,k, i,1,k )
                    d2dy2%cols(ind) = colInd - nx*(ny-1)
                    d2dy2%rows_start(rowInd) = ind
                    ind = ind + 1  ! Increment to next element

                    ! Backwards coupling
                    d2dy2%values(ind) = exchCoeff( problem, i,ny,k, i,ny-1,k )
                    d2dy2%cols(ind) = colInd - nx
                    ind = ind + 1  ! Increment to next element

                    ! Current tile
                    d2dy2%values(ind) = -( exchCoeff( problem, i,ny,k, i,ny-1,k ) + exchCoeff( problem, i,ny,k, i,1,k ) )
                    d2dy2%cols(ind) = colInd
                    d2dy2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                    ind = ind + 1  ! Increment to next element
                    colInd = colInd + 1
                enddo
            else
                do i=1,nx
                   const = exchCoeff( problem, i,ny,k, i,ny-1,k )
                    !lower element
                    d2dy2%values(ind) = const
                    d2dy2%cols(ind) = colInd - nx
                    d2dy2%rows_start(rowInd) = ind
            
                    !increment to next element
                    ind = ind + 1            
            
                    !central element
                    d2dy2%values(ind) = -const
                    d2dy2%cols(ind) = colInd
                    d2dy2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
            
                    ind = ind + 1
                    colInd = colInd + 1
                enddo
            endif
        enddo
                    
        !Multiply by the discretization
        d2dy2%values = d2dy2%values * 1./grid%dy**2
        
        ! Create the sparse matrix for the d^2dy^2
        stat = mkl_sparse_d_create_csr ( d2dy2%A, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, d2dy2%rows_start, d2dy2%rows_end, d2dy2%cols, d2dy2%values)
    endif
    
    !----------------------------------d^2dy^2 ends ----------------------------!
    
    
    !----------------------------------d^2dz^2 begins ----------------------------!
    if ( nz .gt. 1 ) then
        ! Make the d^2/dz^2 matrix
        if ( PBCz ) then
            nNonZero = 3*ntot             ! For each tile, there are 3 non-zero elements corresponding to itself and the two neighbours along z
        else
            nNonZero = 3*ntot - 2*nx*ny   ! With a vacuum boundary, the 2*nx*ny tiles at the two end surfaces only have one neighbour each
        end if
        allocate(d2dz2%values(nNonZero),d2dz2%cols(nNonZero),d2dz2%rows_start(ntot),d2dz2%rows_end(ntot))
    
        ind = 1
        rowInd = 1
        colInd = 1
        ! The bottom boundary
        if ( PBCz ) then
            do j=1,ny
                do i=1,nx
                    ! Ascending column order: the wrap-around neighbour is the topmost tile and
                    ! therefore has the largest column index, so it is emitted last.

                    ! Current tile
                    d2dz2%values(ind) = -( exchCoeff( problem, i,j,1, i,j,nz ) + exchCoeff( problem, i,j,1, i,j,2 ) )
                    d2dz2%cols(ind) = colInd
                    d2dz2%rows_start(rowInd) = ind
                    ind = ind + 1   ! Increment to next element

                    ! Upwards coupling
                    d2dz2%values(ind) = exchCoeff( problem, i,j,1, i,j,2 )
                    d2dz2%cols(ind) = nx*ny + colInd
                    ind = ind + 1   ! Increment to next element

                    ! Downwards coupling (loops around and couples to topmost tile)
                    d2dz2%values(ind) = exchCoeff( problem, i,j,1, i,j,nz )
                    d2dz2%cols(ind) = colInd + ny*nx*(nz-1)
                    d2dz2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                    ind = ind + 1   ! Increment to next element
                    colInd = colInd + 1
                end do
            end do
        else
            do j=1,ny
                do i=1,nx
                    const = exchCoeff( problem, i,j,1, i,j,2 )
                    !central value
                    d2dz2%values(ind) = -const
                    d2dz2%cols(ind) = colInd
                    d2dz2%rows_start(rowInd) = ind
            
                    !increment position
                    ind = ind + 1
            
                    !right-most value
                    d2dz2%values(ind) = const
                    d2dz2%cols(ind) = nx * ny + colInd
                    d2dz2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
            
                    !increment position
                    ind = ind + 1
                    colInd = colInd + 1
                end do
            end do
        endif

        ! Everything in between
        do k=2,nz-1            
            do j=1,ny
                do i=1,nx
                    !left-most value
                    d2dz2%values(ind) = exchCoeff( problem, i,j,k, i,j,k-1 )
                    d2dz2%cols(ind) = colInd - nx * ny
                    d2dz2%rows_start(rowInd) = ind
                    ind = ind + 1   ! Increment to next element
                
                    !central value
                    d2dz2%values(ind) = -( exchCoeff( problem, i,j,k, i,j,k-1 ) + exchCoeff( problem, i,j,k, i,j,k+1 ) )
                    d2dz2%cols(ind) = colInd
                    ind = ind + 1   ! Increment to next element
                
                    !right-most value
                    d2dz2%values(ind) = exchCoeff( problem, i,j,k, i,j,k+1 )
                    d2dz2%cols(ind) = colInd + nx * ny
                    d2dz2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                    ind = ind + 1   ! Increment to next element
                    colInd = colInd + 1
                end do
            end do
        end do

        ! The top boundary
        if ( PBCz ) then
            do j=1,ny
                do i=1,nx
                    ! Ascending column order: the wrap-around neighbour is the bottommost tile and
                    ! therefore has the smallest column index, so it is emitted first.

                    ! Upwards coupling (loops around and couples to bottommost tile)
                    d2dz2%values(ind) = exchCoeff( problem, i,j,nz, i,j,1 )
                    d2dz2%cols(ind) = colInd - nx*ny*(nz - 1)
                    d2dz2%rows_start(rowInd) = ind
                    ind = ind + 1   ! Increment to next element

                    ! Downwards coupling
                    d2dz2%values(ind) = exchCoeff( problem, i,j,nz, i,j,nz-1 )
                    d2dz2%cols(ind) = colInd - nx*ny
                    ind = ind + 1   ! Increment to next element

                    ! Current tile
                    d2dz2%values(ind) = -( exchCoeff( problem, i,j,nz, i,j,nz-1 ) + exchCoeff( problem, i,j,nz, i,j,1 ) )
                    d2dz2%cols(ind) = colInd
                    d2dz2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                    ind = ind + 1   ! Increment to next element
                    colInd = colInd + 1
                end do
            end do
        else
            do j=1,ny
                do i=1,nx
                    const = exchCoeff( problem, i,j,nz, i,j,nz-1 )

                    !left-most value
                    d2dz2%values(ind) = const
                    d2dz2%cols(ind) = colInd - nx * ny
                    d2dz2%rows_start(rowInd) = ind
                
                    !increment position
                    ind = ind + 1
                
                    !central value
                    d2dz2%values(ind) = -const
                    d2dz2%cols(ind) = colInd
                    d2dz2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                
                    !increment position
                    ind = ind + 1
                    colInd = colInd + 1
                end do
            end do
        endif
        
        !Multiply by the discretization
        d2dz2%values = d2dz2%values * 1./grid%dz**2
        
        !Create the sparse matrix for the d^2dz^2
        stat = mkl_sparse_d_create_csr ( d2dz2%A, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, d2dz2%rows_start, d2dz2%rows_end, d2dz2%cols, d2dz2%values)
    endif
    
    !----------------------------------d^2dz^2 ends ----------------------------!
    
            
    !Finally, add up the three diagonals and store in the output sparse matrix, A
    !Store the results temporarily in tmp4
    
    !call writeSparseMatrixToDisk( d2dz2%A, nx*ny*nz, 'd2dz2.dat' )
    !call writeSparseMatrixToDisk( d2dy2%A, nx*ny*nz, 'd2dy2.dat' )
    
    const = 1.
    
    if ( nx .gt. 1) then
        if ( ny .gt. 1) then
            stat = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, d2dx2%A, const, d2dy2%A, tmp)  
            deallocate(d2dy2%values,d2dy2%cols,d2dy2%rows_start,d2dy2%rows_end)
            stat = mkl_sparse_destroy (d2dy2%A)
            
            if ( nz .gt. 1) then
                !nx+ny+nz
                call displayGUIMessage( 'Exchange in x,y,z' )
                stat = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, d2dz2%A, const, tmp, A)       
                deallocate(d2dz2%values,d2dz2%cols,d2dz2%rows_start,d2dz2%rows_end)
                stat = mkl_sparse_destroy (d2dz2%A)
            else
                !nx+ny
                call displayGUIMessage( 'Exchange in x,y' )
                descr%type = SPARSE_MATRIX_TYPE_GENERAL
                descr%mode = SPARSE_FILL_MODE_FULL
                descr%diag = SPARSE_DIAG_NON_UNIT
                stat = mkl_sparse_copy ( tmp, descr, A )
                stat = mkl_sparse_destroy (tmp)
            endif
        else
            if ( nz .gt. 1) then
                !nx+nz
                call displayGUIMessage( 'Exchange in x,z' )
                stat = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, d2dx2%A, const, d2dz2%A, A)  
                deallocate(d2dz2%values,d2dz2%cols,d2dz2%rows_start,d2dz2%rows_end)
                stat = mkl_sparse_destroy (d2dz2%A)
            else
                !nx
                call displayGUIMessage( 'Exchange in x' )
                descr%type = SPARSE_MATRIX_TYPE_GENERAL
                descr%mode = SPARSE_FILL_MODE_FULL
                descr%diag = SPARSE_DIAG_NON_UNIT
                stat = mkl_sparse_copy ( d2dx2%A, descr, A )
            endif
        endif
        deallocate(d2dx2%values,d2dx2%cols,d2dx2%rows_start,d2dx2%rows_end)
        stat = mkl_sparse_destroy (d2dx2%A)
    else
        if ( ny .gt. 1) then
            if ( nz .gt. 1) then
                !ny+nz
                call displayGUIMessage( 'Exchange in y,z' )
                stat = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, d2dy2%A, const, d2dz2%A, A) 
                deallocate(d2dz2%values,d2dz2%cols,d2dz2%rows_start,d2dz2%rows_end)
                stat = mkl_sparse_destroy (d2dz2%A)
            else
                !ny
                call displayGUIMessage( 'Exchange in y' )
                descr%type = SPARSE_MATRIX_TYPE_GENERAL
                descr%mode = SPARSE_FILL_MODE_FULL
                descr%diag = SPARSE_DIAG_NON_UNIT
                stat = mkl_sparse_copy ( d2dy2%A, descr, A )
            endif
            deallocate(d2dy2%values,d2dy2%cols,d2dy2%rows_start,d2dy2%rows_end)
            stat = mkl_sparse_destroy (d2dy2%A)
        else
            if ( nz .gt. 1) then
                !nz
                call displayGUIMessage( 'Exchange in z' )
                descr%type = SPARSE_MATRIX_TYPE_GENERAL
                descr%mode = SPARSE_FILL_MODE_FULL
                descr%diag = SPARSE_DIAG_NON_UNIT
                stat = mkl_sparse_copy ( d2dz2%A, descr, A )
                deallocate(d2dz2%values,d2dz2%cols,d2dz2%rows_start,d2dz2%rows_end)
                stat = mkl_sparse_destroy (d2dz2%A)
            else
                !No exchange terms: a single cell along every direction, i.e. a macrospin. The rest of
                !the code expects an exchange operator, so build one with an explicit zero on every
                !diagonal. This used to fall through with d2dz2 unallocated and die in the deallocate.
                call displayGUIMessage( 'No exchange terms (single cell)' )
                call ComputeZeroExchangeOperator( ntot, A )
                problem%useExchange = .false.
            endif
        endif
    endif
    
    !call writeSparseMatrixToDisk( tmp, nx*ny*nz, 'A_exch.dat' )
    
    !call writeSparseMatrixToDisk( A, nx*ny*nz, 'A_exch.dat' )
    
    call create_COO_values_from_CSR(A,solution%gridinfo)
    

    call trace%end( "ComputeExchangeTerm3D_Uniform", itimer=itimer, verbose=1 )
    end subroutine ComputeExchangeTerm3D_Uniform


    !>-----------------------------------------
    !> @brief
    !> Pairwise exchange coefficient between a cell and one of its neighbours on a uniform
    !> grid, following Heistracher et al., DOI: 10.1016/j.jmmm.2021.168875.
    !>
    !> This is the harmonic mean of the two exchange stiffnesses divided by the stiffness of
    !> the cell owning the matrix row, 2*A_nb/(A_nb + A_self) = A_harm/A_self. The division by
    !> A_self is there because Jfact on a uniform grid is A0/(mu0*Ms) per cell, whereas the
    !> unstructured mesh normalises by the global A0_max instead and so uses A_harm directly.
    !>
    !> The expression used to be written out by hand at 39 places in the routine above, which
    !> made any change to it - such as an interface value specified per pair of materials -
    !> 39 chances to introduce an inconsistency that no test would localise.
    !> @param[in] A_self exchange stiffness of the cell owning the row
    !> @param[in] A_nb exchange stiffness of the neighbouring cell
    !>-----------------------------------------
    !> If the two cells belong to different materials and the user has given an interface
    !> value for that pair, that value replaces the harmonic mean. The coefficient is still
    !> divided by A_self, because Jfact on a uniform grid carries the per-cell A0.
    !>
    !> A cell with A0 = 0 is not magnetic, so a face touching one ignores the interface table
    !> entirely and keeps the ordinary formula. That already leaves the magnetic side
    !> uncoupled, since 2*A_nb/(A_nb + A_self) vanishes when the neighbour is dead, and it is
    !> the same rule the unstructured mesh follows, so the two grid types agree.
    pure function exchCoeff( problem, i_s, j_s, k_s, i_n, j_n, k_n ) result( coeff )
    type(MicroMagProblem),intent(in) :: problem        !> Problem, for A0_map and the phases
    integer,intent(in) :: i_s, j_s, k_s                !> Grid index of the cell owning the row
    integer,intent(in) :: i_n, j_n, k_n                !> Grid index of the neighbouring cell
    real(DP) :: coeff
    real(DP) :: A_self, A_nb, A_face
    integer :: p_self, p_nb

    A_self = problem%A0_map(i_s,j_s,k_s)
    A_nb   = problem%A0_map(i_n,j_n,k_n)

    A_face = -1.0_DP
    if ( problem%n_phase .gt. 1 .and. allocated(problem%A_int) &
         .and. allocated(problem%phase_map) ) then
        p_self = problem%phase_map(i_s,j_s,k_s)
        p_nb   = problem%phase_map(i_n,j_n,k_n)
        !A cell with A0 = 0 is not magnetic, so the override is simply not taken for a face
        !touching one and the ordinary formula below applies unchanged - which already gives
        !no coupling to the magnetic side, since 2*A_nb/(A_nb + A_self) vanishes when A_nb is
        !zero. Declining the override rather than forcing the coefficient to zero is what
        !makes the interface value have no effect at all on such a face, and it mirrors what
        !the unstructured mesh does. It also guarantees A_self > 0 wherever the override is
        !used, so the division below is safe.
        if ( p_self .ne. p_nb .and. A_self .gt. 0.0_DP .and. A_nb .gt. 0.0_DP ) &
            A_face = problem%A_int(p_self,p_nb)
    endif

    if ( A_face .ge. 0.0_DP ) then
        coeff = A_face / A_self
    else if ( ( A_nb + A_self ) .gt. 0.0_DP ) then
        coeff = 2 * A_nb / ( A_nb + A_self )
    else
        !Two neighbouring cells that are both non-magnetic. The harmonic mean is 0/0 there,
        !which with /fpe:0 is an invalid operation rather than a quiet NaN. Two cells with no
        !exchange are not exchange coupled, so the coefficient is zero.
        coeff = 0.0_DP
    endif

    end function exchCoeff


    !>-----------------------------------------
    !> @brief
    !> An exchange operator that is identically zero, for a problem without exchange: A0 = 0
    !> in every tile, or a uniform grid of one cell in total (a macrospin). A problem where
    !> only some tiles have A0 = 0 still has exchange between the others and does not come
    !> here. It is stored with an explicit zero on
    !> every diagonal rather than as an empty matrix, which keeps every row of the CSR arrays
    !> non-empty and gives the returned GridInfo a well-formed operator of the right size.
    !> @param[in] ntot the number of tiles
    !> @param[inout] A the returned sparse matrix
    !---------------------------------------------------------------------------
    subroutine ComputeZeroExchangeOperator( ntot, A )
    integer,intent(in) :: ntot                         !> Number of tiles
    type(sparse_matrix_t),intent(inout) :: A           !> The returned zero operator
    type(MagTenseSparse_d) :: Z                        !> Holds the CSR arrays until A owns a copy
    type(matrix_descr) :: descr
    integer :: i, stat

    allocate(Z%values(ntot),Z%cols(ntot),Z%rows_start(ntot),Z%rows_end(ntot))
    do i = 1, ntot
        Z%values(i) = 0.0_DP
        Z%cols(i) = i
        Z%rows_start(i) = i
        Z%rows_end(i) = i + 1
    end do
    !mkl_sparse_d_create_csr only references the arrays, so A is made as a copy before they go
    stat = mkl_sparse_d_create_csr ( Z%A, SPARSE_INDEX_BASE_ONE, ntot, ntot, Z%rows_start, Z%rows_end, Z%cols, Z%values)
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    descr%mode = SPARSE_FILL_MODE_FULL
    descr%diag = SPARSE_DIAG_NON_UNIT
    stat = mkl_sparse_copy ( Z%A, descr, A )
    stat = mkl_sparse_destroy ( Z%A )
    deallocate(Z%values,Z%cols,Z%rows_start,Z%rows_end)

    end subroutine ComputeZeroExchangeOperator

end module ExchangeUniform
