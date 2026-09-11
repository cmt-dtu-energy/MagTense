module UTIL_MICROMAG

    use MicroMagParameters
    use MKL_SPBLAS
    use IO_GENERAL
    use ISO_C_BINDING
    
    implicit none
    
    contains

    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Converts the dense matrix D (size nx,ny) to a sparse matrix K (size nx,y) )
    !> @params[in] threshold a number specifying the lower limit of values in D that should be considered non-zero
    !> Double precision
    !>-----------------------------------------
    subroutine ConvertDenseToSparse_d( D, K, threshold)
    real(DP),dimension(:,:),intent(in) :: D                 !> Dense input matrix    
    real(SP),intent(in) :: threshold                          !> Values less than this (in absolute) of D are considered zero
    type(MagTenseSparse),intent(inout) :: K                 !> Sparse matrix allocation
    
    integer :: nx,ny, nnonzero
    
    logical,dimension(:,:),allocatable :: mask          !> mask used for finding non-zero values
    integer,dimension(:),allocatable :: colInds         !> Used for storing the values 1...n used for indexing
    integer :: i,j,ind,stat
    
    nx = size(D(:,1))
    ny = size(D(1,:))
    
    allocate(mask(nx,ny))
    
    !find all entries larger than the threshold
    mask = abs(D) .gt. threshold
    nnonzero = count( mask )
        
    allocate( K%values(nnonzero),K%cols(nnonzero))
    allocate( K%rows_start(nx), K%rows_end(nx), colInds(ny) )
    
    do i=1,ny
        colInds(i) = i
    enddo
    
    
    !loop over each row
    ind = 1
    do i=1,nx
        !find all non-zero elements in the i'th row of D
        !starting index of the i'th row
        K%rows_start(i) = ind
        do j=1,ny
            if ( mask(i,j) .eqv. .true. ) then
                K%values( ind ) = D(i,j)
                
                K%cols( ind ) = j
                
                ind = ind + 1
            endif
        enddo                                        
        !ending index of the i'th row
        K%rows_end(i) = ind
    enddo
    
    !make sparse matrix
    stat = mkl_sparse_s_create_csr ( K%A, SPARSE_INDEX_BASE_ONE, nx, ny, K%rows_start, K%rows_end, K%cols, K%values)
    
    !clean up
    deallocate(colInds)
    
    end subroutine ConvertDenseToSparse_d
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Converts the dense matrix D (size nx,ny) to a sparse matrix K (size nx,y) )
    !> @params[in] threshold a number specifying the lower limit of values in D that should be considered non-zero
    !> single precision
    !>-----------------------------------------
    subroutine ConvertDenseToSparse_s( D, K, threshold)
    real(SP),dimension(:,:),intent(in) :: D                 !> Dense input matrix    
    real(SP),intent(in) :: threshold                        !> Values less than this (in absolute) of D are considered zero
    type(MagTenseSparse),intent(inout) :: K                 !> Sparse matrix allocation
    
    integer :: nx,ny, nnonzero
    
    logical,dimension(:,:),allocatable :: mask          !> mask used for finding non-zero values
    integer,dimension(:),allocatable :: colInds         !> Used for storing the values 1...n used for indexing
    integer :: i,j,ind,stat
    character(10) :: prog_str
    
    nx = size(D(:,1))
    ny = size(D(1,:))
    
    allocate(mask(nx,ny))
    
    mask = abs(D) .gt. threshold
    
    nnonzero = count( mask )
    !call displayGUIMessage( 'Number of nonzero demag elements:' )
    !write (prog_str,'(I10.9)') nnonzero
    !call displayGUIMessage( prog_str )
    
    allocate( K%values(nnonzero),K%cols(nnonzero))
    allocate( K%rows_start(nx), K%rows_end(nx), colInds(ny) )
    
    do i=1,ny
        colInds(i) = i
    enddo
    
    
    !loop over each row
    ind = 1
    do i=1,nx
        !find all non-zero elements in the i'th row of D
        !starting index of the i'th row
        K%rows_start(i) = ind
        do j=1,ny
            if ( mask(i,j) .eqv. .true. ) then
                K%values( ind ) = D(i,j)
                
                K%cols( ind ) = j
                
                ind = ind + 1
            endif
        enddo                                        
        !ending index of the i'th row
        K%rows_end(i) = ind
    enddo
    
    !make sparse matrix
    stat = mkl_sparse_s_create_csr ( K%A, SPARSE_INDEX_BASE_ONE, nx, ny, K%rows_start, K%rows_end, K%cols, K%values)
    
    !clean up
    deallocate(colInds)
    
    end subroutine ConvertDenseToSparse_s
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Converts the dense matrix D (size nx,ny) to a sparse matrix K (size nx,y) )
    !> With the matrices being of type complex(kind=4)
    !> @params[in] threshold a number specifying the lower limit of values in D that should be considered non-zero
    !>-----------------------------------------
    subroutine ConvertDenseToSparse_c( D, K, threshold)
    complex(kind=4),dimension(:,:),intent(in) :: D          !> Dense input matrix    
    type(MagTenseSparse_c),intent(inout) :: K                 !> Sparse matrix allocation
    complex(kind=4),intent(in) :: threshold                        !> Values less than this (in absolute) of D are considered zero
    
    
    integer :: nx,ny, nnonzero
    
    logical,dimension(:,:),allocatable :: mask          !> mask used for finding non-zero values
    integer,dimension(:),allocatable :: colInds         !> Used for storing the values 1...n used for indexing
    integer :: i,j,ind,stat
    
    nx = size(D(:,1))
    ny = size(D(1,:))
    
    allocate(mask(nx,ny))
    
    mask = abs(D) .gt. abs(threshold)
    
    nnonzero = count( mask )
    
    allocate( K%values(nnonzero),K%cols(nnonzero))
    allocate( K%rows_start(nx), K%rows_end(nx), colInds(ny) )
    
    do i=1,ny
        colInds(i) = i
    enddo
    
    
    !loop over each row
    ind = 1
    do i=1,nx
        !find all non-zero elements in the i'th row of D
        !starting index of the i'th row
        K%rows_start(i) = ind
        do j=1,ny
            if ( mask(i,j) .eqv. .true. ) then
                K%values( ind ) = D(i,j)
                
                K%cols( ind ) = j
                
                ind = ind + 1
            endif
        enddo                                        
        !ending index of the i'th row
        K%rows_end(i) = ind
    enddo
    
    !make sparse matrix
    stat = mkl_sparse_c_create_csr ( K%A, SPARSE_INDEX_BASE_ONE, nx, ny, K%rows_start, K%rows_end, K%cols, K%values)
    
    !clean up
    deallocate(colInds)
    
    end subroutine ConvertDenseToSparse_c
    
    
    subroutine create_CSR_matrix(rows, columns, values, K, N, eps_criteria, CSR_matrix)
        implicit none
        integer, dimension(:), intent(in) :: rows, columns
        real(dp), dimension(:), intent(in) :: values
        integer, intent(in) :: K, N
        real(dp), intent(in) :: eps_criteria
        type(sparse_matrix_t), intent(out) :: CSR_matrix
        integer :: nnz, stat
        integer, dimension(:), allocatable :: rows_reduced_COO, columns_reduced_COO
        real(dp), dimension(:), allocatable :: values_reduced_COO
        logical, allocatable :: mask1D(:)
        type(sparse_matrix_t) :: COO_matrix
    
        !Find the non-zero values of values
        allocate(mask1D(size(values)))    
        mask1D(:) = .false.
        mask1D = (abs(values) .gt. eps_criteria)
        
        !Then reduce the size of the arrays
        rows_reduced_COO    = pack(rows,mask1D)
        columns_reduced_COO = pack(columns,mask1D)
        values_reduced_COO  = pack(values,mask1D)
        deallocate(mask1D)


    

        ! Create a sparse matrix in COO format from the reduced arrays
        ! See https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-1/mkl-sparse-create-coo.html
        nnz = size(values_reduced_COO)

        stat = mkl_sparse_d_create_coo (COO_matrix, SPARSE_INDEX_BASE_ONE, K, N, nnz, rows_reduced_COO, columns_reduced_COO, values_reduced_COO)

        ! Create a sparse matrix in CSR format from COO format
        ! See https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2025-1/mkl-sparse-convert-csr.html
        stat = mkl_sparse_convert_csr (COO_matrix, SPARSE_OPERATION_NON_TRANSPOSE, CSR_matrix)
    
        stat = mkl_sparse_destroy(COO_matrix)
        
    end subroutine create_CSR_matrix


    subroutine create_COO_values_from_CSR(CSR_matrix, GridInfo)
        implicit none

        type(sparse_matrix_t), intent(inout) :: CSR_matrix
        type(MicroMagGridInfo), intent(inout) :: GridInfo
        type(sparse_matrix_t) :: CSR_copy_matrix
        type(sparse_matrix_t) :: COO_matrix
        integer :: N, K, stat, nnz, indexing
        type(MATRIX_DESCR) :: descr
        type(C_PTR)    :: rows_c, cols_c, values_c
        integer, POINTER :: rows(:), cols(:)
        real(dp), POINTER :: values(:)
        
        descr%type = SPARSE_MATRIX_TYPE_GENERAL 

        !Copy the CSR matrix, convert it to COO, then export this
        stat = mkl_sparse_copy (CSR_matrix, descr, CSR_copy_matrix)
    
        stat = mkl_sparse_convert_coo (CSR_copy_matrix, SPARSE_OPERATION_NON_TRANSPOSE, COO_matrix)
    
        stat = mkl_sparse_destroy (CSR_copy_matrix)
    
        stat = mkl_sparse_d_export_coo (COO_matrix, indexing, K, N, nnz, rows_c, cols_c, values_c)
        
        !   Converting C into Fortran pointers
        call C_F_POINTER(rows_c  , rows  , [nnz])
        call C_F_POINTER(cols_c  , cols  , [nnz])
        call C_F_POINTER(values_c    , values    , [nnz])
    
        ! Save the COO matrix in GridInfo
        allocate(GridInfo%Exch_mat_r(nnz),GridInfo%Exch_mat_c(nnz),GridInfo%Exch_mat_v(nnz))
        GridInfo%Exch_mat_nr = K
        GridInfo%Exch_mat_nc = N
        GridInfo%Exch_mat_ntot = size(rows)
        GridInfo%Exch_mat_r(:) = rows(:)
        GridInfo%Exch_mat_c(:) = cols(:)
        GridInfo%Exch_mat_v(:) = values(:)

        stat = mkl_sparse_destroy (COO_matrix)

            
    end subroutine create_COO_values_from_CSR


    !>-------------------------------------------
    !> Writes the sparse matrix A to disk
    !> @params[in] A sparse matrix according to the Intel MKL format
    !> @params[in] n the no. of rows and columns in A (assumed to be square)
    !>-------------------------------------------
    subroutine writeSparseMatrixToDisk( A, n, filename )
    use, intrinsic :: iso_c_binding
    type(sparse_matrix_t),intent(in) :: A
    integer,intent(in) :: n
    character(*),intent(in) :: filename
    
    real(SP),dimension(:,:),allocatable :: x,y      !> The dense matrix to be created
    real(c_float),dimension(:), allocatable :: xx,yy     
    integer,dimension(1) :: shp1D
    integer,dimension(2) :: shp2D
    integer :: i, stat
    type(matrix_descr) :: descr 
    real(c_float), parameter :: alpha_c = 1.
    real(c_float), parameter :: beta_c = 0. 
    
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    descr%mode = SPARSE_FILL_MODE_FULL
    descr%diag = SPARSE_DIAG_NON_UNIT
    
    allocate( x(n,n), y(n,n) )
    allocate( xx(n*n), yy(n*n) )
    
    x(:,:) = 0.
    xx(:) = 0.
    
    !make unit matrix
    do i=1,n
        x(i,i) = 1.
    enddo
    
    shp1D(1) = n*n
    
    xx = reshape( x, shp1D )
    
    stat = mkl_sparse_s_mm (SPARSE_OPERATION_NON_TRANSPOSE, alpha_c, A, descr, SPARSE_LAYOUT_COLUMN_MAJOR, xx, n, n, beta_c, yy, n)
    
    shp2D(1) = n
    shp2D(2) = n
    y = reshape( yy, shp2D )
    
    open (11, file=filename,	&
			           status='unknown', form='unformatted',	&
			           access='direct', recl=2*n*n)

    write(11,rec=1) y
    
    close(11)
    
    deallocate(x,xx,y,yy)
    
    end subroutine writeSparseMatrixToDisk
    
    
    !>-------------------------------------------
    !> Writes the sparse matrix A to disk
    !> @params[in] A sparse matrix according to the Intel MKL format
    !> @params[in] n the no. of rows and columns in A (assumed to be square)
    !>-------------------------------------------
    subroutine writeSparseMatrixToDisk_c( A, n, filename )
    use, intrinsic :: iso_c_binding
    type(sparse_matrix_t),intent(in) :: A
    integer,intent(in) :: n
    character(*),intent(in) :: filename
    
    complex(kind=4),dimension(:,:),allocatable :: x,y      !> The dense matrix to be created
    complex(c_float_complex),dimension(:), allocatable :: xx,yy     
    integer,dimension(1) :: shp1D
    integer,dimension(2) :: shp2D
    integer :: i, stat
    type(matrix_descr) :: descr 
    complex(c_float_complex), parameter :: alpha_c = cmplx(1.,0.)
    complex(c_float_complex), parameter :: beta_c = cmplx(0.,0.)
    
    
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    descr%mode = SPARSE_FILL_MODE_FULL
    descr%diag = SPARSE_DIAG_NON_UNIT
    
    allocate( x(n,n), y(n,n) )
    allocate( xx(n*n), yy(n*n) )
    
    x(:,:) = cmplx(0.,0.)
    xx(:) = cmplx(0.,0.)
    
    !make unit matrix
    do i=1,n
        x(i,i) = cmplx(1.,0.)
    enddo
    
    shp1D(1) = n*n
    
    xx = reshape( x, shp1D )
    
    stat = mkl_sparse_c_mm (SPARSE_OPERATION_NON_TRANSPOSE, alpha_c, A, descr, SPARSE_LAYOUT_COLUMN_MAJOR, xx, n, n, beta_c, yy, n)
    
    shp2D(1) = n
    shp2D(2) = n
    y = reshape( yy, shp2D )
    
    open (11, file=filename,	&
			           status='unknown', form='unformatted',	&
			           access='direct', recl=1*n*n)

    write(11,rec=1) real(y)
    write(11,rec=2) aimag(y)
    
    close(11)
    
    deallocate(x,xx,y,yy)
    
    end subroutine writeSparseMatrixToDisk_c
    
end module UTIL_MICROMAG