
    module LLODE_Debug
    use MKL_SPBLAS
    implicit none
    
    contains
    
    !>-------------------------------------------
    !> Writes the sparse matrix A to disk
    !> @params[in] A sparse matrix according to the Intel MKL format
    !> @params[in] n the no. of rows and columns in A (assumed to be square)
    !>-------------------------------------------
    subroutine writeSparseMatrixToDisk( A, n, filename )
    type(sparse_matrix_t),intent(in) :: A
    integer,intent(in) :: n
    character(*),intent(in) :: filename
    
    real,dimension(:,:),allocatable :: x,y      !> The dense matrix to be created
    real,dimension(:), allocatable :: xx,yy     
    integer,dimension(1) :: shp1D
    integer,dimension(2) :: shp2D
    integer :: i, stat
    real :: alpha, beta
    type(matrix_descr) :: descr 
    
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
    
    alpha = 1.
    beta = 0.
    
    
    stat = mkl_sparse_d_mm (SPARSE_OPERATION_NON_TRANSPOSE, alpha, A, descr, SPARSE_LAYOUT_COLUMN_MAJOR, xx, n, n, beta, yy, n)
    
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
    
end module LLODE_Debug