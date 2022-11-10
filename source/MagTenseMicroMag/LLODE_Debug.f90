#if USE_MATLAB
#include "fintrf.h"
#endif

    module LLODE_Debug
    use MKL_SPBLAS
    use MicroMagParameters
    implicit none
    
    contains
    
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
    
    
end module LLODE_Debug