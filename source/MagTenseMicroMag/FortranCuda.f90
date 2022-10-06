
!>--------------------------
!> Module for interfacing with CUDA kernels written in C++ via C++ wrapper compiled with icl
!> Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
!>--------------------------
module FortranCuda

use MicroMagParameters
use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_null_ptr, c_associated

implicit none

#if USE_CUDA    

    interface

        !> Kaspar K. Nielsen, kasparkn@gmail.com, DTU / Private, 2019
        !> For de-allocating the gpu arrays
        subroutine cu_icl_destroy(  ) bind(C, name="icl_destroy")
            
        end subroutine 
        
        !> Kaspar K. Nielsen, kasparkn@gmail.com, DTU / Private, 2019
        !> for initialization
        subroutine cu_icl_initDemagMatrices( Kxx, Kxy, Kxz, Kyy, Kyz, Kzz, n ) bind(C,name="icl_initDemagMatrices")
            real*4,dimension(n*n), intent(in) :: Kxx,Kxy,Kxz,Kyy,Kyz,Kzz
            integer*4,intent(in) :: n
        end subroutine cu_icl_initDemagMatrices

        !> Kaspar K. Nielsen, kasparkn@gmail.com, DTU / Private, 2019
        !> doing the actual calculations
        subroutine cu_icl_MVMult_GetH( Mx, My, Mz, Hx, Hy, Hz, n, pref ) bind(C,name="icl_MVMult_GetH")
        real*4,dimension(n),intent(in) :: Mx,My,Mz
        real*4,dimension(n),intent(inout) :: Hx,Hy,Hz
        integer*4,intent(in) :: n
        real*4,intent(in) :: pref    

        end subroutine cu_icl_MVMult_GetH
        
        !> Kaspar K. Nielsen, kasparkn@gmail.com, DTU / Private, 2019
        !> initialization of the sparse solver
        subroutine cu_icl_initDemagMatrices_sparse( n, nnz, mat_no, values, colInds, rowInds ) bind(C,name="icl_initDemagMatrices_sparse")
        integer*4,intent(in) :: n, nnz, mat_no
        real*4,dimension(nnz),intent(in) :: values
        integer*4,dimension(nnz),intent(in) :: colInds
        integer*4,dimension(n+1),intent(in) :: rowInds
        
        end subroutine cu_icl_initDemagMatrices_sparse
        
        !> Kaspar K. Nielsen, kasparkn@gmail.com, DTU / Private, 2019
        !> do the sparse matrix multiplication
        subroutine cu_icl_MVMult_GetH_sparse( Mx, My, Mz, Hx, Hy, Hz, n, pref ) bind(C,name="icl_MVMult_GetH_sparse")
        real*4,dimension(n),intent(in) :: Mx,My,Mz
        real*4,dimension(n),intent(inout) :: Hx,Hy,Hz
        integer*4,intent(in) :: n
        real*4,intent(in) :: pref    
        
        end subroutine cu_icl_MVMult_GetH_sparse
        
    end interface

    contains

    !> Kaspar K. Nielsen, kasparkn@gmail.com, DTU / Private, 2019
    !>
    subroutine cudaDestroy()

        call cu_icl_destroy()


    end subroutine cudaDestroy

    !> Initialization of the sparse CUDA solver
    !> Kaspar K. Nielsen, kasparkn@gmail.com, DTU / Private, 2019
    !> Initializes the CUDA sparse matrix-vector product engine
    !> @param[in] K_in is the array of six elements each being a sparse matrix (Kxx, Kxy, Kxz, Kyy, Kyz and Kzz)
    subroutine cudaInit_sparse( K_in)
    type(MagTenseSparse),dimension(6),intent(in) :: K_in       !> Input sparse matrices

    integer*4,dimension(:),allocatable :: rowInds       !> row indices for the CSR format

    integer*4 :: n,nnz,i                                          !> no. of elements, no. of non-zero elements

    !Total no. of rows / cols in the square matrix. Each of the six matrices have the same size
    n = size(K_in(1)%rows_start)


    allocate(rowInds(n+1) )

    !Call the ICL wrapper for the CUDA code for each of the six demag matrices
        
    !Note that the length of rowInds is nnz+1 as CUDA uses the three-array definition of CSR such that
    !the last element of rowInds is rowInds(end) = nnz + rowInds(1)
    do i=1,6
        !No. of non zero elements
        nnz = size(K_in(i)%values)
        
        

        !copy the row indices
        rowInds(1:n) = K_in(i)%rows_start
        rowInds(n+1) = rowInds(1) + nnz

        call cu_icl_initDemagMatrices_sparse( n, nnz, i, K_in(i)%values, K_in(i)%cols, rowInds )
        
        
    enddo
    deallocate(rowInds)
    end subroutine cudaInit_sparse


    !> Kaspar K. Nielsen, kasparkn@gmail.com, DTU / Private, 2019
    !> Called when doing the sparse matrix-vector product on CUDA
    subroutine cudaMatrVecMult_sparse( Mx_in, My_in, Mz_in, Hx, Hy, Hz, pref )
    real*4,dimension(:),intent(in) :: Mx_in,My_in,Mz_in
    real*4,dimension(:),intent(inout) :: Hx,Hy,Hz
    real*4,intent(in) :: pref
    integer :: i,n

    n = size(Mx_in)

    call cu_icl_MVMult_GetH_sparse( Mx_in, My_in, Mz_in, Hx, Hy, Hz, n, pref )

    end subroutine cudaMatrVecMult_sparse

    !> Kaspar K. Nielsen, kasparkn@gmail.com, DTU / Private, 2019
    !> Called by the main program, converts the 64 bit double's to 32 bit float and
    !> calls the icl wrapper to pass the data to the gpu
    !> in double precision
    subroutine cudaInit_d( Kxx_in, Kxy_in, Kxz_in, Kyy_in, Kyz_in, Kzz_in )
    real*8,dimension(:,:) :: Kxx_in,Kxy_in,Kxz_in,Kyy_in,Kyz_in,Kzz_in

    real*4,dimension(:),allocatable :: Kxx, Kxy, Kxz, Kyy, Kyz, Kzz
    integer*4 :: n,i,j,ind

    n = size(Kxx_in(:,1))

    allocate( Kxx(n*n), Kxy(n*n), Kxz(n*n), Kyy(n*n), Kyz(n*n), Kzz(n*n) )


    !convert the 64 bit vars to 32 bit and make the arrays 1d
    do i=1,n
        do j=1,n
            ind = (i-1)*n + j
            Kxx(ind) = sngl(Kxx_in(i,j))
            Kxy(ind) = sngl(Kxy_in(i,j))
            Kxz(ind) = sngl(Kxz_in(i,j))
            Kyy(ind) = sngl(Kyy_in(i,j))
            Kyz(ind) = sngl(Kyz_in(i,j))
            Kzz(ind) = sngl(Kzz_in(i,j))
        enddo
    enddo


    call cu_icl_initDemagMatrices( Kxx, Kxy, Kxz, Kyy, Kyz, Kzz, n )

    !Clean-up
    deallocate( Kxx, Kxy, Kxz, Kyy, Kyz, Kzz )

    end subroutine cudaInit_d

    !> Kaspar K. Nielsen, kasparkn@gmail.com, DTU / Private, 2019
    !> Called by the main program, converts the 64 bit double's to 32 bit float and
    !> calls the icl wrapper to pass the data to the gpu in single precision
    subroutine cudaInit_s( Kxx_in, Kxy_in, Kxz_in, Kyy_in, Kyz_in, Kzz_in )
    real*4,dimension(:,:) :: Kxx_in,Kxy_in,Kxz_in,Kyy_in,Kyz_in,Kzz_in
        
    integer*4 :: n

    n = size(Kxx_in(:,1))

    call cu_icl_initDemagMatrices( Kxx_in, Kxy_in, Kxz_in, Kyy_in, Kyz_in, Kzz_in, n )



    end subroutine cudaInit_s

    !> Kaspar K. Nielsen, kasparkn@gmail.com, DTU / Private, 2020
    !> Called by the main program, does the dense matrix-vector multiplication on the GPU
    subroutine cudaMatrVecMult( Mx_in, My_in, Mz_in, Hx, Hy, Hz, pref )
    real*4,dimension(:),intent(in) :: Mx_in,My_in,Mz_in
    real*4,dimension(:),intent(inout) :: Hx,Hy,Hz
    real*4,intent(in) :: pref
    integer :: i,n

    n = size(Mx_in)

    !call cuda to do the matrix-vector multiplication
    call cu_icl_MVMult_GetH( Mx_in, My_in, Mz_in, Hx, Hy, Hz, n, pref )
    end subroutine cudaMatrVecMult

#endif    
end module FortranCuda