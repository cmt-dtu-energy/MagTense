
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
        

        subroutine cu_icl_dumpDemagMatrices_dense(filename) bind(C, name="icl_dumpDemagMatrices_dense")
            use, intrinsic :: iso_c_binding, only: c_char
            character(kind=c_char, len=1), intent(in) :: filename(*)
        end subroutine cu_icl_dumpDemagMatrices_dense

        !> The exchange field on the GPU (MagTenseCudaExchange.cu)
        function cu_exchange_init( n, nnz, rowptr, cols, vals ) bind(C, name="cu_exchange_init")
            use, intrinsic :: iso_c_binding, only: c_int, c_double
            integer(c_int) :: cu_exchange_init
            integer(c_int),intent(in) :: n, nnz
            integer(c_int),dimension(*),intent(in) :: rowptr, cols
            real(c_double),dimension(*),intent(in) :: vals
        end function cu_exchange_init

        function cu_exchange_apply( Mx, My, Mz, Hx, Hy, Hz ) bind(C, name="cu_exchange_apply")
            use, intrinsic :: iso_c_binding, only: c_int, c_double
            integer(c_int) :: cu_exchange_apply
            real(c_double),dimension(*),intent(in) :: Mx, My, Mz
            real(c_double),dimension(*),intent(inout) :: Hx, Hy, Hz
        end function cu_exchange_apply

        subroutine cu_exchange_destroy() bind(C, name="cu_exchange_destroy")
        end subroutine cu_exchange_destroy

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


    !> Dump the dense demag matrices currently resident on the GPU.
    !> Writes 6 contiguous blocks (Kxx,Kxy,Kxz,Kyy,Kyz,Kzz) as float32 to the file.
    subroutine cudaDumpDemagDense(filename)
    use, intrinsic :: iso_c_binding, only: c_char, c_null_char
    implicit none
    character(len=*), intent(in) :: filename

    integer :: n, i
    character(kind=c_char, len=1) :: c_filename(len_trim(filename) + 1)

    n = len_trim(filename)

    do i = 1, n
        c_filename(i) = filename(i:i)
    end do
    c_filename(n+1) = c_null_char

    call cu_icl_dumpDemagMatrices_dense(c_filename)
    end subroutine cudaDumpDemagDense


    !> Uploads the exchange operator to the GPU, for cudaExchangeApply. The field there is
    !> H_i = scale_i * sum_j A_ij M_j, so the row scaling is folded into the uploaded values.
    !> MKL's CSR is exported as is (one- or zero-based, rows_start/rows_end) and compacted into the
    !> zero-based three-array form of the kernel.
    !> @param[in] A the exchange operator
    !> @param[in] scale the row scaling, size n
    !> @param[out] status 0: on the GPU, 1: disabled by MAGTENSE_EXCHANGE_GPU=0, <0: failed (the CPU is used)
    !> @param[out] nnz the number of non-zeros of the operator
    subroutine cudaExchangeInit( A, scale, status, nnz )
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_double
    type(sparse_matrix_t),intent(in) :: A
    real(DP),dimension(:),intent(in) :: scale
    integer,intent(out) :: status, nnz

    integer :: stat, indexing, nrows, ncols, base, i, j, k
    type(c_ptr) :: rs_c, re_c, ci_c, v_c
    integer,pointer :: rs(:), re(:), ci(:)
    real(DP),pointer :: v(:)
    integer(c_int),dimension(:),allocatable :: rowptr, cols
    real(c_double),dimension(:),allocatable :: vals

    nnz = 0
    stat = mkl_sparse_d_export_csr( A, indexing, nrows, ncols, rs_c, re_c, ci_c, v_c )
    if ( stat .ne. SPARSE_STATUS_SUCCESS .or. nrows .ne. size(scale) ) then
        status = -2
        return
    endif
    base = merge( 1, 0, indexing .eq. SPARSE_INDEX_BASE_ONE )
    call c_f_pointer( rs_c, rs, [nrows] )
    call c_f_pointer( re_c, re, [nrows] )
    nnz = sum( re - rs )
    call c_f_pointer( ci_c, ci, [maxval(re) - base] )
    call c_f_pointer( v_c, v, [maxval(re) - base] )

    allocate( rowptr(nrows+1), cols(nnz), vals(nnz) )
    k = 0
    rowptr(1) = 0
    do i = 1, nrows
        do j = rs(i) - base + 1, re(i) - base
            k = k + 1
            cols(k) = ci(j) - base
            vals(k) = scale(i) * v(j)
        enddo
        rowptr(i+1) = k
    enddo

    status = cu_exchange_init( nrows, nnz, rowptr, cols, vals )
    deallocate( rowptr, cols, vals )
    end subroutine cudaExchangeInit

    !> H = scale * A_exch * M on the GPU for the three components; status 0 on success
    subroutine cudaExchangeApply( Mx, My, Mz, Hx, Hy, Hz, status )
    real(DP),dimension(:),intent(in) :: Mx, My, Mz
    real(DP),dimension(:),intent(inout) :: Hx, Hy, Hz
    integer,intent(out) :: status

    status = cu_exchange_apply( Mx, My, Mz, Hx, Hy, Hz )
    end subroutine cudaExchangeApply

    subroutine cudaExchangeDestroy()
    call cu_exchange_destroy()
    end subroutine cudaExchangeDestroy



#endif    
end module FortranCuda