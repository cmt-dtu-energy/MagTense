
    !>--------------------------
    !> Module for interfacing with CUDA kernels written in C++ via C++ wrapper compiled with icl
    !> Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !>--------------------------
    module FortranCuda
    
	use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_null_ptr, c_associated


	implicit none

    !Local memory for the M vector arrays to be initialized only once
    real*4,dimension(:),allocatable :: Mx, My, Mz
    
	interface
	
	    !For de-allocating the gpu arrays
	    subroutine cu_icl_destroy(  ) bind(C, name="icl_destroy")
			
        end subroutine 
        !for initialization
        subroutine cu_icl_initDemagMatrices( Kxx, Kxy, Kxz, Kyy, Kyz, Kzz, n ) bind(C,name="icl_initDemagMatrices")
	        real*4,dimension(n*n), intent(in) :: Kxx,Kxy,Kxz,Kyy,Kyz,Kzz
            integer*4,intent(in) :: n
        end subroutine cu_icl_initDemagMatrices
    
        !doing the actual calculations
        subroutine cu_icl_MVMult_GetH( Mx, My, Mz, Hx, Hy, Hz, n, pref ) bind(C,name="icl_MVMult_GetH")
        real*4,dimension(n),intent(in) :: Mx,My,Mz
        real*4,dimension(n),intent(inout) :: Hx,Hy,Hz
        integer*4,intent(in) :: n
        real*4,intent(in) :: pref    
    
        end subroutine cu_icl_MVMult_GetH
        
	end interface

	contains

	subroutine cudaDestroy()
	
		call cu_icl_destroy()
        deallocate(Mx,My,Mz)
	
    end subroutine cudaDestroy

    !> Called by the main program, converts the 64 bit double's to 32 bit float and
    !> calls the icl wrapper to pass the data to the gpu
    subroutine cudaInit( Kxx_in, Kxy_in, Kxz_in, Kyy_in, Kyz_in, Kzz_in )
    real*8,dimension(:,:) :: Kxx_in,Kxy_in,Kxz_in,Kyy_in,Kyz_in,Kzz_in
    
    real*4,dimension(:),allocatable :: Kxx, Kxy, Kxz, Kyy, Kyz, Kzz
    integer*4 :: n,i,j,ind
    
    n = size(Kxx_in(:,1))
    
    allocate( Kxx(n*n), Kxy(n*n), Kxz(n*n), Kyy(n*n), Kyz(n*n), Kzz(n*n) )
    !also initialize the M arrays
    allocate( Mx(n), My(n), Mz(n) )
    
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
    
    end subroutine cudaInit
    
    subroutine cudaMatrVecMult( Mx_in, My_in, Mz_in, Hx, Hy, Hz, pref )
    real*8,dimension(:),intent(in) :: Mx_in,My_in,Mz_in
    real*4,dimension(:),intent(inout) :: Hx,Hy,Hz
    real*4,intent(in) :: pref
    integer :: i,n
    
    n = size(Mx)
    
    !Convert the M arrays to float
    do i=1,n
        Mx(i) = sngl(Mx_in(i))    
        My(i) = sngl(My_in(i))
        Mz(i) = sngl(Mz_in(i))
    enddo
    
    !call cuda to do the matrix-vector multiplication
    call cu_icl_MVMult_GetH( Mx, My, Mz, Hx, Hy, Hz, n, pref )
    end subroutine cudaMatrVecMult
    
    
    end module FortranCuda