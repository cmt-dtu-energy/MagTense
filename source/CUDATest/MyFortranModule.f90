module CudaTest

	use, intrinsic :: iso_c_binding, only : c_float,c_int, c_ptr, c_null_ptr, c_associated


	implicit none

	interface
	
	
	subroutine vecAddTest( val, n ) sbind(C, name="icl_vecAddTest")
	
	real*4,dimension(n) :: val
	integer*4 :: n
	
	

	end subroutine 
	
	subroutine init( val, n ) sbind(C,name="icl_initBuf")
	real*4,dimension(n) :: val
	integer*4 :: n
	
	end subroutine init	
	
	subroutine destroy() sbind(C,name="icl_destroy")
	
	end subroutine destroy
	
	end interface

	contains

	subroutine testCuda( val, n )
	real*4,dimension(n),intent(inout) :: val
	integer*4,intent(in) :: n
	
		write(*,*) 'before', val
		call init( val, n )
		val(:) = 0.
		call vecAddTest(val,n)
		write(*,*) 'after', val
		
		call destroy()
		
	end subroutine testCuda

end module CudaTest