module CudaTest

	use, intrinsic :: iso_c_binding, only : c_float,c_int, c_ptr, c_null_ptr, c_associated


	implicit none

	interface
	
	
		subroutine cuTest( K,n ) sbind(C, name="icl_test")
		use,intrinsic :: iso_c_binding, only : c_float
		real(c_float),dimension(n),intent(in) :: K
		integer*4,intent(in) :: n

		end subroutine 
		
	
	end interface

	contains

	subroutine testCuda( val, n )
	real*4,dimension(n),intent(inout) :: val
	integer*4,intent(in) :: n
		write(*,*) val
		call cuTest( val, n );
		
	end subroutine testCuda

end module CudaTest