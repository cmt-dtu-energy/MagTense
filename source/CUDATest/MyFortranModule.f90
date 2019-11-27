module CudaTest

	use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_null_ptr, c_associated


	implicit none

	interface
	
	
	subroutine vecAddTest(  ) bind(C, name="vecAddTest")
	

	implicit none
	

	end subroutine 
	
	end interface

	contains

	subroutine testCuda()
	
		call vecAddTest()
	
	end subroutine testCuda

end module CudaTest