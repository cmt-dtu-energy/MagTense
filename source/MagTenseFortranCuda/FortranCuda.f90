
    !>--------------------------
    !> Module for interfacing with CUDA kernels written in C++ via C++ wrapper compiled with icl
    !> Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !>--------------------------
    module FortranCuda

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


    
    end module FortranCuda