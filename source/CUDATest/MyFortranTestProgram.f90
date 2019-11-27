program p

	use CudaTest

	use, intrinsic :: iso_fortran_env, only : compiler_version


	print *, "Compiler Version: ", compiler_version()

	call testCuda()

	stop
end program p