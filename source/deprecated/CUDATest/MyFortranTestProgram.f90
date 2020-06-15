program p

	use CudaTest

	use, intrinsic :: iso_fortran_env, only : compiler_version
	real*4,dimension(:),allocatable :: val
	integer :: i,n
	print *, "Compiler Version: ", compiler_version()
	
	n = 10
	
	allocate(val(n))
	
	do i=1,n
		val(i) = i
	enddo
	
	call testCuda(val,n)
	
	
	deallocate(val)
	stop
end program p