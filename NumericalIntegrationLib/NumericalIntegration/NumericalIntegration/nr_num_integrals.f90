module nr_num_integrals
    use integrationDataTypes
implicit none

    contains
    
    !::Functions from Numerical Recipes (2nd edition) 
    
    FUNCTION qromb_mod(func,dat,a,b)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : polint,trapzd
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: a,b
    class(dataCollectionBase), intent(inout),target :: dat
	REAL(DP) :: qromb_mod
	procedure (f_int_dat_vec),intent(in), pointer :: func => null ()
	INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
	REAL(DP), PARAMETER :: EPS=1.0e-6_DP
	REAL(DP), DIMENSION(JMAXP) :: h,s
	REAL(DP) :: dqromb
	INTEGER(I4B) :: j
	h(1)=1.0
	do j=1,JMAX
		call trapzd_mod(func,dat,a,b,s(j),j)
		if (j >= K) then
			call polint(h(j-KM:j),s(j-KM:j),0.0_DP,qromb_mod,dqromb)
			if (abs(dqromb) <= EPS*abs(qromb_mod)) RETURN
		end if
		s(j+1)=s(j)
		h(j+1)=0.25_DP*h(j)
	end do
	call nrerror('qromb: too many steps')
    END FUNCTION qromb_mod
    
    SUBROUTINE trapzd_mod(func,dat,a,b,s,n)
	USE nrtype; USE nrutil, ONLY : arth
	IMPLICIT NONE
    class(dataCollectionBase), intent(inout),target :: dat
	REAL(DP), INTENT(IN) :: a,b
	REAL(DP), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	procedure (f_int_dat_vec),intent(in), pointer :: func => null ()
	REAL(DP) :: del,fsum
    real,dimension(:),allocatable :: xarr
	INTEGER(I4B) :: it
	if (n == 1) then
        allocate(xarr(2))
        xarr(1) = a
        xarr(2) = b
		s=0.5_DP*(b-a)*sum(func( xarr, dat ))
        deallocate(xarr)
    else        
		it=2**(n-2)
		del=(b-a)/it
        allocate(xarr(it))
        xarr = arth(a+0.5_DP*del,del,it)
		fsum=sum(func(xarr,dat))
        deallocate(xarr)
		s=0.5_DP*(s+del*fsum)
	end if
	END SUBROUTINE trapzd_mod

    
    
    end module nr_num_integrals