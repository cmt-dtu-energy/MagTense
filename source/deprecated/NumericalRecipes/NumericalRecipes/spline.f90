	SUBROUTINE spline(x,y,yp1,ypn,y2,n)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : tridag
	IMPLICIT NONE
	INTEGER(I4B) :: n
    REAL(DP), DIMENSION(n), INTENT(IN) :: x,y
	REAL(DP), INTENT(IN) :: yp1,ypn
	REAL(DP), DIMENSION(n), INTENT(OUT) :: y2
	
	REAL(DP), DIMENSION(n) :: a,b,c,r
	!n=assert_eq(size(x),size(y),size(y2),'spline')
	c(1:n-1)=x(2:n)-x(1:n-1)
	r(1:n-1)=6.0_dp*((y(2:n)-y(1:n-1))/c(1:n-1))
	r(2:n-1)=r(2:n-1)-r(1:n-2)
	a(2:n-1)=c(1:n-2)
	b(2:n-1)=2.0_dp*(c(2:n-1)+a(2:n-1))
	b(1)=1.0
	b(n)=1.0
	if (yp1 > 0.99e30_dp) then
		r(1)=0.0
		c(1)=0.0
	else
		r(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
		c(1)=0.5
	end if
	if (ypn > 0.99e30_dp) then
		r(n)=0.0
		a(n)=0.0
	else
		r(n)=(-3.0_dp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
		a(n)=0.5
	end if
	call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
	END SUBROUTINE spline
