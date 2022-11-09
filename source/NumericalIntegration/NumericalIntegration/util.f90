module UTIL_CALL

    implicit none
    contains
subroutine getBilinInterp( table, xval,yval,n1, n2, x, y, res )
integer,intent(in) :: n1,n2
real,intent(in) :: x,y
real,intent(out) :: res
real,dimension(n1,n2),intent(in) :: table
real,dimension(n1),intent(in) :: xval
real,dimension(n2),intent(in) :: yval
integer :: ind_x,ind_y
real :: x_lin,y_lin

    res = 0.

    call locate( xval,n1,x,ind_x )
    
    !::Used if the interpolation comes outside the range
    if ( ind_x .eq. 0 ) return
    
    x_lin = ( x - xval(ind_x) ) / ( xval(ind_x+1) - xval(ind_x) )
    call locate( yval,n2,y,ind_y )
    
    !::Used if the interpolation comes outside the range
    if ( ind_y .eq. 0 ) return
    
    y_lin = (y-yval(ind_y)) / ( yval(ind_y+1)-yval(ind_y))
    
    
    res = (1-y_lin)*( (1.-x_lin) * table(ind_x,ind_y) + x_lin*table(ind_x+1,ind_y) ) &
        +    y_lin *( (1.-x_lin) * table(ind_x,ind_y+1) + x_lin*table(ind_x+1,ind_y+1) ) 


end subroutine getBilinInterp


subroutine interp1_MagTense( x, y, xval, n, yval)
real,dimension(n),intent(in) :: x,y
real,intent(in) :: xval
integer,intent(in) :: n
real,intent(inout) :: yval

integer :: ind_x
real :: x_lin

    yval = 0
    call locate( x, n, xval, ind_x )

    if ( ind_x .eq. 0 ) return
    
    x_lin = ( xval - x(ind_x) ) / ( x(ind_x+1) - x(ind_x) )
    
    yval = ( 1 - x_lin ) * y(ind_x) + x_lin * y(ind_x+1)
    
end subroutine interp1_MagTense


subroutine locate(xx,n,x,j)
    !:: Locate the element j nearest the input-value x in a array named xx of a length n.
    integer :: j,n
    real :: x,xx(n)
    integer :: jl,jm,ju

    jl=0
    ju=n+1
10    if(ju-jl.gt.1)then
      jm=(ju+jl)/2
      if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
        jl=jm
      else
        ju=jm
      endif
    goto 10
    endif
    if(x.eq.xx(1))then
      j=1
    else if(x.eq.xx(n))then
      j=n-1
    else
      j=jl
    endif
    return
end subroutine locate


end module UTIL_CALL