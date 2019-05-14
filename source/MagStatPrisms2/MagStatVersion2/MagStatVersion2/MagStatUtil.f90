module MagStatUtil

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
    
    
    
    
    
end module
    