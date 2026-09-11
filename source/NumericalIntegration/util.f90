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
    !:: Given a set coordintaes x, y and a given x-coordinate, xval, estimate the corresponding y-value.
    !:: Uses a linear interpolation from (x_n, y_n) to (x_{n+1}, y_{n+1}) where x_n and x_{n+1}
    !:: are the x-values closest to xval
real,dimension(n),intent(in) :: x,y ! Coordinates
real,intent(in) :: xval ! x-value to estimate y at
integer,intent(in) :: n ! Number of xy coordinates
real,intent(inout) :: yval  ! Estimated y value

integer :: ind_x
real :: x_lin

    yval = 0
    call locate( x, n, xval, ind_x )

    if ( ind_x .eq. 0 ) return
    
    x_lin = ( xval - x(ind_x) ) / ( x(ind_x+1) - x(ind_x) )

    ! Note how yval = y(ind_x) when xval = x(ind_x) and yval = y(ind_x + 1) when xval = x(ind_x+1)
    ! Linear interpolation applies between these two extremes
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


!>-----------------------------------------
!> @author Rasmus Bj�rk, rabj@dtu.dk, DTU, 2025
!> @brief
!> unique_sort sorts a 1D array and returns only the unique values
!> 
!> @param[in] val - the values that should be sorted and made unique
!> @param[out] unique_val - the unique values from val
!> @param[out] n_unique - the number of unique values in val, i.e. the length of unique_val
!---------------------------------------------------------------------------
subroutine unique_sort(val, unique_val, n_unique)
    implicit none
    integer :: i, min_val, max_val
    integer, dimension(:), intent(in) :: val 
    integer, dimension(:), intent(out) :: unique_val
    integer, intent(out) :: n_unique
    integer, dimension(:), allocatable :: unique
    logical, allocatable :: mask1D(:)
      
    allocate(unique(size(val)))
    min_val = minval(val)-1
    max_val = maxval(val)
    
    allocate(mask1D(size(val)))
    
    i = 0
    do while (min_val < max_val)
        i = i+1
        mask1D = (val > min_val)
        min_val = minval(val, mask1D)
        unique(i) = min_val
    enddo
    unique_val(1:i) = unique(1:i)
    
    n_unique = i
    
    deallocate(unique,mask1D)
    
end subroutine unique_sort


end module UTIL_CALL