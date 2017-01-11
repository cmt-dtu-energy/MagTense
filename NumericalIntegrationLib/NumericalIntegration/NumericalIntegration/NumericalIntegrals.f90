module NumInt
    use QUADPACK
    implicit none

    !real,parameter :: epsabs=0.001,epsrel=0.01
    real :: y1,y2,x,y,epsabs,epsrel
    !real,external :: func
    
    abstract interface
      function func (x,y)
         real :: func
         real, intent (in) :: x,y
      end function func
   end interface

   procedure (func), pointer :: f_ptr => null ()
    
    contains

!::
!::Numerical double integral
!::Performs the numerical double integral of user supplied function fun
!::on the rectangular interval [xl,xh], [yl,yh]
!::returns the result in res (kind = real)
!::The user supplied function fun should take two arguments, x and y: fun=fun(x,y)
subroutine integral2( fun, xl, xh, yl, yh, eps_abs, eps_rel, res, abserr, neval, ier )
real,external :: fun
real,intent(in) :: xl,xh,yl,yh,eps_abs,eps_rel
real,intent(inout) :: res,abserr
integer,intent(inout) :: neval,ier

epsabs = eps_abs
epsrel = eps_rel

!::Set the y-boundaries
y1 = yl
y2 = yh
!::Set the function to be integrated
f_ptr => fun

call qags_x ( h, xl, xh, epsabs, epsrel, res, abserr, neval, ier )

end subroutine integral2

!::Adopted from NR
function h(xx)
real h,xx
!external f
real ss
real :: res,abserr
integer :: neval,ier
! USES g,qgausy,y1,y2
!Called by qgausx. Calls qgausy.

    x = xx
    call qags_y( f, y1, y2, epsabs, epsrel, ss, abserr, neval, ier )
    h = ss
return
end function h

!::adopted from NR
function f( yy )
real f,yy
! USES func
!Called by qgausy. Calls func.

y = yy

f = f_ptr(x,y)

return
end function f

end module NumInt