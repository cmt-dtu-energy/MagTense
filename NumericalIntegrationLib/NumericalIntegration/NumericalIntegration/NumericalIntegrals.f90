module NumInt
    use QUADPACK
    use integrationDataTypes
    implicit none
    
    contains

    !::surf defines the min and max vals of the x, y and z integration limits
    !::dat_arr is an array of 18 elements each being a pointer to a data structure
    !::used by the integration routine
    !::retVector is a 3-element logical vector that determines if the component of the
    !::vector should be computed
    !::handleError is a pointer to a specific subroutine that deals with any errors from the integration
    !::vOut is a 3-element vector that returns the solution of the three vector components
    !::ier is a 2-element array with error-codes from the x- and y-integrations
    !::neval is a 2-element arrant with the maximum number of function evaluations from the x- and y-integrations
    subroutine surface_integral_carth( surf, dat_arr, retVector, handleError, vOut, ier, neval )
    type( surf_carth ), intent(in) :: surf
    type( dat_ptr ), intent(inout), dimension(18) :: dat_arr
    integer,intent(in),dimension(3) :: retVector
    procedure (error_handler), intent(in), pointer :: handleError
    
    real,dimension(3),intent(inout) :: vOut
    
    integer,intent(inout),dimension(2) :: ier,neval
    real,dimension(6,3) :: n_vec
    real,dimension(6,2) :: xp,yp !::defines the plane of integration
    real,dimension(6) :: zp,vecSign !::defines the position of the plane on the third axis
    integer :: i,j,ind
    real :: res
    
    !::Reset input array
    vOut(:) = 0.
    
    neval = 0
      
    !::Setup the six planes to be integrated over
    !::yz- plane
    xp(1,:) = surf%y
    yp(1,:) = surf%z
    zp(1) = surf%x(1)
    vecSign(1) = -1
    n_vec(1,:) = (/ -1, 0, 0 /)
      
    !::yz+ plane
    xp(2,:) = surf%y
    yp(2,:) = surf%z
    zp(2) = surf%x(2)
    vecSign(2) = 1
    n_vec(2,:) = (/ 1, 0, 0 /)
      
    !::xz- plane
    xp(3,:) = surf%x
    yp(3,:) = surf%z
    zp(3) = surf%y(1)
    vecSign(3) = -1
    n_vec(3,:) = (/ 0, -1, 0 /)
      
    !::xz+ plane
    xp(4,:) = surf%x      
    yp(4,:) = surf%z
    zp(4) = surf%y(2)
    vecSign(4) = 1
    n_vec(4,:) = (/ 0, 1, 0 /)
      
    !::xy- plane
    xp(5,:) = surf%x
    yp(5,:) = surf%y
    zp(5) = surf%z(1)
    vecSign(5) = -1
    n_vec(5,:) = (/ 0, 0, -1 /)
      
    !::xy+ plane
    xp(6,:) = surf%x
    yp(6,:) = surf%y
    zp(6) = surf%z(2)
    vecSign(6) = 1
    n_vec(6,:) = (/ 0, 0, 1 /)
    
    !::j loops over the three vector components     
    do j=1,3          
        if ( retVector(j) .eq. 0 ) then       
            !$OMP PARALLEL PRIVATE(i,ind,res)
            do i=1,6
                ind = 6*(j-1) + i
                
                dat_arr(ind)%dat%x1 = xp(i,1)
                dat_arr(ind)%dat%x2 = xp(i,2)
                dat_arr(ind)%dat%y1 = yp(i,1)
                dat_arr(ind)%dat%y2 = yp(i,2)
                
                dat_arr(ind)%dat%z0 = zp(i)
                
                dat_arr(ind)%dat%n_vec = n_vec(i,:)            
                
                call integral2( dat_arr(ind)%dat, res )
            
                call handleError(dat_arr(ind)%dat, dat_arr(ind)%dat%abserr_tot) 
                !$OMP ATOMIC
                vOut(j) = vOut(j) + vecSign(i) * res           
            enddo      
            !$OMP END PARALLEL
        endif        
    enddo            
    neval = 0
    ier = 0
    do i=1,18
        neval(1) = neval(1) + dat_arr(i)%dat%neval_x
        neval(2) = neval(2) + dat_arr(i)%dat%neval_y
        
        if ( dat_arr(i)%dat%ier_x .ne. 0 ) then
            ier(1) = dat_arr(i)%dat%ier_x
        endif
        
        if ( dat_arr(i)%dat%ier_y .ne. 0 ) then
            ier(2) = dat_arr(i)%dat%ier_y
        endif
      
    enddo
    
    end subroutine surface_integral_carth
    
!::
!::Numerical double integral
!::Performs the numerical double integral of user supplied function fun
!::on the rectangular interval [xl,xh], [yl,yh]
!::returns the result in res (kind = real)
!::The user supplied function fun should take two arguments, x and y: fun=fun(x,y)
subroutine integral2( dat, res )
real,intent(inout) :: res
class(dataCollectionBase),intent(inout), target :: dat
procedure (f_int_dat), pointer :: f_ptr => null ()

f_ptr => h

call qags_x ( f_ptr, dat, dat%x1, dat%x2, dat%epsabs, dat%epsrel, res, dat%abserr_x, dat%neval_x, dat%ier_x )

end subroutine integral2

!::Adopted from NR
function h(xx, dat)
real :: h
real,intent(in) :: xx
!external f
real :: ss, res
class(dataCollectionBase), target :: dat
procedure (f_int_dat), pointer :: f_ptr => null ()

    f_ptr => f
    dat%x = xx
    call qags_y( f_ptr, dat, dat%y1, dat%y2, dat%epsabs, dat%epsrel, ss, dat%abserr_y, dat%neval_y, dat%ier_y )
    h = ss
    
return
end function h

!::adopted from NR
function f( yy, dat )
real :: f
real,intent(in) :: yy
class(dataCollectionBase), target :: dat
!Called by qgausy. Calls func.

dat%y = yy

f = dat%f_ptr(dat)

return
end function f

end module NumInt