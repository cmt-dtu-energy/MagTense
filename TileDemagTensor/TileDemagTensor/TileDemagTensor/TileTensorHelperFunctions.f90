module TileTensorHelperFunctions
    use QuadPack
    implicit none

    real,parameter :: pi=3.141592653589793,mu0=4*pi*1e-7
    contains
    
    
    
        !::Returns the L helper function
        function L( x, theta, z )
        real :: L
        real, intent(in) :: x, theta, z
            L = x**2 * ( cos(theta)**2 - 1 ) - z**2
            return
        end function

        function G_p( r, x, theta, z )
        real :: G_p
        real,intent(in) ::r,x,theta,z
            G_p = ( M(r,x,theta,z) + r ) / sqrt(x**2+z**2)
            return
        end function

        function G_m( r, x, theta, z )
        real :: G_m
        real,intent(in) ::r,x,theta,z
            G_m = ( M(r,x,theta,z) - r ) / sqrt(x**2+z**2)
            return
        end function

        !Returns the M helper function
        function M( r, x, theta, z )
        real :: M
        real,intent(in) :: r, x, theta, z
            M = sqrt( r**2 - 2 * x * r * cos(theta) + x**2 + z**2 )
            return
         end function M

    function C_no_sign( theta )
    real :: C_no_sign
    real,intent(in) :: theta
        C_no_sign = asin( cos( theta/2 ) )
        return
    end function
        !Returns the B helper function
    !::The minus sign is due to the convention that in numerical recipes the incomplete elliptic integral of third kind is defined as 1/(1+n).... where in Matlab it is 1/(1-n)....
    function  B( r, x )
    real :: B
    real,intent(in) :: r,x
        if ( abs( r + x  ) .lt. 1e-20 ) then
            !::The limit of B as x -> -r is zero
            B = 0.
        else
            B = -4 * r * x / ( r + x )**2
        endif
        
        if ( B .le. (-1 + 1e-10) ) then
            B = -1 + 1e-10
        endif
            
        return
    end function

        function K( r, x, z )
        real :: K
        real,intent(in) :: r, x, z
            
            !K = 2 * sqrt( r * x / ( ( r + x )**2 + z**2 ) )
            
            !if ( z .lt. 10*tiny(1.) .OR. abs( r + x ) .lt. 10*tiny(1.) ) then
            if ( z .lt. 1e-20 .AND. abs( r + x ) .lt. 1e-20 .OR. abs(r * x) .lt. 1e-20) then
                K = 0.
            else                
                K = 4 * r * x / ( ( r + x )**2 + z**2 )
            endif
            
            !::The implementation used here needs the k parameter and not the m parameter (as used by Matlab). Remember that m = k^2
            !K = K**2
            return
     end function
     function P(x,y,z,xp,yp,zp)
     real, intent(in) :: x,y,z,xp,yp,zp
     real :: P
           P = sqrt( (x-xp)**2 + (y-yp)**2 + (z-zp)**2 )
     end function
    
     
     
    function atan2_custom( y, x )
    real,intent(in) :: x,y
    real :: atan2_custom
    
    !atan2_custom = atan2( y, x )
    
     if ( x .ge. 0 .AND. y .ge. 0 ) then
        !first quadrant
        atan2_custom = atan( y/x )
    elseif ( x .lt. 0 .AND. y .ge. 0 ) then
        !second quadrant
        atan2_custom = pi - atan( abs(y/x) )
    elseif (x .lt. 0 .AND. y .lt. 0 ) then
        !third
        atan2_custom = pi + atan( abs(y/x) )
    elseif ( x .ge. 0 .AND. y .lt. 0 ) then
        !fourth 
        atan2_custom = 2*pi - atan( abs(y/x) )
    endif
    
    end function
     
end module TileTensorHelperFunctions
    
