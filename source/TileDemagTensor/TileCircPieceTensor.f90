    module TileCircPieceTensor
    use QuadPack
    use SpecialFunctions
    use TileTensorHelperFunctions
    
    implicit none
    
    !::Relative distance from the z-axis, xrot/R, below which the six arc integrals are evaluated from
    !::their expansion to first order in xrot (int_axis_dtheta_dz) instead of the closed forms.
    !::The closed forms contain 1/xrot**2 and lose about two digits per decade below xrot ~ 1e-2 R
    !::(relative error ~1e-9 at 1e-3 R, ~5e-8 at 1e-4 R, ~6e-6 at 1e-5 R), while the first-order
    !::expansion has a relative error of about 2 (xrot/R)**2. The two cross at xrot ~ 2e-4 R, where both are about 1e-7.
    real, parameter :: xrot_axis_tol = 2e-4
    
    !::Relative distance from the z-axis, xrot/R, below which the arc parts of the end-surface integrals
    !::capM and capN are evaluated from their expansion to first order in xrot (circPiece_arcMN). The
    !::closed forms contain 1/xrot and lose about one digit per decade below xrot ~ 1e-1 R, while the
    !::expansion has a relative error of about (xrot/R)**2; both are about 1e-10 at the switch.
    real, parameter :: xrot_cap_tol = 1e-5

    !::Relative distance from the z-axis below which the solid angle of the circular sector
    !::(circPiece_omega_sector) is replaced by its on-axis value, with a relative error of about xrot/R.
    real, parameter :: xrot_sector_tol = 1e-10
    
    contains
    
    !Note on the definitions of elliptic integrals in Matlab and in Maple
        !
        !Elliptic integrals of the third kind (Pi)
        !In Maple the arguments are z, nu, k
        !In Matlab the arguments are n, phi, m
        !They correspond like this:
        !z corresponds to acos( phi )
        !nu corresponds to n
        !k**2 corresponds to m

    
        subroutine int_ddx_cos_dtheta_dz( dat, val )
        real,intent(inout) :: val
        class(dataCollectionBase), intent(inout), target :: dat
        real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
        real :: dtheta
        
        val = 0.
           
        call getParameters_rot_trick( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
        
        !::On the z-axis (and close to it) use the elementary small-xrot expansion of the integral
        if ( xrot .lt. xrot_axis_tol * R ) then
            val = int_axis_dtheta_dz( 1, R, xrot, theta1, theta2, z1, z2 )
            return
        endif
            
        if ( theta1 .lt. 0. .AND. theta2 .lt. 0. ) then
            dtheta = theta2-theta1
            theta1 = theta1-2*pi*floor(theta1/(2*pi))
            theta2 = theta1 + dtheta
        else if ( theta2 .gt. 2*pi ) then
            dtheta = theta2 - theta1
            theta2 = theta2-2*pi*floor(theta2/(2*pi))
            theta1 = theta2 - dtheta
        endif
                                                                        
        if ( theta1 .lt. 0. .and. theta2 .gt. 0. ) then
            !first integrate from theta1 to zero
            val = int_ddx_cos_dtheta_dz_fct( 0., z2, R, xrot ) - int_ddx_cos_dtheta_dz_fct( theta1, z2, R, xrot ) - ( int_ddx_cos_dtheta_dz_fct( 0., z1, R, xrot ) - int_ddx_cos_dtheta_dz_fct( theta1, z1, R, xrot ) )
            !then integrate from zero to theta2
            val = val + ( 2 * int_ddx_cos_dtheta_dz_fct( 0., z2, R, xrot ) - int_ddx_cos_dtheta_dz_fct( theta2, z2, R, xrot )) - int_ddx_cos_dtheta_dz_fct( 0., z2, R, xrot ) - ( ( 2 * int_ddx_cos_dtheta_dz_fct( 0., z1, R, xrot ) - int_ddx_cos_dtheta_dz_fct( theta2, z1, R, xrot ) ) - int_ddx_cos_dtheta_dz_fct( 0., z1, R, xrot ) )  
            val = -1*val
        else
            val = (int_ddx_cos_dtheta_dz_fct(theta2,z2,R,xrot) - int_ddx_cos_dtheta_dz_fct(theta1, z2,R,xrot) - ( int_ddx_cos_dtheta_dz_fct(theta2,z1,R,xrot) - int_ddx_cos_dtheta_dz_fct(theta1,z1,R,xrot) ))
        endif                    
        !val = int_ddx_cos_dtheta_dz_fct(theta2,z2,R,xrot)
    end subroutine int_ddx_cos_dtheta_dz
    
    function int_ddx_cos_dtheta_dz_fct( thetap, zp, R, xrot )
    real,intent(in) :: thetap,zp,R,xrot
    real :: int_ddx_cos_dtheta_dz_fct
    real :: elf,elE,elPi
            
        elf = ellf( C_no_sign(thetap), K(R,xrot,zp) )
        elE = elle( C_no_sign(thetap), K(R,xrot,zp) )
        elPi = ellpi( C_no_sign(thetap), B(R,xrot), K(R,xrot,zp) )
        int_ddx_cos_dtheta_dz_fct = zp/(2*R*xrot**2*(R+xrot)*sqrt((R+xrot)**2+zp**2)) * (-(xrot-R)*(R**2+xrot**2) * elPi + (R+xrot)*(((R+xrot)**2+zp**2)*elE - (2*R**2+zp**2)*elF))              
    
    end function int_ddx_cos_dtheta_dz_fct
    
    
        
    subroutine int_ddx_sin_dtheta_dz( dat, val )
    real,intent(inout) :: val
    class(dataCollectionBase), intent(inout), target :: dat
    real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
    real :: dtheta
        
    val = 0.
           
    call getParameters_rot_trick( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
    
    !::On the z-axis (and close to it) use the elementary small-xrot expansion of the integral
    if ( xrot .lt. xrot_axis_tol * R ) then
        val = int_axis_dtheta_dz( 2, R, xrot, theta1, theta2, z1, z2 )
        return
    endif
                                             
            
    val = int_ddx_sin_dtheta_dz_fct(theta2,z2,R,xrot) - int_ddx_sin_dtheta_dz_fct(theta1,z2,R,xrot) - ( int_ddx_sin_dtheta_dz_fct(theta2,z1,R,xrot) - int_ddx_sin_dtheta_dz_fct(theta1,z1,R,xrot) )

    end subroutine
        
    function int_ddx_sin_dtheta_dz_fct(thetap, zp, R, xrot )
    real,intent(in) :: thetap,zp,R,xrot
    real :: int_ddx_sin_dtheta_dz_fct
    
        int_ddx_sin_dtheta_dz_fct = -1/(4*R*xrot**2) * ( (R**2-xrot**2)*log(M(R,xrot,thetap,zp)-zp) + (xrot**2-R**2)*log(M(R,xrot,thetap,zp)+zp) - 2*zp*M(R,xrot,thetap,zp) )
    end function
    
    subroutine int_ddy_cos_dtheta_dz( dat, val )
    real,intent(inout) :: val
    class(dataCollectionBase), intent(inout), target :: dat
    real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
    real :: dtheta
        
    val = 0.
           
    call getParameters_rot_trick( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
    
    !::On the z-axis (and close to it) use the elementary small-xrot expansion of the integral
    if ( xrot .lt. xrot_axis_tol * R ) then
        val = int_axis_dtheta_dz( 3, R, xrot, theta1, theta2, z1, z2 )
        return
    endif
                                                                             
    val = int_ddy_cos_dtheta_dz_fct(theta2,z2,R,xrot) - int_ddy_cos_dtheta_dz_fct(theta1,z2,R,xrot) - ( int_ddy_cos_dtheta_dz_fct(theta2,z1,R,xrot) - int_ddy_cos_dtheta_dz_fct(theta1,z1,R,xrot) )

    end subroutine
    
    function int_ddy_cos_dtheta_dz_fct(thetap, zp, R, xrot )
    real,intent(in) :: thetap,zp,R,xrot
    real :: int_ddy_cos_dtheta_dz_fct
        
            int_ddy_cos_dtheta_dz_fct = -1/(2*R*xrot**2) * ( (R**2+xrot**2)*log(M(R,xrot,thetap,zp)+zp) + zp*M(R,xrot,thetap,zp) )
        
    end function
    
    
    subroutine int_ddy_sin_dtheta_dz( dat, val )
    real,intent(inout) :: val
    class(dataCollectionBase), intent(inout), target :: dat
    real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
    real :: dtheta
        
    val = 0.
           
    call getParameters_rot_trick( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
    
    !::On the z-axis (and close to it) use the elementary small-xrot expansion of the integral
    if ( xrot .lt. xrot_axis_tol * R ) then
        val = int_axis_dtheta_dz( 4, R, xrot, theta1, theta2, z1, z2 )
        return
    endif
        
    if ( theta1 .lt. 0 .AND. theta2 .lt. 0. ) then
        dtheta = theta2-theta1
        theta1 = theta1-2*pi*floor(theta1/(2*pi))
        theta2 = theta1 + dtheta
    else if ( theta2 .gt. 2*pi ) then
        dtheta = theta2 - theta1;
        theta2 = theta2-2*pi*floor(theta2/(2*pi))
        theta1 = theta2 - dtheta
    endif
                                                              
    if ( theta1 .le. 0. .AND. theta2 .gt. 0. ) then
        !first integrate from theta1 to zero
        val = int_ddy_sin_dtheta_dz_fct( 0., z2, R, xrot ) - int_ddy_sin_dtheta_dz_fct( theta1, z2, R, xrot ) - ( int_ddy_sin_dtheta_dz_fct( 0., z1, R, xrot ) - int_ddy_sin_dtheta_dz_fct( theta1, z1, R, xrot ) )
        !then integrate from zero to theta2
        val = val + ( 2 * int_ddy_sin_dtheta_dz_fct( 0., z2, R, xrot ) - int_ddy_sin_dtheta_dz_fct( theta2, z2, R, xrot )) - int_ddy_sin_dtheta_dz_fct( 0., z2, R, xrot ) - ( ( 2 * int_ddy_sin_dtheta_dz_fct( 0., z1, R, xrot ) - int_ddy_sin_dtheta_dz_fct( theta2, z1, R, xrot ) ) - int_ddy_sin_dtheta_dz_fct( 0., z1, R, xrot ) )
        val = -1*val
    else
        val = (int_ddy_sin_dtheta_dz_fct(theta2,z2, R, xrot) - int_ddy_sin_dtheta_dz_fct(theta1, z2, R, xrot) - ( int_ddy_sin_dtheta_dz_fct(theta2,z1, R, xrot) - int_ddy_sin_dtheta_dz_fct(theta1,z1, R, xrot) ))
    endif
            
            
     end subroutine
     
     function int_ddy_sin_dtheta_dz_fct(thetap, zp, R, xrot )
        real,intent(in) :: thetap,zp,R,xrot
        real :: int_ddy_sin_dtheta_dz_fct
        real :: elf,elE,elPi
        
        elf = ellf( C_no_sign(thetap), K(R,xrot,zp) )
        elE = elle( C_no_sign(thetap), K(R,xrot,zp) )
        elPi = ellpi( C_no_sign(thetap), B(R,xrot), K(R,xrot,zp) )
        
            int_ddy_sin_dtheta_dz_fct = zp  / ( 2*R*xrot**2*sqrt((R+xrot)**2+zp**2) ) * ( (xrot-R)**2*elPi + ( (R+xrot)**2 + zp**2 ) * elE - 2* ( xrot**2 + R**2 + zp**2/2) * elF) 
        
     end function
     
     subroutine int_ddz_cos_dtheta_dz( dat, val )
     real,intent(inout) :: val
     class(dataCollectionBase), intent(inout), target :: dat
     real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
     real :: dtheta
        
         val = 0.
           
         call getParameters_rot_trick( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
         
         !::On the z-axis (and close to it) use the elementary small-xrot expansion of the integral
         if ( xrot .lt. xrot_axis_tol * R ) then
             val = int_axis_dtheta_dz( 5, R, xrot, theta1, theta2, z1, z2 )
             return
         endif
            
        if ( theta1 .lt. 0 .AND. theta2 .lt. 0 ) then
            dtheta = theta2-theta1
            theta1 = theta1-2*pi*floor(theta1/(2*pi))
            theta2 = theta1 + dtheta
        elseif ( theta2 .gt. 2*pi ) then
            dtheta = theta2 - theta1
            theta2 = theta2-2*pi*floor(theta2/(2*pi))
            theta1 = theta2 - dtheta
        endif
                                                             
        if ( theta1 .lt. 0 .AND. theta2 .gt. 0. ) then
            !first integrate from theta1 to zero
            val = int_ddz_cos_dtheta_dz_fct( 0., z2, R, xrot ) - int_ddz_cos_dtheta_dz_fct( theta1, z2, R, xrot ) - ( int_ddz_cos_dtheta_dz_fct( 0., z1, R, xrot ) - int_ddz_cos_dtheta_dz_fct( theta1, z1, R, xrot ) )
            !then integrate from zero to theta2
            val = val + ( 2 * int_ddz_cos_dtheta_dz_fct( 0., z2, R, xrot ) - int_ddz_cos_dtheta_dz_fct( theta2, z2, R, xrot )) - int_ddz_cos_dtheta_dz_fct( 0., z2, R, xrot ) - ( ( 2 * int_ddz_cos_dtheta_dz_fct( 0., z1, R, xrot ) - int_ddz_cos_dtheta_dz_fct( theta2, z1, R, xrot ) ) - int_ddz_cos_dtheta_dz_fct( 0., z1, R, xrot ) )
            val = -1*val
        else
            val = (int_ddz_cos_dtheta_dz_fct(theta2,z2, R, xrot) - int_ddz_cos_dtheta_dz_fct(theta1, z2, R, xrot) - ( int_ddz_cos_dtheta_dz_fct(theta2,z1, R, xrot) - int_ddz_cos_dtheta_dz_fct(theta1,z1, R, xrot) ))
        endif
                            
    end subroutine
     
    function int_ddz_cos_dtheta_dz_fct(thetap, zp, R, xrot )
        real,intent(in) :: thetap,zp,R,xrot
        real :: int_ddz_cos_dtheta_dz_fct
        real :: elf,elE
        
        elf = ellf( C_no_sign(thetap), K(R,xrot,zp) )
        elE = elle( C_no_sign(thetap), K(R,xrot,zp) )
        
        int_ddz_cos_dtheta_dz_fct = -1 / (R*xrot*sqrt((R+xrot)**2+zp**2)) * ( ((R+xrot)**2+zp**2) * elE - (R**2 + xrot**2 + zp**2) * elF)
                
    end function
    
    subroutine int_ddz_sin_dtheta_dz( dat, val )
     real,intent(inout) :: val
     class(dataCollectionBase), intent(inout), target :: dat
     real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
     real :: dtheta
        
        val = 0.
           
        call getParameters_rot_trick( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )    
        
        !::On the z-axis (and close to it) use the elementary small-xrot expansion of the integral
        if ( xrot .lt. xrot_axis_tol * R ) then
            val = int_axis_dtheta_dz( 6, R, xrot, theta1, theta2, z1, z2 )
            return
        endif
            
        val = int_ddz_sin_dtheta_dz_fct(theta2,z2,R,xrot) - int_ddz_sin_dtheta_dz_fct(theta1,z2,R,xrot) - ( int_ddz_sin_dtheta_dz_fct(theta2,z1,R,xrot) - int_ddz_sin_dtheta_dz_fct(theta1,z1,R,xrot) )
            
    end subroutine
    
    function int_ddz_sin_dtheta_dz_fct(thetap, zp, R, xrot )
        real,intent(in) :: thetap,zp,R,xrot
        real :: int_ddz_sin_dtheta_dz_fct
        
        int_ddz_sin_dtheta_dz_fct = - M(R,xrot,thetap,zp) / ( R * xrot )
    
    end function
    
    
    subroutine int_ddx_dx_dz( dat, val )
     real,intent(inout) :: val
     class(dataCollectionBase), intent(inout), target :: dat
     real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
     real :: x1,x2,x3,y1,y2,y3
     real :: dtheta,theta0
            
     call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
     call getCorners( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
                                   
    !multiply with sign(cos(theta0)) in order to take the order of
    !integration into account properly (when in 2nd and 3rd
    !quadrants x1 < x3 while in 1st and 4th x3 < x1)
    val = sign(1.,cos(theta0)) * ( int_ddx_dx_dz_fct( x1, z2, y3, x, y, z, theta0 ) - int_ddx_dx_dz_fct( x3, z2, y3, x, y, z, theta0 ) - ( int_ddx_dx_dz_fct(x1,z1, y3, x, y, z, theta0) - int_ddx_dx_dz_fct(x3,z1, y3, x, y, z, theta0) ) )
    
                        
    end subroutine
    !::for the inverted circ piece, i.e. pointing radially inwards
    subroutine int_ddx_dx_dz_inv( dat, val )
     real,intent(inout) :: val
     class(dataCollectionBase), intent(inout), target :: dat
     real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
     real :: x1,x2,x3,y1,y2,y3
     real :: dtheta,theta0
            
     call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
     call getCorners_inv( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
                                   
    !multiply with sign(cos(theta0)) in order to take the order of
    !integration into account properly (when in 2nd and 3rd
    !quadrants x1 < x3 while in 1st and 4th x3 < x1)
    !change sign as the normal vector is now pointing along the positive y-axis
    val = -sign(1.,cos(theta0)) * ( int_ddx_dx_dz_fct( x3, z2, y3, x, y, z, theta0 ) - int_ddx_dx_dz_fct( x2, z2, y3, x, y, z, theta0 ) - ( int_ddx_dx_dz_fct(x3,z1, y3, x, y, z, theta0) - int_ddx_dx_dz_fct(x2,z1, y3, x, y, z, theta0) ) )
                        
    end subroutine
    
    
    function int_ddx_dx_dz_fct( xp, zp, y3, x, y, z, theta0 )
        real,intent(in) :: xp,zp, y3, x, y, z, theta0
        real :: int_ddx_dx_dz_fct,arg
        
        arg = zp-z + P(x,y,z,xp,y3,zp)
        
        if ( arg .lt. 10*tiny(1.) ) then
            arg = 10*tiny(1.)
        endif
             
        int_ddx_dx_dz_fct = -sign(1.,sin(theta0))/ (4*pi) * log( arg )
        
                
    end function
    
    
     subroutine int_ddy_dx_dz(dat, val )
     real,intent(inout) :: val
     class(dataCollectionBase), intent(inout), target :: dat
     real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
     real :: x1,x2,x3,y1,y2,y3
     real :: dtheta,theta0
            
     call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
     call getCorners( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
     
     val = sign(1.,cos(theta0)) * sign(1.,sin(theta0)) * ( int_ddy_dx_dz_fct( x1, z2, y3, x, y, z, theta0 ) - int_ddy_dx_dz_fct( x3, z2, y3, x, y, z, theta0 ) - ( int_ddy_dx_dz_fct(x1,z1, y3, x, y, z, theta0 ) - int_ddy_dx_dz_fct(x3,z1, y3, x, y, z, theta0 )) )
     !val = sign(1.,sin(theta0)) * ( int_ddy_dx_dz_fct( x1, z2, y3, x, y, z, theta0 ) - int_ddy_dx_dz_fct( x3, z2, y3, x, y, z, theta0 ) - ( int_ddy_dx_dz_fct(x1,z1, y3, x, y, z, theta0 ) - int_ddy_dx_dz_fct(x3,z1, y3, x, y, z, theta0 )) )

     end subroutine  
     
     subroutine int_ddy_dx_dz_inv(dat, val )
     real,intent(inout) :: val
     class(dataCollectionBase), intent(inout), target :: dat
     real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
     real :: x1,x2,x3,y1,y2,y3
     real :: dtheta,theta0
            
     call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
     call getCorners_inv( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
     !Change the sign as the normal vector is pointing along the positive y-axis
     val = -sign(1.,cos(theta0)) * sign(1.,sin(theta0)) * ( int_ddy_dx_dz_fct( x3, z2, y3, x, y, z, theta0 ) - int_ddy_dx_dz_fct( x2, z2, y3, x, y, z, theta0 ) - ( int_ddy_dx_dz_fct(x3,z1, y3, x, y, z, theta0 ) - int_ddy_dx_dz_fct(x2,z1, y3, x, y, z, theta0 )) )

     end subroutine  
     
    function int_ddy_dx_dz_fct( xp, zp, y3, x, y, z, theta0 )
        real,intent(in) :: xp,zp, y3, x, y, z, theta0
        real :: int_ddy_dx_dz_fct
        !::Make sure the limit when y-y3 -> 0 is covered. 
        !::In this limit the derivative of the nominator will go to zero (per l'Hopital's rule) 
        !::and the denominator is finite (the derivative of P wrt y is non-zero when x-xp != 0 or z-zp != 0, which is then a requirement
        if ( abs(y-y3) .lt. 10*tiny(1.) ) then
            int_ddy_dx_dz_fct = -1/(4*pi) * atan( huge(1.) )
        else    
            int_ddy_dx_dz_fct = -1/(4*pi) * atan( sign(1.,cos(theta0))*(x-xp)*(z-zp)/ ((y-y3)* P(x,y,z,xp,y3,zp)) )
        endif
        
    end function
    
    
    subroutine int_ddz_dx_dz(dat, val )
     real,intent(inout) :: val
     class(dataCollectionBase), intent(inout), target :: dat
     real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
     real :: x1,x2,x3,y1,y2,y3
     real :: dtheta,theta0
            
     call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
     call getCorners( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
                        
      val = sign(1.,cos(theta0)) * ( int_ddz_dx_dz_fct( x1, z2, y3, x, y, z, theta0 ) - int_ddz_dx_dz_fct( x3, z2, y3, x, y, z, theta0 ) - ( int_ddz_dx_dz_fct(x1,z1, y3, x, y, z, theta0) - int_ddz_dx_dz_fct(x3,z1, y3, x, y, z, theta0) ) )
    end subroutine
    
    subroutine int_ddz_dx_dz_inv(dat, val )
     real,intent(inout) :: val
     class(dataCollectionBase), intent(inout), target :: dat
     real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
     real :: x1,x2,x3,y1,y2,y3
     real :: dtheta,theta0
            
     call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
     call getCorners_inv( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
     !change the sign since the normal vector is now pointing along the positive y-axis
      val = -sign(1.,cos(theta0)) * ( int_ddz_dx_dz_fct( x3, z2, y3, x, y, z, theta0 ) - int_ddz_dx_dz_fct( x2, z2, y3, x, y, z, theta0 ) - ( int_ddz_dx_dz_fct(x3,z1, y3, x, y, z, theta0) - int_ddz_dx_dz_fct(x2,z1, y3, x, y, z, theta0) ) )
    end subroutine
        
    
    function int_ddz_dx_dz_fct( xp, zp, y3, x, y, z, theta0 )
        real,intent(in) :: xp,zp, y3, x, y, z, theta0
        real :: int_ddz_dx_dz_fct
        real :: arg
        
        arg = xp-x + P(x,y,z,xp,y3,zp)
        
        if ( arg .le. 10*tiny(1.) ) then
            arg = 10*tiny(1.)
        endif      
        
        int_ddz_dx_dz_fct =  -sign(1.,sin(theta0)) / (4*pi) * log( arg )
        
        
    end function
        
    subroutine int_ddx_dy_dz(dat, val )
     real,intent(inout) :: val
     class(dataCollectionBase), intent(inout), target :: dat
     real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
     real :: x1,x2,x3,y1,y2,y3
     real :: dtheta,theta0
            
     call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
     call getCorners( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
                                                     
     val = sign(1.,cos(theta0)) * sign(1.,sin(theta0)) * (int_ddx_dy_dz_fct( y2, z2, x3, x, y, z, theta0 ) - int_ddx_dy_dz_fct( y3, z2, x3, x, y, z, theta0 ) - ( int_ddx_dy_dz_fct(y2,z1, x3, x, y, z, theta0) - int_ddx_dy_dz_fct(y3,z1, x3, x, y, z, theta0) ))
     !val = sign(1.,cos(theta0)) * (int_ddx_dy_dz_fct( y2, z2, x3, x, y, z, theta0 ) - int_ddx_dy_dz_fct( y3, z2, x3, x, y, z, theta0 ) - ( int_ddx_dy_dz_fct(y2,z1, x3, x, y, z, theta0) - int_ddx_dy_dz_fct(y3,z1, x3, x, y, z, theta0) ))

    end subroutine
    
    !::inverted version of the circ piece integral
    subroutine int_ddx_dy_dz_inv(dat, val )
     real,intent(inout) :: val
     class(dataCollectionBase), intent(inout), target :: dat
     real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
     real :: x1,x2,x3,y1,y2,y3
     real :: dtheta,theta0
            
     call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
     call getCorners_inv( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
     !sign change as we are now pointing along the positive x-axis with the normal vector  
     val = -sign(1.,cos(theta0)) * sign(1.,sin(theta0)) * (int_ddx_dy_dz_fct( y3, z2, x3, x, y, z, theta0 ) - int_ddx_dy_dz_fct( y1, z2, x3, x, y, z, theta0 ) - ( int_ddx_dy_dz_fct(y3,z1, x3, x, y, z, theta0) - int_ddx_dy_dz_fct(y1,z1, x3, x, y, z, theta0) ))
     !val = int_ddx_dy_dz_fct( y3, z2, x3, x, y, z, theta0 )
    end subroutine
    
    function int_ddx_dy_dz_fct( yp, zp, x3, x, y, z, theta0 )
        real,intent(in) :: yp,zp, x3, x, y, z, theta0
        real :: int_ddx_dy_dz_fct,arg
     
        arg = sign(1.,sin(theta0)) * (y-yp) * (z-zp) / ((x-x3) * P(x,y,z,x3,yp,zp))
        
        
        !::Cover the limit when x-x3 goes to zero
    !    if ( abs( x-x3 ) .lt. 10*tiny(1.) .or. abs(y-yp) .lt. 1e-12 .or. P(x,y,z,x3,yp,zp) .lt. 1e-12) then
    !        int_ddx_dy_dz_fct = -1/(4*pi) * atan(  huge(1.) )
    !    else            
    !        int_ddx_dy_dz_fct =  -1/(4*pi) * atan( arg )
    !    endif
        int_ddx_dy_dz_fct = -1./(4.*pi)*atan( arg )
        
    end function
    
       subroutine int_ddy_dy_dz(dat, val )
         real,intent(inout) :: val
         class(dataCollectionBase), intent(inout), target :: dat
         real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
         real :: x1,x2,x3,y1,y2,y3
         real :: dtheta,theta0
            
         call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
         call getCorners( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
            
         val = -sign(1.,cos(theta0)) * sign(1.,sin(theta0)) * ( int_ddy_dy_dz_fct( y2, z2, x3, x, y, z, theta0 ) - int_ddy_dy_dz_fct( y3, z2, x3, x, y, z, theta0 ) - ( int_ddy_dy_dz_fct(y2,z1, x3, x, y, z, theta0) - int_ddy_dy_dz_fct(y3,z1, x3, x, y, z, theta0) ) )
         !val = -sign(1.,sin(theta0)) *  ( int_ddy_dy_dz_fct( y2, z2, x3, x, y, z, theta0 ) - int_ddy_dy_dz_fct( y3, z2, x3, x, y, z, theta0 ) - ( int_ddy_dy_dz_fct(y2,z1, x3, x, y, z, theta0) - int_ddy_dy_dz_fct(y3,z1, x3, x, y, z, theta0) ) )
         
            
       end subroutine
        
       subroutine int_ddy_dy_dz_inv(dat, val )
         real,intent(inout) :: val
         class(dataCollectionBase), intent(inout), target :: dat
         real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
         real :: x1,x2,x3,y1,y2,y3
         real :: dtheta,theta0
            
         call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
         call getCorners_inv( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
            !change the sign as the normal vector is pointing along the positive y-axis
         val = sign(1.,cos(theta0)) * sign(1.,sin(theta0)) * ( int_ddy_dy_dz_fct( y3, z2, x3, x, y, z, theta0 ) - int_ddy_dy_dz_fct( y1, z2, x3, x, y, z, theta0 ) - ( int_ddy_dy_dz_fct(y3,z1, x3, x, y, z, theta0) - int_ddy_dy_dz_fct(y1,z1, x3, x, y, z, theta0) ) )
            
       end subroutine
       
       function int_ddy_dy_dz_fct( yp, zp, x3, x, y, z, theta0 )
        real,intent(in) :: yp,zp, x3, x, y, z, theta0
        real :: int_ddy_dy_dz_fct,arg
        
        arg = zp-z + P(x,y,z,x3,yp,zp)
        
        if ( arg .le. 10*tiny(1.) ) then
            arg = 10*tiny(1.)
        endif
            
        int_ddy_dy_dz_fct =  1 / (4*pi) * log( arg )
        
        
       end function
    
       subroutine int_ddz_dy_dz(dat, val )
         real,intent(inout) :: val
         class(dataCollectionBase), intent(inout), target :: dat
         real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
         real :: x1,x2,x3,y1,y2,y3
         real :: dtheta,theta0
            
         call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
         call getCorners( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
                                 
         val = sign(1.,sin(theta0)) * ( int_ddz_dy_dz_fct( y2, z2, x3, x, y, z, theta0 ) - int_ddz_dy_dz_fct( y3, z2, x3, x, y, z, theta0 ) - ( int_ddz_dy_dz_fct(y2,z1, x3, x, y, z, theta0) - int_ddz_dy_dz_fct(y3,z1, x3, x, y, z, theta0) ) )
       end subroutine
       
       subroutine int_ddz_dy_dz_inv(dat, val )
         real,intent(inout) :: val
         class(dataCollectionBase), intent(inout), target :: dat
         real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
         real :: x1,x2,x3,y1,y2,y3
         real :: dtheta,theta0
            
         call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
         call getCorners_inv( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
         !change the sign since the normal vector is now pointing along the positive x-axis                  
         val = -sign(1.,sin(theta0)) * ( int_ddz_dy_dz_fct( y3, z2, x3, x, y, z, theta0 ) - int_ddz_dy_dz_fct( y1, z2, x3, x, y, z, theta0 ) - ( int_ddz_dy_dz_fct(y3,z1, x3, x, y, z, theta0) - int_ddz_dy_dz_fct(y1,z1, x3, x, y, z, theta0) ) )
       end subroutine
       
       function int_ddz_dy_dz_fct( yp, zp, x3, x, y, z, theta0 )
        real,intent(in) :: yp,zp, x3, x, y, z, theta0
        real :: int_ddz_dy_dz_fct,arg
        
        arg = yp-y + P(x,y,z,x3,yp,zp)
        
        if ( arg .le. 10*tiny(1.) ) then
            arg = 10*tiny(1.)
        endif
        
        
        int_ddz_dy_dz_fct = -sign(1.,cos(theta0)) / (4*pi) * log( arg )
        
        
       end function
       
    !::------------------------------------------------------------------------------------------------
    !::The integrals over the two end surfaces (z' = z1 and z' = z2) of the circular piece in closed form.
    !::For an end surface at height zp the integrals over the cross section A, which is bounded by the arc,
    !::the vertical chord at x' = xa = R cos(theta2) and the horizontal chord at y' = ya = R sin(theta1), are
    !::   capM(zp) = iint_A dD/dx dA' = [ ln( y'-y + P(xa,y',zp) ) ]_{ya}^{yb} + R [ cos(phi) Et - sin(phi) F ]_{psi1}^{psi2}
    !::   capN(zp) = iint_A dD/dy dA' = [ ln( x'-x + P(x',ya,zp) ) ]_{xa}^{xb} + R [ sin(phi) Et + cos(phi) F ]_{psi1}^{psi2}
    !::   capO(zp) = iint_A dD/dz dA' = sgn(dz) ( Omega_S - Omega_T1 - Omega_T2 )
    !::where D = 1/P is the inverse distance, dz = zp - z, xrot = sqrt(x**2+y**2) and phi = atan2(y,x) are the
    !::distance from the z-axis and the azimuth of the field point, psi = theta' - phi is the rotated angle on
    !::the arc with the limits psi1 = theta1 - phi and psi2 = theta2 - phi, Et and F are the z-derivative
    !::integrals of the arc surface (int_ddz_cos_dtheta_dz_fct with the continuous odd extension
    !::Et(psi) = sgn(psi) ( E(|psi|) - E(0) ) across psi = 0, and int_ddz_sin_dtheta_dz_fct), Omega_S is the solid angle subtended by the circular sector
    !::between theta1 and theta2 and Omega_T1, Omega_T2 are the solid angles of the two triangles (O,1,3)
    !::and (O,3,2) that complete the sector, with the corners 1 = (xb,ya), 2 = (xa,yb) and 3 = (xa,ya).
    !::See "The magnetic field of a homogeneously magnetized circular piece", Sec. 3.3.
    !::The routines return -cap(z1)/(4 pi) and -cap(z2)/(4 pi), so that e.g. Nxz = val2 - val1
    !::in the sign convention of MagTense (H = N.M). For the complementary (inverted) piece the chord
    !::coordinates xa -> xb and ya -> yb are interchanged in the logarithms, the solid angle of A is replaced
    !::by that of the bounding rectangle minus A, and the signs of capM and capN are reversed.
    !::The field point has been reflected into the first quadrant by getN_circPiece, so that
    !::0 < theta1 < theta2 < pi/2, xa < xb and ya < yb.
    !::------------------------------------------------------------------------------------------------
        !::capM and capN at both end surfaces: valx1 = -capM(z1)/(4 pi), valx2 = -capM(z2)/(4 pi), valy1 = -capN(z1)/(4 pi),
        !::valy2 = -capN(z2)/(4 pi), so that Nxz = valx2 - valx1 and Nyz = valy2 - valy1. The arc parts shared by capM
        !::and capN are evaluated once per end surface.
        subroutine int_ddxy_dx_dy(dat, valx1, valx2, valy1, valy2 )
         real,intent(inout) :: valx1, valx2, valy1, valy2
         class(dataCollectionBase), intent(inout), target :: dat
         real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi,capM,capN

         call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )

         call circPiece_capMN( R, theta1, theta2, x, y, z, z1, xrot, phi, .false., capM, capN )
         valx1 = -capM / (4*pi)
         valy1 = -capN / (4*pi)
         call circPiece_capMN( R, theta1, theta2, x, y, z, z2, xrot, phi, .false., capM, capN )
         valx2 = -capM / (4*pi)
         valy2 = -capN / (4*pi)

        end subroutine

        !::for the radially inwards pointing (complementary) circ piece
        subroutine int_ddxy_dx_dy_inv(dat, valx1, valx2, valy1, valy2 )
         real,intent(inout) :: valx1, valx2, valy1, valy2
         class(dataCollectionBase), intent(inout), target :: dat
         real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi,capM,capN

         call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )

         call circPiece_capMN( R, theta1, theta2, x, y, z, z1, xrot, phi, .true., capM, capN )
         valx1 = -capM / (4*pi)
         valy1 = -capN / (4*pi)
         call circPiece_capMN( R, theta1, theta2, x, y, z, z2, xrot, phi, .true., capM, capN )
         valx2 = -capM / (4*pi)
         valy2 = -capN / (4*pi)

        end subroutine

        !::The logarithmic term of capM at (y',z') on the vertical chord x' = x3
        function int_ddx_dx_dy_fct1( yp, zp, x3, x, y, z )
        real,intent(in) :: yp,zp, x3, x, y, z
        real :: int_ddx_dx_dy_fct1,arg

            arg = yp-y + P(x,y,z,x3,yp,zp)
            if ( arg .le. 10*tiny(1.) ) then
                arg = 10*tiny(1.)
            endif

            int_ddx_dx_dy_fct1 = log( arg )

        end function

        !::The logarithmic term of capN at (x',z') on the horizontal chord y' = y3
        function int_ddy_dx_dy_fct1( xp, zp, y3, x, y, z )
        real,intent(in) :: xp,zp, y3, x, y, z
        real :: int_ddy_dx_dy_fct1,arg

            arg = xp-x + P(x,y,z,xp,y3,zp)

            if ( arg .le. 10*tiny(1.) ) then
                arg = 10*tiny(1.)
            endif

            int_ddy_dx_dy_fct1 = log( arg )

        end function

         subroutine int_ddz_dx_dy(dat, val1, val2 )
         real,intent(inout) :: val1,val2
         class(dataCollectionBase), intent(inout), target :: dat
         real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi

         call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )

         val1 = -circPiece_capO( R, theta1, theta2, x, y, z, z1, xrot, phi, .false. ) / (4*pi)
         val2 = -circPiece_capO( R, theta1, theta2, x, y, z, z2, xrot, phi, .false. ) / (4*pi)

         end subroutine

         subroutine int_ddz_dx_dy_inv(dat, val1, val2 )
         real,intent(inout) :: val1,val2
         class(dataCollectionBase), intent(inout), target :: dat
         real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi

         call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )

         val1 = -circPiece_capO( R, theta1, theta2, x, y, z, z1, xrot, phi, .true. ) / (4*pi)
         val2 = -circPiece_capO( R, theta1, theta2, x, y, z, z2, xrot, phi, .true. ) / (4*pi)

         end subroutine

    !::capM(zp) and capN(zp) for the piece (inv = .false.) and for the complementary piece (inv = .true.)
    subroutine circPiece_capMN( R, theta1, theta2, x, y, z, zp, xrot, phi, inv, capM, capN )
    real,intent(in) :: R, theta1, theta2, x, y, z, zp, xrot, phi
    logical,intent(in) :: inv
    real,intent(out) :: capM, capN
    real :: xa, xb, ya, yb, xc, yc, arcM, arcN
        xa = R * cos( theta2 )
        xb = R * cos( theta1 )
        ya = R * sin( theta1 )
        yb = R * sin( theta2 )
        call circPiece_arcMN( R, theta1, theta2, zp - z, xrot, phi, arcM, arcN )
        if ( inv ) then
            xc = xb
            yc = yb
        else
            xc = xa
            yc = ya
        endif
        capM = int_ddx_dx_dy_fct1( yb, zp, xc, x, y, z ) - int_ddx_dx_dy_fct1( ya, zp, xc, x, y, z ) + arcM
        capN = int_ddy_dx_dy_fct1( xb, zp, yc, x, y, z ) - int_ddy_dx_dy_fct1( xa, zp, yc, x, y, z ) + arcN
        if ( inv ) then
            capM = -capM
            capN = -capN
        endif
    end subroutine circPiece_capMN

    !::capO(zp) = sgn(dz) Omega_A for the piece and sgn(dz) ( Omega_rectangle - Omega_A ) for the complementary
    !::piece, with Omega_A = Omega_S - Omega_T1 - Omega_T2 the solid angle subtended by the cross section.
    !::For a field point in the plane of the end surface (dz = 0) the solid angle is discontinuous, and the
    !::mean of the two limits, zero, is returned.
    function circPiece_capO( R, theta1, theta2, x, y, z, zp, xrot, phi, inv )
    real,intent(in) :: R, theta1, theta2, x, y, z, zp, xrot, phi
    logical,intent(in) :: inv
    real :: circPiece_capO
    real :: dz, xa, xb, ya, yb, psi1, psi2, OmS, OmT, OmR
    real,dimension(3) :: pO, p1, p2, p3, p4
        dz = zp - z
        if ( abs( dz ) .lt. 10*tiny(1.) ) then
            circPiece_capO = 0.
            return
        endif
        xa = R * cos( theta2 )
        xb = R * cos( theta1 )
        ya = R * sin( theta1 )
        yb = R * sin( theta2 )
        call circPiece_psi_limits( theta1, theta2, phi, psi1, psi2 )
        OmS = circPiece_omega_sector( R, xrot, psi1, psi2, dz )
        !::positions of the origin and of the corners 1, 2, 3 of the cross section relative to the field point
        pO = (/ -x, -y, dz /)
        p1 = (/ xb-x, ya-y, dz /)
        p2 = (/ xa-x, yb-y, dz /)
        p3 = (/ xa-x, ya-y, dz /)
        OmT = circPiece_omega_triangle( pO, p1, p3 ) + circPiece_omega_triangle( pO, p3, p2 )
        if ( inv ) then
            !::the bounding rectangle with the fourth corner 3' = (xb,yb), split into the triangles (3,1,3') and (3,3',2)
            p4 = (/ xb-x, yb-y, dz /)
            OmR = circPiece_omega_triangle( p3, p1, p4 ) + circPiece_omega_triangle( p3, p4, p2 )
            circPiece_capO = sign( 1., dz ) * ( OmR - ( OmS - OmT ) )
        else
            circPiece_capO = sign( 1., dz ) * ( OmS - OmT )
        endif
    end function circPiece_capO

    !::The limits psi1 = theta1 - phi and psi2 = theta2 - phi of the rotated angle, shifted by a multiple
    !::of 2 pi so that -pi <= psi1 < pi (the integrands are periodic in psi, and the odd extensions used
    !::below are valid for -2 pi < psi < 2 pi).
    subroutine circPiece_psi_limits( theta1, theta2, phi, psi1, psi2 )
    real,intent(in) :: theta1, theta2, phi
    real,intent(out) :: psi1, psi2
    integer :: m
        psi1 = theta1 - phi
        psi2 = theta2 - phi
        m = floor( ( psi1 + pi ) / ( 2*pi ) )
        psi1 = psi1 - 2*pi*m
        psi2 = psi2 - 2*pi*m
    end subroutine circPiece_psi_limits

    !::The arc parts of capM and capN, arcM = R [ cos(phi) Et - sin(phi) F ] and arcN = R [ sin(phi) Et + cos(phi) F ]
    !::between psi1 and psi2, which equal minus the integrals of the inverse distance along the arc,
    !::   arcM = -int_{ya}^{yb} D( sqrt(R**2-y'**2), y', zp ) dy',   arcN = -int_{xa}^{xb} D( x', sqrt(R**2-x'**2), zp ) dx'.
    !::Close to the axis the integrals are expanded to first order in xrot, with A0 = sqrt(R**2+dz**2),
    !::   int D dy' = R/A0 [ sin(psi+phi) ] + R**2 xrot/A0**3 ( cos(phi) [ psi/2 + sin(2psi)/4 ] - sin(phi) [ sin(psi)**2/2 ] ),
    !::   int D dx' = -R/A0 [ cos(psi+phi) ] + R**2 xrot/A0**3 ( sin(phi) [ psi/2 + sin(2psi)/4 ] + cos(phi) [ sin(psi)**2/2 ] ).
    subroutine circPiece_arcMN( R, theta1, theta2, dz, xrot, phi, arcM, arcN )
    real,intent(in) :: R, theta1, theta2, dz, xrot, phi
    real,intent(out) :: arcM, arcN
    real :: psi1, psi2, E0, Et1, Et2, F1, F2, A0, t_c2, t_sc, dS, dC
        call circPiece_psi_limits( theta1, theta2, phi, psi1, psi2 )
        if ( xrot .lt. xrot_cap_tol * R ) then
            A0 = sqrt( R**2 + dz**2 )
            dS = sin( psi2 + phi ) - sin( psi1 + phi )
            dC = cos( psi2 + phi ) - cos( psi1 + phi )
            t_c2 = ( psi2/2 + sin(2*psi2)/4 ) - ( psi1/2 + sin(2*psi1)/4 )
            t_sc = ( sin(psi2)**2 - sin(psi1)**2 ) / 2
            arcM = -( R/A0 * dS + R**2 * xrot / A0**3 * ( cos(phi) * t_c2 - sin(phi) * t_sc ) )
            arcN = -( -R/A0 * dC + R**2 * xrot / A0**3 * ( sin(phi) * t_c2 + cos(phi) * t_sc ) )
        else
            E0 = int_ddz_cos_dtheta_dz_fct( 0., dz, R, xrot )
            Et1 = sign( 1., psi1 ) * ( int_ddz_cos_dtheta_dz_fct( abs(psi1), dz, R, xrot ) - E0 )
            Et2 = sign( 1., psi2 ) * ( int_ddz_cos_dtheta_dz_fct( abs(psi2), dz, R, xrot ) - E0 )
            F1 = int_ddz_sin_dtheta_dz_fct( psi1, dz, R, xrot )
            F2 = int_ddz_sin_dtheta_dz_fct( psi2, dz, R, xrot )
            arcM = R * ( cos(phi) * ( Et2 - Et1 ) - sin(phi) * ( F2 - F1 ) )
            arcN = R * ( sin(phi) * ( Et2 - Et1 ) + cos(phi) * ( F2 - F1 ) )
        endif
    end subroutine circPiece_arcMN

    !::Solid angle of the circular sector 0 <= r' <= R, psi1 <= psi <= psi2 at height dz, seen from the field
    !::point at distance xrot from the axis (manuscript Eq. OmS),
    !::   Omega_S = [ Sc( psi, rho/|dz| ) + Pt( psi ) ]_{psi1}^{psi2},   rho = sqrt( xrot**2 + dz**2 ),
    !::with Sc the continuous antiderivative of the elementary part (circPiece_Sc) and Pt the odd extension
    !::Pt(psi) = sgn(psi) ( Pi(|psi|) - Pi(0) ) of the part containing the elliptic integrals of the third kind,
    !::   Pi(psi) = -|dz| / ( 2 xrot rho sqrt( (R+xrot)**2 + dz**2 ) ) * ( (alpha + beta nu_p) Pi( D, nu_p, k ) - (alpha + beta nu_m) Pi( D, nu_m, k ) )
    !::with D = cos(psi/2), k**2 = 4 R xrot / ( (R+xrot)**2 + dz**2 ), nu_p = 2 xrot/(xrot+rho), nu_m = 2 xrot/(xrot-rho),
    !::alpha = 2 xrot R and beta = -( xrot R + rho**2 ).
    !::Close to the plane of the end surface nu_p -> 1 and nu_m -> -infinity, and the two differences
    !::Pi( D, nu, k ) - Pi( 1, nu, k ) are then formed without cancellation by circPiece_dPi_p and circPiece_dPi_m.
    !::|dz| is kept above 1e-30 xrot so that none of the intermediate quantities underflows or overflows.
    function circPiece_omega_sector( R, xrot, psi1, psi2, dz_in )
    real,intent(in) :: R, xrot, psi1, psi2, dz_in
    real :: circPiece_omega_sector
    real :: dz, rho, c, k2, nu_p, nu_m, nup, one_minus_nu_p, alpha, beta, pref, A, rf1, rjp1, rjm1
    real :: psi, s, cc, q, dPp, dPm
    real,dimension(2) :: Pt
    integer :: i
        if ( xrot .lt. xrot_sector_tol * R ) then
            circPiece_omega_sector = ( psi2 - psi1 ) * ( 1. - abs(dz_in) / sqrt( R**2 + dz_in**2 ) )
            return
        endif
        dz = sign( max( abs( dz_in ), 1e-30 * xrot ), dz_in )
        rho = sqrt( xrot**2 + dz**2 )
        c = rho / abs( dz )
        k2 = 4 * R * xrot / ( ( R + xrot )**2 + dz**2 )
        !::1 - nu_p = (rho - xrot)/(rho + xrot) = dz**2/(rho + xrot)**2 is formed without cancellation
        one_minus_nu_p = dz**2 / ( rho + xrot )**2
        nu_p = 1. - one_minus_nu_p
        nu_m = -2 * xrot * ( rho + xrot ) / dz**2
        !::reciprocal characteristic k**2/nu_m and the coefficient of the arctangent in circPiece_dPi_m
        nup = k2 / nu_m
        A = sqrt( nu_m / ( ( 1. - nu_m ) * ( nu_m - k2 ) ) )
        alpha = 2 * xrot * R
        beta = -( xrot * R + rho**2 )
        pref = -abs( dz ) / ( 2 * xrot * rho * sqrt( ( R + xrot )**2 + dz**2 ) )
        !::complete integrals (D = 1), shared by both limits
        rf1 = rf( 0., 1.-k2, 1. )
        rjp1 = rj( 0., 1.-k2, 1., one_minus_nu_p )
        rjm1 = rj( 0., 1.-k2, 1., 1. - nup )
        do i = 1, 2
            if ( i .eq. 1 ) then
                psi = psi1
            else
                psi = psi2
            endif
            s = cos( abs( psi ) / 2 )
            cc = sin( abs( psi ) / 2 )**2
            q = 1. - k2 * s**2
            dPp = circPiece_dPi_p( s, cc, q, nu_p, one_minus_nu_p, rf1, rjp1 )
            dPm = circPiece_dPi_m( s, cc, q, nup, A, rjm1 )
            Pt(i) = sign( 1., psi ) * pref * ( ( alpha + beta * nu_p ) * dPp - ( alpha + beta * nu_m ) * dPm )
        enddo
        circPiece_omega_sector = ( circPiece_Sc( psi2, c ) - circPiece_Sc( psi1, c ) ) + ( Pt(2) - Pt(1) )
    end function circPiece_omega_sector

    !::Continuous antiderivative of c / ( cos(psi)**2 + c**2 sin(psi)**2 ): Sc = m pi + atan( c tan( psi - m pi ) ), m = floor( psi/pi + 1/2 )
    function circPiece_Sc( psi, c )
    real,intent(in) :: psi, c
    real :: circPiece_Sc
    integer :: m
        m = floor( psi/pi + 0.5 )
        circPiece_Sc = m*pi + atan( c * tan( psi - m*pi ) )
    end function circPiece_Sc

    !::Pi( s, nu, k ) - Pi( 1, nu, k ) for 0 < nu < 1, in the Maple convention (s = sine of the amplitude, the
    !::characteristic entering as 1 - nu t**2, k**2 = k2), from Carlson's functions as in ellpi of SpecialFunctions:
    !::   Pi( s, nu, k ) = s RF( 1-s**2, 1-k2 s**2, 1 ) + nu s**3 RJ( 1-s**2, 1-k2 s**2, 1, 1-nu s**2 ) / 3 .
    !::The last argument of RJ is formed as (1-nu) + nu (1-s**2) from the separately supplied 1 - nu, which keeps
    !::it accurate when nu is close to 1. cc = 1-s**2, q = 1-k2 s**2 and the complete values rf1 = RF( 0, 1-k2, 1 )
    !::and rj1 = RJ( 0, 1-k2, 1, 1-nu ) are supplied by the caller.
    function circPiece_dPi_p( s, cc, q, nu, one_minus_nu, rf1, rj1 )
    real,intent(in) :: s, cc, q, nu, one_minus_nu, rf1, rj1
    real :: circPiece_dPi_p
        circPiece_dPi_p = s * rf( cc, q, 1. ) - rf1 + nu / 3. * ( s**3 * rj( cc, q, 1., one_minus_nu + nu * cc ) - rj1 )
    end function circPiece_dPi_p

    !::Pi( s, nu, k ) - Pi( 1, nu, k ) for nu < 0, by the reciprocal-characteristic transformation
    !::   Pi( s, nu, k ) = F( s, k ) - Pi( s, k2/nu, k ) + A atan( s / ( A sqrt(1-s**2) sqrt(1-k2 s**2) ) ),   A = sqrt( nu / ( (1-nu)(nu-k2) ) ),
    !::in which F - Pi( s, nu', k ) = -nu' s**3 RJ( 1-s**2, 1-k2 s**2, 1, 1-nu' s**2 ) / 3 with the small characteristic
    !::nu' = k2/nu, and in which the difference of the arctangents between s and 1 is formed explicitly. This avoids
    !::the cancellation between F and Pi that occurs for large |nu|, where Pi( s, nu, k ) = O(1/|nu|).
    !::nup = nu', A and the complete value rj1 = RJ( 0, 1-k2, 1, 1-nu' ) are supplied by the caller.
    function circPiece_dPi_m( s, cc, q, nup, A, rj1 )
    real,intent(in) :: s, cc, q, nup, A, rj1
    real :: circPiece_dPi_m
    real :: datan
        if ( s .ge. 0. ) then
            datan = -atan2( A * sqrt( cc ) * sqrt( q ), s )
        else
            datan = -pi + atan2( A * sqrt( cc ) * sqrt( q ), abs( s ) )
        endif
        circPiece_dPi_m = -nup / 3. * ( s**3 * rj( cc, q, 1., 1. - nup * s**2 ) - rj1 ) + A * datan
    end function circPiece_dPi_m

    !::Solid angle of the plane triangle with the vertices p1, p2, p3 relative to the field point
    !::(Van Oosterom and Strackee, IEEE Trans. Biomed. Eng. 30 (1983) 125)
    function circPiece_omega_triangle( p1, p2, p3 )
    real,dimension(3),intent(in) :: p1, p2, p3
    real :: circPiece_omega_triangle
    real :: n1, n2, n3, num, den
    real,dimension(3) :: c23
        n1 = sqrt( sum( p1**2 ) )
        n2 = sqrt( sum( p2**2 ) )
        n3 = sqrt( sum( p3**2 ) )
        c23 = (/ p2(2)*p3(3) - p2(3)*p3(2), p2(3)*p3(1) - p2(1)*p3(3), p2(1)*p3(2) - p2(2)*p3(1) /)
        num = abs( sum( p1 * c23 ) )
        den = n1*n2*n3 + sum( p1*p2 )*n3 + sum( p1*p3 )*n2 + sum( p2*p3 )*n1
        circPiece_omega_triangle = 2 * atan2( num, den )
    end function circPiece_omega_triangle


    !::Small-xrot expansion of the six double integrals over the curved face, to first order in xrot.
    !::With A = R**2 + z**2 the integrands are expanded as 1/M**3 = A**(-3/2) (1 + 3 R xrot cos(t)/A) + O(xrot**2),
    !::after which theta and z separate. With the z-integrals between z1 and z2
    !::   Z3 = [ z / (R**2 sqrt(A)) ],  Z5 = [ z (3R**2 + 2z**2) / (3 R**4 A**(3/2)) ],  Zm1 = [ A**(-1/2) ],  Zm3 = [ A**(-3/2) ]
    !::and the theta-integrals between theta1 and theta2
    !::   c2 = [ t/2 + sin(2t)/4 ] (int cos**2),  s2 = [ t/2 - sin(2t)/4 ] (int sin**2),  sc = [ sin**2/2 ] (int sin cos),
    !::   S = [ sin ] (int cos),  C = [ cos ] (-int sin),  c3 = [ sin - sin**3/3 ] (int cos**3),
    !::   C3 = [ cos**3 ] (-3 int sin cos**2),  S3 = [ sin**3 ] (3 int sin**2 cos),
    !::the integrals are
    !::   1: int cos(t) (xrot - R cos(t)) / M**3  ->  -R c2 Z3 + xrot ( S Z3 - 3 R**2 c3 Z5 )
    !::   2: int sin(t) (xrot - R cos(t)) / M**3  ->  -R sc Z3 + xrot ( -C Z3 + R**2 C3 Z5 )
    !::   3: int cos(t) R sin(t) / M**3           ->   R sc Z3 - xrot R**2 C3 Z5
    !::   4: int sin(t) R sin(t) / M**3           ->   R s2 Z3 + xrot R**2 S3 Z5
    !::   5: int cos(t) z / M**3                  ->  -S Zm1 - xrot R c2 Zm3
    !::   6: int sin(t) z / M**3                  ->   C Zm1 - xrot R sc Zm3
    !::The neglected terms are of relative order 2 (xrot/R)**2. The result does not depend on the branch
    !::of theta, so no wrapping of the angles is needed.
    function int_axis_dtheta_dz( which, R, xrot, theta1, theta2, z1, z2 )
    integer,intent(in) :: which
    real,intent(in) :: R, xrot, theta1, theta2, z1, z2
    real :: int_axis_dtheta_dz
    real :: dZ3, dZ5, dZm1, dZm3, t_c2, t_s2, t_sc, t_s, t_c, t_c3, t_cc3, t_ss3
        dZ3  = z2 / ( R**2 * sqrt( R**2 + z2**2 ) ) - z1 / ( R**2 * sqrt( R**2 + z1**2 ) )
        dZ5  = z2 * ( 3*R**2 + 2*z2**2 ) / ( 3*R**4 * ( R**2 + z2**2 )**1.5 ) - z1 * ( 3*R**2 + 2*z1**2 ) / ( 3*R**4 * ( R**2 + z1**2 )**1.5 )
        dZm1 = 1. / sqrt( R**2 + z2**2 ) - 1. / sqrt( R**2 + z1**2 )
        dZm3 = 1. / ( R**2 + z2**2 )**1.5 - 1. / ( R**2 + z1**2 )**1.5
        t_c2  = ( theta2/2 + sin(2*theta2)/4 ) - ( theta1/2 + sin(2*theta1)/4 )
        t_s2  = ( theta2/2 - sin(2*theta2)/4 ) - ( theta1/2 - sin(2*theta1)/4 )
        t_sc  = ( sin(theta2)**2 - sin(theta1)**2 ) / 2
        t_s   = sin(theta2) - sin(theta1)
        t_c   = cos(theta2) - cos(theta1)
        t_c3  = ( sin(theta2) - sin(theta2)**3/3 ) - ( sin(theta1) - sin(theta1)**3/3 )
        t_cc3 = cos(theta2)**3 - cos(theta1)**3
        t_ss3 = sin(theta2)**3 - sin(theta1)**3
        select case ( which )
        case ( 1 )
            int_axis_dtheta_dz = -R * t_c2 * dZ3 + xrot * ( t_s * dZ3 - 3*R**2 * t_c3 * dZ5 )
        case ( 2 )
            int_axis_dtheta_dz = -R * t_sc * dZ3 + xrot * ( -t_c * dZ3 + R**2 * t_cc3 * dZ5 )
        case ( 3 )
            int_axis_dtheta_dz =  R * t_sc * dZ3 - xrot * R**2 * t_cc3 * dZ5
        case ( 4 )
            int_axis_dtheta_dz =  R * t_s2 * dZ3 + xrot * R**2 * t_ss3 * dZ5
        case ( 5 )
            int_axis_dtheta_dz = -t_s * dZm1 - xrot * R * t_c2 * dZm3
        case ( 6 )
            int_axis_dtheta_dz =  t_c * dZm1 - xrot * R * t_sc * dZm3
        case default
            int_axis_dtheta_dz = 0.
        end select
    end function int_axis_dtheta_dz
    
    
    subroutine getParameters_rot_trick( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
    class(dataCollectionBase), intent(inout), target :: dat
    real,intent(inout) :: x,y,z,R,theta1,theta2,z1,z2, xrot, phi
    
    !radius of the circle piece (distance to the curved face from origo)
    R = dat%r2
    !integration limits of theta
    theta1 = dat%theta1
    theta2 = dat%theta2
    
        
    !integration limits of z
    z1 = dat%z1
    z2 = dat%z2
    
    !coordinates at which the solution is sought
    x = dat%x
    y = dat%y
    z = dat%z
    !do the "change integration limits" trick
    xrot = sqrt( x**2 + y**2 )

    !phi = atan2(y,x)
    phi = atan2_custom( y, x )
        
    theta1 = theta1 - phi
    theta2 = theta2 - phi
        
    !Translate z limits to eliminate the z-coordinate        
    z1 = z1 - z
    z2 = z2 - z            
            
    end subroutine getParameters_rot_trick
    
    
    subroutine getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
    class(dataCollectionBase), intent(inout), target :: dat
    real,intent(inout) :: x,y,z,R,theta1,theta2,z1,z2, xrot, phi
    
    !radius of the circle piece (distance to the curved face from origo)
    R = dat%r2
    !integration limits of theta
    theta1 = dat%theta1
    theta2 = dat%theta2
    
        
    !integration limits of z
    z1 = dat%z1
    z2 = dat%z2
    
    !coordinates at which the solution is sought
    x = dat%x
    y = dat%y
    z = dat%z
    !do the "change integration limits" trick
    xrot = sqrt( x**2 + y**2 )

    phi = atan2_custom( y, x )
        
            
    end subroutine getParameters
    
    subroutine getCorners( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
    real,intent(in) :: R, theta1, theta2
    real,intent(inout) :: theta0, dtheta, x1, x2, x3, y1, y2, y3
    real :: a,b
    dtheta = theta2-theta1
    theta0 = theta1 + dtheta/2
    if ( cos(theta0) .ge. 0 .AND. sin(theta0).ge.0 ) then
        !first quadrant
        a = theta0-dtheta/2
        b = theta0+dtheta/2            
    elseif ( cos(theta0) .lt.0 .AND. sin(theta0) .ge.0 ) then
        !second quadrant
        a = theta0+dtheta/2
        b = theta0-dtheta/2
    elseif ( cos(theta0) .lt.0 .AND. sin(theta0) .lt.0 ) then
        !third
        a = theta0-dtheta/2
        b = theta0+dtheta/2
    elseif ( cos(theta0) .ge.0 .AND. sin(theta0) .lt.0 ) then
        !fourth 
        a = theta0+dtheta/2
        b = theta0-dtheta/2
    endif
        x1 = R * cos( a )
        y1 = R * sin( a )

        x2 = R * cos( b )
        y2 = R * sin( b )

        x3 = x2
        y3 = y1
                         
    end subroutine
    
    !::Returns the corners of the inverted circ piece, i.e. one that is pointing radially inwards
     subroutine getCorners_inv( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
    real,intent(in) :: R, theta1, theta2
    real,intent(inout) :: theta0, dtheta, x1, x2, x3, y1, y2, y3
    real :: a,b
    dtheta = theta2-theta1
    theta0 = theta1 + dtheta/2
    if ( cos(theta0) .ge. 0 .AND. sin(theta0).ge.0 ) then
        !first quadrant
        a = theta0-dtheta/2
        b = theta0+dtheta/2            
    elseif ( cos(theta0) .lt.0 .AND. sin(theta0) .ge.0 ) then
        !second quadrant
        a = theta0+dtheta/2
        b = theta0-dtheta/2
    elseif ( cos(theta0) .lt.0 .AND. sin(theta0) .lt.0 ) then
        !third
        a = theta0-dtheta/2
        b = theta0+dtheta/2
    elseif ( cos(theta0) .ge.0 .AND. sin(theta0) .lt.0 ) then
        !fourth 
        a = theta0+dtheta/2
        b = theta0-dtheta/2
    endif
        x1 = R * cos( a )
        y1 = R * sin( a )

        x2 = R * cos( b )
        y2 = R * sin( b )

        x3 = x1
        y3 = y2
                         
    end subroutine
    
end module TileCircPieceTensor
    