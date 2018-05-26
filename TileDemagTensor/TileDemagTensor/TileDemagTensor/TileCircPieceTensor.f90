module TileCircPieceTensor
    use QuadPack
    use SpecialFunctions
    use TileTensorHelperFunctions
    
    implicit none
    
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
        real :: elf,elE,elPi
        
        elf = ellf( C_no_sign(thetap), K(R,xrot,zp) )
        elE = elle( C_no_sign(thetap), K(R,xrot,zp) )
        elPi = ellpi( C_no_sign(thetap), B(R,xrot), K(R,xrot,zp) )
        
        int_ddz_cos_dtheta_dz_fct = -1 / (R*xrot*sqrt((R+xrot)**2+zp**2)) * ( ((R+xrot)**2+zp**2) * elE - (R**2 + xrot**2 + zp**2) * elF)
                
    end function
    
    subroutine int_ddz_sin_dtheta_dz( dat, val )
     real,intent(inout) :: val
     class(dataCollectionBase), intent(inout), target :: dat
     real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
     real :: dtheta
        
        val = 0.
           
        call getParameters_rot_trick( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )    
            
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

    end subroutine
    
    function int_ddx_dy_dz_fct( yp, zp, x3, x, y, z, theta0 )
        real,intent(in) :: yp,zp, x3, x, y, z, theta0
        real :: int_ddx_dy_dz_fct
     
        !::Cover the limit when x-x3 goes to zero
        if ( abs( x-x3 ) .lt. 10*tiny(1.) ) then
            int_ddx_dy_dz_fct = -1/(4*pi) * atan(  huge(1.) )
        else            
            int_ddx_dy_dz_fct =  -1/(4*pi) * atan( sign(1.,sin(theta0)) * (y-yp) * (z-zp) / ((x-x3) * P(x,y,z,x3,yp,zp)) )
        endif
        
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
       
        subroutine int_ddx_dx_dy(dat, val1, val2 )
         real,intent(inout) :: val1, val2
         class(dataCollectionBase), intent(inout), target :: dat
         real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
         real :: x1,x2,x3,y1,y2,y3
         real :: dtheta,theta0
         real :: val11,val12,val21,val22
         procedure (f_int_dat), pointer :: f_ptr => null ()
            
         call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
         call getCorners( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
                                    
                        
            val11 = int_ddx_dx_dy_fct1(y2,z1, x3, x, y, z) - int_ddx_dx_dy_fct1(y1,z1, x3, x, y, z)
            val12 = int_ddx_dx_dy_fct1(y2,z2, x3, x, y, z) - int_ddx_dx_dy_fct1(y1,z2, x3, x, y, z)
            
            dat%x = x
            dat%y = y
            dat%z = z
            dat%thetas = theta0
            dat%rs = R
            dat%zs = z1
            f_ptr => int_ddx_dx_dy_fct2
            call qags_x ( f_ptr, dat, y1, y2, dat%epsabs, dat%epsrel, val21, dat%abserr_x, dat%neval_x, dat%ier_x )
        
            dat%zs = z2
            call qags_x ( f_ptr, dat, y1, y2, dat%epsabs, dat%epsrel, val22, dat%abserr_x, dat%neval_x, dat%ier_x )
            
            
            val1 = sign(1.,sin(theta0))*sign(1.,cos(theta0))/(4*pi) * ( val21 - val11 )
            
            val2 = sign(1.,sin(theta0))*sign(1.,cos(theta0))/(4*pi) * ( val22 - val12 )
            
        end subroutine
        !::for the radially inwards pointing circ piece
        subroutine int_ddx_dx_dy_inv(dat, val1, val2 )
         real,intent(inout) :: val1, val2
         class(dataCollectionBase), intent(inout), target :: dat
         real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
         real :: x1,x2,x3,y1,y2,y3
         real :: dtheta,theta0
         real :: val11,val12,val21,val22
         procedure (f_int_dat), pointer :: f_ptr => null ()
            
         call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
         call getCorners_inv( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
                                    
                        
          val11 = int_ddx_dx_dy_fct1(y3,z1, x3, x, y, z) - int_ddx_dx_dy_fct1(y1,z1, x3, x, y, z)
          val12 = int_ddx_dx_dy_fct1(y3,z2, x3, x, y, z) - int_ddx_dx_dy_fct1(y1,z2, x3, x, y, z)
            
          dat%x = x
          dat%y = y
          dat%z = z
          dat%thetas = theta0
          dat%rs = R
          dat%zs = z1
          f_ptr => int_ddx_dx_dy_fct2
          call qags_x ( f_ptr, dat, y1, y2, dat%epsabs, dat%epsrel, val21, dat%abserr_x, dat%neval_x, dat%ier_x )
        
          dat%zs = z2
          call qags_x ( f_ptr, dat, y1, y2, dat%epsabs, dat%epsrel, val22, dat%abserr_x, dat%neval_x, dat%ier_x )
            
            
          val1 = sign(1.,sin(theta0))*sign(1.,cos(theta0))/(4*pi) * ( val11 - val21 )
          
          val2 = sign(1.,sin(theta0))*sign(1.,cos(theta0))/(4*pi) * ( val12 - val22 )
            
        end subroutine
        
        
        function int_ddx_dx_dy_fct1( yp, zp, x3, x, y, z )
        real,intent(in) :: yp,zp, x3, x, y, z
        real :: int_ddx_dx_dy_fct1,arg
        
            arg = yp-y + P(x,y,z,x3,yp,zp)
            if ( arg .le. 10*tiny(1.) ) then
                arg = 10*tiny(1.)
            endif     
                
            int_ddx_dx_dy_fct1 = log( arg )
            
            
        end function
        
        
        function int_ddx_dx_dy_fct2( yp, dat )              
        real,intent(in) :: yp
        class(dataCollectionBase), target :: dat
        real :: x,y,z,theta0,R,zp
        real :: int_ddx_dx_dy_fct2
            x = dat%x
            y = dat%y
            z = dat%z
            theta0 = dat%thetas
            R = dat%rs
            zp = dat%zs
            int_ddx_dx_dy_fct2 = 1./ (P(x,y,z,sign(1.,cos(theta0))*sqrt(R**2-yp**2),yp,zp))
        end function 
       
        subroutine int_ddy_dx_dy(dat, val1, val2 )
         real,intent(inout) :: val1, val2
         class(dataCollectionBase), intent(inout), target :: dat
         real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
         real :: x1,x2,x3,y1,y2,y3
         real :: dtheta,theta0
         real :: val11,val12,val21,val22
         procedure (f_int_dat), pointer :: f_ptr => null ()
            
         call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
         call getCorners( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
         
         
            val11 = int_ddy_dx_dy_fct1(x1,z1, y3, x, y, z) - int_ddy_dx_dy_fct1(x3,z1, y3, x, y, z)
            val12 = int_ddy_dx_dy_fct1(x1,z2, y3, x, y, z) - int_ddy_dx_dy_fct1(x3,z2, y3, x, y, z)
            
            dat%x = x
            dat%y = y
            dat%z = z
            dat%thetas = theta0
            dat%rs = R
            dat%zs = z1
            f_ptr => int_ddy_dx_dy_fct2
            call qags_x ( f_ptr, dat, x3, x1, dat%epsabs, dat%epsrel, val21, dat%abserr_x, dat%neval_x, dat%ier_x )
        
            dat%zs = z2
            call qags_x ( f_ptr, dat, x3, x1, dat%epsabs, dat%epsrel, val22, dat%abserr_x, dat%neval_x, dat%ier_x )
            
            
            val1 = sign(1.,sin(theta0))*sign(1.,cos(theta0))/(4*pi) * ( val21 - val11 )
            
            val2 = sign(1.,sin(theta0))*sign(1.,cos(theta0))/(4*pi) * ( val22 - val12 )
         
           
        end subroutine
        
        subroutine int_ddy_dx_dy_inv(dat, val1, val2 )
         real,intent(inout) :: val1, val2
         class(dataCollectionBase), intent(inout), target :: dat
         real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
         real :: x1,x2,x3,y1,y2,y3
         real :: dtheta,theta0
         real :: val11,val12,val21,val22
         procedure (f_int_dat), pointer :: f_ptr => null ()
            
         call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
         call getCorners_inv( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
         
         
            val11 = int_ddy_dx_dy_fct1(x3,z1, y3, x, y, z) - int_ddy_dx_dy_fct1(x2,z1, y3, x, y, z)
            val12 = int_ddy_dx_dy_fct1(x3,z2, y3, x, y, z) - int_ddy_dx_dy_fct1(x2,z2, y3, x, y, z)
            
            dat%x = x
            dat%y = y
            dat%z = z
            dat%thetas = theta0
            dat%rs = R
            dat%zs = z1
            f_ptr => int_ddy_dx_dy_fct2
            call qags_x ( f_ptr, dat, x2, x3, dat%epsabs, dat%epsrel, val21, dat%abserr_x, dat%neval_x, dat%ier_x )
        
            dat%zs = z2
            call qags_x ( f_ptr, dat, x2, x3, dat%epsabs, dat%epsrel, val22, dat%abserr_x, dat%neval_x, dat%ier_x )
            
            
            val1 = sign(1.,sin(theta0))*sign(1.,cos(theta0))/(4*pi) * ( val11 - val21 )
            
            val2 = sign(1.,sin(theta0))*sign(1.,cos(theta0))/(4*pi) * ( val12 - val22 )
         
           
        end subroutine
        
        function int_ddy_dx_dy_fct1( xp, zp, y3, x, y, z )
        real,intent(in) :: xp,zp, y3, x, y, z
        real :: int_ddy_dx_dy_fct1,arg
        
            arg = xp-x + P(x,y,z,xp,y3,zp)
            
            if ( arg .le. 10*tiny(1.) ) then
                arg = 10*tiny(1.)
            endif
            
            int_ddy_dx_dy_fct1 = log( arg )
            
            
        end function
        
        function int_ddy_dx_dy_fct2( xp, dat )              
        real,intent(in) :: xp
        class(dataCollectionBase), target :: dat
        real :: x,y,z,theta0,R,zp
        real :: int_ddy_dx_dy_fct2
            x = dat%x
            y = dat%y
            z = dat%z
            theta0 = dat%thetas
            R = dat%rs
            zp = dat%zs
            int_ddy_dx_dy_fct2 = 1./ (P(x,y,z,xp,sign(1.,sin(theta0))*sqrt(R**2-xp**2),zp))
        end function 
        
         subroutine int_ddz_dx_dy(dat, val1, val2 )
         real,intent(inout) :: val1,val2
         class(dataCollectionBase), intent(inout), target :: dat
         real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
         real :: x1,x2,x3,y1,y2,y3
         real :: dtheta,theta0
         real :: val11,val22,val12,val21
         procedure (f_int_dat), pointer :: f_ptr => null ()
            
         call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
         call getCorners( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
           
            
            val11 = int_ddz_dx_dy_fct1( y2, z1, x3, x, y, z ) - int_ddz_dx_dy_fct1( y1, z1, x3, x, y, z )
            val12 = int_ddz_dx_dy_fct1( y2, z2, x3, x, y, z ) - int_ddz_dx_dy_fct1( y1, z2, x3, x, y, z )
                        
            
            f_ptr => int_ddz_dx_dy_fct2
            
            dat%x = x
            dat%y = y
            dat%z = z
            dat%thetas = theta0
            dat%rs = R
            dat%zs = z1
            
            call qags_x ( f_ptr, dat, y1, y2, dat%epsabs, dat%epsrel, val21, dat%abserr_x, dat%neval_x, dat%ier_x )
            
            dat%zs = z2
            
            call qags_x ( f_ptr, dat, y1, y2, dat%epsabs, dat%epsrel, val22, dat%abserr_x, dat%neval_x, dat%ier_x )
            
            !sign(cos(theta0)) is caused by the order of integration; in
            !the first and fourth quadrants we integrate from x3 to
            !sqrt(R^2-y^2) while in the second and third quadrants we
            !integrate from sqrt(R^2-y^2) to x3
            val1 = -1*sign(1.,sin(theta0))*sign(1.,cos(theta0))/(4*pi) * ( val11 + val21 )  
            val2 = -1*sign(1.,sin(theta0))*sign(1.,cos(theta0))/(4*pi) * ( val12 + val22 )  
            
         end subroutine
         
         subroutine int_ddz_dx_dy_inv(dat, val1, val2 )
         real,intent(inout) :: val1,val2
         class(dataCollectionBase), intent(inout), target :: dat
         real :: x,y,z,R,theta1,theta2,z1,z2,xrot,phi
         real :: x1,x2,x3,y1,y2,y3
         real :: dtheta,theta0
         real :: val11,val22,val12,val21
         procedure (f_int_dat), pointer :: f_ptr => null ()
            
         call getParameters( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
     
         call getCorners_inv( R, theta1, theta2, theta0, dtheta, x1, x2, x3, y1, y2, y3 )
           
            
            val11 = int_ddz_dx_dy_fct1( y3, z1, x3, x, y, z ) - int_ddz_dx_dy_fct1( y1, z1, x3, x, y, z )
            val12 = int_ddz_dx_dy_fct1( y3, z2, x3, x, y, z ) - int_ddz_dx_dy_fct1( y1, z2, x3, x, y, z )
                        
            
            f_ptr => int_ddz_dx_dy_fct2
            
            dat%x = x
            dat%y = y
            dat%z = z
            dat%thetas = theta0
            dat%rs = R
            dat%zs = z1
            
            call qags_x ( f_ptr, dat, y1, y3, dat%epsabs, dat%epsrel, val21, dat%abserr_x, dat%neval_x, dat%ier_x )
            
            dat%zs = z2
            
            call qags_x ( f_ptr, dat, y1, y3, dat%epsabs, dat%epsrel, val22, dat%abserr_x, dat%neval_x, dat%ier_x )
            
            !sign(cos(theta0)) is caused by the order of integration; in
            !the first and fourth quadrants we integrate from x3 to
            !sqrt(R^2-y^2) while in the second and third quadrants we
            !integrate from sqrt(R^2-y^2) to x3
            !change the sign as the x=x3 constant line is to the right of the curved line
            val1 = 1*sign(1.,sin(theta0))*sign(1.,cos(theta0))/(4*pi) * ( val11 + val21 )  
            val2 = 1*sign(1.,sin(theta0))*sign(1.,cos(theta0))/(4*pi) * ( val12 + val22 )  
            
         end subroutine
        
         function int_ddz_dx_dy_fct1( yp, zp, x3, x, y, z )
         real,intent(in) :: yp, zp, x3, x, y, z
         real :: int_ddz_dx_dy_fct1
         
         !::Cover the limit where z-zp goes to zero
         if ( abs ( z-zp ) .lt. 10*tiny(1.) ) then
             int_ddz_dx_dy_fct1 = atan(  huge(1.) )
         else
             
            int_ddz_dx_dy_fct1 =  atan( ((x-x3)*(y-yp) ) / ((z-zp)*P(x,y,z,x3,yp,zp))   )
         endif
                  
         end function
         
         function int_ddz_dx_dy_fct2( yp, dat )         
         real,intent(in) :: yp
         class(dataCollectionBase), target :: dat
         real :: int_ddz_dx_dy_fct2
         real :: zp, x, y, z, theta0,R
         
         x = dat%x
         y = dat%y
         z = dat%z
         theta0 = dat%thetas
         R = dat%rs
         zp = dat%zs
         
         int_ddz_dx_dy_fct2 =  (x-sign(1.,cos(theta0))*sqrt(R**2-yp**2))*(z-zp)/(P(x,y,z,sign(1.,cos(theta0))*sqrt(R**2-yp**2),yp,zp)*((y-yp)**2+(z-zp)**2))
         
         end function
        
        
    subroutine getParameters_rot_trick( dat, x, y, z, R, theta1, theta2, z1, z2, xrot, phi )
    class(dataCollectionBase), intent(in), target :: dat
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
    class(dataCollectionBase), intent(in), target :: dat
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
    