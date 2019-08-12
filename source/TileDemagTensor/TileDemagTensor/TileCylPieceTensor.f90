module TileCylPieceTensor
    use QuadPack
    use TileTensorHelperFunctions
    use SpecialFunctions
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
        !k^2 corresponds to m

    !::integrates over r (analytically) and z (numerically)
    !::val and val2 correspond to different integrations on the surfaces with theta1 and theta2, respectively
    subroutine int_dDdx_dr_dz( dat, val1, val2 )        
    real,intent(inout) :: val1, val2
    procedure (f_int_dat), pointer :: f_ptr => null ()
    class(dataCollectionBase), intent(inout), target :: dat
                
    val1 = 0
    val2 = 0
    f_ptr => int_dDdx_dr_dz_fct

    dat%thetas = dat%theta1
        
    call qags_x ( f_ptr, dat, dat%z1, dat%z2, dat%epsabs, dat%epsrel, val1, dat%abserr_x, dat%neval_x, dat%ier_x )
        
    dat%thetas = dat%theta2
        
    call qags_x ( f_ptr, dat, dat%z1, dat%z2, dat%epsabs, dat%epsrel, val2, dat%abserr_x, dat%neval_x, dat%ier_x )

    end subroutine int_dDdx_dr_dz
        
    function int_dDdx_dr_dz_fct( z,  dat )
    real,intent(in) :: z
    class(dataCollectionBase), intent(inout), target :: dat     
    real :: int_dDdx_dr_dz_fct,r1,r2,thetas,x
        
    r1 = dat%r1
    r2 = dat%r2
    thetas = dat%thetas
    x = dat%x
        
            int_dDdx_dr_dz_fct = ( r2 * x * ( 1 - cos(thetas)**2 ) + z**2 *cos(thetas) ) / ( L(x, thetas,z) * M(r2,x,thetas,z) ) &
        - ( r1 * x * ( 1 - cos(thetas)**2 ) + z**2 *cos(thetas) ) / ( L(x, thetas,z) * M(r1,x,thetas,z) )
              
    end function

    subroutine int_cos_dDdx_dz_dtheta( dat, val1, val2 )
    class(dataCollectionBase), intent(inout), target :: dat
    real,intent(inout) :: val1, val2
    real :: theta1,theta2,dtheta,r1,r2,z1,z2,x
    
    call getParameters( dat, r1, r2, theta1, theta2, z1, z2, x )
    
        !bring theta1 into the range 0 to 2pi and translate theta2
        !accordingly
                        
        if ( theta1 .lt. 0 .and. theta2 .lt. 0 ) then
            dtheta = theta2-theta1
            theta1 = theta1 - 2 * pi * floor(theta1/(2*pi))
            theta2 = theta1 + dtheta
        elseif ( theta2 .gt. 2*pi ) then
            dtheta = theta2 - theta1
            theta2 = theta2-2*pi*floor(theta2/(2*pi))
            theta1 = theta2 - dtheta
        endif
            
            
        if ( theta1 .lt. 0 .AND. theta2 .gt. 0 ) then
            !first integrate from theta1 to zero
            val1 = int_cos_dDdx_dz_dtheta_fct( 0., z2, r1, x ) - int_cos_dDdx_dz_dtheta_fct( theta1, z2, r1, x ) - ( int_cos_dDdx_dz_dtheta_fct( 0., z1, r1, x ) - int_cos_dDdx_dz_dtheta_fct( theta1, z1, r1, x ) )
            !then integrate from zero to theta2
            val1 = val1 + ( int_cos_dDdx_dz_dtheta_fct( 0., z2, r1, x ) - int_cos_dDdx_dz_dtheta_fct( theta2, z2, r1, x ) ) - ( int_cos_dDdx_dz_dtheta_fct( 0., z1, r1, x ) - int_cos_dDdx_dz_dtheta_fct( theta2, z1, r1, x ) )
            !val1 = -val1
        else
            val1 = -(int_cos_dDdx_dz_dtheta_fct( theta2, z2, r1, x ) - int_cos_dDdx_dz_dtheta_fct( theta1, z2, r1, x ) - ( int_cos_dDdx_dz_dtheta_fct( theta2, z1, r1, x ) - int_cos_dDdx_dz_dtheta_fct( theta1, z1, r1, x ) ) )            
            
        endif
        
        if ( theta1 .lt. 0 .AND. theta2 .gt. 0 ) then
            !first integrate from theta1 to zero
            val2 = int_cos_dDdx_dz_dtheta_fct( 0., z2, r2, x ) - int_cos_dDdx_dz_dtheta_fct( theta1, z2, r2, x ) - ( int_cos_dDdx_dz_dtheta_fct( 0., z1, r2, x ) - int_cos_dDdx_dz_dtheta_fct( theta1, z1, r2, x ) )
            !then integrate from zero to theta2
            val2 = val2 + ( int_cos_dDdx_dz_dtheta_fct( 0., z2, r2, x ) - int_cos_dDdx_dz_dtheta_fct( theta2, z2, r2, x ) ) - ( int_cos_dDdx_dz_dtheta_fct( 0., z1, r2, x ) - int_cos_dDdx_dz_dtheta_fct( theta2, z1, r2, x ) )             
            !val2 = -val2
        else
            val2 = -(int_cos_dDdx_dz_dtheta_fct( theta2, z2, r2, x ) - int_cos_dDdx_dz_dtheta_fct( theta1, z2, r2, x ) - ( int_cos_dDdx_dz_dtheta_fct( theta2, z1, r2, x ) - int_cos_dDdx_dz_dtheta_fct( theta1, z1, r2, x ) ) )            
        endif
                
    end subroutine
    function int_cos_dDdx_dz_dtheta_fct( theta, z, rs, x )        
    real :: int_cos_dDdx_dz_dtheta_fct
    real,intent(in) :: theta, z, rs, x
    real :: elF,elE,elPi
    
    elf = ellf( C_no_sign(theta), K(rs,x,z) )
    
    elE = elle( C_no_sign(theta), K(rs,x,z) )
    
    elPi = ellpi( C_no_sign(theta), B(rs,x), K(rs,x,z) )
    
    int_cos_dDdx_dz_dtheta_fct = z / ( 2*x**2 * ( rs + x ) * sqrt((rs+x)**2+z**2 ) )&
                * ((rs-x)*(rs**2+x**2) * elPi + (rs+x)*(((rs+x)**2+z**2)*elE - (2*rs**2 + z**2)*elF ) )
    return
    end function

        !!Returns the definite integral of rs * int(int( sin(theta) * dDdx * dz *
        !!dtheta ) ) from theta1 to theta2 and z1 to z2
        !!Verified numerically
    subroutine int_sin_dDdx_dz_dtheta( dat, val1, val2 )
    class(dataCollectionBase), intent(inout), target :: dat
    real,intent(inout) :: val1,val2
    real :: theta1,theta2,z1,z2,r1,r2,x
    
    call getParameters( dat, r1, r2, theta1, theta2, z1, z2, x )
    
    val1 = int_sin_dDdx_dz_dtheta_fct( theta2, z2, r1, x ) - int_sin_dDdx_dz_dtheta_fct( theta1, z2, r1, x ) - ( int_sin_dDdx_dz_dtheta_fct( theta2, z1, r1, x ) - int_sin_dDdx_dz_dtheta_fct( theta1, z1, r1, x ) )
        
    val2 = int_sin_dDdx_dz_dtheta_fct( theta2, z2, r2, x ) - int_sin_dDdx_dz_dtheta_fct( theta1, z2, r2, x ) - ( int_sin_dDdx_dz_dtheta_fct( theta2, z1, r2, x ) - int_sin_dDdx_dz_dtheta_fct( theta1, z1, r2, x ) )
        
    end subroutine
    function int_sin_dDdx_dz_dtheta_fct( theta, z, rs, x )
    real :: int_sin_dDdx_dz_dtheta_fct
    real,intent(in) :: theta,z,rs,x
        int_sin_dDdx_dz_dtheta_fct = 1./(4.*x**2) * ( ( rs**2 - x**2 ) * ( log( M(rs,x,theta,z) -z) - log( M(rs,x,theta,z) + z ) ) - 2 * z * M(rs,x,theta,z) )
    end function

    subroutine int_r_dDdx_dr_dtheta( dat, val1, val2 )
    class(dataCollectionBase), intent(inout), target :: dat
    real,intent(inout) :: val1,val2
    procedure (f_int_dat), pointer :: f_ptr => null ()               
                
    val1 = 0
    val2 = 0
    f_ptr => int_r_dDdx_dr_dtheta_fct

    dat%zs = dat%z1
        
    call qags_x ( f_ptr, dat, dat%theta1, dat%theta2, dat%epsabs, dat%epsrel, val1, dat%abserr_x, dat%neval_x, dat%ier_x )
        
    dat%zs = dat%z2
        
    call qags_x ( f_ptr, dat, dat%theta1, dat%theta2, dat%epsabs, dat%epsrel, val2, dat%abserr_x, dat%neval_x, dat%ier_x )

        
    end subroutine
    
    function int_r_dDdx_dr_dtheta_fct( theta, dat )
    real :: int_r_dDdx_dr_dtheta_fct
    class(dataCollectionBase), intent(inout), target :: dat
    real,intent(in) :: theta
    real :: r1,r2,x,zs
        
    r1 = dat%r1
    r2 = dat%r2
    zs = dat%zs
    x = dat%x
        
    int_r_dDdx_dr_dtheta_fct =  int_r_dDdx_dr_dtheta_fct1( r2, theta, zs, x ) - int_r_dDdx_dr_dtheta_fct1(r1, theta, zs, x) &
                                + ( int_r_dDdx_dr_dtheta_fct2( r2, theta, zs, x ) - int_r_dDdx_dr_dtheta_fct2( r1, theta, zs, x )) &
                                + ( int_r_dDdx_dr_dtheta_fct3( r2, theta, zs, x) - int_r_dDdx_dr_dtheta_fct3(r1, theta, zs, x) )
    end function
    
    function int_r_dDdx_dr_dtheta_fct1( r, theta, zs, x )
    real :: int_r_dDdx_dr_dtheta_fct1
    real,intent(in) :: r, theta,zs,x
        int_r_dDdx_dr_dtheta_fct1 =  cos(theta) * log( r - x*cos(theta) + M(r,x,theta,zs) )
    end function
    
    function int_r_dDdx_dr_dtheta_fct2( r, theta, zs, x )
    real :: int_r_dDdx_dr_dtheta_fct2
    real,intent(in) :: r, theta,zs,x
        int_r_dDdx_dr_dtheta_fct2 =  -cos(theta) * ( 2*x**2*r*cos(theta)**2 - cos(theta)*(x**3 + x*zs**2) -r*(x**2+zs**2) ) / (M(r,x,theta,zs)*L(x,theta,zs) )
    end function
    
    function int_r_dDdx_dr_dtheta_fct3( r, theta, zs, x )
    real :: int_r_dDdx_dr_dtheta_fct3
    real,intent(in) :: r, theta,zs,x
        int_r_dDdx_dr_dtheta_fct3 = x * ( r * x * cos(theta) - x**2 - zs**2 ) / ( M(r,x,theta,zs) * L(x,theta,zs) )
    end function
        
    !!Returns the definite integral of rs * int( int( sin(theta) * dDdy * dz *
    !!dtheta ) ) evaluated at rs and integrated from theta1 to theta2 and z1 to z2 and with the point of interest, x.
    !Verified (so far) with the unit testing
    subroutine int_sin_dDdy_dz_dtheta( dat, val1, val2 )
    class(dataCollectionBase), intent(inout), target :: dat
    real,intent(inout) :: val1,val2
    real :: theta1,theta2,z1,z2,r1,r2,x,dtheta
    
    call getParameters( dat, r1, r2, theta1, theta2, z1, z2, x )
        !bring theta1 into the range 0 to 2pi and translate theta2
        !accordingly
            
        if ( theta1 .lt. 0 .AND. theta2 .lt. 0 ) then
            dtheta = theta2-theta1
            theta1 = theta1 - 2 * pi * floor(theta1/(2*pi))
            theta2 = theta1 + dtheta
        elseif ( theta2 .gt. 2*pi ) then
            dtheta = theta2 - theta1
            theta2 = theta2 - 2*pi*floor(theta2/(2*pi))
            theta1 = theta2 - dtheta
        endif
            
    
        if ( theta1 .lt. 0 .AND. theta2 .gt. 0 ) then         
            !first integrate from theta1 to zero
            val1 = int_sin_dDdy_dz_dtheta_fct( 0., z2, r1, x ) - int_sin_dDdy_dz_dtheta_fct( theta1, z2, r1, x ) - ( int_sin_dDdy_dz_dtheta_fct( 0., z1, r1, x ) - int_sin_dDdy_dz_dtheta_fct( theta1, z1, r1, x ) )
            !then integrate from zero to theta2
            val1 = val1 + ( int_sin_dDdy_dz_dtheta_fct( 0., z2, r1, x ) - int_sin_dDdy_dz_dtheta_fct( theta2, z2, r1, x ) ) - ( int_sin_dDdy_dz_dtheta_fct( 0., z1, r1, x ) - int_sin_dDdy_dz_dtheta_fct( theta2, z1, r1, x ) ) 
            
            !val1 = -val1
        else
            val1 = -(int_sin_dDdy_dz_dtheta_fct( theta2, z2, r1, x ) - int_sin_dDdy_dz_dtheta_fct( theta1, z2, r1, x )  - ( int_sin_dDdy_dz_dtheta_fct( theta2, z1, r1, x ) - int_sin_dDdy_dz_dtheta_fct( theta1, z1, r1, x ) ) )
        endif
        
        if ( theta1 .lt. 0. .AND. theta2 .gt. 0. ) then        
            !first integrate from theta1 to zero
            val2 = int_sin_dDdy_dz_dtheta_fct( 0., z2, r2, x ) - int_sin_dDdy_dz_dtheta_fct( theta1, z2, r2, x ) - ( int_sin_dDdy_dz_dtheta_fct( 0., z1, r2, x ) - int_sin_dDdy_dz_dtheta_fct( theta1, z1, r2, x ) )
            !then integrate from zero to theta2
            val2 = val2 + ( int_sin_dDdy_dz_dtheta_fct( 0., z2, r2, x ) - int_sin_dDdy_dz_dtheta_fct( theta2, z2, r2, x ) ) - ( int_sin_dDdy_dz_dtheta_fct( 0., z1, r2, x ) - int_sin_dDdy_dz_dtheta_fct( theta2, z1, r2, x ) ) 
            
            !val2 = -val2
        else
            val2 = -(int_sin_dDdy_dz_dtheta_fct( theta2, z2, r2, x ) - int_sin_dDdy_dz_dtheta_fct( theta1, z2, r2, x )  - ( int_sin_dDdy_dz_dtheta_fct( theta2, z1, r2, x ) - int_sin_dDdy_dz_dtheta_fct( theta1, z1, r2, x ) ) )
        endif
        
    end subroutine
    
    function int_sin_dDdy_dz_dtheta_fct( theta, z, rs, x )
    real :: int_sin_dDdy_dz_dtheta_fct
    real :: theta,z,rs,x
    real :: elF,elE,elPi
    
    elf = ellf( C_no_sign(theta), K(rs,x,z) )
    elE = elle( C_no_sign(theta), K(rs,x,z) )
    elPi = ellpi( C_no_sign(theta), B(rs,x), K(rs,x,z) )
    
    int_sin_dDdy_dz_dtheta_fct = -1./2.*z  / ( x**2*sqrt((rs+x)**2+z**2) ) * ( ( rs - x )**2 * elPi + ( ( rs + x )**2 + z**2 ) * elE -2 * ( rs**2 + x**2 +1./2.*z**2 ) * elF )

    return
    end function
    
        !!Returns the definite integral of rs * int(int(cos(theta)*dDdy *
        !!dz*dtheta) from z1 to z2 and theta1 to theta2 with respect to the point
        !!of interest, x
        subroutine int_r_cos_dDdy_dz_dtheta( dat, val1, val2 )
            class(dataCollectionBase), intent(inout), target :: dat
            real,intent(inout) :: val1,val2
            real :: theta1,theta2,z1,z2,r1,r2,x
    
            call getParameters( dat, r1, r2, theta1, theta2, z1, z2, x )
            
            val1 = int_r_cos_dDdy_dz_dtheta_fct( theta2, z2, r1, x ) - int_r_cos_dDdy_dz_dtheta_fct( theta1, z2, r1, x ) - ( int_r_cos_dDdy_dz_dtheta_fct( theta2, z1, r1, x ) - int_r_cos_dDdy_dz_dtheta_fct( theta1, z1, r1, x ) )
            
            val2 = int_r_cos_dDdy_dz_dtheta_fct( theta2, z2, r2, x ) - int_r_cos_dDdy_dz_dtheta_fct( theta1, z2, r2, x ) - ( int_r_cos_dDdy_dz_dtheta_fct( theta2, z1, r2, x ) - int_r_cos_dDdy_dz_dtheta_fct( theta1, z1, r2, x ) )            
        end subroutine 
        
        function int_r_cos_dDdy_dz_dtheta_fct( theta, z, rs, x )
            real :: int_r_cos_dDdy_dz_dtheta_fct
            real,intent(in) :: theta,z,rs,x
            
            int_r_cos_dDdy_dz_dtheta_fct = 1./(4.*x**2) * ( ( rs**2 + x**2 ) * ( log( M(rs,x,theta,z) - z ) - log( M(rs,x,theta,z) + z ) ) - 2 * z * M(rs, x, theta, z) )
        end function
        
        !subroutine int_r_dDdy_dr_dtheta( dat, val1, val2 )
        !    class(dataCollectionBase), intent(inout), target :: dat
        !    real,intent(inout) :: val1,val2
        !    real :: theta1,theta2,z1,z2,r1,r2,x
        !
        !    call getParameters( dat, r1, r2, theta1, theta2, z1, z2, x )
        !    
        !    val1 = int_r_dDdy_dr_dtheta_fct(r2,theta2,z1,x) - int_r_dDdy_dr_dtheta_fct(r1, theta2,z1,x) - ( int_r_dDdy_dr_dtheta_fct(r2,theta1,z1,x) - int_r_dDdy_dr_dtheta_fct(r1,theta1,z1,x) )
        !    
        !    val2 = int_r_dDdy_dr_dtheta_fct(r2,theta2,z2,x) - int_r_dDdy_dr_dtheta_fct(r1, theta2,z2,x) - ( int_r_dDdy_dr_dtheta_fct(r2,theta1,z2,x) - int_r_dDdy_dr_dtheta_fct(r1,theta1,z2,x) )
        !end
        !
        !function int_r_dDdy_dr_dtheta_fct( r, theta, zs, x )
        !real :: int_r_dDdy_dr_dtheta_fct
        !real,intent(in) :: r,theta,zs,x
        !            int_r_dDdy_dr_dtheta_fct =  1./(2.*r*x*sqrt(x**2+zs**2)) * ( 2*r*(x**2+zs**2)*(atanh(G_p(r,x,theta,zs)) - atanh(G_m(r,x,theta,zs))) &
        !        + sqrt(x**2+zs**2)*(-(r**2+x**2+zs**2)*log(r*(r-x*cos(theta)+M(r,x,theta,zs))) &
        !        + (r**2-2*r*x*cos(theta)+x**2+zs**2)*log(r-x*cos(theta)+M(r,x,theta,zs)) &
        !        - 2*r*M(r,x,theta,zs) - (r**2+x**2+zs**2)*log(2.) + 2*r*x*cos(theta) -(r**2+x**2+zs**2) ))
        !    return
        !end function
        
    subroutine int_r_dDdy_dr_dtheta( dat, val1, val2 )
    class(dataCollectionBase), intent(inout), target :: dat
    real,intent(inout) :: val1,val2
    procedure (f_int_dat), pointer :: f_ptr => null ()                    
                
    val1 = 0
    val2 = 0
    f_ptr => int_r_dDdy_dr_dtheta_fct

    dat%zs = dat%z1
        
    call qags_x ( f_ptr, dat, dat%theta1, dat%theta2, dat%epsabs, dat%epsrel, val1, dat%abserr_x, dat%neval_x, dat%ier_x )
        
    dat%zs = dat%z2
        
    call qags_x ( f_ptr, dat, dat%theta1, dat%theta2, dat%epsabs, dat%epsrel, val2, dat%abserr_x, dat%neval_x, dat%ier_x )

        
    end subroutine
    
    function int_r_dDdy_dr_dtheta_fct( theta, dat )
    real :: int_r_dDdy_dr_dtheta_fct
    class(dataCollectionBase), intent(inout), target :: dat
    real,intent(in) :: theta
    real :: r1,r2,x,zs
        
    r1 = dat%r1
    r2 = dat%r2
    zs = dat%zs
    x = dat%x
        
    int_r_dDdy_dr_dtheta_fct =  int_r_dDdy_dr_dtheta_fct1( r2, theta, zs, x ) - int_r_dDdy_dr_dtheta_fct1(r1, theta, zs, x) &
                                + ( int_r_dDdy_dr_dtheta_fct2( r2, theta, zs, x ) - int_r_dDdy_dr_dtheta_fct2( r1, theta, zs, x )) 
    end function
    
    function int_r_dDdy_dr_dtheta_fct1( r, theta, zs, x )
    real :: int_r_dDdy_dr_dtheta_fct1
    real,intent(in) :: r, theta,zs,x
        int_r_dDdy_dr_dtheta_fct1 =  sin(theta) * log( r - x*cos(theta) + M(r,x,theta,zs) )
    end function
    
    function int_r_dDdy_dr_dtheta_fct2( r, theta, zs, x )
    real :: int_r_dDdy_dr_dtheta_fct2
    real,intent(in) :: r, theta,zs,x
        int_r_dDdy_dr_dtheta_fct2 =  -sin(theta) * ( 2*x**2*r*cos(theta)**2 - cos(theta)*(x**3 + x*zs**2) -r*(x**2+zs**2) ) / (M(r,x,theta,zs)*L(x,theta,zs) )
    end function

        !!Returns the definite integral of rs * int( int ( cos(theta) * dDdz *
        !!dz*dtheta ) ) from z1 to z2 and theta1 to theta2
        subroutine int_cos_dDdz_dz_dtheta( dat, val1, val2 )
        class(dataCollectionBase), intent(inout), target :: dat
        real,intent(inout) :: val1,val2
        real :: theta1,theta2,z1,z2,r1,r2,x
        procedure (f_int_dat), pointer :: f_ptr => null ()                    
        call getParameters( dat, r1, r2, theta1, theta2, z1, z2, x )
                                         
        val1 = 0
        val2 = 0
        f_ptr => int_cos_dDdz_dz_dtheta_fct_int

        dat%rs = dat%r1
        
        call qags_x ( f_ptr, dat, dat%theta1, dat%theta2, dat%epsabs, dat%epsrel, val1, dat%abserr_x, dat%neval_x, dat%ier_x )
        
        dat%rs = dat%r2
        
        call qags_x ( f_ptr, dat, dat%theta1, dat%theta2, dat%epsabs, dat%epsrel, val2, dat%abserr_x, dat%neval_x, dat%ier_x )

        end subroutine
        
        function int_cos_dDdz_dz_dtheta_fct_int( theta, dat )
            real :: int_cos_dDdz_dz_dtheta_fct_int
            class(dataCollectionBase), intent(inout), target :: dat
            real, intent(in) :: theta
            real :: z1,z2,rs,x
            
            z1 = dat%z1
            z2 = dat%z2
            
            rs = dat%rs
            x = dat%x
            
            int_cos_dDdz_dz_dtheta_fct_int = int_cos_dDdz_dz_dtheta_fct(theta,z2,rs,x) - int_cos_dDdz_dz_dtheta_fct(theta,z1,rs,x)
        end function

        function int_cos_dDdz_dz_dtheta_fct( theta, z, rs, x )
        real :: int_cos_dDdz_dz_dtheta_fct
        real,intent(in) :: theta,z, rs, x
            
            int_cos_dDdz_dz_dtheta_fct =-rs * cos(theta) / sqrt( rs**2 + x**2 - 2*rs*x*cos(theta) + z**2 )
        end function
        
        
        
        subroutine int_dDdz_dr_dz( dat, val1, val2)
        class(dataCollectionBase), intent(inout), target :: dat
        real,intent(inout) :: val1,val2
        real :: theta1,theta2,z1,z2,r1,r2,x
        
            call getParameters( dat, r1, r2, theta1, theta2, z1, z2, x )
            
            val1 = int_dDdz_dr_dz_fct( r2, z2, theta1, x ) - int_dDdz_dr_dz_fct( r1, z2, theta1, x ) - ( int_dDdz_dr_dz_fct( r2, z1, theta1, x ) - int_dDdz_dr_dz_fct( r1, z1, theta1, x) )
            
            val2 = int_dDdz_dr_dz_fct( r2, z2, theta2, x ) - int_dDdz_dr_dz_fct( r1, z2, theta2, x ) - ( int_dDdz_dr_dz_fct( r2, z1, theta2, x ) - int_dDdz_dr_dz_fct( r1, z1, theta2, x) )
        end subroutine
        
        function int_dDdz_dr_dz_fct( r, z, thetas, x )
        real :: int_dDdz_dr_dz_fct
        real,intent(in) :: r, z, thetas, x
        
        int_dDdz_dr_dz_fct = -log( r - x*cos(thetas) + M(r,x,thetas,z) )
        end function


        !!Returns the definite integral of cos(thetas) * int( int( dDdy * dr * dz )
        !!) integrated from r1 to r2 and z1 to z2 at the point of interest, x.
        !!I could not so far find a closed-form analytical expression for the
        !!second part of the integral so it is done numerically, i.e. the integral
        !!over r has been done analytically while the integral over z is done
        !!numerically.
        subroutine int_dDdy_dr_dz( dat, val1, val2 )
        real,intent(inout) :: val1,val2
        procedure (f_int_dat), pointer :: f_ptr => null ()
        class(dataCollectionBase), intent(inout), target :: dat
                
            val1 = 0
            val2 = 0
            f_ptr => int_dDdy_dr_dz_fct_int

            dat%thetas = dat%theta1
        
            call qags_x ( f_ptr, dat, dat%z1, dat%z2, dat%epsabs, dat%epsrel, val1, dat%abserr_x, dat%neval_x, dat%ier_x )
        
            dat%thetas = dat%theta2
        
            call qags_x ( f_ptr, dat, dat%z1, dat%z2, dat%epsabs, dat%epsrel, val2, dat%abserr_x, dat%neval_x, dat%ier_x )
            
            
        end subroutine
        
        function int_dDdy_dr_dz_fct_int( z, dat )
        real :: int_dDdy_dr_dz_fct_int
        real,intent(in) :: z
        class(dataCollectionBase), intent(inout), target :: dat
        real :: thetas,r1,r2,x
        
            r1 = dat%r1
            r2 = dat%r2
            x = dat%x
            thetas = dat%thetas
                    
            int_dDdy_dr_dz_fct_int =  int_dDdy_dr_dz_fct( r2, z, thetas, x ) - int_dDdy_dr_dz_fct( r1, z, thetas, x )
            return
        end function
        function int_dDdy_dr_dz_fct( r, z, thetas, x )
        real :: int_dDdy_dr_dz_fct
        real,intent(in) :: r,z,thetas,x
            !::necessary to enforce the integrand to be zero when theta is small since there is a discontinuity (zero over zero problem)
            !::
            if ( abs(sin(thetas)) .lt. 1e-14 ) then
                int_dDdy_dr_dz_fct = 0.0
           else            
                int_dDdy_dr_dz_fct = -sin(thetas) * ( r*x*cos(thetas) - x**2 - z**2 ) / ( L( x, thetas, z) * M( r, x, thetas, z ) )
           endif           
           return
        end function
         
        !!Returns the definite integral of int( int( dDdz * r * dr * dtheta ) )
        !!from r1 to r2 and theta1 to theta2 at constant z(=zs) and with respect to
        !!the point of interest, x
        !!This has to be done numerically (over theta given the analytical
        !!expression for the integral over r)
        subroutine int_r_dDdz_dr_dtheta( dat, val1, val2 )
        real,intent(inout) :: val1,val2
        procedure (f_int_dat), pointer :: f_ptr => null ()
        class(dataCollectionBase), intent(inout), target :: dat
                
            val1 = 0
            val2 = 0
            f_ptr => int_r_dDdz_dr_dtheta_fct_int

            dat%zs = dat%z1
        
            call qags_x ( f_ptr, dat, dat%theta1, dat%theta2, dat%epsabs, dat%epsrel, val1, dat%abserr_x, dat%neval_x, dat%ier_x )
        
            dat%zs = dat%z2
        
            call qags_x ( f_ptr, dat, dat%theta1, dat%theta2, dat%epsabs, dat%epsrel, val2, dat%abserr_x, dat%neval_x, dat%ier_x )
            
            
        end subroutine
        
        function int_r_dDdz_dr_dtheta_fct_int( theta, dat )
        real :: int_r_dDdz_dr_dtheta_fct_int
        real,intent(in) :: theta
        class(dataCollectionBase), intent(inout), target :: dat
        real :: r1,r2,zs,x
        r1 = dat%r1
        r2 = dat%r2
        zs = dat%zs
        x = dat%x
        
        int_r_dDdz_dr_dtheta_fct_int =  int_r_dDdz_dr_dtheta_fct(r2,theta,zs,x) - int_r_dDdz_dr_dtheta_fct(r1,theta,zs,x)
        end function
        
        function int_r_dDdz_dr_dtheta_fct( r, theta, zs, x )
        real :: int_r_dDdz_dr_dtheta_fct
        real,intent(in) :: r,theta,zs,x
        int_r_dDdz_dr_dtheta_fct = -zs * (r * x * cos(theta) - x**2 - zs**2 ) / ( M(r, x, theta, zs) * L(x, theta, zs) )
        end function
        
        !::2018-01-09 Kaspar: Something is not correct with this function
        subroutine int_sin_dDdz_dz_dtheta( dat, val1, val2 )
        real,intent(inout) :: val1,val2    
        class(dataCollectionBase), intent(inout), target :: dat
        real :: theta1,theta2,z1,z2,r1,r2,x
        
        call getParameters( dat, r1, r2, theta1, theta2, z1, z2, x )
        

            val1 = int_sin_dDdz_dz_dtheta_fct( theta2, z2, r1, x ) - int_sin_dDdz_dz_dtheta_fct( theta1, z2, r1, x ) - ( int_sin_dDdz_dz_dtheta_fct(theta2, z1, r1, x) - int_sin_dDdz_dz_dtheta_fct(theta1, z1, r1, x) )
            
            val2 = int_sin_dDdz_dz_dtheta_fct( theta2, z2, r2, x ) - int_sin_dDdz_dz_dtheta_fct( theta1, z2, r2, x ) - ( int_sin_dDdz_dz_dtheta_fct(theta2, z1, r2, x) - int_sin_dDdz_dz_dtheta_fct(theta1, z1, r2, x) )
        end subroutine
        
        function int_sin_dDdz_dz_dtheta_fct( theta, z, rs, x )
        real :: int_sin_dDdz_dz_dtheta_fct
        real,intent(in) :: theta,z,rs,x
        int_sin_dDdz_dz_dtheta_fct = - sqrt(rs**2 - 2*rs*x*cos(theta) + x**2 + z**2) / x
        return
        end function
        
        subroutine getParameters( dat, r1, r2, theta1, theta2, z1, z2, x )
        class(dataCollectionBase), intent(in), target :: dat
        real,intent(inout) :: r1, r2,theta1,theta2,z1,z2,x
        r1 = dat%r1
        r2 = dat%r2
    
        theta1 = dat%theta1
        theta2 = dat%theta2
    
        z1 = dat%z1
        z2 = dat%z2
    
        x = dat%x
        end subroutine
    
end module TileCylPieceTensor
    
