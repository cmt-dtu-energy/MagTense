module TileNComponents
    use TileCylPieceTensor
    use TileRectanagularPrismTensor
    use TileCircPieceTensor
    implicit none
    
    !::General base-type for alle the different tile types
    type MagTile
        !::Specific for a cylindrical tile piece
        real :: r0, theta0, z0, dr, dtheta, dz
        
        !::Specific for a rectangular prism
        real :: a, b, c
        
        !::Generel variables, shared among all tile types
        real,dimension(3) :: M
        real,dimension(3) :: u_ea,u_oa1,u_oa2    
        real :: mu_r_ea,mu_r_oa,Mrem
        integer :: tileType        !::defines whether the tile is cylindrical, a prism, an ellipsoid and so on
        real,dimension(3) :: offset !::the centre coordinates relative to the global coordinate system
        real,dimension(3) :: rotAngles !:: rotation angles (phi_x, phi_y, phi_z) about the principle axes of the tile with respect to the centre of the tile
        integer :: magnetType !::defines whether the tile is a hard or soft magnet
        integer :: stateFunctionIndex !::index matching an entry into an array of type MagStatStateFunction. Used by soft ferromagnets (when interpolation on an M vs H curve is necessary)
        
        !::internal variables that are used for averaging the internal field
        integer,dimension(3) :: n_ave
        real,dimension(:,:),allocatable :: H_ave_pts,H_ave
        real,dimension(:,:,:,:),allocatable :: N_ave_pts
        logical :: isIterating
        integer :: fieldEvaluation
    end type MagTile
    
    integer,parameter :: tileTypeCylPiece=1,tileTypePrism=2,tileTypeCircPiece=3,tileTypeEllipsoid=4
    integer,parameter :: magnetTypeHard=1,magnetTypeSoft=2
    integer,parameter :: fieldEvaluationCentre=1,fieldEvaluationAverage=2
    
    contains
    
    
    subroutine  getN_CylPiece( cylP, x, N )
    type(MagTile),intent(in) :: cylP
    real,intent(in) :: x
    class( dataCollectionBase ), pointer :: dat
    real :: theta1, theta2
    real :: int_dDdx_dr_dz1,int_dDdx_dr_dz2
    real :: int_cos_dDdx_dz_dtheta1, int_cos_dDdx_dz_dtheta2
    real :: int_sin_dDdx_dz_dtheta1, int_sin_dDdx_dz_dtheta2
    real :: int_r_dDdx_dr_dtheta1, int_r_dDdx_dr_dtheta2
    real :: int_sin_dDdy_dz_dtheta1, int_sin_dDdy_dz_dtheta2
    real :: int_r_cos_dDdy_dz_dtheta1,int_r_cos_dDdy_dz_dtheta2
    real :: int_r_dDdy_dr_dtheta1, int_r_dDdy_dr_dtheta2
    real :: int_cos_dDdz_dz_dtheta1,int_cos_dDdz_dz_dtheta2
    real :: int_dDdz_dr_dz1, int_dDdz_dr_dz2
    real :: int_dDdy_dr_dz1,int_dDdy_dr_dz2
    real :: int_r_dDdz_dr_dtheta1, int_r_dDdz_dr_dtheta2
    real :: int_sin_dDdz_dz_dtheta1,int_sin_dDdz_dz_dtheta2
    real,dimension(3,3), intent(inout) :: N
            
        N(:,:) = 0
        allocate(dat)

        !get the integration limits
        dat%r1 = cylP%r0 - cylP%dr/2
        dat%r2 = cylP%r0 + cylP%dr/2

        dat%theta1 = cylP%theta0 - cylP%dtheta/2
        dat%theta2 = cylP%theta0 + cylP%dtheta/2

        dat%z1 = cylP%z0 - cylP%dz/2
        dat%z2 = cylP%z0 + cylP%dz/2                
         
        dat%x = x
        
        call int_dDdx_dr_dz( dat, int_dDdx_dr_dz1, int_dDdx_dr_dz2 )
        
        call int_cos_dDdx_dz_dtheta( dat, int_cos_dDdx_dz_dtheta1, int_cos_dDdx_dz_dtheta2 )
        
        call int_sin_dDdx_dz_dtheta( dat, int_sin_dDdx_dz_dtheta1, int_sin_dDdx_dz_dtheta2 )
        
        call int_r_dDdx_dr_dtheta( dat, int_r_dDdx_dr_dtheta1, int_r_dDdx_dr_dtheta2 )
        
        call int_sin_dDdy_dz_dtheta( dat, int_sin_dDdy_dz_dtheta1, int_sin_dDdy_dz_dtheta2 )
        
        call int_r_cos_dDdy_dz_dtheta( dat, int_r_cos_dDdy_dz_dtheta1, int_r_cos_dDdy_dz_dtheta2 )
        
        call int_r_dDdy_dr_dtheta( dat, int_r_dDdy_dr_dtheta1, int_r_dDdy_dr_dtheta2 )
        
        call int_cos_dDdz_dz_dtheta( dat, int_cos_dDdz_dz_dtheta1, int_cos_dDdz_dz_dtheta2 )
        
        call int_dDdz_dr_dz( dat, int_dDdz_dr_dz1, int_dDdz_dr_dz2)
        
        call int_dDdy_dr_dz( dat, int_dDdy_dr_dz1, int_dDdy_dr_dz2 )
        
        call int_r_dDdz_dr_dtheta( dat, int_r_dDdz_dr_dtheta1, int_r_dDdz_dr_dtheta2 )
        
       call int_sin_dDdz_dz_dtheta( dat, int_sin_dDdz_dz_dtheta1, int_sin_dDdz_dz_dtheta2 )
        
        
        theta1 = dat%theta1
        theta2 = dat%theta2
        
        
        N(1,1) = int_sin_dDdy_dz_dtheta2 - int_sin_dDdy_dz_dtheta1 &
            + cos(theta2) * int_dDdy_dr_dz2 - cos(theta1) * int_dDdy_dr_dz1 &
            + int_r_dDdz_dr_dtheta2 - int_r_dDdz_dr_dtheta1

        N(1,2) = int_r_cos_dDdy_dz_dtheta1 - int_r_cos_dDdy_dz_dtheta2 &
            + sin(theta2) * int_dDdy_dr_dz2 - sin(theta1) * int_dDdy_dr_dz1
        
        !N(1,3) = int_cos_dDdz_dz_dtheta1 - int_cos_dDdz_dz_dtheta2 &
        !    + sin(theta2) * int_dDdz_dr_dz2 - sin(theta1) * int_dDdz_dr_dz1
        
        N(2,1) = int_sin_dDdx_dz_dtheta1 - int_sin_dDdx_dz_dtheta2 &
            + cos(theta1) * int_dDdx_dr_dz1 - cos(theta2) * int_dDdx_dr_dz2
        
        
        N(2,2) = int_cos_dDdx_dz_dtheta2 - int_cos_dDdx_dz_dtheta1 &
                + sin(theta1) * int_dDdx_dr_dz1 - sin(theta2) * int_dDdx_dr_dz2 &
                + int_r_dDdz_dr_dtheta2 - int_r_dDdz_dr_dtheta1
        
        N(3,2) = int_r_dDdy_dr_dtheta1 - int_r_dDdy_dr_dtheta2
        
        !N(2,3) = int_sin_dDdz_dz_dtheta1 - int_sin_dDdz_dz_dtheta2 &
        !    + cos(theta1) * int_dDdz_dr_dz1 - cos(theta2) * int_dDdz_dr_dz2
        N(2,3) = N(3,2)
        N(3,1) = int_r_dDdx_dr_dtheta1 - int_r_dDdx_dr_dtheta2
       
        N(1,3) = N(3,1)        
        
               
        
        
        !::2018-01-08. Kaspar makes this hack to ensure consistency - there is an unknown problem with this component; it should be symmetric with 3,2 and appears not to be
        !::It is therefore forced to be symmetric for the time being until the problem has been figured out.
        !::2018-01-09. Kaspar figured out the probrlem in one of the tensor-functions (there was a division with r which should not be there). It is now fixed
        !N(2,3) = N(3,2)
        
        N(3,3) = int_sin_dDdy_dz_dtheta2 + int_cos_dDdx_dz_dtheta2 &
            - int_sin_dDdy_dz_dtheta1 - int_cos_dDdx_dz_dtheta1 &
            + cos(theta2) * int_dDdy_dr_dz2 - sin(theta2) * int_dDdx_dr_dz2 &
            -(cos(theta1) * int_dDdy_dr_dz1 - sin(theta1) * int_dDdx_dr_dz1)
    deallocate(dat)
    end subroutine getN_CylPiece

    subroutine getN_circPiece( tile, pos, Nout )
    type(MagTile),intent(in) :: tile
    class( dataCollectionBase ), pointer :: dat
    real,dimension(3),intent(in) :: pos
    real,dimension(3,3),intent(inout) :: Nout
    real :: int_ddx_dy_dz_val,int_ddx_dx_dz_val,int_ddx_cos_dtheta_dz_val,int_ddx_sin_dtheta_dz_val
    real :: int_ddx_dx_dy_val1,int_ddx_dx_dy_val2,int_ddy_dy_dz_val, int_ddy_dx_dz_val,int_ddy_cos_dtheta_dz_val
    real :: int_ddy_sin_dtheta_dz_val, int_ddy_dx_dy_val1,int_ddy_dx_dy_val2,int_ddz_dy_dz_val
    real :: int_ddz_dx_dz_val,int_ddz_cos_dtheta_dz_val,int_ddz_sin_dtheta_dz_val,int_ddz_dx_dy_val1,int_ddz_dx_dy_val2
    real :: phi,Nxx,Nyy,Nzz,Nxy,Nyx,Nxz,Nzx,Nyz,Nzy
    real :: R,x,y
            Nout(:,:) = 0.            

            !get the integration limits            
             allocate(dat)

            !get the integration limits
            dat%r1 = tile%r0 - tile%dr/2
            dat%r2 = tile%r0 + tile%dr/2

            dat%theta1 = tile%theta0 - tile%dtheta/2
            dat%theta2 = tile%theta0 + tile%dtheta/2
            
            dat%z1 = tile%z0 - tile%dz/2
            dat%z2 = tile%z0 + tile%dz/2    
            
            dat%x = pos(1)
            dat%y = pos(2)
            dat%z = pos(3)
            
            x = pos(1)
            y = pos(2)
            R = tile%r0 + tile%dr/2
           
            !x-component of the field
            !yz-plane (Mx)
            call int_ddx_dy_dz( dat, int_ddx_dy_dz_val )
                        
            !xz-plane (My)
            call int_ddx_dx_dz( dat, int_ddx_dx_dz_val )
                        
            !theta-z-plane (Mr=cos(theta0)*Mx + sin(theta0)*My)
            call int_ddx_cos_dtheta_dz( dat, int_ddx_cos_dtheta_dz_val )            
            
            call int_ddx_sin_dtheta_dz( dat, int_ddx_sin_dtheta_dz_val )            
            
            !xy-plane (Mz) (UP,positive)
            call int_ddx_dx_dy( dat, int_ddx_dx_dy_val1, int_ddx_dx_dy_val2 )            
            
            
            !y-component of the field
            !yz-plane (Mx)
            call int_ddy_dy_dz( dat, int_ddy_dy_dz_val )
                        
            !xz-plane (My)
            call int_ddy_dx_dz( dat, int_ddy_dx_dz_val )
                        
            !theta-z-plane (Mr=cos(theta0)*Mx + sin(theta0)*My)
            call int_ddy_cos_dtheta_dz( dat, int_ddy_cos_dtheta_dz_val )
            call int_ddy_sin_dtheta_dz( dat, int_ddy_sin_dtheta_dz_val )
            
            
            !xy-plane (Mz) (UP,positive, DOWN negative)
            call int_ddy_dx_dy( dat, int_ddy_dx_dy_val1, int_ddy_dx_dy_val2 )
                        
            !z-component of the field
            !yz-plane (Mx)
            call int_ddz_dy_dz( dat, int_ddz_dy_dz_val )            
            
            !xz-plane (My)
            call int_ddz_dx_dz( dat, int_ddz_dx_dz_val )            
            
            !theta-z-plane (Mr=cos(theta0)*Mx + sin(theta0)*My)
            call int_ddz_cos_dtheta_dz(dat, int_ddz_cos_dtheta_dz_val )
            call int_ddz_sin_dtheta_dz(dat, int_ddz_sin_dtheta_dz_val )
            
            
            !xy-plane (Mz) (UP,positive)
            call int_ddz_dx_dy ( dat, int_ddz_dx_dy_val1, int_ddz_dx_dy_val2 )
                        
            !define the rotation-trick-variable
            phi = atan2(y,x)
            

            
            Nxx = int_ddx_dy_dz_val + R/(4*pi) * &
                ( cos(phi) * ( cos(phi) * int_ddx_cos_dtheta_dz_val - sin(phi) * int_ddx_sin_dtheta_dz_val ) &
                + sin(phi) * ( cos(phi) * int_ddy_cos_dtheta_dz_val - sin(phi) * int_ddy_sin_dtheta_dz_val ) )
            
            !Kaspar: Something is not correct with int_ddx_sin_dtheta_dz. 
            !Kaspar: Fixed!
            Nxy = int_ddx_dx_dz_val + R/(4*pi) * &
                ( cos(phi) * ( sin(phi) * int_ddx_cos_dtheta_dz_val + cos(phi) * int_ddx_sin_dtheta_dz_val ) &
                + sin(phi) * ( sin(phi) * int_ddy_cos_dtheta_dz_val + cos(phi) * int_ddy_sin_dtheta_dz_val ) )
            
            
            Nxz = int_ddx_dx_dy_val2 - int_ddx_dx_dy_val1
                                    
            Nyx = int_ddy_dy_dz_val + R/(4*pi) * &
                ( sin(phi) * ( cos(phi) * int_ddx_cos_dtheta_dz_val - sin(phi) * int_ddx_sin_dtheta_dz_val ) &
                - cos(phi) * ( cos(phi) * int_ddy_cos_dtheta_dz_val - sin(phi) * int_ddy_sin_dtheta_dz_val ) )
            
            Nyy = int_ddy_dx_dz_val + R/(4*pi) * &
                ( sin(phi) * ( sin(phi) * int_ddx_cos_dtheta_dz_val + cos(phi) * int_ddx_sin_dtheta_dz_val ) &
                - cos(phi) * ( sin(phi) * int_ddy_cos_dtheta_dz_val + cos(phi) * int_ddy_sin_dtheta_dz_val ) )
            
            Nyz = int_ddy_dx_dy_val2 - int_ddy_dx_dy_val1
            
            Nzx = int_ddz_dy_dz_val - R/(4*pi) * ( cos(phi) * int_ddz_cos_dtheta_dz_val - sin(phi) * int_ddz_sin_dtheta_dz_val )
            
            Nzy = int_ddz_dx_dz_val - R/(4*pi) * ( sin(phi) * int_ddz_cos_dtheta_dz_val + cos(phi) * int_ddz_sin_dtheta_dz_val )
            
            Nzz = int_ddz_dx_dy_val2 - int_ddz_dx_dy_val1
            
            Nout(1,1) = Nxx;
            Nout(1,2) = Nxy;
            Nout(1,3) = Nxz;
            Nout(2,1) = Nyx;
            Nout(2,2) = Nyy;
            Nout(2,3) = Nyz;
            Nout(3,1) = Nzx;
            Nout(3,2) = Nzy;
            Nout(3,3) = Nzz;

        deallocate(dat)
    
    end subroutine getN_circPiece
    
    
    !::Calculates N from the analytical expression in 3D
    !::Given the prism tile (prism) and the position vector to it (pos = (x,y,z) )
    !::Returns a (3,3) array N_out
    subroutine getN_prism_3D( prism, pos, N_out )
    type(MagTile),intent(in) :: prism
    real,intent(in),dimension(3) :: pos
    real,intent(out),dimension(3,3) :: N_out
    
    real :: a,b,c,x,y,z
    real :: nom,denom,nom_l,nom_h,denom_l,denom_h    
    real,parameter :: lim_scl_h=1.01,lim_scl_l=0.98
    real :: xl,xh,yl,yh,zl,zh,al,ah,bl,bh,cl,ch,lim
    
    a = prism%a
    b = prism%b
    c = prism%c
    
    x = pos(1)
    y = pos(2)
    z = pos(3)
    
    !::Diagonal elements
    N_out(1,1) = 1./(4.*pi) * ( atan(f_3D(a,b,c,x,y,z))   + atan(f_3D(a,b,c,-x,y,z))  + atan(f_3D(a,b,c,x,-y,z)) + &
                                atan(f_3D(a,b,c,x,y,-z))  + atan(f_3D(a,b,c,-x,-y,z)) + atan(f_3D(a,b,c,x,-y,-z)) + &
                                atan(f_3D(a,b,c,-x,y,-z)) + atan(f_3D(a,b,c,-x,-y,-z)) )



    N_out(2,2) = 1./(4.*pi) * ( atan(g_3D(a,b,c,x,y,z))   + atan(g_3D(a,b,c,-x,y,z))  + atan(g_3D(a,b,c,x,-y,z)) + &
                                atan(g_3D(a,b,c,x,y,-z))  + atan(g_3D(a,b,c,-x,-y,z)) + atan(g_3D(a,b,c,x,-y,-z)) + &
                                atan(g_3D(a,b,c,-x,y,-z)) + atan(g_3D(a,b,c,-x,-y,-z)) )
                            
                            

    N_out(3,3) = 1./(4.*pi) * ( atan(h_3D(a,b,c,x,y,z))   + atan(h_3D(a,b,c,-x,y,z))  + atan(h_3D(a,b,c,x,-y,z)) + &
                                atan(h_3D(a,b,c,x,y,-z))  + atan(h_3D(a,b,c,-x,-y,z)) + atan(h_3D(a,b,c,x,-y,-z)) + &
                                atan(h_3D(a,b,c,-x,y,-z)) + atan(h_3D(a,b,c,-x,-y,-z)) )                            
                            

    !::Off-diagonal elements
    nom = FF_3D(a,b,c,x,y,z)  * FF_3D(-a,-b,c,x,y,z) * FF_3D(a,-b,-c,x,y,z) * FF_3D(-a,b,-c,x,y,z)
    denom = FF_3D(a,-b,c,x,y,z) * FF_3D(-a,b,c,x,y,z)  * FF_3D(a,b,-c,x,y,z)  * FF_3D(-a,-b,-c,x,y,z)

    if ( denom .eq. 0 .or. nom .eq. 0 ) then
        !Find the limit
        lim = getF_limit(a,b,c,x,y,z,FF_3D)
        N_out(1,2) = -1./(4.*pi) * log( lim )        
    else
        N_out(1,2) = -1./(4.*pi) * log( nom / denom )
    endif


    !::the tensor is symmetric
    N_out(2,1) = N_out(1,2)

    nom = GG_3D(a,b,c,x,y,z)  * GG_3D(-a,-b,c,x,y,z) * GG_3D(a,-b,-c,x,y,z) * GG_3D(-a,b,-c,x,y,z)
    denom = GG_3D(a,-b,c,x,y,z) * GG_3D(-a,b,c,x,y,z)  * GG_3D(a,b,-c,x,y,z)  * GG_3D(-a,-b,-c,x,y,z)

    if ( denom .eq. 0 .or. nom .eq. 0 ) then
    !Find the limit
      
        lim = getF_limit(a,b,c,x,y,z,GG_3D)       
    
        N_out(2,3) = -1./(4.*pi) * log( lim )
    else
        N_out(2,3) = -1./(4.*pi) * log( nom/denom )
    endif


    N_out(3,2) = N_out(2,3)

    nom = HH_3D(a,b,c,x,y,z)  * HH_3D(-a,-b,c,x,y,z) * HH_3D(a,-b,-c,x,y,z) * HH_3D(-a,b,-c,x,y,z)
    denom = HH_3D(a,-b,c,x,y,z) * HH_3D(-a,b,c,x,y,z)  * HH_3D(a,b,-c,x,y,z)  * HH_3D(-a,-b,-c,x,y,z)

    if ( denom .eq. 0 .or. nom .eq. 0 ) then
    
        lim = getF_limit(a,b,c,x,y,z,HH_3D)           
    
        N_out(1,3) = -1./(4.*pi) * log( lim )
    else
        N_out(1,3) = -1./(4.*pi) * log( nom / denom )
    endif


    N_out(3,1) = N_out(1,3)
    
    end subroutine
        
end module TileNComponents
    
    