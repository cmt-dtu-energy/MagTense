module TileNComponents
    use TileCylPieceTensor
    use TileRectanagularPrismTensor
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
    end type MagTile
    
    integer,parameter :: tileTypeCylPiece=1,tileTypePrism=2,tileTypeEllipsoid=3
    integer,parameter :: magnetTypeHard=1,magnetTypeSoft=2
    
    
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
    
    