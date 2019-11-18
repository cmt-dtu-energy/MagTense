module TileNComponents
    use TileCylPieceTensor
    use TileRectanagularPrismTensor
    use TileCircPieceTensor
    use TilePlanarCoilTensor
    use TileTriangle
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
        real,dimension(3) :: color !! color rgb triplet
        integer :: magnetType !::defines whether the tile is a hard or soft magnet
        integer :: stateFunctionIndex !::index matching an entry into an array of type MagStateFunction. Used by soft ferromagnets (when interpolation on an M vs H curve is necessary)
        integer :: includeInIteration,exploitSymmetry
        real,dimension(3) :: symmetryOps !! 1 for symmetry, -1 for anti-symmetry ((1) for about xy plane, (2) for about (xz) plane and (3) for about yz plane)
        real :: Mrel !! the current relative change of the magnetization (i.e. abs((M1-M2)/M2 ) where M1 is the latest magnetization norm and M2 is the previous one
        real,dimension(3,4) :: vert !Vertices, used e.g. by a tetrahedron
        
        !::internal variables that are used for averaging the internal field
        integer,dimension(3) :: n_ave
        real,dimension(:,:),allocatable :: H_ave_pts,H_ave
        real,dimension(:,:,:,:),allocatable :: N_ave_pts
        logical :: isIterating
        integer :: fieldEvaluation
        
        
    end type MagTile
    
    
    !> Interface for the call to a specific implementation of the demagnetization tensor
    !! @param prismTile is the tile from which the tensor is wanted
    !! @param diffPos is a (3,1) vector with the relative position of the point of interest (relative to the tile's offset)
    !! @param N is a (3,3) matrix where the tensor is returned to
    abstract interface
      
      subroutine N_tensor_subroutine ( prismTile, diffPos, N )
         import MagTile
         type(MagTile),intent(in) :: prismTile
         real,dimension(3),intent(in) :: diffPos
         real,dimension(3,3),intent(inout) :: N
      end subroutine N_tensor_subroutine
    end interface
    
    integer,parameter :: tileTypeCylPiece=1,tileTypePrism=2,tileTypeCircPiece=3,tileTypeCircPieceInverted=4,tileTypeTetrahedron=5,tileTypeEllipsoid=10,tileTypePlanarCoil=101
    integer,parameter :: magnetTypeHard=1,magnetTypeSoft=2
    integer,parameter :: fieldEvaluationCentre=1,fieldEvaluationAverage=2
    
    contains
    
    subroutine setupEvaluationPoints( cylTile )
    type(MagTile),intent(inout) :: cylTile
    real,dimension(:),allocatable :: r, theta, z
    real :: dr,dtheta,dz
    integer :: n
    integer :: cnt,j,k,l
    
        cyltile%fieldevaluation  = fieldevaluationaverage
        cyltile%n_ave(1) = 1
        cyltile%n_ave(2) = 1
        cyltile%n_ave(3) = 1
        n = cyltile%n_ave(1)*cyltile%n_ave(2)*cyltile%n_ave(3)
        allocate(cyltile%H_ave_pts(n,3),cyltile%H_ave(n,3),r(n),theta(n),z(n))                
        dr = cyltile%dr / cyltile%n_ave(1)
        dtheta = cyltile%dtheta / cyltile%n_ave(2)
        dz = cyltile%dz / cyltile%n_ave(3)
                
        !::set as default value so other users of the code don't have to update.
        cylTile%isIterating = .false.
                
        cnt = 1
        do j=1,cyltile%n_ave(1)
            do k=1,cyltile%n_ave(2)
                do l=1,cyltile%n_ave(3)
                    r(cnt) = dr/2 + (j-1)*dr + cylTile%r0 - cylTile%dr / 2
                    theta(cnt) = dtheta/2 + (k-1)*dtheta + cylTile%theta0 - cylTile%dtheta / 2
                    z(cnt) = dz/2 + (l-1)*dz + cylTile%z0 - cylTile%dz / 2
                    cnt = cnt + 1
                enddo
            enddo
        enddo
        cyltile%h_ave_pts(:,1) = r * cos( theta ) + cyltile%offset(1)
        cyltile%h_ave_pts(:,2) = r * sin( theta ) + cyltile%offset(2)
        cyltile%h_ave_pts(:,3) = z + cyltile%offset(3)
        cylTile%H_ave(:,:) = 0.
        deallocate(r,theta,z)            
        
    end subroutine setupEvaluationPoints
    
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
    real :: x_
            
        N(:,:) = 0
        allocate(dat)
        
        
        x_ = x
        
        if ( abs(x) .le. 1e-15 .AND. x .ge. 0. ) then
            x_ = 1e-15
        elseif ( abs(x) .le. 1e-15 .AND. x .lt. 0. ) then
            x_ = -1e-15
        endif
        
        

        !get the integration limits
        dat%r1 = cylP%r0 - cylP%dr/2
        dat%r2 = cylP%r0 + cylP%dr/2

        dat%theta1 = cylP%theta0 - cylP%dtheta/2
        dat%theta2 = cylP%theta0 + cylP%dtheta/2

        dat%z1 = cylP%z0 - cylP%dz/2
        dat%z2 = cylP%z0 + cylP%dz/2                
         
        dat%x = x_
        
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
        
        N = N * 1./(4.*pi)
        
    deallocate(dat)
    end subroutine getN_CylPiece

    !::Get the demag tensor for a piece of a circle defined with the angular extent dtheta, radius R and angle through the middle theta0.
    !::The piece is then delimited by to orthogonal cords that are parallel to a primary axis each, respectively.
    !::This object is thus defined by three points in the xy-plane: (x1,y1), (x2,y2) and (x3,y3) given by
    !::x1 = R * cos( theta - dtheta/2 ), y1 = R * sin( theta - theta/2 )
    !::x2 = R * cos( theta + dtheta/2 ), y2 = R * sin( theta + dtheta/2 )
    !::x3 = R * cos( theta + dtheta/2 ), y3 = R * sin( theta - dtheta/2 )
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
    real :: R,x,y,theta0,dtheta
    real,dimension(3,3) :: RefMatr,SymmMatr
    real,dimension(3) :: pos_
    integer :: i
            Nout(:,:) = 0.            

            !::Hack to circumvent problem with x=0, y=0 and z=0
            pos_ = pos
            do i=1,3
                if ( abs(pos(i)) .le. 1e-15 .AND. pos(i) .ge. 0. ) then
                    pos_(i) = 1e-15
                elseif ( abs(pos(i)) .le. 1e-15 .AND. pos(i) .lt. 0. ) then
                    pos_(i) = -1e-15
                endif
            enddo
            
            
            
            !get the integration limits            
             allocate(dat)

            !get the integration limits
            dat%r1 = tile%r0 - tile%dr/2
            dat%r2 = tile%r0 + tile%dr/2

            theta0 = tile%theta0
            dtheta = tile%dtheta
                        
            
            dat%z1 = tile%z0 - tile%dz/2
            dat%z2 = tile%z0 + tile%dz/2    
            
            dat%x = pos_(1)
            dat%y = pos_(2)
            dat%z = pos_(3)
            
            !::Rotate to 1st quadrant and mirror the x- and y-coordinates accordingly
            dat%x = dat%x * sign(1.,cos(theta0))
            dat%y = dat%y * sign(1.,sin(theta0))
            
            x = dat%x
            y = dat%y
            R = tile%r0 + tile%dr/2

            
            
            !::Reflection matrix for correctly mirroring the N-tensor
            !::Note that the reflection is on the field and not the tensor components
            SymmMatr(:,:) = 1.
            RefMatr(:,:) = 0.
            RefMatr(1,1) = 1.
            RefMatr(2,2) = 1.
            RefMatr(3,3) = 1.
            
            !define the rotation-trick-variable            
            phi = atan2_custom( y, x )
            
            if ( sign(1.,cos(tile%theta0)) .lt. 0. .AND. sign(1.,sin(tile%theta0)) .ge. 0. ) then
                !Second quadrant                
                theta0 = pi - theta0
                !::Reflection about the y-axis (i.e. Mx and Hx changes sign)
                SymmMatr(:,1) = -1. 
                RefMatr(1,1) = -1.
            elseif ( sign(1.,cos(tile%theta0)) .lt. 0. .AND. sign(1.,sin(tile%theta0)) .lt. 0. ) then
                !Third quadrant
                theta0 = theta0 - pi
                !::Reflection about the y-axis (Mx and Hx change sign)
                SymmMatr(:,1) = -1.
                RefMatr(2,2) = -1.
                !::reflection about the x-axis (My and Hy change sign)
                SymmMatr(:,2) = -1.                                
                RefMatr(1,1) = -1.
            elseif ( sign(1.,cos(tile%theta0)) .ge. 0. .AND. sign(1.,sin(tile%theta0)) .lt. 0. ) then
                !Fourth quadrant
                theta0 = 2*pi - theta0
                
                !::Reflection about the x-axis (My and Hy change sign)
                SymmMatr(:,2) = -1.                
                RefMatr(2,2) = -1.
            
            endif
                            
            
            dat%theta1 = theta0 - tile%dtheta/2
            dat%theta2 = theta0 + tile%dtheta/2
            
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
            
            Nout(1,1) = Nxx
            Nout(1,2) = Nxy
            Nout(1,3) = Nxz
            Nout(2,1) = Nyx
            Nout(2,2) = Nyy
            Nout(2,3) = Nyz
            Nout(3,1) = Nzx
            Nout(3,2) = Nzy
            Nout(3,3) = Nzz
            
            !::Mirror the solution (change sign of  Mx if mirror in y-axis, sign of My if mirror in x-axis)
            Nout = matmul( RefMatr, SymmMatr * Nout )

        deallocate(dat)
    
    end subroutine getN_circPiece
    
    
    !::Get the demag tensor for a piece of a circle defined with the angular extent dtheta, radius R and angle through the middle theta0.
    !::The piece is then delimited by to mutually orthogonal lines that are parallel to a primary axis each, respectively and extent on the external part of the circle
    !::This object is thus defined by three points in the xy-plane: (x1,y1), (x2,y2) and (x3,y3) given by
    !::x1 = R * cos( theta - dtheta/2 ), y1 = R * sin( theta - theta/2 )
    !::x2 = R * cos( theta + dtheta/2 ), y2 = R * sin( theta + dtheta/2 )
    !::x3 = R * cos( theta - dtheta/2 ), y3 = R * sin( theta + dtheta/2 )
    !::It then extends in the z-direction by dz centered at z0
    subroutine getN_circPiece_Inv( tile, pos, Nout )
    type(MagTile),intent(in) :: tile
    class( dataCollectionBase ), pointer :: dat
    real,dimension(3),intent(in) :: pos
    real,dimension(3,3),intent(inout) :: Nout
    real :: int_ddx_dy_dz_val,int_ddx_dx_dz_val,int_ddx_cos_dtheta_dz_val,int_ddx_sin_dtheta_dz_val
    real :: int_ddx_dx_dy_val1,int_ddx_dx_dy_val2,int_ddy_dy_dz_val, int_ddy_dx_dz_val,int_ddy_cos_dtheta_dz_val
    real :: int_ddy_sin_dtheta_dz_val, int_ddy_dx_dy_val1,int_ddy_dx_dy_val2,int_ddz_dy_dz_val
    real :: int_ddz_dx_dz_val,int_ddz_cos_dtheta_dz_val,int_ddz_sin_dtheta_dz_val,int_ddz_dx_dy_val1,int_ddz_dx_dy_val2
    real :: phi,Nxx,Nyy,Nzz,Nxy,Nyx,Nxz,Nzx,Nyz,Nzy
    real :: R,x,y,theta0,dtheta
    real,dimension(3,3) :: RefMatr,SymmMatr
    real,dimension(3) :: pos_
    integer :: i
            Nout(:,:) = 0.            

            !::Hack to circumvent problem with x=0, y=0 and z=0
            !::15-6-2018: Note from Kaspar: Hacks like this WILL come back and bite you in the ass
            !::Changed 1e-10 to 1e-20 - otherwise the orders when trying to do anything serious get messed up
            pos_ = pos
            do i=1,3
                if ( abs(pos(i)) .le. 1e-15 .AND. pos(i) .ge. 0. ) then
                    pos_(i) = 1e-15
                elseif ( abs(pos(i)) .le. 1e-15 .AND. pos(i) .lt. 0. ) then
                    pos_(i) = -1e-15
                endif
            enddo
            
            
            
            !get the integration limits            
             allocate(dat)

            !get the integration limits
            dat%r1 = tile%r0 - tile%dr/2
            dat%r2 = tile%r0 + tile%dr/2

            theta0 = tile%theta0
            dtheta = tile%dtheta
                        
            
            dat%z1 = tile%z0 - tile%dz/2
            dat%z2 = tile%z0 + tile%dz/2    
            
            dat%x = pos_(1)
            dat%y = pos_(2)
            dat%z = pos_(3)
            
            !::Rotate to 1st quadrant and mirror the x- and y-coordinates accordingly
            dat%x = dat%x * sign(1.,cos(theta0))
            dat%y = dat%y * sign(1.,sin(theta0))
            
            x = dat%x
            y = dat%y
            R = tile%r0 + tile%dr/2

            
            
            !::Reflection matrix for correctly mirroring the N-tensor
            !::Note that the reflection is on the field and not the tensor components
            SymmMatr(:,:) = 1.
            RefMatr(:,:) = 0.
            RefMatr(1,1) = 1.
            RefMatr(2,2) = 1.
            RefMatr(3,3) = 1.
            
            !define the rotation-trick-variable            
            phi = atan2_custom( y, x )
            
            if ( sign(1.,cos(tile%theta0)) .lt. 0. .AND. sign(1.,sin(tile%theta0)) .ge. 0. ) then
                !Second quadrant                
                theta0 = pi - theta0
                !::Reflection about the y-axis (i.e. Mx and Hx changes sign)
                SymmMatr(:,1) = -1. 
                RefMatr(1,1) = -1.
            elseif ( sign(1.,cos(tile%theta0)) .lt. 0. .AND. sign(1.,sin(tile%theta0)) .lt. 0. ) then
                !Third quadrant
                theta0 = theta0 - pi
                !::Reflection about the y-axis (Mx and Hx change sign)
                SymmMatr(:,1) = -1.
                RefMatr(2,2) = -1.
                !::reflection about the x-axis (My and Hy change sign)
                SymmMatr(:,2) = -1.                                
                RefMatr(1,1) = -1.
            elseif ( sign(1.,cos(tile%theta0)) .ge. 0. .AND. sign(1.,sin(tile%theta0)) .lt. 0. ) then
                !Fourth quadrant
                theta0 = 2*pi - theta0
                
                !::Reflection about the x-axis (My and Hy change sign)
                SymmMatr(:,2) = -1.                
                RefMatr(2,2) = -1.
            
            endif
                            
            
            dat%theta1 = theta0 - tile%dtheta/2
            dat%theta2 = theta0 + tile%dtheta/2
            
            !x-component of the field
            !yz-plane (Mx)
            call int_ddx_dy_dz_inv( dat, int_ddx_dy_dz_val )
                        
            !xz-plane (My)
            call int_ddx_dx_dz_inv( dat, int_ddx_dx_dz_val )
                        
            !theta-z-plane (Mr=cos(theta0)*Mx + sin(theta0)*My)
            call int_ddx_cos_dtheta_dz( dat, int_ddx_cos_dtheta_dz_val )            
            !::Change the sign since the inverted circ piece has the opposite normal vectors
            int_ddx_cos_dtheta_dz_val = -int_ddx_cos_dtheta_dz_val
            
            call int_ddx_sin_dtheta_dz( dat, int_ddx_sin_dtheta_dz_val )            
            int_ddx_sin_dtheta_dz_val = -int_ddx_sin_dtheta_dz_val
            
            !xy-plane (Mz) (UP,positive)
            call int_ddx_dx_dy_inv( dat, int_ddx_dx_dy_val1, int_ddx_dx_dy_val2 )            
            
            
            !y-component of the field
            !yz-plane (Mx)
            call int_ddy_dy_dz_inv( dat, int_ddy_dy_dz_val )
                        
            !xz-plane (My)
            call int_ddy_dx_dz_inv( dat, int_ddy_dx_dz_val )
                        
            !theta-z-plane (Mr=cos(theta0)*Mx + sin(theta0)*My)
            call int_ddy_cos_dtheta_dz( dat, int_ddy_cos_dtheta_dz_val )
            int_ddy_cos_dtheta_dz_val = -int_ddy_cos_dtheta_dz_val
            call int_ddy_sin_dtheta_dz( dat, int_ddy_sin_dtheta_dz_val )
            int_ddy_sin_dtheta_dz_val = -int_ddy_sin_dtheta_dz_val
            
            
            !xy-plane (Mz) (UP,positive, DOWN negative)
            call int_ddy_dx_dy_inv( dat, int_ddy_dx_dy_val1, int_ddy_dx_dy_val2 )
                        
            !z-component of the field
            !yz-plane (Mx)
            call int_ddz_dy_dz_inv( dat, int_ddz_dy_dz_val )            
            
            !xz-plane (My)
            call int_ddz_dx_dz_inv( dat, int_ddz_dx_dz_val )            
            
            !theta-z-plane (Mr=cos(theta0)*Mx + sin(theta0)*My)
            call int_ddz_cos_dtheta_dz(dat, int_ddz_cos_dtheta_dz_val )
            int_ddz_cos_dtheta_dz_val = -int_ddz_cos_dtheta_dz_val
            call int_ddz_sin_dtheta_dz(dat, int_ddz_sin_dtheta_dz_val )
            int_ddz_sin_dtheta_dz_val = -int_ddz_sin_dtheta_dz_val
            
            
            !xy-plane (Mz) (UP,positive)
            call int_ddz_dx_dy_inv ( dat, int_ddz_dx_dy_val1, int_ddz_dx_dy_val2 )
                        
            
            
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
            
            Nout(1,1) = Nxx
            Nout(1,2) = Nxy
            Nout(1,3) = Nxz
            Nout(2,1) = Nyx
            Nout(2,2) = Nyy
            Nout(2,3) = Nyz
            Nout(3,1) = Nzx
            Nout(3,2) = Nzy
            Nout(3,3) = Nzz
            
            !::Mirror the solution (change sign of  Mx if mirror in y-axis, sign of My if mirror in x-axis)
            Nout = matmul( RefMatr, SymmMatr * Nout )

        deallocate(dat)
    
    end subroutine getN_circPiece_Inv
    
    !::Calculates N from the analytical expression in 3D
    !::Given the prism tile (prism) and the position vector to it (pos = (x,y,z) )
    !::Returns a (3,3) array N_out
    subroutine getN_prism_3D( prism, pos, N_out )
    type(MagTile),intent(in) :: prism
    real,intent(in),dimension(3) :: pos
    real,intent(inout),dimension(3,3) :: N_out
    
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
    
    !!Change the sign so that the output tensor follows the same definition as all the other tensors
    N_out = -1.* N_out
    
    end subroutine
        
    
    subroutine getN_tensor_tetrahedron ( tile, r, N )         
         type(MagTile),intent(in) :: tile
         real,dimension(3),intent(in) :: r
         real,dimension(3,3),intent(inout) :: N
         
         real,dimension(3,3) :: N_tmp,P
         integer :: i
         
         N(:,:) = 0.
         
         !Loop through all the four triangular faces
         do i=1,4
            call getN_Triangle( cshift( tile%vert, (i-1), 2 ), r, N_tmp, P )
            N = N + N_tmp
         enddo
         
         
      end subroutine getN_tensor_tetrahedron
    
    
    subroutine getN_dummy( k )
    integer :: k
    
    end subroutine getN_dummy
    
    !::
    !::Returns the demag tensor for a planar coil with inner radius Rin, outer radius Rout
    !::and NL windings between these two radii. Each winding is assumed to carry the same current, I, which
    !::is considered to be equivalent to the magnetization of the "tile"
    !::It is assumed tha NL = 100 yields a good approximation to a continuous coil. Then the field strength is governed by the current, I
    subroutine getN_PlanarCoil( tile, pos, Nout )
    type(MagTile),intent(in) :: tile
    class( dataCollectionBase ), pointer :: dat
    real,dimension(3),intent(in) :: pos
    real,dimension(3,3),intent(inout) :: Nout
    
    real,dimension(2) :: Npol,tmp
    real :: r,z,theta,Rin,Rout,R_loop
    integer :: i
    integer,parameter :: NL=100
    
    Npol(:) = 0.
    tmp(:) = 0.
    Nout(:,:) = 0.
    
    allocate(dat)
    
    
    
    !::Get the loop dimensions
    Rin = tile%a
    Rout = tile%b
    
    !::Convert to polar coordinates
    r = sqrt( pos(1)**2 + pos(2)**2 )
    z = pos(3)
    theta = atan2( pos(2), pos(1) )
    
    do i=1,NL
        R_loop = Rin + (Rout - Rin) * (i-1.)/(NL-1.)
        call getTensorFromSingleLoop( dat, R_loop, r , z, tmp )
        Npol = Npol + tmp
    enddo
    
    !::Convert from polar to Cartesian coordinates
    !::The resulting tensor is diagonal as there are no cross-terms
    !::due to the symmetry of the coil
    Nout(1,1) = cos( theta ) * Npol(1)
    Nout(2,2) = sin( theta ) * Npol(1)
    Nout(3,3) = Npol(2)
    
    deallocate(dat)
    end subroutine getN_PlanarCoil
    
    
    
end module TileNComponents
    
    