
    module MagStat2GetSolution

    use TileNComponents
    use MagTileIO
    implicit none
    
    
    
    contains
    
    !::
    !::General function to call given a number of (different) tiles    
    !::Input arguments
    !::tiles: Array of type MagTile, size n_tiles
    !::H the output magnetic field, size [n_ele,3]
    !::pts the points at which the field should be evaluated, size [n_ele,3]
    !::integer n_tiles, the number of tiles
    !::integer n_ele, the number of points at which to evaluate the field
    subroutine getFieldFromTiles( tiles, H, pts, n_tiles, n_ele )
    type(MagTile),intent(inout),dimension(n_tiles) :: tiles
    real,dimension(n_ele,3),intent(inout) :: H
    real,dimension(n_ele,3),intent(in) :: pts
    integer,intent(in) :: n_tiles,n_ele
    
    integer :: i,prgCnt,tid,prog,OMP_GET_THREAD_NUM
    real,dimension(:,:),allocatable :: H_tmp
    real,dimension(:,:,:),allocatable :: N_out
    integer,parameter :: cbCnt = 10
    
    
    H(:,:) = 0.
    
    prgCnt = 0
    prog = 0
    
  
    !$OMP PARALLEL DO PRIVATE(i,H_tmp,N_out)
    do i=1,n_tiles
        
        !Make sure to allocate H_tmp on the heap and for each thread
        !$OMP CRITICAL
        allocate( H_tmp(n_ele,3),N_out(3,3,n_ele) )
        H_tmp(:,:) = 0.
        N_out(:,:,:) = 0.
        !$OMP END CRITICAL
        !::Here a selection of which subroutine to use should be done, i.e. whether the tile
        !:: is cylindrical, a prism or an ellipsoid
        select case (tiles(i)%tileType )
        case (tileTypeCylPiece)
            call getFieldFromCylTile( tiles(i), H_tmp, pts, n_ele, N_out )    
        case (tileTypePrism)
            call getFieldFromRectangularPrismTile( tiles(i), H_tmp, pts, n_ele, N_out )
        case (tileTypeEllipsoid)
        case default        
            
        end select
        
        
        !$OMP CRITICAL
        H = H + H_tmp
        !prgCnt = prgCnt + 1
        !tid = omp_get_thread_num()
        !if ( tid .eq. 0 ) then
        !    if ( floor( 100. * prgCnt / (1.*n_tiles) ) .ge. cbCnt ) then
        !        prgCnt = 0
        !        prog = prog + 10                
        !        call displayProgress( prog )
        !    endif
        !endif        
        deallocate(H_tmp,N_out)
       !$OMP END CRITICAL
       
        
        
        
        
    enddo
    !$OMP END PARALLEL DO
    !::subtract M of a tile in points that are inside that tile in order to actually get H (only for CylindricalTiles as these actually calculate the B-field (divided by mu0)
    call SubtractMFromCylindricalTiles( H, tiles, pts, n_tiles, n_ele)
        
    end subroutine getFieldFromTiles
    
    
    subroutine SubtractMFromCylindricalTiles( H, tiles, pts, n_tiles, n_pts)
    real,dimension(n_pts,3),intent(inout) :: H
    type(MagTile),dimension(n_tiles),intent(in) :: tiles
    real,dimension(n_pts,3),intent(in) :: pts
    integer,intent(in) :: n_tiles,n_pts
    
    real,dimension(:),allocatable :: r,theta,z
    real :: rmin,rmax,thetamin,thetamax,zmin,zmax
    integer :: i
    
    allocate( r(n_pts), theta(n_pts), z(n_pts) )
    !::Convert from Cartesian to cylindrical coordinates
    r = sqrt( pts(:,1)**2 + pts(:,2)**2 )
    theta = acos( pts(:,1) / r )
    !Correct for being in either the third or foruth quadrants
    where ( pts(:,2) .lt. 0 )
        theta = 2 * pi - theta    
    endwhere
    
    z = pts(:,3)
    
    !::loop over each tile
    do i=1,n_tiles
        if ( tiles(i)%tileType .eq. tileTypeCylPiece ) then
            rmin = tiles(i)%r0 - tiles(i)%dr/2
            rmax = tiles(i)%r0 + tiles(i)%dr/2
        
            thetamin = tiles(i)%theta0 - tiles(i)%dtheta/2
            thetamax = tiles(i)%theta0 + tiles(i)%dtheta/2
        
            zmin = tiles(i)%z0 - tiles(i)%dz/2
            zmax = tiles(i)%z0 + tiles(i)%dz/2
        
            !::then the point is inside the tile and the M-vector of that tile should be subtracted in order to get H
            where( rmin .le. r .AND. r .le. rmax .AND. thetamin .le. theta .AND. theta .le. thetamax .AND. zmin .le. z .AND. z .le. zmax )
                H(:,1) = H(:,1) - tiles(i)%M(1)
                H(:,2) = H(:,2) - tiles(i)%M(2)
                H(:,3) = H(:,3) - tiles(i)%M(3)
            endwhere
        endif
        
    enddo
    deallocate(r,theta,z)
    end subroutine SubtractMFromCylindricalTiles
    
    !::
    !::Specific implementation for a single cylindrical tile
    !::
    subroutine getFieldFromCylTile( cylTile, H, pts, n_ele, N_out )
    type(MagTile), intent(inout) :: cylTile
    real,dimension(n_ele,3),intent(inout) :: H
    real,dimension(n_ele,3) :: pts
    integer,intent(in) :: n_ele
    real,dimension(3,3,n_ele),intent(inout) :: N_out
    real,dimension(:),allocatable :: r,x,phi
    real,dimension(:,:),allocatable :: pts_local
    real :: phi_orig,z_orig
    real,dimension(3) :: M_orig,M_tmp
    real,dimension(3,3) :: N,Rz
    integer :: i
    
      !::Run the calculation
      allocate( r(n_ele), x(n_ele), phi(n_ele), pts_local(n_ele,3) )
      
      !::Include the offset between the global coordinate system and the tile's coordinate system
      !::the pts array is always in global coordinates
      pts_local(:,1) = pts(:,1) - cylTile%offset(1)
      pts_local(:,2) = pts(:,2) - cylTile%offset(2)
      pts_local(:,3) = pts(:,3) - cylTile%offset(3)
      
      x(:) = 0.
      phi(:) = 0.
      H(:,:) = 0.
      !the length of the radius vector
      r = sqrt( pts_local(:,1)**2 + pts_local(:,2)**2 )
      !Find the rotation angle. It is assumed that r != 0 in all points since the tensor-field diverges here      
      phi = acos( pts_local(:,1) / r )
      !correct for being in the third or foruth quadrants
      where ( pts_local(:,2) .lt. 0 )
          phi = 2 * pi - phi !phi = -phi
      endwhere
      
      !save the original orientation of the tile
      phi_orig = cylTile%theta0
      z_orig = cylTile%z0      
      M_orig = cylTile%M
      do i=1,n_ele
          !Offset the angle
          cylTile%theta0 = phi_orig - phi(i)
          !Offset the z-coordinate
          cylTile%z0 = z_orig - pts_local(i,3)
          call getN_CylPiece( cylTile, r(i), N )
          N_out(:,:,i) = N
          !Get the rotation vector
          call getRotZ( -phi(i), Rz )
          !Rotate the magnetization vector
          M_tmp = matmul( Rz, M_orig )
          !find the field at the rotated point
          H(i,:) = 1./(4.*pi) * matmul( N, M_tmp )
          !get the negative rotation
          call getRotZ( phi(i), Rz )
          !rotate the field back to the original orientation
          H(i,:) = matmul( Rz, H(i,:) )
          
          !Revert the changes in angle and z position
          cylTile%theta0 = phi_orig
          cylTile%z0 = z_orig
      enddo
      
      deallocate(r,x,phi,pts_local)
    
    end subroutine getFieldFromCylTile
    
    !::
    !::Returns the magnetic field from a rectangular prism
    !::
    subroutine getFieldFromRectangularPrismTile( prismTile, H, pts, n_ele, N_out )
    type(MagTile),intent(in) :: prismTile
    real,dimension(n_ele,3),intent(inout) :: H
    real,dimension(n_ele,3) :: pts
    integer,intent(in) :: n_ele
    real,dimension(3,3,n_ele),intent(inout) :: N_out
    real,dimension(3) :: diffPos,dotProd
    real,dimension(3,3) :: rotMat,rotMatInv    
    integer :: i
    
    !::get the rotation matrices
    call getRotationMatrices( prismTile, rotMat, rotMatInv)
    
    do i=1,n_ele
        !::1. The relative position vector between the origo of the current tile and the point at which the tensor is required
        diffPos = pts(i,:) - prismTile%offset
        
        !::2. rotate the position vector according to the rotation of the prism
        diffPos = matmul( rotMat, diffPos )
        
        !::3. get the demag tensor
        call getN_prism_3D( prismTile, diffPos, N_out(:,:,i) )
        
        !::4. rotate the magnetization vector from the global system to the rotated frame and get the field (dotProd)
        call getDotProd( N_out(:,:,i), matmul( rotMat, prismTile%M ), dotProd )
        
        !::5. Rotate the resulting field back to the global coordinate system
        dotProd = matmul( rotMatInv, dotProd )        
        
        !::. Update the solution. The minus sign comes from the definition of the demag tensor (the demagfield is assumed negative)
        H(i,:) = -dotProd
    enddo
    end subroutine getFieldFromRectangularPrismTile
    

!---------------------------------------------------------------------------------------!
    !::Helper routines

    !::Dot product between (3,3) matrix and (1,3) vector
    subroutine getDotProd( N, M, dot_prod )
        real,intent(in),dimension(3,3) :: N
        real,intent(in),dimension(3) :: M
        real,intent(inout),dimension(3) :: dot_prod

        dot_prod(1) = sum(N(1,:) * M(:))
        dot_prod(2) = sum(N(2,:) * M(:))
        dot_prod(3) = sum(N(3,:) * M(:))


    end subroutine getDotProd
    
    
    !::
    !::Calculates the rotation matrix and its inverse for a given tile
    !::
    subroutine getRotationMatrices( tile, rotMat, rotMatInv)
    type(MagTile),intent(in) :: tile    
    real,intent(inout),dimension(3,3) :: rotMat,rotMatInv
    
    real,dimension(3,3) :: RotX,RotY,RotZ


    !::The minus sign is important since a rotated prism can be represented with a rotation about the given axis in the opposite direction
    call getRotX( -tile%rotAngles(1), RotX )
    call getRotY( -tile%rotAngles(2), RotY )
    call getRotZ( -tile%rotAngles(3), RotZ )

    !::Find the rotation matrix as defined by yaw, pitch and roll
    !::Rot = RotZ * RotY * RotX
    rotMat = matmul( matmul( RotZ, RotY ), RotX )
    !::As the rotation matrices are orthogonal we have that înv(R) = transp(R)
    !:: And thus inv(R1R2) = transp(R2)transp(R1) (note the change in multiplication order)
    !::and that transp( R1R2R3 ) = transp(R3)transp(R1R2) = transp(R3) transp(R2) transp(R1)
    !:: inv(Rot) = inv(RotZ * RotY * RotX ) =>
    !:: inv(Rot) = transp(Rot) = transp( RotZ * RotY * RotX )
    !::                        = transp( RotX ) * transp( RotZ * RotY )
    !::                        = transp( RotX ) * transp( RotY ) * transp( RotZ )
    rotMatInv = matmul( matmul( transpose(RotX), transpose(RotY) ), transpose( RotZ ) )
    

    end subroutine getRotationMatrices

    !::
    !::Returns the rotation matrix for a rotation of angle radians about the x-axis
    !::
    subroutine getRotX( angle, rotMat )
    real,intent(in) :: angle
    real,intent(inout),dimension(3,3) :: rotMat

    !::fortran matrices are (row,col)
    rotMat(1,1) = 1
    !::Top row
    rotMat(1,2:3) = 0
    !::left column
    rotMat(2:3,1) = 0
    rotMat(2,2) = cos(angle)
    rotMat(3,3) = cos(angle)
    rotMat(2,3) = -sin(angle)
    rotMat(3,2) = sin(angle)

    end subroutine getRotX

    !::
    !::Returns the rotation matrix for a rotation of angle radians about the y-axis
    !::
    subroutine getRotY( angle, rotMat )
    real,intent(in) :: angle
    real,intent(inout),dimension(3,3) :: rotMat

    !::fortran matrices are (row,col)
    !::top row
    rotMat(1,1) = cos(angle)
    rotMat(1,2) = 0
    rotMat(1,3) = sin(angle)
    !::middle row
    rotMat(2,1) = 0
    rotMat(2,2) = 1
    rotMat(2,3) = 0
    !::bottom row
    rotMat(3,1) = -sin(angle)
    rotMat(3,2) = 0
    rotMat(3,3) = cos(angle)

    end subroutine getRotY
    
    !::Get rotation matrix for rotation about the z-axis
    subroutine getRotZ( phi, Rz )
        real*8,intent(in) :: phi
        real*8,dimension(3,3),intent(inout) :: Rz
        Rz(:,:) = 0
    
        Rz(1,1) = cos(phi)
        Rz(1,2) = -sin(phi)
        Rz(2,1) = sin(phi)
        Rz(2,2) = cos(phi)
        
        Rz(3,3) = 1.
    
    end subroutine getRotZ
    
    
    end module MagStat2GetSolution
    
  