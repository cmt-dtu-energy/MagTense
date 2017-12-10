module MagStat2GetSolution
    use TileNComponents

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
    
    integer :: i
    real,dimension(:,:),allocatable :: H_tmp
    
    
    allocate(H_tmp(n_ele,3))
    H(:,:) = 0.
    
    !$OMP PARALLEL DO PRIVATE(i,H_tmp)
    do i=1,n_tiles
        H_tmp(:,:) = 0.
        
        !::Here a selection of which subroutine to use should be done, i.e. whether the tile
        !:: is cylindrical, a prism or an ellipsoid
        call getFieldFromCylTile( tiles(i), H_tmp, pts, n_ele )
        !$OMP CRITICAL
        H = H + H_tmp
        !$OMP END CRITICAL
    enddo
    !$OMP END PARALLEL DO
    
    deallocate(H_tmp)
    end subroutine getFieldFromTiles
    
    !::
    !::Specific implementation for a single cylindrical tile
    !::
    subroutine getFieldFromCylTile( cylTile, H, pts, n_ele )
    type(MagTile), intent(inout) :: cylTile
    real,dimension(n_ele,3),intent(inout) :: H
    real,dimension(n_ele,3) :: pts
    integer,intent(in) :: n_ele
    real,dimension(:),allocatable :: r,x,phi
    real :: phi_orig,z_orig
    real,dimension(3) :: M_orig,M_tmp
    real,dimension(3,3) :: N,Rz
    integer :: i
    
      !::Run the calculation
      allocate( r(n_ele), x(n_ele), phi(n_ele))
      
      x(:) = 0.
      phi(:) = 0.
      H(:,:) = 0.
      !the length of the radius vector
      r = sqrt( pts(:,1)**2 + pts(:,2)**2 )
      !Find the rotation angle. It is assumed that r != 0 in all points since the tensor-field diverges here      
      phi = acos( pts(:,1) / r )
      !Change the sign if y is negative
      where ( pts(:,2) .lt. 0 )
          phi = -phi
      endwhere
      
      !save the original orientation of the tile
      phi_orig = cylTile%theta0
      z_orig = cylTile%z0      
      M_orig = cylTile%M
      do i=1,n_ele
          !Offset the angle
          cylTile%theta0 = phi_orig - phi(i)
          !Offset the z-coordinate
          cylTile%z0 = z_orig - pts(i,3)
          call getN( cylTile, r(i), N )
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
      
      deallocate(r,x,phi)
    
    end subroutine getFieldFromCylTile
    
    subroutine getRotZ( phi, Rz )
        real*8,intent(in) :: phi
        real*8,dimension(3,3),intent(inout) :: Rz
        Rz(:,:) = 0
    
        Rz(1,1) = cos(phi)
        Rz(1,2) = -sin(phi)
        Rz(2,1) = sin(phi)
        Rz(2,2) = cos(phi)
    
        end subroutine getRotZ
    
    end module MagStat2GetSolution
    
  