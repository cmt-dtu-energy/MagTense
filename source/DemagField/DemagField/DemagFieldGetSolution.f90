
    module DemagFieldGetSolution
    
    use TileNComponents

    implicit none
      
    contains
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !!General function to call given a number of (different) tiles
    !!The demagnetization tensor Nout can be provided. If this is the case, the function uses this for calculations. Other it is computed.
    !!Input arguments:
    !!@param tiles: Array of type MagTile, size n_tiles
    !!@param H the output magnetic field, size [n_ele,3]
    !!@param pts the points at which the field should be evaluated, size [n_ele,3]
    !!@param n_tiles, the number of tiles
    !!@param n_ele, the number of points at which to evaluate the field
    !!@param Nout the demag tensor calculated by this function (size (n_tiles,n_pts,3,3) )
    !!
    subroutine getFieldFromTiles( tiles, H, pts, n_tiles, n_ele, Nout )
        type(MagTile),intent(inout),dimension(n_tiles) :: tiles
        real,dimension(n_ele,3),intent(inout) :: H
        real,dimension(n_ele,3),intent(in) :: pts
        integer,intent(in) :: n_tiles,n_ele
        real,dimension(:,:,:,:),allocatable,optional :: Nout
    
        integer :: i,prgCnt,tid,prog,OMP_GET_THREAD_NUM
        real,dimension(:,:),allocatable :: H_tmp
        integer,parameter :: cbCnt = 10
        logical :: useStoredN
    
        !!If Nout is provided, the function uses this for calculations
        if ( present( Nout ) ) then 
            if ( .NOT. allocated(Nout) ) then
                allocate(Nout(n_tiles,n_ele,3,3))
                Nout(:,:,:,:) = 0.
                useStoredN = .false.
            else
                useStoredN = .true.  
            endif
        else
            useStoredN = .false.
        endif
                
        allocate(H_tmp(n_ele,3))
        H(:,:) = 0.
    
        prgCnt = 0
        prog = 0
    
        ! $OMP PARALLEL DO PRIVATE(i,H_tmp)    
        do i=1,n_tiles
        
            !Make sure to allocate H_tmp on the heap and for each thread
            ! $OMP CRITICAL
            H_tmp(:,:) = 0.        
        
            ! $OMP END CRITICAL
            !! Here a selection of which subroutine to use should be done, i.e. whether the tile
            !! is cylindrical, a prism or an ellipsoid or another geometry
            select case (tiles(i)%tileType )
            case (tileTypeCylPiece)
                if ( present(Nout) ) then
                    call getFieldFromCylTile( tiles(i), H_tmp, pts, n_ele, Nout(i,:,:,:), useStoredN )    
                else
                    call getFieldFromCylTile( tiles(i), H_tmp, pts, n_ele )
                endif
            
            case (tileTypePrism)
                if ( present(Nout) ) then
                    call getFieldFromRectangularPrismTile( tiles(i), H_tmp, pts, n_ele, Nout(i,:,:,:), useStoredN )
                else
                    call getFieldFromRectangularPrismTile( tiles(i), H_tmp, pts, n_ele )
                endif
            case (tileTypeCircPiece )
                 if ( present(Nout) ) then
                    call getFieldFromCircPieceTile( tiles(i), H_tmp, pts, n_ele, Nout(i,:,:,:), useStoredN )
                else
                    call getFieldFromCircPieceTile( tiles(i), H_tmp, pts, n_ele )
                endif
            case (tileTypeCircPieceInverted )
                 if ( present(Nout) ) then
                    call getFieldFromCircPieceInvertedTile( tiles(i), H_tmp, pts, n_ele, Nout(i,:,:,:), useStoredN )
                else
                    call getFieldFromCircPieceInvertedTile( tiles(i), H_tmp, pts, n_ele )
                endif
            case (tileTypeEllipsoid)
                !!@todo add the existing code for spheroids, and correct Ellipsoids to Spheroids
            case (tileTypePlanarCoil )
                if ( present(Nout) ) then
                    call getFieldFromPlanarCoilTile( tiles(i), H_tmp, pts, n_ele, Nout(i,:,:,:), useStoredN )
                else
                    call getFieldFromPlanarCoilTile( tiles(i), H_tmp, pts, n_ele )            
                endif
            
            case default        
            
            end select
        
            !!@todo Can this code be removed as well as the variables prgCnt and prog?
            ! $OMP CRITICAL
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
           ! $OMP END CRITICAL
       
        
        
        enddo
        ! $OMP END PARALLEL DO
    
        deallocate(H_tmp)
    
        !!Subtract M of a tile in points that are inside that tile in order to actually get H (only for CylindricalTiles as these actually calculate the B-field (divided by mu0)
        !!@todo Should this function not only be called if tileTypeCylPiece?
        call SubtractMFromCylindricalTiles( H, tiles, pts, n_tiles, n_ele)
    
    end subroutine getFieldFromTiles
    
    
    
!---------------------------------------------------------------------------------------!
!------------------------ Specific tile geometries -------------------------------------!
!---------------------------------------------------------------------------------------!    

    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !!Returns the magnetic field from a rectangular prism
    !!
    subroutine getFieldFromRectangularPrismTile( prismTile, H, pts, n_ele, N_out, useStoredN )
        type(MagTile),intent(in) :: prismTile
        real,dimension(n_ele,3),intent(inout) :: H
        real,dimension(n_ele,3) :: pts
        integer,intent(in) :: n_ele
        real,dimension(n_ele,3,3),intent(inout),optional :: N_out
        logical,intent(in),optional :: useStoredN
    
        procedure (N_tensor_subroutine), pointer :: N_tensor => null ()
    
        N_tensor => getN_prism_3D
    
        !! Check to see if we should use symmetry
        if ( prismTile%exploitSymmetry .eq. 1 ) then        
            call getFieldFromTile_symm(prismTile, H, pts, n_ele, N_tensor, N_out, useStoredN )        
        else        
            call getFieldFromTile( prismTile, H, pts, n_ele, N_tensor, N_out, useStoredN )       
        endif
    
        !!@todo Can this be removed?
        !! The minus sign comes from the definition of the demag tensor (the demagfield is assumed negative)
        !! Change in the tensor subroutine in order to make the behavior of the tensor components of the various geometries conform
        !H = -1. * H
    
    end subroutine getFieldFromRectangularPrismTile      
        
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !!Returns the magnetic field from a cylindrical tile
    !!
    subroutine getFieldFromCylTile( cylTile, H, pts, n_ele, N_out, useStoredN )
        type(MagTile), intent(inout) :: cylTile
        real,dimension(n_ele,3),intent(inout) :: H
        real,dimension(n_ele,3) :: pts
        integer,intent(in) :: n_ele
        real,dimension(n_ele,3,3),intent(inout),optional :: N_out
        real,dimension(:),allocatable :: r,x,phi
        real,dimension(:,:),allocatable :: pts_local
        real :: phi_orig,z_orig
        real,dimension(3) :: M_orig !,M_tmp
        real,dimension(3,3) :: N,Rz,Rz_inv
        integer :: i
        logical,intent(in),optional :: useStoredN
    
    
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
      
        r = sqrt( pts_local(:,1)**2 + pts_local(:,2)**2 ) !< The length of the radius vector
      
        !!Find the rotation angle. It is assumed that r != 0 in all points since the tensor-field diverges here      
        phi = acos( pts_local(:,1) / r )
      
        !!Correct for being in the third or fourth quadrants
        where ( pts_local(:,2) .lt. 0 )
            !phi = 2 * pi - phi
            phi = -phi
        endwhere
      
        !!Save the original orientation of the tile
        phi_orig = cylTile%theta0
        z_orig = cylTile%z0      
        M_orig = cylTile%M
        do i=1,n_ele
          
            cylTile%theta0 = phi_orig - phi(i) !< Offset the angle
            cylTile%z0 = z_orig - pts_local(i,3) !< Offset the z-coordinate
          
            if ( present( useStoredN ) .eqv. .true. ) then
                if ( useStoredN .eqv. .false. ) then
                    call getN_CylPiece( cylTile, r(i), N )
              
                    !<Change basis from rotation trick coordinate system to the local system of the tile
                    !!The magnetic field in the global coordinate system is given by:
                    !!H = Rot(phi) * ( N *  (Rot(-phi) * M) )
                    !!Get the rotation matrix for rotating the magnetization with the trick-angle
                    call getRotZ( -phi(i), Rz_inv )
                    !!Get the rotation matrix for rotating the magnetic field back to the original coordinate system
                    call getRotZ( phi(i), Rz )
                    N = matmul( Rz, N )
                    N = matmul( N, Rz_inv )
                  
                    N_out(i,:,:) = N
              
                else
                    N = N_out(i,:,:)
                endif
            else
                call getN_CylPiece( cylTile, r(i), N )
                !<Change basis from rotation trick coordinate system to the local system of the tile
                !!The magnetic field in the global coordinate system is given by:
                !!H = Rot(phi) * ( N *  (Rot(-phi) * M) )
                !!Get the rotation matrix for rotating the magnetization with the trick-angle
                call getRotZ( -phi(i), Rz_inv )
                !!Get the rotation matrix for rotating the magnetic field back to the original coordinate system
                call getRotZ( phi(i), Rz )
                N = matmul( Rz, N )
                N = matmul( N, Rz_inv )
                  
                N_out(i,:,:) = N
            endif
          
            !!@todo Can this code be removed?
            !Rotate the magnetization vector
            !M_tmp = matmul( Rz, M_orig )
            !find the field at the rotated point
            !H(i,:) = matmul( N, M_tmp )
            !get the negative rotation
            !call getRotZ( phi(i), Rz )
            !rotate the field back to the original orientation
            !H(i,:) = matmul( Rz, H(i,:) )
            H(i,:) = matmul( N, M_orig )
          
            !!Revert the changes in angle and z position
            cylTile%theta0 = phi_orig
            cylTile%z0 = z_orig
        enddo
      
        deallocate(r,x,phi,pts_local)
    
    end subroutine getFieldFromCylTile
    
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !!Function to calcuate H within a cylindrical tile that is rotated, as for CylindricalTiles it is actually the B-field divided by mu0 that is calculated
    !!
    subroutine SubtractMFromCylindricalTiles( H, tiles, pts, n_tiles, n_pts)
        real,dimension(n_pts,3),intent(inout) :: H
        type(MagTile),dimension(n_tiles),intent(in) :: tiles
        real,dimension(n_pts,3),intent(in) :: pts
        integer,intent(in) :: n_tiles,n_pts
    
        real,dimension(:),allocatable :: r,theta,z
        real,dimension(:,:),allocatable :: pts_local
        real :: rmin,rmax,thetamin,thetamax,zmin,zmax
        integer :: i
    
        allocate( r(n_pts), theta(n_pts), z(n_pts),  pts_local(n_pts,3) )
              
        !::loop over each tile
        do i=1,n_tiles
            if ( tiles(i)%tileType .eq. tileTypeCylPiece ) then
            
                !::Include the offset between the global coordinate system and the tile's coordinate system
                !::the pts array is always in global coordinates
                pts_local(:,1) = pts(:,1) - tiles(i)%offset(1)
                pts_local(:,2) = pts(:,2) - tiles(i)%offset(2)
                pts_local(:,3) = pts(:,3) - tiles(i)%offset(3)
                !::Convert from Cartesian to cylindrical coordinates
                r = sqrt( pts_local(:,1)**2 + pts_local(:,2)**2 )
                theta = acos( pts_local(:,1) / r )
                !Correct for being in either the third or foruth quadrants
                where ( pts_local(:,2) .lt. 0 )
                    theta = 2 * pi - theta    
                endwhere
    
                z = pts_local(:,3)
            
                rmin = tiles(i)%r0 - tiles(i)%dr/2
                rmax = tiles(i)%r0 + tiles(i)%dr/2
        
                thetamin = tiles(i)%theta0 - tiles(i)%dtheta/2
                thetamax = tiles(i)%theta0 + tiles(i)%dtheta/2                   
            
                zmin = tiles(i)%z0 - tiles(i)%dz/2
                zmax = tiles(i)%z0 + tiles(i)%dz/2
        
                where( rmin .le. r .AND. r .lt. rmax .AND. zmin .le. z .AND. z .lt. zmax  .AND. thetamin .le. theta .AND. theta .lt. thetamax .OR. &
                       rmin .le. r .AND. r .lt. rmax .AND. zmin .le. z .AND. z .lt. zmax  .AND. thetamin - 2*pi .le. theta .AND. theta .lt. thetamax - 2*pi .OR. &
                       rmin .le. r .AND. r .lt. rmax .AND. zmin .le. z .AND. z .lt. zmax  .AND. thetamin + 2*pi .le. theta .AND. theta .lt. thetamax + 2*pi )
                    H(:,1) = H(:,1) - tiles(i)%M(1)
                    H(:,2) = H(:,2) - tiles(i)%M(2)
                    H(:,3) = H(:,3) - tiles(i)%M(3)
                endwhere
            endif
        
        enddo
        deallocate(r,theta,z,pts_local)
    end subroutine SubtractMFromCylindricalTiles
    
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !!Returns the magnetic field from a tile that is a piece of circle
    !!
    subroutine getFieldFromCircPieceTile( circTile, H, pts, n_ele, N_out, useStoredN )
        type(MagTile),intent(in) :: circTile
        real,dimension(n_ele,3),intent(inout) :: H
        real,dimension(n_ele,3) :: pts
        integer,intent(in) :: n_ele
        real,dimension(n_ele,3,3),intent(inout),optional :: N_out
        logical,intent(in),optional :: useStoredN
    
        procedure (N_tensor_subroutine), pointer :: N_tensor => null ()
    
        N_tensor => getN_circPiece
        !! check to see if we should use symmetry
        if ( circTile%exploitSymmetry .eq. 1 ) then
        
            call getFieldFromTile_symm( circTile, H, pts, n_ele, N_tensor, N_out, useStoredN )
        else
            call getFieldFromTile( circTile, H, pts, n_ele, N_tensor, N_out, useStoredN )
            !!@todo Can this code be removed?
            !!::get the rotation matrices
            !call getRotationMatrices( circTile, rotMat, rotMatInv)
            !
            !do i=1,n_ele
            !    !::1. The relative position vector between the origo of the current tile and the point at which the tensor is required
            !    diffPos = pts(i,:) - circTile%offset
            !
            !    !::2. rotate the position vector according to the rotation of the tile, i.e. rotate
            !    !::the position vector to align with the local coordinate system of the tile
            !    diffPos = matmul( rotMat, diffPos )
            !
            !    !::3. get the demag tensor                
            !    if ( present( useStoredN ) .eq. .true. ) then
            !        if ( useStoredN .eq. .false. ) then
            !             call getN_circPiece( circTile, diffPos, N )
            !             N_out(i,:,:) = N
            !        else
            !            N = N_out(i,:,:)
            !        endif
            !    else
            !        call getN_circPiece( circTile, diffPos, N )
            !    endif
            !
            !    !::4. rotate the magnetization vector from the global system to the rotated frame and get the field (dotProd)
            !    call getDotProd( N, matmul( rotMat, circTile%M ), dotProd )
            !
            !    !::5. Rotate the resulting field back to the global coordinate system
            !    dotProd = matmul( rotMatInv, dotProd )        
            !
            !    !::. Update the solution. 
            !    H(i,:) = dotProd
            !enddo
        endif
    
    end subroutine getFieldFromCircPieceTile
    
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !!Returns the magnetic field from a piece of circle that is inverted, i.e. pointing inwards
    !!
    subroutine getFieldFromCircPieceInvertedTile( circTile, H, pts, n_ele, N_out, useStoredN )
        type(MagTile),intent(in) :: circTile
        real,dimension(n_ele,3),intent(inout) :: H
        real,dimension(n_ele,3) :: pts
        integer,intent(in) :: n_ele
        real,dimension(n_ele,3,3),intent(inout),optional :: N_out
        logical,intent(in),optional :: useStoredN  
    
        procedure (N_tensor_subroutine), pointer :: N_tensor => null ()
    
        N_tensor => getN_circPiece_Inv
        !! check to see if we should use symmetry
        if ( circTile%exploitSymmetry .eq. 1 ) then
        
            call getFieldFromTile_symm( circTile, H, pts, n_ele, N_tensor, N_out, useStoredN )
        else
            call getFieldFromTile( circTile, H, pts, n_ele, N_tensor, N_out, useStoredN )
            !!@todo Can this code be removed?
            !
            !!::get the rotation matrices
            !call getRotationMatrices( circTile, rotMat, rotMatInv)
            !
            !do i=1,n_ele
            !    !::1. The relative position vector between the origo of the current tile and the point at which the tensor is required
            !    diffPos = pts(i,:) - circTile%offset
            !
            !    !::2. rotate the position vector according to the rotation of the prism
            !    diffPos = matmul( rotMat, diffPos )
            !
            !    !::3. get the demag tensor                
            !    if ( present( useStoredN ) .eq. .true. ) then
            !        if ( useStoredN .eq. .false. ) then
            !             call getN_circPiece_Inv( circTile, diffPos, N )
            !             N_out(i,:,:) = N
            !        else
            !            N = N_out(i,:,:)
            !        endif
            !    else
            !        call getN_circPiece_Inv( circTile, diffPos, N )
            !    endif
            !
            !    !::4. rotate the magnetization vector from the global system to the rotated frame and get the field (dotProd)
            !    call getDotProd( N, matmul( rotMat, circTile%M ), dotProd )
            !
            !    !::5. Rotate the resulting field back to the global coordinate system
            !    dotProd = matmul( rotMatInv, dotProd )        
            !
            !    !::. Update the solution. 
            !    H(i,:) = dotProd
            !enddo
        endif
    
    end subroutine getFieldFromCircPieceInvertedTile
    
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !!Returns the field from a planar coil with N windings and the current I running through it
    !!tile%M is assumed to contain the current as tile%M = [I,0,0] in the unit of Amps
    !!@todo Rotation has not been implemented yet for this tile type
    !!
    subroutine getFieldFromPlanarCoilTile( tile, H, pts, n_ele, N_out, useStoredN )
        type(MagTile),intent(in) :: tile
        real,dimension(n_ele,3),intent(inout) :: H
        real,dimension(n_ele,3) :: pts
        integer,intent(in) :: n_ele
        real,dimension(n_ele,3,3),intent(inout),optional :: N_out
        logical,intent(in),optional :: useStoredN
        real,dimension(3) :: diffPos,dotProd
        real,dimension(3,3) :: N
        !real,dimension(3,3) :: rotMat,rotMatInv
        integer :: i
    
        !::get the rotation matrices
        !call getRotationMatrices( tile, rotMat, rotMatInv)
    
        do i=1,n_ele
            !! The relative position vector between the origo of the current tile and the point at which the tensor is required
            diffPos = pts(i,:) - tile%offset
        
            !! rotate the position vector according to the rotation of the prism
           ! diffPos = matmul( rotMat, diffPos )
        
            !! Get the demag tensor                
            if ( present( useStoredN ) .eqv. .true. ) then
                if ( useStoredN .eqv. .false. ) then
                    call getN_PlanarCoil( tile, diffPos, N )
                 
                    N_out(i,:,:) = N
                else
                    N = N_out(i,:,:)
                endif
            else
                call getN_PlanarCoil( tile, diffPos, N )
            endif
        
            !! Rotate the magnetization vector from the global system to the rotated frame and get the field (dotProd)
            !call getDotProd( N, matmul( rotMat, tile%M ), dotProd )
            call getDotProd( N, tile%M, dotProd )
        
            !! Rotate the resulting field back to the global coordinate system
            !dotProd = matmul( rotMatInv, dotProd )        
        
            !! Update the solution. 
            H(i,:) = dotProd
        enddo
    end subroutine getFieldFromPlanarCoilTile
    

    
!---------------------------------------------------------------------------------------!
!------------------------ General tile routines ----------------------------------------!
!---------------------------------------------------------------------------------------!
    
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !!Calculates the actual field from a tile of a given geometry
    !!
    subroutine getFieldFromTile( tile, H, pts, n_ele, N_tensor, N_out, useStoredN )
        type(MagTile),intent(in) :: tile
        real,dimension(n_ele,3),intent(inout) :: H
        real,dimension(n_ele,3),intent(in) :: pts
        integer,intent(in) :: n_ele
        procedure (N_tensor_subroutine), pointer, intent(in) :: N_tensor
        real,dimension(n_ele,3,3),intent(inout),optional :: N_out
        logical,intent(in),optional:: useStoredN
            
        real,dimension(3,3) :: rotMat,rotMatInv,N
        integer :: i
        real,dimension(3) :: diffPos,dotProd        
    
        !! get the rotation matrices
        call getRotationMatrices( tile, rotMat, rotMatInv)
    
        do i=1,n_ele
            !! The relative position vector between the origo of the current tile and the point at which the tensor is required
            diffPos = pts(i,:) - tile%offset
        
            !! Rotate the position vector according to the rotation of the prism
            diffPos = matmul( rotMat, diffPos )
        
            !! Get the demag tensor                
            if ( present( useStoredN ) .eqv. .true. ) then
                if ( useStoredN .eqv. .false. ) then
                     call N_tensor( tile, diffPos, N )
                     N_out(i,:,:) = N
                else
                    N = N_out(i,:,:)
                endif
            else
                call N_tensor( tile, diffPos, N )
            endif
        
            !! Rotate the magnetization vector from the global system to the rotated frame and get the field (dotProd)
            call getDotProd( N, matmul( rotMat, tile%M ), dotProd )
        
            !! Rotate the resulting field back to the global coordinate system
            dotProd = matmul( rotMatInv, dotProd )        
        
            !! Update the solution.
            H(i,:) = dotProd
        enddo
    
    end subroutine getFieldFromTile
    
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !> This function assumes 8-fold symmetry in the tile, i.e. that it has siblings at 7 other locations according to the conventional
    !! mirror operations about the global coordinate planes
    !! @param prismTile is the tile in question
    !! @param H the field vector to be updated by this calculation
    !! @param pts the points at which to find the solution (n_ele,3)
    !! @param n_ele integer, no. of points at which to evaluate
    !! @param N_tensor pointer to the subroutine that calculates the N-tensor
    !! @N_out the resulting demag tensor for each of the points in question 
    !! @ useStoredN logical. If true then use the values in N_out else calculate the tensor
    !!
    subroutine getFieldFromTile_symm(tile, H, pts, n_ele, N_tensor, N_out, useStoredN )
        type(MagTile),intent(in) :: tile
        real,dimension(n_ele,3),intent(inout) :: H
        real,dimension(n_ele,3),intent(in) :: pts
        integer,intent(in) :: n_ele
        procedure (N_tensor_subroutine), pointer, intent(in) :: N_tensor
        real,dimension(n_ele,3,3),intent(inout),optional :: N_out
        logical,intent(in),optional:: useStoredN
        
        real,dimension(:,:,:),allocatable :: N_tmp    
        real,dimension(3,3) :: rotMat,rotMatInv,N
        integer :: i,j
        real,dimension(3) :: diffPos,dotProd
        real,dimension(8,3,3) :: symm_op_M, symm_op_H
        !!@todo Why is this temporary variable used instead of just useStoredN
        logical :: useStoredN_tmp
    
        !! get the rotation matrices for the tile
        call getRotationMatrices( tile, rotMat, rotMatInv)
    
        !! get the symmetry operation matrices
        !! symm_op_M is to be applied to the magnetization vector while
        !! symm_op_H is to be applied to the full expression and the diffPos
        call getSymmOpMatrices( symm_op_M, symm_op_H, tile%symmetryOps )
    
        allocate( N_tmp(n_ele,3,3) )
        N_tmp(:,:,:) = 0.
    
        if ( present( useStoredN ) .eqv. .false. ) then
            useStoredN_tmp = .false.
        else
            useStoredN_tmp = useStoredN
        endif
    
        !! go through each point
        do i=1,n_ele                  
        
            if ( useStoredN_tmp .eqv. .false. ) then
            
                !! loop over each symmetry operation
                do j=1,8                        
                
                    !! The relative position vector between the origo of the current tile and the point at which the tensor is required
                    !! apply symmetry operation to the point of interest (and leave the tile intact)
                    diffPos = matmul( symm_op_H(j,:,:), pts(i,:) ) - tile%offset
                
                    !! rotate the position vector according to the rotation of the tile
                    diffPos = matmul( rotMat, diffPos )            
                
                    !! do the calculation
                    call N_tensor( tile, diffPos, N )
                    
                    !! Change the sign of N according to the mirror symmetry (or anti-symmetry)                
                    N = matmul( symm_op_H(j,:,:), matmul( N, symm_op_M(j,:,:) ) )
                    !! remember to add the solution as we are exploting 8-fold symmetry
                    N_tmp(i,:,:) = N + N_tmp(i,:,:)
                enddo                                     
            else                
                N_tmp(i,:,:) = N_out(i,:,:)
            endif
        
            !! rotate the magnetization vector from the global system to the rotated frame and get the field (dotProd)
            call getDotProd( N_tmp(i,:,:), matmul( rotMat, tile%M ), dotProd )
        
            !! Rotate the resulting field back to the global coordinate system
            dotProd = matmul( rotMatInv, dotProd )        
        
            !! Update the solution. 
            H(i,:) = dotProd
        
        enddo
    
        if ( present(N_out) .eqv. .true. ) then
            N_out = N_tmp
        endif
    
    
        deallocate(N_tmp)
    end subroutine getFieldFromTile_symm
    
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !> Returns the 8 symmetry operation matrices for
    !! conventional mirroring about the three principal planes and their combinations
    !! @param symm_op_M array of 8 3x3 matrices containing the symmetry ops for the magnetization
    !! @param symm_op_H correspondingly but for the field and the point of interest
    !! @param symm_ops is a (3,1) array with the symmetry operations to perform along the three principal planes (xy, xz and yz)
    !! if @param symm_ops is 1 the operation is symmetric, if it is -1 the operation is anti-symmetric
    !!
    subroutine getSymmOpMatrices( symm_M, symm_H, symm_ops )
        real,dimension(8,3,3),intent(inout) :: symm_M, symm_H
        real,dimension(3),intent(in) :: symm_ops
        real,dimension(3,3) :: SymmX,SymmY,SymmZ    
        real :: theta_x,theta_y,theta_z
    
        SymmX(:,:) = 0.
        SymmX(1,1) = -1.
        SymmX(2,2) = 1.
        SymmX(3,3) = 1.
    
        SymmY = SymmX
        SymmY(1,1) = 1.
        SymmY(2,2) = -1.
    
        SymmZ = SymmY
        SymmZ(2,2) = 1.
        SymmZ(3,3) = -1.  
    
        !! reset
        symm_M(:,:,:) = 0.    
    
        !! insert unit matrices
        symm_M(:,1,1) = 1.
        symm_M(:,2,2) = 1.
        symm_M(:,3,3) = 1.
    
        symm_H = symm_M
    
        !!the octants are defined so that 1-4 are the same as the normal quadrants and 5-8 likewise but for negative z
    
        !! first operation is unity, i.e. nothing happens
    
        !! from 1st to 2nd octant, i.e. about yz plane    
        !! apply the symmetric / anti-symmetric op and mirror back to the original position
        symm_M(2,:,:) = matmul( SymmX(:,:) * symm_ops(3), SymmX(:,:) )
    
        symm_H(2,:,:) = SymmX
    
        !!from the 1st to the 4th octant, i.e. about the xz-plane    
        symm_M(3,:,:) = matmul( SymmY(:,:) * symm_ops(2), SymmY(:,:) )
    
        symm_H(3,:,:) = SymmY
    
        !!from the 1st to the 5th octant, i.e. about the xy-plane    
        symm_M(4,:,:) = matmul( SymmZ(:,:) * symm_ops(1), SymmZ(:,:) )
    
        symm_H(4,:,:) = SymmZ
    
        !! from the 1st to the 3rd octant (through the 2nd octant), i.e. first about yz and then about xz
        symm_M(5,:,:) = matmul( symm_M(3,:,:) , symm_M(2,:,:) )
        symm_H(5,:,:) = matmul( symm_H(3,:,:) , symm_H(2,:,:) )
    
        !! from the 1st to the 6th octant, i.e. first about yz and the about xy
        symm_M(6,:,:) = matmul( symm_M(4,:,:), symm_M(2,:,:) )
        symm_H(6,:,:) = matmul( symm_H(4,:,:), symm_H(2,:,:) )
    
        !! from the 1st to the 8th octants i.e. first about xy and then about xz
        symm_M(7,:,:) = matmul( symm_M(4,:,:), symm_M(3,:,:) )
        symm_H(7,:,:) = matmul( symm_H(4,:,:), symm_H(3,:,:) )
    
        !! from the 1st to the 7th octant i.e. firsth through yz, then through xy and then through xz
        symm_M(8,:,:) = matmul( symm_M(2,:,:), matmul( symm_M(3,:,:), symm_M(4,:,:) ) )
        symm_H(8,:,:) = matmul( symm_H(2,:,:), matmul( symm_H(3,:,:), symm_H(4,:,:) ) )
    
    end subroutine getSymmOpMatrices
    

!---------------------------------------------------------------------------------------!
!------------------------------ Helper routines ----------------------------------------!
!---------------------------------------------------------------------------------------!
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Dot product between (3,3) matrix and (1,3) vector
    !!
    subroutine getDotProd( N, M, dot_prod )
        real,intent(in),dimension(3,3) :: N
        real,intent(in),dimension(3) :: M
        real,intent(inout),dimension(3) :: dot_prod

        dot_prod(1) = sum(N(1,:) * M(:))
        dot_prod(2) = sum(N(2,:) * M(:))
        dot_prod(3) = sum(N(3,:) * M(:))

    end subroutine getDotProd
    
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Calculates the rotation matrix and its inverse for a given tile
    !!
    subroutine getRotationMatrices( tile, rotMat, rotMatInv)
        type(MagTile),intent(in) :: tile    
        real,intent(inout),dimension(3,3) :: rotMat,rotMatInv
    
        real,dimension(3,3) :: RotX,RotY,RotZ

        !! The minus sign is important since a rotated prism can be represented with a rotation about the given axis in the opposite direction
        call getRotX( -tile%rotAngles(1), RotX )
        call getRotY( -tile%rotAngles(2), RotY )
        call getRotZ( -tile%rotAngles(3), RotZ )

        !::Find the rotation matrix as defined by yaw, pitch and roll
        !::Rot = RotZ * RotY * RotX
        rotMat = matmul( matmul( RotZ, RotY ), RotX )
        !::As the rotation matrices are orthogonal we have that inv(R) = transp(R)
        !:: And thus inv(R1R2) = transp(R2)transp(R1) (note the change in multiplication order)
        !::and that transp( R1R2R3 ) = transp(R3)transp(R1R2) = transp(R3) transp(R2) transp(R1)
        !:: inv(Rot) = inv(RotZ * RotY * RotX ) =>
        !:: inv(Rot) = transp(Rot) = transp( RotZ * RotY * RotX )
        !::                        = transp( RotX ) * transp( RotZ * RotY )
        !::                        = transp( RotX ) * transp( RotY ) * transp( RotZ )
        rotMatInv = matmul( matmul( transpose(RotX), transpose(RotY) ), transpose( RotZ ) )
    
    end subroutine getRotationMatrices

    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Returns the rotation matrix for a rotation of angle radians about the x-axis
    !!
    subroutine getRotX( angle, rotMat )
        real,intent(in) :: angle
        real,intent(inout),dimension(3,3) :: rotMat

        !! Fortran matrices are (row,col)
        rotMat(1,1) = 1
        !! Top row
        rotMat(1,2:3) = 0
        !! Left column
        rotMat(2:3,1) = 0
        rotMat(2,2) = cos(angle)
        rotMat(3,3) = cos(angle)
        rotMat(2,3) = -sin(angle)
        rotMat(3,2) = sin(angle)

    end subroutine getRotX

    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Returns the rotation matrix for a rotation of angle radians about the y-axis
    !!
    subroutine getRotY( angle, rotMat )
        real,intent(in) :: angle
        real,intent(inout),dimension(3,3) :: rotMat

        !! Fortran matrices are (row,col)
        !! Top row
        rotMat(1,1) = cos(angle)
        rotMat(1,2) = 0
        rotMat(1,3) = sin(angle)
        !! Middle row
        rotMat(2,1) = 0
        rotMat(2,2) = 1
        rotMat(2,3) = 0
        !! Bottom row
        rotMat(3,1) = -sin(angle)
        rotMat(3,2) = 0
        rotMat(3,3) = cos(angle)

    end subroutine getRotY
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Get rotation matrix for rotation about the z-axis
    !!
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
    
    end module DemagFieldGetSolution
    
  
