
    module IterateMagnetSolution

    use DemagFieldGetSolution
    use MagParameters
    use spline

    

    implicit none

     INTERFACE
                FUNCTION displayIteration_fct(err,err_max)
                REAL, INTENT(IN) :: err,err_max
                integer :: displayIteration_fct
                END FUNCTION displayIteration_fct
     END INTERFACE
    
    contains   
    
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Iterates over all tiles to find the self-consistent solution and returns the tiles with updated magnetization vectors
    !! @param tiles array of tiles each define a part of the geometry
    !! @param n the number of tiles
    !! @param stateFunction array of type statefunction defining the look-up table for soft materials
    !! @param n_stf size of the statefunction array
    !! @param T the temperature to use when looking up into the statefunction table
    !! @param err_max maximum allowed relative error for the iterations
    !! @param max_ite max. no. of iterations to be done
    !! @param dispIteFct pointer to a subroutine that can display how the iteration is going
    !! @param resumeIteration real defining whether the first iteration should assume a finite magnetization in the input tiles and use that as the starting point
    !! if @param resumeIteration is zero no resume is done. Else resume is assumed and the value of resumeIteration is used as the starting value of lambda (should thus be positive)
    !! Error state is returned in i_err:
    !! 0 no error
    !! -1 max number of iterations reached without satisfying the minimum error
    subroutine iterateMagnetization( tiles, n, stateFunction, n_stf, T, err_max, max_ite, dispIteFct, resumeIteration )
        type(MagTile),dimension(n),intent(inout) :: tiles
        integer,intent(in) :: n
        type(MagStateFunction),intent(in),dimension(n_stf) :: stateFunction    
        integer,intent(in) :: n_stf
        procedure(displayIteration_fct),pointer,intent(in) :: dispIteFct
        real,intent(in) :: resumeIteration 
        real,intent(in) :: T    
        real,optional :: err_max
        integer,optional :: max_ite 
    
        !!@todo This parameter is not used. Remove?
        integer,parameter :: cnt_max = 200
        
        logical :: done
        integer :: i,j,cnt,i_err,lambdaCnt
        real,dimension(3) :: M
        real,dimension(:,:),allocatable :: H,H_old
        real,dimension(:),allocatable :: Hnorm,Hnorm_old,err_val,Mnorm,Mnorm_old
        type(NStoreArr),dimension(:),allocatable :: Nstore
        real,dimension(3) :: pts    
        real :: H_par,H_trans_1,H_trans_2,err,lambda,tmp
        real,dimension(4) :: maxRelDiffArr
        real,dimension(3,3) :: rotMat,rotMatInv
        logical :: lCh
    
        !::set default error margin if nothing has been specified by the user
        if ( .not. present(err_max) ) err_max = 1e-10
        if ( .not. present(max_ite) ) max_ite = 100
        i_err = 0
    
        allocate(H(n,3),H_old(n,3),Hnorm(n),Hnorm_old(n),err_val(n),Mnorm(n),Mnorm_old(n))
        allocate(Nstore(n))
        
        !! Initialization
        H(:,:) = 0
        maxRelDiffArr(:) = 0.
        cnt = 0
        Mnorm(:) = 0.
        Mnorm_old(:) = 0.
    
        !! Relaxation parameter
        lambda = 1.
        lambdaCnt = 0
    
        done = .false.
        
        !!Iteration loop
        do
            
            if ( done .eqv. .true. ) then
                exit
            endif
        
            !! Keep the magnetization from the current iteration and carry it to the next as the old value
            Mnorm_old = Mnorm
        
            !! If this is the first iteration and we should resume the iteration (resumeIteration != 0) then don't find the magnetization
            if ( cnt .eq. 0 .AND. resumeIteration .ne. 0. ) then
                lambda = resumeIteration
                !! set the old M value to the initial values from the input tiles
                do i=1,n
                    Mnorm_old(i) = sqrt( sum( tiles(i)%M**2 ) )
                enddo
            
            else
                !! Loop over each prism and set their respective magnetization vectors    
                do i=1,n            
                    if ( tiles(i)%tileType .lt. 100 .AND. tiles(i)%includeInIteration .ne. 0 ) then
                        select case ( tiles(i)%magnetType )
                        case ( MagnetTypeHard )
                            !!@todo Why is the CylPiece handled in a special way?
                            if ( tiles(i)%tileType .eq. tileTypeCylPiece ) then
                                call getM_HardMagnet_n( tiles(i)%H_ave, tiles(i), tiles(i)%n_ave(1)*tiles(i)%n_ave(2)*tiles(i)%n_ave(3) )
                            else
                                call getM_HardMagnet_1( H(i,:), tiles(i) )
                            endif
                
                        case ( MagnetTypeSoft )
                            call getM_SoftMagnet( T, H(i,:), tiles(i), stateFunction, n_stf)
                            
                        case ( MagnetTypeSoftConstPerm )
                            call getM_SoftConstMur( tiles(i), H(i,:) )
                            
                        case default
                        end select
            
                        tiles(i)%isIterating = .true.
                    endif
                enddo
                
                    
                !! set the Mnorm array
                do i=1,n
                    Mnorm(i) = sqrt( sum( tiles(i)%M**2 ) )
                enddo
                !! Ensure that the "old" magnetization is not zero initially
                if ( cnt .eq. 0 ) then
                    Mnorm_old = Mnorm
                endif
            
            endif        
            
            
            H_old = H   !< Reset H array
            H(:,:) = 0  !< Make sure to reset the H field
        
            !! Get the field at the center at each tile from all other tiles        
            do i=1,n
                if ( tiles(i)%includeInIteration .ne. 0 ) then
                    !! Consider the different possible tiles
                    select case (tiles(i)%tileType )
                    case (tileTypeCylPiece)                
                        !! Find the contribution to the internal field of the i'th tile from all tiles           
                        pts(1) = cos(tiles(i)%theta0) * tiles(i)%r0 + tiles(i)%offset(1)
                        pts(2) = sin(tiles(i)%theta0) * tiles(i)%r0 + tiles(i)%offset(2)
                        pts(3) = tiles(i)%z0 + tiles(i)%offset(3)
                        if ( tiles(i)%fieldEvaluation == fieldEvaluationCentre ) then
                            !! Do nothing as the field is already know at the tile centers
                    
                        elseif ( tiles(i)%fieldEvaluation == fieldEvaluationAverage ) then                    
                            !::Get the field in the pre-defined points. This is used to find the average magnetization (a few lines up in the call to getM_HardMagnet_n)
                            tiles(i)%H_ave(:,:) = 0.
                            call getFieldFromTiles( tiles, tiles(i)%H_ave(:,:), tiles(i)%H_ave_pts, n, tiles(i)%n_ave(1) * tiles(i)%n_ave(2) * tiles(i)%n_ave(3), tiles(i)%N_ave_pts )
                    
                        endif
                    
                    case(tileTypePrism)
                        !! No rotation is needed as the offset of the prism is with respect to the center of the prism
                        pts = tiles(i)%offset      
                        
                    case(tileTypeSphere)
                        !! No rotation is possible, as the tile is a sphere
                        pts = tiles(i)%offset 
            
                    case(tileTypeSpheroid)
                        !! No rotation is needed
                        pts = tiles(i)%offset 
                        
                    case(tileTypeCircPiece)
                        !! Find the center point of the circ piece
                        !! First find the center point of the piece in the reference frame of the piece itself
                        pts(1) = 0.5 * ( tiles(i)%r0 + tiles(i)%dr/2 ) * ( cos( tiles(i)%theta0 + tiles(i)%dtheta/2 ) + cos( tiles(i)%theta0 - tiles(i)%dtheta/2 ) )
                        pts(2) = 0.5 * ( tiles(i)%r0 + tiles(i)%dr/2 ) * ( sin( tiles(i)%theta0 + tiles(i)%dtheta/2 ) + sin( tiles(i)%theta0 - tiles(i)%dtheta/2 ) )
                        pts(3) = tiles(i)%z0
                
                        !! Then rotate this point about the three principal axes
                        !! Note that the rotMat gives the rotation from the global coordinate system to the reference system of the tile
                        !! The rotMatInv matrix gives the rotation from the tile's local coordinate system to the global one and is just the inverse rotation of the former.
                        call getRotationMatrices( tiles(i), rotMat, rotMatInv)
                        
                        pts = matmul( rotMatInv, pts )  !< We require the point in the global coordinate system.
                        
                        pts = pts + tiles(i)%offset     !< Finally add the offset to translate the circ piece
                                
            
                    case(tileTypeCircPieceInverted)
                        !! Find the center point of the circ piece
                        !! First find the center point of the piece in the reference frame of the piece itself
                        !pts(1) = 0.5 * ( tiles(i)%r0 + tiles(i)%dr/2 ) * ( cos( tiles(i)%theta0 + tiles(i)%dtheta/2 ) + cos( tiles(i)%theta0 - tiles(i)%dtheta/2 ) )
                        !pts(2) = 0.5 * ( tiles(i)%r0 + tiles(i)%dr/2 ) * ( sin( tiles(i)%theta0 + tiles(i)%dtheta/2 ) + sin( tiles(i)%theta0 - tiles(i)%dtheta/2 ) )
                        !!@Why the factor of 1.01?
                        pts(1) = 1.01*( tiles(i)%r0 + tiles(i)%dr/2 ) * cos( tiles(i)%theta0 )
                        pts(2) = 1.01*( tiles(i)%r0 + tiles(i)%dr/2 ) * sin( tiles(i)%theta0 )
                        pts(3) = tiles(i)%z0
                
                        !! Then rotate this point about the three principal axes
                        !! Note that the rotMat gives the rotation from the global coordinate system to the reference system of the tile
                        !! The rotMatInv matrix gives the rotation from the tile's local coordinate system to the global one and is just the inverse rotation of the former.
                        call getRotationMatrices( tiles(i), rotMat, rotMatInv)
                        
                        pts = matmul( rotMatInv, pts )  !< We require the point in the global coordinate system.
                        
                        pts = pts + tiles(i)%offset     !< Finally add the offset to translate the circ piece
                            
            
                    case default
                
                    end select
                    !! Tiles with type >100 are special tiles like a coil that are not iterated over
                    if ( tiles(i)%tileType .lt. 100 ) then
                        
                        !>If the magnet is assumed to have a constant permeability, then its own field is not included in this summation but rather calculated explicitly (see a few lines down)
                        if ( tiles(i)%magnetType .eq. magnetTypeSoftConstPerm ) then
                            tiles(i)%excludeFromSummation = .true.
                        endif
                        
                        call getFieldFromTiles( tiles, H(i,:), pts, n, 1, Nstore(i)%N )     !< Get the field in the i'th tile from all tiles           
                    
                        
                        !! When lambda == 1 then the new solution dominates. As lambda is decreased, the step towards the new solution is dampened
                        H(i,:) = H_old(i,:) + lambda * ( H(i,:) - H_old(i,:) )
                        
                        !>Set the flag back to false such that the tile is included in the next round of summation
                        if ( tiles(i)%magnetType .eq. magnetTypeSoftConstPerm ) then
                            tiles(i)%excludeFromSummation = .false.
                        endif
                    endif
                endif            
            enddo
        
            
            cnt = cnt + 1   !< Update iteration step
            
            !call saveMagnetizationDebug( tiles, n, cnt )   !< Only for debugging purposes, saves the magnetization out to a binary file
            
            if ( cnt .gt. 2 ) then
                !!@todo Remove?
                !Hnorm = sqrt( H(:,1)**2 + H(:,2)**2 + H(:,3)**2 )
                !Hnorm_old = sqrt( H_old(:,1)**2 + H_old(:,2)**2 + H_old(:,3)**2 )
            
                err_val = 0.
                !!@todo Remove?
                !where ( Hnorm_old .ne. 0. )
                !    err_val = abs( (Hnorm-Hnorm_old)/Hnorm_old )
                !endwhere
                where ( Mnorm_old .ne. 0. )
                    err_val = abs( (Mnorm-Mnorm_old)/Mnorm_old )
                endwhere
            
                !! Update the errorval in the tiles so that this can be used for re-iterating the model
                do i=1,n
                    tiles(i)%Mrel = err_val(i)
                enddo
            
                err = maxval(err_val)
            
                !! Update maxRelDiff array            
                maxRelDiffArr = cshift( maxRelDiffArr, 1 )
                maxRelDiffArr(1) = err
            
                if ( err .lt. err_max * lambda ) then
                    done = .true.

                else if ( cnt .gt. max_ite ) then     
                    done = .true.
                    i_err = -1
                endif    
            endif   
            
            !! Update lambda
            if ( cnt .gt. 4 .AND. lambdaCnt .gt. 4 ) then            
                call updateLambda( lambda, maxRelDiffArr, lCh )
                if ( lCh .eqv. .true. ) lambdaCnt = 1
            endif
        
            lambdaCnt = lambdaCnt + 1

             !call displayIteration( err, err_max * lambda )
             tmp = dispIteFct( err, err_max * lambda )
        
        enddo    
        deallocate(H,H_old,Hnorm,Hnorm_old,Nstore,err_val,Mnorm,Mnorm_old)
    end subroutine IterateMagnetization
    

    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !<
    !! Reduces lambda and sets lCh to true, if maxDiff oscillates.
    subroutine updateLambda( lambda, maxDiff, lCh )
        real,intent(inout) :: lambda
        real,intent(in),dimension(4) :: maxDiff
        logical,intent(inout) :: lCh

        lCh = .false.

        if ( (maxDiff(1) .gt. maxDiff(2) .AND. maxDiff(2) .lt. maxDiff(3) .AND. maxDiff(3) .gt. maxDiff(4)) .OR. &
             maxDiff(4) .gt. maxDiff(3) .AND. maxDiff(3) .gt. maxDiff(2) .AND. maxDiff(2) .gt. maxDiff(1) .OR. &
             maxDiff(4) .gt. maxDiff(3) .AND. maxDiff(3) .gt. maxDiff(2) .AND. maxDiff(2) .lt. maxDiff(1) ) then
            lambda = lambda * 0.5
            lCh = .true.
        endif
    end subroutine updateLambda
    
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !<
    !! Updates the magnetization vector of a tile based on the total magnetic field at it's center for a single tile
    !!@param H: a vector with three elements containing the magnetic field
    !!@param tile: type MagTile and is the tile under consideration
    !!
    !!@todo Why does this function exist in a _1 and an _n version? The _1 version seems only to be called if tiles(i)%tileType .eq. tileTypeCylPiece.
    subroutine getM_HardMagnet_1( H, tile )
        real,dimension(3), intent(in) :: H
        type(MagTile), intent(inout) :: tile

        real,dimension(3) :: M,un_x,un_y,un_z
        real :: H_par,H_trans_1,H_trans_2
    
        un_x(1) = 1
        un_x(2:3) = 0
    
        un_y(1) = 0
        un_y(2) = 1
        un_y(3) = 0
    
        un_z(1:2) = 0
        un_z(3) = 1
    
        !! Find the M vector in the i'th tile given the total internal field at the i'th tile's center
        !! Magnetization in the orthonormal basis of the easy axis
        M = 0
        
        H_par = dot_product( H, tile%u_ea )             !< Magnetic field along the easy axis
        M(1) = tile%Mrem + (tile%mu_r_ea - 1) * H_par   !< Magnetization along the easy axis   
        
        H_trans_1 = dot_product( H, tile%u_oa1 )        !< Magnetic field orthogonal to the easy axis
        M(2) = (tile%mu_r_oa - 1) * H_trans_1           !< Magnetization orthogonal to the easy axis       
                
        H_trans_2 = dot_product( H, tile%u_oa2 )        !< The other magnetic field component orthogonal to the easy axis
        M(3) = (tile%mu_r_oa - 1) * H_trans_2           !< Magnetization orthogonal to the easy axis       
        
        !! Magnetization of the i'th tile (defined with respect to the global coordinate system)
        tile%M(1) = dot_product( (M(1) * tile%u_ea + M(2) * tile%u_oa1 + M(3) * tile%u_oa2), un_x )
        tile%M(2) = dot_product( (M(1) * tile%u_ea + M(2) * tile%u_oa1 + M(3) * tile%u_oa2), un_y )
        tile%M(3) = dot_product( (M(1) * tile%u_ea + M(2) * tile%u_oa1 + M(3) * tile%u_oa2), un_z ) 
    
    end subroutine getM_HardMagnet_1

    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !<
    !! Updates the magnetization vector of a tile based on the total magnetic field at it's center
    !!@param H: a vector with three elements containing the magnetic field at the n tiles
    !!@param tile: type MagTile and is the tile under consideration
    !!@param n: the number of field evaluations
    subroutine getM_HardMagnet_n( H, tile, n )
        real,dimension(n,3), intent(in) :: H
        type(MagTile), intent(inout) :: tile
        integer,intent(in) :: n

        real,dimension(3) :: M,un_x,un_y,un_z
        real :: H_par,H_trans_1,H_trans_2
        integer :: i
    
        un_x(1) = 1
        un_x(2:3) = 0
    
        un_y(1) = 0
        un_y(2) = 1
        un_y(3) = 0
    
        un_z(1:2) = 0
        un_z(3) = 1
    
        !Find the M vector in the i'th tile given the total internal field at the i'th tile's center
        !Magnetization in the orthonormal basis of the easy axis
        M = 0
    
        do i=1,n
            H_par = dot_product( H(i,:), tile%u_ea )                !< Magnetic field along the easy axis           
            M(1) = M(1) + tile%Mrem + (tile%mu_r_ea - 1) * H_par    !< Magnetization along the easy axis       
            
            H_trans_1 = dot_product( H(i,:), tile%u_oa1 )           !< Magnetic field orthogonal to the easy axis
            M(2) = M(2) + (tile%mu_r_oa - 1) * H_trans_1            !< Magnetization orthogonal to the easy axis       
            
            H_trans_2 = dot_product( H(i,:), tile%u_oa2 )           !< The other magnetic field component orthogonal to the easy axis
            M(3) = M(3) + (tile%mu_r_oa - 1) * H_trans_2            !< Magnetization orthogonal to the easy axis       
        
        enddo
        
        !! Get the average value
        M = M / n
        
        !! Magnetization of the i'th tile (defined with respect to the global coordinate system)
        tile%M(1) = dot_product( (M(1) * tile%u_ea + M(2) * tile%u_oa1 + M(3) * tile%u_oa2), un_x )
        tile%M(2) = dot_product( (M(1) * tile%u_ea + M(2) * tile%u_oa1 + M(3) * tile%u_oa2), un_y )
        tile%M(3) = dot_product( (M(1) * tile%u_ea + M(2) * tile%u_oa1 + M(3) * tile%u_oa2), un_z ) 
    
    end subroutine getM_HardMagnet_n


    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !<
    !! Derives the magnetization of a soft ferromagnet given the input field and the state function
    !!
    subroutine getM_SoftMagnet( T, H, tile, stateFunction, n_stf)
        real,intent(in) :: T
        real,dimension(3), intent(in) :: H
        type(MagTile), intent(inout) :: tile
        type(MagStateFunction),intent(in),dimension(n_stf) :: stateFunction
        integer,intent(in) :: n_stf
    
        integer :: index    
        real :: Hnorm,Mnorm
    
        Hnorm = sqrt( H(1)**2 + H(2)**2 + H(3)**2 )
    
        index = tile%stateFunctionIndex
        !call splint( stateFunction(index)%H, stateFunction(index)%M(1,:), stateFunction(index)%y2a, Hnorm, Mnorm, stateFunction(index)%nH )
         call spline_b_val ( stateFunction(index)%nH, stateFunction(index)%H, stateFunction(index)%M(1,:), Hnorm, Mnorm )
         if ( Hnorm .ne. 0 ) then
             tile%M = Mnorm * H / Hnorm  !< Assume that the M vector is along the H vector
         else
             tile%M(:) = 0
         endif
     
    end subroutine getM_SoftMagnet
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !< Function called to get the magntization of a soft ferromagnet with constant permeability (isotropic)
    !! @param tile the tile of which the magnetization is wanted
    !! @param H the magnetic field at the tile
    !!
    subroutine getM_SoftConstMur( tile, H )
        type(MagTile),intent(inout) :: tile
        real,dimension(3),intent(in) :: H
        
        !Assuming B = mur * mu0 * (H + M) = mur * mu0 * H => M = H(mur-1)
        
        tile%M = H * ( tile%mu_r_ea - 1. )
    end subroutine getM_SoftConstMur
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !< Function called to dump each iteration's magnetization vector for each tile for debugging purposes
    !! @param tiles is an array of tiles who's magnetization vectors are to be saved
    !! @param n the number of tiles
    !! @param snaps the snapshot ID
    !!
    subroutine saveMagnetizationDebug( tiles, n, snaps )
        type(MagTile),dimension(n),intent(in) :: tiles
        integer,intent(in) :: n,snaps
    
        real,dimension(:,:),allocatable :: M
        integer :: i
    
        allocate(M(n,3))
    
        M(:,:) = 0.
    
        do i=1,n
            M(i,:) = tiles(i)%M    
        enddo
    
        open (14, file='M_debug.dat',	&
			               status='unknown', form='unformatted',action='write',	&
			               access='direct', recl=2*3*n)
        write(14,rec=snaps) M
    
        close(14)

        deallocate(M)
    
    end subroutine saveMagnetizationDebug

    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !< Function to load tiles from python module to the corresponding one fortran
    !! @param tileType - defines whether the tile is cylindrical, a prism, an ellipsoid and so on
    !! @param offset - the centre coordinates relative to the global coordinate system
    !! @param rotAngles - rotation angles (phi_x, phi_y, phi_z) about the principle axes of the tile with respect to the centre of the tile
    !! @param color - color rgb triplet
    !! @param magnetType - defines whether the tile is a hard or soft magnet
    !! @param stateFunctionIndex - index matching an entry into an array of type MagStateFunction. Used by soft ferromagnets (when interpolation on an M vs H curve is necessary)
    !! @param symmetryOps - 1 for symmetry, -1 for anti-symmetry ((1) for about xy plane, (2) for about (xz) plane and (3) for about yz plane)
    !! @param Mrel - the current relative change of the magnetization (i.e. abs((M1-M2)/M2 ) where M1 is the latest magnetization norm and M2 is the previous one
    !!
    subroutine loadTiles( centerPos, dev_center, tile_size, vertices, Mag, u_ea, u_oa1, u_oa2, &
        mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
        includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, tiles)

        !::Specific for a cylindrical tile piece
        real(8),dimension(n_tiles,3),intent(in) :: centerPos
        real(8),dimension(n_tiles,3),intent(in) :: dev_center
            
        !::Specific for a rectangular prism
        real(8),dimension(n_tiles,3),intent(in) :: tile_size

        !::Specific for a tetrahedron
        real(8),dimension(n_tiles,3,4),intent(in) :: vertices
            
        !::Generel variables, shared among all tile types
        real(8),dimension(n_tiles,3),intent(in) :: Mag
        real(8),dimension(n_tiles,3),intent(in) :: u_ea,u_oa1,u_oa2    
        real(8),dimension(n_tiles),intent(in) :: mu_r_ea,mu_r_oa,Mrem
        integer(4),dimension(n_tiles),intent(in) :: tileType
        real(8),dimension(n_tiles,3),intent(in) :: offset
        real(8),dimension(n_tiles,3),intent(in) :: rotAngles
        real(8),dimension(n_tiles,3),intent(in) :: color
        integer(4),dimension(n_tiles),intent(in) :: magnetType
        integer(4),dimension(n_tiles),intent(in) :: stateFunctionIndex
        integer(4),dimension(n_tiles),intent(in) :: includeInIteration,exploitSymmetry
        real(8),dimension(n_tiles,3),intent(in) :: symmetryOps
        real(8),dimension(n_tiles),intent(in) :: Mrel

        type(MagTile),dimension(n_tiles),intent(inout) :: tiles
        integer(4),intent(in) :: n_tiles
        integer :: i

        do i=1,n_tiles
            tiles(i)%r0 = centerPos(i,1)
            tiles(i)%theta0 = centerPos(i,2)
            tiles(i)%z0 = centerPos(i,3)
            tiles(i)%dr = dev_center(i,1)
            tiles(i)%dtheta = dev_center(i,2)
            tiles(i)%dz = dev_center(i,3)
            tiles(i)%a = tile_size(i,1)
            tiles(i)%b = tile_size(i,2)
            tiles(i)%c = tile_size(i,3)
            tiles(i)%vert(:,1) = vertices(i,:,1)
            tiles(i)%vert(:,2) = vertices(i,:,2)
            tiles(i)%vert(:,3) = vertices(i,:,3)
            tiles(i)%vert(:,4) = vertices(i,:,4)
            tiles(i)%M = Mag(i,:)
            tiles(i)%u_ea = u_ea(i,:)
            tiles(i)%u_oa1 = u_oa1(i,:)
            tiles(i)%u_oa2 = u_oa2(i,:)
            tiles(i)%mu_r_ea = mu_r_ea(i)
            tiles(i)%mu_r_oa = mu_r_oa(i)
            tiles(i)%Mrem = Mrem(i)
            tiles(i)%tileType = tileType(i)
            tiles(i)%offset = offset(i,:)
            tiles(i)%rotAngles = rotAngles(i,:)
            tiles(i)%color = color(i,:)
            tiles(i)%magnetType = magnetType(i)
            tiles(i)%stateFunctionIndex = stateFunctionIndex(i)
            tiles(i)%includeInIteration = includeInIteration(i)
            tiles(i)%exploitSymmetry = exploitSymmetry(i)
            tiles(i)%symmetryOps = symmetryOps(i,:)
            tiles(i)%Mrel = Mrel(i)

            if ( tiles(i)%tileType ==  tileTypeCylPiece ) then
                call setupEvaluationPoints( tiles(i) )
            endif

            !> Set the default behavior to false such that the tile by default is included in the calculation
            tiles(i)%excludeFromSummation = .false.

        enddo

    end subroutine loadTiles

    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !< Function to load the state function for each tile from python to fortran
    !!
    subroutine loadStateFunction( nT, nH, stateFcn, data_stateFcn, n_stateFcn )
        
        integer,intent(in) :: nT,nH
        type(MagStateFunction),dimension(n_stateFcn),intent(inout) :: stateFcn
        real,dimension(nH,nT),intent(in) :: data_stateFcn
        integer,intent(in) :: n_stateFcn
        integer :: i, j, k

        do i=1,n_stateFcn
            stateFcn(i)%nT = nT-1
            stateFcn(i)%nH = nH-1
            allocate( stateFcn(i)%M(stateFcn(i)%nT,stateFcn(i)%nH) )
            allocate( stateFcn(i)%T(stateFcn(i)%nT) )
            allocate( stateFcn(i)%H(stateFcn(i)%nH) )
            allocate( stateFcn(i)%y2a(stateFcn(i)%nH) )
        
            stateFcn(i)%T(1) = data_stateFcn(1,2)
            stateFcn(i)%T(2) = data_stateFcn(1,3)
            stateFcn(i)%T(3) = data_stateFcn(1,4)

            do j=2,nH
                k = j - 1
                stateFcn(i)%H(k) = data_stateFcn(j,1)
                stateFcn(i)%M(1,k) = data_stateFcn(j,2)
                stateFcn(i)%M(2,k) = data_stateFcn(j,3)
                stateFcn(i)%M(3,k) = data_stateFcn(j,4)
                
                !! make the spline derivatives for later interpolation
                !call splie2( sngl(stateFunction(i)%T), sngl(stateFunction(i)%H), sngl(stateFunction(i)%M), sngl(stateFunction(i)%nT), sngl(stateFunction(i)%nH), sngl(stateFunction(i)%y2a) )
                !!@todo If this code is deprecated, the spline.f90 code can be removed as it is only called here.
                !call spline( stateFunction(i)%H, stateFunction(i)%M(1,:), 1e30, 1e30, stateFunction(i)%y2a, stateFunction(i)%nH )
            
            enddo
        enddo

    end subroutine loadStateFunction  
    
    end module IterateMagnetSolution
    
