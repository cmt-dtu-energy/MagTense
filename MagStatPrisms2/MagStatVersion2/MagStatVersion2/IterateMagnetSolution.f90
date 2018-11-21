
    module IterateMagnetSolution

    use MagStat2GetSolution
    use MagStatParameters
    use MagStatUtil
    !use MagTileIO
    implicit none
    
    contains
    
    
    
!::
!::Iterates over all tiles to find the self-consistent solution and returns the tiles with updated magnetization vectors
!::This subroutine is only relevant for permanent magnets that may be modeled with a permeability vector
!::
!::Error state is returned in i_err:
!::0 no error
!:: -1 max number of iterations reached without satisfying the minimum error
subroutine iterateMagnetization( tiles, n, stateFunction, n_stf, T, err_max, max_ite, displayIteration_fct )
    type(MagTile),dimension(n),intent(inout) :: tiles
    integer,intent(in) :: n
    type(MagStatStateFunction),intent(in),dimension(n_stf) :: stateFunction    
    integer,intent(in) :: n_stf

    real,intent(in) :: T    
    real,optional :: err_max
    INTERFACE
		FUNCTION displayIteration_fct(err,err_max)		
		REAL, INTENT(IN) :: err,err_max		
		END FUNCTION displayIteration_fct
	END INTERFACE
    
    
    integer,parameter :: cnt_max = 200
        
    logical :: done
    integer :: i,j,cnt,i_err,lambdaCnt
    real,dimension(3) :: M
    real,dimension(:,:),allocatable :: H,H_old
    real,dimension(:),allocatable :: Hnorm,Hnorm_old,err_val
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
    
    allocate(H(n,3),H_old(n,3),Hnorm(n),Hnorm_old(n),err_val(n))
    allocate(Nstore(n))
    H(:,:) = 0
    maxRelDiffArr(:) = 0.
    cnt = 0
    
    !::Relaxation parameter
    lambda = 1.
    lambdaCnt = 0
    
    done = .false.
    !::Iteration loop
    do
        if ( done .eq. .true. ) then
            exit
        endif
        
        !::Loop over each prism and set their respective magnetization vectors
        do i=1,n
            if ( tiles(i)%tileType .lt. 100 .AND. tiles(i)%includeInIteration .ne. 0 ) then
                select case ( tiles(i)%magnetType )
                case ( MagnetTypeHard )
                    if ( tiles(i)%tileType .eq. tileTypeCylPiece ) then
                        call getM_HardMagnet_n( tiles(i)%H_ave, tiles(i), tiles(i).n_ave(1)*tiles(i).n_ave(2)*tiles(i).n_ave(3) )
                    else
                        call getM_HardMagnet_1( H(i,:), tiles(i) )
                    endif
                
                case ( MagnetTypeSoft )
                    call getM_SoftMagnet( T, H(i,:), tiles(i), stateFunction, n_stf)
                case default
                end select
            
                tiles(i)%isIterating = .true.
            endif            
        enddo
        !reset H array
        H_old = H
        !make sure to reset the H field
        H(:,:) = 0
        !::Get the field at the center at each tile from all other tiles
        
        do i=1,n
            if ( tiles(i)%includeInIteration .ne. 0 ) then
                select case (tiles(i)%tileType )
                case (tileTypeCylPiece)                
                    !find the contribution to the internal field of the i'th tile from all tiles           
                    pts(1) = cos(tiles(i)%theta0) * tiles(i)%r0 + tiles(i)%offset(1)
                    pts(2) = sin(tiles(i)%theta0) * tiles(i)%r0 + tiles(i)%offset(2)
                    pts(3) = tiles(i)%z0 + tiles(i)%offset(3)
                    if ( tiles(i)%fieldEvaluation == fieldEvaluationCentre ) then
                    
            
                    
                    elseif ( tiles(i)%fieldEvaluation == fieldEvaluationAverage ) then                    
                        !::Get the field in the pre-defined points. This is used to find the average magnetization (a few lines up in the call to getM_HardMagnet_n)
                        tiles(i)%H_ave(:,:) = 0.
                        call getFieldFromTiles( tiles, tiles(i)%H_ave(:,:), tiles(i)%H_ave_pts, n, tiles(i)%n_ave(1) * tiles(i)%n_ave(2) * tiles(i)%n_ave(3), tiles(i)%N_ave_pts )
                    
                    endif
                    
                case(tileTypePrism)
                    !::no rotation is needed as the offset of the prism is with respect to the center of the prism
                    pts = tiles(i)%offset                
            
                case(tileTypeCircPiece)
                    !find the center point of the circ piece
                    !::First find the center point of the piece in the reference frame of the piece itself
                    pts(1) = 0.5 * ( tiles(i)%r0 + tiles(i)%dr/2 ) * ( cos( tiles(i)%theta0 + tiles(i)%dtheta/2 ) + cos( tiles(i)%theta0 - tiles(i)%dtheta/2 ) )
                    pts(2) = 0.5 * ( tiles(i)%r0 + tiles(i)%dr/2 ) * ( sin( tiles(i)%theta0 + tiles(i)%dtheta/2 ) + sin( tiles(i)%theta0 - tiles(i)%dtheta/2 ) )
                    pts(3) = tiles(i)%z0
                
                    !::Then rotate this point about the three principal axes
                    !::Note that the rotMat gives the rotation from the global coordinate system to the 
                    !::reference system of the tile
                    !::The rotMatInv matrix gives the rotation from the tile's local coordinate system to the global one
                    !::and is just the inverse rotation of the former.
                    call getRotationMatrices( tiles(i), rotMat, rotMatInv)
                    !::We require the point in the global coordinate system.
                    pts = matmul( rotMatInv, pts )
                    !::Finally add the offset to translate the circ piece
                    pts = pts + tiles(i)%offset
                                
            
                case(tileTypeCircPieceInverted)
                    !find the center point of the circ piece
                    !::First find the center point of the piece in the reference frame of the piece itself
                    !pts(1) = 0.5 * ( tiles(i)%r0 + tiles(i)%dr/2 ) * ( cos( tiles(i)%theta0 + tiles(i)%dtheta/2 ) + cos( tiles(i)%theta0 - tiles(i)%dtheta/2 ) )
                    !pts(2) = 0.5 * ( tiles(i)%r0 + tiles(i)%dr/2 ) * ( sin( tiles(i)%theta0 + tiles(i)%dtheta/2 ) + sin( tiles(i)%theta0 - tiles(i)%dtheta/2 ) )
                    pts(1) = 1.01*( tiles(i)%r0 + tiles(i)%dr/2 ) * cos( tiles(i)%theta0 )
                    pts(2) = 1.01*( tiles(i)%r0 + tiles(i)%dr/2 ) * sin( tiles(i)%theta0 )
                    pts(3) = tiles(i)%z0
                
                    !::Then rotate this point about the three principal axes
                    !::Note that the rotMat gives the rotation from the global coordinate system to the 
                    !::reference system of the tile
                    !::The rotMatInv matrix gives the rotation from the tile's local coordinate system to the global one
                    !::and is just the inverse rotation of the former.
                    call getRotationMatrices( tiles(i), rotMat, rotMatInv)
                    !::We require the point in the global coordinate system.
                    pts = matmul( rotMatInv, pts )
                    !::Finally add the offset to translate the circ piece
                    pts = pts + tiles(i)%offset
                
                            
            
                case default
                
                end select
                !::Tiles with type >100 are special tiles like a coil that are not iterated over
                if ( tiles(i)%tileType .lt. 100 ) then
                    !::Get the field in the i'th tile from all tiles           
                    
                    call getFieldFromTiles( tiles, H(i,:), pts, n, 1, Nstore(i)%N )
                    
                   
                    !::When lambda == 1 then the new solution dominates. As lambda is decreased, the step towards the new solution is slowed
           
                    H(i,:) = H_old(i,:) + lambda * ( H(i,:) - H_old(i,:) )
                endif
            endif            
        enddo
        
        !::Update iteration step
        cnt = cnt + 1
        if ( cnt .gt. 1 ) then
            Hnorm = sqrt( H(:,1)**2 + H(:,2)**2 + H(:,3)**2 )
            Hnorm_old = sqrt( H_old(:,1)**2 + H_old(:,2)**2 + H_old(:,3)**2 )
            err_val = 0.
            where ( Hnorm_old .ne. 0. )
                err_val = abs( (Hnorm-Hnorm_old)/Hnorm_old )
            endwhere
            err = maxval(err_val)
            
            
            !::Update maxRelDiff array            
            maxRelDiffArr = cshift( maxRelDiffArr, 1 )
            maxRelDiffArr(1) = err
            
            if ( err .lt. err_max * lambda ) then
                done = .true.
                

            else if ( cnt .gt. max_ite ) then     
                done = .true.
                i_err = -1
            endif    
        endif   
        !::update lambda
        if ( cnt .gt. 4 .AND. lambdaCnt .gt. 4 ) then            
            call updateLambda( lambda, maxRelDiffArr, lCh )
            if ( lCh .eq. .true. ) lambdaCnt = 1
        endif
        
        lambdaCnt = lambdaCnt + 1

         !call displayIteration( err, err_max * lambda )
         tmp = displayIteration_fct( err, err_max * lambda )
        
    enddo    
    deallocate(H,H_old,Hnorm,Hnorm_old,Nstore,err_val)
end subroutine IterateMagnetization
    

!:: Reduces lambda and sets lCh to true, if maxDiff oscillates.
subroutine updateLambda( lambda, maxDiff, lCh )
real,intent(inout) :: lambda
real,intent(in),dimension(4) :: maxDiff
logical,intent(inout) :: lCh

lCh = .false.

if ( (maxDiff(1) .gt. maxDiff(2) .AND. maxDiff(2) .lt. maxDiff(3) .AND. maxDiff(3) .gt. maxDiff(4)) .OR. &

    lambda = lambda * 0.5
    lCh = .true.
endif
end subroutine updateLambda
!::
!::Updates the magnetization vector of a tile based on the total magnetic field at it's center
!::H is a vector with three elements containing the magnetic field
!::tile is of type MagTile and is the tile under consideration
!::
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
    
    !Find the M vector in the i'th tile given the total internal field
    !at the i'th tile's center
    !Magnetization in the orthonormal basis of the easy axis
    M = 0
    !Magnetic field along the easy axis
    H_par = dot_product( H, tile%u_ea )
           
    !Magnetization along the easy axis       
    M(1) = tile%Mrem + (tile%mu_r_ea - 1) * H_par

    !Magnetic field orthogonal to the easy axis
    H_trans_1 = dot_product( H, tile%u_oa1 )
    !Magnetization orthogonal to the easy axis       
    M(2) = (tile%mu_r_oa - 1) * H_trans_1

    !the other magnetic field component orthogonal to the easy axis
    H_trans_2 = dot_product( H, tile%u_oa2 )
    !Magnetization orthogonal to the easy axis       
    M(3) = (tile%mu_r_oa - 1) * H_trans_2
    !magnetization of the i'th tile (defined with respect to the global
    !coordinate system)
    tile%M(1) = dot_product( (M(1) * tile%u_ea + M(2) * tile%u_oa1 + M(3) * tile%u_oa2), un_x )
    tile%M(2) = dot_product( (M(1) * tile%u_ea + M(2) * tile%u_oa1 + M(3) * tile%u_oa2), un_y )
    tile%M(3) = dot_product( (M(1) * tile%u_ea + M(2) * tile%u_oa1 + M(3) * tile%u_oa2), un_z ) 
    
end subroutine getM_HardMagnet_1

!::
!::n is the number of field evaluations
!::
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
    
    !Find the M vector in the i'th tile given the total internal field
    !at the i'th tile's center
    !Magnetization in the orthonormal basis of the easy axis
    M = 0
    
    do i=1,n
    
        !Magnetic field along the easy axis
        H_par = dot_product( H(i,:), tile%u_ea )
           
        !Magnetization along the easy axis       
        M(1) = M(1) + tile%Mrem + (tile%mu_r_ea - 1) * H_par

        !Magnetic field orthogonal to the easy axis
        H_trans_1 = dot_product( H(i,:), tile%u_oa1 )
        !Magnetization orthogonal to the easy axis       
        M(2) = M(2) + (tile%mu_r_oa - 1) * H_trans_1

        !the other magnetic field component orthogonal to the easy axis
        H_trans_2 = dot_product( H(i,:), tile%u_oa2 )
        !Magnetization orthogonal to the easy axis       
        M(3) = M(3) + (tile%mu_r_oa - 1) * H_trans_2
        
    enddo
    !::get the average value
    M = M / n
    !magnetization of the i'th tile (defined with respect to the global
    !coordinate system)
    tile%M(1) = dot_product( (M(1) * tile%u_ea + M(2) * tile%u_oa1 + M(3) * tile%u_oa2), un_x )
    tile%M(2) = dot_product( (M(1) * tile%u_ea + M(2) * tile%u_oa1 + M(3) * tile%u_oa2), un_y )
    tile%M(3) = dot_product( (M(1) * tile%u_ea + M(2) * tile%u_oa1 + M(3) * tile%u_oa2), un_z ) 
    
    
end subroutine getM_HardMagnet_n


!::
!::Derives the magnetization of a soft ferromagnet given the input field and the state function
!::
subroutine getM_SoftMagnet( T, H, tile, stateFunction, n_stf)
    real,intent(in) :: T
    real,dimension(3), intent(in) :: H
    type(MagTile), intent(inout) :: tile
    type(MagStatStateFunction),intent(in),dimension(n_stf) :: stateFunction
    integer,intent(in) :: n_stf
    
    real :: Hnorm,Mnorm
    
    Hnorm = sqrt( H(1)**2 + H(2)**2 + H(3)**2 )
    
    call getBilinInterp( stateFunction(tile%stateFunctionIndex)%M, stateFunction(tile%stateFunctionIndex)%T, &
                     stateFunction(tile%stateFunctionIndex)%H, stateFunction(tile%stateFunctionIndex)%nT, &
                     stateFunction(tile%stateFunctionIndex)%nH, T, Hnorm, Mnorm )   
     
     if ( Hnorm .ne. 0 ) then
        !::assumes that the M vector is along the H vector
        tile%M = Mnorm * H / Hnorm
     else
         tile%M(:) = 0
     endif
     

end subroutine getM_SoftMagnet
    
end module IterateMagnetSolution
    