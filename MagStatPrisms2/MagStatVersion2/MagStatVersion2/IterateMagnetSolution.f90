
    module IterateMagnetSolution

    use MagStat2GetSolution
    use MagStatParameters
    use MagStatUtil
    use MagTileIO
    implicit none
    
    
    contains
    
    
    
!::
!::Iterates over all tiles to find the self-consistent solution and returns the tiles with updated magnetization vectors
!::This subroutine is only relevant for permanent magnets that may be modeled with a permeability vector
!::
!::Error state is returned in i_err:
!::0 no error
!:: -1 max number of iterations reached without satisfying the minimum error
subroutine iterateMagnetization( tiles, n, stateFunction, n_stf, T, err_max )
    type(MagTile),dimension(n),intent(inout) :: tiles
    integer,intent(in) :: n
    type(MagStatStateFunction),intent(in),dimension(n_stf) :: stateFunction    
    integer,intent(in) :: n_stf
    real,intent(in) :: T
    real,optional :: err_max
    
    
    integer,parameter :: cnt_max = 1000
        
    logical :: done
    integer :: i,j,cnt,i_err,lambdaCnt
    real,dimension(3) :: M
    real,dimension(:,:),allocatable :: H,H_old
    real,dimension(:),allocatable :: Hnorm,Hnorm_old
    real,dimension(3) :: pts    
    real :: H_par,H_trans_1,H_trans_2,err,lambda    
    real,dimension(4) :: maxRelDiffArr
    logical :: lCh
    
    !::set default error margin if nothing has been specified by the user
    if ( .not. present(err_max) ) err_max = 1e-10
    
    i_err = 0
    
    allocate(H(n,3),H_old(n,3),Hnorm(n),Hnorm_old(n))
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
            
            select case ( tiles(i)%magnetType )
            case ( MagnetTypeHard )
                call getM_HardMagnet( H(i,:), tiles(i) )
            case ( MagnetTypeSoft )
                call getM_SoftMagnet( T, H(i,:), tiles(i), stateFunction, n_stf)
            case default
            end select
            
            
        enddo
        !reset H array
        H_old = H
        !make sure to reset the H field
        H(:,:) = 0
        !::Get the field at the center at each tile from all other tiles
        do i=1,n
            
            select case (tiles(i)%tileType )
            case (tileTypeCylPiece)
                !find the contribution to the internal field of the i'th tile from all tiles           
                pts(1) = cos(tiles(i)%theta0) * tiles(i)%r0
                pts(2) = sin(tiles(i)%theta0) * tiles(i)%r0
                pts(3) = tiles(i)%z0
            case(tileTypePrism)
                pts = tiles(i)%offset
            case default
                
            end select
            
           
           !::Get the field in the i'th tile from all tiles           
           call getFieldFromTiles( tiles, H(i,:), pts, n, 1 )           
           
           !subtract the magnetization of the i'th tile to get H
           
           !::When lambda == 1 then the new solution dominates. As lambda is decreased, the step towards the new solution is slowed
           !::If the tile is a cylinder piece then M should be subtracted due to the way the tensor was derived. If the tile is a prism them M should not be subtracted. Really idiotic (Kaspar, 21-12-2017)
           !select case (tiles(i)%tileType )
           !case (tileTypeCylPiece)
           !    H(i,:) = H_old(i,:) + lambda * ( H(i,:) - tiles(i)%M - H_old(i,:) )
           !case(tileTypePrism)
           !    H(i,:) = H_old(i,:) + lambda * ( H(i,:) - H_old(i,:) )
           !case default
           !     
           ! end select
           H(i,:) = H_old(i,:) + lambda * ( H(i,:) - H_old(i,:) )
           
        enddo
        
        !::Update iteration step
        cnt = cnt + 1
        if ( cnt .gt. 1 ) then
            Hnorm = sqrt( H(:,1)**2 + H(:,2)**2 + H(:,3)**2 )
            Hnorm_old = sqrt( H_old(:,1)**2 + H_old(:,2)**2 + H_old(:,3)**2 )
            err = maxval(abs( (Hnorm-Hnorm_old)/Hnorm_old ))
            
            !::Update maxRelDiff array            
            maxRelDiffArr = cshift( maxRelDiffArr, 1 )
            maxRelDiffArr(1) = err
            
            if ( err .lt. err_max * lambda ) then
                done = .true.
                
            else if ( cnt .gt. cnt_max ) then     
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
        call displayIteration( err, err_max * lambda )
                      
        
    enddo    
    deallocate(H,H_old,Hnorm,Hnorm_old)    
end subroutine IterateMagnetization
    

!:: Reduces lambda and sets lCh to true, if maxDiff oscillates.
subroutine updateLambda( lambda, maxDiff, lCh )
real,intent(inout) :: lambda
real,intent(in),dimension(4) :: maxDiff
logical,intent(inout) :: lCh

lCh = .false.

if ( (maxDiff(1) .gt. maxDiff(2) .AND. maxDiff(2) .lt. maxDiff(3) .AND. maxDiff(3) .gt. maxDiff(4)) .OR. &
     (maxDiff(1) .lt. maxDiff(2) .AND. maxDiff(2) .gt. maxDiff(3) .AND. maxDiff(3) .lt. maxDiff(4))   ) then

    lambda = lambda * 0.5
    lCh = .true.
endif
end subroutine updateLambda
!::
!::Updates the magnetization vector of a tile based on the total magnetic field at it's center
!::H is a vector with three elements containing the magnetic field
!::tile is of type MagTile and is the tile under consideration
!::
subroutine getM_HardMagnet( H, tile )
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
    
end subroutine getM_HardMagnet

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
    