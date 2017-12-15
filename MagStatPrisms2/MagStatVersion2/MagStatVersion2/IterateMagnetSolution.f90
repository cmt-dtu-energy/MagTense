module IterateMagnetSolution
use MagStat2GetSolution

    implicit none
    
    
    contains
    
    
    
!::
!::Iterates over all tiles to find the self-consistent solution and returns the tiles with updated magnetization vectors
!::This subroutine is only relevant for permanent magnets that may be modeled with a permeability vector
!::
!::Error state is returned in i_err:
!::0 no error
!:: -1 max number of iterations reached without satisfying the minimum error
subroutine iterateMagnetization( tiles, n, err_max )
    type(MagTile),dimension(n),intent(inout) :: tiles
    integer,intent(in) :: n    
    real,optional :: err_max
    
    integer,parameter :: cnt_max = 30
        
    logical :: done
    integer :: i,j,cnt,i_err
    real,dimension(3) :: M
    real,dimension(:,:),allocatable :: H,H_old
    real,dimension(:),allocatable :: Hnorm,Hnorm_old
    real,dimension(3) :: un_x,un_y,un_z,pts    
    real :: H_par,H_trans_1,H_trans_2,err
    
    !::set default error margin if nothing has been specified by the user
    if ( .not. present(err_max) ) err_max = 1e-10
    
    i_err = 0
    
    
    un_x(1) = 1
    un_x(2:3) = 0
    
    un_y(1) = 0
    un_y(2) = 1
    un_y(3) = 0
    
    un_z(1:2) = 0
    un_z(3) = 1
    
    allocate(H(n,3),H_old(n,3),Hnorm(n),Hnorm_old(n))
    H(:,:) = 0
    
    cnt = 0
    
    done = .false.
    
    do
        if ( done .eq. .true. ) then
            exit
        endif
        
        
        do i=1,n

           !Find the M vector in the i'th tile given the total internal field
           !at the i'th tile's center
           !Magnetization in the orthonormal basis of the easy axis
           M(:) = 0
           !Magnetic field along the easy axis
           H_par = dot_product( H(i,:), tiles(i)%u_ea )
           
           !Magnetization along the easy axis       
           M(1) = tiles(i)%Mrem + (tiles(i)%mu_r_ea - 1) * H_par

           !Magnetic field orthogonal to the easy axis
           H_trans_1 = dot_product( H(i,:), tiles(i)%u_oa1 )
           !Magnetization orthogonal to the easy axis       
           M(2) = (tiles(i)%mu_r_oa - 1) * H_trans_1

           !the other magnetic field component orthogonal to the easy axis
           H_trans_2 = dot_product( H(i,:), tiles(i)%u_oa2 )
           !Magnetization orthogonal to the easy axis       
           M(3) = (tiles(i)%mu_r_oa - 1) * H_trans_2
           !magnetization of the i'th tile (defined with respect to the global
           !coordinate system)
           tiles(i)%M(1) = dot_product( (M(1) * tiles(i)%u_ea + M(2) * tiles(i)%u_oa1 + M(3) * tiles(i)%u_oa2), un_x )
           tiles(i)%M(2) = dot_product( (M(1) * tiles(i)%u_ea + M(2) * tiles(i)%u_oa1 + M(3) * tiles(i)%u_oa2), un_y )
           tiles(i)%M(3) = dot_product( (M(1) * tiles(i)%u_ea + M(2) * tiles(i)%u_oa1 + M(3) * tiles(i)%u_oa2), un_z ) 
        enddo
        !reset H array
        H_old = H
        !make sure to reset the H field
        H(:,:) = 0
        
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
            H(i,:) = H(i,:) - tiles(i)%M
        enddo
        
        cnt = cnt + 1
        if ( cnt .gt. 1 ) then
            Hnorm = sqrt( H(:,1)**2 + H(:,2)**2 + H(:,3)**2 )
            Hnorm_old = sqrt( H_old(:,1)**2 + H_old(:,2)**2 + H_old(:,3)**2 )
            err = maxval(abs( (Hnorm-Hnorm_old)/Hnorm_old ))
            
            if ( err .lt. err_max ) then
                done = .true.
                
            else if ( cnt .gt. cnt_max ) then     
                done = .true.
                i_err = -1
            endif
        endif
           
    enddo    
    deallocate(H,H_old,Hnorm,Hnorm_old)    
    end subroutine IterateMagnetization
    
    
end module IterateMagnetSolution
    