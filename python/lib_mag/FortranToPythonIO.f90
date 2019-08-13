module FortranToPythonIO

use DemagFieldGetSolution
use IterateMagnetSolution
implicit none

contains

!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!>
!> function for displaying output (to the terminal)
!! @param err the current relative error
!! @param err_max the current maximum allowed error
function dispIte_fct( err, err_max )
    real,intent(in) :: err,err_max
    integer :: dispIte_fct
   
    write(*,*) err,err_max
  
    dispIte_fct = 0
    
end function dispIte_fct


subroutine getNFromTiles( centerPos, dev_center, rect_size, Mag, u_ea, u_oa1, u_oa2, mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
    includeInIteration, exploitSymmetry, symmetryOps, Mrel, pts, n_tiles, n_pts, N )
    !::Specific for a cylindrical tile piece
    real(8),dimension(n_tiles,3),intent(in) :: centerPos
    real(8),dimension(n_tiles,3),intent(in) :: dev_center
        
    !::Specific for a rectangular prism
    real(8),dimension(n_tiles,3),intent(in) :: rect_size
        
    !::Generel variables, shared among all tile types
    real(8),dimension(n_tiles,3),intent(in) :: Mag
    real(8),dimension(n_tiles,3),intent(in) :: u_ea,u_oa1,u_oa2    
    real(8),dimension(n_tiles),intent(in) :: mu_r_ea,mu_r_oa,Mrem
    integer(4),dimension(n_tiles),intent(in) :: tileType        !::defines whether the tile is cylindrical, a prism, an ellipsoid and so on
    real(8),dimension(n_tiles,3),intent(in) :: offset !::the centre coordinates relative to the global coordinate system
    real(8),dimension(n_tiles,3),intent(in) :: rotAngles !:: rotation angles (phi_x, phi_y, phi_z) about the principle axes of the tile with respect to the centre of the tile
    real(8),dimension(n_tiles,3),intent(in) :: color !! color rgb triplet
    integer(4),dimension(n_tiles),intent(in) :: magnetType !::defines whether the tile is a hard or soft magnet
    integer(4),dimension(n_tiles),intent(in) :: stateFunctionIndex !::index matching an entry into an array of type MagStateFunction. Used by soft ferromagnets (when interpolation on an M vs H curve is necessary)
    integer(4),dimension(n_tiles),intent(in) :: includeInIteration,exploitSymmetry
    real(8),dimension(n_tiles,3),intent(in) :: symmetryOps !! 1 for symmetry, -1 for anti-symmetry ((1) for about xy plane, (2) for about (xz) plane and (3) for about yz plane)
    real(8),dimension(n_tiles),intent(in) :: Mrel !! the current relative change of the magnetization (i.e. abs((M1-M2)/M2 ) where M1 is the latest magnetization norm and M2 is the previous one


    real(8),dimension(n_pts,3),intent(in) :: pts
    real(8),dimension(n_pts,3) :: H    
    real(8),dimension(n_tiles,n_pts,3,3),intent(out) :: N
    type(MagTile),dimension(n_tiles) :: tiles
    integer(4),intent(in) :: n_tiles
    integer(4),intent(in) :: n_pts
    integer :: i

    !::initialise MagTile with specified parameters
    call loadTiles( centerPos, dev_center, rect_size, Mag, u_ea, u_oa1, u_oa2, mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
        includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, tiles )

    !::do the calculation
    N(:,:,:,:) = 0
    H(:,:) = 0
    do i=1,n_tiles
        if ( tiles(i)%tileType .eq. tileTypeCylPiece ) then
            call getFieldFromCylTile( tiles(i), H, pts, n_pts, N(i,:,:,:), .false. )
        else if (tiles(i)%tileType .eq. tileTypeCircPiece ) then          
            call getFieldFromCircPieceTile( tiles(i), H, pts, n_pts, N(i,:,:,:), .false. )
        else if (tiles(i)%tileType .eq. tileTypeCircPieceInverted ) then          
            call getFieldFromCircPieceInvertedTile( tiles(i), H, pts, n_pts, N(i,:,:,:), .false. )
        else if (tiles(i)%tileType .eq. tileTypePrism ) then    
            call getFieldFromRectangularPrismTile( tiles(i), H, pts, n_pts, N(i,:,:,:), .false. )
        endif
    enddo

end subroutine getNFromTiles


subroutine getHFromTiles( centerPos, dev_center, rect_size, Mag, u_ea, u_oa1, u_oa2, mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
    includeInIteration, exploitSymmetry, symmetryOps, Mrel, pts, n_tiles, n_pts, H, N, useProvidedN )
    !::Specific for a cylindrical tile piece
    real(8),dimension(n_tiles,3),intent(in) :: centerPos
    real(8),dimension(n_tiles,3),intent(in) :: dev_center
        
    !::Specific for a rectangular prism
    real(8),dimension(n_tiles,3),intent(in) :: rect_size
        
    !::Generel variables, shared among all tile types
    real(8),dimension(n_tiles,3),intent(in) :: Mag
    real(8),dimension(n_tiles,3),intent(in) :: u_ea,u_oa1,u_oa2    
    real(8),dimension(n_tiles),intent(in) :: mu_r_ea,mu_r_oa,Mrem
    integer(4),dimension(n_tiles),intent(in) :: tileType        !::defines whether the tile is cylindrical, a prism, an ellipsoid and so on
    real(8),dimension(n_tiles,3),intent(in) :: offset !::the centre coordinates relative to the global coordinate system
    real(8),dimension(n_tiles,3),intent(in) :: rotAngles !:: rotation angles (phi_x, phi_y, phi_z) about the principle axes of the tile with respect to the centre of the tile
    real(8),dimension(n_tiles,3),intent(in) :: color !! color rgb triplet
    integer(4),dimension(n_tiles),intent(in) :: magnetType !::defines whether the tile is a hard or soft magnet
    integer(4),dimension(n_tiles),intent(in) :: stateFunctionIndex !::index matching an entry into an array of type MagStateFunction. Used by soft ferromagnets (when interpolation on an M vs H curve is necessary)
    integer(4),dimension(n_tiles),intent(in) :: includeInIteration,exploitSymmetry
    real(8),dimension(n_tiles,3),intent(in) :: symmetryOps !! 1 for symmetry, -1 for anti-symmetry ((1) for about xy plane, (2) for about (xz) plane and (3) for about yz plane)
    real(8),dimension(n_tiles),intent(in) :: Mrel !! the current relative change of the magnetization (i.e. abs((M1-M2)/M2 ) where M1 is the latest magnetization norm and M2 is the previous one


    real(8),dimension(n_pts,3),intent(in) :: pts
    real(8),dimension(n_pts,3),intent(out) :: H
    real(8),dimension(n_pts,3) :: H_tmp
    real(8),dimension(n_tiles,n_pts,3,3),intent(in) :: N
    real(8),dimension(n_tiles,n_pts,3,3) :: N_tmp
    logical,intent(in) :: useProvidedN

    type(MagTile),dimension(n_tiles) :: tiles
    integer(4),intent(in) :: n_tiles
    integer(4),intent(in) :: n_pts
    integer :: i

    !::initialise MagTile with specified parameters
    call loadTiles( centerPos, dev_center, rect_size, Mag, u_ea, u_oa1, u_oa2, mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
        includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, tiles )

    if ( useProvidedN .eqv. .true. ) then
        write(*,*) 'Finding solution at requested points with precalculated N'

        N_tmp = N
        H(:,:) = 0.
    
        do i=1,n_tiles
        
            H_tmp(:,:) = 0.        

            select case (tiles(i)%tileType )

            case (tileTypeCylPiece)
                call getFieldFromCylTile( tiles(i), H_tmp, pts, n_pts, N_tmp(i,:,:,:), useProvidedN )

            case (tileTypePrism)
                call getFieldFromRectangularPrismTile( tiles(i), H_tmp, pts, n_pts, N_tmp(i,:,:,:), useProvidedN )

            case (tileTypeCircPiece )
                 call getFieldFromCircPieceTile( tiles(i), H_tmp, pts, n_pts, N_tmp(i,:,:,:), useProvidedN )

            case (tileTypeCircPieceInverted )
                 call getFieldFromCircPieceInvertedTile( tiles(i), H_tmp, pts, n_pts, N_tmp(i,:,:,:), useProvidedN )

            case (tileTypeEllipsoid)
                !!@todo add the existing code for spheroids, and correct Ellipsoids to Spheroids

            case (tileTypePlanarCoil )
                call getFieldFromPlanarCoilTile( tiles(i), H_tmp, pts, n_pts, N_tmp(i,:,:,:), useProvidedN )
            
            case default        
            
            end select
        
            H = H + H_tmp

        enddo
    
        call SubtractMFromCylindricalTiles( H, tiles, pts, n_tiles, n_pts)

    else
        write(*,*) 'Finding solution at requested points'
        call getFieldFromTiles( tiles, H, pts, n_tiles, n_pts )
    endif

end subroutine getHFromTiles


subroutine IterateTiles( centerPos, dev_center, rect_size, Mag, u_ea, u_oa1, u_oa2, mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
    includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, nT, nH, n_stateFcn, data_stateFcn, T, maxErr, nIteMax, Mag_out, Mrel_out )
    !::Specific for a cylindrical tile piece
    real(8),dimension(n_tiles,3),intent(in) :: centerPos
    real(8),dimension(n_tiles,3),intent(in) :: dev_center
        
    !::Specific for a rectangular prism
    real(8),dimension(n_tiles,3),intent(in) :: rect_size
        
    !::Generel variables, shared among all tile types
    real(8),dimension(n_tiles,3),intent(in) :: Mag
    real(8),dimension(n_tiles,3),intent(out) :: Mag_out
    real(8),dimension(n_tiles,3),intent(in) :: u_ea,u_oa1,u_oa2    
    real(8),dimension(n_tiles),intent(in) :: mu_r_ea,mu_r_oa,Mrem
    integer(4),dimension(n_tiles),intent(in) :: tileType        !::defines whether the tile is cylindrical, a prism, an ellipsoid and so on
    real(8),dimension(n_tiles,3),intent(in) :: offset !::the centre coordinates relative to the global coordinate system
    real(8),dimension(n_tiles,3),intent(in) :: rotAngles !:: rotation angles (phi_x, phi_y, phi_z) about the principle axes of the tile with respect to the centre of the tile
    real(8),dimension(n_tiles,3),intent(in) :: color !! color rgb triplet
    integer(4),dimension(n_tiles),intent(in) :: magnetType !::defines whether the tile is a hard or soft magnet
    integer(4),dimension(n_tiles),intent(in) :: stateFunctionIndex !::index matching an entry into an array of type MagStateFunction. Used by soft ferromagnets (when interpolation on an M vs H curve is necessary)
    integer(4),dimension(n_tiles),intent(in) :: includeInIteration,exploitSymmetry
    real(8),dimension(n_tiles,3),intent(in) :: symmetryOps !! 1 for symmetry, -1 for anti-symmetry ((1) for about xy plane, (2) for about (xz) plane and (3) for about yz plane)
    real(8),dimension(n_tiles),intent(in) :: Mrel !! the current relative change of the magnetization (i.e. abs((M1-M2)/M2 ) where M1 is the latest magnetization norm and M2 is the previous one
    real(8),dimension(n_tiles),intent(out) :: Mrel_out

    real(8) :: maxErr, T
    integer(4) :: nIteMax
    type(MagStateFunction),dimension(n_stateFcn) :: stateFcn
    integer(4),intent(in) :: n_stateFcn
    integer(4),intent(in) :: nT,nH
    real(8),dimension(nH,nT),intent(in) :: data_stateFcn
    real(8) :: resumeIteration
    procedure(displayIteration_fct),pointer :: disp_fct => null()

    type(MagTile),dimension(n_tiles) :: tiles
    integer(4),intent(in) :: n_tiles
    integer :: i

    disp_fct => dispIte_fct

    !::initialise MagTile with specified parameters
    call loadTiles( centerPos, dev_center, rect_size, Mag, u_ea, u_oa1, u_oa2, mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
        includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, tiles )
    
    !::load state function from table
    call loadMagStateFunction( nT, nH, stateFcn, data_stateFcn, n_stateFcn)
    
    write(*,*) 'Doing iteration'
    call iterateMagnetization( tiles, n_tiles, stateFcn, n_stateFcn, T, maxErr, nIteMax, disp_fct, resumeIteration )

    do i=1,n_tiles
        Mag_out(i,:) = tiles(i)%M
        Mrel_out(i) = tiles(i)%Mrel
    enddo

end subroutine IterateTiles


subroutine runSimulation( centerPos, dev_center, rect_size, Mag, u_ea, u_oa1, u_oa2, mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
    includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, n_stateFcn, nT, nH, data_stateFcn, T, maxErr, nIteMax, iterateSolution, returnSolution, n_pts, pts, H, Mag_out, Mrel_out )
    !::Specific for a cylindrical tile piece
    real(8),dimension(n_tiles,3),intent(in) :: centerPos
    real(8),dimension(n_tiles,3),intent(in) :: dev_center
        
    !::Specific for a rectangular prism
    real(8),dimension(n_tiles,3),intent(in) :: rect_size
        
    !::Generel variables, shared among all tile types
    real(8),dimension(n_tiles,3),intent(in) :: Mag
    real(8),dimension(n_tiles,3),intent(out) :: Mag_out
    real(8),dimension(n_tiles,3),intent(in) :: u_ea,u_oa1,u_oa2
    !! real(8),dimension(n_tiles,3),intent(out) :: u_ea_out,u_oa1_out,u_oa2_out 
    real(8),dimension(n_tiles),intent(in) :: mu_r_ea,mu_r_oa,Mrem
    integer(4),dimension(n_tiles),intent(in) :: tileType        !::defines whether the tile is cylindrical, a prism, an ellipsoid and so on
    real(8),dimension(n_tiles,3),intent(in) :: offset !::the centre coordinates relative to the global coordinate system
    real(8),dimension(n_tiles,3),intent(in) :: rotAngles !:: rotation angles (phi_x, phi_y, phi_z) about the principle axes of the tile with respect to the centre of the tile
    real(8),dimension(n_tiles,3),intent(in) :: color !! color rgb triplet
    integer(4),dimension(n_tiles),intent(in) :: magnetType !::defines whether the tile is a hard or soft magnet
    integer(4),dimension(n_tiles),intent(in) :: stateFunctionIndex !::index matching an entry into an array of type MagStateFunction. Used by soft ferromagnets (when interpolation on an M vs H curve is necessary)
    integer(4),dimension(n_tiles),intent(in) :: includeInIteration,exploitSymmetry
    real(8),dimension(n_tiles,3),intent(in) :: symmetryOps !! 1 for symmetry, -1 for anti-symmetry ((1) for about xy plane, (2) for about (xz) plane and (3) for about yz plane)
    real(8),dimension(n_tiles),intent(in) :: Mrel !! the current relative change of the magnetization (i.e. abs((M1-M2)/M2 ) where M1 is the latest magnetization norm and M2 is the previous one
    real(8),dimension(n_tiles),intent(out) :: Mrel_out

    real(8) :: maxErr, T
    integer(4) :: nIteMax
    logical,intent(in) :: iterateSolution, returnSolution
    type(MagStateFunction),dimension(n_stateFcn) :: stateFcn
    integer(4),intent(in) :: n_stateFcn !! Currently just one state function supported, can be extendend to input variable
    integer(4),intent(in):: nT,nH
    real(8),dimension(nH,nT),intent(in) :: data_stateFcn
    real(8) :: start,finish,resumeIteration
    procedure(displayIteration_fct),pointer :: disp_fct => null()

    real(8),dimension(n_pts,3),intent(in) :: pts
    real(8),dimension(n_pts,3),intent(out) :: H
    type(MagTile),dimension(n_tiles) :: tiles
    integer(4),intent(in) :: n_tiles
    integer(4),intent(in) :: n_pts
    integer :: i

    call cpu_time(start)

    disp_fct => dispIte_fct
    
    !!@todo no support for resuming iterations at the moment
    resumeIteration = 0

    !::initialise MagTile with specified parameters
    call loadTiles( centerPos, dev_center, rect_size, Mag, u_ea, u_oa1, u_oa2, mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
        includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, tiles )

    !::load state function from table
    call loadMagStateFunction( nT, nH, stateFcn, data_stateFcn, n_stateFcn )
    
    if ( iterateSolution .eqv. .true. ) then
        write(*,*) 'Doing iteration'
        call iterateMagnetization( tiles, n_tiles, stateFcn, n_stateFcn, T, maxErr, nIteMax, disp_fct, resumeIteration )

        do i=1,n_tiles
            ! centerPos(i,1) = tiles(i)%r0
            ! centerPos(i,2) = tiles(i)%theta0
            ! centerPos(i,3) = tiles(i)%z0
            ! dev_center(i,1) = tiles(i)%dr
            ! dev_center(i,2) = tiles(i)%dtheta
            ! dev_center(i,3) = tiles(i)%dz
            ! rect_size(i,1) = tiles(i)%a
            ! rect_size(i,2) = tiles(i)%b
            ! rect_size(i,3) = tiles(i)%c
            Mag_out(i,:) = tiles(i)%M
            ! u_ea_out(i,:) = tiles(i)%u_ea
            ! u_oa1_out(i,:) = tiles(i)%u_oa1
            ! u_oa2_out(i,:) = tiles(i)%u_oa2
            ! mu_r_ea(i) = tiles(i)%mu_r_ea
            ! mu_r_oa(i) = tiles(i)%mu_r_oa
            ! Mrem(i) = tiles(i)%Mrem
            ! tileType(i) = tiles(i)%tileType
            ! offset(i,:) = tiles(i)%offset
            ! rotAngles(i,:) = tiles(i)%rotAngles
            ! color(i,:) = tiles(i)%color
            ! magnetType(i) = tiles(i)%magnetType
            ! stateFunctionIndex(i) = tiles(i)%stateFunctionIndex
            ! includeInIteration(i) = tiles(i)%includeInIteration
            ! exploitSymmetry(i) = tiles(i)%exploitSymmetry
            ! symmetryOps(i,:) = tiles(i)%symmetryOps
            Mrel_out(i) = tiles(i)%Mrel
        enddo
    endif
    
    if ( returnSolution .eqv. .true. ) then
        write(*,*) 'Finding solution at requested points'
        call getFieldFromTiles( tiles, H, pts, n_tiles, n_pts )
    endif

    call cpu_time(finish)
    
    write(*,*) 'Elapsed time', finish-start

end subroutine runSimulation


end module FortranToPythonIO
