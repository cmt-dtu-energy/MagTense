module FortranToPythonIO

use DemagFieldGetSolution
use IterateMagnetSolution
use LandauLifshitzSolution
use IntegrationDataTypes
use MagTenseMicroMagIO
implicit none

contains

!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!> function for displaying output (to the terminal)
!! @param err the current relative error
!! @param err_max the current maximum allowed error
!!
function dispIte_fct( err, err_max )
    real(8),intent(in) :: err,err_max
    integer(4) :: dispIte_fct
   
    write(*,*) err,err_max
  
    dispIte_fct = 0
    
end function dispIte_fct

function dispIte_fct_no_output( err, err_max )
    real(8),intent(in) :: err,err_max
    integer(4) :: dispIte_fct_no_output
  
    dispIte_fct_no_output = 0
    
end function dispIte_fct_no_output


!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!< Function getNFromTiles
!! @param tileType - defines whether the tile is cylindrical, a prism, an ellipsoid and so on
!! @param offset - the centre coordinates relative to the global coordinate system
!! @param rotAngles - rotation angles (phi_x, phi_y, phi_z) about the principle axes of the tile with respect to the centre of the tile
!! @param color - color rgb triplet
!! @param magnetType - defines whether the tile is a hard or soft magnet
!! @param stateFunctionIndex - index matching an entry into an array of type MagStateFunction. Used by soft ferromagnets (when interpolation on an M vs H curve is necessary)
!! @param symmetryOps - 1 for symmetry, -1 for anti-symmetry ((1) for about xy plane, (2) for about (xz) plane and (3) for about yz plane)
!! @param Mrel - the current relative change of the magnetization (i.e. abs((M1-M2)/M2 ) where M1 is the latest magnetization norm and M2 is the previous one
!!
subroutine getNFromTiles( centerPos, dev_center, tile_size, vertices, Mag, u_ea, u_oa1, u_oa2, &
    mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
    includeInIteration, exploitSymmetry, symmetryOps, Mrel, pts, n_tiles, n_pts, N )
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

    real(8),dimension(n_pts,3),intent(in) :: pts
    real(8),dimension(n_pts,3) :: H    
    real(8),dimension(n_tiles,n_pts,3,3),intent(out) :: N
    type(MagTile),dimension(n_tiles) :: tiles
    integer(4),intent(in) :: n_tiles
    integer(4),intent(in) :: n_pts
    integer :: i

    !::initialise MagTile with specified parameters
    call loadTiles( centerPos, dev_center, tile_size, vertices, Mag, u_ea, u_oa1, u_oa2, &
        mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
        includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, tiles )

    !::do the calculation
    N(:,:,:,:) = 0
    H(:,:) = 0
    ! $OMP PARALLEL DO PRIVATE(i)
    do i=1,n_tiles
        select case ( tiles(i)%tileType )
        case ( tileTypeCylPiece )
            call getFieldFromCylTile( tiles(i), H, pts, n_pts, N(i,:,:,:), .false. )
        case ( tileTypePrism )
            call getFieldFromRectangularPrismTile( tiles(i), H, pts, n_pts, N(i,:,:,:), .false. )
        case ( tileTypeSphere )
            call getFieldFromSphereTile( tiles(i), H, pts, n_pts, N(i,:,:,:), .false. )
        case ( tileTypeSpheroid )
            call getFieldFromSpheroidTile( tiles(i), H, pts, n_pts, N(i,:,:,:), .false. )
        case ( tileTypeCircPiece )
            call getFieldFromCircPieceTile( tiles(i), H, pts, n_pts, N(i,:,:,:), .false. )
        case ( tileTypeCircPieceInverted )
            call getFieldFromCircPieceInvertedTile( tiles(i), H, pts, n_pts, N(i,:,:,:), .false. )
        case ( tileTypeTetrahedron )
            call getFieldFromTetrahedronTile( tiles(i), H, pts, n_pts, N(i,:,:,:), .false. )
        case (tileTypePlanarCoil )
            call getFieldFromPlanarCoilTile( tiles(i), H, pts, n_pts, N(i,:,:,:), .false. )
        case default
        end select
    enddo
    ! $OMP END PARALLEL DO

end subroutine getNFromTiles

!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!< Function getHFromTiles
!! @param tileType - defines whether the tile is cylindrical, a prism, an ellipsoid and so on
!! @param offset - the centre coordinates relative to the global coordinate system
!! @param rotAngles - rotation angles (phi_x, phi_y, phi_z) about the principle axes of the tile with respect to the centre of the tile
!! @param color - color rgb triplet
!! @param magnetType - defines whether the tile is a hard or soft magnet
!! @param stateFunctionIndex - index matching an entry into an array of type MagStateFunction. Used by soft ferromagnets (when interpolation on an M vs H curve is necessary)
!! @param symmetryOps - 1 for symmetry, -1 for anti-symmetry ((1) for about xy plane, (2) for about (xz) plane and (3) for about yz plane)
!! @param Mrel - the current relative change of the magnetization (i.e. abs((M1-M2)/M2 ) where M1 is the latest magnetization norm and M2 is the previous one
!!
subroutine getHFromTiles( centerPos, dev_center, tile_size, vertices, Mag, u_ea, u_oa1, u_oa2, &
    mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
    includeInIteration, exploitSymmetry, symmetryOps, Mrel, pts, n_tiles, n_pts, H, N, useStoredN )
    
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

    real(8),dimension(n_pts,3),intent(in) :: pts
    real(8),dimension(n_pts,3),intent(out) :: H
    real(8),dimension(n_pts,3) :: H_tmp
    real(8),dimension(n_tiles,n_pts,3,3),intent(in) :: N
    real(8),dimension(n_tiles,n_pts,3,3) :: N_out
    logical,intent(in) :: useStoredN

    type(MagTile),dimension(n_tiles) :: tiles
    integer(4),intent(in) :: n_tiles
    integer(4),intent(in) :: n_pts
    integer :: i

    !::initialise MagTile with specified parameters
    call loadTiles( centerPos, dev_center, tile_size, vertices, Mag, u_ea, u_oa1, u_oa2, &
        mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
        includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, tiles )

    N_out = N
    H(:,:) = 0.

    ! $OMP PARALLEL DO PRIVATE(i,H_tmp)    
    do i=1,n_tiles

        !Make sure to allocate H_tmp on the heap and for each thread
        ! $OMP CRITICAL
        H_tmp(:,:) = 0.

        ! $OMP END CRITICAL
        !! Here a selection of which subroutine to use should be done, i.e. whether the tile
        !! is cylindrical, a prism or an ellipsoid or another geometry
        select case ( tiles(i)%tileType )
        case ( tileTypeCylPiece )
            if ( useStoredN .eqv. .true. ) then
                call getFieldFromCylTile( tiles(i), H_tmp, pts, n_pts, N_out(i,:,:,:), useStoredN )    
            else
                call getFieldFromCylTile( tiles(i), H_tmp, pts, n_pts )
            endif
        case ( tileTypePrism )
            if ( useStoredN .eqv. .true. ) then
                call getFieldFromRectangularPrismTile( tiles(i), H_tmp, pts, n_pts, N_out(i,:,:,:), useStoredN )
            else
                call getFieldFromRectangularPrismTile( tiles(i), H_tmp, pts, n_pts )
            endif
        case ( tileTypeSphere )
            if ( useStoredN .eqv. .true. ) then
                call getFieldFromSphereTile( tiles(i), H_tmp, pts, n_pts, N_out(i,:,:,:), useStoredN )
            else
                call getFieldFromSphereTile( tiles(i), H_tmp, pts, n_pts )
            endif
        case ( tileTypeSpheroid )
            if ( useStoredN .eqv. .true. ) then
                call getFieldFromSpheroidTile( tiles(i), H_tmp, pts, n_pts, N_out(i,:,:,:), useStoredN )
            else
                call getFieldFromSpheroidTile( tiles(i), H_tmp, pts, n_pts )
            endif
        case ( tileTypeCircPiece )
            if ( useStoredN .eqv. .true. ) then
                call getFieldFromCircPieceTile( tiles(i), H_tmp, pts, n_pts, N_out(i,:,:,:), useStoredN )
            else
                call getFieldFromCircPieceTile( tiles(i), H_tmp, pts, n_pts )
            endif
        case ( tileTypeCircPieceInverted )
            if ( useStoredN .eqv. .true. ) then
                call getFieldFromCircPieceInvertedTile( tiles(i), H_tmp, pts, n_pts, N_out(i,:,:,:), useStoredN )
            else
                call getFieldFromCircPieceInvertedTile( tiles(i), H_tmp, pts, n_pts )
            endif
        case ( tileTypeTetrahedron )
            if ( useStoredN .eqv. .true. ) then
                call getFieldFromTetrahedronTile( tiles(i), H_tmp, pts, n_pts, N_out(i,:,:,:), useStoredN )
            else
                call getFieldFromTetrahedronTile( tiles(i), H_tmp, pts, n_pts )
            endif
        case ( tileTypePlanarCoil )
            if ( useStoredN .eqv. .true. ) then
                call getFieldFromPlanarCoilTile( tiles(i), H_tmp, pts, n_pts, N_out(i,:,:,:), useStoredN )
            else
                call getFieldFromPlanarCoilTile( tiles(i), H_tmp, pts, n_pts )            
            endif

        case default

        end select

        H = H + H_tmp

    enddo
    ! $OMP END PARALLEL DO

    call SubtractMFromCylindricalTiles( H, tiles, pts, n_tiles, n_pts )

end subroutine getHFromTiles

!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!< Function IterateTiles
!! @param tileType - defines whether the tile is cylindrical, a prism, an ellipsoid and so on
!! @param offset - the centre coordinates relative to the global coordinate system
!! @param rotAngles - rotation angles (phi_x, phi_y, phi_z) about the principle axes of the tile with respect to the centre of the tile
!! @param color - color rgb triplet
!! @param magnetType - defines whether the tile is a hard or soft magnet
!! @param stateFunctionIndex - index matching an entry into an array of type MagStateFunction. Used by soft ferromagnets (when interpolation on an M vs H curve is necessary)
!! @param symmetryOps - 1 for symmetry, -1 for anti-symmetry ((1) for about xy plane, (2) for about (xz) plane and (3) for about yz plane)
!! @param Mrel - the current relative change of the magnetization (i.e. abs((M1-M2)/M2 ) where M1 is the latest magnetization norm and M2 is the previous one
!!
subroutine IterateTiles( centerPos, dev_center, tile_size, vertices, Mag, u_ea, u_oa1, u_oa2, &
    mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
    includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, nT, nH, n_stateFcn, &
    data_stateFcn, T, maxErr, nIteMax, Mag_out, Mrel_out )
    
    !::Specific for a cylindrical tile piece
    real(8),dimension(n_tiles,3),intent(in) :: centerPos
    real(8),dimension(n_tiles,3),intent(in) :: dev_center
        
    !::Specific for a rectangular prism
    real(8),dimension(n_tiles,3),intent(in) :: tile_size
    
    !::Specific for a tetrahedron
    real(8),dimension(n_tiles,3,4),intent(in) :: vertices

    !::Generel variables, shared among all tile types
    real(8),dimension(n_tiles,3),intent(in) :: Mag
    real(8),dimension(n_tiles,3),intent(out) :: Mag_out
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
    real(8),dimension(n_tiles),intent(out) :: Mrel_out

    real(8),intent(in) :: maxErr, T
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

    !! default value is zero, i.e. do not resume iteration
    resumeIteration = 0.
    
    disp_fct => dispIte_fct_no_output

    !::initialise MagTile with specified parameters
    call loadTiles( centerPos, dev_center, tile_size, vertices, Mag, u_ea, u_oa1, u_oa2, &
        mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
        includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, tiles )
    
    !::load state function from table
    call loadStateFunction( nT, nH, stateFcn, data_stateFcn, n_stateFcn)
    call iterateMagnetization( tiles, n_tiles, stateFcn, n_stateFcn, T, maxErr, nIteMax, disp_fct, resumeIteration )

    do i=1,n_tiles
        Mag_out(i,:) = tiles(i)%M
        Mrel_out(i) = tiles(i)%Mrel
    enddo

end subroutine IterateTiles

!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!< Function runSimulation
!! @param tileType - defines whether the tile is cylindrical, a prism, an ellipsoid and so on
!! @param offset - the centre coordinates relative to the global coordinate system
!! @param rotAngles - rotation angles (phi_x, phi_y, phi_z) about the principle axes of the tile with respect to the centre of the tile
!! @param color - color rgb triplet
!! @param magnetType - defines whether the tile is a hard or soft magnet
!! @param stateFunctionIndex - index matching an entry into an array of type MagStateFunction. Used by soft ferromagnets (when interpolation on an M vs H curve is necessary)
!! @param symmetryOps - 1 for symmetry, -1 for anti-symmetry ((1) for about xy plane, (2) for about (xz) plane and (3) for about yz plane)
!! @param Mrel - the current relative change of the magnetization (i.e. abs((M1-M2)/M2 ) where M1 is the latest magnetization norm and M2 is the previous one
!!
subroutine runSimulation( centerPos, dev_center, tile_size, vertices, Mag, u_ea, u_oa1, u_oa2, &
    mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
    includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, n_stateFcn, nT, nH, &
    data_stateFcn, T, maxErr, nIteMax, iterateSolution, returnSolution, n_pts, pts, H, Mag_out, Mrel_out, console )
    
    !::Specific for a cylindrical tile piece
    real(8),dimension(n_tiles,3),intent(in) :: centerPos
    real(8),dimension(n_tiles,3),intent(in) :: dev_center
        
    !::Specific for a rectangular prism
    real(8),dimension(n_tiles,3),intent(in) :: tile_size

    !::Specific for a tetrahedron
    real(8),dimension(n_tiles,3,4),intent(in) :: vertices
        
    !::Generel variables, shared among all tile types
    real(8),dimension(n_tiles,3),intent(in) :: Mag
    real(8),dimension(n_tiles,3),intent(out) :: Mag_out
    real(8),dimension(n_tiles,3),intent(in) :: u_ea,u_oa1,u_oa2
    !! real(8),dimension(n_tiles,3),intent(out) :: u_ea_out,u_oa1_out,u_oa2_out 
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
    logical, intent(in) :: console

    real(8),dimension(n_pts,3),intent(in) :: pts
    real(8),dimension(n_pts,3),intent(out) :: H
    real(8),dimension(n_pts,3) :: H_tmp
    type(MagTile),dimension(n_tiles) :: tiles
    integer(4),intent(in) :: n_tiles
    integer(4),intent(in) :: n_pts
    integer :: i

    call cpu_time(start)

    if (console) then
        disp_fct => dispIte_fct
    else
        disp_fct => dispIte_fct_no_output
    endif
    
    !!@todo no support for resuming iterations at the moment
    resumeIteration = 0

    !::initialise MagTile with specified parameters
    call loadTiles( centerPos, dev_center, tile_size, vertices, Mag, u_ea, u_oa1, u_oa2, &
        mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
        includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, tiles )

    !::load state function from table
    call loadStateFunction( nT, nH, stateFcn, data_stateFcn, n_stateFcn )
    
    if ( iterateSolution .eqv. .true. ) then
        if (console) then
            write(*,*) 'Doing iteration'
        endif
        call iterateMagnetization( tiles, n_tiles, stateFcn, n_stateFcn, T, maxErr, nIteMax, disp_fct, resumeIteration )

        do i=1,n_tiles
            Mag_out(i,:) = tiles(i)%M
            Mrel_out(i) = tiles(i)%Mrel
        enddo
    endif
    
    if ( returnSolution .eqv. .true. ) then
        if (console) then
            write(*,*) 'Finding solution at requested points'
        endif

        H(:,:) = 0.
    
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
                call getFieldFromCylTile( tiles(i), H_tmp, pts, n_pts )
            case (tileTypePrism)
                call getFieldFromRectangularPrismTile( tiles(i), H_tmp, pts, n_pts )
            case (tileTypeSphere)
                call getFieldFromSphereTile( tiles(i), H_tmp, pts, n_pts )
            case (tileTypeSpheroid)
                call getFieldFromSpheroidTile( tiles(i), H_tmp, pts, n_pts )
            case (tileTypeCircPiece )
                call getFieldFromCircPieceTile( tiles(i), H_tmp, pts, n_pts )
            case (tileTypeCircPieceInverted )
                call getFieldFromCircPieceInvertedTile( tiles(i), H_tmp, pts, n_pts )
            case (tileTypeTetrahedron )
                call getFieldFromTetrahedronTile( tiles(i), H_tmp, pts, n_pts )
            case (tileTypePlanarCoil )
                call getFieldFromPlanarCoilTile( tiles(i), H_tmp, pts, n_pts )
            case default
            end select
    
            H = H + H_tmp
    
        enddo
        ! $OMP END PARALLEL DO
    
        call SubtractMFromCylindricalTiles( H, tiles, pts, n_tiles, n_pts )
    endif

    call cpu_time(finish)
    
    if (console) then
        write(*,*) 'Elapsed time', finish-start
    endif

end subroutine runSimulation


subroutine RunMicroMagSimulation( ntot, grid_n, grid_L, grid_type, u_ea, ProblemMode, solver, A0, Ms, K0, &
    gamma, alpha, MaxT0, nt_Hext, Hext, nt, t, m0, dem_thres, useCuda, dem_appr, N_ret, N_file_out, &
    N_load, N_file_in, setTimeDis, nt_alpha, alphat, tol, thres, useCVODE, nt_conv, t_conv, &
    conv_tol, grid_pts, grid_ele, grid_nod, grid_nnod, exch_nval, exch_nrow, exch_val, exch_rows, &
    exch_rowe, exch_col, grid_abc )

    type(MicroMagProblem) :: problem
    type(MicroMagSolution) :: solution
    real(8),dimension(nt),intent(out) :: t
    real(8),dimension(nt,ntot,ndim,3),intent(out) :: M

    real(8),dimension(n_pts,3),intent(out) :: pts, H_exc, H_ext, H_dem, H_ani


    call loadMicroMagProblemPy( ntot, grid_n, grid_L, grid_type, u_ea, ProblemMode, solver, A0, Ms, K0, &
        gamma, alpha, MaxT0, nt_Hext, Hext, nt, t, m0, dem_thres, useCuda, dem_appr, N_ret, N_file_out, &
        N_load, N_file_in, setTimeDis, nt_alpha, alphat, tol, thres, useCVODE, nt_conv, t_conv, &
        conv_tol, grid_pts, grid_ele, grid_nod, grid_nnod, exch_nval, exch_nrow, exch_val, exch_rows, &
        exch_rowe, exch_col, grid_abc, problem )

    call SolveLandauLifshitzEquation( problem, solution )

    call returnMicroMagSolutionPy( solution, nt, ntot, ndim, t, M, pts, H_exc, H_ext, H_dem, H_ani )

end subroutine RunMicroMagSimulation


end module FortranToPythonIO
