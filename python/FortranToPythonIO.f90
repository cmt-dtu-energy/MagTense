module FortranToPythonIO

    use MagParameters
    use DemagFieldGetSolution
    use IterateMagnetSolution
#if USE_MICROMAG
    use MagTenseMicroMagPyIO
    use LandauLifshitzSolution
#endif
    implicit none

    contains

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

        integer(4),intent(in) :: n_tiles, n_pts

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

        integer(4),intent(in) :: n_tiles, n_pts

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
        real(8),dimension(n_tiles,n_pts,3,3),intent(inout) :: N
        logical,intent(in) :: useStoredN

        type(MagTile),dimension(n_tiles) :: tiles
        integer :: i

        !::initialise MagTile with specified parameters
        call loadTiles( centerPos, dev_center, tile_size, vertices, Mag, u_ea, u_oa1, u_oa2, &
            mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
            includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, tiles )

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
                    call getFieldFromCylTile( tiles(i), H_tmp, pts, n_pts, N(i,:,:,:), useStoredN )
                else
                    call getFieldFromCylTile( tiles(i), H_tmp, pts, n_pts )
                endif
            case ( tileTypePrism )
                if ( useStoredN .eqv. .true. ) then
                    call getFieldFromRectangularPrismTile( tiles(i), H_tmp, pts, n_pts, N(i,:,:,:), useStoredN )
                else
                    call getFieldFromRectangularPrismTile( tiles(i), H_tmp, pts, n_pts )
                endif
            case ( tileTypeSphere )
                if ( useStoredN .eqv. .true. ) then
                    call getFieldFromSphereTile( tiles(i), H_tmp, pts, n_pts, N(i,:,:,:), useStoredN )
                else
                    call getFieldFromSphereTile( tiles(i), H_tmp, pts, n_pts )
                endif
            case ( tileTypeSpheroid )
                if ( useStoredN .eqv. .true. ) then
                    call getFieldFromSpheroidTile( tiles(i), H_tmp, pts, n_pts, N(i,:,:,:), useStoredN )
                else
                    call getFieldFromSpheroidTile( tiles(i), H_tmp, pts, n_pts )
                endif
            case ( tileTypeCircPiece )
                if ( useStoredN .eqv. .true. ) then
                    call getFieldFromCircPieceTile( tiles(i), H_tmp, pts, n_pts, N(i,:,:,:), useStoredN )
                else
                    call getFieldFromCircPieceTile( tiles(i), H_tmp, pts, n_pts )
                endif
            case ( tileTypeCircPieceInverted )
                if ( useStoredN .eqv. .true. ) then
                    call getFieldFromCircPieceInvertedTile( tiles(i), H_tmp, pts, n_pts, N(i,:,:,:), useStoredN )
                else
                    call getFieldFromCircPieceInvertedTile( tiles(i), H_tmp, pts, n_pts )
                endif
            case ( tileTypeTetrahedron )
                if ( useStoredN .eqv. .true. ) then
                    call getFieldFromTetrahedronTile( tiles(i), H_tmp, pts, n_pts, N(i,:,:,:), useStoredN )
                else
                    call getFieldFromTetrahedronTile( tiles(i), H_tmp, pts, n_pts )
                endif
            case ( tileTypePlanarCoil )
                if ( useStoredN .eqv. .true. ) then
                    call getFieldFromPlanarCoilTile( tiles(i), H_tmp, pts, n_pts, N(i,:,:,:), useStoredN )
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
        
        integer(4),intent(in) :: n_tiles, n_stateFcn, nT, nH
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
        real(8),dimension(nH,nT),intent(in) :: data_stateFcn
        real(8) :: resumeIteration

        type(MagTile),dimension(n_tiles) :: tiles
        integer :: i

        !! default value is zero, i.e. do not resume iteration
        resumeIteration = 0.

        !::initialise MagTile with specified parameters
        call loadTiles( centerPos, dev_center, tile_size, vertices, Mag, u_ea, u_oa1, u_oa2, &
            mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
            includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, tiles )

        !::load state function from table
        call loadStateFunction( nT, nH, stateFcn, data_stateFcn, n_stateFcn)
        call iterateMagnetization( tiles, n_tiles, stateFcn, n_stateFcn, T, maxErr, nIteMax, resumeIteration )

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
        data_stateFcn, T, maxErr, nIteMax, iterateSolution, returnSolution, n_pts, pts, H, Mag_out, Mrel_out )

        integer(4),intent(in) :: n_tiles, n_pts, n_stateFcn, nT, nH

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
        real(8),dimension(nH,nT),intent(in) :: data_stateFcn
        real(8) :: start,finish,resumeIteration

        real(8),dimension(n_pts,3),intent(in) :: pts
        real(8),dimension(n_pts,3),intent(out) :: H
        real(8),dimension(n_pts,3) :: H_tmp
        type(MagTile),dimension(n_tiles) :: tiles
        integer :: i

        call cpu_time(start)

        !!@todo no support for resuming iterations at the moment
        resumeIteration = 0

        !::initialise MagTile with specified parameters
        call loadTiles( centerPos, dev_center, tile_size, vertices, Mag, u_ea, u_oa1, u_oa2, &
            mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
            includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, tiles )

        !::load state function from table
        call loadStateFunction( nT, nH, stateFcn, data_stateFcn, n_stateFcn )

        if ( iterateSolution .eqv. .true. ) then
            write(*,*) 'Doing iteration'
            call iterateMagnetization( tiles, n_tiles, stateFcn, n_stateFcn, T, maxErr, nIteMax, resumeIteration )

            do i=1,n_tiles
                Mag_out(i,:) = tiles(i)%M
                Mrel_out(i) = tiles(i)%Mrel
            enddo
        endif

        if ( returnSolution .eqv. .true. ) then
            write(*,*) 'Finding solution at requested points'

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
                
                ! $OMP CRITICAL
                H = H + H_tmp
                ! $OMP END CRITICAL

            enddo
            ! $OMP END PARALLEL DO

            call SubtractMFromCylindricalTiles( H, tiles, pts, n_tiles, n_pts )
        endif

        call cpu_time(finish)

        write(*,*) 'Elapsed time', finish-start

    end subroutine runSimulation


    subroutine RunMicroMagSimulation( ntot, grid_n, grid_L, grid_type, u_ea, ProblemMode, solver, A0, Ms, K0, &
        gamma, alpha_mm, MaxT0, nt_Hext, Hext, nt, t, m0, dem_thres, useCuda, dem_appr, N_ret, N_file_out, &
        N_load, N_file_in, setTimeDis, nt_alpha, alphat, tol, thres, useCVODE, nt_conv, t_conv, &
        conv_tol, grid_pts, grid_ele, grid_nod, grid_nnod, exch_nval, exch_nrow, exch_val, exch_rows, &
        exch_rowe, exch_col, grid_abc, usePrecision, nThreadsMatlab, N_ave, t_out, M_mm, pts, H_exc, H_ext, H_dem, H_ani)

        integer(4),intent(in) :: ntot, nt_Hext, nt, nt_alpha, nt_conv, grid_nnod, exch_nval, exch_nrow
        integer(4),dimension(3),intent(in) :: grid_n, N_ave
        real(8),dimension(3),intent(in) :: grid_L
        real(8),dimension(ntot,3),intent(in) :: grid_pts
        integer(4),dimension(4,ntot),intent(in) :: grid_ele
        real(8),dimension(grid_nnod,3),intent(in) :: grid_nod
        real(8),dimension(ntot, 3),intent(in) :: grid_abc, u_ea
        real(8),dimension(nt_Hext, 4),intent(in) :: Hext
        real(8),dimension(3*ntot),intent(in) :: m0
        real(8),dimension(nt_alpha,2),intent(in) :: alphat
        integer(4),dimension(exch_nval),intent(in) :: exch_val, exch_col
        integer(4),dimension(exch_nrow),intent(in) :: exch_rows, exch_rowe
        real(8),dimension(nt_conv),intent(in) :: t_conv
        integer(4),intent(in) :: ProblemMode, solver, grid_type, dem_appr, usePrecision, nThreadsMatlab
        integer(4),intent(in) :: N_ret, N_load, setTimeDis, useCuda, useCVODE
        real(8),intent(in) :: A0, Ms, K0, gamma, alpha_mm, MaxT0, tol, thres, conv_tol, dem_thres
        character*256,intent(in) :: N_file_in, N_file_out

        real(8),dimension(nt),intent(in) :: t
        real(8),dimension(nt),intent(out) :: t_out
        real(8),dimension(nt,ntot,1,3),intent(out) :: M_mm
        real(8),dimension(nt,ntot,1,3),intent(out) :: H_exc, H_ext, H_dem, H_ani
        real(8),dimension(ntot,3),intent(out) :: pts

#if USE_MICROMAG
        type(MicroMagProblem) :: problem
        type(MicroMagSolution) :: solution

        call loadMicroMagProblem( ntot, grid_n, grid_L, grid_type, u_ea, ProblemMode, solver, A0, Ms, K0, &
            gamma, alpha_mm, MaxT0, nt_Hext, Hext, nt, t, m0, dem_thres, useCuda, dem_appr, N_ret, N_file_out, &
            N_load, N_file_in, setTimeDis, nt_alpha, alphat, tol, thres, useCVODE, nt_conv, t_conv, &
            conv_tol, grid_pts, grid_ele, grid_nod, grid_nnod, exch_nval, exch_nrow, exch_val, exch_rows, &
            exch_rowe, exch_col, grid_abc, usePrecision, nThreadsMatlab, N_ave, problem )

        call SolveLandauLifshitzEquation( problem, solution )

        t_out = solution%t_out
        M_mm = solution%M_out
        pts = solution%pts
        H_exc = solution%H_exc
        H_ext = solution%H_ext
        H_dem = solution%H_dem
        H_ani = solution%H_ani
#else
        write(*,*) 'Compiled without micromagnetic part. Returning zeros.'
        t_out = 0.
        M_mm(:,:,:,:) = 0.
        pts(:,:) = 0.
        H_exc(:,:,:,:) = 0.
        H_ext(:,:,:,:) = 0.
        H_dem(:,:,:,:) = 0.
        H_ani(:,:,:,:) = 0.
#endif

    end subroutine RunMicroMagSimulation

end module FortranToPythonIO
