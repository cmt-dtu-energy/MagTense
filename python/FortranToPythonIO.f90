module FortranToPythonIO
    use auxInit_mod
    use timer_mod
    use trace_mod
    use TileNComponents
    use TileTetrahedronPairTensor
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
        includeInIteration, exploitSymmetry, symmetryOps, Mrel, pts, n_tiles, n_pts, N , Obs_size)

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
        real(8),dimension(n_pts,3),intent(in), optional :: Obs_size
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
            case ( tileTypeAvgPrism )
                if (present(Obs_size)) then
                    call getAvgFieldFromRPT(tiles(i), H, pts, n_pts, N(i,:,:,:), .false., Obs_size)
                else
                    call getAvgFieldFromRPT(tiles(i), H, pts, n_pts, N(i,:,:,:), .false.)
                end if
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
            case (tileTypeUniformField )
                call getFieldFromUniformFieldTile( tiles(i), H, pts, n_pts, N(i,:,:,:), .false. )
            case default
            end select
        enddo
        ! $OMP END PARALLEL DO

    end subroutine getNFromTiles

    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !< Function getNTensorTetrahedronPair
    !! Analytically exact, target-volume-averaged demagnetization tensor between a
    !! uniformly magnetized SOURCE tetrahedron and a TARGET tetrahedron, i.e. the
    !! tensor N with < H >_{V_target} = N M. This is the finite-target counterpart
    !! of the tetrahedron-source to POINT-target tensor reached through
    !! getNFromTiles with tileType = 5, which is unchanged.
    !!
    !! @param source_vertices - (3,4) vertices of the source tetrahedron, column k being vertex k
    !! @param target_vertices - (3,4) vertices of the receiving tetrahedron, same layout and frame
    !! @param N - the resulting (3,3) tensor
    !! @param ierr - 0 on success, non-zero for a degenerate tetrahedron or a non-finite result
    !!
    subroutine getNTensorTetrahedronPair( source_vertices, target_vertices, N, ierr )

        real(8),dimension(3,4),intent(in) :: source_vertices
        real(8),dimension(3,4),intent(in) :: target_vertices
        real(8),dimension(3,3),intent(out) :: N
        integer(4),intent(out) :: ierr

        integer :: status

        call getN_tensor_tetrahedron_tetrahedron( source_vertices, target_vertices, N, status )

        ierr = status

    end subroutine getNTensorTetrahedronPair

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
        includeInIteration, exploitSymmetry, symmetryOps, Mrel, pts, n_tiles, n_pts, H, N, useStoredN, Obs_size )

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
        real(8),intent(in),dimension(n_pts,3), optional :: Obs_size

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
            case ( tileTypeAvgPrism )
                    if (useStoredN) then
                        if (present(Obs_size)) then
                            call getAvgFieldFromRPT( tiles(i), H_tmp, pts, n_pts, N(i,:,:,:), useStoredN=useStoredN, Obs_size=Obs_size )
                        else
                            call getAvgFieldFromRPT( tiles(i), H_tmp, pts, n_pts, N(i,:,:,:), useStoredN=useStoredN )
                        end if
                    else
                        if (present(Obs_size)) then
                            call getAvgFieldFromRPT( tiles(i), H_tmp, pts, n_pts, Obs_size=Obs_size )
                        else
                            call getAvgFieldFromRPT( tiles(i), H_tmp, pts, n_pts )
                        end if
                    end if
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
            case ( tileTypeUniformField )
                !Not a geometry: a uniform applied field, the same at every point
                if ( useStoredN .eqv. .true. ) then
                    call getFieldFromUniformFieldTile( tiles(i), H_tmp, pts, n_pts, N(i,:,:,:), useStoredN )
                else
                    call getFieldFromUniformFieldTile( tiles(i), H_tmp, pts, n_pts )
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
        data_stateFcn, T, maxErr, nIteMax, iterateSolution, returnSolution, n_pts, pts, H, Mag_out, Mrel_out, Obs_size )
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
        real(8),dimension(n_pts,3),intent(in), optional :: Obs_size
        type(MagTile),dimension(n_tiles) :: tiles
        integer :: i

        call cpu_time(start)

        !!@todo no support for resuming iterations at the moment
        resumeIteration = 0
        !::initialise MagTile with specified parameters
        call loadTiles( centerPos, dev_center, tile_size, vertices, Mag, u_ea, u_oa1, u_oa2, &
            mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
            includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, tiles )
        !print *, "After loadtiles, M:", tiles(1)%M
        !print *, "Mrem: ", Mrem
        !::load state function from table
        call loadStateFunction( nT, nH, stateFcn, data_stateFcn, n_stateFcn )
        !print *, "After loadstatefunction, M:", tiles(1)%M
        !print *, "Mrem: ", Mrem

        if ( iterateSolution .eqv. .true. ) then
            call iterateMagnetization( tiles, n_tiles, stateFcn, n_stateFcn, T, maxErr, nIteMax, resumeIteration )

            do i=1,n_tiles
                Mag_out(i,:) = tiles(i)%M
                Mrel_out(i) = tiles(i)%Mrel
            enddo
        endif
        !print *, "After iterateMagnetization, M:", tiles(1)%M
        !print *, "Mrem: ", Mrem

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
                case (tileTypeAvgPrism)
                    !write(*,*)  "TileType avgPrism detected"
                    if (present(Obs_size)) then
                        !write(*,*)  "This is what we want!"
                       ! print *, "Obs_size = ", Obs_size
                        call getAvgFieldFromRPT(tiles(i), H_tmp, pts, n_pts, Obs_size=Obs_size)
                    else
                        call getAvgFieldFromRPT(tiles(i), H_tmp, pts, n_pts)
                    end if
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
                case (tileTypeUniformField )
                    !Not a geometry: a uniform applied field, the same at every point
                    call getFieldFromUniformFieldTile( tiles(i), H_tmp, pts, n_pts )
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
        K1, K2, K0_arr, CrysAxis, gamma, alpha_mm, temperature, nt_Hext, nt_Hext_out, Hext, nt, t, m0, dem_thres, useCuda, dem_appr, N_ret, N_file_out, &
        N_load, N_file_in, setTimeDis, nt_alpha, alphat, tol, thres, useCVODE, nt_conv, t_conv, &
        conv_tol, grid_pts, grid_ele, grid_nod, grid_nnod, exch_nval, exch_nrow, exch_val, exch_rows, &
        exch_cols, grid_abc, N_ave, CV, useReturnHall, useAvgN, & 
		exch_weigh, exch_meth, exch_intpn, passExch, exch_ncols, exch_presize, &
        n_macro, shiftVec, macroShape, sampleShape, exchPBC, hysteresis_solver, &
        H_start, H_end, dH_initial, dH_min, dH_max, maxHextSteps, dM_min, dM_target, dM_reject, dH_grow, dH_shrink, switch_refine_dH, use_switch_refine, &
        min_tol, min_maxiter, min_maxrot, min_fallback, min_saddle_check, min_predictor, &
        t_out, M_mm, pts, H_exc, H_ext, H_dem, H_ani, n_Hext_accepted, &
		n_tot_Exch, ExchMat_r, ExchMat_c, ExchMat_v, ExchMat_nr, ExchMat_nc, dummy_run, &
        useCDFMM, cdfmm_order, cdfmm_depth, cdfmm_basis, cdfmm_precision, &
        log_dir,timer_log_file, trace_log_file, window_enabled, window_interval, trace_enabled, flush_each, trace_verbose, timer_enabled, useDemag, rng_seed, &
        n_phase, phase_id, A_int, &
        E_out, n_feval, min_iter, min_torque, min_status, min_eig )

        !nt_Hext is the number of rows in the Hext array; nt_Hext_out is the third extent of the
        !returned M and H arrays. There used to be a third, n_Hext, which was accepted and declared
        !but never referenced anywhere in this routine - the number of distinct applied fields is
        !taken from size(problem%Hext,1) inside the solver. It has been removed rather than left as
        !a parameter that looks meaningful from the Python side.
        integer(4), intent(in) :: ntot, nt_conv, grid_type, nt_Hext, nt_alpha, nt, grid_nnod, exch_nval, exch_nrow, exch_ncols, exch_presize
        integer(4), intent(in) :: nt_Hext_out
        integer(4),dimension(3),intent(in) :: grid_n, N_ave
        real(8),dimension(3),intent(in) :: grid_L
        real(8),dimension(3),intent(in) :: H_start, H_end
        real(8),intent(in) :: dH_initial, dH_min, dH_max, dM_min, dM_target, dM_reject, dH_grow, dH_shrink, switch_refine_dH
        integer(4),intent(in) :: hysteresis_solver, maxHextSteps, use_switch_refine
        !> Energy minimizer settings (solver = 3), see MicroMagProblem
        real(8),intent(in) :: min_tol, min_maxrot
        integer(4),intent(in) :: min_maxiter, min_fallback, min_saddle_check, min_predictor
        real(8),dimension(ntot,3),intent(in) :: grid_pts
        integer(4),dimension(4,ntot),intent(in) :: grid_ele
        real(8),dimension(grid_nnod,3),intent(in) :: grid_nod
        real(8),dimension(ntot, 3),intent(in) :: grid_abc, u_ea
        real(8),dimension(nt_Hext, 4),intent(in) :: Hext
        real(8),dimension(3*ntot),intent(in) :: m0
        real(8),dimension(nt_alpha,2),intent(in) :: alphat
        !The exchange operator values are real. Declaring them integer(4) here made f2py generate
        !an 'i' argument, and it force-casts the float64 array from Python to int32 without an
        !error - which destroys entries of order 1/dx**2. rows/cols stay integer.
        real(8),dimension(exch_nval),intent(in) :: exch_val
        integer(4),dimension(exch_nval),intent(in) :: exch_cols, exch_rows
        real(8),dimension(nt_conv),intent(in) :: t_conv
		integer(4),intent(in) :: ProblemMode, solver, useCuda, dem_appr, useAvgN
		integer(4),intent(in) :: N_ret, N_load, setTimeDis, useCVODE, useReturnHall, exch_meth, exch_intpn, passExch, useDemag
        real(8),intent(in) :: gamma, alpha_mm, tol, thres, conv_tol, dem_thres
		real(8),dimension(ntot),intent(in) :: A0, Ms, K0, K1, K2, temperature
        !> Optional interface exchange between materials. n_phase = 1 disables it.
        integer(4),intent(in) :: n_phase
        integer(4),dimension(ntot),intent(in) :: phase_id
        real(8),dimension(n_phase,n_phase),intent(in) :: A_int
        real(8),dimension(ntot,6,3),intent(in) :: K0_arr
        real(8),dimension(ntot,3,3),intent(in):: CrysAxis
        character*256,intent(in) :: N_file_in, N_file_out
		real(8), intent(in) :: CV, exch_weigh

        integer(4),dimension(3),intent(in) :: n_macro
        real(8),dimension(3),intent(in) :: shiftVec, macroShape, sampleShape
        integer(4),dimension(3),intent(in) :: exchPBC

        real(8),dimension(nt),intent(in) :: t
        real(8),dimension(nt),intent(out) :: t_out
        real(8),dimension(nt,ntot,nt_Hext_out,3),intent(out) :: M_mm
        real(8),dimension(nt,ntot,nt_Hext_out,3),intent(out) :: H_exc, H_ext, H_dem, H_ani
        integer(4),intent(out) :: n_Hext_accepted
        real(8),dimension(ntot,3),intent(out) :: pts
        !> Energies (exchange, external, demag, anisotropy) [J] at the output times and applied fields,
        !> and the relaxation diagnostics per applied field: field evaluations spent, minimizer
        !> iterations, final max relative torque and status (-1 LL integration, 0 converged, 1 converged
        !> after an LL fallback, 2 not converged)
        real(8),dimension(nt,nt_Hext_out,4),intent(out) :: E_out
        integer(4),dimension(nt_Hext_out),intent(out) :: n_feval, min_iter, min_status
        real(8),dimension(nt_Hext_out),intent(out) :: min_torque, min_eig
        integer :: n_copy
		
		integer,intent(out) :: n_tot_Exch
		integer,dimension(exch_presize*ntot),intent(out)  :: ExchMat_r
		integer,dimension(exch_presize*ntot),intent(out)  :: ExchMat_c
		real(8),dimension(exch_presize*ntot),intent(out)  :: ExchMat_v
		integer,intent(out) :: ExchMat_nr,ExchMat_nc
        
        integer, intent(in) :: dummy_run
        integer(4), intent(in) :: useCDFMM
        integer(4), intent(in) :: cdfmm_order, cdfmm_depth, cdfmm_basis
        integer(4), intent(in) :: cdfmm_precision
        !> Seed for the stochastic thermal field: 0 keeps the compiler default sequence (identical
        !> on every run), a positive value seeds deterministically, a negative value seeds from the
        !> clock so that repeated runs are independent Monte-Carlo samples.
        integer(4), intent(in) :: rng_seed

        !-------------------- timer and trace modules --------------------------------------
        character*256,intent(in) :: timer_log_file, trace_log_file, log_dir
        integer, intent(in) :: window_enabled, trace_enabled, flush_each
        integer, intent(in) :: timer_enabled                 !> 1 writes the timing log file
        real(8), intent(in) :: window_interval
        integer, intent(in) :: trace_verbose
        !-----------------------------------------------------------------------------------
        
        logical :: window_enabled_l, trace_enabled_l, flush_each_l
        logical :: use_cdfmm
        !! Local auxiliary instance - avoids referencing the auxInit_mod module
        !! global across the f2py name-mangling boundary on Windows.
        type(auxInit_t) :: auxInit_local

#if USE_MICROMAG
        type(MicroMagProblem) :: problem
        type(MicroMagSolution) :: solution



        !---------------------- initiaize auxiliary modules -----------------------------
        call initAux(auxInit_local, log_dir, timer_log_file, trace_log_file, window_enabled, &
            window_interval, trace_enabled, flush_each, trace_verbose, timer_enabled)
        !---------------------------------------------------------------------------------




        !exchPBC is an integer flag per direction at both ends of this call - loadMicroMagProblem
        !takes integer(4) and stores it in problem%macrogrid%exchPBC, which is integer - so it is
        !passed straight through. It used to be routed through a local declared integer but
        !assigned from merge(.true., .false., ...), which is a logical-to-integer assignment that
        !only compiles as an Intel extension and that gfortran rejects outright.
        use_cdfmm = merge(.true., .false., useCDFMM /= 0)
        call loadMicroMagProblem( ntot, grid_n, grid_L, grid_type, u_ea, ProblemMode, solver, A0, Ms, K0, &
            gamma, alpha_mm, temperature, nt_Hext, Hext, nt, t, m0, dem_thres, useCuda, dem_appr, N_ret, N_file_out, &
            N_load, N_file_in, setTimeDis, nt_alpha, alphat, tol, thres, useCVODE, nt_conv, t_conv, &
            conv_tol, grid_pts, grid_ele, grid_nod, grid_nnod, exch_nval, exch_nrow, exch_val, exch_rows, &
            exch_cols, grid_abc, N_ave, &
            CV, useReturnHall, useAvgN, exch_weigh, exch_meth, exch_intpn, &
            n_macro, shiftVec, macroShape, sampleShape, exchPBC, &
            passExch, exch_ncols, CrysAxis, K0_arr, K1, K2, n_phase, phase_id, A_int, problem, dummy_run, &
            use_cdfmm, cdfmm_order, cdfmm_depth, cdfmm_basis, cdfmm_precision, &
            useDemag, rng_seed)

        if (hysteresis_solver .eq. 2) then
            problem%adaptiveHext = .true.
            problem%maxHextSteps = maxHextSteps
            problem%nHextAccepted = 0
            problem%H_start = H_start
            problem%H_end = H_end
            problem%dH_initial = dH_initial
            problem%dH_min = dH_min
            problem%dH_max = dH_max
            problem%dH_grow = dH_grow
            problem%dH_shrink = dH_shrink
            problem%dM_min = dM_min
            problem%dM_target = dM_target
            problem%dM_reject = dM_reject
            problem%switch_refine_dH = switch_refine_dH
            problem%use_switch_refine = (use_switch_refine /= 0)
        else if (hysteresis_solver .eq. 1) then
            problem%adaptiveHext = .false.
        else
            error stop 'hysteresis_solver must be 1 (static) or 2 (adaptive)'
        end if

        problem%min_tol = min_tol
        problem%min_maxiter = min_maxiter
        problem%min_maxrot = min_maxrot
        problem%min_fallback = min_fallback
        problem%min_saddle_check = min_saddle_check
        problem%min_predictor = min_predictor

        call SolveLandauLifshitzEquation( problem, solution )

        !The solver sizes these by the number of applied fields it allocated, which matches
        !nt_Hext_out for every calling mode; the min() only guards against a caller that sizes the
        !output differently.
        E_out = 0.
        n_feval = 0
        min_iter = 0
        min_torque = 0.
        min_status = 0
        min_eig = 0.
        n_copy = min( size(solution%n_feval), nt_Hext_out )
        E_out(:,1:n_copy,:) = solution%E_out(:,1:n_copy,:)
        n_feval(1:n_copy) = solution%n_feval(1:n_copy)
        min_iter(1:n_copy) = solution%min_iter(1:n_copy)
        min_torque(1:n_copy) = solution%min_torque(1:n_copy)
        min_status(1:n_copy) = solution%min_status(1:n_copy)
        min_eig(1:n_copy) = solution%min_eig(1:n_copy)


        t_out = solution%t_out
        M_mm = solution%M_out
        pts = solution%pts

        !solution%H_* are only allocated at the full (nt, ntot, nt_Hext, 3) size when
        !useReturnHall is set; otherwise they are a (1,1,1,3) placeholder. Copying the placeholder
        !into these full-size intent(out) arrays is a non-conforming array assignment, and it left
        !the caller holding uninitialised memory - values of order 1e+60 and 1e-227 were coming
        !back. Return explicit zeros in that case instead.
        if ( useReturnHall .eq. 1 ) then
            H_exc = solution%H_exc
            H_ext = solution%H_ext
            H_dem = solution%H_dem
            H_ani = solution%H_ani
        else
            H_exc = 0.
            H_ext = 0.
            H_dem = 0.
            H_ani = 0.
        endif
		n_Hext_accepted = problem%nHextAccepted
        !The exchange matrix in COO form only exists when the solver built the exchange operator
        !itself. With passExch the matrix comes from the caller and these arrays stay unallocated,
        !so return an empty matrix instead of reading them (the bounds-checked build stops here).
        if (.not. allocated(solution%gridinfo%Exch_mat_r)) then
            n_tot_Exch = 0
            ExchMat_r = 0
            ExchMat_c = 0
            ExchMat_v = 0.
            ExchMat_nr = 0
            ExchMat_nc = 0
        else
        n_tot_Exch = solution%gridinfo%Exch_mat_ntot

		if (exch_presize*ntot < n_tot_Exch) then
            write(*,*) 'ExchMat_presize is too small to copy all exchange matrix values. It is set to ', exch_presize*ntot, ' but the exchange matrix has ', n_tot_Exch, ' entries.'
            write(*,*) 'Please increase the value of exch_presize to at least ', n_tot_Exch/ntot, ' in the Python script.'
            write(*,*) 'Returning zeros for exchange matrix.'
            ExchMat_r = 0
            ExchMat_c = 0
            ExchMat_v = 0.
        else
            ExchMat_r(1:n_tot_Exch) = solution%gridinfo%Exch_mat_r
            ExchMat_c(1:n_tot_Exch) = solution%gridinfo%Exch_mat_c
            ExchMat_v(1:n_tot_Exch) = solution%gridinfo%Exch_mat_v
        end if

        ExchMat_nr = solution%gridinfo%Exch_mat_nr
		ExchMat_nc = solution%gridinfo%Exch_mat_nc
        end if
#else
        write(*,*) 'Compiled without micromagnetic part. Returning zeros.'
        n_tot_Exch = 0
        ExchMat_r = 0
        ExchMat_c = 0
        ExchMat_v = 0.
        ExchMat_nr = 0
        ExchMat_nc = 0

        t_out = 0.
        M_mm(:,:,:,:) = 0.
        pts(:,:) = 0.
        H_exc(:,:,:,:) = 0.
        H_ext(:,:,:,:) = 0.
        H_dem(:,:,:,:) = 0.
        H_ani(:,:,:,:) = 0.
        n_Hext_accepted = 0
        E_out = 0.
        n_feval = 0
        min_iter = 0
        min_torque = 0.
        min_status = 0
        min_eig = 0.
#endif


    !---------------------- finalize auxiliary modules -----------------------------
    call trace%finalize()
    call timer%log_finalize()
    !---------------------------------------------------------------------------------

    end subroutine RunMicroMagSimulation

end module FortranToPythonIO
