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

 

!----------------------------------------------------------------------------------
! FMM-backed H-field using FMM3D Laplace dipoles
! Approximates each tile as a single dipole located at its centre, with
! dipole moment = Mag(i,:) * Volume(i).
! Currently volume is estimated from tile_size for prisms; other shapes use a crude
! fallback volume from the geometric mean of tile_size. Refine as needed.
!
subroutine getHFromTilesFMM( centerPos, dev_center, tile_size, vertices, Mag, u_ea, u_oa1, u_oa2, &
    mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
    includeInIteration, exploitSymmetry, symmetryOps, Mrel, pts, n_tiles, n_pts, H, eps )
    implicit none
    integer(4),intent(in) :: n_tiles, n_pts
    real(8),dimension(n_tiles,3),intent(in) :: centerPos, dev_center
    real(8),dimension(n_tiles,3),intent(in) :: tile_size
    real(8),dimension(n_tiles,3,4),intent(in) :: vertices
    real(8),dimension(n_tiles,3),intent(in) :: Mag, u_ea, u_oa1, u_oa2
    real(8),dimension(n_tiles),intent(in) :: mu_r_ea, mu_r_oa, Mrem
    integer(4),dimension(n_tiles),intent(in) :: tileType, magnetType, stateFunctionIndex
    integer(4),dimension(n_tiles),intent(in) :: includeInIteration, exploitSymmetry
    real(8),dimension(n_tiles,3),intent(in) :: offset, rotAngles, color, symmetryOps
    real(8),dimension(n_tiles),intent(in) :: Mrel
    real(8),dimension(n_pts,3),intent(in) :: pts
    real(8),dimension(n_pts,3),intent(out) :: H
    real(8),intent(in),optional :: eps
    real(8) :: fourpi

#if USE_FMM3D
    ! ---------- real implementation when USE_FMM3D=1 ----------
    integer(8) :: nd, nsrc8, ntgt8, ier
    real(8) :: fmm_eps, vol_i
    real(8),allocatable :: source(:,:), targ(:,:), dipvec(:,:,:), pottarg(:,:), gradtarg(:,:,:)
    integer :: i, j
    !-------------------------- Interface to FMM3D dipole FMM --------------------------
    ! Mainly included to allow f2py to parse the subroutine signature, potentially catching errors at compile time
    interface
      subroutine lfmm3d_t_d_g_vec(nd,eps,nsource,source,dipvec,ntarg,targ,pottarg,gradtarg,ier)
        implicit none
        integer(8), intent(in) :: nd, nsource, ntarg
        real(8),    intent(in) :: eps
        real(8) :: source(3,nsource), dipvec(nd,3,nsource), targ(3,ntarg)
        real(8) :: pottarg(nd,ntarg), gradtarg(nd,3,ntarg)
        integer(8) :: ier
      end subroutine lfmm3d_t_d_g_vec
    end interface
    !---------------------------------------------------------------------------------
    !------------------- define 4pi  -------------------------------------------------
    fourpi = 12.566370614359172d0
    !---------------------------------------------------------------------------------
    !---------------- set FMM precision - 1e-6 if not provided -----------------------
    fmm_eps = merge(eps, 1.0d-6, present(eps))
    !---------------------------------------------------------------------------------
    !---------------- allocate tmp arrays for FMM3D call -----------------------------
    allocate(source(3,n_tiles), dipvec(1,3,n_tiles))
    allocate(targ(3,n_pts), pottarg(1,n_pts), gradtarg(1,3,n_pts))
    !---------------------------------------------------------------------------------
    !---------------- rotate targets from (n_pts,3) to (3,n_pts) ---------------------
    do j = 1, n_pts
      targ(1,j)=pts(j,1)
      targ(2,j)=pts(j,2)
      targ(3,j)=pts(j,3)
    end do
    !------------------------------------------------------------------------------
    !--------------- iterate over source tiles to get position and dipole moment -----
    do i = 1, n_tiles
      !---------------- get and rotate source position ----------------
      source(1,i) = offset(i,1)
      source(2,i) = offset(i,2)
      source(3,i) = offset(i,3)
      ! TODO - maybe take into account centerPos and dev_center as well?
      !----------------------------------------------------------------
      !------------------- get volume of tile -------------------------------
      ! TODO - refine with proper volume calculations for other shapes
      if (tileType(i) == 2) then
       vol_i = max(0d0,tile_size(i,1))*max(0d0,tile_size(i,2))*max(0d0,tile_size(i,3))
      else
       vol_i = (max(0d0,tile_size(i,1))*max(0d0,tile_size(i,2))*max(0d0,tile_size(i,3)))**(1d0/3d0)
       vol_i = vol_i**3
      end if
      !------------------------------------------------------------------------------------
      !------------- convert magnetization to dipole moment --------------
      dipvec(1,1,i) = Mag(i,1)*vol_i !* mu0
      dipvec(1,2,i) = Mag(i,2)*vol_i !* mu0
      dipvec(1,3,i) = Mag(i,3)*vol_i !* mu0
      !--------------------------------------------------------------------
    end do
    !----------------------------------------------------------------------------
    !-------------------  set integer variables for FMM3D call -------------------
    nd=1_8
    nsrc8=int(n_tiles,8)
    ntgt8=int(n_pts,8)
    !----------------------------------------------------------------------------
    !-------------------- call FMM3D Laplace dipole FMM ------------------------
    call lfmm3d_t_d_g_vec(nd,fmm_eps,nsrc8,source,dipvec,ntgt8,targ,pottarg,gradtarg,ier)
    !----------------------------------------------------------------------------
    !------------------ rotate gradient, scale and negate to get H-field in (n_pts,3) -------------------
    do j = 1, n_pts
      H(j,1) = -gradtarg(1,1,j) / fourpi
      H(j,2) = -gradtarg(1,2,j) / fourpi
      H(j,3) = -gradtarg(1,3,j) / fourpi
    end do
    !---------------------------------------------------------------------------------------------------
    deallocate(source, dipvec, targ, pottarg, gradtarg)
#else 
    !------------------ Fallback implementation when USE_FMM3D=0 ------------------
    print *, "WARNING: getHFromTilesFMM called but MagTense built without FMM3D support. - Returning zero field." 
    H(:,:) = 0.0d0
    return
    !----------------------------------------------------------------------------
#endif

end subroutine getHFromTilesFMM



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
        K1, K2, K0_arr, CrysAxis, gamma, alpha_mm, MaxT0, nt_Hext, Hext, nt, t, m0, dem_thres, useCuda, dem_appr, N_ret, N_file_out, &
        N_load, N_file_in, setTimeDis, nt_alpha, alphat, tol, thres, useCVODE, nt_conv, t_conv, &
        conv_tol, grid_pts, grid_ele, grid_nod, grid_nnod, exch_nval, exch_nrow, exch_val, exch_rows, &
        exch_rowe, exch_col, grid_abc, usePrecision, nThreadsMatlab, N_ave, CV, useReturnHall, demigstp, & 
		exch_weigh, exch_meth, exch_intpn, passExch, exch_ncols, exch_presize, &
        t_out, M_mm, pts, H_exc, H_ext, H_dem, H_ani, &
		n_tot_Exch, ExchMat_r, ExchMat_c, ExchMat_v, ExchMat_nr, ExchMat_nc)

        integer(4), intent(in) :: ntot, nt_conv, grid_type, nt_Hext, nt_alpha, nt, grid_nnod, exch_nval, exch_nrow, exch_ncols, exch_presize
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
        integer(4),dimension(exch_nval),intent(in) :: exch_rows
		integer(4),dimension(exch_nrow),intent(in) :: exch_rowe
        real(8),dimension(nt_conv),intent(in) :: t_conv
		integer(4),intent(in) :: ProblemMode, solver, useCuda, dem_appr, usePrecision, nThreadsMatlab
		integer(4),intent(in) :: N_ret, N_load, setTimeDis, useCVODE, useReturnHall, demigstp, exch_meth, exch_intpn, passExch
        real(8),intent(in) :: gamma, alpha_mm, MaxT0, tol, thres, conv_tol, dem_thres
		real(8),dimension(ntot),intent(in) :: A0, Ms, K0, K1, K2
        real(8),dimension(ntot,6,3),intent(in) :: K0_arr
        real(8),dimension(ntot,3,3),intent(in):: CrysAxis
        character*256,intent(in) :: N_file_in, N_file_out
		real(8), intent(in) :: CV, exch_weigh

        real(8),dimension(nt),intent(in) :: t
        real(8),dimension(nt),intent(out) :: t_out
        real(8),dimension(nt,ntot,1,3),intent(out) :: M_mm
        real(8),dimension(nt,ntot,1,3),intent(out) :: H_exc, H_ext, H_dem, H_ani
        real(8),dimension(ntot,3),intent(out) :: pts
		
		integer,intent(out) :: n_tot_Exch
		integer,dimension(exch_presize*ntot),intent(out)  :: ExchMat_r
		integer,dimension(exch_presize*ntot),intent(out)  :: ExchMat_c
		real(8),dimension(exch_presize*ntot),intent(out)  :: ExchMat_v
		integer,intent(out) :: ExchMat_nr,ExchMat_nc

#if USE_MICROMAG
        type(MicroMagProblem) :: problem
        type(MicroMagSolution) :: solution

        call loadMicroMagProblem( ntot, grid_n, grid_L, grid_type, u_ea, ProblemMode, solver, A0, Ms, K0, &
            gamma, alpha_mm, MaxT0, nt_Hext, Hext, nt, t, m0, dem_thres, useCuda, dem_appr, N_ret, N_file_out, &
            N_load, N_file_in, setTimeDis, nt_alpha, alphat, tol, thres, useCVODE, nt_conv, t_conv, &
            conv_tol, grid_pts, grid_ele, grid_nod, grid_nnod, exch_nval, exch_nrow, exch_val, exch_rows, &
            exch_rowe, exch_col, grid_abc, usePrecision, nThreadsMatlab, N_ave, &
			CV, useReturnHall, demigstp, exch_weigh, exch_meth, exch_intpn,	passExch, exch_ncols, &
            CrysAxis, K0_arr, K1, K2, problem )

        call SolveLandauLifshitzEquation( problem, solution )

        t_out = solution%t_out
        M_mm = solution%M_out
        pts = solution%pts
        H_exc = solution%H_exc
        H_ext = solution%H_ext
        H_dem = solution%H_dem
        H_ani = solution%H_ani
		
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
		
		write(*,*) 'Done saving returned variables'
		
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
