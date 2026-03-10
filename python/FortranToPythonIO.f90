module FortranToPythonIO
    use omp_mod
    use timer_mod
    use trace_mod
    use MagParameters
    use DemagFieldGetSolution
    use IterateMagnetSolution
#if USE_MICROMAG
    use MagTenseMicroMagPyIO
    use LandauLifshitzSolution
#endif
#if USE_FMM3D
    use fmm_nbor_tensor_mod
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

!#if USE_FMM3D
!    ! ---------- real implementation when USE_FMM3D=1 ----------
!    integer(8) :: nd, nsrc8, ntgt8, ier
!    real(8) :: fmm_eps, vol_i
!    real(8),allocatable :: source(:,:), targ(:,:), dipvec(:,:,:), pottarg(:,:), gradtarg(:,:,:)
!    integer :: i, j
!    !-------------------------- Interface to FMM3D dipole FMM --------------------------
!    ! Mainly included to allow f2py to parse the subroutine signature, potentially catching errors at compile time
!    interface
!      subroutine lfmm3d_t_d_g_vec(nd,eps,nsource,source,dipvec,ntarg,targ,pottarg,gradtarg,ier)
!        implicit none
!        integer(8), intent(in) :: nd, nsource, ntarg
!        real(8),    intent(in) :: eps
!        real(8) :: source(3,nsource), dipvec(nd,3,nsource), targ(3,ntarg)
!        real(8) :: pottarg(nd,ntarg), gradtarg(nd,3,ntarg)
!        integer(8) :: ier
!      end subroutine lfmm3d_t_d_g_vec
!    end interface
!    !---------------------------------------------------------------------------------
!    !------------------- define 4pi  -------------------------------------------------
!    fourpi = 12.566370614359172d0
!    !---------------------------------------------------------------------------------
!    !---------------- set FMM precision - 1e-6 if not provided -----------------------
!    fmm_eps = merge(eps, 1.0d-6, present(eps))
!    !---------------------------------------------------------------------------------
!    !---------------- allocate tmp arrays for FMM3D call -----------------------------
!    allocate(source(3,n_tiles), dipvec(1,3,n_tiles))
!    allocate(targ(3,n_pts), pottarg(1,n_pts), gradtarg(1,3,n_pts))
!    !---------------------------------------------------------------------------------
!    !---------------- rotate targets from (n_pts,3) to (3,n_pts) ---------------------
!    do j = 1, n_pts
!      targ(1,j)=pts(j,1)
!      targ(2,j)=pts(j,2)
!      targ(3,j)=pts(j,3)
!    end do
!    !------------------------------------------------------------------------------
!    !--------------- iterate over source tiles to get position and dipole moment -----
!    do i = 1, n_tiles
!      !---------------- get and rotate source position ----------------
!      source(1,i) = offset(i,1)
!      source(2,i) = offset(i,2)
!      source(3,i) = offset(i,3)
!      ! TODO - maybe take into account centerPos and dev_center as well?
!      !----------------------------------------------------------------
!      !------------------- get volume of tile -------------------------------
!      ! TODO - refine with proper volume calculations for other shapes
!      if (tileType(i) == 2) then
!       vol_i = max(0d0,tile_size(i,1))*max(0d0,tile_size(i,2))*max(0d0,tile_size(i,3))
!      else
!       vol_i = (max(0d0,tile_size(i,1))*max(0d0,tile_size(i,2))*max(0d0,tile_size(i,3)))**(1d0/3d0)
!       vol_i = vol_i**3
!      end if
!      !------------------------------------------------------------------------------------
!      !------------- convert magnetization to dipole moment --------------
!      dipvec(1,1,i) = Mag(i,1)*vol_i !* mu0
!      dipvec(1,2,i) = Mag(i,2)*vol_i !* mu0
!      dipvec(1,3,i) = Mag(i,3)*vol_i !* mu0
!      !--------------------------------------------------------------------
!    end do
!    !----------------------------------------------------------------------------
!    !-------------------  set integer variables for FMM3D call -------------------
!    nd=1_8
!    nsrc8=int(n_tiles,8)
!    ntgt8=int(n_pts,8)
!    !----------------------------------------------------------------------------
!    !-------------------- call FMM3D Laplace dipole FMM ------------------------
!    call lfmm3d_t_d_g_vec(nd,fmm_eps,nsrc8,source,dipvec,ntgt8,targ,pottarg,gradtarg,ier)
!    !----------------------------------------------------------------------------
!    !------------------ rotate gradient, scale and negate to get H-field in (n_pts,3) -------------------
!    do j = 1, n_pts
!      H(j,1) = -gradtarg(1,1,j) / fourpi
!      H(j,2) = -gradtarg(1,2,j) / fourpi
!      H(j,3) = -gradtarg(1,3,j) / fourpi
!    end do
!    !---------------------------------------------------------------------------------------------------
!    deallocate(source, dipvec, targ, pottarg, gradtarg)
!#else 
    !------------------ Fallback implementation when USE_FMM3D=0 ------------------
    print *, "WARNING: getHFromTilesFMM called but MagTense built without FMM3D support. - Returning zero field." 
    H(:,:) = 0.0d0
    return
    !----------------------------------------------------------------------------
!#endif

end subroutine getHFromTilesFMM

!----------------------------------------------------------------------------------
! FMM-backed H-field AT SOURCE POSITIONS using FMM3D (dipoles-only, src->src)
! Each tile i is a dipole at its centre with moment m_i = M_i * Volume_i.
! Returns H on the tile centres (n_tiles,3).
!
subroutine getHOnSourcesFMM( centerPos, dev_center, tile_size, vertices, Mag, u_ea, u_oa1, u_oa2, &
    mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType, stateFunctionIndex, &
    includeInIteration, exploitSymmetry, symmetryOps, Mrel, n_tiles, H_src, eps )
  implicit none
  !---------------- args ----------------
  integer(4),intent(in) :: n_tiles
  real(8),dimension(n_tiles,3),intent(in) :: centerPos, dev_center
  real(8),dimension(n_tiles,3),intent(in) :: tile_size
  real(8),dimension(n_tiles,3,4),intent(in) :: vertices
  real(8),dimension(n_tiles,3),intent(in) :: Mag, u_ea, u_oa1, u_oa2
  real(8),dimension(n_tiles),intent(in) :: mu_r_ea, mu_r_oa, Mrem
  integer(4),dimension(n_tiles),intent(in) :: tileType, magnetType, stateFunctionIndex
  integer(4),dimension(n_tiles),intent(in) :: includeInIteration, exploitSymmetry
  real(8),dimension(n_tiles,3),intent(in) :: offset, rotAngles, color, symmetryOps
  real(8),dimension(n_tiles),intent(in) :: Mrel
  real(8),dimension(n_tiles,3),intent(out) :: H_src
  real(8),intent(in),optional :: eps

  real(8) :: fourpi
!#if USE_FMM3D
!  !-------------- locals --------------
!  integer :: i, ier4
!  integer :: nd4, nsrc4
!  real(8) :: fmm_eps, vol_i
!  real(8),allocatable :: source(:,:)          ! (3, n_tiles)
!  real(8),allocatable :: dipvec(:,:,:)        ! (1, 3, n_tiles)
!  real(8),allocatable :: pot(:,:)             ! (1, n_tiles)
!  real(8),allocatable :: grad(:,:,:)          ! (1, 3, n_tiles)
!
!  interface
!    subroutine lfmm3d_s_d_g_vec(nd,eps,nsource,source,dipvec,pot,grad,ier)
!      implicit none
!      integer, intent(in) :: nd, nsource
!      double precision, intent(in) :: eps
!      double precision :: source(3,nsource)
!      double precision :: dipvec(nd,3,nsource)
!      double precision :: pot(nd,nsource)
!      double precision :: grad(nd,3,nsource)
!      integer :: ier
!    end subroutine lfmm3d_s_d_g_vec
!  end interface
!
!  fourpi  = 12.566370614359172d0
!  fmm_eps = merge(eps, 1.0d-6, present(eps))
!
!  allocate(source(3,n_tiles), dipvec(1,3,n_tiles))
!  allocate(pot(1,n_tiles), grad(1,3,n_tiles))
!
!  ! Pack sources (tile centres) and dipole moments (M*V)
!  do i = 1, n_tiles
!    ! Use offset as the global tile centre; adjust if your pipeline prefers centerPos+dev_center
!    source(1,i) = offset(i,1)
!    source(2,i) = offset(i,2)
!    source(3,i) = offset(i,3)
!
!    ! Volume estimate: exact for prisms; crude fallback for others
!    if (tileType(i) == 2) then
!      vol_i = max(0d0,tile_size(i,1))*max(0d0,tile_size(i,2))*max(0d0,tile_size(i,3))
!    else
!      vol_i = (max(0d0,tile_size(i,1))*max(0d0,tile_size(i,2))*max(0d0,tile_size(i,3)))**(1d0/3d0)
!      vol_i = vol_i**3
!    end if
!
!    dipvec(1,1,i) = Mag(i,1) * vol_i
!    dipvec(1,2,i) = Mag(i,2) * vol_i
!    dipvec(1,3,i) = Mag(i,3) * vol_i
!  end do
!
!  nd4   = 1
!  nsrc4 = n_tiles
!  call lfmm3d_s_d_g_vec(nd4, fmm_eps, nsrc4, source, dipvec, pot, grad, ier4)
!  if (ier4 /= 0) then
!    write(*,*) 'FMM3D error in getHOnSourcesFMM: ier =', ier4
!  end if
!
!  ! Map ∇u -> H: H = -grad / (4π)
!  do i = 1, n_tiles
!    H_src(i,1) = -grad(1,1,i) / fourpi
!    H_src(i,2) = -grad(1,2,i) / fourpi
!    H_src(i,3) = -grad(1,3,i) / fourpi
!  end do
!
!
!  call add_self_field(  H_src, centerPos, dev_center, tile_size, vertices, &
!                        Mag, u_ea, u_oa1, u_oa2, mu_r_ea, mu_r_oa, Mrem, tileType, &
!                        offset, rotAngles, color, magnetType, stateFunctionIndex,   &
!                        includeInIteration, exploitSymmetry, symmetryOps, Mrel)
!
!
!    call add_neighbour_corrections(H_src, centerPos, dev_center, tile_size, vertices, Mag, u_ea, u_oa1, u_oa2, &
!                                mu_r_ea, mu_r_oa, Mrem, tileType, offset, rotAngles, color, magnetType,     &
!                                stateFunctionIndex, includeInIteration, exploitSymmetry, symmetryOps, Mrel, &
!                                radius_cells = 2)
!
!  deallocate(source, dipvec, pot, grad)
!#else
!  fourpi = 12.566370614359172d0
  H_src(:,:) = 0.0d0
  print *, "WARNING: getHOnSourcesFMM called without FMM3D support; returning zeros."
!#endif
end subroutine getHOnSourcesFMM


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
        K1, K2, K0_arr, CrysAxis, gamma, alpha_mm, MaxT0, nt_Hext, nt_Hext_out, Hext, nt, t, m0, dem_thres, useCuda, dem_appr, N_ret, N_file_out, &
        N_load, N_file_in, setTimeDis, nt_alpha, alphat, tol, thres, useCVODE, nt_conv, t_conv, &
        conv_tol, grid_pts, grid_ele, grid_nod, grid_nnod, exch_nval, exch_nrow, exch_val, exch_rows, &
        exch_rowe, exch_col, grid_abc, usePrecision, nThreadsMatlab, N_ave, CV, useReturnHall, demigstp, & 
		exch_weigh, exch_meth, exch_intpn, passExch, exch_ncols, exch_presize, &
        t_out, M_mm, pts, H_exc, H_ext, H_dem, H_ani, &
		n_tot_Exch, ExchMat_r, ExchMat_c, ExchMat_v, ExchMat_nr, ExchMat_nc, dummy_run, fmm_cells_per_node, eps_fmm, ifunif, nlmin, nlmax, allow_fmm_short_circuit, fmm_min_n, &
        log_dir,timer_log_file, trace_log_file, window_enabled, window_interval, trace_enabled, flush_each, trace_verbose )

        integer(4), intent(in) :: ntot, nt_conv, grid_type, nt_Hext, nt_alpha, nt, grid_nnod, exch_nval, exch_nrow, exch_ncols, exch_presize
        integer(4), intent(in) :: nt_Hext_out
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
        real(8),dimension(nt,ntot,nt_Hext_out,3),intent(out) :: M_mm
        real(8),dimension(nt,ntot,nt_Hext_out,3),intent(out) :: H_exc, H_ext, H_dem, H_ani
        real(8),dimension(ntot,3),intent(out) :: pts
		
		integer,intent(out) :: n_tot_Exch
		integer,dimension(exch_presize*ntot),intent(out)  :: ExchMat_r
		integer,dimension(exch_presize*ntot),intent(out)  :: ExchMat_c
		real(8),dimension(exch_presize*ntot),intent(out)  :: ExchMat_v
		integer,intent(out) :: ExchMat_nr,ExchMat_nc
        
        integer, intent(in) :: dummy_run
        integer(4), intent(in) :: fmm_cells_per_node
        real(8), intent(in) :: eps_fmm
        integer(4), intent(in) :: ifunif
        integer(4), intent(in) :: nlmin
        integer(4), intent(in) :: nlmax
        integer(4), intent(in) :: allow_fmm_short_circuit
        integer(4), intent(in) :: fmm_min_n

        !-------------------- timer and trace modules --------------------------------------
        character*256,intent(in) :: timer_log_file, trace_log_file, log_dir
        integer, intent(in) :: window_enabled, trace_enabled, flush_each
        real(8), intent(in) :: window_interval
        integer, intent(in) :: trace_verbose
        
        logical :: window_enabled_l, trace_enabled_l, flush_each_l

        !-----------------------------------------------------------------------------------

#if USE_MICROMAG
        type(MicroMagProblem) :: problem
        type(MicroMagSolution) :: solution



        !---------------------- initiaize auxiliary modules -----------------------------
        call omp%init()
        call omp%info()


        window_enabled_l = merge(.true., .false., window_enabled /= 0)
        trace_enabled_l = merge(.true., .false., trace_enabled /= 0)
        flush_each_l = merge(.true., .false., flush_each /= 0)

        call timer%log_init(trim(log_dir), trim(timer_log_file), window_enabled=window_enabled_l, window_interval=window_interval)
        call trace%trace_init(trim(log_dir), trim(trace_log_file), enabled=trace_enabled_l, unit=97, flush_each=flush_each_l, verbose=trace_verbose)
        
        !call timer%log_init("timing.log", window_enabled=.true., window_interval=30.0d0)
        !call trace%init("trace.log", enabled=.false., unit=97, flush_each=.true.)
        !---------------------------------------------------------------------------------



        call loadMicroMagProblem( ntot, grid_n, grid_L, grid_type, u_ea, ProblemMode, solver, A0, Ms, K0, &
            gamma, alpha_mm, MaxT0, nt_Hext, Hext, nt, t, m0, dem_thres, useCuda, dem_appr, N_ret, N_file_out, &
            N_load, N_file_in, setTimeDis, nt_alpha, alphat, tol, thres, useCVODE, nt_conv, t_conv, &
            conv_tol, grid_pts, grid_ele, grid_nod, grid_nnod, exch_nval, exch_nrow, exch_val, exch_rows, &
            exch_rowe, exch_col, grid_abc, usePrecision, nThreadsMatlab, N_ave, &
			CV, useReturnHall, demigstp, exch_weigh, exch_meth, exch_intpn,	passExch, exch_ncols, &
            CrysAxis, K0_arr, K1, K2, problem, dummy_run, fmm_cells_per_node, eps_fmm, ifunif, nlmin, nlmax, allow_fmm_short_circuit, fmm_min_n)

        print *, " starting SolveLandauLifshitzEquation "
        call SolveLandauLifshitzEquation( problem, solution )
        print *, " finished SolveLandauLifshitzEquation "


        t_out = solution%t_out
        print *, " M_mm", shape(M_mm)
        print *, " M_mm size", shape(solution%M_out)
        M_mm = solution%M_out
        print *, " pts ", shape(pts)
        print *, " pts size", shape(solution%pts)
        pts = solution%pts
        print *, " H_exc size", shape(solution%H_exc)
        H_exc = solution%H_exc
        print *, " H_ext size", shape(solution%H_ext)
        H_ext = solution%H_ext
        print *, " H_dem size", shape(solution%H_dem)
        H_dem = solution%H_dem
        print *, " H_ani size", shape(solution%H_ani)
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
#endif


    !---------------------- finalize auxiliary modules -----------------------------
    call trace%finalize()
    call timer%log_finalize()
    !---------------------------------------------------------------------------------

    end subroutine RunMicroMagSimulation













    !================== Helper functions for FMM correction ==========================================

        !---------------- for getting self field correction ----------------------------------------

!    subroutine add_self_field(H_fmm, centerPos, dev_center, tile_size, vertices, &
!                          Mag, u_ea, u_oa1, u_oa2, mu_r_ea, mu_r_oa, Mrem, tileType, &
!                          offset, rotAngles, color, magnetType, stateFunctionIndex, &
!                          includeInIteration, exploitSymmetry, symmetryOps, Mrel)
!        implicit none
!        ! ---- dummy arguments (assumed-shape everywhere) ----
!        real(8), intent(inout) :: H_fmm(:,:)               ! (N,3)
!        real(8), intent(in)    :: centerPos(:,:), dev_center(:,:), tile_size(:,:), vertices(:,:,:)
!        real(8), intent(in)    :: Mag(:,:), u_ea(:,:), u_oa1(:,:), u_oa2(:,:)
!        real(8), intent(in)    :: mu_r_ea(:), mu_r_oa(:), Mrem(:)
!        integer(4), intent(in) :: tileType(:), magnetType(:), stateFunctionIndex(:)
!        real(8), intent(in)    :: offset(:,:), rotAngles(:,:), color(:,:), symmetryOps(:,:), Mrel(:)
!        integer(4), intent(in) :: includeInIteration(:), exploitSymmetry(:)
!
!        ! ---- locals ----
!        integer(4) :: n_all, t
!        real(8)    :: pts(1,3)
!        real(8)    :: Ktt(1,1,3,3)
!        real(8)    :: Kloc(3,3)
!
!        ! optional sanity (defensive): ensure second dim==3 for relevant arrays
!        ! if (size(H_fmm,2) /= 3) stop 'add_self_field: H_fmm second dim must be 3'
!        ! if (size(Mag,2)   /= 3) stop 'add_self_field: Mag second dim must be 3'
!
!        n_all = size(H_fmm,1)
!
!        do t = 1, n_all
!            pts(1,:) = offset(t,:)  ! evaluate at tile t centre
!
!            call getNFromTiles( centerPos(t:t,:), dev_center(t:t,:), tile_size(t:t,:), vertices(t:t,:,:), &
!                                Mag(t:t,:), u_ea(t:t,:), u_oa1(t:t,:), u_oa2(t:t,:),                       &
!                                mu_r_ea(t:t), mu_r_oa(t:t), Mrem(t:t), tileType(t:t),                      &
!                                offset(t:t,:), rotAngles(t:t,:), color(t:t,:), magnetType(t:t),            &
!                                stateFunctionIndex(t:t), includeInIteration(t:t),                          &
!                                exploitSymmetry(t:t), symmetryOps(t:t,:), Mrel(t:t),                       &
!                                pts, n_tiles=1, n_pts=1, N=Ktt )
!
!            ! Slice directly: Ktt(1,1,:,:) is already 3x3 (no RESHAPE needed)
!            Kloc      = Ktt(1,1,:,:)
!            H_fmm(t,:) = H_fmm(t,:) + matmul(Kloc, Mag(t,:))
!        end do
!    end subroutine add_self_field
!        !------------------------------------------------------------------------------------------------
!
!subroutine add_neighbour_corrections(H_fmm, centerPos, dev_center, tile_size, vertices, &
!                                     Mag, u_ea, u_oa1, u_oa2, mu_r_ea, mu_r_oa, Mrem, tileType, &
!                                     offset, rotAngles, color, magnetType, stateFunctionIndex,   &
!                                     includeInIteration, exploitSymmetry, symmetryOps, Mrel,     &
!                                     radius_cells)
!  implicit none
!  ! in/out
!  real(8), intent(inout) :: H_fmm(:,:)                 ! (N,3)
!
!  ! inputs
!  real(8), intent(in)    :: centerPos(:,:), dev_center(:,:), tile_size(:,:), vertices(:,:,:)
!  real(8), intent(in)    :: Mag(:,:), u_ea(:,:), u_oa1(:,:), u_oa2(:,:)
!  real(8), intent(in)    :: mu_r_ea(:), mu_r_oa(:), Mrem(:)
!  integer(4), intent(in) :: tileType(:), magnetType(:), stateFunctionIndex(:)
!  real(8), intent(in)    :: offset(:,:), rotAngles(:,:), color(:,:), symmetryOps(:,:), Mrel(:)
!  integer(4), intent(in) :: includeInIteration(:), exploitSymmetry(:)
!  integer(4), intent(in) :: radius_cells
!
!  ! locals
!  integer(4) :: N, t, j, jj, n_nbr
!  integer(4), allocatable :: nbr_idx(:)
!  real(8) :: pitch(3), pts(1,3)
!  real(8), allocatable :: Ktn(:,:,:,:)   ! (1, n_nbr, 3, 3)
!  real(8) :: Rvec(3), Kdip(3,3), Kloc(3,3), Vol_j
!  real(8) :: dH(3)
!#if USE_FMM3D
!
!  N = size(H_fmm,1)
!
!  ! simple pitch estimate (mean size per axis)
!  pitch(1) = sum(tile_size(:,1))/real(N,8)
!  pitch(2) = sum(tile_size(:,2))/real(N,8)
!  pitch(3) = sum(tile_size(:,3))/real(N,8)
!
!  do t = 1, N
!    ! ---- pass 1: count neighbours ----
!    n_nbr = 0
!    do j = 1, N
!      if ( is_neighbour(t, j, offset, tile_size, pitch, radius_cells) ) n_nbr = n_nbr + 1
!    end do
!    if (n_nbr == 0) cycle
!
!    ! ---- pass 2: allocate and fill index list exactly ----
!    allocate(nbr_idx(n_nbr))
!    n_nbr = 0
!    do j = 1, N
!      if ( is_neighbour(t, j, offset, tile_size, pitch, radius_cells) ) then
!        n_nbr = n_nbr + 1
!        nbr_idx(n_nbr) = j
!      end if
!    end do
!
!    allocate(Ktn(1, n_nbr, 3, 3))
!
!    ! exact prism tensors for these neighbours at target centre
!    pts(1,:) = offset(t,:)
!    call getNFromTiles( centerPos(nbr_idx,:), dev_center(nbr_idx,:), tile_size(nbr_idx,:), vertices(nbr_idx,:,:), &
!                        Mag(nbr_idx,:), u_ea(nbr_idx,:), u_oa1(nbr_idx,:), u_oa2(nbr_idx,:),                       &
!                        mu_r_ea(nbr_idx), mu_r_oa(nbr_idx), Mrem(nbr_idx), tileType(nbr_idx),                      &
!                        offset(nbr_idx,:), rotAngles(nbr_idx,:), color(nbr_idx,:), magnetType(nbr_idx),            &
!                        stateFunctionIndex(nbr_idx), includeInIteration(nbr_idx),                                   &
!                        exploitSymmetry(nbr_idx), symmetryOps(nbr_idx,:), Mrel(nbr_idx),                           &
!                        pts, n_tiles = n_nbr, n_pts = 1, N = Ktn )
!
!    ! accumulate correction ΔH_t = Σ (K_prism - K_dip*V) @ M_j
!    dH = 0.0d0
!    do jj = 1, n_nbr
!      j      = nbr_idx(jj)
!      Rvec   = pts(1,:) - offset(j,:)                   ! source -> target
!      call dipole_tensor_3x3(Rvec, Kdip)                ! per-unit-M
!      Vol_j  = tile_size(j,1) * tile_size(j,2) * tile_size(j,3)
!      Kdip   = Kdip * Vol_j                             ! FMM used m = M * V
!      Kloc   = Ktn(1, jj, :, :)                         ! exact prism tensor
!      dH     = dH + matmul( Kloc - Kdip, Mag(j,:) )
!    end do
!
!    H_fmm(t,:) = H_fmm(t,:) + dH
!
!    deallocate(Ktn)
!    deallocate(nbr_idx)
!  end do
!#endif
!end subroutine add_neighbour_corrections

        !------------------------------------------------------------------------------------------------
    !============================================================================================================

end module FortranToPythonIO
