module FortranToPythonIO
    use auxInit_mod
    use timer_mod
    use trace_mod
    use TileNComponents
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
    includeInIteration, exploitSymmetry, symmetryOps, Mrel, pts, n_tiles, n_pts, & 
    fmm_eps, fmm_nterms_in, fmm_cells_per_node, fmm_nlmin, fmm_nlmax, fmm_ifunif,  &
    do_target, do_fi, &
    H )
    !------------ IMPORTANT NOTE ------------
    ! n_pts MUST match the number of output points. 
    ! if do_target then it is len(pts) - if !do_target then it is n_tiles (for direct eval of sources as targets)
    !---------------------------------------------

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
    real(8),intent(in) :: fmm_eps
    integer(4),intent(in) :: fmm_nterms_in, fmm_cells_per_node, fmm_nlmin, fmm_nlmax, fmm_ifunif
    integer(4) :: do_target, do_FI
    real(8) :: fourpi


#if USE_FMM3D
    real(8) , contiguous, pointer:: source(:,:), dipvec(:,:,:), grad(:,:,:)
    real(8) , contiguous, pointer:: dipvec_s(:,:), grad_s(:,:), pot(:), pottarg(:), gradtarg(:,:),   targ(:,:)
    class(FMM3DTree), pointer :: fmm_tree
    integer :: nd, ier
    integer :: i, j
    integer :: nterms_in, cells_per_node, nlmin, nlmax, ifunif
    real(8) :: eps
    real(8) :: vol_i
    integer(8) :: nsource, ier8, ntarg

    interface 
    subroutine lfmm3d_t_d_g(eps,nsource,source,dipvec,ntarg,targ,pottarg,gradtarg,ier)
        implicit none
        ! FMM3D is built without /assume:underscore /names:lowercase (see
        ! external/FMM3D/make.inc.windows), so on Windows it exports the
        ! default upper-case name while this file, compiled with those
        ! options, would otherwise look for lfmm3d_t_d_g_. On Linux both
        ! sides already agree on the default lower-case form.
#ifdef _WIN32
        !DEC$ ATTRIBUTES ALIAS:"LFMM3D_T_D_G" :: lfmm3d_t_d_g
#endif
        integer(8), intent(in) :: nsource, ntarg
        real(8),    intent(in) :: eps
        real(8) :: source(3,nsource), dipvec(3,nsource), targ(3,ntarg)
        real(8) :: pottarg(ntarg), gradtarg(3,ntarg)
        integer(8) :: ier
    end subroutine lfmm3d_t_d_g
    subroutine lfmm3d_s_d_g(eps,nsource,source,dipvec,pot,grad,ier)
        implicit none
#ifdef _WIN32
        !DEC$ ATTRIBUTES ALIAS:"LFMM3D_S_D_G" :: lfmm3d_s_d_g
#endif
        integer(8), intent(in) :: nsource
        real(8),    intent(in) :: eps
        real(8) :: source(3,nsource), dipvec(3,nsource)
        real(8) :: pot(nsource), grad(3,nsource)
        integer(8) :: ier
    end subroutine lfmm3d_s_d_g
    end interface


    ier = 0
    ier8 = 0
    nsource = n_tiles
    ntarg = n_pts

    eps = fmm_eps
    nterms_in = fmm_nterms_in
    cells_per_node = fmm_cells_per_node
    nlmin = fmm_nlmin
    nlmax = fmm_nlmax
    ifunif = fmm_ifunif
   !------------------- define 4pi  -------------------------------------------------
   fourpi = 12.566370614359172d0
   !---------------------------------------------------------------------------------
   !---------------- allocate tmp arrays for FMM3D call -----------------------------
   allocate(source(3,n_tiles), dipvec(1,3,n_tiles), grad(1,3,n_tiles))
   allocate(targ(3,n_pts), pottarg(n_pts), gradtarg(3,n_pts), pot(n_tiles))
   dipvec_s => dipvec(1,:,:)
   grad_s => grad(1,:,:)

   grad (:,:,:) = 0.0d0
   !---------------------------------------------------------------------------------


   H(:,:) = 0.0d0


   do i = 1, n_tiles
      !---------------- get and rotate source position ----------------
      source(1,i) = offset(i,1)
      source(2,i) = offset(i,2)
      source(3,i) = offset(i,3)
      ! TODO - maybe take into account centerPos and dev_center as well?
      !----------------------------------------------------------------
      !------------------- get volume of tile -------------------------------
        vol_i = real(tile_size(i,1), DP) * &
                real(tile_size(i,2), DP) * &
                real(tile_size(i,3), DP)
 
      !------------------------------------------------------------------------------------
      !------------- convert magnetization to dipole moment --------------
      dipvec(1,1,i) = Mag(i,1) * vol_i !* Mrem(i)
      dipvec(1,2,i) = Mag(i,2) * vol_i !* Mrem(i)
      dipvec(1,3,i) = Mag(i,3) * vol_i !* Mrem(i)
   end do

   do i = 1, n_pts
      targ(1,i)=pts(i,1)
      targ(2,i)=pts(i,2)
      targ(3,i)=pts(i,3)
   end do


   if (do_FI .eq. 1) then
    if (do_target .eq. 1) then
        call lfmm3d_t_d_g(eps, nsource, source, dipvec_s, ntarg, targ, pottarg, gradtarg, ier8)
    else
        call lfmm3d_s_d_g(eps, nsource, source, dipvec_s, pot, grad_s, ier8)
    end if
   else
        allocate(fmm_tree)
        nd = 1
        fmm_tree%nterms_in = nterms_in
        call fmm_tree%build_tree( source, eps, cells_per_node, ier, ifunif, nlmin, nlmax)
        call fmm_tree%make_and_eval(dipvec, grad)
        !missing direvt eval
        !dealloc releases the source array as well, since build_tree takes ownership of it, so
        !drop our own reference to it here rather than deallocating it a second time below.
        call fmm_tree%dealloc()
        deallocate(fmm_tree)
        nullify(source)
   end if



   if (do_FI .eq. 1 .and. do_target .eq. 1)then 
        do j = 1, n_pts
            H(j,1) = - gradtarg(1,j) / fourpi
            H(j,2) = - gradtarg(2,j) / fourpi
            H(j,3) = - gradtarg(3,j) / fourpi
        end do
   else
        do j = 1, n_pts
            H(j,1) = grad(1,1,j) / fourpi
            H(j,2) = grad(1,2,j) / fourpi
            H(j,3) = grad(1,3,j) / fourpi
        end do
    end if

    if (associated(source)) deallocate(source)
    deallocate(dipvec, grad)
    deallocate(targ, pottarg, gradtarg, pot)
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
        K1, K2, K0_arr, CrysAxis, gamma, alpha_mm, temperature, MaxT0, nt_Hext, nt_Hext_out, Hext, nt, t, m0, dem_thres, useCuda, dem_appr, N_ret, N_file_out, &
        N_load, N_file_in, setTimeDis, nt_alpha, alphat, tol, thres, useCVODE, nt_conv, t_conv, &
        conv_tol, grid_pts, grid_ele, grid_nod, grid_nnod, exch_nval, exch_nrow, exch_val, exch_rows, &
        exch_cols, grid_abc, usePrecision, nThreadsMatlab, N_ave, CV, useReturnHall, useAvgN, demigstp, & 
		exch_weigh, exch_meth, exch_intpn, passExch, exch_ncols, exch_presize, &
        n_macro, shiftVec, macroShape, sampleShape, exchPBC, hysteresis_solver, &
        H_start, H_end, dH_initial, dH_min, dH_max, maxHextSteps, dM_min, dM_target, dM_reject, dH_grow, dH_shrink, switch_refine_dH, use_switch_refine, &
        t_out, M_mm, pts, H_exc, H_ext, H_dem, H_ani, n_Hext_accepted, &
		n_tot_Exch, ExchMat_r, ExchMat_c, ExchMat_v, ExchMat_nr, ExchMat_nc, dummy_run, fmm_cells_per_node, eps_fmm, ifunif, nlmin, nlmax, allow_fmm_short_circuit, fmm_min_n, fmm_nterms, useFMM, &
        log_dir,timer_log_file, trace_log_file, window_enabled, window_interval, trace_enabled, flush_each, trace_verbose, useDemag, rng_seed )

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
		integer(4),intent(in) :: ProblemMode, solver, useCuda, dem_appr, usePrecision, nThreadsMatlab, useAvgN
		integer(4),intent(in) :: N_ret, N_load, setTimeDis, useCVODE, useReturnHall, demigstp, exch_meth, exch_intpn, passExch, useDemag
        real(8),intent(in) :: gamma, alpha_mm, MaxT0, tol, thres, conv_tol, dem_thres
		real(8),dimension(ntot),intent(in) :: A0, Ms, K0, K1, K2, temperature
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
        integer(4), intent(in) :: fmm_nterms
        integer(4), intent(in) :: useFMM
        !> Seed for the stochastic thermal field: 0 keeps the compiler default sequence (identical
        !> on every run), a positive value seeds deterministically, a negative value seeds from the
        !> clock so that repeated runs are independent Monte-Carlo samples.
        integer(4), intent(in) :: rng_seed

        !-------------------- timer and trace modules --------------------------------------
        character*256,intent(in) :: timer_log_file, trace_log_file, log_dir
        integer, intent(in) :: window_enabled, trace_enabled, flush_each
        real(8), intent(in) :: window_interval
        integer, intent(in) :: trace_verbose
        !-----------------------------------------------------------------------------------
        
        logical :: window_enabled_l, trace_enabled_l, flush_each_l
        logical :: use_fmm
        !! Local auxiliary instance - avoids referencing the auxInit_mod module
        !! global across the f2py name-mangling boundary on Windows.
        type(auxInit_t) :: auxInit_local

#if USE_MICROMAG
        type(MicroMagProblem) :: problem
        type(MicroMagSolution) :: solution



        !---------------------- initiaize auxiliary modules -----------------------------
        call initAux(auxInit_local, log_dir, timer_log_file, trace_log_file, window_enabled, &
            window_interval, trace_enabled, flush_each, trace_verbose)
        !---------------------------------------------------------------------------------




        !exchPBC is an integer flag per direction at both ends of this call - loadMicroMagProblem
        !takes integer(4) and stores it in problem%macrogrid%exchPBC, which is integer - so it is
        !passed straight through. It used to be routed through a local declared integer but
        !assigned from merge(.true., .false., ...), which is a logical-to-integer assignment that
        !only compiles as an Intel extension and that gfortran rejects outright.
        use_fmm = merge(.true., .false., useFMM /= 0)
        call loadMicroMagProblem( ntot, grid_n, grid_L, grid_type, u_ea, ProblemMode, solver, A0, Ms, K0, &
            gamma, alpha_mm, temperature, MaxT0, nt_Hext, Hext, nt, t, m0, dem_thres, useCuda, dem_appr, N_ret, N_file_out, &
            N_load, N_file_in, setTimeDis, nt_alpha, alphat, tol, thres, useCVODE, nt_conv, t_conv, &
            conv_tol, grid_pts, grid_ele, grid_nod, grid_nnod, exch_nval, exch_nrow, exch_val, exch_rows, &
            exch_cols, grid_abc, usePrecision, nThreadsMatlab, N_ave, &
            CV, useReturnHall, useAvgN, demigstp, exch_weigh, exch_meth, exch_intpn, &
            n_macro, shiftVec, macroShape, sampleShape, exchPBC, &
            passExch, exch_ncols, CrysAxis, K0_arr, K1, K2, problem, dummy_run, fmm_cells_per_node, eps_fmm, ifunif, nlmin, nlmax, allow_fmm_short_circuit, fmm_min_n, fmm_nterms, use_fmm, &
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

        call SolveLandauLifshitzEquation( problem, solution )


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
        n_Hext_accepted = 0
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
