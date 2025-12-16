module fmm_nbor_tensor_mod
    use MKL_SPBLAS
    use MKL_DFTI


    use MicroMagParameters
    use fmm3d_tree_mod
    use TileNComponents
    use DemagFieldGetSolution
  implicit none
contains 

        !------------------- for getting near-field correction ------------------------------------------
    pure subroutine dipole_tensor_3x3(Rvec, Kdip)
        implicit none
        ! Kdip = (1/(4π r^3)) * (3 r̂ r̂^T - I), mapping M -> H (SI)
        real(8), intent(in)  :: Rvec(3)
        real(8), intent(out) :: Kdip(3,3)
        real(8), parameter   :: fourpi = 12.56637061435917295385d0
        real(8), parameter   :: epsR   = 1.0d-15
        real(8) :: r2, rmag, invR3, ux, uy, uz

        r2  = Rvec(1)*Rvec(1) + Rvec(2)*Rvec(2) + Rvec(3)*Rvec(3)
        rmag = sqrt(r2)

        if (rmag <= epsR) then
            Kdip = 0.0d0
            return
        end if

        invR3 = 1.0d0 / (fourpi * r2 * rmag)
        ux = Rvec(1)/rmag
        uy = Rvec(2)/rmag
        uz = Rvec(3)/rmag

        Kdip(1,1) = (3.0d0*ux*ux - 1.0d0) * invR3
        Kdip(1,2) = (3.0d0*ux*uy)         * invR3
        Kdip(1,3) = (3.0d0*ux*uz)         * invR3
        Kdip(2,1) = Kdip(1,2)
        Kdip(2,2) = (3.0d0*uy*uy - 1.0d0) * invR3
        Kdip(2,3) = (3.0d0*uy*uz)         * invR3
        Kdip(3,1) = Kdip(1,3)
        Kdip(3,2) = Kdip(2,3)
        Kdip(3,3) = (3.0d0*uz*uz - 1.0d0) * invR3
    end subroutine dipole_tensor_3x3

    pure logical function is_neighbour(idx_t, idx_s, offset, size, pitch, rad_cells) result(mask)
        implicit none
        integer(4), intent(in) :: idx_t, idx_s, rad_cells
        real(8),   intent(in) :: offset(:,:), size(:,:), pitch(3)
        real(8) :: dx, dy, dz, tx, ty, tz

        if (idx_t == idx_s) then
            mask = .false.
            return
        end if

        dx = abs(offset(idx_s,1) - offset(idx_t,1))
        dy = abs(offset(idx_s,2) - offset(idx_t,2))
        dz = abs(offset(idx_s,3) - offset(idx_t,3))

        tx = rad_cells * pitch(1) + 0.5d0*(size(idx_s,1) + size(idx_t,1))
        ty = rad_cells * pitch(2) + 0.5d0*(size(idx_s,2) + size(idx_t,2))
        tz = rad_cells * pitch(3) + 0.5d0*(size(idx_s,3) + size(idx_t,3))

        mask = (dx <= tx) .and. (dy <= ty) .and. (dz <= tz)
    end function is_neighbour

subroutine BuildNeighbourList_FromTree(problem, tree)
  implicit none
  type(MicroMagProblem), intent(inout) :: problem
  type(FMM3DTree),       intent(in)    :: tree

  integer :: ntot, nneigh_max
  integer :: ilev, ibox, i, jbox
  integer :: istarts, iends, jstart, jend
  integer :: p_t, p_s, t_idx, j_idx
  integer :: nb, cnt

  ! --- basic sanity ---
  if (.not. tree%is_built) stop "BuildNeighbourList_FromTree: tree not built"
  if (.not. associated(tree%laddr))  stop "BuildNeighbourList_FromTree: tree%laddr not associated"
  if (.not. associated(tree%isrcse)) stop "BuildNeighbourList_FromTree: tree%isrcse not associated"
  if (.not. associated(tree%isrc))   stop "BuildNeighbourList_FromTree: tree%isrc not associated"
  if (.not. associated(tree%nlist1)) stop "BuildNeighbourList_FromTree: tree%nlist1 not associated"
  if (.not. associated(tree%list1))  stop "BuildNeighbourList_FromTree: tree%list1 not associated"

  ntot = tree%nsource   ! should match size(problem%grid%pts,1) in your usage

  ! Allocate / reset n_nbors
  if (associated(problem%n_nbors)) nullify(problem%n_nbors)
  allocate(problem%n_nbors(ntot))
  problem%n_nbors = 0

  !===========================
  ! Pass 1: count neighbours
  !===========================
  do ilev = 0, tree%nlevels
    do ibox = tree%laddr(1,ilev), tree%laddr(2,ilev)

      istarts = tree%isrcse(1,ibox)
      iends   = tree%isrcse(2,ibox)
      if (iends < istarts) cycle

      ! number of sources in all list1 boxes (same for all targets in this ibox),
      ! minus 1 for self if ibox is included in its own list1.
      nb = 0
      do i = 1, tree%nlist1(ibox)
        jbox   = tree%list1(i,ibox)
        jstart = tree%isrcse(1,jbox)
        jend   = tree%isrcse(2,jbox)
        if (jend >= jstart) nb = nb + (jend - jstart + 1)
        if (jbox == ibox)   nb = nb - 1
      end do
      if (nb < 0) nb = 0

      do p_t = istarts, iends
        t_idx = tree%isrc(p_t)          ! original index for this target
        problem%n_nbors(t_idx) = nb
      end do

    end do
  end do

  nneigh_max = maxval(problem%n_nbors) 
  nneigh_max = nneigh_max + 1 ! To include self 

  ! Allocate / reset nbr_idx
  if (associated(problem%nbr_idx)) nullify(problem%nbr_idx)
  allocate(problem%nbr_idx(ntot, nneigh_max))
  !problem%nbr_idx = -1

  !===========================
  ! Pass 2: fill neighbour ids
  !===========================
  do ilev = 0, tree%nlevels
    do ibox = tree%laddr(1,ilev), tree%laddr(2,ilev)

      istarts = tree%isrcse(1,ibox)
      iends   = tree%isrcse(2,ibox)
      if (iends < istarts) cycle

      do p_t = istarts, iends
        t_idx = tree%isrc(p_t)
        cnt   = 0

        do i = 1, tree%nlist1(ibox)
          jbox   = tree%list1(i,ibox)
          jstart = tree%isrcse(1,jbox)
          jend   = tree%isrcse(2,jbox)
          if (jend < jstart) cycle

          do p_s = jstart, jend
            j_idx = tree%isrc(p_s)
            !if (j_idx == t_idx) cycle   ! exclude self
            cnt = cnt + 1
            if (cnt <= nneigh_max) then
              problem%nbr_idx(t_idx, cnt) = j_idx
            else
              stop "BuildNeighbourList_FromTree: cnt exceeded nneigh_max (unexpected)"
            end if
          end do
        end do

        ! overwrite with exact count (robust, avoids relying on the "-1 if jbox==ibox" assumption)
        problem%n_nbors(t_idx) = cnt

      end do 
    end do
  end do

end subroutine BuildNeighbourList_FromTree



subroutine BuildNeighbourDemagTensor(problem)
  !-----------------------------------------------------------------
  ! Construct exact prism demag tensors for near neighbours of
  ! every cell, for a (possibly) non-uniform grid.
  !
  ! Neighbours are defined geometrically via is_neighbour(), which
  ! uses:
  !   - cell centres: offset(ℓ,1:3)  = problem%grid%pts(ℓ,1:3)
  !   - cell sizes:   size_cell(ℓ,1:3) = problem%grid%abc(ℓ,1:3)
  !                   (for uniform grids, we fall back to dx,dy,dz)
  !
  ! Inputs:
  !   problem        : MicroMagProblem (owns grid & output storage)
  !
  ! Outputs (attached to problem):
  !   problem%nbr_idx(ntot, nneigh_max)
  !   problem%Nnbr   (ntot, nneigh_max, 3, 3)
  !
  ! Neighbour relation is symmetric by construction.
  !-----------------------------------------------------------------
  implicit none

  type(MicroMagProblem), intent(inout) :: problem

  ! Local aliases / pointers to problem fields
  integer,  dimension(:,:),   pointer :: nbr_idx_p
  real(SP), dimension(:,:,:,:), pointer :: Nnbr_p

  ! Grid-related locals
  integer :: ntot
  real(DP), contiguous, pointer :: offset(:,:)      ! (ntot,3)
  real(DP), contiguous, pointer :: size_cell(:,:)   ! (ntot,3)
  real(DP) :: pitch(3)
  integer  :: nneigh_max

  ! Temp neighbour list
  integer, contiguous, pointer  :: nbr_idx_loc(:,:)
  integer, contiguous, pointer  :: n_nbors(:)

  ! Demag construction locals
  integer :: t, m, q, s
  type(MagTile), allocatable :: tiles(:)
  real(DP), dimension(1,3) :: pts_t
  real(DP), allocatable :: H_dummy(:,:)
  real(DP), allocatable :: Nout(:,:,:,:)  ! (m, 1, 3, 3)
  integer :: d


  !-----------
  integer :: ier, i 
  type(FMM3DTree), pointer :: fmm_tree
  real(DP), contiguous, pointer :: sources(:,:)
  real(DP) :: eps_fmm

  !---------------------------------------
  ! 1. Determine ntot from pts array
  !---------------------------------------
  ntot = size(problem%grid%pts, dim=1)



  !offset(:,:) => problem%grid%pts(1:ntot, 1:3)
  !size_cell(:,:) => problem%grid%abc(1:ntot, 1:3)


  allocate(offset(ntot,3))
  allocate(size_cell(ntot,3))
  offset(:,:) = problem%grid%pts(:,:)
  size_cell(:,:) = problem%grid%abc(:,:)

  ! Nominal pitch: for uniform, dx,dy,dz; for non-uniform, use mean cell size
  do d = 1, 3
    pitch(d) = sum(size_cell(:,d)) / real(ntot, DP)
  end do

  !---------------------------------------
  ! 2. Clean up any existing neighbour data
  !---------------------------------------
  if (associated(problem%nbr_idx)) deallocate(problem%nbr_idx)
  if (associated(problem%Nnbr))    deallocate(problem%Nnbr)

  !----------- build neighbour list --------------
  !call BuildNeighbourList(offset, size_cell, pitch, radius_cells, &
  !                        nbr_idx_loc, nneigh_max)



  !------------------- built neighbour list based on list1 in tree --------------
  !        first builts a temporary tree to get the list
  !        then stores all points in list 1 for each tile as a neighbour list
  allocate(fmm_tree)
  allocate(sources(3,ntot))
  do i = 1, ntot
      !------------------ positions -> FMM (3, N) ------------------
      sources(1,i) = real(problem%grid%pts(i,1), DP)
      sources(2,i) = real(problem%grid%pts(i,2), DP)
      sources(3,i) = real(problem%grid%pts(i,3), DP)
      !------------------------------------------------------------
  end do
  !eps_fmm = 1.0e-4_DP
  call fmm_tree%build_tree( sources, problem%fmm_eps, problem%fmm_cells_per_node , ier)
  call BuildNeighbourList_FromTree(problem, fmm_tree)
  nneigh_max = maxval(problem%n_nbors) 
  print *, " nbor list built  from tree has ", sum(problem%n_nbors), " total nbor"
  print *, " max nbor = ", maxval(problem%n_nbors)
  deallocate(fmm_tree)
  deallocate(sources)
  !---------------------------------------------------------------------------------


  !------------- built neighbour list based on geometry ----------------------------
  !            first finds all "touching" tiles
  !            then iteratively adds all nbors nbors 'radius_cell'-times to get full nbor list
!  call BuildNeighbourList(offset, size_cell, radius_cells, &
!                           nbr_idx_loc, nneigh_max, n_nbors)
!   print *, " nbor list built from geometric has ", sum(n_nbors), " total nbor"
!   print *, " max nbor = ", nneigh_max
  ! problem%n_nbors => n_nbors
  ! problem%nbr_idx => nbr_idx_loc
  !--------------------------------------------------------------------------------
  nbr_idx_p => problem%nbr_idx

  ! If no neighbours at all, allocate empty Nnbr and bail
  if (nneigh_max <= 0) then
    allocate(problem%Nnbr(ntot, 0, 3, 3))
    problem%Nnbr = 0.0_SP
    deallocate(offset, size_cell)!, nbr_idx_loc)
    return
  end if


  !print *, " allocate neighbour demag tensors"
  !---------------------------------------
  ! 4. Allocate Nnbr pointer and attach to problem
  !---------------------------------------
  allocate(Nnbr_p(ntot, nneigh_max, 3, 3))
  Nnbr_p = 0.0_SP
  problem%Nnbr => Nnbr_p

  ! H_dummy is just a 1x3 placeholder for getFieldFromTiles
  allocate(H_dummy(1,3))
  H_dummy = 0.0_DP

  !print *, " looping over target cells to build demag tensors"
  !---------------------------------------
  ! 5. Loop over target cells and build demag tensors
  !---------------------------------------
 !$omp parallel do default(none) &
 !$omp shared(ntot, problem, nbr_idx_p, Nnbr_p, offset, size_cell, H_dummy) &
 !$omp private(m, t, q, s, tiles, pts_t, Nout) schedule(static)
  do t = 1, ntot

    ! Number of neighbours for this target
    m = problem%n_nbors(t)
    !m = count(nbr_idx_p(t,:) > 0)   ! since we use -1 for "no neighbour"

    if (m <= 0) cycle

    ! Build MagTile array for these m neighbours
    allocate(tiles(m))
    allocate(Nout(m,1,3,3))
    Nout = 0.0_DP

    do q = 1, m
      s = nbr_idx_p(t, q)

      tiles(q)%tileType        = 2           ! prism
      tiles(q)%a               = size_cell(s,1)
      tiles(q)%b               = size_cell(s,2)
      tiles(q)%c               = size_cell(s,3)
      tiles(q)%exploitSymmetry = 0
      tiles(q)%rotAngles(:)    = 0.0_DP
      tiles(q)%M(:)            = 0.0_DP

      tiles(q)%offset(1)       = offset(s,1)
      tiles(q)%offset(2)       = offset(s,2)
      tiles(q)%offset(3)       = offset(s,3)
    end do

    ! Target point is the centre of the target tile
    pts_t(1,1) = offset(t,1)
    pts_t(1,2) = offset(t,2)
    pts_t(1,3) = offset(t,3)

    ! Compute demag tensors from all neighbour tiles at this target
    call getFieldFromTiles( tiles, H_dummy, pts_t, m, 1, Nout, .false. )

    ! Store tensors in (t, neighbour-slot, 3, 3)
    Nnbr_p(t, 1:m, :, :) = real(Nout(:,1,:,:), SP)

    deallocate(Nout)
    deallocate(tiles)

  end do
  !$omp end parallel do

  !print *, " cleanup temporaries"
  !---------------------------------------
  ! 6. Clean up temporaries
  !---------------------------------------
  deallocate(H_dummy)
  deallocate(offset, size_cell)
  !deallocate(nbr_idx_loc)
  !print *, "done building nonuniform demag tensor"


  !------------ convert neighbour tensor to diff between exact and dipole -------------
  !call convert_Nnbr_to_diffTens(problem) ! no conversion needed currently as
  ! conversion needed if we allow FMM to first calculate dipole and then add correction
  ! but currently we dispable FMM direct evaluation and use only neighbour tensor
  problem%diffTens => problem%Nnbr
  !------------------------------------------------------------------------------------



  call build_FMM_sparse_nborTensor(problem)

end subroutine BuildNeighbourDemagTensor

subroutine convert_Nnbr_to_diffTens(problem)
  implicit none
  type(MicroMagProblem),  intent(inout)    :: problem
  !---------------------------------------
  integer :: ntot, t, m, jidx
  real(DP) :: dx, dy, dz, volj
  real(DP) :: xt, yt, zt, xj, yj, zj
  real(DP) :: Rvec(3), Kdip(3,3), Kloc(3,3)
  real(SP), contiguous, pointer :: diffTens(:,:,:,:)
  !---------------------------------------

  ntot = size(problem%grid%pts,1)
  dx = real(problem%grid%dx, DP)
  dy = real(problem%grid%dy, DP)
  dz = real(problem%grid%dz, DP)


  allocate(diffTens(size(problem%Nnbr,1), size(problem%Nnbr,2), size(problem%Nnbr,3), size(problem%Nnbr,4)))
  diffTens(:,:,:,:) = 0.0_SP
  problem%diffTens => diffTens

  do t = 1, ntot 
    xt = real(problem%grid%pts(t,1),DP)
    yt = real(problem%grid%pts(t,2),DP)
    zt = real(problem%grid%pts(t,3),DP)
    do m=1, problem%n_nbors(t)
      jidx = problem%nbr_idx(t,m)
      if (jidx < 0) cycle   ! empty slot (boundary)
      xj = real(problem%grid%pts(jidx,1),DP)
      yj = real(problem%grid%pts(jidx,2),DP)
      zj = real(problem%grid%pts(jidx,3),DP)


      Kloc(:,:) = real(problem%Nnbr(t,m,:,:), DP)
      Rvec(1) = xt - xj
      Rvec(2) = yt - yj
      Rvec(3) = zt - zj
      call dipole_tensor_3x3(Rvec, Kdip)
      if (problem%grid%gridType .eq. gridTypeUniform) then
          volj = dx * dy * dz
      else
          volj = problem%grid%abc(jidx,1) * problem%grid%abc(jidx,2) * problem%grid%abc(jidx,3)
      endif

      !problem%diffTens(t, m, :,:) = Kloc - Kdip * volj  ! correction compared to dipole 
      problem%diffTens(t, m, :,:) = Kloc  ! full tensor - if eval_direct is never called 
    end do
  end do
end subroutine convert_Nnbr_to_diffTens


subroutine build_FMM_sparse_nborTensor(problem)
  implicit none
  type(MicroMagProblem), intent(inout) :: problem

  integer :: ntot, t, m, j, row
  integer :: nb_valid, nnz_total
  integer :: pos
  integer :: stat, idx

  real(SP) :: v_xx, v_xy, v_xz, v_yy, v_yz, v_zz

  ! Guards
  if (.not. associated(problem%nbr_idx))   stop "build_FMM_sparse_nborTensor: nbr_idx not associated"
  if (.not. associated(problem%n_nbors))  stop "build_FMM_sparse_nborTensor: n_nbors not associated"
  if (.not. associated(problem%diffTens)) stop "build_FMM_sparse_nborTensor: diffTens not associated"

  ntot = size(problem%grid%pts, dim=1)

  ! Destroy if already built
  if (problem%K_fmm_built) call destroy_nbrcorr_sparse(problem)

  ! Count nnz: one scalar entry per valid neighbour link
  nnz_total = 0
  do t = 1, ntot
    do m = 1, problem%n_nbors(t)
      j = problem%nbr_idx(t,m)
      if (j < 0) cycle
      nnz_total = nnz_total + 1
    end do
  end do

  ! Allocate all 6 matrices
  do idx = 1, 6
    problem%K_fmm_s(idx)%nrows   = ntot
    problem%K_fmm_s(idx)%ncols   = ntot
    problem%K_fmm_s(idx)%nvalues = nnz_total

    allocate(problem%K_fmm_s(idx)%rows_start(ntot))
    allocate(problem%K_fmm_s(idx)%rows_end(ntot))
    allocate(problem%K_fmm_s(idx)%cols(nnz_total))
    allocate(problem%K_fmm_s(idx)%values(nnz_total))

    problem%K_fmm_s(idx)%rows_start = 0
    problem%K_fmm_s(idx)%rows_end   = 0
    problem%K_fmm_s(idx)%cols       = 0
    problem%K_fmm_s(idx)%values     = 0.0_SP
    !problem%K_fmm_s(idx)%A          = SPARSE_MATRIX_T_NULL
  end do

  ! Build row pointers (rows_start/rows_end) – identical for all 6
  pos = 1
  do t = 1, ntot
    nb_valid = 0
    do m = 1, problem%n_nbors(t)
      if (problem%nbr_idx(t,m) >= 0) nb_valid = nb_valid + 1
    end do

    ! row t has nb_valid entries
    do idx = 1, 6
      problem%K_fmm_s(idx)%rows_start(t) = pos
      problem%K_fmm_s(idx)%rows_end(t)   = pos + nb_valid   ! one-past-end
    end do

    pos = pos + nb_valid
  end do

  if (pos /= nnz_total + 1) stop "build_FMM_sparse_nborTensor: nnz mismatch"

  ! Fill cols + values
  do t = 1, ntot
    pos = problem%K_fmm_s(1)%rows_start(t)

    do m = 1, problem%n_nbors(t)
      j = problem%nbr_idx(t,m)
      if (j < 0) cycle

      ! Column index is source cell index
      do idx = 1, 6
        problem%K_fmm_s(idx)%cols(pos) = j
      end do

      ! Extract the 6 unique components from diffTens(t,m,:,:)
      v_xx = problem%diffTens(t,m,1,1)
      v_xy = problem%diffTens(t,m,1,2)
      v_xz = problem%diffTens(t,m,1,3)
      v_yy = problem%diffTens(t,m,2,2)
      v_yz = problem%diffTens(t,m,2,3)
      v_zz = problem%diffTens(t,m,3,3)

      problem%K_fmm_s(1)%values(pos) = v_xx
      problem%K_fmm_s(2)%values(pos) = v_xy
      problem%K_fmm_s(3)%values(pos) = v_xz
      problem%K_fmm_s(4)%values(pos) = v_yy
      problem%K_fmm_s(5)%values(pos) = v_yz
      problem%K_fmm_s(6)%values(pos) = v_zz

      pos = pos + 1
    end do

    if (pos /= problem%K_fmm_s(1)%rows_end(t)) stop "build_FMM_sparse_nborTensor: row fill mismatch"
  end do

  ! Create MKL handles
  do idx = 1, 6
    stat = mkl_sparse_s_create_csr( problem%K_fmm_s(idx)%A, SPARSE_INDEX_BASE_ONE, &
                                   ntot, ntot, &
                                   problem%K_fmm_s(idx)%rows_start, problem%K_fmm_s(idx)%rows_end, &
                                   problem%K_fmm_s(idx)%cols, problem%K_fmm_s(idx)%values )
    if (stat /= SPARSE_STATUS_SUCCESS) stop "build_FMM_sparse_nborTensor: mkl_sparse_s_create_csr failed"

    stat = mkl_sparse_optimize(problem%K_fmm_s(idx)%A)
    if (stat /= SPARSE_STATUS_SUCCESS) stop "build_FMM_sparse_nborTensor: mkl_sparse_optimize failed"
  end do

  ! Descriptor (GENERAL) – store once in problem for reuse
  problem%K_fmm_descr_s%type = SPARSE_MATRIX_TYPE_GENERAL

  problem%K_fmm_built = .true.

end subroutine build_FMM_sparse_nborTensor


subroutine destroy_nbrcorr_sparse(problem)
  implicit none
  type(MicroMagProblem), intent(inout) :: problem

  integer :: idx, stat

  do idx = 1, 6
    ! if (problem%K_fmm_s(idx)%A .ne. SPARSE_MATRIX_T_NULL) then
    !   stat = mkl_sparse_destroy(problem%K_fmm_s(idx)%A)
    !   problem%K_fmm_s(idx)%A = SPARSE_MATRIX_T_NULL
    ! end if

    if (allocated(problem%K_fmm_s(idx)%rows_start)) deallocate(problem%K_fmm_s(idx)%rows_start)
    if (allocated(problem%K_fmm_s(idx)%rows_end))   deallocate(problem%K_fmm_s(idx)%rows_end)
    if (allocated(problem%K_fmm_s(idx)%cols))       deallocate(problem%K_fmm_s(idx)%cols)
    if (allocated(problem%K_fmm_s(idx)%values))     deallocate(problem%K_fmm_s(idx)%values)

    problem%K_fmm_s(idx)%nrows   = 0
    problem%K_fmm_s(idx)%ncols   = 0
    problem%K_fmm_s(idx)%nvalues = 0
  end do

  problem%K_fmm_built = .false.

end subroutine destroy_nbrcorr_sparse

end module fmm_nbor_tensor_mod