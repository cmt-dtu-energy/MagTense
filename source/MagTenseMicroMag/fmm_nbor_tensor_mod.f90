

module fmm_nbor_tensor_mod
    use MKL_SPBLAS
    use MKL_DFTI
    use MicroMagParameters
    use fmm3d_tree_mod
    use TileNComponents
    use DemagFieldGetSolution
    use IO_GENERAL
    use trace_mod
  implicit none
contains 

!=====================================================================
!> Compute dipole tensor Kdip for vector Rvec
!!  Kdip = (1/(4π r^3)) * (3 r̂ r̂^T - I), mapping M -> H (SI)
!! 
!!  arguments:
!!   Rvec  - (3) vector from source to target
!!   Kdip  - (3,3) dipole tensor output
!=====================================================================
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



!=====================================================================
!> Check if two cells are neighbours based on their positions and sizes
!!  NOTE - not currently used - kept for reference. Replaced by FMM tree-based neighbour list.
!!
!!  arguments:
!!   idx_t     - target cell index
!!   idx_s     - source cell index
!!   offset    - (ntot,3) array of cell positions
!!   size      - (ntot,3) array of cell sizes
!!   pitch     - (3) array of cell pitches
!!   rad_cells - radius [in cells] to consider neighbours. fx. 1 = immediate neighbours, 2 = next-nearest, etc.
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




!=====================================================================
!> Build neighbour list from FMM tree's list1 (near-field) boxes
!!
!!  arguments:
!!   problem  - MicroMagProblem object to fill n_nbors and nbr_idx
!!   tree     - FMM3DTree object from which to extract neighbour info
!=====================================================================
subroutine BuildNeighbourList_FromTree(problem, tree)
  implicit none
  type(MicroMagProblem), intent(inout) :: problem
  type(FMM3DTree),       intent(in)    :: tree

  integer :: ntot, nneigh_max
  integer :: ilev, ibox, i, jbox
  integer :: istarts, iends, jstart, jend
  integer :: p_t, p_s, t_idx, j_idx
  integer :: nb, cnt
  integer, save :: itimer = 0

  call trace%begin( "BuildNeighbourList_FromTree", itimer=itimer, verbose=1 )

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

  print *, " counting neighbours from FMM tree ..."
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


  print *, " filling neighbour list from FMM tree ..."
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


  call trace%end( "BuildNeighbourList_FromTree", itimer=itimer, verbose=1 )

end subroutine BuildNeighbourList_FromTree


!=====================================================================
!> Build dense neighbour list (all-to-all)
!! Not currently used - for testing purposes only.
!!
!!  arguments:
!!   problem  - MicroMagProblem object to fill n_nbors and nbr_idx
!=====================================================================
subroutine BuildNeighbourList_Dense(problem)
  implicit none
  type(MicroMagProblem), intent(inout) :: problem
  !------------------------------------------------------
  integer :: i,j, ntot
  !------------------------------------------------------

  ntot = size(problem%grid%pts, dim=1)

  if (associated(problem%n_nbors)) nullify(problem%n_nbors)
  allocate(problem%n_nbors(ntot))
  problem%n_nbors = ntot

  if (associated(problem%nbr_idx)) nullify(problem%nbr_idx)
  allocate(problem%nbr_idx(ntot, ntot))
  do i = 1, ntot
      do j = 1, ntot
          problem%nbr_idx(i,j) = j
      end do
  end do

end subroutine BuildNeighbourList_Dense

!=====================================================================
!> Build neighbour list from FMM tree's list1 (near-field) boxes
!!
!!  arguments:
!!   problem  - MicroMagProblem object to fill n_nbors and nbr_idx
!=====================================================================
subroutine BuildNeighbourDemagTensor(problem)
  !-----------------------------------------------------------------
  ! Construct exact prism demag tensors for near neighbours of
  ! every cell, for a (possibly) non-uniform grid.
  !
  ! This version:
  !  - Keeps the OpenMP directive right on the loop over t
  !  - Allocates ALL workspaces OUTSIDE the parallel loop
  !  - Calls getFieldFromTiles once per neighbour (m calls per target)
  !    using per-thread allocatable "whole arrays" to satisfy an
  !    allocatable dummy arg interface for Nout.
  !
  ! Neighbour-slot semantics are preserved: result goes into (t,q,:,:).
  !-----------------------------------------------------------------
  use omp_lib
  implicit none

  type(MicroMagProblem), intent(inout) :: problem

  ! Local aliases / pointers to problem fields
  integer,  dimension(:,:),       pointer :: nbr_idx_p
  real(SP), dimension(:,:,:,:),   pointer :: Nnbr_p

  ! Grid-related locals
  integer :: ntot
  real(DP), allocatable :: offset(:,:)      ! (ntot,3)
  real(DP), allocatable :: size_cell(:,:)   ! (ntot,3)
  real(DP) :: pitch(3)
  integer  :: nneigh_max
  integer  :: d

  ! Neighbour / tensor loop locals (declared here; user said they moved defs to top)
  integer :: t, m, q, s
  real(DP) :: pts_t(1,3)

  ! FMM-tree neighbour construction locals
  integer :: ier, i
  type(FMM3DTree), pointer :: fmm_tree
  real(DP), contiguous, pointer :: sources(:,:)

  !-------------------------
  ! Per-thread workspaces
  !-------------------------
  type :: TileWS
    type(MagTile), allocatable :: a(:)      ! (1)
  end type TileWS

  type :: NoutWS
    real(DP), allocatable :: a(:,:,:,:)     ! (1,1,3,3)
  end type NoutWS

  type :: HdummyWS
    real(DP), allocatable :: a(:,:)         ! (1,3)
  end type HdummyWS

  integer :: nthreads, tid
  type(TileWS),  allocatable :: tiles_thr(:)
  type(NoutWS),  allocatable :: nout_thr(:)
  type(HdummyWS),allocatable :: hdum_thr(:)

  integer, save :: itimer = 0
  call trace%begin( "BuildNeighbourDemagTensor", itimer=itimer, verbose=2 )


  ier = 0


  !---------------------------------------
  ! 1. Determine ntot from pts array
  !---------------------------------------
  ntot = size(problem%grid%pts, dim=1)

  allocate(offset(ntot,3))
  allocate(size_cell(ntot,3))
  offset(:,:)    = problem%grid%pts(:,:)
  size_cell(:,:) = problem%grid%abc(:,:)

  do d = 1, 3
    pitch(d) = sum(size_cell(:,d)) / real(ntot, DP)
  end do

  !---------------------------------------
  ! 2. Clean up any existing neighbour data
  !---------------------------------------
  if (associated(problem%nbr_idx)) deallocate(problem%nbr_idx)
  if (associated(problem%Nnbr))    deallocate(problem%Nnbr)

  !---------------------------------------
  ! 3. Build neighbour list from FMM tree lists
  !---------------------------------------
  allocate(fmm_tree)
  allocate(sources(3,ntot))

  do i = 1, ntot
    sources(1,i) = real(problem%grid%pts(i,1), DP)
    sources(2,i) = real(problem%grid%pts(i,2), DP)
    sources(3,i) = real(problem%grid%pts(i,3), DP)
  end do

  call fmm_tree%build_tree(sources, problem%fmm_eps, problem%fmm_cells_per_node, ier, problem%ifunif, problem%nlmin, problem%nlmax)
  if (ier /= 0) stop "BuildNeighbourDemagTensor: fmm_tree%build_tree failed"

  call BuildNeighbourList_FromTree(problem, fmm_tree)

  nneigh_max = maxval(problem%n_nbors)
  print *, " nbor list built from tree has ", sum(problem%n_nbors), " total nbor"
  print *, " max nbor = ", nneigh_max

  deallocate(fmm_tree)
  deallocate(sources)


  if (problem%allow_fmm_short_circuit .eq. 1 .and. ntot < problem%fmm_min_n) then
    print *, " Short circuiting FMM - disabling FMM"
    call dealloc_fmm_arrays(problem)
    problem%use_fmm = .false.
    return
  end if 




  nbr_idx_p => problem%nbr_idx

  ! If no neighbours at all, allocate empty Nnbr and bail
  if (nneigh_max <= 0) then
    allocate(problem%Nnbr(ntot, 0, 3, 3))
    problem%Nnbr = 0.0_SP
    deallocate(offset, size_cell)
    return
  end if

  !---------------------------------------
  ! 4. Allocate Nnbr pointer and attach to problem
  !---------------------------------------
  allocate(Nnbr_p(ntot, nneigh_max, 3, 3))
  Nnbr_p = 0.0_SP
  problem%Nnbr => Nnbr_p

  call displayGUIMessage(" build demag tensor from fmm nbors (per-neighbour calls)")

  !---------------------------------------
  ! 5. Allocate per-thread workspaces OUTSIDE the parallel loop
  !    Determine number of threads via a parallel singleton.
  !---------------------------------------
  nthreads = 1
  !$omp parallel default(none) shared(nthreads)
    !$omp single
      nthreads = omp_get_num_threads()
    !$omp end single
  !$omp end parallel

  allocate(tiles_thr(nthreads))
  allocate(nout_thr(nthreads))
  allocate(hdum_thr(nthreads))

  do tid = 1, nthreads
    allocate(tiles_thr(tid)%a(1))
    allocate(nout_thr(tid)%a(1,1,3,3))
    allocate(hdum_thr(tid)%a(1,3))
    nout_thr(tid)%a = 0.0_DP
    hdum_thr(tid)%a = 0.0_DP
  end do

  !---------------------------------------
  ! 6. Loop over target cells and build demag tensors
  !    (OpenMP directive stays right at the loop)
  !---------------------------------------
  !$omp parallel do default(none) &
  !$omp shared(ntot, problem, nbr_idx_p, Nnbr_p, offset, size_cell, tiles_thr, nout_thr, hdum_thr) &
  !$omp private(t, m, q, s, pts_t, tid) schedule(static)
  do t = 1, ntot

    tid = omp_get_thread_num() + 1

    m = problem%n_nbors(t)
    if (m <= 0) cycle

    pts_t(1,1) = offset(t,1)
    pts_t(1,2) = offset(t,2)
    pts_t(1,3) = offset(t,3)

    do q = 1, m
      s = nbr_idx_p(t,q)
      if (s < 1) cycle   ! skip invalid slots (-1/0). Adjust if your valid range differs.

      tiles_thr(tid)%a(1)%tileType        = 2
      tiles_thr(tid)%a(1)%a               = size_cell(s,1)
      tiles_thr(tid)%a(1)%b               = size_cell(s,2)
      tiles_thr(tid)%a(1)%c               = size_cell(s,3)
      tiles_thr(tid)%a(1)%exploitSymmetry = 0
      tiles_thr(tid)%a(1)%rotAngles(:)    = 0.0_DP
      tiles_thr(tid)%a(1)%M(:)            = 0.0_DP
      tiles_thr(tid)%a(1)%offset(1)       = offset(s,1)
      tiles_thr(tid)%a(1)%offset(2)       = offset(s,2)
      tiles_thr(tid)%a(1)%offset(3)       = offset(s,3)

      nout_thr(tid)%a = 0.0_DP

      ! Whole allocatable actual argument (OK with allocatable dummy)
      call getFieldFromTiles( tiles_thr(tid)%a, hdum_thr(tid)%a, pts_t, 1, 1, nout_thr(tid)%a, .false. )

      Nnbr_p(t, q, :, :) = sngl(nout_thr(tid)%a(1,1,:,:))
    end do

  end do
  !$omp end parallel do

  call displayGUIMessage(" done building demag tensor from fmm nbors")

  !---------------------------------------
  ! 7. Clean up temporaries
  !---------------------------------------
  do tid = 1, nthreads
    if (allocated(hdum_thr(tid)%a))  deallocate(hdum_thr(tid)%a)
    if (allocated(nout_thr(tid)%a))  deallocate(nout_thr(tid)%a)
    if (allocated(tiles_thr(tid)%a)) deallocate(tiles_thr(tid)%a)
  end do
  deallocate(hdum_thr, nout_thr, tiles_thr)

  deallocate(offset, size_cell)

  !---------------------------------------
  ! 8. Attach diffTens and build CSR sparse matrices (opt)
  !---------------------------------------
  problem%diffTens => problem%Nnbr
  call build_FMM_sparse_nborTensor_opt(problem)


  call trace%end( "BuildNeighbourDemagTensor", itimer=itimer, verbose=2 )
end subroutine BuildNeighbourDemagTensor




!=====================================================================
!> Convert Nnbr to diffTens by subtracting dipole tensor contribution
!!  Note - not currently used - kept for reference.
!! 
!!
!!  arguments:
!!   problem  - MicroMagProblem object containing Nnbr; diffTens will be allocated here
subroutine convert_Nnbr_to_diffTens(problem)
  implicit none
  type(MicroMagProblem),  intent(inout)    :: problem
  !---------------------------------------
  integer :: ntot, t, m, jidx
  real(DP) :: dx, dy, dz, volj
  real(DP) :: xt, yt, zt, xj, yj, zj
  real(DP) :: Rvec(3), Kdip(3,3), Kloc(3,3)
  real(SP), contiguous, pointer :: diffTens(:,:,:,:)
  integer, save :: itimer = 0
  !---------------------------------------
  call trace%begin( "convert_Nnbr_to_diffTens", itimer=itimer, verbose=2 )

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

      problem%diffTens(t, m, :,:) = Kloc - Kdip * volj  ! correction compared to dipole 
      !problem%diffTens(t, m, :,:) = Kloc  ! full tensor - if eval_direct is never called 
    end do
  end do

  call trace%end( "convert_Nnbr_to_diffTens", itimer=itimer, verbose=2 )
end subroutine convert_Nnbr_to_diffTens



!=====================================================================
!> Build sparse neighbour correction tensors in CSR format
!! NOTE - optimized version - faster than build_FMM_sparse_nborTensor, but more complex
!! 
!!  arguments:
!!   problem  - MicroMagProblem object containing nbr_idx, n_nbors, diffTens
subroutine build_fmm_sparse_nborTensor_opt(problem)
  implicit none
  type(MicroMagProblem), intent(inout) :: problem

  integer :: ntot, t, m, j, idx
  integer :: nnz_total
  integer :: pos, k
  integer :: stat

  integer, allocatable :: row_nnz(:)

  integer :: clk_start, clk_end, clk_rate
  real(SP) :: t_sec
  character(len=128) :: msg
  integer, save :: itimer = 0
  call trace%begin( "build_fmm_sparse_nborTensor_opt", itimer=itimer, verbose=2 )


  !-------- optimized not fully tested - use with caution.

  call system_clock(count_rate=clk_rate)

  ! Guards
  if (.not. associated(problem%nbr_idx))   stop "build_fmm_sparse_nborTensor_opt: nbr_idx not associated"
  if (.not. associated(problem%n_nbors))  stop "build_fmm_sparse_nborTensor_opt: n_nbors not associated"
  if (.not. associated(problem%diffTens)) stop "build_fmm_sparse_nborTensor_opt: diffTens not associated"

  call displayGUIMessage(" Starting FMM sparse nbor build (opt)")

  ntot = size(problem%grid%pts, dim=1)

  ! Destroy if already built
  if (problem%K_fmm_built) call destroy_nbrcorr_sparse(problem)

  allocate(row_nnz(ntot))
  row_nnz = 0

  !--------------------------------------------------------------------------
  ! Pass 1: Count valid neighbours per row
  !--------------------------------------------------------------------------
  call displayGUIMessage(" Pass 1/4: Counting valid neighbours per row")
  call system_clock(clk_start)

  !$omp parallel do default(shared) private(t,m,j) schedule(static)
  do t = 1, ntot
    do m = 1, problem%n_nbors(t)
      j = problem%nbr_idx(t,m)
      if (j >= 0) row_nnz(t) = row_nnz(t) + 1
    end do
  end do
  !$omp end parallel do

  call system_clock(clk_end)
  t_sec = real(clk_end - clk_start, SP) / real(clk_rate, SP)
  write(msg,'(A,F8.3,A)') "  -> Pass 1 time: ", t_sec, " s"
  call displayGUIMessage(msg)

  nnz_total = sum(row_nnz)

  !--------------------------------------------------------------------------
  ! Pass 2: Allocate CSR storage
  !--------------------------------------------------------------------------
  call displayGUIMessage(" Pass 2/4: Allocating CSR storage")
  call system_clock(clk_start)

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
  end do

  call system_clock(clk_end)
  t_sec = real(clk_end - clk_start, SP) / real(clk_rate, SP)
  write(msg,'(A,F8.3,A)') "  -> Pass 2 time: ", t_sec, " s"
  call displayGUIMessage(msg)

  !--------------------------------------------------------------------------
  ! Pass 3: Build CSR row pointers (prefix sum)
  !--------------------------------------------------------------------------
  call displayGUIMessage(" Pass 3/4: Building CSR row pointers")
  call system_clock(clk_start)

  pos = 1
  do t = 1, ntot
    do idx = 1, 6
      problem%K_fmm_s(idx)%rows_start(t) = pos
      problem%K_fmm_s(idx)%rows_end(t)   = pos + row_nnz(t)
    end do
    pos = pos + row_nnz(t)
  end do

  if (pos /= nnz_total + 1) stop "build_fmm_sparse_nborTensor_opt: nnz mismatch"

  call system_clock(clk_end)
  t_sec = real(clk_end - clk_start, SP) / real(clk_rate, SP)
  write(msg,'(A,F8.3,A)') "  -> Pass 3 time: ", t_sec, " s"
  call displayGUIMessage(msg)

  !--------------------------------------------------------------------------
  ! Pass 4: Fill cols + values
  !--------------------------------------------------------------------------
  call displayGUIMessage(" Pass 4/4: Filling CSR cols/values")
  call system_clock(clk_start)

  !$omp parallel do default(shared) private(t,m,j,pos,k) schedule(static)
  do t = 1, ntot
    pos = problem%K_fmm_s(1)%rows_start(t)
    k   = 0

    do m = 1, problem%n_nbors(t)
      j = problem%nbr_idx(t,m)
      if (j < 0) cycle
      k = k + 1

      problem%K_fmm_s(1)%cols(pos+k-1) = j
      problem%K_fmm_s(2)%cols(pos+k-1) = j
      problem%K_fmm_s(3)%cols(pos+k-1) = j
      problem%K_fmm_s(4)%cols(pos+k-1) = j
      problem%K_fmm_s(5)%cols(pos+k-1) = j
      problem%K_fmm_s(6)%cols(pos+k-1) = j

      problem%K_fmm_s(1)%values(pos+k-1) = problem%diffTens(t,m,1,1)
      problem%K_fmm_s(2)%values(pos+k-1) = problem%diffTens(t,m,1,2)
      problem%K_fmm_s(3)%values(pos+k-1) = problem%diffTens(t,m,1,3)
      problem%K_fmm_s(4)%values(pos+k-1) = problem%diffTens(t,m,2,2)
      problem%K_fmm_s(5)%values(pos+k-1) = problem%diffTens(t,m,2,3)
      problem%K_fmm_s(6)%values(pos+k-1) = problem%diffTens(t,m,3,3)
    end do

    if (k /= row_nnz(t)) stop "build_fmm_sparse_nborTensor_opt: row fill mismatch"
  end do
  !$omp end parallel do

  call system_clock(clk_end)
  t_sec = real(clk_end - clk_start, SP) / real(clk_rate, SP)
  write(msg,'(A,F8.3,A)') "  -> Pass 4 time: ", t_sec, " s"
  call displayGUIMessage(msg)

  deallocate(row_nnz)

  !--------------------------------------------------------------------------
  ! Create MKL handles
  !--------------------------------------------------------------------------
  call displayGUIMessage(" Creating MKL CSR handles")
  call system_clock(clk_start)

  do idx = 1, 6
    stat = mkl_sparse_s_create_csr( problem%K_fmm_s(idx)%A, SPARSE_INDEX_BASE_ONE, &
                                   ntot, ntot, &
                                   problem%K_fmm_s(idx)%rows_start, problem%K_fmm_s(idx)%rows_end, &
                                   problem%K_fmm_s(idx)%cols, problem%K_fmm_s(idx)%values )
    if (stat /= SPARSE_STATUS_SUCCESS) stop "build_fmm_sparse_nborTensor_opt: mkl_sparse_s_create_csr failed"
  end do

  call system_clock(clk_end)
  t_sec = real(clk_end - clk_start, SP) / real(clk_rate, SP)
  write(msg,'(A,F8.3,A)') "  -> MKL handle creation time: ", t_sec, " s"
  call displayGUIMessage(msg)

  problem%K_fmm_descr_s%type = SPARSE_MATRIX_TYPE_GENERAL
  problem%K_fmm_built = .true.

  call displayGUIMessage(" Finished FMM sparse nbor build (opt)")


  call trace%end( "build_fmm_sparse_nborTensor_opt", itimer=itimer, verbose=2 )

end subroutine build_fmm_sparse_nborTensor_opt



!=====================================================================
!> Build sparse neighbour correction tensors in CSR format
!! NOTE - semi-deprecated - use build_fmm_sparse_nborTensor_opt instead
!!        code is kept as it is simpler to read/understand - but is slower. 
!!
!!  arguments:
!!   problem  - MicroMagProblem object containing nbr_idx, n_nbors, diffTens
!=====================================================================
subroutine build_FMM_sparse_nborTensor(problem)
  implicit none
  type(MicroMagProblem), intent(inout) :: problem

  integer :: ntot, t, m, j, row
  integer :: nb_valid, nnz_total
  integer :: pos
  integer :: stat, idx

  real(SP) :: v_xx, v_xy, v_xz, v_yy, v_yz, v_zz
  integer, save :: itimer = 0
  call trace%begin( "build_FMM_sparse_nborTensor", itimer=itimer, verbose=2 )

  !TODO - this could be optimized by parralizing

  call displayGUIMessage( " Starting FMM sparse nbor built" )


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


    call displayGUIMessage( " row pointers built" )

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

    call displayGUIMessage( " K_fmm_s fully built" )


  ! Create MKL handles
  do idx = 1, 6
    stat = mkl_sparse_s_create_csr( problem%K_fmm_s(idx)%A, SPARSE_INDEX_BASE_ONE, &
                                   ntot, ntot, &
                                   problem%K_fmm_s(idx)%rows_start, problem%K_fmm_s(idx)%rows_end, &
                                   problem%K_fmm_s(idx)%cols, problem%K_fmm_s(idx)%values )
    if (stat /= SPARSE_STATUS_SUCCESS) stop "build_FMM_sparse_nborTensor: mkl_sparse_s_create_csr failed"

    !stat = mkl_sparse_optimize(problem%K_fmm_s(idx)%A)
    !if (stat /= SPARSE_STATUS_SUCCESS) stop "build_FMM_sparse_nborTensor: mkl_sparse_optimize failed"
  end do

      call displayGUIMessage( " mkl_sparse matrix of FMM nbor built" )


  ! Descriptor (GENERAL) – store once in problem for reuse
  problem%K_fmm_descr_s%type = SPARSE_MATRIX_TYPE_GENERAL

  problem%K_fmm_built = .true.

  call trace%end( "build_FMM_sparse_nborTensor", itimer=itimer, verbose=2 )

end subroutine build_FMM_sparse_nborTensor


!=====================================================================
!> Destroy FMM neighbour correction sparse matrices
!!
!!  arguments:
!!   problem  - MicroMagProblem object to deallocate FMM sparse matrices from
!=====================================================================
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


!=====================================================================
!> Deallocate FMM neighbour tensor arrays
!!
!!  arguments:
!!   problem  - MicroMagProblem object to deallocate FMM arrays from
!=====================================================================
subroutine dealloc_fmm_arrays(problem)
  implicit none
  type(MicroMagProblem), intent(inout) :: problem

  ! Deallocate neighbour tensor arrays
  if (associated(problem%Nnbr))    deallocate(problem%Nnbr)
  !if (associated(problem%diffTens)) deallocate(problem%diffTens) !currently just points to Nnbr
  if (associated(problem%nbr_idx))  deallocate(problem%nbr_idx)
  if (associated(problem%n_nbors))  deallocate(problem%n_nbors)

  ! Deallocate sparse FMM neighbour correction matrices
  if (problem%K_fmm_built) then
    call destroy_nbrcorr_sparse(problem)
  end if

end subroutine dealloc_fmm_arrays
end module fmm_nbor_tensor_mod
