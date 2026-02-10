!================================================
! timer_mod.f90  (OpenMP-capable, no MPI)
!================================================
MODULE timer_mod
  USE omp_timer_mod, only : wallclock, omp_timer
  USE omp_lock_mod,  only : omp_lock_t
  USE omp_mod,       only : omp
  implicit none
  private

  integer, parameter :: NAME_LEN  = 128
  integer, parameter :: MAX_DEPTH = 64

  !---------------------- Per-thread timer state ------------------------------
  type :: timer_thread_t
    integer(8), allocatable :: calls(:)            ! number of calls per timer id
    real(8),    allocatable :: total(:)            ! accumulated exclusive time per timer id

    integer :: depth = 0                           ! nesting depth (per thread)
    integer :: stack_id(MAX_DEPTH) = 0             ! active timer ids (per thread)
    real(8) :: stack_t0(MAX_DEPTH) = 0.0d0         ! start time for each stack level
    real(8) :: child_acc(MAX_DEPTH) = 0.0d0        ! accumulated child time for level
  end type timer_thread_t
  !----------------------------------------------------------------------------

  !---------------------- Global registry (shared) -----------------------------
  integer, save :: ntimer = 0                               ! number of registered timers
  character(len=NAME_LEN), allocatable, save :: names(:)    ! timer names (1..capacity)
  !-----------------------------------------------------------------------------

  !------------------------ Registry lock (protects registration) --------------
  type(omp_lock_t), save :: reg_lock                         ! lock for registry updates
  logical, save :: initialized = .false.                     ! has init() run?
  !-----------------------------------------------------------------------------

  !-------------------- Lock protecting timing log output ----------------------
  type(omp_lock_t), save :: log_lock
  !-----------------------------------------------------------------------------

  !-------------------- Window snapshot (for delta printing) -------------------
  integer(8), allocatable, save :: win_prev_calls(:)         ! previous snapshot calls
  real(8),    allocatable, save :: win_prev_total(:)         ! previous snapshot totals
  logical, save :: win_snap_initialized = .false.             ! snapshot initialised?
  !-----------------------------------------------------------------------------

  !-------------------- Per-thread timer storage -------------------------------
  type(timer_thread_t), allocatable, save :: timers(:)       ! timers(tid)
  !-----------------------------------------------------------------------------

  !----------------------------- Public timer object ---------------------------
  type, public :: timer_t
    logical :: enabled  = .true.          ! global on/off switch for begin/end accounting
    integer :: capacity = 0               ! max number of registered timers

    logical :: log_enabled    = .false.   ! enable writing to timing.log
    integer :: log_unit       = -1        ! unit for timing.log
    logical :: log_flush_each = .false.   ! flush after each dump (slow but safe)

    logical :: window_enabled  = .false.  ! enable periodic window dumps
    real(8) :: window_interval = 30.0d0   ! seconds between dumps
    real(8) :: next_dump_time  = 0.0d0    ! wallclock time threshold for next dump

    real(8) :: last_dump_time  = 0.0d0    ! wallclock of last window dump

  contains
    procedure, nopass :: timer_init
    procedure, nopass :: begin
    procedure, nopass :: end
    procedure, nopass :: reset
    procedure, nopass :: print
    procedure, nopass :: set_enabled

    procedure, nopass :: log_init
    procedure, nopass :: log_finalize
  end type timer_t
  !-----------------------------------------------------------------------------

  !---------------------- Global singleton timer object ------------------------
  type(timer_t), public, save :: timer
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
!> Initialise timer module (allocates thread storage, registry, locks)
!=============================================================================
  subroutine timer_init()
    integer :: nthreads, tid

    if (initialized) return

    !---------------------- Initialise locks -----------------------------------
    call reg_lock%init()
    call log_lock%init()
    !---------------------------------------------------------------------------

    !---------------------- Choose registry capacity ---------------------------
    timer%capacity = 4096
    !---------------------------------------------------------------------------

    !---------------------- Allocate registry name table -----------------------
    allocate(names(timer%capacity))
    names = ""
    ntimer = 0
    !---------------------------------------------------------------------------

    !---------------------- Allocate per-thread timer arrays -------------------
    nthreads = omp%max_threads
    allocate(timers(0:nthreads-1))

    !$omp parallel private(tid)
      tid = omp%thread_id()
      allocate(timers(tid)%calls(timer%capacity))
      allocate(timers(tid)%total(timer%capacity))
      timers(tid)%calls = 0_8
      timers(tid)%total = 0.0d0
      timers(tid)%depth = 0
      timers(tid)%stack_id = 0
      timers(tid)%stack_t0 = 0.0d0
      timers(tid)%child_acc = 0.0d0
    !$omp end parallel
    !---------------------------------------------------------------------------

    !-------------------- Reset window snapshot state --------------------------
    win_snap_initialized = .false.
    if (allocated(win_prev_calls)) deallocate(win_prev_calls)
    if (allocated(win_prev_total)) deallocate(win_prev_total)
    !---------------------------------------------------------------------------

    initialized = .true.

  end subroutine timer_init


!=============================================================================
!> Enable/disable timing globally
!> enabled : logical
!=============================================================================
  subroutine set_enabled(enabled)
    logical, intent(in) :: enabled
    timer%enabled = enabled
  end subroutine set_enabled


!=============================================================================
!> Reset all timer counters and stacks (keeps registry)
!=============================================================================
  subroutine reset()
    integer :: tid
    real(8) :: tnow

    if (.not. allocated(timers)) call timer%timer_init()

    !-------------------- Reset per-thread state -------------------------------
    if (omp%in_parallel()) then
      tid = omp%thread_id()
      call reset_thread(tid)
    else
      !$omp parallel private(tid)
        tid = omp%thread_id()
        call reset_thread(tid)
      !$omp end parallel
    end if
    !---------------------------------------------------------------------------

    !-------------------- Reset window snapshot/time base ----------------------
    call log_lock%lock()

      tnow = wallclock()

      !----------------- Initialise delta baseline at t = now ------------------
      ! This makes the first dump print "last window_interval seconds" instead
      ! of using the first trigger to initialise the snapshot.
      if (.not. allocated(win_prev_calls)) then
        allocate(win_prev_calls(timer%capacity), win_prev_total(timer%capacity))
      end if
      win_prev_calls = 0_8
      win_prev_total = 0.0d0
      win_snap_initialized = .true.
      !---------------------------------------------------------------------------

      timer%last_dump_time = tnow
      timer%next_dump_time = tnow + timer%window_interval

    call log_lock%unlock()
    !---------------------------------------------------------------------------

  end subroutine reset



!=============================================================================
!> Reset single thread timer state
!> tid : thread id
!=============================================================================
  subroutine reset_thread(tid)
    integer, intent(in) :: tid

    timers(tid)%calls = 0_8
    timers(tid)%total = 0.0d0

    timers(tid)%depth = 0
    timers(tid)%stack_id = 0
    timers(tid)%stack_t0 = 0.0d0
    timers(tid)%child_acc = 0.0d0

  end subroutine reset_thread


!=============================================================================
!> Begin timer region (exclusive time)
!> name   : timer label
!> itimer : timer id (0 => auto-register)
!=============================================================================
  subroutine begin(name, itimer)
    character(len=*), intent(in) :: name
    integer,          intent(inout) :: itimer

    integer :: tid, id
    real(8) :: tnow

    if (.not. timer%enabled) return
    if (.not. allocated(timers)) call timer%timer_init()

    tid = omp%thread_id()

    !---------------------- Register timer name if needed ----------------------
    if (itimer <= 0) then
      id = get_timer_id(name)
      itimer = id
    else
      id = itimer
    end if
    !---------------------------------------------------------------------------

    if (id < 1 .or. id > timer%capacity) error stop "timer_mod: invalid timer id."

    tnow = wallclock()

    !---------------------- Push on per-thread stack ---------------------------
    timers(tid)%depth = timers(tid)%depth + 1
    if (timers(tid)%depth > MAX_DEPTH) error stop "timer_mod: MAX_DEPTH exceeded."

    timers(tid)%stack_id(timers(tid)%depth)  = id
    timers(tid)%stack_t0(timers(tid)%depth)  = tnow
    timers(tid)%child_acc(timers(tid)%depth) = 0.0d0
    !---------------------------------------------------------------------------

  end subroutine begin


!=============================================================================
!> End timer region (exclusive time, proper nesting enforced with warnings)
!> itimer : timer id
!=============================================================================
  subroutine end(itimer)
    integer, intent(in) :: itimer

    integer :: tid, depth, id
    real(8) :: tnow, dt, excl

    if (.not. timer%enabled) return
    if (.not. allocated(timers)) call timer%timer_init()

    tid = omp%thread_id()
    id  = itimer

    if (id < 1 .or. id > timer%capacity) then
      call warn_mismatch("end() called with invalid itimer id.")
      return
    end if

    depth = timers(tid)%depth
    if (depth <= 0) then
      call warn_mismatch("end() called but stack empty.")
      return
    end if

    if (timers(tid)%stack_id(depth) /= id) then
      call warn_mismatch("mismatched nesting: end(id) != last begin(id).")
      return
    end if

    tnow = wallclock()
    dt   = tnow - timers(tid)%stack_t0(depth)
    excl = dt - timers(tid)%child_acc(depth)
    if (excl < 0.0d0) excl = 0.0d0

    !---------------------- Accumulate exclusive time --------------------------
    timers(tid)%calls(id) = timers(tid)%calls(id) + 1_8
    timers(tid)%total(id) = timers(tid)%total(id) + excl
    !---------------------------------------------------------------------------

    !---------------------- Pop and add child time to parent -------------------
    timers(tid)%depth = depth - 1
    if (timers(tid)%depth > 0) then
      timers(tid)%child_acc(timers(tid)%depth) = timers(tid)%child_acc(timers(tid)%depth) + dt
    end if
    !---------------------------------------------------------------------------

    !---------------------- Periodic window dump (tid=0 trigger) ---------------
    call maybe_window_dump()
    !---------------------------------------------------------------------------

  end subroutine end


!=============================================================================
!> Initialise timing log file (fresh) + reset omp_timer offset for this run
!> filename       : log filename (default "timing.log")
!> unit           : Fortran unit (default 98)
!> flush_each     : flush after each dump (slow)
!> window_enabled : enable periodic window dumps
!> window_interval: seconds between window dumps
!=============================================================================
  subroutine log_init(log_dir, filename, unit, flush_each, window_enabled, window_interval)
    !DEC$ ATTRIBUTES ALIAS:"log_init_" :: log_init
    character(len=*), intent(in), optional :: log_dir
    character(len=*), intent(in), optional :: filename
    integer,          intent(in), optional :: unit
    logical,          intent(in), optional :: flush_each
    logical,          intent(in), optional :: window_enabled
    real(8),          intent(in), optional :: window_interval
    character(len=256) :: dir
    character(len=256) :: full_fn
    character(len=256) :: fn
    logical :: log_unit_open

    if (.not. allocated(timers)) call timer%timer_init()

    fn = "timing.log"
    if (present(filename)) fn = filename

    timer%log_unit = 98
    if (present(unit)) timer%log_unit = unit

    timer%log_flush_each = .false.
    if (present(flush_each)) timer%log_flush_each = flush_each

    timer%window_enabled = .false.
    if (present(window_enabled)) timer%window_enabled = window_enabled

    timer%window_interval = 30.0d0
    if (present(window_interval)) timer%window_interval = window_interval

    dir = "logs"
    if (present(log_dir)) dir = log_dir

    !---------------------------- Ensure log dir --------------------------------
    call ensure_dir_exists(trim(dir))
    !---------------------------------------------------------------------------

    full_fn = trim(dir)//"/"//trim(fn)



    !------------------- Reset wallclock zero point for this run ---------------
    call omp_timer%init()
    !---------------------------------------------------------------------------

    !----------------------- Open timing log (fresh) ---------------------------
    call log_lock%lock()

      inquire(unit=timer%log_unit, opened=log_unit_open)
      if (log_unit_open) close(timer%log_unit)

      open(timer%log_unit, file=trim(full_fn), status="replace", action="write")
      timer%log_enabled = .true.

      win_snap_initialized = .false.
      if (allocated(win_prev_calls)) deallocate(win_prev_calls)
      if (allocated(win_prev_total)) deallocate(win_prev_total)
      allocate(win_prev_calls(timer%capacity), win_prev_total(timer%capacity))
      win_prev_calls = 0_8
      win_prev_total = 0.0d0

      write(timer%log_unit,'(a)') "============================================================"
      write(timer%log_unit,'(a)') "TIMING LOG START"
      write(timer%log_unit,'(a)') "============================================================"
      if (timer%log_flush_each) flush(timer%log_unit)

    call log_lock%unlock()
    !---------------------------------------------------------------------------

    !------------------------ Reset counters + window base ---------------------
    call timer%reset()
    !---------------------------------------------------------------------------

  end subroutine log_init


subroutine ensure_dir_exists(d)
  use iso_fortran_env, only: error_unit
  implicit none
  character(len=*), intent(in) :: d
  integer :: stat

  if (len_trim(d) == 0) return

  call execute_command_line('mkdir -p "'//trim(d)//'"', exitstat=stat)
  if (stat /= 0) then
    write(error_unit,'(a,i0)') 'TRACE: failed to create log directory, exitstat=', stat
  end if
end subroutine ensure_dir_exists


!=============================================================================
!> Finalise timing log (final summary + close)
!=============================================================================
  subroutine log_finalize()
    !DEC$ ATTRIBUTES ALIAS:"log_finalize_" :: log_finalize
    if (.not. timer%log_enabled) return

    call log_lock%lock()
      write(timer%log_unit,'(a)') ""
      write(timer%log_unit,'(a)') "--------------------- FINAL TIMER SUMMARY -------------------"
    call log_lock%unlock()

    call timer%print(unit=timer%log_unit, tid=0)

    call log_lock%lock()
      write(timer%log_unit,'(a)') "============================================================"
      write(timer%log_unit,'(a)') "TIMING LOG END"
      write(timer%log_unit,'(a)') "============================================================"
      if (timer%log_flush_each) flush(timer%log_unit)
      close(timer%log_unit)
      timer%log_enabled = .false.
      timer%log_unit    = -1
    call log_lock%unlock()

  end subroutine log_finalize


!=============================================================================
!> Check if we should dump a timing window; only tid=0 performs dump
!=============================================================================
  subroutine maybe_window_dump()
    real(8) :: tnow
    integer :: tid

    if (.not. timer%log_enabled) return
    if (.not. timer%window_enabled) return

    tid = omp%thread_id()
    if (tid /= 0) return

    tnow = wallclock()
    if (tnow < timer%next_dump_time) return

    timer%next_dump_time = tnow + timer%window_interval
    call dump_window(tnow)

  end subroutine maybe_window_dump


!=============================================================================
!> Dump delta timing summary for tid=0 since previous snapshot
!> tnow : current wallclock time
!=============================================================================
  subroutine dump_window(tnow)
    real(8), intent(in) :: tnow

    integer :: i
    integer(8), allocatable :: cur_calls(:), del_calls(:)
    real(8),    allocatable :: cur_total(:), del_total(:)
    real(8) :: t_rel, t_rel_r, dt, dt_r, other_time
    character(len=128) :: title

    if (.not. timer%log_enabled) return

    allocate(cur_calls(timer%capacity), cur_total(timer%capacity))
    allocate(del_calls(timer%capacity), del_total(timer%capacity))

    cur_calls = 0_8
    cur_total = 0.0d0

    !---------------------- Snapshot current totals (tid=0) --------------------
    do i = 1, ntimer
      cur_calls(i) = timers(0)%calls(i)
      cur_total(i) = timers(0)%total(i)
    end do
    !---------------------------------------------------------------------------

    !---------------------- Compute delta since previous snapshot --------------
    call log_lock%lock()

      if (.not. win_snap_initialized) then
        win_prev_calls = cur_calls
        win_prev_total = cur_total
        win_snap_initialized = .true.

        call log_lock%unlock()
        deallocate(cur_calls, cur_total, del_calls, del_total)
        timer%last_dump_time = tnow
        return
      end if

      do i = 1, ntimer
        del_calls(i) = cur_calls(i) - win_prev_calls(i)
        del_total(i) = cur_total(i) - win_prev_total(i)
      end do

      win_prev_calls = cur_calls
      win_prev_total = cur_total

    call log_lock%unlock()
    !---------------------------------------------------------------------------

    !---------------------- Compute window times -------------------------------
    t_rel   = tnow
    t_rel_r = dble(nint(10.0d0*t_rel))/10.0d0

    dt   = tnow - timer%last_dump_time
    dt_r = dble(nint(10.0d0*dt))/10.0d0

    timer%last_dump_time = tnow

    other_time = max(0.0d0, dt - sum(del_total(1:ntimer)))

    write(title,'(a,f0.1,a)') "Window timer summary for last ", dt_r, " seconds"
    !---------------------------------------------------------------------------

    !---------------------- Print delta summary (single ==== region) -----------
    call log_lock%lock()
      write(timer%log_unit,'(a)') ""
      write(timer%log_unit,'(a)') "============================================================"
      write(timer%log_unit,'(a,f0.1,a)') "Time = ", t_rel_r, " s"
      write(timer%log_unit,'(a)') trim(title)
      write(timer%log_unit,'(a)') "------------------------------------------------------------"
    call log_lock%unlock()

    call print_from_arrays(unit=timer%log_unit, calls=del_calls, total=del_total, title="", total_ref=dt, other_time=other_time, show_total=.false.)
    !---------------------------------------------------------------------------

    deallocate(cur_calls, cur_total, del_calls, del_total)
  end subroutine dump_window


!=============================================================================
!> Print cumulative timer summary
!> unit : output unit (default = stdout)
!> tid  : if present -> print this tid only, else SUM over threads
!=============================================================================
  subroutine print(unit, tid)
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: tid

    integer :: u, i, nthreads, t, tsel
    integer(8), allocatable :: calls(:)
    real(8),    allocatable :: total(:)

    u = 6
    if (present(unit)) u = unit
    if (.not. allocated(timers)) call timer%timer_init()

    nthreads = size(timers)

    allocate(calls(timer%capacity), total(timer%capacity))
    calls = 0_8
    total = 0.0d0

    !------------------ Select aggregation mode --------------------------------
    if (present(tid)) then
      tsel = tid
      if (tsel < 0 .or. tsel > nthreads-1) then
        error stop "timer_mod: print() called with invalid tid."
      end if

      do i = 1, ntimer
        calls(i) = timers(tsel)%calls(i)
        total(i) = timers(tsel)%total(i)
      end do

      call print_from_arrays(unit=u, calls=calls, total=total, title="Timer summary (exclusive wall time, tid = 0)", show_total=.true.)
    else
      do t = 0, nthreads-1
        do i = 1, ntimer
          calls(i) = calls(i) + timers(t)%calls(i)
          total(i) = total(i) + timers(t)%total(i)
        end do
      end do

      call print_from_arrays(unit=u, calls=calls, total=total, title="Timer summary (exclusive wall time, SUM over threads)", show_total=.true.)
    end if
    !---------------------------------------------------------------------------

    deallocate(calls, total)
  end subroutine print


!=============================================================================
!> Register timer name, returning id (1..capacity)
!=============================================================================
  integer function get_timer_id(name) result(id)
    character(len=*), intent(in) :: name
    integer :: i

    call reg_lock%lock()

    !---------------------- Lookup existing name ------------------------------
    do i = 1, ntimer
      if (trim(names(i)) == trim(name)) then
        id = i
        call reg_lock%unlock()
        return
      end if
    end do
    !---------------------------------------------------------------------------

    !---------------------- Register new timer --------------------------------
    ntimer = ntimer + 1
    if (ntimer > timer%capacity) then
      call reg_lock%unlock()
      error stop "timer_mod: capacity exceeded."
    end if

    names(ntimer) = trim(name)
    id = ntimer
    !---------------------------------------------------------------------------

    call reg_lock%unlock()

  end function get_timer_id


!=============================================================================
!> Print timer summary (filtered/sorted) from explicit arrays
!> unit       : output unit
!> calls      : calls per timer (1..capacity)
!> total      : total time per timer (1..capacity)
!> title      : header title (optional; empty string => no title line)
!> total_ref  : optional reference total (e.g. window wall time)
!> other_time : optional "other" time to include as separate row
!> show_total : print "total[s]" line (default true)
!=============================================================================
  subroutine print_from_arrays(unit, calls, total, title, total_ref, other_time, show_total)
    integer, intent(in) :: unit
    integer(8), intent(in) :: calls(:)
    real(8),    intent(in) :: total(:)
    character(len=*), intent(in) :: title
    real(8),    intent(in), optional :: total_ref
    real(8),    intent(in), optional :: other_time
    logical,    intent(in), optional :: show_total

    integer :: i, k, nkeep, nitem, i_other
    real(8) :: total_sum, avg, oth
    real(8), allocatable :: pct(:), total_w(:)
    integer(8), allocatable :: calls_w(:)
    integer, allocatable :: idx(:), keep(:)
    character(len=NAME_LEN), allocatable :: names_w(:)
    logical :: do_total

    oth = 0.0d0
    if (present(other_time)) oth = max(0.0d0, other_time)

    nitem = ntimer
    if (oth > 0.0d0) nitem = ntimer + 1

    do_total = .true.
    if (present(show_total)) do_total = show_total

    allocate(pct(nitem), idx(nitem), keep(nitem), calls_w(nitem), total_w(nitem), names_w(nitem))
    pct     = 0.0d0
    idx     = 0
    keep    = 0
    calls_w = 0_8
    total_w = 0.0d0
    names_w = ""

    !---------------------- Build working arrays (optionally add "other") ------
    do i = 1, ntimer
      idx(i)        = i
      calls_w(i)    = calls(i)
      total_w(i)    = total(i)
      names_w(i)    = names(i)
    end do

    if (oth > 0.0d0) then
      i_other          = ntimer + 1
      idx(i_other)     = i_other
      calls_w(i_other) = 1_8
      total_w(i_other) = oth
      names_w(i_other) = "other"
    end if
    !---------------------------------------------------------------------------

    !---------------------- Compute total sum ----------------------------------
    if (present(total_ref)) then
      total_sum = max(0.0d0, total_ref)
    else
      total_sum = 0.0d0
      do i = 1, nitem
        if (calls_w(i) > 0_8) total_sum = total_sum + total_w(i)
      end do
    end if
    !---------------------------------------------------------------------------

    call log_lock%lock()

      if (len_trim(title) > 0) then
        write(unit,'(a)') "============================================================"
        write(unit,'(a)') trim(title)
      end if

      if (do_total) then
        write(unit,'(a,es14.6)') "  total[s] = ", total_sum
      end if

      write(unit,'(a)') "------------------------------------------------------------"
      write(unit,'(a)') "  name                               total[s]       avg[s]      %     calls"
      write(unit,'(a)') "------------------------------------------------------------"

    call log_lock%unlock()

    if (total_sum <= 0.0d0) then
      call log_lock%lock()
      write(unit,'(a)') "  (no timers recorded in this window)"
      write(unit,'(a)') "============================================================"
      call log_lock%unlock()
      deallocate(pct, idx, keep, calls_w, total_w, names_w)
      return
    end if

    !---------------------- Percentages ----------------------------------------
    do i = 1, nitem
      if (calls_w(i) > 0_8) then
        pct(i) = 100.0d0 * total_w(i) / total_sum
      else
        pct(i) = 0.0d0
      end if
    end do
    !---------------------------------------------------------------------------

    !---------------------- Sort by percentage (descending) --------------------
    call sort_desc_by_real(idx, pct)
    !---------------------------------------------------------------------------

    !---------------------- Filter < 0.5% -------------------------------------
    nkeep = 0
    do k = 1, nitem
      i = idx(k)
      if (calls_w(i) <= 0_8) cycle
      if (pct(i) < 0.5d0) cycle
      nkeep = nkeep + 1
      keep(nkeep) = i
    end do
    !---------------------------------------------------------------------------

    !---------------------- Print table ----------------------------------------
    do k = 1, nkeep
      i = keep(k)
      if (calls_w(i) > 0_8) then
        avg = total_w(i) / dble(calls_w(i))
      else
        avg = 0.0d0
      end if

      call log_lock%lock()
      write(unit,'(1x,a35,1x,es12.4,1x,es10.3,1x,f6.2,1x,i10)') &
        trim(names_w(i)), total_w(i), avg, pct(i), calls_w(i)
      call log_lock%unlock()
    end do

    call log_lock%lock()
      write(unit,'(a)') "============================================================"
      if (timer%log_flush_each) flush(unit)
    call log_lock%unlock()
    !---------------------------------------------------------------------------

    deallocate(pct, idx, keep, calls_w, total_w, names_w)
  end subroutine print_from_arrays


!=============================================================================
!> Warning printer for nesting mismatch (does not abort)
!=============================================================================
  subroutine warn_mismatch(msg)
    character(len=*), intent(in) :: msg
    call log_lock%lock()
      write(*,'(a)') "TIMER WARNING: "//trim(msg)
    call log_lock%unlock()
  end subroutine warn_mismatch


!=============================================================================
!> Sort idx by val(idx(i)) descending (insertion sort; small arrays)
!=============================================================================
  subroutine sort_desc_by_real(idx, val)
    integer, intent(inout) :: idx(:)
    real(8), intent(in)    :: val(:)

    integer :: i, j, key
    real(8) :: keyv

    !---------------------- Insertion sort (small ntimer) ----------------------
    do i = 2, size(idx)
      key  = idx(i)
      keyv = val(key)

      j = i - 1
      do while (j >= 1 .and. val(idx(j)) < keyv)
        idx(j+1) = idx(j)
        j = j - 1
      end do
      idx(j+1) = key
    end do
    !---------------------------------------------------------------------------

  end subroutine sort_desc_by_real

END MODULE timer_mod
