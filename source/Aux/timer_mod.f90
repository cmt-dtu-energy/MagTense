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

  !-------------------------- Thread-local timer state -------------------------
  type :: timer_thread_t
    integer(8), allocatable :: calls(:)
    real(8),    allocatable :: start(:)
    real(8),    allocatable :: total(:)
    integer,    allocatable :: stack(:)
    integer :: istack = 0
  end type timer_thread_t
  !-----------------------------------------------------------------------------

  !---------------------- Public timer object interface ------------------------
  type, public :: timer_t
    logical :: enabled  = .true.
    integer :: capacity = 0
  contains
    procedure, nopass :: init
    procedure, nopass :: begin
    procedure, nopass :: end
    procedure, nopass :: reset
    procedure, nopass :: print
    procedure, nopass :: set_enabled
  end type timer_t
  !-----------------------------------------------------------------------------

  !---------------------- Global singleton timer object ------------------------
  type(timer_t), public, save :: timer
  !-----------------------------------------------------------------------------

  !--------------------------- Shared timer registry ---------------------------
  integer, save :: ntimer = 0
  character(len=NAME_LEN), allocatable, save :: names(:)
  !-----------------------------------------------------------------------------

  !-------------------------- Per-thread timer storage -------------------------
  type(timer_thread_t), allocatable, save :: timers(:)
  !-----------------------------------------------------------------------------

  !------------------- Lock protecting registry (ntimer/names) -----------------
  type(omp_lock_t), save :: reg_lock
  !-----------------------------------------------------------------------------

  !-------------------------- Module initialisation flag -----------------------
  logical, save :: initialized = .false.
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
!> Initialise timer module
!> capacity        : maximum number of timers
!> enabled         : enable/disable timing
!> sync_omp_offset : sync omp_timer offset across threads
!=============================================================================
  subroutine init(capacity, enabled, sync_omp_offset)
    integer, intent(in), optional :: capacity
    logical, intent(in), optional :: enabled
    logical, intent(in), optional :: sync_omp_offset
    integer :: cap

    cap = 512
    if (present(capacity)) cap = capacity

    timer%capacity = cap
    timer%enabled  = .true.
    if (present(enabled)) timer%enabled = enabled

    !------------------- Initialise omp timer source ---------------------------
    if (present(sync_omp_offset)) omp_timer%sync = sync_omp_offset
    call omp_timer%init()
    !---------------------------------------------------------------------------

    !-------- Initialise registry lock once (safe in parallel) -----------------
    !$omp critical(timer_mod_init_lock)
      call reg_lock%init()
      initialized = .true.
    !$omp end critical(timer_mod_init_lock)
    !---------------------------------------------------------------------------

    !--------------------- Allocate shared registry ----------------------------
    if (allocated(names)) deallocate(names)
    allocate(names(cap))
    names  = ""
    ntimer = 0
    !---------------------------------------------------------------------------

    !------------------- Allocate per-thread storage ---------------------------
    if (allocated(timers)) deallocate(timers)

    if (omp%in_parallel()) then
      call init_parallel_storage(cap)
    else
      !$omp parallel
        call init_parallel_storage(cap)
      !$omp end parallel
    end if
    !---------------------------------------------------------------------------

  end subroutine init


!=============================================================================
!> Allocate and clear per-thread timer arrays
!> cap : timer capacity
!=============================================================================
  subroutine init_parallel_storage(cap)
    integer, intent(in) :: cap
    integer :: nthreads, tid

    !---------------------- Allocate timers array ------------------------------
    !$omp single
      nthreads = omp%numthreads()
      allocate(timers(0:nthreads-1))
    !$omp end single
    !---------------------------------------------------------------------------

    tid = omp%thread_id()

    !------------------ Allocate thread-local arrays ---------------------------
    if (allocated(timers(tid)%calls)) then
      deallocate(timers(tid)%calls, timers(tid)%start, timers(tid)%total, timers(tid)%stack)
    end if

    allocate(timers(tid)%calls(cap), timers(tid)%start(cap), timers(tid)%total(cap))
    allocate(timers(tid)%stack(MAX_DEPTH))

    timers(tid)%calls  = 0_8
    timers(tid)%start  = 0.0d0
    timers(tid)%total  = 0.0d0
    timers(tid)%stack  = 0
    timers(tid)%istack = 0
    !---------------------------------------------------------------------------

  end subroutine init_parallel_storage


!=============================================================================
!> Enable or disable timing
!> enabled : logical flag
!=============================================================================
  subroutine set_enabled(enabled)
    logical, intent(in) :: enabled
    timer%enabled = enabled
  end subroutine set_enabled


!=============================================================================
!> Reset all timers (keeps registry)
!=============================================================================
  subroutine reset()
    integer :: tid

    if (.not. allocated(timers)) call timer%init()

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

  end subroutine reset


!=============================================================================
!> Reset single thread timer state
!> tid : thread id
!=============================================================================
  subroutine reset_thread(tid)
    integer, intent(in) :: tid
    timers(tid)%calls  = 0_8
    timers(tid)%start  = 0.0d0
    timers(tid)%total  = 0.0d0
    timers(tid)%stack  = 0
    timers(tid)%istack = 0
  end subroutine reset_thread


!=============================================================================
!> Register a timer name and return its id
!> label : timer name
!=============================================================================
  integer function register_timer(label) result(itimer)
    character(len=*), intent(in) :: label
    integer :: i
    character(len=NAME_LEN) :: lab

    if (.not. initialized .or. .not. allocated(names) .or. .not. allocated(timers)) call timer%init()

    lab = ""
    lab(1:min(len_trim(label), NAME_LEN)) = label(1:min(len_trim(label), NAME_LEN))

    !------------------- Protected registry access -----------------------------
    call reg_lock%lock()

    do i = 1, ntimer
      if (trim(names(i)) == trim(lab)) then
        itimer = i
        call reg_lock%unlock()
        return
      end if
    end do

    if (ntimer >= timer%capacity) then
      call reg_lock%unlock()
      error stop "timer_mod: timer registry capacity exceeded."
    end if

    ntimer = ntimer + 1
    names(ntimer) = lab
    itimer = ntimer

    call reg_lock%unlock()
    !---------------------------------------------------------------------------

  end function register_timer


!=============================================================================
!> Begin timing a region
!> label  : timer name
!> itimer : timer id (0 = auto-register)
!=============================================================================
  subroutine begin(label, itimer)
    character(len=*), intent(in) :: label
    integer, intent(inout)       :: itimer
    integer :: tid, otimer
    real(8) :: t

    if (.not. allocated(timers)) call timer%init()
    if (itimer <= 0) itimer = 0
    if (itimer == 0) itimer = register_timer(label)
    if (.not. timer%enabled) return

    tid = 0
    if (omp%in_parallel()) tid = omp%thread_id()

    t = wallclock()

    !---------------- Pause current active timer -------------------------------
    if (timers(tid)%istack > 0) then
      otimer = timers(tid)%stack(timers(tid)%istack)
      timers(tid)%total(otimer) = timers(tid)%total(otimer) + (t - timers(tid)%start(otimer))
    end if
    !---------------------------------------------------------------------------

    !---------------- Push new timer on stack ----------------------------------
    if (timers(tid)%istack >= MAX_DEPTH) error stop "timer_mod: MAX_DEPTH exceeded."
    timers(tid)%istack = timers(tid)%istack + 1
    timers(tid)%stack(timers(tid)%istack) = itimer

    timers(tid)%start(itimer) = t
    timers(tid)%calls(itimer) = timers(tid)%calls(itimer) + 1_8
    !---------------------------------------------------------------------------

  end subroutine begin


!=============================================================================
!> End timing a region
!> itimer : timer id
!=============================================================================
  subroutine end(itimer)
    integer, intent(in) :: itimer
    integer :: tid, top, parent
    real(8) :: t

    if (.not. allocated(timers)) call timer%init()
    if (.not. timer%enabled) return

    tid = 0
    if (omp%in_parallel()) tid = omp%thread_id()

    if (timers(tid)%istack <= 0) error stop "timer_mod: end() with empty stack."

    top = timers(tid)%stack(timers(tid)%istack)
    if (top /= itimer) error stop "timer_mod: end() not matching top of stack."

    t = wallclock()

    !------------------ Accumulate timer --------------------------------------
    timers(tid)%total(itimer) = timers(tid)%total(itimer) + (t - timers(tid)%start(itimer))
    !---------------------------------------------------------------------------

    !------------------ Pop stack ---------------------------------------------
    timers(tid)%istack = timers(tid)%istack - 1
    !---------------------------------------------------------------------------

    !------------------ Resume parent timer -----------------------------------
    if (timers(tid)%istack > 0) then
      parent = timers(tid)%stack(timers(tid)%istack)
      timers(tid)%start(parent) = t
    end if
    !---------------------------------------------------------------------------

  end subroutine end


!=============================================================================
!> Print timer summary
!> unit : output unit (default = stdout)
!> tid  : thread id to print (optional; if absent -> aggregate over threads)
!=============================================================================
  subroutine print(unit, tid)
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: tid
    integer :: u, i, tsel, nthreads, t
    integer(8), allocatable :: sum_calls(:)
    real(8),    allocatable :: sum_total(:)
    real(8) :: avg
    logical :: one_thread

    u = 6
    if (present(unit)) u = unit
    if (.not. allocated(timers)) call timer%init()

    nthreads = size(timers)

    allocate(sum_calls(timer%capacity), sum_total(timer%capacity))
    sum_calls = 0_8
    sum_total = 0.0d0

    one_thread = .false.

    !------------------ Select aggregation mode --------------------------------
    if (present(tid)) then
      one_thread = .true.
      tsel = tid
      if (tsel < 0 .or. tsel > nthreads-1) then
        error stop "timer_mod: print() called with invalid tid."
      end if

      do i = 1, ntimer
        sum_calls(i) = timers(tsel)%calls(i)
        sum_total(i) = timers(tsel)%total(i)
      end do

    else
      !---------------- Aggregate SUM over threads -----------------------------
      do t = 0, nthreads-1
        do i = 1, ntimer
          sum_calls(i) = sum_calls(i) + timers(t)%calls(i)
          sum_total(i) = sum_total(i) + timers(t)%total(i)
        end do
      end do
      !------------------------------------------------------------------------
    end if
    !---------------------------------------------------------------------------

    write(u,'(a)') "============================================================"
    if (one_thread) then
      write(u,'(a,i0)') "Timer summary (exclusive wall time, thread tid = ", tsel
    else
      write(u,'(a)') "Timer summary (exclusive wall time, aggregated SUM over threads)"
      write(u,'(a,i0)') "  threads = ", nthreads
    end if
    write(u,'(a)') "------------------------------------------------------------"
    write(u,'(a)') "  id   calls        total[s]       avg[s]    name"
    write(u,'(a)') "------------------------------------------------------------"

    do i = 1, ntimer
      if (sum_calls(i) > 0_8) then
        avg = sum_total(i) / dble(sum_calls(i))
      else
        avg = 0.0d0
      end if
      write(u,'(i4,1x,i10,1x,es14.6,1x,es12.4,2x,a)') &
        i, sum_calls(i), sum_total(i), avg, trim(names(i))
    end do

    write(u,'(a)') "============================================================"

    deallocate(sum_calls, sum_total)
  end subroutine print

END MODULE timer_mod
