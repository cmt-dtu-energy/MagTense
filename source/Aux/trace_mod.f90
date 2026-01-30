!================================================
! trace_mod.f90  (OpenMP-capable trace logging, no MPI)
!================================================
MODULE trace_mod
  USE omp_timer_mod, only : wallclock
  USE omp_lock_mod,  only : omp_lock_t
  USE omp_mod,       only : omp
  USE timer_mod,     only : timer
  implicit none
  private

  integer, parameter :: ID_LEN    = 128
  integer, parameter :: MAX_DEPTH = 64

  !---------------------------- Thread-local trace state -----------------------
  type :: trace_thread_t
    integer :: level = 0                         ! indentation depth for this thread
    character(len=ID_LEN) :: stack(MAX_DEPTH)    ! active begin stack for this thread
  end type trace_thread_t
  !-----------------------------------------------------------------------------

  !----------------------------- Public trace object ---------------------------
  type, public :: trace_t
    logical :: enabled    = .false.              ! master on/off switch
    integer :: unit       = -1                   ! file unit used for trace output
    logical :: flush_each = .false.              ! flush after each write
    integer :: verbose   = 0                   ! verbosity level (not used yet)

    integer :: nthreads = 1                      ! number of threads (cached)
    type(trace_thread_t), allocatable :: t(:)    ! per-thread state (0..nthreads-1)

    type(omp_lock_t) :: io_lock                  ! lock protecting file writes
    logical :: initialized = .false.             ! has init() been called?
  contains
    procedure, nopass :: init
    procedure, nopass :: finalize
    procedure, nopass :: set_enabled
    procedure, nopass :: begin
    procedure, nopass :: end
    procedure, nopass :: tag
  end type trace_t
  !-----------------------------------------------------------------------------

  !---------------------- Global singleton trace object ------------------------
  type(trace_t), public, save :: trace
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
!> Initialise trace logging
!> filename   : output filename (default "trace.log")
!> enabled    : enable/disable tracing
!> unit       : Fortran unit (default 97)
!> flush_each : flush file after each line (slow but safe)
!=============================================================================
  subroutine init(log_dir, filename, enabled, unit, flush_each, verbose)
    character(len=*), intent(in), optional :: log_dir
    character(len=*), intent(in), optional :: filename
    logical,          intent(in), optional :: enabled
    integer,          intent(in), optional :: unit
    logical,          intent(in), optional :: flush_each
    integer,          intent(in), optional :: verbose
    character(len=256) :: dir
    character(len=256) :: full_fn
    character(len=256) :: fn

    trace%enabled = .false.
    if (present(enabled)) trace%enabled = enabled

    trace%flush_each = .false.
    if (present(flush_each)) trace%flush_each = flush_each

    trace%unit = 97
    if (present(unit)) trace%unit = unit

    fn = "trace.log"
    if (present(filename)) fn = filename

    if (present(verbose)) trace%verbose = verbose

    dir = "logs"
    if (present(log_dir)) dir = log_dir
    full_fn = trim(dir)//"/"//trim(fn)

    !------------------------- Initialise locks once ---------------------------
    call trace%io_lock%init()
    !---------------------------------------------------------------------------

    !-------------------- Allocate per-thread state ----------------------------
    if (allocated(trace%t)) deallocate(trace%t)

    trace%nthreads = omp%max_threads
    allocate(trace%t(0:trace%nthreads-1))
    trace%t(:)%level = 0
    !---------------------------------------------------------------------------

    !---------------------------- Open log file --------------------------------
    if (trace%enabled) then
      open(trace%unit, file=trim(full_fn), status="replace", action="write", position="append")

      call trace_write_line("INFO ", "TRACE START", 0, force_master=.true.)

      if (trace%flush_each) flush(trace%unit)
    end if
    !---------------------------------------------------------------------------

    trace%initialized = .true.
  end subroutine init


!=============================================================================
!> Finalise trace logging (writes timer summary and closes file)
!=============================================================================
  subroutine finalize()
    if (.not. trace%initialized) return

    !---------------------- Write timer summary -------------------------------
    if (trace%enabled .and. trace%unit > 0) then
      call trace_write_line("INFO ", "TIMER SUMMARY", 0, force_master=.true.)
      call timer%print(unit=trace%unit, tid=0)
      if (trace%flush_each) flush(trace%unit)
    end if
    !---------------------------------------------------------------------------

    !-------------------------- Close down safely ------------------------------
    if (trace%enabled .and. trace%unit > 0) then
      call trace_write_line("INFO ", "TRACE END", 0, force_master=.true.)
      if (trace%flush_each) flush(trace%unit)
      close(trace%unit)
    end if
    !---------------------------------------------------------------------------

    trace%enabled      = .false.
    trace%unit         = -1

    if (allocated(trace%t)) deallocate(trace%t)

    call trace%io_lock%destroy()

    trace%nthreads     = 1
    trace%initialized = .false.

  end subroutine finalize


!=============================================================================
!> Enable/disable tracing (does not open/close files)
!> enabled : logical flag
!=============================================================================
  subroutine set_enabled(enabled)
    logical, intent(in) :: enabled
    trace%enabled = enabled
  end subroutine set_enabled


!=============================================================================
!> Begin a traced region (thread-safe)
!> id     : region label
!> itimer : optional timer id for aggregation (0 = auto-register)
!> verbose : optional verbose level
!=============================================================================
  subroutine begin(id, itimer, verbose)
    character(len=*), intent(in)              :: id
    integer,          intent(inout), optional :: itimer
    integer,          intent(in), optional    :: verbose
    integer :: verbose_in

    integer :: tid
    character(len=ID_LEN) :: lab

    if (.not. trace%initialized) call trace%init(enabled=.false.)
    
    !--------------------------- Verbosity check ------------------------------
    ! only proceed if the message verbosity is less than or equal to the trace verbosity
    if (present(verbose)) then
      verbose_in = verbose
    else
      verbose_in = 1
    end if
    if (verbose_in > trace%verbose) return
    !---------------------------------------------------------------------------


    tid = omp%thread_id()

    lab = ""
    lab(1:min(len_trim(id), ID_LEN)) = id(1:min(len_trim(id), ID_LEN))

    !----------------------------- Timer coupling ------------------------------
    if (present(itimer)) then
      call timer%begin(lab, itimer)
    end if
    !---------------------------------------------------------------------------

    !---------------------- Update per-thread stack ----------------------------
    trace%t(tid)%level = trace%t(tid)%level + 1
    if (trace%t(tid)%level > MAX_DEPTH) error stop "trace_mod: MAX_DEPTH exceeded."
    trace%t(tid)%stack(trace%t(tid)%level) = lab
    !---------------------------------------------------------------------------

    if (.not. trace%enabled) return

    ! Indent is one level less than current depth (children are deeper).
    call trace_write_line("BEGIN", lab, trace%t(tid)%level-1)

  end subroutine begin


!=============================================================================
!> End a traced region (thread-safe, enforces proper nesting)
!> id     : region label (must match last begin on this thread)
!> itimer : optional timer id for aggregation
!> verbose : optional verbose level
!=============================================================================
  subroutine end(id, itimer, verbose)
    character(len=*), intent(in)           :: id
    integer,          intent(in), optional :: itimer
    integer,          intent(in), optional :: verbose
    integer :: verbose_in
    integer :: tid
    character(len=ID_LEN) :: lab, want

    if (.not. trace%initialized) call trace%init(enabled=.false.)

   !--------------------------- Verbosity check ------------------------------
    ! only proceed if the message verbosity is less than or equal to the trace verbosity
    if (present(verbose)) then
      verbose_in = verbose
    else
      verbose_in = 1
    end if
    if (verbose_in > trace%verbose) return

    !---------------------------------------------------------------------------


    tid = omp%thread_id()

    lab = ""
    lab(1:min(len_trim(id), ID_LEN)) = id(1:min(len_trim(id), ID_LEN))

    !----------------------------- Timer coupling ------------------------------
    if (present(itimer)) then
      call timer%end(itimer)
    end if
    !---------------------------------------------------------------------------

    !------------------------ Check nesting correctness ------------------------
    if (trace%t(tid)%level <= 0) then
      error stop "trace_mod: end() called but stack empty."
    end if

    want = trace%t(tid)%stack(trace%t(tid)%level)
    if (trim(want) /= trim(lab)) then
      error stop "trace_mod: mismatched end() (id != last begin())."
    end if
    !---------------------------------------------------------------------------

    if (trace%enabled) then
      call trace_write_line("END  ", lab, trace%t(tid)%level-1)
    end if

    trace%t(tid)%level = trace%t(tid)%level - 1

  end subroutine end


!=============================================================================
!> Write a tag message at current indentation level (thread-safe)
!> message : freeform message
!=============================================================================
  subroutine tag(message)
    character(len=*), intent(in) :: message

    integer :: tid
    character(len=ID_LEN) :: msg

    if (.not. trace%initialized) call trace%init(enabled=.false.)
    if (.not. trace%enabled) return

    tid = omp%thread_id()

    msg = ""
    msg(1:min(len_trim(message), ID_LEN)) = message(1:min(len_trim(message), ID_LEN))

    call trace_write_line("TAG  ", msg, trace%t(tid)%level)

  end subroutine tag


!=============================================================================
!> Write one trace line (serialised with io_lock)
!> kind         : 5-char code ("BEGIN","END  ","TAG  ","INFO ")
!> label        : message text (not indented)
!> level        : indentation level (2 spaces per level)
!> force_master : if true, use tid=0 in output
!=============================================================================
  subroutine trace_write_line(kind, label, level, force_master)
    character(len=*), intent(in) :: kind
    character(len=*), intent(in) :: label
    integer,          intent(in) :: level
    logical,          intent(in), optional :: force_master

    integer :: tid
    logical :: fm

    if (.not. trace%enabled) return

    fm = .false.
    if (present(force_master)) fm = force_master

    tid = omp%thread_id()
    if (fm) tid = 0

    !------------------------ Serialize file output ----------------------------
    call trace%io_lock%lock()

      write(trace%unit,'(es14.6,1x,"tid=",i0,3x,a,a,1x,a)') &
          wallclock(), tid, repeat("  ", max(0,level)), adjustl(kind), trim(label)

      if (trace%flush_each) flush(trace%unit)

    call trace%io_lock%unlock()
    !---------------------------------------------------------------------------

  end subroutine trace_write_line

END MODULE trace_mod
