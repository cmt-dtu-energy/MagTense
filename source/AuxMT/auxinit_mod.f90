MODULE auxInit_mod
  USE omp_mod,   only : omp
  USE timer_mod, only : timer
  USE trace_mod, only : trace
  implicit none
  private

  !------------------------- Public auxiliary wrapper -------------------------
  type, public :: auxInit_t
    character(len=256) :: log_dir        = "logs"        ! directory for auxiliary log files
    character(len=256) :: timer_log_file = "timer.log"   ! timer log file name
    character(len=256) :: trace_log_file = "trace.log"   ! trace log file name

    logical :: window_enabled = .true.                   ! enable periodic timer windows
    logical :: trace_enabled  = .false.                  ! enable trace logging
    logical :: flush_each     = .true.                   ! flush log output after each write

    real(8) :: window_interval = 30.0d0                  ! seconds between timer windows
    integer :: trace_verbose  = 1                        ! trace verbosity level
  end type auxInit_t
  !----------------------------------------------------------------------------

  public :: initAux
  !----------------------------------------------------------------------------

  !---------------------- Global singleton auxiliary object --------------------
  type(auxInit_t), public, save :: auxInit
  !----------------------------------------------------------------------------

CONTAINS

!=============================================================================
!> Initialise the auxiliary OpenMP, timer and trace modules
!> log_dir         : optional directory for auxiliary log files
!> timer_log_file  : optional timer log file name
!> trace_log_file  : optional trace log file name
!> window_enabled  : optional integer flag for periodic timer windows
!> window_interval : optional seconds between timer windows
!> trace_enabled   : optional integer flag for trace logging
!> flush_each      : optional integer flag for flushing after each write
!> trace_verbose   : optional trace verbosity level
!=============================================================================
  subroutine initAux(self, log_dir, timer_log_file, trace_log_file, window_enabled, &
      window_interval, trace_enabled, flush_each, trace_verbose)
    !DEC$ ATTRIBUTES ALIAS:"initaux_" :: initAux
    type(auxInit_t), intent(inout) :: self
    character(len=*), intent(in), optional :: log_dir, timer_log_file, trace_log_file
    integer, intent(in), optional :: window_enabled, trace_enabled, flush_each
    real(8), intent(in), optional :: window_interval
    integer, intent(in), optional :: trace_verbose

    !---------------------- Store the requested input parameters ---------------
    if (present(log_dir)) self%log_dir = log_dir
    if (present(timer_log_file)) self%timer_log_file = timer_log_file
    if (present(trace_log_file)) self%trace_log_file = trace_log_file
    if (present(window_enabled)) self%window_enabled = merge(.true., .false., window_enabled /= 0)
    if (present(trace_enabled)) self%trace_enabled = merge(.true., .false., trace_enabled /= 0)
    if (present(flush_each)) self%flush_each = merge(.true., .false., flush_each /= 0)
    if (present(window_interval)) self%window_interval = window_interval
    if (present(trace_verbose)) self%trace_verbose = trace_verbose
    !---------------------------------------------------------------------------

    !---------------------- Initialise auxiliary modules -----------------------
    call omp%init()
    if (self%trace_enabled) call omp%info()
    call timer%log_init(trim(self%log_dir), trim(self%timer_log_file), &
      flush_each=self%flush_each, window_enabled=self%window_enabled, &
      window_interval=self%window_interval)
    call trace%trace_init(trim(self%log_dir), trim(self%trace_log_file), &
      enabled=self%trace_enabled, unit=97, flush_each=self%flush_each, &
      verbose=self%trace_verbose)
    !---------------------------------------------------------------------------

  end subroutine initAux

END MODULE auxInit_mod
