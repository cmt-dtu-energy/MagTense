!================================================
! Copied from ITA-Solar/Dispatch-Bifrost/ - Michael Haahr
!================================================

!*******************************************************************************
!> Support tic/toc timing, as in MATLAB, and accurate wallclock() function.
!> The timing is generally much more accurate if this module is compiled with
!> OMP active.
!*******************************************************************************
MODULE omp_timer_mod
  USE omp_mod
  implicit none
  private
  real(8), external:: omp_get_wtime
  type omp_timer_t
    logical:: sync=.true.
  contains
    procedure:: init
    procedure, nopass:: get
    procedure, nopass:: set
  end type
  type(omp_timer_t), public:: omp_timer
  real(8), save:: offset=0d0
  !omp threadprivate (offset)
  public wallclock
CONTAINS

!===============================================================================
FUNCTION wallclock() result (time)
  real(8):: time
  !.............................................................................
  time = omp_get_wtime()
  if (offset == 0d0) then
    offset=time
  end if
  time = time-offset
END FUNCTION wallclock

!===============================================================================
!> Iniitallize offset for all threads
!===============================================================================
SUBROUTINE init (self)
  class(omp_timer_t):: self
  !-----------------------------------------------------------------------------
  !$omp parallel
  offset = omp_get_wtime()
  !$omp end parallel
END SUBROUTINE init

!===============================================================================
!> Get offset
!===============================================================================
SUBROUTINE get (wc)
  real(8):: wc
  !.............................................................................
  wc = offset
END SUBROUTINE get

!===============================================================================
!> Set offset
!===============================================================================
SUBROUTINE set (wc)
  real(8):: wc
  !.............................................................................
  if (omp_timer%sync) then
    offset = wc
  end if
END SUBROUTINE set

END MODULE omp_timer_mod