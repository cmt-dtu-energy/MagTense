MODULE omp_mod
  USE omp_lib
  implicit none
  private

  public :: omp_get_thread_num, omp_get_num_threads
  public :: omp_in_parallel

  type, public :: omp_t
    ! Per-thread
    integer :: nt     = 1
    integer :: tid    = 0
    logical :: master = .false.

    ! Tier 1
    integer :: max_threads       = -1
    integer :: num_procs         = -1
    logical :: dynamic_teams     = .false.
    integer :: max_active_levels = -1
    logical :: nested_parallel   = .false.   ! derived from max_active_levels
    integer :: proc_bind         = -1
    integer :: sched_kind        = -1
    integer :: sched_chunk       = -1
    integer :: openmp_version    = -1   ! from _OPENMP macro

    ! Tier 2 (places)
    integer :: place_num  = -1
    integer :: num_places = -1

  contains
    procedure, nopass :: init
    procedure, nopass :: info
    procedure, nopass :: in_parallel
    procedure, nopass :: is_master
    procedure, nopass :: numthreads
    procedure, nopass :: thread_id
  end type omp_t

  type(omp_t), public, save :: omp
  !$omp threadprivate(omp)

CONTAINS

  subroutine init()
    !DEC$ ATTRIBUTES ALIAS:"init_" :: init
    integer :: nthreads, tid_local
    integer :: k, c

#ifndef _OPENMP
    error stop "This code must be compiled with OpenMP enabled."
#endif
#if _OPENMP < 201811
    error stop "OpenMP 5.0+ required (_OPENMP >= 201811)."
#endif


    ! Set default schedule if needed - Should be called before any parallel region
    ! For more dynamic control use the "export OMP_SCHEDULE="static,4" " environment variable 
    !call set_schedule("static", 4)  ! default schedule

    !$omp parallel private(nthreads, tid_local, k, c) default(none)
      tid_local = omp_get_thread_num()
      nthreads  = omp_get_num_threads()

      omp%tid    = tid_local
      omp%nt     = nthreads
      omp%master = (tid_local == 0)

      ! Tier 2 per-thread placement (place)
      omp%place_num = omp_get_place_num()

      !$omp master
        omp%openmp_version    = _OPENMP
        omp%max_threads       = omp_get_max_threads()
        omp%num_procs         = omp_get_num_procs()
        omp%dynamic_teams     = (omp_get_dynamic() /= 0)

        omp%max_active_levels = omp_get_max_active_levels()
        omp%nested_parallel   = (omp%max_active_levels > 1)

        omp%proc_bind         = omp_get_proc_bind()

        call omp_get_schedule(k, c)
        omp%sched_kind  = k
        omp%sched_chunk = c

        omp%num_places = omp_get_num_places()
      !$omp end master
    !$omp end parallel
  end subroutine init

  subroutine info(unit)
    !DEC$ ATTRIBUTES ALIAS:"info_" :: info
    integer, intent(in), optional :: unit
    integer :: u
    u = 6
    if (present(unit)) u = unit

#ifndef _OPENMP
    error stop "This code must be compiled with OpenMP enabled."
#endif
#if _OPENMP < 201811
    error stop "OpenMP 5.0+ required (_OPENMP >= 201811)."
#endif

    if (omp_in_parallel()) then
      call info_parallel(u)
    else
      !$omp parallel default(none) shared(u)
        call info_parallel(u)
      !$omp end parallel
    end if
  end subroutine info

  subroutine info_parallel(u)
    integer, intent(in) :: u
    integer :: i, nt_local, tid_local
    integer :: mal

    tid_local = omp_get_thread_num()
    nt_local  = omp_get_num_threads()

    omp%tid    = tid_local
    omp%nt     = nt_local
    omp%master = (tid_local == 0)

    omp%place_num = omp_get_place_num()

    !$omp master
      mal = omp_get_max_active_levels()
      write(u,'(a)') "============================================================"
      write(u,'(a)') "OpenMP runtime info"
      write(u,'(a,i0)') "  _OPENMP (macro)        : ", _OPENMP
      write(u,'(a,a)')  "  OpenMP spec (mapped)   : ", trim(openmp_spec_string(_OPENMP))
      write(u,'(a,i0)') "  num_threads (team)     : ", nt_local
      write(u,'(a,i0)') "  max_threads            : ", omp_get_max_threads()
      write(u,'(a,i0)') "  num_procs              : ", omp_get_num_procs()
      write(u,'(a,l1)') "  dynamic teams          : ", (omp_get_dynamic() /= 0)
      write(u,'(a,i0)') "  max active levels      : ", mal
      write(u,'(a,l1)') "  nested parallel        : ", (mal > 1)
      write(u,'(a,a)')  "  proc_bind              : ", proc_bind_to_string(omp_get_proc_bind())
      call print_schedule(u)
      write(u,'(a,i0)') "  num_places             : ", omp_get_num_places()
      write(u,'(a)') "Per-thread placement:"
      write(u,'(a)') "------------------------------------------------------------"
    !$omp end master

    do i = 0, nt_local-1
      !$omp barrier
      if (tid_local == i) then
        write(u,'(a,i0,a,i0)') "  tid=", tid_local, "  place=", omp%place_num
      end if
      !$omp barrier
    end do

    !$omp master
      write(u,'(a)') "============================================================"
    !$omp end master
  end subroutine info_parallel

  LOGICAL FUNCTION in_parallel()
    in_parallel = omp_in_parallel()
  END FUNCTION in_parallel

  LOGICAL FUNCTION is_master()
    is_master = omp%master
  END FUNCTION is_master

  INTEGER FUNCTION numthreads()
    if (omp_in_parallel()) then
      numthreads = omp_get_num_threads()
    else
      numthreads = omp%nt
    end if
  END FUNCTION numthreads

  INTEGER FUNCTION thread_id()
    if (omp_in_parallel()) then
      thread_id = omp_get_thread_num()
    else
      thread_id = omp%tid
    end if
  END FUNCTION thread_id

  pure function openmp_spec_string(v) result(s)
    integer, intent(in) :: v
    character(len=16) :: s
    select case (v)
    case (201811); s = "5.0"
    case (202011); s = "5.1"
    case (202111); s = "5.2"
    case default
      write(s,'("unknown (",i0,")")') v
    end select
  end function openmp_spec_string

  pure function proc_bind_to_string(pb) result(s)
    integer, intent(in) :: pb
    character(len=16) :: s
    select case (pb)
    case (omp_proc_bind_false)  ; s = "false"
    case (omp_proc_bind_true)   ; s = "true"
    case (omp_proc_bind_master) ; s = "master"
    case (omp_proc_bind_close)  ; s = "close"
    case (omp_proc_bind_spread) ; s = "spread"
    case default                ; s = "unknown"
    end select
  end function proc_bind_to_string

  pure function schedule_to_string(k) result(s)
    integer, intent(in) :: k
    character(len=16) :: s
    select case (k)
    case (omp_sched_static)  ; s = "static, chunk="
    case (omp_sched_dynamic) ; s = "dynamic, chunk="
    case (omp_sched_guided)  ; s = "guided, chunk="
    case (omp_sched_auto)    ; s = "auto, chunk="
    case default             ; s = "unknown, chunk="
    end select
  end function schedule_to_string

  subroutine print_schedule(u)
    integer, intent(in) :: u
    integer :: k, c
    call omp_get_schedule(k, c)
    write(u,'(a,a,i0)') "  schedule               : ", schedule_to_string(k), c
  end subroutine print_schedule

  subroutine set_schedule(kind, chunk)
  character(len=*), intent(in) :: kind
  integer, intent(in), optional :: chunk
  integer :: c

  c = 0
  if (present(chunk)) c = chunk

  select case (kind)
  case ("static")
    call omp_set_schedule(omp_sched_static, c)
  case ("dynamic")
    call omp_set_schedule(omp_sched_dynamic, c)
  case ("guided")
    call omp_set_schedule(omp_sched_guided, c)
  case ("auto")
    call omp_set_schedule(omp_sched_auto, c)
  case default
    error stop "set_schedule: unknown kind (use static/dynamic/guided/auto)"
  end select
end subroutine set_schedule


END MODULE omp_mod
