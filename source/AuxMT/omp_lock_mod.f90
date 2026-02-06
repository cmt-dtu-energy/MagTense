MODULE omp_lock_mod
  USE omp_lib
  implicit none
  private

  !-----------------------------------------------------------------------------
  !> OpenMP lock wrapper
  !-----------------------------------------------------------------------------
  type, public :: omp_lock_t
    integer(omp_lock_kind) :: lk                 ! underlying OpenMP lock
    logical :: inited = .false.                  ! guard: init called?
  contains
    procedure :: init
    procedure :: destroy
    procedure :: lock
    procedure :: unlock
    procedure :: trylock
  end type omp_lock_t

CONTAINS

  !=============================================================================
  !> Initialise lock (idempotent)
  !=============================================================================
  subroutine init(self)
    class(omp_lock_t), intent(inout) :: self
    if (self%inited) return
    call omp_init_lock(self%lk)
    self%inited = .true.
  end subroutine init

  !=============================================================================
  !> Destroy lock (idempotent)
  !=============================================================================
  subroutine destroy(self)
    class(omp_lock_t), intent(inout) :: self
    if (.not. self%inited) return
    call omp_destroy_lock(self%lk)
    self%inited = .false.
  end subroutine destroy

  !=============================================================================
  !> Acquire lock (blocks until acquired)
  !=============================================================================
  subroutine lock(self)
    class(omp_lock_t), intent(inout) :: self
    if (.not. self%inited) call self%init()
    call omp_set_lock(self%lk)
  end subroutine lock

  !=============================================================================
  !> Release lock
  !=============================================================================
  subroutine unlock(self)
    class(omp_lock_t), intent(inout) :: self
    if (.not. self%inited) error stop "omp_lock_mod: unlock called on uninitialised lock."
    call omp_unset_lock(self%lk)
  end subroutine unlock

  !=============================================================================
  !> Try to acquire lock (non-blocking)
  !> Returns .true. if acquired, .false. otherwise.
  !=============================================================================
  logical function trylock(self)
    class(omp_lock_t), intent(inout) :: self
    if (.not. self%inited) call self%init()
    trylock = (omp_test_lock(self%lk) /= 0)
  end function trylock

END MODULE omp_lock_mod
