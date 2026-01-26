module sort_mod
    use, intrinsic :: iso_fortran_env, only: int32, real32, real64
    implicit none
    private

    !===============================================================================================
    ! Public API
    !===============================================================================================
    public :: sort
    public :: sort_r
    public :: argsort
    public :: argsort_r
    public :: apply_perm

    !-----------------------------------------------------------------------------------------------
    ! Module parameters
    !-----------------------------------------------------------------------------------------------
    integer, parameter :: INSERTION_CUTOFF_DEFAULT = 24
    !-----------------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------------
    ! Generic interfaces
    !-----------------------------------------------------------------------------------------------
    interface sort
        module procedure sort_inplace_i4
        module procedure sort_inplace_f4
        module procedure sort_inplace_f8
        module procedure sort_copy_i4
        module procedure sort_copy_f4
        module procedure sort_copy_f8
    end interface sort

    interface sort_r
        module procedure sort_inplace_i4_r
        module procedure sort_inplace_f4_r
        module procedure sort_inplace_f8_r
        module procedure sort_copy_i4_r
        module procedure sort_copy_f4_r
        module procedure sort_copy_f8_r
    end interface sort_r

    interface argsort
        module procedure argsort_i4
        module procedure argsort_f4
        module procedure argsort_f8
    end interface argsort

    interface argsort_r
        module procedure argsort_i4_r
        module procedure argsort_f4_r
        module procedure argsort_f8_r
    end interface argsort_r

    interface apply_perm
        module procedure apply_perm_i4
        module procedure apply_perm_f4
        module procedure apply_perm_f8
    end interface apply_perm

contains

    !===============================================================================================
    !> Sort values in-place (ascending).
    !!
    !!  Arguments:
    !!    values  [inout] Array to be sorted in-place.
    !!    algo    [in]    Optional sorting algorithm selector (default "quicksort").
    !!    cutoff  [in]    Optional insertion-sort cutoff for small partitions.
    !===============================================================================================
    subroutine sort_inplace_i4(values, algo, cutoff)
        implicit none
        integer(int32), intent(inout)           :: values(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n
        integer :: insertion_cutoff

        !-------------------------------------------------------------------------------------------
        ! Validate algo
        !-------------------------------------------------------------------------------------------
        call require_quicksort(algo, "sort_inplace_i4")
        !-------------------------------------------------------------------------------------------

        !-------------------------------------------------------------------------------------------
        ! Early return
        !-------------------------------------------------------------------------------------------
        n = size(values)
        if (n <= 1) return
        !-------------------------------------------------------------------------------------------

        !-------------------------------------------------------------------------------------------
        ! Cutoff
        !-------------------------------------------------------------------------------------------
        insertion_cutoff = INSERTION_CUTOFF_DEFAULT
        if (present(cutoff)) insertion_cutoff = cutoff
        !-------------------------------------------------------------------------------------------

        call qsort_values_i4(values, 1, n, insertion_cutoff)
    end subroutine sort_inplace_i4


    subroutine sort_inplace_f4(values, algo, cutoff)
        implicit none
        real(real32), intent(inout)             :: values(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n
        integer :: insertion_cutoff

        call require_quicksort(algo, "sort_inplace_f4")

        n = size(values)
        if (n <= 1) return

        insertion_cutoff = INSERTION_CUTOFF_DEFAULT
        if (present(cutoff)) insertion_cutoff = cutoff

        call qsort_values_f4(values, 1, n, insertion_cutoff)
    end subroutine sort_inplace_f4


    subroutine sort_inplace_f8(values, algo, cutoff)
        implicit none
        real(real64), intent(inout)             :: values(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n
        integer :: insertion_cutoff

        call require_quicksort(algo, "sort_inplace_f8")

        n = size(values)
        if (n <= 1) return

        insertion_cutoff = INSERTION_CUTOFF_DEFAULT
        if (present(cutoff)) insertion_cutoff = cutoff

        call qsort_values_f8(values, 1, n, insertion_cutoff)
    end subroutine sort_inplace_f8


    !===============================================================================================
    !> Sort values in-place (descending).
    !!
    !!  Arguments:
    !!    values  [inout] Array to be sorted in-place.
    !!    algo    [in]    Optional sorting algorithm selector (default "quicksort").
    !!    cutoff  [in]    Optional insertion-sort cutoff for small partitions.
    !===============================================================================================
    subroutine sort_inplace_i4_r(values, algo, cutoff)
        implicit none
        integer(int32), intent(inout)           :: values(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n
        integer :: insertion_cutoff

        call require_quicksort(algo, "sort_inplace_i4_r")

        n = size(values)
        if (n <= 1) return

        insertion_cutoff = INSERTION_CUTOFF_DEFAULT
        if (present(cutoff)) insertion_cutoff = cutoff

        call qsort_values_i4_r(values, 1, n, insertion_cutoff)
    end subroutine sort_inplace_i4_r


    subroutine sort_inplace_f4_r(values, algo, cutoff)
        implicit none
        real(real32), intent(inout)             :: values(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n
        integer :: insertion_cutoff

        call require_quicksort(algo, "sort_inplace_f4_r")

        n = size(values)
        if (n <= 1) return

        insertion_cutoff = INSERTION_CUTOFF_DEFAULT
        if (present(cutoff)) insertion_cutoff = cutoff

        call qsort_values_f4_r(values, 1, n, insertion_cutoff)
    end subroutine sort_inplace_f4_r


    subroutine sort_inplace_f8_r(values, algo, cutoff)
        implicit none
        real(real64), intent(inout)             :: values(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n
        integer :: insertion_cutoff

        call require_quicksort(algo, "sort_inplace_f8_r")

        n = size(values)
        if (n <= 1) return

        insertion_cutoff = INSERTION_CUTOFF_DEFAULT
        if (present(cutoff)) insertion_cutoff = cutoff

        call qsort_values_f8_r(values, 1, n, insertion_cutoff)
    end subroutine sort_inplace_f8_r


    !===============================================================================================
    !> Sort input_values into sorted_values (ascending).
    !!
    !!  Arguments:
    !!    input_values   [in]   Input array (not modified).
    !!    sorted_values  [out]  Sorted output (same size as input_values).
    !!    algo           [in]   Optional sorting algorithm selector (default "quicksort").
    !!    cutoff         [in]   Optional insertion-sort cutoff for small partitions.
    !===============================================================================================
    subroutine sort_copy_i4(input_values, sorted_values, algo, cutoff)
        implicit none
        integer(int32), intent(in)              :: input_values(:)
        integer(int32), intent(out)             :: sorted_values(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n

        call require_quicksort(algo, "sort_copy_i4")

        n = size(input_values)
        if (size(sorted_values) /= n) error stop "sort_copy_i4: size mismatch"

        sorted_values = input_values
        call sort_inplace_i4(sorted_values, algo="quicksort", cutoff=cutoff)
    end subroutine sort_copy_i4


    subroutine sort_copy_f4(input_values, sorted_values, algo, cutoff)
        implicit none
        real(real32), intent(in)                :: input_values(:)
        real(real32), intent(out)               :: sorted_values(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n

        call require_quicksort(algo, "sort_copy_f4")

        n = size(input_values)
        if (size(sorted_values) /= n) error stop "sort_copy_f4: size mismatch"

        sorted_values = input_values
        call sort_inplace_f4(sorted_values, algo="quicksort", cutoff=cutoff)
    end subroutine sort_copy_f4


    subroutine sort_copy_f8(input_values, sorted_values, algo, cutoff)
        implicit none
        real(real64), intent(in)                :: input_values(:)
        real(real64), intent(out)               :: sorted_values(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n

        call require_quicksort(algo, "sort_copy_f8")

        n = size(input_values)
        if (size(sorted_values) /= n) error stop "sort_copy_f8: size mismatch"

        sorted_values = input_values
        call sort_inplace_f8(sorted_values, algo="quicksort", cutoff=cutoff)
    end subroutine sort_copy_f8


    !===============================================================================================
    !> Sort input_values into sorted_values (descending).
    !!
    !!  Arguments:
    !!    input_values   [in]   Input array (not modified).
    !!    sorted_values  [out]  Sorted output (same size as input_values).
    !!    algo           [in]   Optional sorting algorithm selector (default "quicksort").
    !!    cutoff         [in]   Optional insertion-sort cutoff for small partitions.
    !===============================================================================================
    subroutine sort_copy_i4_r(input_values, sorted_values, algo, cutoff)
        implicit none
        integer(int32), intent(in)              :: input_values(:)
        integer(int32), intent(out)             :: sorted_values(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n

        call require_quicksort(algo, "sort_copy_i4_r")

        n = size(input_values)
        if (size(sorted_values) /= n) error stop "sort_copy_i4_r: size mismatch"

        sorted_values = input_values
        call sort_inplace_i4_r(sorted_values, algo="quicksort", cutoff=cutoff)
    end subroutine sort_copy_i4_r


    subroutine sort_copy_f4_r(input_values, sorted_values, algo, cutoff)
        implicit none
        real(real32), intent(in)                :: input_values(:)
        real(real32), intent(out)               :: sorted_values(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n

        call require_quicksort(algo, "sort_copy_f4_r")

        n = size(input_values)
        if (size(sorted_values) /= n) error stop "sort_copy_f4_r: size mismatch"

        sorted_values = input_values
        call sort_inplace_f4_r(sorted_values, algo="quicksort", cutoff=cutoff)
    end subroutine sort_copy_f4_r


    subroutine sort_copy_f8_r(input_values, sorted_values, algo, cutoff)
        implicit none
        real(real64), intent(in)                :: input_values(:)
        real(real64), intent(out)               :: sorted_values(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n

        call require_quicksort(algo, "sort_copy_f8_r")

        n = size(input_values)
        if (size(sorted_values) /= n) error stop "sort_copy_f8_r: size mismatch"

        sorted_values = input_values
        call sort_inplace_f8_r(sorted_values, algo="quicksort", cutoff=cutoff)
    end subroutine sort_copy_f8_r


    !===============================================================================================
    !> Argsort: sorted_indices such that keys(sorted_indices) is ascending.
    !!
    !!  Arguments:
    !!    keys           [in]   Values to be sorted (not modified).
    !!    sorted_indices [out]  Permutation indices (1-based), size(sorted_indices) == size(keys).
    !!    algo           [in]   Optional sorting algorithm selector (default "quicksort").
    !!    cutoff         [in]   Optional insertion-sort cutoff for small partitions.
    !===============================================================================================
    subroutine argsort_i4(keys, sorted_indices, algo, cutoff)
        implicit none
        integer(int32), intent(in)              :: keys(:)
        integer,        intent(out)             :: sorted_indices(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n, i
        integer :: insertion_cutoff

        call require_quicksort(algo, "argsort_i4")

        n = size(keys)
        if (size(sorted_indices) /= n) error stop "argsort_i4: size mismatch"

        if (n <= 1) then
            if (n == 1) sorted_indices(1) = 1
            return
        end if

        do i = 1, n
            sorted_indices(i) = i
        end do

        insertion_cutoff = INSERTION_CUTOFF_DEFAULT
        if (present(cutoff)) insertion_cutoff = cutoff

        call qsort_indices_i4(keys, sorted_indices, 1, n, insertion_cutoff)
    end subroutine argsort_i4


    subroutine argsort_f4(keys, sorted_indices, algo, cutoff)
        implicit none
        real(real32), intent(in)                :: keys(:)
        integer,      intent(out)               :: sorted_indices(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n, i
        integer :: insertion_cutoff

        call require_quicksort(algo, "argsort_f4")

        n = size(keys)
        if (size(sorted_indices) /= n) error stop "argsort_f4: size mismatch"

        if (n <= 1) then
            if (n == 1) sorted_indices(1) = 1
            return
        end if

        do i = 1, n
            sorted_indices(i) = i
        end do

        insertion_cutoff = INSERTION_CUTOFF_DEFAULT
        if (present(cutoff)) insertion_cutoff = cutoff

        call qsort_indices_f4(keys, sorted_indices, 1, n, insertion_cutoff)
    end subroutine argsort_f4


    subroutine argsort_f8(keys, sorted_indices, algo, cutoff)
        implicit none
        real(real64), intent(in)                :: keys(:)
        integer,      intent(out)               :: sorted_indices(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n, i
        integer :: insertion_cutoff

        call require_quicksort(algo, "argsort_f8")

        n = size(keys)
        if (size(sorted_indices) /= n) error stop "argsort_f8: size mismatch"

        if (n <= 1) then
            if (n == 1) sorted_indices(1) = 1
            return
        end if

        do i = 1, n
            sorted_indices(i) = i
        end do

        insertion_cutoff = INSERTION_CUTOFF_DEFAULT
        if (present(cutoff)) insertion_cutoff = cutoff

        call qsort_indices_f8(keys, sorted_indices, 1, n, insertion_cutoff)
    end subroutine argsort_f8


    !===============================================================================================
    !> Argsort: sorted_indices such that keys(sorted_indices) is descending.
    !!
    !!  Arguments:
    !!    keys           [in]   Values to be sorted (not modified).
    !!    sorted_indices [out]  Permutation indices (1-based), size(sorted_indices) == size(keys).
    !!    algo           [in]   Optional sorting algorithm selector (default "quicksort").
    !!    cutoff         [in]   Optional insertion-sort cutoff for small partitions.
    !===============================================================================================
    subroutine argsort_i4_r(keys, sorted_indices, algo, cutoff)
        implicit none
        integer(int32), intent(in)              :: keys(:)
        integer,        intent(out)             :: sorted_indices(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n, i
        integer :: insertion_cutoff

        call require_quicksort(algo, "argsort_i4_r")

        n = size(keys)
        if (size(sorted_indices) /= n) error stop "argsort_i4_r: size mismatch"

        if (n <= 1) then
            if (n == 1) sorted_indices(1) = 1
            return
        end if

        do i = 1, n
            sorted_indices(i) = i
        end do

        insertion_cutoff = INSERTION_CUTOFF_DEFAULT
        if (present(cutoff)) insertion_cutoff = cutoff

        call qsort_indices_i4_r(keys, sorted_indices, 1, n, insertion_cutoff)
    end subroutine argsort_i4_r


    subroutine argsort_f4_r(keys, sorted_indices, algo, cutoff)
        implicit none
        real(real32), intent(in)                :: keys(:)
        integer,      intent(out)               :: sorted_indices(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n, i
        integer :: insertion_cutoff

        call require_quicksort(algo, "argsort_f4_r")

        n = size(keys)
        if (size(sorted_indices) /= n) error stop "argsort_f4_r: size mismatch"

        if (n <= 1) then
            if (n == 1) sorted_indices(1) = 1
            return
        end if

        do i = 1, n
            sorted_indices(i) = i
        end do

        insertion_cutoff = INSERTION_CUTOFF_DEFAULT
        if (present(cutoff)) insertion_cutoff = cutoff

        call qsort_indices_f4_r(keys, sorted_indices, 1, n, insertion_cutoff)
    end subroutine argsort_f4_r


    subroutine argsort_f8_r(keys, sorted_indices, algo, cutoff)
        implicit none
        real(real64), intent(in)                :: keys(:)
        integer,      intent(out)               :: sorted_indices(:)
        character(len=*), intent(in), optional  :: algo
        integer,          intent(in), optional  :: cutoff

        integer :: n, i
        integer :: insertion_cutoff

        call require_quicksort(algo, "argsort_f8_r")

        n = size(keys)
        if (size(sorted_indices) /= n) error stop "argsort_f8_r: size mismatch"

        if (n <= 1) then
            if (n == 1) sorted_indices(1) = 1
            return
        end if

        do i = 1, n
            sorted_indices(i) = i
        end do

        insertion_cutoff = INSERTION_CUTOFF_DEFAULT
        if (present(cutoff)) insertion_cutoff = cutoff

        call qsort_indices_f8_r(keys, sorted_indices, 1, n, insertion_cutoff)
    end subroutine argsort_f8_r


    !===============================================================================================
    !> Apply a permutation: output_values(i) = input_values(sorted_indices(i)).
    !!
    !!  Arguments:
    !!    input_values   [in]   Input array.
    !!    sorted_indices [in]   Permutation indices (1-based).
    !!    output_values  [out]  Permuted output array.
    !===============================================================================================
    subroutine apply_perm_i4(input_values, sorted_indices, output_values)
        implicit none
        integer(int32), intent(in)  :: input_values(:)
        integer,        intent(in)  :: sorted_indices(:)
        integer(int32), intent(out) :: output_values(:)

        integer :: n
        integer :: i

        n = size(sorted_indices)
        if (size(input_values)  /= n) error stop "apply_perm_i4: size mismatch"
        if (size(output_values) /= n) error stop "apply_perm_i4: size mismatch"

        do i = 1, n
            output_values(i) = input_values(sorted_indices(i))
        end do
    end subroutine apply_perm_i4


    subroutine apply_perm_f4(input_values, sorted_indices, output_values)
        implicit none
        real(real32), intent(in)  :: input_values(:)
        integer,      intent(in)  :: sorted_indices(:)
        real(real32), intent(out) :: output_values(:)

        integer :: n
        integer :: i

        n = size(sorted_indices)
        if (size(input_values)  /= n) error stop "apply_perm_f4: size mismatch"
        if (size(output_values) /= n) error stop "apply_perm_f4: size mismatch"

        do i = 1, n
            output_values(i) = input_values(sorted_indices(i))
        end do
    end subroutine apply_perm_f4


    subroutine apply_perm_f8(input_values, sorted_indices, output_values)
        implicit none
        real(real64), intent(in)  :: input_values(:)
        integer,      intent(in)  :: sorted_indices(:)
        real(real64), intent(out) :: output_values(:)

        integer :: n
        integer :: i

        n = size(sorted_indices)
        if (size(input_values)  /= n) error stop "apply_perm_f8: size mismatch"
        if (size(output_values) /= n) error stop "apply_perm_f8: size mismatch"

        do i = 1, n
            output_values(i) = input_values(sorted_indices(i))
        end do
    end subroutine apply_perm_f8


    !-----------------------------------------------------------------------------------------------
    ! Algo guard: only "quicksort" supported for now.
    !-----------------------------------------------------------------------------------------------
    subroutine require_quicksort(algo, caller_name)
        implicit none
        character(len=*), intent(in), optional :: algo
        character(len=*), intent(in)           :: caller_name

        character(len=:), allocatable :: algo_local

        if (.not. present(algo)) return

        algo_local = adjustl(algo)

        if (trim(algo_local) /= "quicksort") then
            error stop caller_name // ": unsupported algo='" // trim(algo_local) // "' (only 'quicksort')"
        end if
    end subroutine require_quicksort


    !===============================================================================================
    ! Internal sorting kernels: values (ascending / descending)
    !===============================================================================================

    recursive subroutine qsort_values_i4(values, lo, hi, cutoff)
        implicit none
        integer(int32), intent(inout) :: values(:)
        integer,        intent(in)    :: lo, hi, cutoff

        integer :: i, j, mid
        integer(int32) :: pivot, tmp

        if (hi - lo + 1 <= cutoff) then
            call insertion_values_i4(values, lo, hi)
            return
        end if

        mid   = (lo + hi) / 2
        pivot = values(mid)

        i = lo
        j = hi
        do
            do while (values(i) < pivot); i = i + 1; end do
            do while (values(j) > pivot); j = j - 1; end do

            if (i <= j) then
                tmp = values(i); values(i) = values(j); values(j) = tmp
                i = i + 1
                j = j - 1
            end if

            if (i > j) exit
        end do

        if (lo < j) call qsort_values_i4(values, lo, j, cutoff)
        if (i  < hi) call qsort_values_i4(values, i,  hi, cutoff)
    end subroutine qsort_values_i4


    subroutine insertion_values_i4(values, lo, hi)
        implicit none
        integer(int32), intent(inout) :: values(:)
        integer,        intent(in)    :: lo, hi

        integer :: i, j
        integer(int32) :: held

        do i = lo + 1, hi
            held = values(i)
            j = i - 1
            do while (j >= lo .and. values(j) > held)
                values(j+1) = values(j)
                j = j - 1
            end do
            values(j+1) = held
        end do
    end subroutine insertion_values_i4


    recursive subroutine qsort_values_i4_r(values, lo, hi, cutoff)
        implicit none
        integer(int32), intent(inout) :: values(:)
        integer,        intent(in)    :: lo, hi, cutoff

        integer :: i, j, mid
        integer(int32) :: pivot, tmp

        if (hi - lo + 1 <= cutoff) then
            call insertion_values_i4_r(values, lo, hi)
            return
        end if

        mid   = (lo + hi) / 2
        pivot = values(mid)

        i = lo
        j = hi
        do
            do while (values(i) > pivot); i = i + 1; end do
            do while (values(j) < pivot); j = j - 1; end do

            if (i <= j) then
                tmp = values(i); values(i) = values(j); values(j) = tmp
                i = i + 1
                j = j - 1
            end if

            if (i > j) exit
        end do

        if (lo < j) call qsort_values_i4_r(values, lo, j, cutoff)
        if (i  < hi) call qsort_values_i4_r(values, i,  hi, cutoff)
    end subroutine qsort_values_i4_r


    subroutine insertion_values_i4_r(values, lo, hi)
        implicit none
        integer(int32), intent(inout) :: values(:)
        integer,        intent(in)    :: lo, hi

        integer :: i, j
        integer(int32) :: held

        do i = lo + 1, hi
            held = values(i)
            j = i - 1
            do while (j >= lo .and. values(j) < held)
                values(j+1) = values(j)
                j = j - 1
            end do
            values(j+1) = held
        end do
    end subroutine insertion_values_i4_r


    recursive subroutine qsort_values_f4(values, lo, hi, cutoff)
        implicit none
        real(real32), intent(inout) :: values(:)
        integer,      intent(in)    :: lo, hi, cutoff

        integer :: i, j, mid
        real(real32) :: pivot, tmp

        if (hi - lo + 1 <= cutoff) then
            call insertion_values_f4(values, lo, hi)
            return
        end if

        mid   = (lo + hi) / 2
        pivot = values(mid)

        i = lo
        j = hi
        do
            do while (values(i) < pivot); i = i + 1; end do
            do while (values(j) > pivot); j = j - 1; end do

            if (i <= j) then
                tmp = values(i); values(i) = values(j); values(j) = tmp
                i = i + 1
                j = j - 1
            end if

            if (i > j) exit
        end do

        if (lo < j) call qsort_values_f4(values, lo, j, cutoff)
        if (i  < hi) call qsort_values_f4(values, i,  hi, cutoff)
    end subroutine qsort_values_f4


    subroutine insertion_values_f4(values, lo, hi)
        implicit none
        real(real32), intent(inout) :: values(:)
        integer,      intent(in)    :: lo, hi

        integer :: i, j
        real(real32) :: held

        do i = lo + 1, hi
            held = values(i)
            j = i - 1
            do while (j >= lo .and. values(j) > held)
                values(j+1) = values(j)
                j = j - 1
            end do
            values(j+1) = held
        end do
    end subroutine insertion_values_f4


    recursive subroutine qsort_values_f4_r(values, lo, hi, cutoff)
        implicit none
        real(real32), intent(inout) :: values(:)
        integer,      intent(in)    :: lo, hi, cutoff

        integer :: i, j, mid
        real(real32) :: pivot, tmp

        if (hi - lo + 1 <= cutoff) then
            call insertion_values_f4_r(values, lo, hi)
            return
        end if

        mid   = (lo + hi) / 2
        pivot = values(mid)

        i = lo
        j = hi
        do
            do while (values(i) > pivot); i = i + 1; end do
            do while (values(j) < pivot); j = j - 1; end do

            if (i <= j) then
                tmp = values(i); values(i) = values(j); values(j) = tmp
                i = i + 1
                j = j - 1
            end if

            if (i > j) exit
        end do

        if (lo < j) call qsort_values_f4_r(values, lo, j, cutoff)
        if (i  < hi) call qsort_values_f4_r(values, i,  hi, cutoff)
    end subroutine qsort_values_f4_r


    subroutine insertion_values_f4_r(values, lo, hi)
        implicit none
        real(real32), intent(inout) :: values(:)
        integer,      intent(in)    :: lo, hi

        integer :: i, j
        real(real32) :: held

        do i = lo + 1, hi
            held = values(i)
            j = i - 1
            do while (j >= lo .and. values(j) < held)
                values(j+1) = values(j)
                j = j - 1
            end do
            values(j+1) = held
        end do
    end subroutine insertion_values_f4_r


    recursive subroutine qsort_values_f8(values, lo, hi, cutoff)
        implicit none
        real(real64), intent(inout) :: values(:)
        integer,      intent(in)    :: lo, hi, cutoff

        integer :: i, j, mid
        real(real64) :: pivot, tmp

        if (hi - lo + 1 <= cutoff) then
            call insertion_values_f8(values, lo, hi)
            return
        end if

        mid   = (lo + hi) / 2
        pivot = values(mid)

        i = lo
        j = hi
        do
            do while (values(i) < pivot); i = i + 1; end do
            do while (values(j) > pivot); j = j - 1; end do

            if (i <= j) then
                tmp = values(i); values(i) = values(j); values(j) = tmp
                i = i + 1
                j = j - 1
            end if

            if (i > j) exit
        end do

        if (lo < j) call qsort_values_f8(values, lo, j, cutoff)
        if (i  < hi) call qsort_values_f8(values, i,  hi, cutoff)
    end subroutine qsort_values_f8


    subroutine insertion_values_f8(values, lo, hi)
        implicit none
        real(real64), intent(inout) :: values(:)
        integer,      intent(in)    :: lo, hi

        integer :: i, j
        real(real64) :: held

        do i = lo + 1, hi
            held = values(i)
            j = i - 1
            do while (j >= lo .and. values(j) > held)
                values(j+1) = values(j)
                j = j - 1
            end do
            values(j+1) = held
        end do
    end subroutine insertion_values_f8


    recursive subroutine qsort_values_f8_r(values, lo, hi, cutoff)
        implicit none
        real(real64), intent(inout) :: values(:)
        integer,      intent(in)    :: lo, hi, cutoff

        integer :: i, j, mid
        real(real64) :: pivot, tmp

        if (hi - lo + 1 <= cutoff) then
            call insertion_values_f8_r(values, lo, hi)
            return
        end if

        mid   = (lo + hi) / 2
        pivot = values(mid)

        i = lo
        j = hi
        do
            do while (values(i) > pivot); i = i + 1; end do
            do while (values(j) < pivot); j = j - 1; end do

            if (i <= j) then
                tmp = values(i); values(i) = values(j); values(j) = tmp
                i = i + 1
                j = j - 1
            end if

            if (i > j) exit
        end do

        if (lo < j) call qsort_values_f8_r(values, lo, j, cutoff)
        if (i  < hi) call qsort_values_f8_r(values, i,  hi, cutoff)
    end subroutine qsort_values_f8_r


    subroutine insertion_values_f8_r(values, lo, hi)
        implicit none
        real(real64), intent(inout) :: values(:)
        integer,      intent(in)    :: lo, hi

        integer :: i, j
        real(real64) :: held

        do i = lo + 1, hi
            held = values(i)
            j = i - 1
            do while (j >= lo .and. values(j) < held)
                values(j+1) = values(j)
                j = j - 1
            end do
            values(j+1) = held
        end do
    end subroutine insertion_values_f8_r


    !===============================================================================================
    ! Internal sorting kernels: argsort indices (ascending / descending)
    !===============================================================================================

    recursive subroutine qsort_indices_i4(keys, sorted_indices, lo, hi, cutoff)
        implicit none
        integer(int32), intent(in)    :: keys(:)
        integer,        intent(inout) :: sorted_indices(:)
        integer,        intent(in)    :: lo, hi, cutoff

        integer :: i, j, mid, tmp_index
        integer(int32) :: pivot_key

        if (hi - lo + 1 <= cutoff) then
            call insertion_indices_i4(keys, sorted_indices, lo, hi)
            return
        end if

        mid       = (lo + hi) / 2
        pivot_key = keys(sorted_indices(mid))

        i = lo
        j = hi
        do
            do while (keys(sorted_indices(i)) < pivot_key); i = i + 1; end do
            do while (keys(sorted_indices(j)) > pivot_key); j = j - 1; end do

            if (i <= j) then
                tmp_index = sorted_indices(i); sorted_indices(i) = sorted_indices(j); sorted_indices(j) = tmp_index
                i = i + 1
                j = j - 1
            end if

            if (i > j) exit
        end do

        if (lo < j) call qsort_indices_i4(keys, sorted_indices, lo, j, cutoff)
        if (i  < hi) call qsort_indices_i4(keys, sorted_indices, i,  hi, cutoff)
    end subroutine qsort_indices_i4


    subroutine insertion_indices_i4(keys, sorted_indices, lo, hi)
        use, intrinsic :: iso_fortran_env, only: int32
        implicit none

        integer(int32), intent(in)    :: keys(:)
        integer,        intent(inout) :: sorted_indices(:)
        integer,        intent(in)    :: lo, hi

        integer :: i
        integer :: j
        integer :: held_index

        do i = lo + 1, hi
            held_index = sorted_indices(i)
            j = i - 1

            do while (j >= lo)
                if (keys(sorted_indices(j)) <= keys(held_index)) exit
                sorted_indices(j+1) = sorted_indices(j)
                j = j - 1
            end do

            sorted_indices(j+1) = held_index
        end do
    end subroutine insertion_indices_i4



    recursive subroutine qsort_indices_i4_r(keys, sorted_indices, lo, hi, cutoff)
        implicit none
        integer(int32), intent(in)    :: keys(:)
        integer,        intent(inout) :: sorted_indices(:)
        integer,        intent(in)    :: lo, hi, cutoff

        integer :: i, j, mid, tmp_index
        integer(int32) :: pivot_key

        if (hi - lo + 1 <= cutoff) then
            call insertion_indices_i4_r(keys, sorted_indices, lo, hi)
            return
        end if

        mid       = (lo + hi) / 2
        pivot_key = keys(sorted_indices(mid))

        i = lo
        j = hi
        do
            do while (keys(sorted_indices(i)) > pivot_key); i = i + 1; end do
            do while (keys(sorted_indices(j)) < pivot_key); j = j - 1; end do

            if (i <= j) then
                tmp_index = sorted_indices(i); sorted_indices(i) = sorted_indices(j); sorted_indices(j) = tmp_index
                i = i + 1
                j = j - 1
            end if

            if (i > j) exit
        end do

        if (lo < j) call qsort_indices_i4_r(keys, sorted_indices, lo, j, cutoff)
        if (i  < hi) call qsort_indices_i4_r(keys, sorted_indices, i,  hi, cutoff)
    end subroutine qsort_indices_i4_r


    subroutine insertion_indices_i4_r(keys, sorted_indices, lo, hi)
        use, intrinsic :: iso_fortran_env, only: int32
        implicit none

        integer(int32), intent(in)    :: keys(:)
        integer,        intent(inout) :: sorted_indices(:)
        integer,        intent(in)    :: lo, hi

        integer :: i
        integer :: j
        integer :: held_index

        do i = lo + 1, hi
            held_index = sorted_indices(i)
            j = i - 1

            do while (j >= lo)
                if (keys(sorted_indices(j)) >= keys(held_index)) exit
                sorted_indices(j+1) = sorted_indices(j)
                j = j - 1
            end do

            sorted_indices(j+1) = held_index
        end do
    end subroutine insertion_indices_i4_r


    recursive subroutine qsort_indices_f4(keys, sorted_indices, lo, hi, cutoff)
        implicit none
        real(real32), intent(in)      :: keys(:)
        integer,      intent(inout)   :: sorted_indices(:)
        integer,      intent(in)      :: lo, hi, cutoff

        integer :: i, j, mid, tmp_index
        real(real32) :: pivot_key

        if (hi - lo + 1 <= cutoff) then
            call insertion_indices_f4(keys, sorted_indices, lo, hi)
            return
        end if

        mid       = (lo + hi) / 2
        pivot_key = keys(sorted_indices(mid))

        i = lo
        j = hi
        do
            do while (keys(sorted_indices(i)) < pivot_key); i = i + 1; end do
            do while (keys(sorted_indices(j)) > pivot_key); j = j - 1; end do

            if (i <= j) then
                tmp_index = sorted_indices(i); sorted_indices(i) = sorted_indices(j); sorted_indices(j) = tmp_index
                i = i + 1
                j = j - 1
            end if

            if (i > j) exit
        end do

        if (lo < j) call qsort_indices_f4(keys, sorted_indices, lo, j, cutoff)
        if (i  < hi) call qsort_indices_f4(keys, sorted_indices, i,  hi, cutoff)
    end subroutine qsort_indices_f4


    subroutine insertion_indices_f4(keys, sorted_indices, lo, hi)
        use, intrinsic :: iso_fortran_env, only: real32
        implicit none

        real(real32), intent(in)      :: keys(:)
        integer,      intent(inout)   :: sorted_indices(:)
        integer,      intent(in)      :: lo, hi

        integer :: i
        integer :: j
        integer :: held_index

        do i = lo + 1, hi
            held_index = sorted_indices(i)
            j = i - 1

            do while (j >= lo)
                if (keys(sorted_indices(j)) <= keys(held_index)) exit
                sorted_indices(j+1) = sorted_indices(j)
                j = j - 1
            end do

            sorted_indices(j+1) = held_index
        end do
    end subroutine insertion_indices_f4



    recursive subroutine qsort_indices_f4_r(keys, sorted_indices, lo, hi, cutoff)
        implicit none
        real(real32), intent(in)      :: keys(:)
        integer,      intent(inout)   :: sorted_indices(:)
        integer,      intent(in)      :: lo, hi, cutoff

        integer :: i, j, mid, tmp_index
        real(real32) :: pivot_key

        if (hi - lo + 1 <= cutoff) then
            call insertion_indices_f4_r(keys, sorted_indices, lo, hi)
            return
        end if

        mid       = (lo + hi) / 2
        pivot_key = keys(sorted_indices(mid))

        i = lo
        j = hi
        do
            do while (keys(sorted_indices(i)) > pivot_key); i = i + 1; end do
            do while (keys(sorted_indices(j)) < pivot_key); j = j - 1; end do

            if (i <= j) then
                tmp_index = sorted_indices(i); sorted_indices(i) = sorted_indices(j); sorted_indices(j) = tmp_index
                i = i + 1
                j = j - 1
            end if

            if (i > j) exit
        end do

        if (lo < j) call qsort_indices_f4_r(keys, sorted_indices, lo, j, cutoff)
        if (i  < hi) call qsort_indices_f4_r(keys, sorted_indices, i,  hi, cutoff)
    end subroutine qsort_indices_f4_r


    subroutine insertion_indices_f4_r(keys, sorted_indices, lo, hi)
        use, intrinsic :: iso_fortran_env, only: real32
        implicit none

        real(real32), intent(in)      :: keys(:)
        integer,      intent(inout)   :: sorted_indices(:)
        integer,      intent(in)      :: lo, hi

        integer :: i
        integer :: j
        integer :: held_index

        do i = lo + 1, hi
            held_index = sorted_indices(i)
            j = i - 1

            do while (j >= lo)
                if (keys(sorted_indices(j)) >= keys(held_index)) exit
                sorted_indices(j+1) = sorted_indices(j)
                j = j - 1
            end do

            sorted_indices(j+1) = held_index
        end do
    end subroutine insertion_indices_f4_r



    recursive subroutine qsort_indices_f8(keys, sorted_indices, lo, hi, cutoff)
        implicit none
        real(real64), intent(in)      :: keys(:)
        integer,      intent(inout)   :: sorted_indices(:)
        integer,      intent(in)      :: lo, hi, cutoff

        integer :: i, j, mid, tmp_index
        real(real64) :: pivot_key

        if (hi - lo + 1 <= cutoff) then
            call insertion_indices_f8(keys, sorted_indices, lo, hi)
            return
        end if

        mid       = (lo + hi) / 2
        pivot_key = keys(sorted_indices(mid))

        i = lo
        j = hi
        do
            do while (keys(sorted_indices(i)) < pivot_key); i = i + 1; end do
            do while (keys(sorted_indices(j)) > pivot_key); j = j - 1; end do

            if (i <= j) then
                tmp_index = sorted_indices(i); sorted_indices(i) = sorted_indices(j); sorted_indices(j) = tmp_index
                i = i + 1
                j = j - 1
            end if

            if (i > j) exit
        end do

        if (lo < j) call qsort_indices_f8(keys, sorted_indices, lo, j, cutoff)
        if (i  < hi) call qsort_indices_f8(keys, sorted_indices, i,  hi, cutoff)
    end subroutine qsort_indices_f8


    subroutine insertion_indices_f8(keys, sorted_indices, lo, hi)
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none

        real(real64), intent(in)      :: keys(:)
        integer,      intent(inout)   :: sorted_indices(:)
        integer,      intent(in)      :: lo, hi

        integer :: i
        integer :: j
        integer :: held_index

        do i = lo + 1, hi
            held_index = sorted_indices(i)
            j = i - 1

            do while (j >= lo)
                if (keys(sorted_indices(j)) <= keys(held_index)) exit
                sorted_indices(j+1) = sorted_indices(j)
                j = j - 1
            end do

            sorted_indices(j+1) = held_index
        end do
    end subroutine insertion_indices_f8



    
    recursive subroutine qsort_indices_f8_r(keys, sorted_indices, lo, hi, cutoff)
        implicit none
        real(real64), intent(in)      :: keys(:)
        integer,      intent(inout)   :: sorted_indices(:)
        integer,      intent(in)      :: lo, hi, cutoff

        integer :: i, j, mid, tmp_index
        real(real64) :: pivot_key

        if (hi - lo + 1 <= cutoff) then
            call insertion_indices_f8_r(keys, sorted_indices, lo, hi)
            return
        end if

        mid       = (lo + hi) / 2
        pivot_key = keys(sorted_indices(mid))

        i = lo
        j = hi
        do
            do while (keys(sorted_indices(i)) > pivot_key); i = i + 1; end do
            do while (keys(sorted_indices(j)) < pivot_key); j = j - 1; end do

            if (i <= j) then
                tmp_index = sorted_indices(i); sorted_indices(i) = sorted_indices(j); sorted_indices(j) = tmp_index
                i = i + 1
                j = j - 1
            end if

            if (i > j) exit
        end do

        if (lo < j) call qsort_indices_f8_r(keys, sorted_indices, lo, j, cutoff)
        if (i  < hi) call qsort_indices_f8_r(keys, sorted_indices, i,  hi, cutoff)
    end subroutine qsort_indices_f8_r


    subroutine insertion_indices_f8_r(keys, sorted_indices, lo, hi)
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none

        real(real64), intent(in)      :: keys(:)
        integer,      intent(inout)   :: sorted_indices(:)
        integer,      intent(in)      :: lo, hi

        integer :: i
        integer :: j
        integer :: held_index

        do i = lo + 1, hi
            held_index = sorted_indices(i)
            j = i - 1

            do while (j >= lo)
                if (keys(sorted_indices(j)) >= keys(held_index)) exit
                sorted_indices(j+1) = sorted_indices(j)
                j = j - 1
            end do

            sorted_indices(j+1) = held_index
        end do
    end subroutine insertion_indices_f8_r


end module sort_mod
