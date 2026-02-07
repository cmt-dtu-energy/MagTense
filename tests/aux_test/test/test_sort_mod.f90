module test_quicksort_mod
    use, intrinsic :: iso_fortran_env, only: int32, real32, real64
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use sort_mod,  only: sort, sort_r, argsort, argsort_r, apply_perm
    implicit none
    private

    public :: collect_quicksort_tests

    integer, parameter :: N_RANDOM = 1000

contains

    !===============================================================================================
    !> Collect all quicksort unit tests.
    !===============================================================================================
    subroutine collect_quicksort_tests(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            !---------------------------------------------------------------------------------------
            ! sort: integer(int32)
            !---------------------------------------------------------------------------------------
            new_unittest("quicksort: sort i4 ascending (default cutoff)",    test_sort_i4_asc_default), &
            new_unittest("quicksort: sort i4 ascending (cutoff=1)",         test_sort_i4_asc_cut1), &
            new_unittest("quicksort: sort i4 ascending (cutoff=8)",         test_sort_i4_asc_cut8), &
            new_unittest("quicksort: sort i4 ascending (cutoff=64)",        test_sort_i4_asc_cut64), &
            new_unittest("quicksort: sort i4 descending (default cutoff)",  test_sort_i4_desc_default), &
            new_unittest("quicksort: sort i4 descending (cutoff=8)",        test_sort_i4_desc_cut8), &
            new_unittest("quicksort: sort i4 descending (cutoff=64)",       test_sort_i4_desc_cut64), &
            !---------------------------------------------------------------------------------------
            ! argsort: integer(int32)
            !---------------------------------------------------------------------------------------
            new_unittest("quicksort: argsort i4 ascending (default cutoff)", test_argsort_i4_asc_default), &
            new_unittest("quicksort: argsort i4 ascending (cutoff=1)",      test_argsort_i4_asc_cut1), &
            new_unittest("quicksort: argsort i4 ascending (cutoff=64)",     test_argsort_i4_asc_cut64), &
            new_unittest("quicksort: argsort i4 descending (default cutoff)",test_argsort_i4_desc_default), &
            new_unittest("quicksort: argsort i4 descending (cutoff=8)",     test_argsort_i4_desc_cut8), &
            !---------------------------------------------------------------------------------------
            ! sort: real(real32)
            !---------------------------------------------------------------------------------------
            new_unittest("quicksort: sort f4 ascending (default cutoff)",    test_sort_f4_asc_default), &
            new_unittest("quicksort: sort f4 ascending (cutoff=8)",         test_sort_f4_asc_cut8), &
            new_unittest("quicksort: sort f4 descending (cutoff=1)",        test_sort_f4_desc_cut1), &
            new_unittest("quicksort: sort f4 descending (cutoff=64)",       test_sort_f4_desc_cut64), &
            !---------------------------------------------------------------------------------------
            ! argsort: real(real32)
            !---------------------------------------------------------------------------------------
            new_unittest("quicksort: argsort f4 ascending (cutoff=1)",       test_argsort_f4_asc_cut1), &
            new_unittest("quicksort: argsort f4 descending (default cutoff)",test_argsort_f4_desc_default), &
            !---------------------------------------------------------------------------------------
            ! sort: real(real64)
            !---------------------------------------------------------------------------------------
            new_unittest("quicksort: sort f8 ascending (default cutoff)",    test_sort_f8_asc_default), &
            new_unittest("quicksort: sort f8 ascending (cutoff=64)",        test_sort_f8_asc_cut64), &
            new_unittest("quicksort: sort f8 descending (cutoff=8)",        test_sort_f8_desc_cut8), &
            new_unittest("quicksort: sort f8 descending (cutoff=64)",       test_sort_f8_desc_cut64), &
            !---------------------------------------------------------------------------------------
            ! argsort: real(real64)
            !---------------------------------------------------------------------------------------
            new_unittest("quicksort: argsort f8 ascending (default cutoff)", test_argsort_f8_asc_default), &
            new_unittest("quicksort: argsort f8 descending (cutoff=1)",      test_argsort_f8_desc_cut1) &
        ]
    end subroutine collect_quicksort_tests


    !===============================================================================================
    ! Tests: sort integer(int32)
    !===============================================================================================
    subroutine test_sort_i4_asc_default(error)
        type(error_type), allocatable, intent(out) :: error
        call run_sort_i4_case(error, cutoff_present=.false., cutoff_value=0, descending=.false., seed_tag=11001)
    end subroutine test_sort_i4_asc_default

    subroutine test_sort_i4_asc_cut1(error)
        type(error_type), allocatable, intent(out) :: error
        call run_sort_i4_case(error, cutoff_present=.true., cutoff_value=1, descending=.false., seed_tag=11002)
    end subroutine test_sort_i4_asc_cut1

    subroutine test_sort_i4_asc_cut8(error)
        type(error_type), allocatable, intent(out) :: error
        call run_sort_i4_case(error, cutoff_present=.true., cutoff_value=8, descending=.false., seed_tag=11003)
    end subroutine test_sort_i4_asc_cut8

    subroutine test_sort_i4_asc_cut64(error)
        type(error_type), allocatable, intent(out) :: error
        call run_sort_i4_case(error, cutoff_present=.true., cutoff_value=64, descending=.false., seed_tag=11004)
    end subroutine test_sort_i4_asc_cut64

    subroutine test_sort_i4_desc_default(error)
        type(error_type), allocatable, intent(out) :: error
        call run_sort_i4_case(error, cutoff_present=.false., cutoff_value=0, descending=.true., seed_tag=11005)
    end subroutine test_sort_i4_desc_default

    subroutine test_sort_i4_desc_cut8(error)
        type(error_type), allocatable, intent(out) :: error
        call run_sort_i4_case(error, cutoff_present=.true., cutoff_value=8, descending=.true., seed_tag=11006)
    end subroutine test_sort_i4_desc_cut8

    subroutine test_sort_i4_desc_cut64(error)
        type(error_type), allocatable, intent(out) :: error
        call run_sort_i4_case(error, cutoff_present=.true., cutoff_value=64, descending=.true., seed_tag=11007)
    end subroutine test_sort_i4_desc_cut64


    !===============================================================================================
    ! Tests: argsort integer(int32)
    !===============================================================================================
    subroutine test_argsort_i4_asc_default(error)
        type(error_type), allocatable, intent(out) :: error
        call run_argsort_i4_case(error, cutoff_present=.false., cutoff_value=0, descending=.false., seed_tag=12001)
    end subroutine test_argsort_i4_asc_default

    subroutine test_argsort_i4_asc_cut1(error)
        type(error_type), allocatable, intent(out) :: error
        call run_argsort_i4_case(error, cutoff_present=.true., cutoff_value=1, descending=.false., seed_tag=12002)
    end subroutine test_argsort_i4_asc_cut1

    subroutine test_argsort_i4_asc_cut64(error)
        type(error_type), allocatable, intent(out) :: error
        call run_argsort_i4_case(error, cutoff_present=.true., cutoff_value=64, descending=.false., seed_tag=12003)
    end subroutine test_argsort_i4_asc_cut64

    subroutine test_argsort_i4_desc_default(error)
        type(error_type), allocatable, intent(out) :: error
        call run_argsort_i4_case(error, cutoff_present=.false., cutoff_value=0, descending=.true., seed_tag=12004)
    end subroutine test_argsort_i4_desc_default

    subroutine test_argsort_i4_desc_cut8(error)
        type(error_type), allocatable, intent(out) :: error
        call run_argsort_i4_case(error, cutoff_present=.true., cutoff_value=8, descending=.true., seed_tag=12005)
    end subroutine test_argsort_i4_desc_cut8


    !===============================================================================================
    ! Tests: sort real(real32)
    !===============================================================================================
    subroutine test_sort_f4_asc_default(error)
        type(error_type), allocatable, intent(out) :: error
        call run_sort_f4_case(error, cutoff_present=.false., cutoff_value=0, descending=.false., seed_tag=21001)
    end subroutine test_sort_f4_asc_default

    subroutine test_sort_f4_asc_cut8(error)
        type(error_type), allocatable, intent(out) :: error
        call run_sort_f4_case(error, cutoff_present=.true., cutoff_value=8, descending=.false., seed_tag=21002)
    end subroutine test_sort_f4_asc_cut8

    subroutine test_sort_f4_desc_cut1(error)
        type(error_type), allocatable, intent(out) :: error
        call run_sort_f4_case(error, cutoff_present=.true., cutoff_value=1, descending=.true., seed_tag=21003)
    end subroutine test_sort_f4_desc_cut1

    subroutine test_sort_f4_desc_cut64(error)
        type(error_type), allocatable, intent(out) :: error
        call run_sort_f4_case(error, cutoff_present=.true., cutoff_value=64, descending=.true., seed_tag=21004)
    end subroutine test_sort_f4_desc_cut64


    !===============================================================================================
    ! Tests: argsort real(real32)
    !===============================================================================================
    subroutine test_argsort_f4_asc_cut1(error)
        type(error_type), allocatable, intent(out) :: error
        call run_argsort_f4_case(error, cutoff_present=.true., cutoff_value=1, descending=.false., seed_tag=22001)
    end subroutine test_argsort_f4_asc_cut1

    subroutine test_argsort_f4_desc_default(error)
        type(error_type), allocatable, intent(out) :: error
        call run_argsort_f4_case(error, cutoff_present=.false., cutoff_value=0, descending=.true., seed_tag=22002)
    end subroutine test_argsort_f4_desc_default


    !===============================================================================================
    ! Tests: sort real(real64)
    !===============================================================================================
    subroutine test_sort_f8_asc_default(error)
        type(error_type), allocatable, intent(out) :: error
        call run_sort_f8_case(error, cutoff_present=.false., cutoff_value=0, descending=.false., seed_tag=31001)
    end subroutine test_sort_f8_asc_default

    subroutine test_sort_f8_asc_cut64(error)
        type(error_type), allocatable, intent(out) :: error
        call run_sort_f8_case(error, cutoff_present=.true., cutoff_value=64, descending=.false., seed_tag=31002)
    end subroutine test_sort_f8_asc_cut64

    subroutine test_sort_f8_desc_cut8(error)
        type(error_type), allocatable, intent(out) :: error
        call run_sort_f8_case(error, cutoff_present=.true., cutoff_value=8, descending=.true., seed_tag=31003)
    end subroutine test_sort_f8_desc_cut8

    subroutine test_sort_f8_desc_cut64(error)
        type(error_type), allocatable, intent(out) :: error
        call run_sort_f8_case(error, cutoff_present=.true., cutoff_value=64, descending=.true., seed_tag=31004)
    end subroutine test_sort_f8_desc_cut64


    !===============================================================================================
    ! Tests: argsort real(real64)
    !===============================================================================================
    subroutine test_argsort_f8_asc_default(error)
        type(error_type), allocatable, intent(out) :: error
        call run_argsort_f8_case(error, cutoff_present=.false., cutoff_value=0, descending=.false., seed_tag=32001)
    end subroutine test_argsort_f8_asc_default

    subroutine test_argsort_f8_desc_cut1(error)
        type(error_type), allocatable, intent(out) :: error
        call run_argsort_f8_case(error, cutoff_present=.true., cutoff_value=1, descending=.true., seed_tag=32002)
    end subroutine test_argsort_f8_desc_cut1


    !===============================================================================================
    ! RNG utilities (deterministic per test)
    !===============================================================================================
    subroutine seed_rng(tag)
        integer, intent(in) :: tag

        integer :: n_seed
        integer, allocatable :: seed(:)
        integer :: i
        integer :: x

        call random_seed(size=n_seed)
        allocate(seed(n_seed))

        x = tag
        do i = 1, n_seed
            ! Simple integer LCG step; keep positive.
            x = iand(1103515245 * x + 12345 + 97*i, huge(1))
            if (x <= 0) x = -x + 1
            seed(i) = x
        end do

        call random_seed(put=seed)
    end subroutine seed_rng


    subroutine fill_random_i4(values, seed_tag)
        integer(int32), intent(out) :: values(:)
        integer,        intent(in)  :: seed_tag

        real(real64) :: r(size(values))
        integer :: i

        call seed_rng(seed_tag)
        call random_number(r)

        do i = 1, size(values)
            ! Range approx [-1_000_000, 1_000_000]
            values(i) = int(2.0e6_real64 * (r(i) - 0.5_real64), int32)
        end do
    end subroutine fill_random_i4


    subroutine fill_random_f4(values, seed_tag)
        real(real32), intent(out) :: values(:)
        integer,      intent(in)  :: seed_tag

        call seed_rng(seed_tag)
        call random_number(values)

        values = 1.0e3_real32 * (values - 0.5_real32)
    end subroutine fill_random_f4


    subroutine fill_random_f8(values, seed_tag)
        real(real64), intent(out) :: values(:)
        integer,      intent(in)  :: seed_tag

        call seed_rng(seed_tag)
        call random_number(values)

        values = 1.0e3_real64 * (values - 0.5_real64)
    end subroutine fill_random_f8


    !===============================================================================================
    ! Core test runners: integer(int32)
    !===============================================================================================
    subroutine run_sort_i4_case(error, cutoff_present, cutoff_value, descending, seed_tag)
        type(error_type), allocatable, intent(out) :: error
        logical, intent(in) :: cutoff_present
        integer, intent(in) :: cutoff_value
        logical, intent(in) :: descending
        integer, intent(in) :: seed_tag

        integer(int32) :: values(N_RANDOM)

        call fill_random_i4(values, seed_tag)

        if (.not. descending) then
            if (cutoff_present) then
                call sort(values, algo="quicksort", cutoff=cutoff_value)
            else
                call sort(values, algo="quicksort")
            end if
            call assert_sorted_i4_asc(values, error)
        else
            if (cutoff_present) then
                call sort_r(values, algo="quicksort", cutoff=cutoff_value)
            else
                call sort_r(values, algo="quicksort")
            end if
            call assert_sorted_i4_desc(values, error)
        end if
    end subroutine run_sort_i4_case


    subroutine run_argsort_i4_case(error, cutoff_present, cutoff_value, descending, seed_tag)
        type(error_type), allocatable, intent(out) :: error
        logical, intent(in) :: cutoff_present
        integer, intent(in) :: cutoff_value
        logical, intent(in) :: descending
        integer, intent(in) :: seed_tag

        integer(int32) :: keys(N_RANDOM)
        integer        :: sorted_indices(N_RANDOM)
        integer(int32) :: keys_via_perm(N_RANDOM)
        integer(int32) :: keys_via_sort(N_RANDOM)

        call fill_random_i4(keys, seed_tag)

        if (.not. descending) then
            if (cutoff_present) then
                call argsort(keys, sorted_indices, algo="quicksort", cutoff=cutoff_value)
                call sort(keys, keys_via_sort, algo="quicksort", cutoff=cutoff_value)
            else
                call argsort(keys, sorted_indices, algo="quicksort")
                call sort(keys, keys_via_sort, algo="quicksort")
            end if
            call apply_perm(keys, sorted_indices, keys_via_perm)
            call assert_sorted_i4_asc(keys_via_perm, error)
            if (allocated(error)) return
            call assert_equal_i4(keys_via_perm, keys_via_sort, error)
        else
            if (cutoff_present) then
                call argsort_r(keys, sorted_indices, algo="quicksort", cutoff=cutoff_value)
                call sort_r(keys, keys_via_sort, algo="quicksort", cutoff=cutoff_value)
            else
                call argsort_r(keys, sorted_indices, algo="quicksort")
                call sort_r(keys, keys_via_sort, algo="quicksort")
            end if
            call apply_perm(keys, sorted_indices, keys_via_perm)
            call assert_sorted_i4_desc(keys_via_perm, error)
            if (allocated(error)) return
            call assert_equal_i4(keys_via_perm, keys_via_sort, error)
        end if

        if (allocated(error)) return
        call assert_is_permutation(sorted_indices, N_RANDOM, error)
    end subroutine run_argsort_i4_case


    subroutine assert_sorted_i4_asc(values, error)
        integer(int32), intent(in) :: values(:)
        type(error_type), allocatable, intent(out) :: error
        integer :: i

        do i = 2, size(values)
            call check(error, values(i-1) <= values(i), "not sorted ascending")
            if (allocated(error)) return
        end do
    end subroutine assert_sorted_i4_asc


    subroutine assert_sorted_i4_desc(values, error)
        integer(int32), intent(in) :: values(:)
        type(error_type), allocatable, intent(out) :: error
        integer :: i

        do i = 2, size(values)
            call check(error, values(i-1) >= values(i), "not sorted descending")
            if (allocated(error)) return
        end do
    end subroutine assert_sorted_i4_desc


    subroutine assert_equal_i4(a, b, error)
        integer(int32), intent(in) :: a(:), b(:)
        type(error_type), allocatable, intent(out) :: error

        call check(error, size(a) == size(b), "size mismatch")
        if (allocated(error)) return

        call check(error, all(a == b), "arrays differ")
    end subroutine assert_equal_i4


    !===============================================================================================
    ! Core test runners: real(real32)
    !===============================================================================================
    subroutine run_sort_f4_case(error, cutoff_present, cutoff_value, descending, seed_tag)
        type(error_type), allocatable, intent(out) :: error
        logical, intent(in) :: cutoff_present
        integer, intent(in) :: cutoff_value
        logical, intent(in) :: descending
        integer, intent(in) :: seed_tag

        real(real32) :: values(N_RANDOM)

        call fill_random_f4(values, seed_tag)

        if (.not. descending) then
            if (cutoff_present) then
                call sort(values, algo="quicksort", cutoff=cutoff_value)
            else
                call sort(values, algo="quicksort")
            end if
            call assert_sorted_f4_asc(values, error)
        else
            if (cutoff_present) then
                call sort_r(values, algo="quicksort", cutoff=cutoff_value)
            else
                call sort_r(values, algo="quicksort")
            end if
            call assert_sorted_f4_desc(values, error)
        end if
    end subroutine run_sort_f4_case


    subroutine run_argsort_f4_case(error, cutoff_present, cutoff_value, descending, seed_tag)
        type(error_type), allocatable, intent(out) :: error
        logical, intent(in) :: cutoff_present
        integer, intent(in) :: cutoff_value
        logical, intent(in) :: descending
        integer, intent(in) :: seed_tag

        real(real32) :: keys(N_RANDOM)
        integer      :: sorted_indices(N_RANDOM)
        real(real32) :: keys_via_perm(N_RANDOM)
        real(real32) :: keys_via_sort(N_RANDOM)

        call fill_random_f4(keys, seed_tag)

        if (.not. descending) then
            if (cutoff_present) then
                call argsort(keys, sorted_indices, algo="quicksort", cutoff=cutoff_value)
                call sort(keys, keys_via_sort, algo="quicksort", cutoff=cutoff_value)
            else
                call argsort(keys, sorted_indices, algo="quicksort")
                call sort(keys, keys_via_sort, algo="quicksort")
            end if
            call apply_perm(keys, sorted_indices, keys_via_perm)
            call assert_sorted_f4_asc(keys_via_perm, error)
            if (allocated(error)) return
            call assert_equal_f4(keys_via_perm, keys_via_sort, error)
        else
            if (cutoff_present) then
                call argsort_r(keys, sorted_indices, algo="quicksort", cutoff=cutoff_value)
                call sort_r(keys, keys_via_sort, algo="quicksort", cutoff=cutoff_value)
            else
                call argsort_r(keys, sorted_indices, algo="quicksort")
                call sort_r(keys, keys_via_sort, algo="quicksort")
            end if
            call apply_perm(keys, sorted_indices, keys_via_perm)
            call assert_sorted_f4_desc(keys_via_perm, error)
            if (allocated(error)) return
            call assert_equal_f4(keys_via_perm, keys_via_sort, error)
        end if

        if (allocated(error)) return
        call assert_is_permutation(sorted_indices, N_RANDOM, error)
    end subroutine run_argsort_f4_case


    subroutine assert_sorted_f4_asc(values, error)
        real(real32), intent(in) :: values(:)
        type(error_type), allocatable, intent(out) :: error
        integer :: i

        do i = 2, size(values)
            call check(error, values(i-1) <= values(i), "not sorted ascending")
            if (allocated(error)) return
        end do
    end subroutine assert_sorted_f4_asc


    subroutine assert_sorted_f4_desc(values, error)
        real(real32), intent(in) :: values(:)
        type(error_type), allocatable, intent(out) :: error
        integer :: i

        do i = 2, size(values)
            call check(error, values(i-1) >= values(i), "not sorted descending")
            if (allocated(error)) return
        end do
    end subroutine assert_sorted_f4_desc


    subroutine assert_equal_f4(a, b, error)
        real(real32), intent(in) :: a(:), b(:)
        type(error_type), allocatable, intent(out) :: error

        call check(error, size(a) == size(b), "size mismatch")
        if (allocated(error)) return

        call check(error, all(a == b), "arrays differ")
    end subroutine assert_equal_f4


    !===============================================================================================
    ! Core test runners: real(real64)
    !===============================================================================================
    subroutine run_sort_f8_case(error, cutoff_present, cutoff_value, descending, seed_tag)
        type(error_type), allocatable, intent(out) :: error
        logical, intent(in) :: cutoff_present
        integer, intent(in) :: cutoff_value
        logical, intent(in) :: descending
        integer, intent(in) :: seed_tag

        real(real64) :: values(N_RANDOM)

        call fill_random_f8(values, seed_tag)

        if (.not. descending) then
            if (cutoff_present) then
                call sort(values, algo="quicksort", cutoff=cutoff_value)
            else
                call sort(values, algo="quicksort")
            end if
            call assert_sorted_f8_asc(values, error)
        else
            if (cutoff_present) then
                call sort_r(values, algo="quicksort", cutoff=cutoff_value)
            else
                call sort_r(values, algo="quicksort")
            end if
            call assert_sorted_f8_desc(values, error)
        end if
    end subroutine run_sort_f8_case


    subroutine run_argsort_f8_case(error, cutoff_present, cutoff_value, descending, seed_tag)
        type(error_type), allocatable, intent(out) :: error
        logical, intent(in) :: cutoff_present
        integer, intent(in) :: cutoff_value
        logical, intent(in) :: descending
        integer, intent(in) :: seed_tag

        real(real64) :: keys(N_RANDOM)
        integer      :: sorted_indices(N_RANDOM)
        real(real64) :: keys_via_perm(N_RANDOM)
        real(real64) :: keys_via_sort(N_RANDOM)

        call fill_random_f8(keys, seed_tag)

        if (.not. descending) then
            if (cutoff_present) then
                call argsort(keys, sorted_indices, algo="quicksort", cutoff=cutoff_value)
                call sort(keys, keys_via_sort, algo="quicksort", cutoff=cutoff_value)
            else
                call argsort(keys, sorted_indices, algo="quicksort")
                call sort(keys, keys_via_sort, algo="quicksort")
            end if
            call apply_perm(keys, sorted_indices, keys_via_perm)
            call assert_sorted_f8_asc(keys_via_perm, error)
            if (allocated(error)) return
            call assert_equal_f8(keys_via_perm, keys_via_sort, error)
        else
            if (cutoff_present) then
                call argsort_r(keys, sorted_indices, algo="quicksort", cutoff=cutoff_value)
                call sort_r(keys, keys_via_sort, algo="quicksort", cutoff=cutoff_value)
            else
                call argsort_r(keys, sorted_indices, algo="quicksort")
                call sort_r(keys, keys_via_sort, algo="quicksort")
            end if
            call apply_perm(keys, sorted_indices, keys_via_perm)
            call assert_sorted_f8_desc(keys_via_perm, error)
            if (allocated(error)) return
            call assert_equal_f8(keys_via_perm, keys_via_sort, error)
        end if

        if (allocated(error)) return
        call assert_is_permutation(sorted_indices, N_RANDOM, error)
    end subroutine run_argsort_f8_case


    subroutine assert_sorted_f8_asc(values, error)
        real(real64), intent(in) :: values(:)
        type(error_type), allocatable, intent(out) :: error
        integer :: i

        do i = 2, size(values)
            call check(error, values(i-1) <= values(i), "not sorted ascending")
            if (allocated(error)) return
        end do
    end subroutine assert_sorted_f8_asc


    subroutine assert_sorted_f8_desc(values, error)
        real(real64), intent(in) :: values(:)
        type(error_type), allocatable, intent(out) :: error
        integer :: i

        do i = 2, size(values)
            call check(error, values(i-1) >= values(i), "not sorted descending")
            if (allocated(error)) return
        end do
    end subroutine assert_sorted_f8_desc


    subroutine assert_equal_f8(a, b, error)
        real(real64), intent(in) :: a(:), b(:)
        type(error_type), allocatable, intent(out) :: error

        call check(error, size(a) == size(b), "size mismatch")
        if (allocated(error)) return

        call check(error, all(a == b), "arrays differ")
    end subroutine assert_equal_f8


    !===============================================================================================
    ! Helper: permutation check for sorted_indices
    !===============================================================================================
    subroutine assert_is_permutation(indices, n, error)
        integer, intent(in) :: indices(:)
        integer, intent(in) :: n
        type(error_type), allocatable, intent(out) :: error

        logical, allocatable :: seen(:)
        integer :: i
        integer :: v

        call check(error, size(indices) == n, "permutation size mismatch")
        if (allocated(error)) return

        allocate(seen(n))
        seen = .false.

        do i = 1, n
            v = indices(i)

            call check(error, v >= 1 .and. v <= n, "index out of bounds in permutation")
            if (allocated(error)) return

            call check(error, .not. seen(v), "duplicate index in permutation")
            if (allocated(error)) return

            seen(v) = .true.
        end do

        call check(error, all(seen), "permutation missing values")
    end subroutine assert_is_permutation

end module test_quicksort_mod


program tester
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type
    use test_quicksort_mod, only: collect_quicksort_tests
    implicit none

    integer :: stat
    type(testsuite_type), allocatable :: testsuites(:)

    stat = 0

    testsuites = [ &
        new_testsuite("aux/quicksort", collect_quicksort_tests) &
    ]

    call run_testsuite(testsuites(1)%collect, error_unit, stat)

    if (stat > 0) then
        write(error_unit, '(i0,1x,a)') stat, "test(s) failed!"
        error stop 1
    end if
end program tester
