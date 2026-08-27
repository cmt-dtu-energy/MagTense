!> Adapter between MagTense micromagnetics and dip-fmm's persistent Fortran API.
module DipFmmDemag
    use, intrinsic :: iso_c_binding, only: c_double, c_float, c_int
    use MicroMagParameters
#if USE_CDFMM
    use cdfmm_fortran
#endif
    implicit none
    private

    public :: initialiseDipFmm, evaluateDipFmm, destroyDipFmm
    public :: dipFmmIsActive

#if USE_CDFMM
    type(cdfmm_plan_t), save :: dip_fmm_plan
    real(c_float), dimension(:), allocatable, save :: moment_x, moment_y, moment_z
#endif
    logical, save :: dip_fmm_active = .false.

contains

    !> Construct the fixed-geometry cuboid plan and persistent moment buffers.
    subroutine initialiseDipFmm(problem, ierr, message)
        type(MicroMagProblem), intent(in) :: problem
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message
#if USE_CDFMM
        type(cdfmm_options_t) :: options
        real(c_double) :: cell_size(3)
        integer(c_int) :: cdfmm_ierr
        integer :: ntot
#endif

        call destroyDipFmm()
        ierr = 0
        message = ''

#if USE_CDFMM
        if (problem%grid%gridType /= gridTypeUniform) then
            ierr = 1
            message = 'dip-fmm currently supports only uniform MagTense grids'
            return
        end if
        if (problem%cdfmm_precision /= CDFMM_PRECISION_FLOAT32) then
            ierr = 1
            message = 'MagTense demag arrays are FP32; cdfmm_precision must be FLOAT32 (0)'
            return
        end if

        ntot = size(problem%grid%pts, dim=1)
        cell_size = [real(problem%grid%dx, c_double), &
                     real(problem%grid%dy, c_double), &
                     real(problem%grid%dz, c_double)]

        options = cdfmm_options()
        options%order = problem%cdfmm_order
        options%depth = problem%cdfmm_depth
        options%basis = problem%cdfmm_basis
        options%precision = problem%cdfmm_precision
        options%backend = CDFMM_BACKEND_CPU_STATIC
        if (cdfmm_one_mkl_available()) then
            options%static_matrix_backend = CDFMM_STATIC_MATRIX_ONE_MKL
        else
            options%static_matrix_backend = CDFMM_STATIC_MATRIX_PORTABLE
        end if
#if USE_CUDA
        if (problem%useCuda == useCudaTrue) then
            options%backend = CDFMM_BACKEND_CUDA_FULL
        end if
#endif

        call cdfmm_create_uniform_cuboids( &
            dip_fmm_plan, &
            problem%grid%pts(:, 1), &
            problem%grid%pts(:, 2), &
            problem%grid%pts(:, 3), &
            cell_size, options, cdfmm_ierr)
        if (cdfmm_ierr /= CDFMM_SUCCESS) then
            ierr = int(cdfmm_ierr)
            message = cdfmm_last_error()
            return
        end if

        allocate(moment_x(ntot), moment_y(ntot), moment_z(ntot))
        moment_x = 0.0_c_float
        moment_y = 0.0_c_float
        moment_z = 0.0_c_float
        dip_fmm_active = .true.
#else
        ierr = 1
        message = 'MagTense was compiled without USE_CDFMM=1'
#endif
    end subroutine initialiseDipFmm

    !> Evaluate H for the current normalised magnetisation without rebuilding geometry.
    subroutine evaluateDipFmm(problem, solution, ierr, message)
        type(MicroMagProblem), intent(in) :: problem
        type(MicroMagSolution), intent(inout) :: solution
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message
#if USE_CDFMM
        integer(c_int) :: cdfmm_ierr
        real(c_double) :: cell_volume
#endif

        ierr = 0
        message = ''

        if (problem%useDemag == useDemagFalse) then
            solution%HmX = 0.0_SP
            solution%HmY = 0.0_SP
            solution%HmZ = 0.0_SP
            return
        end if

#if USE_CDFMM
        if (.not. dip_fmm_active) then
            ierr = 1
            message = 'dip-fmm evaluation requested without an active plan'
            return
        end if

        ! dip-fmm consumes total magnetic moment m_i = dV * Ms_i * mhat_i.
        cell_volume = real(problem%grid%dx, c_double) * &
                      real(problem%grid%dy, c_double) * &
                      real(problem%grid%dz, c_double)
        moment_x = real(solution%Mx * problem%Ms * cell_volume, c_float)
        moment_y = real(solution%My * problem%Ms * cell_volume, c_float)
        moment_z = real(solution%Mz * problem%Ms * cell_volume, c_float)

        ! cdfmm_evaluate already returns H=-grad(phi) with G=1/(4*pi*r).
        call cdfmm_evaluate( &
            dip_fmm_plan, moment_x, moment_y, moment_z, &
            solution%HmX, solution%HmY, solution%HmZ, cdfmm_ierr)
        if (cdfmm_ierr /= CDFMM_SUCCESS) then
            ierr = int(cdfmm_ierr)
            message = cdfmm_last_error()
        end if
#else
        ierr = 1
        message = 'MagTense was compiled without USE_CDFMM=1'
#endif
    end subroutine evaluateDipFmm

    !> Release all dip-fmm state. Safe to call when no plan is active.
    subroutine destroyDipFmm()
#if USE_CDFMM
        call dip_fmm_plan%destroy()
        if (allocated(moment_x)) deallocate(moment_x)
        if (allocated(moment_y)) deallocate(moment_y)
        if (allocated(moment_z)) deallocate(moment_z)
#endif
        dip_fmm_active = .false.
    end subroutine destroyDipFmm

    logical function dipFmmIsActive()
        dipFmmIsActive = dip_fmm_active
    end function dipFmmIsActive

end module DipFmmDemag
