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
    !> Volume of every cell or tile, for the moment m_i = V_i * Ms_i * mhat_i.
    real(c_double), dimension(:), allocatable, save :: cell_volume
#endif
    logical, save :: dip_fmm_active = .false.

contains

    !> Construct the fixed-geometry cuboid plan and persistent moment buffers.
    !>
    !> Uniform grid with the uniform tree and no root: one common cell size on the bounding
    !> box (the original path). Otherwise every cell passes its own full side lengths:
    !> grid%abc on an unstructuredPrisms (octree) grid, or the common size repeated on a
    !> uniform grid. problem%cdfmm_tree selects the fixed-depth uniform tree (0) or dip-fmm's
    !> capacity-driven adaptive tree (1); with problem%cdfmm_root the uniform tree is every
    !> occupied box of that root split down to cdfmm_depth (adaptive builder, capacity 1).
    subroutine initialiseDipFmm(problem, ierr, message)
        type(MicroMagProblem), intent(in) :: problem
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message
#if USE_CDFMM
        type(cdfmm_options_t) :: options
        real(c_double) :: cell_size(3)
        real(c_double), dimension(:), allocatable :: hx, hy, hz, xs, ys, zs
        real(c_double), dimension(:,:), allocatable :: vertices
        real(c_double) :: root_centre(3), origin(3), unit, corner(3, 4), centroid(3)
        integer(c_int) :: cdfmm_ierr
        integer :: ntot, i, k
        logical :: uniform_cells, tetrahedra
#endif

        call destroyDipFmm()
        ierr = 0
        message = ''

#if USE_CDFMM
        if (problem%cdfmm_precision /= CDFMM_PRECISION_FLOAT32) then
            ierr = 1
            message = 'MagTense demag arrays are FP32; cdfmm_precision must be FLOAT32 (0)'
            return
        end if
        if (problem%cdfmm_tree /= 0 .and. problem%cdfmm_tree /= 1) then
            ierr = 1
            message = 'cdfmm_tree must be 0 (uniform tree) or 1 (adaptive tree)'
            return
        end if

        ntot = size(problem%grid%pts, dim=1)
        tetrahedra = .false.
        select case (problem%grid%gridType)
        case (gridTypeUniform)
            uniform_cells = .true.
        case (gridTypeTetrahedron)
            uniform_cells = .false.
            tetrahedra = .true.
            if (.not. (allocated(problem%grid%nodes) .and. allocated(problem%grid%elements))) then
                ierr = 1
                message = 'dip-fmm on a tetrahedron grid needs grid%nodes and grid%elements'
                return
            end if
            if (size(problem%grid%elements, dim=1) < 4 .or. size(problem%grid%elements, dim=2) /= ntot) then
                ierr = 1
                message = 'grid%elements must hold the four corner nodes of every tetrahedron'
                return
            end if
        case (gridTypeUnstructuredPrisms)
            uniform_cells = .false.
            if (.not. allocated(problem%grid%abc)) then
                ierr = 1
                message = 'dip-fmm on an unstructuredPrisms grid needs the tile sizes grid%abc'
                return
            end if
            if (size(problem%grid%abc, dim=1) /= ntot) then
                ierr = 1
                message = 'grid%abc must hold one row of side lengths per tile'
                return
            end if
        case default
            ierr = 1
            message = 'dip-fmm supports uniform, unstructuredPrisms and tetrahedron MagTense grids'
            return
        end select
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

        allocate(cell_volume(ntot))
        if (tetrahedra) then
            ! Exact tetrahedra: each element is passed as its centroid plus the offsets of
            ! its four corners. As for the prisms below, the geometry is given in a mesh unit
            ! (here the edge of the cube of mean element volume) measured from the lower
            ! corner of the nodes (or of the root), and the moments use the volumes in the
            ! same unit, so the field is unchanged. Irregular tetrahedra have no exactly
            ! equal pairs, so nothing is snapped.
            allocate(xs(ntot), ys(ntot), zs(ntot), vertices(12, ntot))
            do i = 1, ntot
                do k = 1, 4
                    corner(:, k) = real(problem%grid%nodes(:, problem%grid%elements(k, i)), c_double)
                end do
                centroid = 0.25_c_double * sum(corner, dim=2)
                xs(i) = centroid(1)
                ys(i) = centroid(2)
                zs(i) = centroid(3)
                do k = 1, 4
                    vertices(3 * k - 2:3 * k, i) = corner(:, k) - centroid
                end do
                cell_volume(i) = abs(dot_product(corner(:, 4) - corner(:, 1), &
                    cross3(corner(:, 2) - corner(:, 1), corner(:, 3) - corner(:, 1)))) / 6.0_c_double
            end do
            unit = (sum(cell_volume) / real(ntot, c_double))**(1.0_c_double / 3.0_c_double)
            if (.not. (minval(cell_volume) > 0.0_c_double)) then
                ierr = 1
                message = 'dip-fmm needs tetrahedra of positive volume'
                deallocate(cell_volume, xs, ys, zs, vertices)
                return
            end if
            if (problem%cdfmm_root(4) > 0.0_DP) then
                origin = real(problem%cdfmm_root(1:3) - problem%cdfmm_root(4), c_double)
            else
                origin = real(minval(problem%grid%nodes(1:3, :), dim=2), c_double)
            end if
            xs = (xs - origin(1)) / unit
            ys = (ys - origin(2)) / unit
            zs = (zs - origin(3)) / unit
            vertices = vertices / unit
            cell_volume = cell_volume / unit**3
            if (problem%cdfmm_root(4) > 0.0_DP) then
                root_centre = (real(problem%cdfmm_root(1:3), c_double) - origin) / unit
                if (problem%cdfmm_tree == 1) then
                    call cdfmm_create_adaptive_tetrahedra( &
                        dip_fmm_plan, xs, ys, zs, vertices, problem%cdfmm_capacity, &
                        problem%cdfmm_max_depth, options, cdfmm_ierr, root_centre, &
                        real(problem%cdfmm_root(4), c_double) / unit)
                else
                    ! Uniform tree of cdfmm_depth on the given root (capacity 1), as for prisms.
                    call cdfmm_create_adaptive_tetrahedra( &
                        dip_fmm_plan, xs, ys, zs, vertices, 1, problem%cdfmm_depth, options, &
                        cdfmm_ierr, root_centre, real(problem%cdfmm_root(4), c_double) / unit)
                end if
            else if (problem%cdfmm_tree == 1) then
                call cdfmm_create_adaptive_tetrahedra( &
                    dip_fmm_plan, xs, ys, zs, vertices, problem%cdfmm_capacity, &
                    problem%cdfmm_max_depth, options, cdfmm_ierr)
            else
                call cdfmm_create_tetrahedra(dip_fmm_plan, xs, ys, zs, vertices, options, cdfmm_ierr)
            end if
            deallocate(xs, ys, zs, vertices)
        else if (uniform_cells .and. problem%cdfmm_tree == 0 .and. .not. (problem%cdfmm_root(4) > 0.0_DP)) then
            ! The original path: one common cuboid size on the fixed-depth tree, whose root
            ! dip-fmm takes as the bounding box of the cells.
            cell_volume = cell_size(1) * cell_size(2) * cell_size(3)
            call cdfmm_create_uniform_cuboids( &
                dip_fmm_plan, &
                problem%grid%pts(:, 1), &
                problem%grid%pts(:, 2), &
                problem%grid%pts(:, 3), &
                cell_size, options, cdfmm_ierr)
        else
            allocate(hx(ntot), hy(ntot), hz(ntot))
            if (uniform_cells) then
                hx = cell_size(1)
                hy = cell_size(2)
                hz = cell_size(3)
            else
                hx = real(problem%grid%abc(:, 1), c_double)
                hy = real(problem%grid%abc(:, 2), c_double)
                hz = real(problem%grid%abc(:, 3), c_double)
            end if
            ! dip-fmm builds one exact near-field operator per class of pairs whose
            ! displacement and sizes agree bit for bit. Mesh coordinates in metres are not
            ! exactly representable (2.5e-9 m), so equal pairs of an octree mesh would all
            ! differ in their last bits and every pair would be built separately (290 s
            ! instead of 15 s for 324k tiles). The plan therefore gets the geometry in
            ! units of the smallest tile edge, measured from the lower corner of the root
            ! (or of the tiles), snapped to 2^-20 of that unit: a dyadic mesh becomes exact
            ! and nothing moves by more than 5e-7 of a tile. The field is unchanged, since
            ! the moments use the volumes in the same unit (m / r^3 is scale invariant).
            unit = minval(min(hx, min(hy, hz)))
            if (.not. (unit > 0.0_c_double)) then
                ierr = 1
                message = 'dip-fmm needs positive tile sizes'
                deallocate(cell_volume, hx, hy, hz)
                return
            end if
            if (problem%cdfmm_root(4) > 0.0_DP) then
                origin = real(problem%cdfmm_root(1:3) - problem%cdfmm_root(4), c_double)
            else
                origin = [minval(problem%grid%pts(:, 1) - 0.5_c_double * hx), &
                          minval(problem%grid%pts(:, 2) - 0.5_c_double * hy), &
                          minval(problem%grid%pts(:, 3) - 0.5_c_double * hz)]
            end if
            allocate(xs(ntot), ys(ntot), zs(ntot))
            xs = snap((problem%grid%pts(:, 1) - origin(1)) / unit)
            ys = snap((problem%grid%pts(:, 2) - origin(2)) / unit)
            zs = snap((problem%grid%pts(:, 3) - origin(3)) / unit)
            hx = snap(hx / unit)
            hy = snap(hy / unit)
            hz = snap(hz / unit)
            cell_volume = hx * hy * hz
            if (problem%cdfmm_tree == 0 .and. problem%cdfmm_root(4) > 0.0_DP) then
                ! Uniform tree on a given (aligned) root: the adaptive builder with capacity 1
                ! splits every occupied box down to cdfmm_depth, i.e. the uniform tree of that
                ! depth on this root, with the empty boxes left out. Its leaves hold whole cells
                ! when the root is aligned, which the bounding box of the cells does not ensure.
                root_centre = snap((real(problem%cdfmm_root(1:3), c_double) - origin) / unit)
                call cdfmm_create_adaptive_variable_cuboids( &
                    dip_fmm_plan, xs, ys, zs, hx, hy, hz, 1, problem%cdfmm_depth, options, cdfmm_ierr, &
                    root_centre, snap(real(problem%cdfmm_root(4), c_double) / unit))
            else if (problem%cdfmm_tree == 1) then
                if (problem%cdfmm_root(4) > 0.0_DP) then
                    root_centre = snap((real(problem%cdfmm_root(1:3), c_double) - origin) / unit)
                    call cdfmm_create_adaptive_variable_cuboids( &
                        dip_fmm_plan, xs, ys, zs, hx, hy, hz, problem%cdfmm_capacity, &
                        problem%cdfmm_max_depth, options, cdfmm_ierr, root_centre, &
                        snap(real(problem%cdfmm_root(4), c_double) / unit))
                else
                    call cdfmm_create_adaptive_variable_cuboids( &
                        dip_fmm_plan, xs, ys, zs, hx, hy, hz, problem%cdfmm_capacity, &
                        problem%cdfmm_max_depth, options, cdfmm_ierr)
                end if
            else
                call cdfmm_create_variable_cuboids( &
                    dip_fmm_plan, xs, ys, zs, hx, hy, hz, options, cdfmm_ierr)
            end if
            deallocate(hx, hy, hz, xs, ys, zs)
        end if
        if (cdfmm_ierr /= CDFMM_SUCCESS) then
            ierr = int(cdfmm_ierr)
            message = cdfmm_last_error()
            deallocate(cell_volume)
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

        ! dip-fmm consumes total magnetic moment m_i = V_i * Ms_i * mhat_i (per-cell
        ! volumes, set when the plan was created).
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
        if (allocated(cell_volume)) deallocate(cell_volume)
#endif
        dip_fmm_active = .false.
    end subroutine destroyDipFmm

    logical function dipFmmIsActive()
        dipFmmIsActive = dip_fmm_active
    end function dipFmmIsActive

#if USE_CDFMM
    !> Round to the nearest multiple of 2^-20 (exact in double for |x| < 2^32).
    elemental function snap(x) result(y)
        real(c_double), intent(in) :: x
        real(c_double) :: y
        real(c_double), parameter :: grid = 1048576.0_c_double
        y = anint(x * grid) / grid
    end function snap

    pure function cross3(a, b) result(c)
        real(c_double), intent(in) :: a(3), b(3)
        real(c_double) :: c(3)
        c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
    end function cross3
#endif

end module DipFmmDemag
