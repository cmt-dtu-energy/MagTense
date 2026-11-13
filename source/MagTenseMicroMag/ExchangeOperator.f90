!>-----------------------------------------
!> Exchange operator: the single entry point of the solver for the exchange matrix.
!>
!>   prepareExchangeMesh    mesh analysis of an unstructured grid (or a passed exchange matrix)
!>   ComputeExchangeTerm3D  builds problem%A_exch with the operator selected by exch_method
!>
!> Operators, one module each:
!>   uniform grid                       ExchangeUniform          7-point Laplacian
!>   MicroMagExchMethodDirectLaplacianNeumann, MicroMagExchMethodGGNeumann
!>                                      ExchangeDirectLaplacian  Direct Laplacian / Green-Gauss
!>   MicroMagExchMethodCellHessian      ExchangeCellHessian      cell Hessian (least squares)
!>   MicroMagExchMethodSushi            ExchangeSushi            symmetric, energy based (SUSHI-g)
!> The cell Hessian and SUSHI-g operators take their neighbours and faces from
!> ExchangeMeshTopology; the Direct Laplacian uses GridInfo from the mesh analysis modules.
!>
!> A new operator for unstructured grids: a module with a routine of the signature
!>   ( grid, pbc, A0n, phase, Aint, GridInfo, Aexch_matrix )
!> (see computeExchangeSushi), a MicroMagExchMethod* constant in MicroMagParameters, a branch in
!> ComputeExchangeTerm3D below, the file in the Makefile and MagTenseMicroMag.vfproj, and the
!> name in the exch_meth table of python/src/magtense/micromag.py.
!---------------------------------------------------------------------------
module ExchangeOperator
    use MKL_SPBLAS
    use MicroMagParameters
    use IO_GENERAL
    use UTIL_MICROMAG
    use trace_mod
    use UnstructuredMeshAnalysis
    use TetrahedralMeshAnalysis
    use ExchangeDirectLaplacian
    use ExchangeCellHessian
    use ExchangeSushi
    use ExchangeUniform
    implicit none
    private

    public :: prepareExchangeMesh, ComputeExchangeTerm3D, passDifferentialOperators

    contains

    !>-----------------------------------------
    !> Prepares an unstructured grid for the exchange operator before the interaction matrices
    !> are built: an exchange matrix passed from outside is converted to CSR, otherwise the mesh
    !> analysis fills GridInfo (faces, normals, areas, volumes) for the Direct Laplacian and
    !> Green-Gauss operators. Uniform grids need nothing here.
    !> @param[inout] problem the problem data structure
    !> @param[inout] solution the solution data structure, whose GridInfo is filled
    !---------------------------------------------------------------------------
    subroutine prepareExchangeMesh( problem, solution )
    type(MicroMagProblem),intent(inout) :: problem      !> Problem data structure
    type(MicroMagSolution),intent(inout) :: solution    !> Solution data structure

        if (( problem%grid%gridType .eq. gridTypeUnstructuredPrisms ) .or. &
            ( problem%grid%gridType .eq. gridTypeTetrahedron )) then
            if ( problem%passExch .eq. passExchTrue) then
                call displayGUIMessage( 'Passing exchange matrix' )
                call passDifferentialOperators(problem)
            elseif ( problem%grid%gridType .eq. gridTypeUnstructuredPrisms ) then
                call CartesianUnstructuredMeshAnalysis(problem%grid%pts, problem%grid%abc, solution%gridinfo, &
                    problem%macrogrid%exchPBC)
            else
                !A tetrahedral mesh is given by its nodes and connectivity rather than by positions
                !and sizes, so it has its own analysis. It fills the same GridInfo, which means the
                !exchange operator below is built in exactly the same way as for the prisms.
                call TetrahedralUnstructuredMeshAnalysis(problem%grid%nodes, problem%grid%elements, &
                    solution%gridinfo, problem%macrogrid%exchPBC)
            endif
        endif

    end subroutine prepareExchangeMesh


    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Calculates the exhange term matrix
    !> which means it produces the differential operator d^2/dx^2 + d^2/dy^2 + d^2/dz^2 and returns this in the sparse matrix A
    !---------------------------------------------------------------------------   
    subroutine ComputeExchangeTerm3D( grid, A, problem, solution )
    type(MicroMagGrid),intent(in) :: grid               !> Struct containing the grid information    
    type(sparse_matrix_t),intent(inout) :: A            !> The returned matrix from the sparse matrix creator
    type(MicroMagProblem),intent(inout) :: problem      !> Problem data structure    
    type(MicroMagSolution),intent(inout) :: solution    !> Solution data structure
    real(DP),dimension(:),allocatable :: A0_normalized  !> Normalized A0 for uneven anisotropy
    real(DP),dimension(:,:),allocatable :: A_int_normalized  !> Interface exchange, same scaling
    integer,dimension(:),allocatable :: phase_local     !> Material index of each tile
    integer :: dlLinLen, dlLinStat                      !> Length and status of MAGTENSE_DL_LINEAR


    if ( problem%passExch .eq. passExchTrue) return ! Skip this if the exchange has already been passed from outside

    !With A0 = 0 in every tile there is no exchange at all, and the operator is identically zero.
    !It has to be caught here: the normalisation below is then 0/0, which is NaN inside the
    !python and MATLAB extensions (where /fpe:0 does not trap), and on an unstructured mesh the
    !NaN zeroes every interpolation weight, all sparse matrices come out empty and the solver
    !segfaults. The uniform grid survived only because it never reads A0_normalized.
    if ( .not. any( problem%A0 .gt. 0.0_DP ) ) then
        call displayGUIMessage( 'No exchange terms (A0 = 0 in every tile)' )
        call ComputeZeroExchangeOperator( grid%nx * grid%ny * grid%nz, A )
        call create_COO_values_from_CSR( A, solution%gridinfo )
        problem%useExchange = .false.
        return
    endif

    allocate( A0_normalized(size(problem%A0)) )
    A0_normalized = problem%A0 / ( maxval(problem%A0) )   ! Normalized by the largest exchange factor

    !The interface values have to be scaled exactly as A0 is, since they replace a harmonic
    !mean of the scaled values. The negative entries are sentinels meaning 'use the harmonic
    !mean' and are passed through untouched so that they stay recognisable after scaling.
    if ( allocated(problem%A_int) ) then
        allocate( A_int_normalized(size(problem%A_int,1),size(problem%A_int,2)) )
        where ( problem%A_int .ge. 0.0_DP )
            A_int_normalized = problem%A_int / ( maxval(problem%A0) )
        elsewhere
            A_int_normalized = -1.0_DP
        end where
    else
        allocate( A_int_normalized(1,1) )
        A_int_normalized = -1.0_DP
    endif

    if ( allocated(problem%phase_id) ) then
        allocate( phase_local(size(problem%phase_id)) )
        phase_local = problem%phase_id
    else
        allocate( phase_local(size(problem%A0)) )
        phase_local = 1
    endif
        
    if ( grid%gridType .eq. gridTypeUniform ) then
        call ComputeExchangeTerm3D_Uniform( grid, A, problem, solution )
    elseif (( grid%gridType .eq. gridTypeTetrahedron ) .or. (grid%gridType .eq. gridTypeUnstructuredPrisms)) then
    
         if ( problem%macrogrid%exchPBC(1) .or. problem%macrogrid%exchPBC(2) .or. problem%macrogrid%exchPBC(3) ) then
            !The tiles at the two ends of each periodic direction have been linked together in the mesh
            !analysis, i.e. they share a face and enter each others interpolation stencils, so the
            !exchange operator below is computed in exactly the same way as without periodic boundaries.
            !Both mesh analyses do this, the prisms geometrically and the tetrahedra by identifying the
            !nodes on the two boundary planes. An exchange matrix passed from the outside carries its
            !own boundary conditions and has already returned above.
            call displayGUIMessage( 'Using periodic exchange boundary conditions' )
        endif
        if ( problem%exch_method .eq. MicroMagExchMethodCellHessian ) then
            !Cell-centred Hessian operator: consistent on irregular tetrahedra and on graded and
            !refined prism grids, free surfaces imposed in the fit, periodic only along exchPBC.
            call computeExchangeCellHessian( problem%grid, problem%macrogrid%exchPBC .ne. 0, A0_normalized, &
                                             phase_local, A_int_normalized, solution%gridinfo, problem%A_exch )
        elseif ( problem%exch_method .eq. MicroMagExchMethodSushi ) then
            !Symmetric, energy based operator (SUSHI-g): the field is exactly -dE/dm, so the exchange
            !energy is conserved without damping and decreases with damping; free surfaces natural,
            !periodic only along exchPBC.
            call computeExchangeSushi( problem%grid, problem%macrogrid%exchPBC .ne. 0, A0_normalized, &
                                       phase_local, A_int_normalized, solution%gridinfo, problem%A_exch )
        else
            !Tetrahedral meshes use a quadratic face fit by default; prism grids use the linear fit of
            !Poulsen et al. 2022 with all gradient components, which is exact for linear fields on
            !uniform, graded and locally refined grids. Setting the environment variable
            !MAGTENSE_DL_LINEAR (to anything) selects that linear fit on tetrahedra too, i.e. the
            !method as published, for comparison.
            call get_environment_variable( 'MAGTENSE_DL_LINEAR', length=dlLinLen, status=dlLinStat )
            call computeDifferentialOperatorsFromMesh_DirectLap(solution%gridinfo, problem%exch_interpn, problem%exch_weight, problem%exch_method, A0_normalized, phase_local, A_int_normalized, problem%A_exch, &
                quadFit = ( problem%grid%gridType .eq. gridTypeTetrahedron ) .and. .not. ( dlLinStat .eq. 0 .and. dlLinLen .gt. 0 ))
        endif
    endif

    
    end subroutine ComputeExchangeTerm3D


    subroutine passDifferentialOperators(problem)
        implicit none
        type(MicroMagProblem),intent(inout) :: problem     !> The problem data structure
        real(dp) :: eps_criteria
    
        eps_criteria = 1.0e-12
        call displayGUIMessage( 'Creating Exchange matrix start' )
        call create_CSR_matrix(problem%grid%A_exch_load%rows, problem%grid%A_exch_load%cols, problem%grid%A_exch_load%values, problem%grid%A_exch_load%nrows, problem%grid%A_exch_load%ncols, eps_criteria, problem%A_exch)
        call displayGUIMessage( 'Creating Exchange matrix end' )
    
    end subroutine passDifferentialOperators

end module ExchangeOperator
