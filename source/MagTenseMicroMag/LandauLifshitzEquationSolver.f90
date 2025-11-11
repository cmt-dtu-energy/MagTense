! #if USE_MATLAB
! include 'mkl_blas.f90'
! #endif

!include 'mkl_vml.f90'
    module LandauLifshitzSolution
    use ODE_Solvers
    use integrationDataTypes
    use MKL_SPBLAS
    use MKL_DFTI
 !   use MKL_VML
    use BLAS95
#if USE_MATLAB        
    use MagTenseMicroMagIO
#else
    use MagTenseMicroMagPyIO
#endif
    use MicroMagParameters
    use UTIL_CALL
    use UTIL_MICROMAG
    use DemagAuxFunctions
    use DemagFieldGetSolution
    use FortranCuda
    !use, intrinsic :: omp_lib
    use omp_lib
    use IO_GENERAL
    use iso_fortran_env
    use UnstructuredMeshAnalysis
    use DifferentialOperators
    implicit none
       
    
    !>Module variables
    type(MicroMagSolution) :: gb_solution
    type(MicroMagProblem) :: gb_problem
    
    real(DP),dimension(:),allocatable :: crossX,crossY,crossZ   !>Cross product terms
    real(DP),dimension(:),allocatable :: HeffX,HeffY,HeffZ      !>Effective fields
    real(DP),dimension(:),allocatable :: HeffX2,HeffY2,HeffZ2      !>Effective fields
    
    private :: gb_solution,gb_problem,crossX,crossY,crossZ,HeffX,HeffY,HeffZ,HeffX2,HeffY2,HeffZ2
    
    contains
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> @param[inout] prob data structure containing the problem to be solved
    !> @param[inout] sol data structure containing the solution    
    !>-----------------------------------------
    subroutine SolveLandauLifshitzEquation( prob, sol )
    !DEC$ ATTRIBUTES ALIAS:"solvelandaulifshitzequation_" :: SolveLandauLifshitzEquation
    type(MicroMagProblem),intent(inout) :: prob     !> The problem data structure
    type(MicroMagSolution),intent(inout) :: sol     !> The solution data structure    
    type(MicroMagGridInfo) :: gridinfo              !> The grid information structure
    integer :: ntot,i,j,k,ind,nt,nt_Hext,stat       !> total no. of tiles
    procedure(dydt_fct), pointer :: fct             !> Input function pointer for the function to be integrated
    procedure(callback_fct),pointer :: cb_fct       !> Callback function for displaying progress
    real(DP),dimension(:,:,:),allocatable :: M_out        !> Internal buffer for the solution (M) on the form (3*ntot,nt)
    character*(100) :: prog_str 
    
    !Save internal representation of the problem and the solution
    gb_solution = sol
    gb_problem = prob
    
    ntot = gb_problem%grid%nx * gb_problem%grid%ny * gb_problem%grid%nz
    
    !Analyze the mesh, if needed
    if ( gb_problem%grid%gridType .eq. gridTypeUnstructuredPrisms ) then
        if ( gb_problem%passExch .eq. passExchTrue) then
            call displayGUIMessage( 'Passing exchange matrix' )
            call passDifferentialOperators(gb_problem)
        else    
            call CartesianUnstructuredMeshAnalysis(gb_problem%grid%pts, gb_problem%grid%abc, gridinfo)
            call computeDifferentialOperatorsFromMesh_DirectLap(gridinfo, gb_problem%exch_interpn, gb_problem%exch_weight, gb_problem%exch_method, gb_problem%Jfact, gb_problem%A_exch)
            gb_solution%gridinfo = gridinfo
        endif
    endif
      
    call displayGUIMessage( 'Initializing matrices' )
    !Calculate the interaction matrices
    call initializeInteractionMatrices( gb_problem )
    
    !write(prog_str,'(A37, F8.4, A5)') 'Demagnetization tensor memory usage: ', 6*storage_size(gb_problem%Kxx)*ntot/(8*2**30), ' gigabytes'
    !call displayGUIMessage( trim(prog_str) )    
    
    !Copy the demag tensor to CUDA
    if ( gb_problem%useCuda .eq. useCudaTrue ) then
        call displayGUIMessage( 'Copying to CUDA' )
#if USE_CUDA
        !Initialize the Cuda arrays and load the demag tensors into the GPU memory
        if ( ( gb_problem%demag_approximation .eq. DemagApproximationThreshold ) .or. ( gb_problem%demag_approximation .eq. DemagApproximationThresholdFraction ) ) then
           
            !If the matrices are sparse
            call cudaInit_sparse( gb_problem%K_s )        
            
        else
            !if the matrices are dense
            call cudaInit_s( gb_problem%Kxx, gb_problem%Kxy, gb_problem%Kxz, gb_problem%Kyy, gb_problem%Kyz, gb_problem%Kzz )
            
        endif
#else
        call displayGUIMessage( 'MagTense not compiled with CUDA - exiting!' )
        stop
#endif
    endif
   
    allocate( gb_solution%pts(ntot,3) )
    if ( gb_problem%grid%gridType .eq. gridTypeUniform ) then   
        do k=1,gb_problem%grid%nz
            do j=1,gb_problem%grid%ny            
                do i=1,gb_problem%grid%nx
                    ind = i + (j-1) * gb_problem%grid%nx + (k-1) * gb_problem%grid%nx * gb_problem%grid%ny
                    gb_solution%pts(ind,1) = gb_problem%grid%x(i,j,k)
                    gb_solution%pts(ind,2) = gb_problem%grid%y(i,j,k)
                    gb_solution%pts(ind,3) = gb_problem%grid%z(i,j,k)
                enddo
            enddo
        enddo
    endif
    
    if (( gb_problem%grid%gridType .eq. gridTypeTetrahedron ) .or. (gb_problem%grid%gridType .eq. gridTypeUnstructuredPrisms)) then
        do i=1,gb_problem%grid%nx
            gb_solution%pts( i, 1 ) = gb_problem%grid%pts( i, 1 )
            gb_solution%pts( i, 2 ) = gb_problem%grid%pts( i, 2 )
            gb_solution%pts( i, 3 ) = gb_problem%grid%pts( i, 3 )
        enddo
    endif    
    
    
    call displayGUIMessage( 'Initializing solution' )
    !Initialize the solution, i.e. allocate various arrays
    call initializeSolution( gb_problem, gb_solution )
         
    !Set the initial values for m (remember that M is organized such that mx = m(1:ntot), my = m(ntot+1:2*ntot), mz = m(2*ntot+1:3*ntot)
    allocate(gb_solution%t_out(size(gb_problem%t)))
        
    
    call displayGUIMessage( 'Running solution' )
    
    !Do the solution
    fct => dmdt_fct
    cb_fct => displayGUIProgressMessage
    
    gb_solution%HextInd = 1;
    
    !Go through a range of applied fields and find the equilibrium solution for each of them
    !The no. of applied fields to consider
    nt = size( gb_problem%t ) 
        
    if ( gb_problem%solver .eq. MicroMagSolverExplicit ) then
        !Run several different applied fields
        nt_Hext = size(gb_problem%Hext, 1) 
    else if ( gb_problem%solver .eq. MicroMagSolverDynamic ) then
        !Simply do a time evolution as specified in the problem  
        nt_Hext = 1
    endif
    
    allocate(M_out(3*ntot,nt,nt_Hext))   
    allocate(gb_solution%M_out(size(gb_problem%t),ntot,nt_Hext,3))
    !Allocate the arrays for the different fields
    !Only if these are to be returned are they saved at every time step
    if (gb_problem%useReturnHall .eq. useReturnHallTrue) then
        allocate(gb_solution%H_exc(size(gb_problem%t),ntot,nt_Hext,3))
        allocate(gb_solution%H_ext(size(gb_problem%t),ntot,nt_Hext,3))
        allocate(gb_solution%H_dem(size(gb_problem%t),ntot,nt_Hext,3))
        allocate(gb_solution%H_ani(size(gb_problem%t),ntot,nt_Hext,3))  
    else
        allocate(gb_solution%H_exc(1,1,1,3))
        allocate(gb_solution%H_ext(1,1,1,3))
        allocate(gb_solution%H_dem(1,1,1,3))
        allocate(gb_solution%H_ani(1,1,1,3))
        gb_solution%H_exc(:,:,:,:) = 0
        gb_solution%H_ext(:,:,:,:) = 0
        gb_solution%H_dem(:,:,:,:) = 0
        gb_solution%H_ani(:,:,:,:) = 0
    endif

    !loop over the range of applied fields
    do i=1,nt_Hext
        !Applied field
        gb_solution%HextInd = i
        
        if (gb_problem%solver .eq. MicroMagSolverExplicit) then               
            write(prog_str,'(A20, I5.2, A8, I5.2, A6, F6.2, A7)') 'External Field nr.: ', i, ' out of ', nt_Hext, ' i.e. ', real(i)/real(nt_Hext)*100,'% done'
            call displayGUIMessage( trim(prog_str) )
        endif

        call MagTense_ODE( fct, gb_problem%t, gb_problem%m0, gb_solution%t_out, M_out(:,:,i), cb_fct, gb_problem%setTimeDisplay, gb_problem%tol, gb_problem%thres_value, gb_problem%useCVODE, gb_problem%t_conv, gb_problem%conv_tol )  
        
        !The initial state of the next solution is the previous solution result
        gb_problem%m0 = M_out(:,nt,i)
        
        !Store the solution
        gb_solution%M_out(:,:,i,1) =  transpose( M_out(1:ntot,:,i) )
        gb_solution%M_out(:,:,i,2) =  transpose( M_out((ntot+1):2*ntot,:,i) )
        gb_solution%M_out(:,:,i,3) =  transpose( M_out((2*ntot+1):3*ntot,:,i)  )
            
        call StoreHeffComponents ( gb_problem, gb_solution )
            
    enddo

    !clean up
    deallocate(crossX,crossY,crossZ,HeffX,HeffY,HeffZ,HeffX2,HeffY2,HeffZ2, M_out)
    
    !clean up
    if (gb_problem%CV > 0) then
        deallocate(gb_solution%u1, gb_solution%u2, gb_solution%u3, gb_solution%u4, gb_solution%u5, gb_solution%u6)
    endif
    
    !clean-up
    stat = DftiFreeDescriptor(gb_problem%desc_hndl_FFT_M_H)
   
#if USE_CUDA
    if ( gb_problem%useCuda .eqv. useCudaTrue ) then
        call cudaDestroy()
    endif
#endif
    
    !Return the correct state
    sol = gb_solution
    prob = gb_problem
    
    end subroutine SolveLandauLifshitzEquation

    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Defines the function that gives the derivative dmdt that is to be integrated
    !> in the Landau-Lifshitz equation. This function only works on a uniform grid
    !> @param[in] t the time at which the derivative is requested
    !> @param[in] m array size n holding the m_i values corresponding to the time t
    !> @param[inout] dmdt array size n for the derivatives at the time t
    !---------------------------------------------------------------------------    
    subroutine dmdt_fct ( t, m, dmdt )
        !---- inputs/outputs ----
        real(DP),intent(in)              :: t
        real(DP),dimension(:),intent(in) :: m
        real(DP),dimension(:),intent(inout) :: dmdt

        !---- locals ----
        integer :: ntot
        
        !---- persistent counters/accumulators ----
#if USE_TIMING
        real(DP) :: t0, t1
        integer,          save :: call_count = 0
        integer, parameter     :: NPRINT = 200   ! print every N calls (tune as you like)

        real(DP), save :: acc_total    = 0.0_DP
        real(DP), save :: acc_exch     = 0.0_DP
        real(DP), save :: acc_ext      = 0.0_DP
        real(DP), save :: acc_aniso    = 0.0_DP
        real(DP), save :: acc_demag    = 0.0_DP
        real(DP), save :: acc_heff     = 0.0_DP
        real(DP), save :: acc_cross    = 0.0_DP
        real(DP), save :: acc_double   = 0.0_DP
        real(DP), save :: acc_assemble = 0.0_DP
#endif

        !------------------------------------------
        ntot = gb_problem%grid%nx * gb_problem%grid%ny * gb_problem%grid%nz

        if ( .not. allocated(crossX) ) then
            allocate( crossX(ntot), crossY(ntot), crossZ(ntot) )
            allocate( HeffX(ntot), HeffY(ntot), HeffZ(ntot) )
            allocate( HeffX2(ntot), HeffY2(ntot), HeffZ2(ntot) )
        endif
#if USE_TIMING
        t0 = walltime()
        ! Update magnetisation (load from m)
        t1 = walltime()
#endif
        gb_solution%Mx = m(1:ntot)
        gb_solution%My = m(ntot+1:2*ntot)
        gb_solution%Mz = m(2*ntot+1:3*ntot)
        
        ! Exchange term
#if USE_TIMING
        acc_assemble = acc_assemble + (walltime() - t1)

        t1 = walltime()
#endif
        call updateExchangeTerms( gb_problem, gb_solution )

#if USE_TIMING

        acc_exch = acc_exch + (walltime() - t1)

        ! External field
        t1 = walltime()
#endif
        call updateExternalField( gb_problem, gb_solution, t )
#if USE_TIMING

        acc_ext = acc_ext + (walltime() - t1)

        ! Anisotropy term
        t1 = walltime()
    
#endif
        call updateAnisotropy(  gb_problem, gb_solution )
#if USE_TIMING
        acc_aniso = acc_aniso + (walltime() - t1)

        ! Demagnetising field (FMM)
        t1 = walltime()
#endif
        call updateDemagfieldFMM( gb_problem, gb_solution )
        ! call updateDemagfield( gb_problem, gb_solution )
#if USE_TIMING
        acc_demag = acc_demag + (walltime() - t1)

        ! Effective field combine
        t1 = walltime()
#endif
        HeffX = gb_solution%HhX + gb_solution%HjX + gb_solution%HmX + gb_solution%HkX
        HeffY = gb_solution%HhY + gb_solution%HjY + gb_solution%HmY + gb_solution%HkY
        HeffZ = gb_solution%HhZ + gb_solution%HjZ + gb_solution%HmZ + gb_solution%HkZ
#if USE_TIMING
        acc_heff = acc_heff + (walltime() - t1)

        ! Precession term: m x Heff
        t1 = walltime()
#endif
        crossX = -1.0_DP * ( gb_solution%My * HeffZ - gb_solution%Mz * HeffY )
        crossY = -1.0_DP * ( gb_solution%Mz * HeffX - gb_solution%Mx * HeffZ )
        crossZ = -1.0_DP * ( gb_solution%Mx * HeffY - gb_solution%My * HeffX )
#if USE_TIMING
        acc_cross = acc_cross + (walltime() - t1)

        ! Damping term: m x (m x Heff)
        t1 = walltime()
#endif
        HeffX2 = gb_solution%My * crossZ - gb_solution%Mz * crossY
        HeffY2 = gb_solution%Mz * crossX - gb_solution%Mx * crossZ
        HeffZ2 = gb_solution%Mx * crossY - gb_solution%My * crossX
#if USE_TIMING
        acc_double = acc_double + (walltime() - t1)

        ! Assemble dmdt
        t1 = walltime()
#endif
        dmdt(1:ntot)               = -gb_problem%gamma * crossX - alpha(t,gb_problem) * HeffX2
        dmdt(ntot+1:2*ntot)        = -gb_problem%gamma * crossY - alpha(t,gb_problem) * HeffY2
        dmdt(2*ntot+1:3*ntot)      = -gb_problem%gamma * crossZ - alpha(t,gb_problem) * HeffZ2
#if USE_TIMING
        acc_assemble = acc_assemble + (walltime() - t1)

        acc_total = acc_total + (walltime() - t0)
        call_count = call_count + 1
        if (mod(call_count, NPRINT) == 0) then
            ! Compact one-liner per category (seconds over last NPRINT calls)
            write(*,'(A, I0, A, 1X, A, F10.6, 1X, A, F10.6, 1X, A, F10.6, 1X, A, F10.6, 1X, A, F10.6, 1X, A, F10.6, 1X, A, F10.6, 1X, A, F10.6, 1X, A, F10.6)') &
                'dmdt_fct timing (last ', NPRINT, ' calls)  ', &
                'total=',    acc_total,    ' exch=',   acc_exch,   ' ext=',     acc_ext, &
                ' aniso=',   acc_aniso,    ' demag=',  acc_demag,  ' heff=',    acc_heff, &
                ' mxh=',     acc_cross,    ' m×(mxh)=',acc_double, ' asmb=',    acc_assemble

            ! reset accumulators for next block
            acc_total    = 0.0_DP
            acc_exch     = 0.0_DP
            acc_ext      = 0.0_DP
            acc_aniso    = 0.0_DP
            acc_demag    = 0.0_DP
            acc_heff     = 0.0_DP
            acc_cross    = 0.0_DP
            acc_double   = 0.0_DP
            acc_assemble = 0.0_DP
        end if
#endif
    end subroutine dmdt_fct

    real(8) function walltime()
        use omp_lib, only: omp_get_wtime
        walltime = omp_get_wtime()
    end function walltime

    
    !>--------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Calulates the alpha-coefficient for the Landau-Lifshitz equation (can be time dependent)
    !> @param[in] t the time at which to evaluate alpha
    !> @param[in] problem the problem on which the solution is solved
    function alpha( t, problem )
    real(DP) :: alpha
    real(DP),intent(in) :: t
    type(MicroMagProblem),intent(in) :: problem
    
    if ( problem%alpha0 .eq. 0 ) then
        !Interpolate to get the applied field at time t
        call interp1_MagTense( problem%alpha(:,1), problem%alpha(:,2), t, size(problem%alpha(:,1)), alpha )
        
    else
        alpha = problem%alpha0
    endif
    
    end function alpha
    
    
    
    !>-----------------------------------------
    !> @author Rasmus Bjoerk, rabj@dtu.dk, DTU, 2020
    !> @brief
    !> Defines the function that calculates and stores the individual terms
    !> of the effective magnetic field 
    !> @param[in] problem the specification of the problem
    !> @param[in] solution the solution structure where the fields are stored
    !---------------------------------------------------------------------------    
    subroutine StoreHeffComponents ( problem, solution )
    type(MicroMagProblem),intent(in) :: problem
    type(MicroMagSolution),intent(inout) :: solution
    integer :: i,j,nt
    
    i = gb_solution%HextInd       
    nt = size( gb_problem%t ) 
    
    do j=1,nt
        !Calculate the effective field for the values where the magnetization is known
        gb_solution%Mx = gb_solution%M_out(j,:,i,1)
        gb_solution%My = gb_solution%M_out(j,:,i,2)
        gb_solution%Mz = gb_solution%M_out(j,:,i,3)
    
        !Exchange term    
        call updateExchangeTerms( gb_problem, gb_solution )
        !External field
        call updateExternalField( gb_problem, gb_solution, gb_problem%t(j) )
        !Anisotropy term
        call updateAnisotropy(  gb_problem, gb_solution )
        !Demag. field
        !call updateDemagfield( gb_problem, gb_solution )
        call updateDemagfieldFMM( gb_problem, gb_solution )

    
        if (gb_problem%useReturnHall .eq. useReturnHallTrue) then
            !Store the components of the effective field
            gb_solution%H_exc(j,:,i,1) = gb_solution%HjX
            gb_solution%H_exc(j,:,i,2) = gb_solution%HjY
            gb_solution%H_exc(j,:,i,3) = gb_solution%HjZ
        
            gb_solution%H_ext(j,:,i,1) = gb_solution%HhX
            gb_solution%H_ext(j,:,i,2) = gb_solution%HhY
            gb_solution%H_ext(j,:,i,3) = gb_solution%HhZ
        
            gb_solution%H_dem(j,:,i,1) = gb_solution%HmX
            gb_solution%H_dem(j,:,i,2) = gb_solution%HmY
            gb_solution%H_dem(j,:,i,3) = gb_solution%HmZ
        
            gb_solution%H_ani(j,:,i,1) = gb_solution%HkX
            gb_solution%H_ani(j,:,i,2) = gb_solution%HkY
            gb_solution%H_ani(j,:,i,3) = gb_solution%HkZ
        endif
    enddo
    
    end subroutine StoreHeffComponents
    
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Initializes the solution arrays
    !> @param[in] problem the problem from which the solution is found
    !> @param[inout] solution the solution data structure
    subroutine initializeSolution( problem, solution )
    type(MicroMagProblem),intent(in) :: problem
    type(MicroMagSolution),intent(inout) :: solution
    character*(100) :: prog_str 
    
    integer :: ntot, i
    !character(50) :: prog_str
    
    if ( problem%problemMode .eq. ProblemModeNew ) then
        !No. of grid points
        ntot = problem%grid%nx * problem%grid%ny * problem%grid%nz
        !Magnetization
        allocate( solution%Mx(ntot), solution%My(ntot), solution%Mz(ntot) )
        solution%Mx(:) = 0.
        solution%My(:) = 0.
        solution%Mz(:) = 0.
        
        !Magnetization in single variables
        allocate( solution%Mx_s(ntot), solution%My_s(ntot), solution%Mz_s(ntot) )
        solution%Mx_s(:) = 0.
        solution%My_s(:) = 0.
        solution%Mz_s(:) = 0.
        
        !Exchange effective field
        allocate( solution%HjX(ntot), solution%HjY(ntot), solution%HjZ(ntot) )
        solution%HjX(:) = 0.
        solution%HjY(:) = 0.
        solution%HjZ(:) = 0.
        
        !External effective field
        allocate( solution%HhX(ntot), solution%HhY(ntot), solution%HhZ(ntot) )
        solution%HhX(:) = 0.
        solution%HhY(:) = 0.
        solution%HhZ(:) = 0.
        
        !Anisotropy
        allocate( solution%HkX(ntot), solution%HkY(ntot), solution%HkZ(ntot) )
        solution%HkX(:) = 0.
        solution%HkY(:) = 0.
        solution%HkZ(:) = 0.
        
        !Demag field
        allocate( solution%HmX(ntot), solution%HmY(ntot), solution%HmZ(ntot) )
        solution%HmX(:) = 0.
        solution%HmY(:) = 0.
        solution%HmZ(:) = 0.
        
        if ( ( problem%demag_approximation .eq. DemagApproximationFFTThreshold ) .or. ( problem%demag_approximation .eq. DemagApproximationFFTThresholdFraction ) ) then
            !allocate the Fourier Transform of the magnetization
            allocate( solution%Mx_FT(ntot), solution%My_FT(ntot), solution%Mz_FT(ntot) )
            allocate( solution%HmX_c(ntot), solution%HmY_c(ntot), solution%HmZ_c(ntot) )
        endif
        
        
    endif
    
    !If a random noise is present initialize the random vectors
    if (problem%CV > 0) then
        call AddUncertaintyToDemagField( problem, solution)
    endif
    
    
    end subroutine initializeSolution
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Calculates and returns the exchange terms
    !> @param[in] problem, the struct containing the current problem
    !> @param[inout] solution, struct containing the current solution    
    !>-----------------------------------------
    subroutine updateExchangeTerms( problem, solution )
    type(MicroMagProblem),intent(in) :: problem
    type(MicroMagSolution),intent(inout) :: solution
    
    integer :: stat, ntot
    type(MATRIX_DESCR) :: descr
    real(DP) :: alpha, beta
    real(DP), dimension(:), allocatable :: temp
    
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    descr%mode = SPARSE_FILL_MODE_FULL
    descr%diag = SPARSE_DIAG_NON_UNIT
    
    alpha = -2.! * solution%Jfact
    !alpha = -2. * problem%A0 / ( mu0 * 8.0e5 )
    beta = 0.
    
    ntot = problem%grid%nx * problem%grid%ny * problem%grid%nz
    allocate(temp(ntot))
    !Effective field in the X-direction. Note that the scalar alpha is multiplied on from the left, such that
    !y = alpha * (A_exch * Mx )
    stat = mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%A_exch, descr, solution%Mx, beta, temp )
    !If we have an unstructued grid, the exchange matrix already includes Jfact
    if ( gb_problem%grid%gridType .eq. gridTypeUniform ) then  
        solution%HjX = temp * problem%Jfact
    else
        solution%HjX = temp
    endif
    
    !Effective field in the Y-direction
    stat = mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%A_exch, descr, solution%My, beta, temp )
    if ( gb_problem%grid%gridType .eq. gridTypeUniform ) then 
        solution%HjY = temp * problem%Jfact
    else
        solution%HjY = temp
    endif
    
    !Effective field in the Z-direction
    stat = mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%A_exch, descr, solution%Mz, beta, temp )
    if ( gb_problem%grid%gridType .eq. gridTypeUniform ) then 
        solution%HjZ = temp * problem%Jfact
    else
        solution%HjZ = temp
    endif
    
    deallocate(temp)
    
    end subroutine updateExchangeTerms

    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Calculates and returns the external field
    !> @param[in] problem, the struct containing the current problem
    !> @param[inout] solution, struct containing the current solution    
    !> @param[in] t the current time
    !>-----------------------------------------
    subroutine updateExternalField( problem, solution, t )
    type(MicroMagProblem),intent(in) :: problem         !> Problem data structure    
    type(MicroMagSolution),intent(inout) :: solution    !> Solution data structure
    real(DP),intent(in) :: t
    
    real(DP) :: HextX,HextY,HextZ
    
    if ( problem%solver .eq. MicroMagSolverExplicit ) then
         !Assume the field to be constant in time (we are finding the equilibrium solution at a given applied field)
        solution%HhX = -problem%Hext(solution%HextInd,2)
        solution%HhY = -problem%Hext(solution%HextInd,3)
        solution%HhZ = -problem%Hext(solution%HextInd,4)

    elseif ( problem%solver .eq. MicroMagSolverDynamic ) then
        
        !Interpolate to get the applied field at time t
        call interp1_MagTense( problem%Hext(:,1), problem%Hext(:,2), t, size(problem%Hext(:,1)), HextX )
        call interp1_MagTense( problem%Hext(:,1), problem%Hext(:,3), t, size(problem%Hext(:,1)), HextY )
        call interp1_MagTense( problem%Hext(:,1), problem%Hext(:,4), t, size(problem%Hext(:,1)), HextZ )
        
        solution%HhX = -HextX
        solution%HhY = -HextY
        solution%HhZ = -HextZ
        
    elseif ( problem%solver .eq. MicroMagSolverImplicit ) then
        !not implemented yet
    endif
    
    
    end subroutine updateExternalField
    
    !>-----------------------------------------
    !> @author Rasmus Bjørk, rabj@dtu.dk, DTU, 2025
    !> @brief
    !> Calculates and returns the effective field from the anisotropy for a general crystal anisotropy
    !> @param[in] problem, the struct containing the current problem
    !> @param[inout] solution, struct containing the current solution        
    !>-----------------------------------------
    subroutine updateAnisotropy( problem, solution)
    type(MicroMagProblem),intent(in) :: problem         !> Problem data structure    
    type(MicroMagSolution),intent(inout) :: solution    !> Solution data structure
    
    !real :: prefact                                    !> Multiplicative scalar factor
    type(MATRIX_DESCR) :: descr                         !>descriptor for the sparse matrix-vector multiplication
    real(DP),dimension(:),allocatable   :: Mx_rot, My_rot, Mz_rot, Hkx_rot, Hky_rot, Hkz_rot
    
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    descr%mode = SPARSE_FILL_MODE_FULL
    descr%diag = SPARSE_DIAG_NON_UNIT
       
    ! 3x3 matrix specifying the local coordinate system
    !                             [v1_x v2_x v3_x]
    !problem%CrystalAxis(i,:,:) = [v1_y v2_y v3_y]
    !                             [v1_z v2_z v3_z]
        
    ! The anisotropy matrix is given in the local coordinate system, i.e. the crystal axis,
    ! and is specified according to Eq. (2) in https://doi.org/10.1088/1361-665X/aafff8 but generalized to
    ! independent coordinates:
    !  problem%Kfact_arr(i,:,:) = [alpha1_x   alpha1_y   alpha1_z  ]
    !                             [alpha11_x  alpha11_x  alpha11_x ]
    !                             [alpha12_x  alpha12_y  alpha12_z ]   
    !                             [alpha111_x alpha111_y alpha111_z]   
    !                             [alpha112_x alpha112_y alpha112_z]   
    !                             [alpha123   0          0         ]
    !
    ! In the above matrix notation, a uniaxial anisotropy in the z-direction is given by:
    ! Uniaxial in z = [0 0 Kfact]
    !                 [0 0 0 ]
    !                 [0 0 0 ]
    !                 [0 0 0 ]
    !                 [0 0 0 ]
    !                 [0 0 0 ]
    !
    ! In the above matrix notation, a cubic anisotropy is given by (http://wpage.unina.it/mdaquino/PhD_thesis/main/node13.html):
    ! Cubic in cartesian = [0 0 0]
    !                      [0 0 0]
    !                      [Kfact1 Kfact1 Kfact1]  
    !                      [0 0 0]
    !                      [0 0 0]
    !                      [Kfact2 0 0]
    
    allocate(Mx_rot(size(solution%Mx)), My_rot(size(solution%My)), Mz_rot(size(solution%Mz)))
    Mx_rot = problem%CrystalAxis(:,1,1)*solution%Mx(:) + problem%CrystalAxis(:,1,2)*solution%My(:) + problem%CrystalAxis(:,1,3)*solution%Mz(:)
    My_rot = problem%CrystalAxis(:,2,1)*solution%Mx(:) + problem%CrystalAxis(:,2,2)*solution%My(:) + problem%CrystalAxis(:,2,3)*solution%Mz(:)
    Mz_rot = problem%CrystalAxis(:,3,1)*solution%Mx(:) + problem%CrystalAxis(:,3,2)*solution%My(:) + problem%CrystalAxis(:,3,3)*solution%Mz(:)
        
    allocate(Hkx_rot(size(Mx_rot)), Hky_rot(size(My_rot)), Hkz_rot(size(Mz_rot)))
    Hkx_rot = -2*Mx_rot(:)*(problem%Kfact_arr(:,1,1) + 2*problem%Kfact_arr(:,2,1)*Mx_rot(:)**2 + problem%Kfact_arr(:,3,3)*My_rot(:)**2 + problem%Kfact_arr(:,3,2)*Mz_rot(:)**2 + 3*problem%Kfact_arr(:,4,1)*Mx_rot(:)**4 + problem%Kfact_arr(:,5,2)*My_rot(:)**4 + problem%Kfact_arr(:,5,3)*Mz_rot(:)**4 + 2*problem%Kfact_arr(:,5,1)*Mx_rot(:)**2*My_rot(:)**2 + 2*problem%Kfact_arr(:,5,1)*Mx_rot(:)**2*Mz_rot(:)**2 + problem%Kfact_arr(:,6,1)*My_rot(:)**2*Mz_rot(:)**2 )
    Hky_rot = -2*My_rot(:)*(problem%Kfact_arr(:,1,2) + 2*problem%Kfact_arr(:,2,2)*My_rot(:)**2 + problem%Kfact_arr(:,3,3)*Mx_rot(:)**2 + problem%Kfact_arr(:,3,1)*Mz_rot(:)**2 + 3*problem%Kfact_arr(:,4,2)*My_rot(:)**4 + problem%Kfact_arr(:,5,1)*Mx_rot(:)**4 + problem%Kfact_arr(:,5,3)*Mz_rot(:)**4 + 2*problem%Kfact_arr(:,5,2)*Mx_rot(:)**2*My_rot(:)**2 + 2*problem%Kfact_arr(:,5,2)*My_rot(:)**2*Mz_rot(:)**2 + problem%Kfact_arr(:,6,1)*Mx_rot(:)**2*Mz_rot(:)**2 )
    Hkz_rot = -2*Mz_rot(:)*(problem%Kfact_arr(:,1,3) + 2*problem%Kfact_arr(:,2,3)*Mz_rot(:)**2 + problem%Kfact_arr(:,3,2)*Mx_rot(:)**2 + problem%Kfact_arr(:,3,1)*My_rot(:)**2 + 3*problem%Kfact_arr(:,4,3)*Mz_rot(:)**4 + problem%Kfact_arr(:,5,1)*Mx_rot(:)**4 + problem%Kfact_arr(:,5,2)*My_rot(:)**4 + 2*problem%Kfact_arr(:,5,3)*Mx_rot(:)**2*Mz_rot(:)**2 + 2*problem%Kfact_arr(:,5,3)*My_rot(:)**2*Mz_rot(:)**2 + problem%Kfact_arr(:,6,1)*Mx_rot(:)**2*My_rot(:)**2 )
    
    solution%Hkx(:) = problem%CrystalAxis(:,1,1)*Hkx_rot(:) + problem%CrystalAxis(:,2,1)*Hky_rot(:) + problem%CrystalAxis(:,3,1)*Hkz_rot(:)
    solution%Hky(:) = problem%CrystalAxis(:,1,2)*Hkx_rot(:) + problem%CrystalAxis(:,2,2)*Hky_rot(:) + problem%CrystalAxis(:,3,2)*Hkz_rot(:)
    solution%Hkz(:) = problem%CrystalAxis(:,1,3)*Hkx_rot(:) + problem%CrystalAxis(:,2,3)*Hky_rot(:) + problem%CrystalAxis(:,3,3)*Hkz_rot(:)
    
    deallocate(Mx_rot, My_rot, Mz_rot, Hkx_rot, Hky_rot, Hkz_rot)

    end subroutine updateAnisotropy    




   !===============================================================
! FMM-backed demag field update (sources->sources)
! Uses FMM3D Laplace kernel with dipoles only:
!   u(x) = - v_j · ∇(1/|x - x_j|)
!   H(x) = -∇u(x) / (4π)
!===============================================================
! subroutine updateDemagfieldFMM(problem, solution)
!       implicit none
!   !------------------ Arguments ------------------
!   type(MicroMagProblem),  intent(inout) :: problem
!   type(MicroMagSolution), intent(inout) :: solution
! #if USE_FMM3D
!   !------------------ Locals ---------------------
!   integer :: ntot, i
!   real(DP) :: fourpi
!   ! FMM arrays (Laplace dipole, source->source)
!   integer :: nd, ier
!   double precision , allocatable :: source(:,:), dipvec(:,:,:)
!   double precision , allocatable :: pot(:,:), grad(:,:,:)
!   double precision :: vol_i
!   double precision :: mx, my, mz
!   double precision :: eps_fmm
!   double precision :: avg
!   !------------------------------------------------
!   interface
!     subroutine lfmm3d_s_d_g_vec(nd,eps,nsource,source,dipvec,pot,grad,ier)
!       implicit none
!       integer, intent(in) :: nd, nsource
!       double precision, intent(in) :: eps   
!       double precision :: source(3,nsource)
!       double precision :: dipvec(nd,3,nsource)
!       double precision :: pot(nd,nsource)
!       double precision :: grad(nd,3,nsource)
!       integer :: ier
!     end subroutine lfmm3d_s_d_g_vec
!   end interface
!   ! Dummy charge (all zeros; we only use dipoles)
!   !------------------ Setup ----------------------
!   fourpi = 12.566370614359172D0
!   ! Pick a tolerance. If you have a field in problem with a natural tol, use that.
!   eps_fmm = 1.0D-6

!   ! Number of cells
!   ntot = problem%grid%nx * problem%grid%ny * problem%grid%nz

!   ! Ensure single->double copies of M exist (you already do this in your current routine)
!   solution%Mx_s = real(solution%Mx, SP)
!   solution%My_s = real(solution%My, SP)
!   solution%Mz_s = real(solution%Mz, SP)

!   !------------------ Allocate FMM work arrays ------------------
!   allocate(source(3, ntot))
!   allocate(dipvec(1, 3, ntot))
!   allocate(pot(1, ntot))
!   allocate(grad(1, 3, ntot))

!   !print *, 'FMM3D: Calculating demag field with FMM for', ntot, 'cells.'

!   !------------------ Pack sources and dipoles ------------------
!   ! Assumptions (adjust if your field names differ):
!   !   - problem%grid%pts(ntot,3) holds cell-centre coordinates (x,y,z)
!   !   - problem%grid%abc(ntot,3) holds cell edge lengths (dx,dy,dz) per cell
!   !eps_fmm = 1e6



!   do i = 1, ntot
!     ! positions -> FMM (3, N)
!     source(1,i) = real(problem%grid%pts(i,1), DP)
!     source(2,i) = real(problem%grid%pts(i,2), DP)
!     source(3,i) = real(problem%grid%pts(i,3), DP)


!     vol_i = real(problem%grid%dx, DP) * &
!             real(problem%grid%dy, DP) * &
!             real(problem%grid%dz, DP)


!     ! dipole moment m = M * ΔV  (units consistent with your M)
!     mx = real(solution%Mx_s(i), DP) * vol_i * problem%Ms(i)
!     my = real(solution%My_s(i), DP) * vol_i * problem%Ms(i)
!     mz = real(solution%Mz_s(i), DP) * vol_i * problem%Ms(i)


!     dipvec(1,1,i) = mx
!     dipvec(1,2,i) = my
!     dipvec(1,3,i) = mz

!   end do


! eps_fmm = 1.0D-6


!   !------------------ Call FMM (sources->sources) ------------------
!   nd   = 1
  
!   !print *, "source = ", source
!   !print*, "dipvec = ", dipvec
  
!   !print *, " Calling FMM3D backend"
!   call lfmm3d_s_d_g_vec( nd, eps_fmm, ntot, &
!        source, dipvec, pot, grad, ier )
! !print *, " FMM3D backend returned"

!   !grad = grad * avg

!   if (ier /= 0) then
!     write(*,*) 'FMM3D returned error code in updateDemagfieldFMM: ier =', ier
!     stop " FMM3D error"
!   end if    

!   !------------------ Map grad -> H (and to single) ------------------
!   ! For Laplace kernel used by FMM3D examples:
!   !   u(x) = - m · ∇(1/r),  grad = ∇u
!   !   H = -grad / (4π)
!   ! NOTE: If your existing matrix-based routes embed an extra scaling
!   !       via problem%Mfact, consider multiplying by that here for
!   !       exact equivalence:
!   !       H *= problem%Mfact
!   do i = 1, ntot
!     solution%HmX(i) = real( grad(1,1,i) / fourpi, SP )
!     solution%HmY(i) = real( grad(1,2,i) / fourpi, SP )
!     solution%HmZ(i) = real( grad(1,3,i) / fourpi, SP )
!   end do


!   !--------------- account for self field correction --------------
!   call add_self_correction(problem, solution)
!   !-----------------------------------------------------------------
!   !--------------- account for neighbour correction --------------
!     call add_neighbour_correction(problem, solution)
!   !-----------------------------------------------------------------


!   !print *, " finished FMM demag field calculation"

!   !------------------ Cleanup ------------------
!   deallocate(source, dipvec, pot, grad)




! #else
!   integer :: ntot
!   ntot = problem%grid%nx * problem%grid%ny * problem%grid%nz
!   solution%HmX(1:ntot) = 0.0_SP
!   solution%HmY(1:ntot) = 0.0_SP
!   solution%HmZ(1:ntot) = 0.0_SP
!   print *, 'WARNING: updateDemagfieldFMM called but built without FMM3D; returning zero field.'
! #endif


!     if (problem%CV > 0) then
!         !print *, 'Adding uncertainty to demag field with CV = ', problem%CV
!         solution%HmX = solution%HmX + solution%HmX*problem%CV*sqrt(-2d0*log(solution%u1))*cos(2d0*pi*solution%u2)
!         solution%HmY = solution%HmY + solution%HmY*problem%CV*sqrt(-2d0*log(solution%u3))*cos(2d0*pi*solution%u4)
!         solution%HmZ = solution%HmZ + solution%HmZ*problem%CV*sqrt(-2d0*log(solution%u5))*cos(2d0*pi*solution%u6)
!     endif

! end subroutine updateDemagfieldFMM

subroutine updateDemagfieldFMM(problem, solution)
  implicit none
  !------------------ Arguments ------------------
  type(MicroMagProblem),  intent(inout) :: problem
  type(MicroMagSolution), intent(inout) :: solution

#if USE_FMM3D
  !------------------ Timing ---------------------
#if USE_TIMING
  integer,          save :: call_count_fmm = 0
  integer, parameter     :: NPRINT = 200    ! keep same cadence as dmdt_fct
  real(DP), save :: acc_total_fmm = 0.0_DP
  real(DP), save :: acc_cast      = 0.0_DP
  real(DP), save :: acc_alloc     = 0.0_DP
  real(DP), save :: acc_pack      = 0.0_DP
  real(DP), save :: acc_fmm       = 0.0_DP
  real(DP), save :: acc_map       = 0.0_DP
  real(DP), save :: acc_self      = 0.0_DP
  real(DP), save :: acc_neigh     = 0.0_DP
  real(DP), save :: acc_cleanup   = 0.0_DP
  real(DP), save :: acc_noise     = 0.0_DP

  real(DP) :: t0, t1
#endif
  !------------------ Locals ---------------------
  integer :: ntot, i
  real(DP) :: fourpi
  integer :: nd, ier
  double precision , allocatable :: source(:,:), dipvec(:,:,:)
  double precision , allocatable :: pot(:,:), grad(:,:,:)
  double precision :: vol_i
  double precision :: mx, my, mz
  double precision :: eps_fmm
  !------------------------------------------------
  interface
    subroutine lfmm3d_s_d_g_vec(nd,eps,nsource,source,dipvec,pot,grad,ier)
      implicit none
      integer, intent(in) :: nd, nsource
      double precision, intent(in) :: eps
      double precision :: source(3,nsource)
      double precision :: dipvec(nd,3,nsource)
      double precision :: pot(nd,nsource)
      double precision :: grad(nd,3,nsource)
      integer :: ier
    end subroutine lfmm3d_s_d_g_vec
  end interface

#if USE_TIMING
  t0 = walltime()
#endif
  fourpi  = 12.566370614359172D0
  eps_fmm = 1.0D-4

  ntot = problem%grid%nx * problem%grid%ny * problem%grid%nz

  ! Cast M (DP->SP) if your fields expect it

#if USE_TIMING
  t1 = walltime()
#endif
  solution%Mx_s = real(solution%Mx, SP)
  solution%My_s = real(solution%My, SP)
  solution%Mz_s = real(solution%Mz, SP)
#if USE_TIMING
  acc_cast = acc_cast + (walltime() - t1)

  !------------------ Allocate FMM work arrays ------------------
  t1 = walltime()
#endif
  allocate(source(3, ntot))
  allocate(dipvec(1, 3, ntot))
  allocate(pot(1, ntot))
  allocate(grad(1, 3, ntot))
#if USE_TIMING
  acc_alloc = acc_alloc + (walltime() - t1)

  !------------------ Pack sources and dipoles ------------------
  t1 = walltime()
#endif
  do i = 1, ntot
    source(1,i) = real(problem%grid%pts(i,1), DP)
    source(2,i) = real(problem%grid%pts(i,2), DP)
    source(3,i) = real(problem%grid%pts(i,3), DP)

    vol_i = real(problem%grid%dx, DP) * &
            real(problem%grid%dy, DP) * &
            real(problem%grid%dz, DP)

    mx = real(solution%Mx_s(i), DP) * vol_i * problem%Ms(i)
    my = real(solution%My_s(i), DP) * vol_i * problem%Ms(i)
    mz = real(solution%Mz_s(i), DP) * vol_i * problem%Ms(i)

    dipvec(1,1,i) = mx
    dipvec(1,2,i) = my
    dipvec(1,3,i) = mz
  end do
#if USE_TIMING
  acc_pack = acc_pack + (walltime() - t1)

  t1 = walltime()
#endif
  !------------------ Call FMM (sources->sources) ------------------
  nd = 1
  call lfmm3d_s_d_g_vec( nd, eps_fmm, ntot, source, dipvec, pot, grad, ier )
  
  if (ier /= 0) then
    write(*,*) 'FMM3D returned error code in updateDemagfieldFMM: ier =', ier
    stop " FMM3D error"
  end if
#if USE_TIMING
  acc_fmm = acc_fmm + (walltime() - t1)
  !------------------ Map grad -> H (and to single) ----------------
  t1 = walltime()
#endif
  do i = 1, ntot
    solution%HmX(i) = real( grad(1,1,i) / fourpi, SP )
    solution%HmY(i) = real( grad(1,2,i) / fourpi, SP )
    solution%HmZ(i) = real( grad(1,3,i) / fourpi, SP )
  end do
#if USE_TIMING
  acc_map = acc_map + (walltime() - t1)

  !--------------- Self & neighbour corrections -------------------
  t1 = walltime()
#endif
  call add_self_correction(problem, solution)
#if USE_TIMING
  acc_self = acc_self + (walltime() - t1)

  t1 = walltime()
#endif
  call add_neighbour_correction(problem, solution)
#if USE_TIMING
  acc_neigh = acc_neigh + (walltime() - t1)
#endif

  !------------------ Cleanup -------------------------------------
#if USE_TIMING
  t1 = walltime()
#endif
  deallocate(source, dipvec, pot, grad)
#if USE_TIMING
  acc_cleanup = acc_cleanup + (walltime() - t1)
#endif

#else
  integer :: ntot
  ntot = problem%grid%nx * problem%grid%ny * problem%grid%nz
  solution%HmX(1:ntot) = 0.0_SP
  solution%HmY(1:ntot) = 0.0_SP
  solution%HmZ(1:ntot) = 0.0_SP
  print *, 'WARNING: updateDemagfieldFMM called but built without FMM3D; returning zero field.'
#endif

  !------------------ Optional field noise (CV) --------------------
#if USE_TIMING
  t1 = walltime()
#endif
  if (problem%CV > 0) then
      solution%HmX = solution%HmX + solution%HmX*problem%CV*sqrt(-2d0*log(solution%u1))*cos(2d0*pi*solution%u2)
      solution%HmY = solution%HmY + solution%HmY*problem%CV*sqrt(-2d0*log(solution%u3))*cos(2d0*pi*solution%u4)
      solution%HmZ = solution%HmZ + solution%HmZ*problem%CV*sqrt(-2d0*log(solution%u5))*cos(2d0*pi*solution%u6)
  end if
#if USE_TIMING
  acc_noise = acc_noise + (walltime() - t1)

  !------------------ Final accounting/print -----------------------
  acc_total_fmm = acc_total_fmm + (walltime() - t0)
  call_count_fmm = call_count_fmm + 1

  if (mod(call_count_fmm, NPRINT) == 0) then
      write(*,'(A,I0,A,1X, A,F10.6,1X,A,F10.6,1X,A,F10.6,1X,A,F10.6,1X,A,F10.6,1X,A,F10.6,1X,A,F10.6,1X,A,F10.6,1X,A,F10.6,1X,A,F10.6)') &
           'updateDemagfieldFMM timing (last ', NPRINT, ' calls) ', &
           'total=',  acc_total_fmm, ' cast=',   acc_cast,    ' alloc=',  acc_alloc,  ' pack=',  acc_pack, &
           ' fmm=',   acc_fmm,       ' map=',    acc_map,     ' self=',   acc_self,   ' neigh=', acc_neigh, &
           ' free=',  acc_cleanup,   ' noise=',  acc_noise

      acc_total_fmm = 0.0_DP
      acc_cast      = 0.0_DP
      acc_alloc     = 0.0_DP
      acc_pack      = 0.0_DP
      acc_fmm       = 0.0_DP
      acc_map       = 0.0_DP
      acc_self      = 0.0_DP
      acc_neigh     = 0.0_DP
      acc_cleanup   = 0.0_DP
      acc_noise     = 0.0_DP
  end if
#endif

end subroutine updateDemagfieldFMM


    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Calculates and returns the effective demag field
    !> @param[in] problem, the struct containing the current problem
    !> @param[inout] solution, struct containing the current solution        
    !>-----------------------------------------
    subroutine updateDemagfield( problem, solution)
    type(MicroMagProblem),intent(inout) :: problem         !> Problem data structure    
    type(MicroMagSolution),intent(inout) :: solution    !> Solution data structure
    integer :: stat,ntot,i
    type(matrix_descr) :: descr
    real(SP) :: pref,alpha,beta
    complex(kind=4) :: alpha_c, beta_c
    real(SP), dimension(:), allocatable :: temp
    character*(100) :: prog_str 
    
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    descr%mode = SPARSE_FILL_MODE_FULL
    descr%diag = SPARSE_DIAG_NON_UNIT
    
    ntot = problem%grid%nx * problem%grid%ny * problem%grid%nz
    allocate(temp(ntot))
    
    ! Convert the magnetization to single before the demag calculation
    solution%Mx_s = real(solution%Mx, SP)
    solution%My_s = real(solution%My, SP)
    solution%Mz_s = real(solution%Mz, SP)
    
    
    if ( ( problem%demag_approximation .eq. DemagApproximationThreshold ) .or. ( problem%demag_approximation .eq. DemagApproximationThresholdFraction ) ) then
        if ( problem%useCuda .eq. useCudaFalse ) then
            !Do the matrix multiplications using sparse matrices
            alpha = -1.0
            beta = 0.
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(1)%A, descr, solution%Mx_s, beta, temp )
            beta = 1.0
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(2)%A, descr, solution%My_s, beta, temp )
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(3)%A, descr, solution%Mz_s, beta, temp )
        
            solution%HmX = temp * ( problem%Mfact )
            !ntot = problem%grid%nx * problem%grid%ny * problem%grid%nz
            !call vsmul( ntot, solution%HmX, -problem%Mfact, solution%HmX )
            
            beta = 0.
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(2)%A, descr, solution%Mx_s, beta, temp )
            beta = 1.0
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(4)%A, descr, solution%My_s, beta, temp )
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(5)%A, descr, solution%Mz_s, beta, temp )
        
            solution%HmY = temp * ( problem%Mfact )
            !call vsmul( ntot, solution%HmY, -problem%Mfact, solution%HmY )
          
            beta = 0.
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(3)%A, descr, solution%Mx_s, beta, temp )
            beta = 1.0
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(5)%A, descr, solution%My_s, beta, temp )
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(6)%A, descr, solution%Mz_s, beta, temp )
        
            solution%HmZ = temp * ( problem%Mfact )
            !call vsmul( ntot, solution%HmZ, -problem%Mfact, solution%HmZ )
        else
#if USE_CUDA
            !Do the sparse matrix multiplication using CUDA
            pref = sngl(-1 )!* problem%Mfact)                                
            call cudaMatrVecMult_sparse( solution%Mx_s, solution%My_s, solution%Mz_s, solution%HmX, solution%HmY, solution%HmZ, pref )
            temp = solution%HmX * problem%Mfact
            solution%HmX = temp
            temp = solution%HmY * problem%Mfact
            solution%HmY = temp
            temp = solution%HmZ * problem%Mfact
            solution%HmZ = temp
#endif
        endif
        
    elseif ( ( problem%demag_approximation .eq. DemagApproximationFFTThreshold ) .or. ( problem%demag_approximation .eq. DemagApproximationFFTThresholdFraction ) ) then
        
        if ( problem%useCuda .eq. useCudaFalse ) then
            !fourier transform Mx, My and Mz
            ntot = problem%grid%nx * problem%grid%ny * problem%grid%nz
                
            !Convert to complex format
            do i=1,ntot
                solution%Mx_FT(i) = cmplx( solution%Mx_s(i), 0. )
                solution%My_FT(i) = cmplx( solution%My_s(i), 0. )
                solution%Mz_FT(i) = cmplx( solution%Mz_s(i), 0. )
            enddo
        
            stat = DftiComputeForward( problem%desc_hndl_FFT_M_H, solution%Mx_FT )
            !normalization
            solution%Mx_FT = solution%Mx_FT / ntot
                
            stat = DftiComputeForward( problem%desc_hndl_FFT_M_H, solution%My_FT )
            !normalization
            solution%My_FT = solution%My_FT / ntot
        
            stat = DftiComputeForward( problem%desc_hndl_FFT_M_H, solution%Mz_FT )
            !normalization
            solution%Mz_FT = solution%Mz_FT / ntot
            
            !sparse matrix multiplication with the demag matrices in fourier space...                
            ! use problem%K_s(1..6) with cuda or MKL to do the sparse matrix-vector product with the FFT(M) and subsequently the IFT on the whole thing to get H
        
            !First Hx = Kxx * Mx + Kxy * My + Kxz * Mz
        
            alpha_c = cmplx(1.,0)
            beta_c = cmplx(0.,0.)
            stat = mkl_sparse_c_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha_c, problem%K_s_c(1)%A, descr, solution%Mx_FT, beta_c, solution%HmX_c )
            beta_c = cmplx(1.0,0.)
            stat = mkl_sparse_c_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha_c, problem%K_s_c(2)%A, descr, solution%My_FT, beta_c, solution%HmX_c )
            stat = mkl_sparse_c_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha_c, problem%K_s_c(3)%A, descr, solution%Mz_FT, beta_c, solution%HmX_c )
        
        
            !Fourier transform backwards to get the field
            stat = DftiComputeBackward( problem%desc_hndl_FFT_M_H, solution%HmX_C )
        
            !Get the field
            solution%HmX = -problem%Mfact * real(solution%HmX_c)
            !call vsmul( ntot, real(solution%HmX_c), -problem%Mfact, solution%HmX )
        
        
            !Second Hy = Kyx * Mx + Kyy * My + Kyz * Mz        
            beta_c = cmplx(0.,0.)
            stat = mkl_sparse_c_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha_c, problem%K_s_c(2)%A, descr, solution%Mx_FT, beta_c, solution%HmY_c )
            beta_c = cmplx(1.0,0.)
            stat = mkl_sparse_c_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha_c, problem%K_s_c(4)%A, descr, solution%My_FT, beta_c, solution%HmY_c )
            stat = mkl_sparse_c_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha_c, problem%K_s_c(5)%A, descr, solution%Mz_FT, beta_c, solution%HmY_c )
        
            !Fourier transform backwards to get the field
            stat = DftiComputeBackward( problem%desc_hndl_FFT_M_H, solution%HmY_c )
        
            !Get the field        
            solution%HmY = -problem%Mfact * real(solution%HmY_c)
            !call vsmul( ntot, real(solution%HmY_c), -problem%Mfact, solution%HmY )
        
        
            !Third Hz = Kzx * Mx + Kzy * My + Kzz * Mz        
            beta_c = cmplx(0.,0.)
            stat = mkl_sparse_c_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha_c, problem%K_s_c(3)%A, descr, solution%Mx_FT, beta_c, solution%HmZ_c )
            beta_c = cmplx(1.0,0.)
            stat = mkl_sparse_c_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha_c, problem%K_s_c(5)%A, descr, solution%My_FT, beta_c, solution%HmZ_c )
            stat = mkl_sparse_c_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha_c, problem%K_s_c(6)%A, descr, solution%Mz_FT, beta_c, solution%HmZ_c )
        
            !Fourier transform backwards to get the field
            stat = DftiComputeBackward( problem%desc_hndl_FFT_M_H, solution%HmZ_c )
            !finally, get the field out        
            solution%HmZ = -problem%Mfact * real(solution%HmZ_c)
            !call vsmul( ntot, real(solution%HmZ_c), -problem%Mfact, solution%HmZ )
        
        else
            !!No CUDA support for this part yet
            call displayGUIMessage( 'CUDA not currently supported with FFT - exiting!' )
            stop
        endif
        
    else        
        !Default way of doing the problem
        if ( problem%useCuda .eq. useCudaFalse ) then
            !Needs to be checked for proper matrix calculation (Kxx is an n x n matrix while Mx should be n x 1 column vector and the result an n x 1 column vector)
            !Note that the demag tensor is symmetric such that Kxy = Kyx and we only store what is needed.
            !solution%HmX = - problem%Mfact * ( matmul( problem%Kxx, solution%Mx ) + matmul( problem%Kxy, solution%My ) + matmul( problem%Kxz, solution%Mz ) )
            !solution%HmY = - problem%Mfact * ( matmul( problem%Kxy, solution%Mx ) + matmul( problem%Kyy, solution%My ) + matmul( problem%Kyz, solution%Mz ) )
            !solution%HmZ = - problem%Mfact * ( matmul( problem%Kxz, solution%Mx ) + matmul( problem%Kyz, solution%My ) + matmul( problem%Kzz, solution%Mz ) )
            
            alpha = -1.! * problem%Mfact
            beta = 0.0
            !Hmx = Kxx * Mx
            call gemv( problem%Kxx, solution%Mx_s, solution%HmX, alpha, beta )
                                   
            beta = 1.0
            !Hmx = Hmx + Kxy * My
            call gemv( problem%Kxy, solution%My_s, solution%HmX, alpha, beta )
            
            !Hmx = Hmx + Kxz * Mz
            call gemv( problem%Kxz, solution%Mz_s, solution%HmX, alpha, beta )
            
            beta = 0.0
            !HmY = Kyx * Mx
            call gemv( problem%Kxy, solution%Mx_s, solution%HmY, alpha, beta )
            
            beta = 1.0
            !HmY = HmY + Kyy * My
            call gemv( problem%Kyy, solution%My_s, solution%HmY, alpha, beta )
            
            !Hmy = Hmy + Kyz * Mz
            call gemv( problem%Kyz, solution%Mz_s, solution%HmY, alpha, beta )
            
            beta = 0.0
            !HmZ = Kzx * Mx
            call gemv( problem%Kxz, solution%Mx_s, solution%HmZ, alpha, beta )
            
            beta = 1.0
            !HmZ = HmZ + Kzy * My
            call gemv( problem%Kyz, solution%My_s, solution%HmZ, alpha, beta )
            
            !HmZ = HmZ + Kzz * Mz
            call gemv( problem%Kzz, solution%Mz_s, solution%HmZ, alpha, beta )
            
            temp = solution%HmX * problem%Mfact
            solution%HmX = temp
            temp = solution%HmY * problem%Mfact
            solution%HmY = temp
            temp = solution%HmZ * problem%Mfact
            solution%HmZ = temp
            
        else
            pref = sngl(-1)! * problem%Mfact)
#if USE_CUDA
            call cudaMatrVecMult( solution%Mx_s, solution%My_s, solution%Mz_s, solution%HmX, solution%HmY, solution%HmZ, pref )
            temp = solution%HmX * problem%Mfact
            solution%HmX = temp
            temp = solution%HmY * problem%Mfact
            solution%HmY = temp
            temp = solution%HmZ * problem%Mfact
            solution%HmZ = temp
#endif
        endif 
    endif
        
       
    deallocate(temp)
    
    if (problem%CV > 0) then
        solution%HmX = solution%HmX + solution%HmX*problem%CV*sqrt(-2d0*log(solution%u1))*cos(2d0*pi*solution%u2)
        solution%HmY = solution%HmY + solution%HmY*problem%CV*sqrt(-2d0*log(solution%u3))*cos(2d0*pi*solution%u4)
        solution%HmZ = solution%HmZ + solution%HmZ*problem%CV*sqrt(-2d0*log(solution%u5))*cos(2d0*pi*solution%u6)
    endif
    
    end subroutine updateDemagfield
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Initializes the interaction matrices
    !> @param[inout] problem struct containing the problem information    
    !---------------------------------------------------------------------------   
    subroutine initializeInteractionMatrices( problem )
    type(MicroMagProblem), intent(inout) :: problem         !> Struct containing the grid information
    
    
    !Setup the grid
    call setupGrid( problem%grid )
    
    !Demagnetization tensor matrix
    !call ComputeDemagfieldTensor( problem )
    call BuildNeighbourDemagTensor( problem, 2) ! hardcode 2 levels of neighbour cells for now

    !Anisotropy matrix
    call ComputeAnisotropyTerm3D( problem )
    
    !Exhange term matrix
    call ComputeExchangeTerm3D( problem%grid, problem%A_exch )
    
    
    end subroutine initializeInteractionMatrices
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Sets up the grid for the model. Only uniform grid is supported currently, but this will change in due time
    !> @param[inout] grid the struct that contains the information about the current grid
    !>-----------------------------------------
    
    subroutine setupGrid( grid )
    type(MicroMagGrid),intent(inout) :: grid            !> Grid information to be generated
    integer :: i,j,k,ind
    
    !Setup the grid depending on which type of grid it is
    if ( grid%gridType .eq. gridTypeUniform ) then
        
        !Allocate the grid
        allocate( grid%x(grid%nx,grid%ny,grid%nz),grid%y(grid%nx,grid%ny,grid%nz),grid%z(grid%nx,grid%ny,grid%nz) )
        allocate( grid%dV( grid%nx * grid%ny * grid%nz ) )
        allocate( grid%pts(grid%nx*grid%ny*grid%nz,3) )
        grid%pts(:,:) = 0.
        
        if ( grid%nx .gt. 1 ) then
            grid%dx = grid%Lx / grid%nx
            do i=1,grid%nx
                grid%x(i,:,:) = -grid%Lx/2 + (i-1) * grid%dx + grid%dx/2
            enddo
            
        else
            grid%x(:,:,:) = 0.
            grid%dx = grid%Lx
        endif
    
        if ( grid%ny .gt. 1 ) then
            grid%dy = grid%Ly / grid%ny
            do i=1,grid%ny
                grid%y(:,i,:) = -grid%Ly/2 + (i-1) * grid%dy + grid%dy/2
            enddo            
        else
            grid%y(:,:,:) = 0.
            grid%dy = grid%Ly
        endif

        if ( grid%nz .gt. 1 ) then
            grid%dz = grid%Lz / grid%nz
            do i=1,grid%nz
                grid%z(:,:,i) = -grid%Lz/2 + (i-1) * grid%dz + grid%dz/2
            enddo            
        else
            grid%z(:,:,:) = 0.
            grid%dz = grid%Lz
        endif

        grid%dV(:) = grid%dx * grid%dy * grid%dz
        do k=1,grid%nz
            do j=1,grid%ny
                do i=1,grid%nx
                    ind  = i + (j-1) * grid%nx + (k-1) * grid%ny * grid%nx
                    grid%pts( ind, 1 ) = grid%x(i,j,k)
                    grid%pts( ind, 2 ) = grid%y(i,j,k)
                    grid%pts( ind, 3 ) = grid%z(i,j,k)
                enddo
            enddo
        enddo
    endif
    
    end subroutine setupGrid
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2022
    !> @brief
    !> Calculates and returns the demag field tensor
    !> The calculations
    !> @param[inout] problem, the struct containing the problem
    
    !>-----------------------------------------
    subroutine ComputeDemagfieldTensor( problem )
    type(MicroMagProblem),intent(inout) :: problem                !> Grid data structure    
    
    type(MagTile),dimension(1) :: tile                            !> Tile representing the current tile under consideration
    real(DP),dimension(:,:),allocatable :: H                      !> The field and the corresponding evaluation point arrays
    integer :: i,j,k,nx,ny,nz,ntot,ind                            !> Internal counters and index variables
    integer :: i_a,j_a,k_a,nx_ave,ny_ave,nz_ave                   !> Internal counters and index variables for avering the demag tensor over the recieving tile
    real(DP),dimension(:),allocatable :: dx,dy,dz
    real(DP), dimension(:,:),allocatable :: pts_arr
    real(DP),dimension(:,:,:,:),allocatable :: Nout,Noutave             !> Temporary storage for the demag tensor            
    integer,dimension(4) :: indx_ele
    real :: rate
    integer :: c1,c2,cr,cm
    character(10) :: prog_str
    
    ! First initialize the system_clock
    call system_clock(count_rate=cr)
    call system_clock(count_max=cm)
    rate = REAL(cr)
    
    nx = problem%grid%nx
    ny = problem%grid%ny
    nz = problem%grid%nz
    ntot = nx * ny * nz
    
    !The number of elements to average the receiving tile over
    nx_ave = problem%N_ave(1)
    ny_ave = problem%N_ave(2)
    nz_ave = problem%N_ave(3)
    
    !Demag tensor components
    allocate( problem%Kxx(ntot,ntot), problem%Kxy(ntot,ntot), problem%Kxz(ntot,ntot) )
    allocate( problem%Kyy(ntot,ntot), problem%Kyz(ntot,ntot) )
    allocate( problem%Kzz(ntot,ntot) )
        
        
    if ( problem%demagTensorLoadState .gt. DemagTensorReturnMemory ) then
        call loadDemagTensorFromDisk( problem )
    else

        CALL SYSTEM_CLOCK(c1)
 
        !call mkl_set_num_threads(problem%nThreadsMatlab)
        !call omp_set_num_threads(problem%nThreadsMatlab)
        call omp_set_num_threads(1)
               
        if ( problem%grid%gridType .eq. gridTypeUniform ) then
            
            if (nx_ave*ny_ave*nz_ave > 1) then
                call displayGUIMessage( 'Averaging the N_tensor not supported for this tile type' )
            endif
        
            !$OMP PARALLEL DO SHARED(problem) PRIVATE(ind, tile, H, Nout, k, j, i)
        
            !for each element find the tensor for all evaluation points (i.e. all elements)
            do k=1,nz
                do j=1,ny                
                    do i=1,nx
                        !Setup template tile
                        tile(1)%tileType = 2 !(for prism)
                        !dimensions of the tile
                        tile(1)%a = problem%grid%dx
                        tile(1)%b = problem%grid%dy
                        tile(1)%c = problem%grid%dz
                        tile(1)%exploitSymmetry = 0 !0 for no and this is important
                        tile(1)%rotAngles(:) = 0. !ensure that these are indeed zero
                        tile(1)%M(:) = 0.
                        
                        !Set the center of the tile to be the current point
                        tile(1)%offset(1) = problem%grid%x(i,j,k)
                        tile(1)%offset(2) = problem%grid%y(i,j,k)
                        tile(1)%offset(3) = problem%grid%z(i,j,k)
                        
                        allocate(Nout(1,ntot,3,3))
                        allocate(H(ntot,3))
                        
                        call getFieldFromTiles( tile, H, problem%grid%pts, 1, ntot, Nout, .false. )
                        
                        !Copy Nout into the proper structure used by the micro mag model
                        ind = (k-1) * nx * ny + (j-1) * nx + i
                    
                        problem%Kxx(:,ind) = sngl(Nout(1,:,1,1))
                        problem%Kxy(:,ind) = sngl(Nout(1,:,1,2))
                        problem%Kxz(:,ind) = sngl(Nout(1,:,1,3))
                    
                        !Not stored due to symmetry  (Kxy = Kyx)
                        !Kyx(ind,:) = Nout(1,:,2,1)
                        problem%Kyy(:,ind) = sngl(Nout(1,:,2,2))
                        problem%Kyz(:,ind) = sngl(Nout(1,:,2,3))
                    
                        !Not stored due to symmetry (Kxz = Kzx)
                        !Kzx(ind,:) = Nout(1,:,3,1)
                        !Not stored due to symmetry (Kyz = Kzy)
                        !Kzy(ind,:) = Nout(1,:,3,2)
                        problem%Kzz(:,ind) = sngl(Nout(1,:,3,3))
                    
                        !Clean up
                        deallocate(Nout)
                        deallocate(H)
                    enddo
                enddo
            enddo
            
            !$OMP END PARALLEL DO
            
        elseif ( problem%grid%gridType .eq. gridTypeTetrahedron ) then
        
            if (nx_ave*ny_ave*nz_ave > 1) then
                call displayGUIMessage( 'Averaging the N_tensor not supported for this tile type' )
            endif
    
            !$OMP PARALLEL DO SHARED(problem) PRIVATE(ind, indx_ele, tile, H, Nout, pts_arr, i)
                        
            !for each element find the tensor for all evaluation points (i.e. all elements)
            do i=1,nx
                !Setup template tile
                tile(1)%tileType = 5 !(for tetrahedron)
                tile(1)%exploitSymmetry = 0 !0 for no and this is important
                tile(1)%rotAngles(:) = 0. !ensure that these are indeed zero
                tile(1)%M(:) = 0.
            
                indx_ele = problem%grid%elements(:,i)
                tile(1)%vert(:,:) = problem%grid%nodes(:,indx_ele)   
                
                allocate(Nout(1,ntot,3,3))
                allocate(H(ntot,3))
                
                allocate(pts_arr(ntot,3))
                pts_arr(:,1) =  problem%grid%pts(:,1)
                pts_arr(:,2) =  problem%grid%pts(:,2)
                pts_arr(:,3) =  problem%grid%pts(:,3)
                
                !call getFieldFromTiles( tile, H, problem%grid%pts, 1, ntot, Nout, .false. )
                call getFieldFromTiles( tile, H, pts_arr, 1, ntot, Nout, .false. )
                    
                !Copy Nout into the proper structure used by the micro mag model
                ind = i
                    
                problem%Kxx(:,ind) = sngl(Nout(1,:,1,1))
                problem%Kxy(:,ind) = sngl(Nout(1,:,1,2))
                problem%Kxz(:,ind) = sngl(Nout(1,:,1,3))
                    
                !Not stored due to symmetry  (Kxy = Kyx)
                !Kyx(ind,:) = Nout(1,:,2,1)
                problem%Kyy(:,ind) = sngl(Nout(1,:,2,2))
                problem%Kyz(:,ind) = sngl(Nout(1,:,2,3))
                    
                !Not stored due to symmetry (Kxz = Kzx)
                !Kzx(ind,:) = Nout(1,:,3,1)
                !Not stored due to symmetry (Kyz = Kzy)
                !Kzy(ind,:) = Nout(1,:,3,2)
                problem%Kzz(:,ind) = sngl(Nout(1,:,3,3))
                
                !Clean up
                deallocate(pts_arr)
                deallocate(Nout)
                deallocate(H)
            enddo
            
            !$OMP END PARALLEL DO
            
        elseif ( problem%grid%gridType .eq. gridTypeUnstructuredPrisms ) then
            !call displayGUIMessage( 'Constructing the Tensormap' )
            !call ConstructDemagTensorMap( problem )
            !call displayGUIMessage( 'Done constructing the Tensormap' )
            
            !$OMP PARALLEL DO SHARED(problem) PRIVATE(ind, tile, H, Nout,Noutave,dx,dy,dz,pts_arr)
            
            !for each element find the tensor for all evaluation points (i.e. all elements)
            do i=1,ntot
                !Setup template tile
                tile(1)%tileType = 2 !(for prism)
                tile(1)%exploitSymmetry = 0 !0 for no and this is important
                tile(1)%rotAngles(:) = 0. !ensure that these are indeed zero
                tile(1)%M(:) = 0.
            
                !dimensions of the tile
                tile(1)%a = problem%grid%abc(i,1)
                tile(1)%b = problem%grid%abc(i,2)
                tile(1)%c = problem%grid%abc(i,3)
                
                !Set the center of the tile to be the current point
                tile(1)%offset(1) = problem%grid%pts(i,1)
                tile(1)%offset(2) = problem%grid%pts(i,2)
                tile(1)%offset(3) = problem%grid%pts(i,3)
                
                allocate(Nout(1,ntot,3,3))
                allocate(H(ntot,3))
                
                allocate(Noutave(1,ntot,3,3))
                allocate(dx(ntot))
                allocate(dy(ntot))
                allocate(dz(ntot))
                allocate(pts_arr(ntot,3))
                
                Noutave(1,:,:,:) = 0;

                !Calculate the spacing between the points to do average in
                !Note that abc is the full side length of the tile - it is divided with 1/2 in the demag tensor calculation
                !to make it compatible to the expression in Smith_2010
                dx = problem%grid%abc(:,1)/(nx_ave+1)
                dy = problem%grid%abc(:,2)/(ny_ave+1)
                dz = problem%grid%abc(:,3)/(nz_ave+1)
                do k_a=1,nz_ave
                    do j_a=1,ny_ave                
                        do i_a=1,nx_ave
                            !x = -2; a = 6; N = 4; dx = a/(N+1); figure; hold all; plot(x,0,'kd'); plot(x-a/2,0,'bd'); plot(x+a/2,0,'bd'); plot((x-a/2)+(1:N)*dx,0,'k*');
                            pts_arr(:,1) =  (problem%grid%pts(:,1)-problem%grid%abc(:,1)/2)+dx(:)*i_a
                            pts_arr(:,2) =  (problem%grid%pts(:,2)-problem%grid%abc(:,2)/2)+dy(:)*j_a
                            pts_arr(:,3) =  (problem%grid%pts(:,3)-problem%grid%abc(:,3)/2)+dz(:)*k_a
                            !call getFieldFromTiles( tile, H, problem%grid%pts, 1, ntot, Nout, .false. )
                            call getFieldFromTiles( tile, H, pts_arr, 1, ntot, Nout, .false. )
                            
                            Noutave = Noutave+Nout
                            
                        enddo
                    enddo
                enddo
                
                Nout = Noutave/(nx_ave*ny_ave*nz_ave)
                
                !Copy Nout into the proper structure used by the micro mag model
                ind = i
                    
                problem%Kxx(:,ind) = sngl(Nout(1,:,1,1))
                problem%Kxy(:,ind) = sngl(Nout(1,:,1,2))
                problem%Kxz(:,ind) = sngl(Nout(1,:,1,3))
                    
                !Not stored due to symmetry  (Kxy = Kyx)
                !Kyx(ind,:) = Nout(1,:,2,1)
                problem%Kyy(:,ind) = sngl(Nout(1,:,2,2))
                problem%Kyz(:,ind) = sngl(Nout(1,:,2,3))
                    
                !Not stored due to symmetry (Kxz = Kzx)
                !Kzx(ind,:) = Nout(1,:,3,1)
                !Not stored due to symmetry (Kyz = Kzy)
                !Kzy(ind,:) = Nout(1,:,3,2)
                problem%Kzz(:,ind) = sngl(Nout(1,:,3,3))
                
                !Clean up
                deallocate(Nout)
                deallocate(H)
                deallocate(Noutave)
                deallocate(dx)
                deallocate(dy)
                deallocate(dz)
                deallocate(pts_arr)
            enddo
        
            !$OMP END PARALLEL DO
            
        endif
        
        CALL SYSTEM_CLOCK(c2)

        !Display the time to compute the demag tensor and its first entry
        !call displayGUIMessage( 'Time demag tensor:' )
        !write (prog_str,'(f10.3)') (c2 - c1)/rate
        !call displayGUIMessage( prog_str )
        
        !call displayGUIMessage( 'Kxx(1,1):' )
        !write (prog_str,'(f10.3)') problem%Kxx(1,1)
        !call displayGUIMessage( prog_str )      
        
        !Write the demag tensors to disk if asked to do so            
        if ( problem%demagTensorReturnState .gt. DemagTensorReturnMemory ) then
            call writeDemagTensorToDisk( problem )
        endif            
        
    endif
    
        
    !Make a sparse matrix out of the dense matrices by specifying a threshold
    if ( ( problem%demag_approximation .eq. DemagApproximationThreshold ) .or. ( problem%demag_approximation .eq. DemagApproximationThresholdFraction ) ) then     
        call ApplyThresholdDense( problem )
    elseif ( ( problem%demag_approximation .eq. DemagApproximationFFTThreshold ) .or. ( problem%demag_approximation .eq. DemagApproximationFFTThresholdFraction ) ) then
        call ApplyThresholdFFT( problem )
    endif
    
    
    end subroutine ComputeDemagfieldTensor
    
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Calculates the exhange term matrix
    !> which means it produces the differential operator d^2/dx^2 + d^2/dy^2 + d^2/dz^2 and returns this in the sparse matrix A
    !---------------------------------------------------------------------------   
    subroutine ComputeExchangeTerm3D( grid, A )
    type(MicroMagGrid),intent(in) :: grid             !> Struct containing the grid information    
    type(sparse_matrix_t),intent(inout) :: A          !> The returned matrix from the sparse matrix creator
    
    if ( grid%gridType .eq. gridTypeUniform ) then
        call ComputeExchangeTerm3D_Uniform( grid, A )
    elseif (( grid%gridType .eq. gridTypeTetrahedron ) .or. (grid%gridType .eq. gridTypeUnstructuredPrisms)) then
        !call ConvertExchangeTerm3D_NonUniform( grid, A )
    endif

    
    end subroutine ComputeExchangeTerm3D
    
        
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Calculates the exhange term matrix on a uniform grid
    !> which means it produces the differential operator d^2/dx^2 + d^2/dy^2 + d^2/dz^2 and returns this in the sparse matrix A
    !---------------------------------------------------------------------------   
    subroutine ComputeExchangeTerm3D_Uniform( grid, A )
    type(MicroMagGrid),intent(in) :: grid             !> Struct containing the grid information    
    type(sparse_matrix_t),intent(inout) :: A           !> The returned matrix from the sparse matrix creator
        
    integer :: stat                                   !> Status value for the various sparse matrix operations        
    type(MagTenseSparse_d) :: d2dx2, d2dy2, d2dz2       !> Sparse matrices for the double derivatives with respect to x, y and z, respectively.
    type(sparse_matrix_t) :: tmp                      !> Temporary sparse matrices used for internal calculations
    integer :: ind, ntot,colInd,rowInd                !> Internal counter for indexing, the total no. of elements in the current sparse matrix being manipulated    
    integer :: i,j,k,nx,ny,nz                         !> For-loop counters
    type(matrix_descr) :: descr                         !> Describes a sparse matrix operation
    real(DP) :: const
    
    !Find the three sparse matrices for the the individual directions, respectively. Then add them to get the total matrix
    !It is assumed that the magnetization vector to operate on is in fact a single column of Mx, My and Mz respectively.
    
    nx = grid%nx
    ny = grid%ny
    nz = grid%nz
    
    !----------------------------------d^2dx^2 begins -----------------------------!
    !Make the d^2/dx^2 matrix. The no. of non-zero elements is 3 * nx*ny*nz - 2 * ny * nz
    ntot = nz * ( ny * 4 + ny * (nx-2)*3 )
    allocate(d2dx2%values(ntot),d2dx2%cols(ntot),d2dx2%rows_start(nx*ny*nz),d2dx2%rows_end(nx*ny*nz))
    
    ind = 1
    rowInd = 1
    colInd = 1
    
    do k=1,nz
        do j=1,ny
            
            !The left boundary
            d2dx2%values(ind) = -1.
            d2dx2%cols(ind) = colInd
            d2dx2%rows_start(rowInd) = ind            
            ind = ind + 1
            
            d2dx2%values(ind) = 1.
            d2dx2%cols(ind) = colInd + 1
            d2dx2%rows_end(rowInd) = ind+1
            rowInd = rowInd + 1
            ind = ind + 1
            
            
            !Go through one row at a time
            do i=2,nx-1
                            
                !Left-most point
                d2dx2%values( ind ) = 1.
                d2dx2%cols(ind) = colInd
                !update where the row starts
                d2dx2%rows_start(rowInd) = ind                           
                ind = ind + 1
                
                !Center point
                d2dx2%values( ind ) = -2.
                d2dx2%cols(ind) = colInd + 1
                ind = ind + 1
                
                !Right-most point
                d2dx2%values( ind ) = 1.
                d2dx2%cols(ind) = colInd + 2
                d2dx2%rows_end(rowInd) = ind+1
                rowInd = rowInd + 1
                ind = ind + 1
                
                colInd = colInd + 1
            enddo            
            !The right boundary
            d2dx2%values(ind) = 1.
            d2dx2%cols(ind) = colInd
            !update where the row starts
            d2dx2%rows_start(rowInd) = ind            
            ind = ind + 1
            
            d2dx2%values(ind) = -1.
            d2dx2%cols(ind) = colInd+1
            d2dx2%rows_end(rowInd) = ind+1
            rowInd = rowInd + 1
            ind = ind + 1     
            
            colInd = colInd + 2
        enddo
    enddo
    
    !Multiply by the discretization
    d2dx2%values = d2dx2%values * 1./grid%dx**2
        
    !Create the sparse matrix for the d^2dx^2
    stat = mkl_sparse_d_create_csr ( d2dx2%A, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, d2dx2%rows_start, d2dx2%rows_end, d2dx2%cols, d2dx2%values)
    
    
    !----------------------------------d^2dx^2 ends -----------------------------!
    
    
    !----------------------------------d^2dy^2 begins ----------------------------!
    ntot = nz * ( nx * 2 + (ny-2) * nx * 3 + nx * 2 )
    !Make the d^2/dy^2 matrix. The no. of non-zero elements is 3 * nx*ny*nz - 2 * ny * nz just as for d^2dx^2
    allocate(d2dy2%values(ntot),d2dy2%cols(ntot),d2dy2%rows_start(nx*ny*nz),d2dy2%rows_end(nx*ny*nz))
    
    ind = 1
    rowInd = 1
    colInd = 1
    do k=1,nz
        !The bottom boundary
        do i=1,nx
            d2dy2%values(ind) = -1.
            d2dy2%cols(ind) = colInd
            d2dy2%rows_start(rowInd) = ind
            
            !increment to next element
            ind = ind + 1
            
            d2dy2%values(ind) = 1.
            d2dy2%cols(ind) = colInd + nx
            d2dy2%rows_end(rowInd) = ind+1
            rowInd = rowInd + 1
            
            !increment to next element
            ind = ind + 1
            
            colInd = colInd + 1
        enddo
        
        !Everything in between
        do j=2,ny-1
                    
            
            do i=1,nx
                
                !lower value
                d2dy2%values(ind) = 1.
                d2dy2%cols(ind) = colInd-nx
                d2dy2%rows_start(rowInd) = ind
                !increment to next element
                ind = ind + 1  
                
                !central value
                d2dy2%values(ind) = -2.
                d2dy2%cols(ind) = colInd
            
                !increment to next element
                ind = ind + 1
            
                !upper value
                d2dy2%values(ind) = 1.
                d2dy2%cols(ind) = colInd + nx
                d2dy2%rows_end(rowInd) = ind+1
                rowInd = rowInd + 1
                
                !increment to next element
                ind = ind + 1
                colInd = colInd + 1
            enddo
                        
        
        enddo
        
        !The top boundary    
        do i=1,nx
            !lower element
            d2dy2%values(ind) = 1.
            d2dy2%cols(ind) = colInd - nx
            d2dy2%rows_start(rowInd) = ind
            
            !increment to next element
            ind = ind + 1            
            
            !central element
            d2dy2%values(ind) = -1.
            d2dy2%cols(ind) = colInd
            d2dy2%rows_end(rowInd) = ind + 1
            rowInd = rowInd + 1
            
            ind = ind + 1
            colInd = colInd + 1
            
        enddo
    enddo
    
    !Multiply by the discretization
    d2dy2%values = d2dy2%values * 1./grid%dy**2
        
    !Create the sparse matrix for the d^2dy^2
    stat = mkl_sparse_d_create_csr ( d2dy2%A, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, d2dy2%rows_start, d2dy2%rows_end, d2dy2%cols, d2dy2%values)
    
    !----------------------------------d^2dy^2 ends ----------------------------!
    
    
    !----------------------------------d^2dz^2 begins ----------------------------!
    if ( nz .gt. 1 ) then
        !Make the d^2/dz^2 matrix. The no. of non-zero elements is 3 * nx*ny*nz - 2 * ny * nz just as for d^2dx^2 and d^2dy^2
        ntot = 2 * 2 * nx * ny + 3 * (nz-2) * nx * ny
        !Make the d^2/dy^2 matrix. The no. of non-zero elements is 3 * nx*ny*nz - 2 * ny * nz just as for d^2dx^2
        allocate(d2dz2%values(ntot),d2dz2%cols(ntot),d2dz2%rows_start(nx*ny*nz),d2dz2%rows_end(nx*ny*nz))
    
        ind = 1
        rowInd = 1
        colInd = 1
        !The z=1 face        
        do j=1,ny
            do i=1,nx
                !central value
                d2dz2%values(ind) = -1.
                d2dz2%cols(ind) = colInd
                d2dz2%rows_start(rowInd) = ind
            
                !increment position
                ind = ind + 1
            
                !right-most value
                d2dz2%values(ind) = 1.
                d2dz2%cols(ind) = nx * ny + colInd
                d2dz2%rows_end(rowInd) = ind + 1
                rowInd = rowInd + 1
            
                !increment position
                ind = ind + 1
                colInd = colInd + 1
            enddo
        enddo
        !Everything in between
        do k=2,nz-1            
            do j=1,ny
                do i=1,nx
                    !left-most value
                    d2dz2%values(ind) = 1.
                    d2dz2%cols(ind) = colInd - nx * ny
                    d2dz2%rows_start(rowInd) = ind
            
                    !increment position
                    ind = ind + 1
                
                    !central value
                    d2dz2%values(ind) = -2.
                    d2dz2%cols(ind) = colInd
                            
                    !increment position
                    ind = ind + 1
                
                    !right-most value
                    d2dz2%values(ind) = 1.
                    d2dz2%cols(ind) = colInd + nx * ny
                    d2dz2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                
                    !increment position
                    ind = ind + 1
                    colInd = colInd + 1
                enddo
            enddo
        
        enddo
    
    
        !The z=nz face        
        do j=1,ny
            do i=1,nx
                
                !left-most value
                d2dz2%values(ind) = 1.
                d2dz2%cols(ind) = colInd - nx * ny
                d2dz2%rows_start(rowInd) = ind
            
                !increment position
                ind = ind + 1
            
                !central value
                d2dz2%values(ind) = -1.
                d2dz2%cols(ind) = colInd
                d2dz2%rows_end(rowInd) = ind + 1
                rowInd = rowInd + 1
            
                !increment position
                ind = ind + 1
                colInd = colInd + 1
            
            enddo
        enddo
    
    
        !Multiply by the discretization
        d2dz2%values = d2dz2%values * 1./grid%dz**2
        
        !Create the sparse matrix for the d^2dz^2
        stat = mkl_sparse_d_create_csr ( d2dz2%A, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, d2dz2%rows_start, d2dz2%rows_end, d2dz2%cols, d2dz2%values)
    endif
    
    !----------------------------------d^2dz^2 ends ----------------------------!
    
            
    !Finally, add up the three diagonals and store in the output sparse matrix, A
    !store the results temporarily in tmp4
    
    !call writeSparseMatrixToDisk( d2dz2%A, nx*ny*nz, 'd2dz2.dat' )
    !call writeSparseMatrixToDisk( d2dy2%A, nx*ny*nz, 'd2dy2.dat' )
    
    const = 1.
    stat = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, d2dx2%A, const, d2dy2%A, tmp)    
    
    !call writeSparseMatrixToDisk( tmp, nx*ny*nz, 'A_exch.dat' )
    
    if ( nz .gt. 1 ) then    
        stat = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, d2dz2%A, const, tmp, A)
        !clean up        
        deallocate(d2dz2%values,d2dz2%cols,d2dz2%rows_start,d2dz2%rows_end)
        stat = mkl_sparse_destroy (d2dz2%A)
    else
        descr%type = SPARSE_MATRIX_TYPE_GENERAL
        descr%mode = SPARSE_FILL_MODE_FULL
        descr%diag = SPARSE_DIAG_NON_UNIT
        stat = mkl_sparse_copy ( tmp, descr, A )
        stat = mkl_sparse_destroy (tmp)
    endif
    
    !call writeSparseMatrixToDisk( A, nx*ny*nz, 'A_total.dat' )
    
    !clean up
    deallocate(d2dx2%values,d2dx2%cols,d2dx2%rows_start,d2dx2%rows_end)
    deallocate(d2dy2%values,d2dy2%cols,d2dy2%rows_start,d2dy2%rows_end)
    
    stat = mkl_sparse_destroy (d2dx2%A)
    stat = mkl_sparse_destroy (d2dy2%A)
    
    !call writeSparseMatrixToDisk( A, nx*ny*nz, 'A_exch.dat' )
    
    end subroutine ComputeExchangeTerm3D_Uniform
       
    !>-----------------------------------------
    !> @author Rasmus Bjørk, rabj@dtu.dk, DTU, 2020
    !> @brief
    !> Converts the loaded information from Matlab in CSR 
    !> format to a CSR MKL type
    !---------------------------------------------------------------------------   
    subroutine ConvertExchangeTerm3D_NonUniform(grid, A)
    type(MicroMagGrid),intent(in) :: grid             !> Struct containing the grid information    
    type(sparse_matrix_t),intent(out) :: A            !> The returned matrix from the sparse matrix creator
                
    integer :: stat                                   !> Status value for the various sparse matrix operations        
    integer :: nx,ny,nz                               !> Dimensions
    type(matrix_descr) :: descr                       !> Describes a sparse matrix operation
    
    nx = grid%nx
    ny = grid%ny
    nz = grid%nz
    
    stat = mkl_sparse_d_create_csr ( A, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, grid%A_exch_load%rows_start, grid%A_exch_load%rows_end, grid%A_exch_load%cols, grid%A_exch_load%values)
        
    end subroutine ConvertExchangeTerm3D_NonUniform
    
    
    !>-----------------------------------------
    !> @author Rasmus Bjørk, rabj@dtu.dk, DTU, 2020
    !> @brief
    !> Calculates the anisotropy term sparse matrix assuming the effective field anisotropy is linear in m    
    !> @param[inout] problem the data structure containing the problem
    !---------------------------------------------------------------------------   
    subroutine ComputeAnisotropyTerm3D( problem )
    type(MicroMagProblem),intent(inout) :: problem       !> Struct containing the problem
    
    call ComputeAnisotropyTerm3D_General( problem )
    
    end subroutine ComputeAnisotropyTerm3D
    
    !>-----------------------------------------
    !> @author Rasmus Bjørk, rabj@dtu.dk, DTU, 2020
    !> @brief
    !> Calculates the anisotropy term matrix on any grid  
    !> @param[inout] problem the data structure containing the problem
    !---------------------------------------------------------------------------   
    subroutine ComputeAnisotropyTerm3D_General( problem )
    type(MicroMagProblem),intent(inout) :: problem             !> Struct containing the problem
    
    integer :: nx,ny,nz,ntot, i
    
    !--- We use the general formulation introduced in updateAnisotropy
    !--- If the user has specified a value for the uniaxial anisotropy or the cubic anisotropy, we use transform
    !--- those to the general matrix formulation
    allocate(problem%Kfact_arr(size(problem%Ms),6,3))
    problem%Kfact_arr(:,:,:) = 0.0
    
    if (any(problem%K0 .ne. 0)) then !Uniaxial anisotropy
        call displayGUIMessage( 'Assuming uniaxial anisotropy' )
        do i = 1,size(problem%Ms)
            problem%Kfact_arr(i,1,3) = problem%K0(i) / ( mu0 * problem%Ms(i) )
        enddo
        
        !Set the crystal axis as the user only specified u_ea
        problem%CrystalAxis(:,1,3) =  problem%u_ea(:,1)
        problem%CrystalAxis(:,2,3) =  problem%u_ea(:,2)
        problem%CrystalAxis(:,3,3) =  problem%u_ea(:,3)
        
    elseif (any(problem%K1 .ne. 0)) then !Cubic anisotropy
        call displayGUIMessage( 'Assuming cubic anisotropy' )
        do i = 1,size(problem%Ms)
            problem%Kfact_arr(i,3,1) = problem%K1(i) / ( mu0 * problem%Ms(i) )
            problem%Kfact_arr(i,3,2) = problem%K1(i) / ( mu0 * problem%Ms(i) )
            problem%Kfact_arr(i,3,3) = problem%K1(i) / ( mu0 * problem%Ms(i) )
            problem%Kfact_arr(i,6,1) = problem%K2(i) / ( mu0 * problem%Ms(i) )
        enddo
    else !General anisotropy
        do i = 1,size(problem%Ms)
            problem%Kfact_arr(i,:,:) = problem%K0_arr(i,:,:) / ( mu0 * problem%Ms(i) )
        enddo
    endif
    
    end subroutine ComputeAnisotropyTerm3D_General

    


    !================== Helper functions for FMM correction ==========================================

        !---------------- for getting self field correction ----------------------------------------

    !-------------------------------
    !< Function add_self_correction_fmm
    !< Add the self-demag correction to the FMM-calculated demag field
    !< assumes uniform cartesian grid
    subroutine add_self_correction(problem, solution)
        implicit none
        type(MicroMagProblem),  intent(in)    :: problem
        type(MicroMagSolution), intent(inout) :: solution

        ! --- locals ---
        integer :: ntot, i
        type(MagTile), dimension(1) :: tile
        real(DP), dimension(1,3) :: pts
        real(DP), dimension(1,3) :: H_dummy
        !real(DP), dimension(1,1,3,3) :: Nout
        real(DP),dimension(:,:,:,:),allocatable :: Nout
        real(DP) :: Nloc(3,3)
        real(DP) :: dHx, dHy, dHz, dH(3)
        real(DP) :: dx, dy, dz

        ! --- grid sizes & count ---
        ntot = problem%grid%nx * problem%grid%ny * problem%grid%nz
        dx = real(problem%grid%dx, DP)
        dy = real(problem%grid%dy, DP)
        dz = real(problem%grid%dz, DP)

        ! --- build a single prism tile identical to a voxel ---
        tile(1)%tileType         = 2                 ! prism
        tile(1)%a                = dx                ! full side lengths (as in your tensor builder)
        tile(1)%b                = dy
        tile(1)%c                = dz
        tile(1)%exploitSymmetry  = 0                 ! critical: no symmetry tricks
        tile(1)%rotAngles(:)     = 0.0_DP            ! axis-aligned
        tile(1)%M(:)             = 0.0_DP            ! not used; we just want the tensor
        tile(1)%offset(1)        = 0.0_DP            ! centre arbitrarily at origin
        tile(1)%offset(2)        = 0.0_DP
        tile(1)%offset(3)        = 0.0_DP

        ! evaluate at the tile's centre (same point)
        pts(1,1) = 0.0_DP
        pts(1,2) = 0.0_DP
        pts(1,3) = 0.0_DP

       ! print *, " Adding self-demag correction to FMM demag field..."

        ! --- get the 3x3 self demag tensor once ---
        call getFieldFromTiles( tile, H_dummy, pts, 1, 1, Nout, .false. )
        Nloc(:,:) = Nout(1,1,:,:)

        ! --- apply H += -Nloc * M for every cell (in place) ---
        !$omp parallel do private(i,dH) shared(problem,solution,Nloc,ntot) default(none)
        do i = 1, ntot
            !dHx = -( Nloc(1,1)*real(solution%Mx_s(i),DP) + Nloc(1,2)*real(solution%My_s(i),DP) + Nloc(1,3)*real(solution%Mz_s(i),DP) )
            !dHy = -( Nloc(2,1)*real(solution%Mx_s(i),DP) + Nloc(2,2)*real(solution%My_s(i),DP) + Nloc(2,3)*real(solution%Mz_s(i),DP) )
            !dHz = -( Nloc(3,1)*real(solution%Mx_s(i),DP) + Nloc(3,2)*real(solution%My_s(i),DP) + Nloc(3,3)*real(solution%Mz_s(i),DP) )

            dH = matmul(Nloc, [solution%Mx_s(i), solution%My_s(i),solution%Mz_s(i)] ) * problem%Ms(i)

            solution%HmX(i) = solution%HmX(i) - real(dH(1), SP)
            solution%HmY(i) = solution%HmY(i) - real(dH(2), SP)
            solution%HmZ(i) = solution%HmZ(i) - real(dH(3), SP)
        end do
        !$omp end parallel do
    end subroutine add_self_correction



        !------------------------------------------------------------------------------------------------

        !------------------- for getting near-field correction ------------------------------------------
    pure subroutine dipole_tensor_3x3(Rvec, Kdip)
        implicit none
        ! Kdip = (1/(4π r^3)) * (3 r̂ r̂^T - I), mapping M -> H (SI)
        real(8), intent(in)  :: Rvec(3)
        real(8), intent(out) :: Kdip(3,3)
        real(8), parameter   :: fourpi = 12.56637061435917295385d0
        real(8), parameter   :: epsR   = 1.0d-15
        real(8) :: r2, rmag, invR3, ux, uy, uz

        r2  = Rvec(1)*Rvec(1) + Rvec(2)*Rvec(2) + Rvec(3)*Rvec(3)
        rmag = sqrt(r2)

        if (rmag <= epsR) then
            Kdip = 0.0d0
            return
        end if

        invR3 = 1.0d0 / (fourpi * r2 * rmag)
        ux = Rvec(1)/rmag
        uy = Rvec(2)/rmag
        uz = Rvec(3)/rmag

        Kdip(1,1) = (3.0d0*ux*ux - 1.0d0) * invR3
        Kdip(1,2) = (3.0d0*ux*uy)         * invR3
        Kdip(1,3) = (3.0d0*ux*uz)         * invR3
        Kdip(2,1) = Kdip(1,2)
        Kdip(2,2) = (3.0d0*uy*uy - 1.0d0) * invR3
        Kdip(2,3) = (3.0d0*uy*uz)         * invR3
        Kdip(3,1) = Kdip(1,3)
        Kdip(3,2) = Kdip(2,3)
        Kdip(3,3) = (3.0d0*uz*uz - 1.0d0) * invR3
    end subroutine dipole_tensor_3x3

    pure logical function is_neighbour(idx_t, idx_s, offset, size, pitch, rad_cells) result(mask)
        implicit none
        integer(4), intent(in) :: idx_t, idx_s, rad_cells
        real(8),   intent(in) :: offset(:,:), size(:,:), pitch(3)
        real(8) :: dx, dy, dz, tx, ty, tz

        if (idx_t == idx_s) then
            mask = .false.
            return
        end if

        dx = abs(offset(idx_s,1) - offset(idx_t,1))
        dy = abs(offset(idx_s,2) - offset(idx_t,2))
        dz = abs(offset(idx_s,3) - offset(idx_t,3))

        tx = rad_cells * pitch(1) + 0.5d0*(size(idx_s,1) + size(idx_t,1))
        ty = rad_cells * pitch(2) + 0.5d0*(size(idx_s,2) + size(idx_t,2))
        tz = rad_cells * pitch(3) + 0.5d0*(size(idx_s,3) + size(idx_t,3))

        mask = (dx <= tx) .and. (dy <= ty) .and. (dz <= tz)
    end function is_neighbour






subroutine BuildNeighbourDemagTensor(problem, radius_cells)
  ! Construct exact prism demag tensors for near neighbours of every cell.
  ! Outputs:
  !   nbr_idx(ntot, nneigh_max)  = linear indices of neighbours (−1 if none)
  !   Nnbr   (ntot, nneigh_max, 3,3) = demag tensor at target-centre from each neighbour tile
  implicit none
  type(MicroMagProblem), intent(inout)  :: problem
  integer,               intent(in)  :: radius_cells
  integer, dimension(:,:), pointer :: nbr_idx(:,:)
  real(SP),dimension(:,:,:,:), pointer:: Nnbr(:,:,:,:)
  ! --- locals ---
  integer :: nx, ny, nz, ntot, r, nneigh_max
  integer :: i, j, k, ii, jj, kk, di, dj, dk, m, lin_t, lin_s
  type(MagTile), allocatable :: tiles(:)
  real(DP) :: dx, dy, dz
  real(DP), dimension(1,3) :: pts_t
  real(DP), allocatable :: H_dummy(:,:)
  real(DP), allocatable :: Nout(:,:,:,:)  ! (n_tiles, 1, 3, 3)
 integer :: mj
  ! --- only uniform prism grids in this first pass ---
  if (problem%grid%gridType /= gridTypeUniform) then
    stop 'BuildNeighbourDemagTensor: only gridTypeUniform supported in this routine.'
  end if




  !print *, " starting BuildNeighbourDemagTensor with radius_cells = ", radius_cells

  nx   = problem%grid%nx
  ny   = problem%grid%ny
  nz   = problem%grid%nz
  ntot = nx * ny * nz
  r    = radius_cells

  dx = real(problem%grid%dx, DP)
  dy = real(problem%grid%dy, DP)
  dz = real(problem%grid%dz, DP)

  ! Maximum neighbours in a 3D Chebyshev ball of radius r (remove self)
  nneigh_max = (2*r + 1)**3 - 1

  ! Allocate outputs
  if (associated(problem%nbr_idx)) deallocate(problem%nbr_idx)
  if (associated(problem%Nnbr))    deallocate(problem%Nnbr)
  allocate(nbr_idx(ntot, nneigh_max))
  allocate(Nnbr(   ntot, nneigh_max, 3, 3))

  problem%nbr_idx => nbr_idx
problem%Nnbr    => Nnbr



  nbr_idx = -1
  Nnbr    = 0.0_DP
  
  allocate(H_dummy(1,3))

  ! OpenMP-friendly outer loop over target cells
  !$omp parallel do collapse(3) default(none) &
  !$omp shared(nx,ny,nz,ntot,r,nneigh_max,problem,nbr_idx,Nnbr,dx,dy,dz) &
  !$omp private(i,j,k,lin_t,ii,jj,kk,di,dj,dk,m,lin_s,tiles,pts_t,H_dummy,Nout) schedule(static)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        ! Linear index of target tile (column-major flattening)
        lin_t = (k-1)*nx*ny + (j-1)*nx + i

        ! Enumerate neighbours in Chebyshev ball of radius r
        m = 0
        do dk = -r, r
          kk = k + dk
          if (kk < 1 .or. kk > nz) cycle
          do dj = -r, r
            jj = j + dj
            if (jj < 1 .or. jj > ny) cycle
            do di = -r, r
              ii = i + di
              if (ii < 1 .or. ii > nx) cycle
              if (di == 0 .and. dj == 0 .and. dk == 0) cycle  ! skip self
              m = m + 1
              nbr_idx(lin_t, m) = (kk-1)*nx*ny + (jj-1)*nx + ii
            end do
          end do
        end do

        if (m > 0) then
          ! Build MagTile array for these m neighbours
          allocate(tiles(m))

          do mj = 1, m
            tiles(mj)%tileType        = 2           ! prism
            tiles(mj)%a               = dx
            tiles(mj)%b               = dy
            tiles(mj)%c               = dz
            tiles(mj)%exploitSymmetry = 0
            tiles(mj)%rotAngles(:)    = 0.0
            tiles(mj)%M(:)            = 0.0
          end do

          ! Set offsets to neighbour centres using grid coordinate arrays
          ! problem%grid%x(i,j,k), y(i,j,k), z(i,j,k)
          do jj = 1, m
            lin_s = nbr_idx(lin_t, jj)
            ! Recover (is, js, ks) for lin_s if you don’t have direct xyz arrays:
            ! But since you have problem%grid%pts(lin, :), we can use that directly:
            tiles(jj)%offset(1) = problem%grid%pts(lin_s, 1)
            tiles(jj)%offset(2) = problem%grid%pts(lin_s, 2)
            tiles(jj)%offset(3) = problem%grid%pts(lin_s, 3)
          end do

          ! Target point is the centre of the target tile
          pts_t(1,1) = problem%grid%pts(lin_t, 1)
          pts_t(1,2) = problem%grid%pts(lin_t, 2)
          pts_t(1,3) = problem%grid%pts(lin_t, 3)

          ! Small temporaries for the call
          allocate(Nout(m,1,3,3))

          call getFieldFromTiles( tiles, H_dummy, pts_t, m, 1, Nout, .false. )

          ! Store tensors in (t, neighbour-slot, 3, 3)
          Nnbr(lin_t, 1:m, :, :) = Nout(:,1,:,:)

          deallocate(Nout)
          deallocate(tiles)
        end if
            
            ! Unused neighbour slots (m+1:nneigh_max) remain {idx=-1, tensor=0}
        end do
    end do
    end do
    !$omp end parallel do
    deallocate(H_dummy)

      !print *, " finished BuildNeighbourDemagTensor with radius_cells = ", radius_cells


    !problem%nbr_idx => nbr_idx
    !problem%Nnbr => Nnbr
end subroutine BuildNeighbourDemagTensor


subroutine add_neighbour_correction(problem, solution)
  ! Finite-size neighbour patch for FMM demag field:
  ! H += sum_m ( Nnbr(t,m) - Kdip(R_tj)*V ) * M_j
  implicit none
  type(MicroMagProblem),  intent(in)    :: problem
  type(MicroMagSolution), intent(inout) :: solution

  integer  :: ntot, nneigh_max, t, m, jidx
  real(DP) :: dx, dy, dz, volj
  real(DP) :: xt, yt, zt, xj, yj, zj
  real(DP) :: Rvec(3), Kdip(3,3), Kloc(3,3)
  real(DP) :: dH(3)
  real(DP) :: Mj(3)

  ! --- guards ---
  if (.not. associated(problem%nbr_idx)) stop 'add_neighbour_correction: nbr_idx not associated'
  if (.not. associated(problem%Nnbr))    stop 'add_neighbour_correction: Nnbr not associated'


  !print *, " Adding neighbour-demag correction to FMM demag field..."

  ntot       = problem%grid%nx * problem%grid%ny * problem%grid%nz
  nneigh_max = size(problem%nbr_idx, 2)

  dx = real(problem%grid%dx, DP)
  dy = real(problem%grid%dy, DP)
  dz = real(problem%grid%dz, DP)
  volj = dx * dy * dz

  !$omp parallel do default(none) &
  !$omp shared(solution, problem, nneigh_max, volj, ntot)  &
  !$omp private(t,m,jidx,xt,yt,zt,xj,yj,zj,Rvec,Kdip,Kloc,dH,Mj)
  do t = 1, ntot
    dH = 0.0_DP

    ! target centre
    xt = real(problem%grid%pts(t,1), DP)
    yt = real(problem%grid%pts(t,2), DP)
    zt = real(problem%grid%pts(t,3), DP)

    do m = 1, nneigh_max
      jidx = problem%nbr_idx(t, m)
      if (jidx < 0) cycle   ! empty slot (boundary)

      ! neighbour j magnetisation (single -> double)
      Mj(1) = real(solution%Mx_s(jidx), DP)
      Mj(2) = real(solution%My_s(jidx), DP)
      Mj(3) = real(solution%Mz_s(jidx), DP)

      ! exact prism tensor for (j -> t), SP -> DP
      Kloc(:,:) = real(problem%Nnbr(t, m, :, :), DP)

      ! dipole tensor at target from source centre (per-unit M), then scale by V
      xj = real(problem%grid%pts(jidx,1), DP)
      yj = real(problem%grid%pts(jidx,2), DP)
      zj = real(problem%grid%pts(jidx,3), DP)
      Rvec(1) = xt - xj
      Rvec(2) = yt - yj
      Rvec(3) = zt - zj
      call dipole_tensor_3x3(Rvec, Kdip)
      Kdip = Kdip * volj

      ! accumulate (Kloc - Kdip) * Mj
      dH = dH + matmul( Kloc - Kdip, Mj * problem%Ms(jidx) )
    end do

    ! add correction in place (double -> single)
    solution%HmX(t) = solution%HmX(t) - real(dH(1), kind=SP)
    solution%HmY(t) = solution%HmY(t) - real(dH(2), kind=SP)
    solution%HmZ(t) = solution%HmZ(t) - real(dH(3), kind=SP)
  end do
  !$omp end parallel do

end subroutine add_neighbour_correction



        !------------------------------------------------------------------------------------------------
    !============================================================================================================




    
end module LandauLifshitzSolution
    