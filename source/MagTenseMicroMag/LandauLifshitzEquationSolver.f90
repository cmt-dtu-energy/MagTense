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

#if USE_FMM3D
    use fmm3d_tree_mod
    use fmm_nbor_tensor_mod
#endif

    use trace_mod
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
    type(MicroMagProblem),intent(inout) :: prob          !> The problem data structure
    type(MicroMagSolution),intent(inout) :: sol          !> The solution data structure
    type(MicroMagGridInfo) :: gridinfo                   !> The grid information structure
    integer :: ntot,i,j,k,ind,nt,nt_Hext,stat            !> total no. of tiles (Note from F. Durhuus: n_Hext would make more sense than nt_Hext)
    integer :: n_reject_max
    procedure(dydt_fct), pointer :: fct                  !> Input function pointer for the function to be integrated
    procedure(no_argument_fct), pointer :: fct_thermal   !> Function pointer for the function that updates the stochastic thermal field
    procedure(callback_fct),pointer :: cb_fct            !> Callback function for displaying progress
    real(DP),dimension(:,:,:),allocatable :: M_out       !> Internal buffer for the solution (M) on the form (3*ntot,nt)
    real(DP) :: A0_max
    real(DP) :: t_step, gamma, kB, alpha0, alphaGilbert  !> Parameters needed for thermal factor
    real(DP),dimension(:),allocatable :: volCells        !> Cell volumes
    character*(100) :: prog_str 
    real :: rate
    integer :: c1,c2,cr,cm 
    integer, save :: itimer = 0
    call trace%begin( "SolveLandauLifshitzEquation", itimer=itimer, verbose=1 )
    
    ! First initialize the system_clock
    call system_clock(count_rate=cr)
    call system_clock(count_max=cm)
    rate = REAL(cr)
    
    write(prog_str,'(A37, A8, A12)') 'MagTense version 2.3.0, compiled on: ', __TIME__, __DATE__ 
    call displayGUIMessage( trim(prog_str) ) 
    
    !Save internal representation of the problem and the solution
    gb_solution = sol
    gb_problem = prob
    
    !Calculate the local scaled coefficients for the LLG equation
    !"J" : exchange term
    if (( gb_problem%grid%gridType .eq. gridTypeTetrahedron ) .or. (gb_problem%grid%gridType .eq. gridTypeUnstructuredPrisms)) then
        A0_max = maxval( gb_problem%A0 )
        gb_problem%Jfact = A0_max / ( mu0 * gb_problem%Ms )   ! Normalized by the largest exchange factor, as this is needed in the exchange calculation for the unstructured meshed
    else
        gb_problem%Jfact = gb_problem%A0 / ( mu0 * gb_problem%Ms )
    endif
    
    !"M" : demagnetization term
    gb_problem%Mfact = gb_problem%Ms
    !"K" : anisotropy term
    !TODO - is this correct or should it be allocated first????
    gb_problem%Kfact = gb_problem%K0 / ( mu0 * gb_problem%Ms )
    
    ntot = gb_problem%grid%nx * gb_problem%grid%ny * gb_problem%grid%nz
    allocate( gb_solution%pts(ntot,3) )
    
    !Analyze the mesh, if needed
    if ( gb_problem%grid%gridType .eq. gridTypeUnstructuredPrisms ) then
        if ( gb_problem%passExch .eq. passExchTrue) then
            call displayGUIMessage( 'Passing exchange matrix' )
            call passDifferentialOperators(gb_problem)
        else    
            call CartesianUnstructuredMeshAnalysis(gb_problem%grid%pts, gb_problem%grid%abc, gb_solution%gridinfo)
        endif
    endif

    if (( gb_problem%grid%gridType .eq. gridTypeTetrahedron ) .or. (gb_problem%grid%gridType .eq. gridTypeUnstructuredPrisms)) then
        do i=1,gb_problem%grid%nx
            gb_solution%pts( i, 1 ) = gb_problem%grid%pts( i, 1 )
            gb_solution%pts( i, 2 ) = gb_problem%grid%pts( i, 2 )
            gb_solution%pts( i, 3 ) = gb_problem%grid%pts( i, 3 )
        end do
    endif   
    
    if ( gb_problem%grid%gridType .eq. gridTypeUniform ) then
        !Setup the grid
        call setupGrid( gb_problem%grid )
        
        !Set up the mesh for a uniform grid
        allocate( gb_problem%A0_map(gb_problem%grid%nx,gb_problem%grid%ny,gb_problem%grid%nz) )
        if ( gb_problem%grid%gridType .eq. gridTypeUniform ) then   
            do k=1,gb_problem%grid%nz
                do j=1,gb_problem%grid%ny            
                    do i=1,gb_problem%grid%nx
                        ind = i + (j-1) * gb_problem%grid%nx + (k-1) * gb_problem%grid%nx * gb_problem%grid%ny
                        gb_solution%pts(ind,1) = gb_problem%grid%x(i,j,k)
                        gb_solution%pts(ind,2) = gb_problem%grid%y(i,j,k)
                        gb_solution%pts(ind,3) = gb_problem%grid%z(i,j,k)
                    
                        gb_problem%A0_map(i,j,k) = gb_problem%A0(ind)
                    end do
                end do
            end do
        endif
    endif
        
    
    call displayGUIMessage( 'Initializing matrices' )
    !Calculate the interaction matrices
    call initializeInteractionMatrices( gb_problem, gb_solution )
    
    !write(prog_str,'(A37, F8.4, A5)') 'Demagnetization tensor memory usage: ', 6*storage_size(gb_problem%Kxx)*ntot/(8*2**30), ' gigabytes'
    !call displayGUIMessage( trim(prog_str) )    
    
    if ( gb_problem%useDemag .eq. useDemagTrue ) then
        !Copy the demag tensor to CUDA
        if ( gb_problem%useCuda .eq. useCudaTrue ) then
            call displayGUIMessage( 'Copying to CUDA' ) 
#if USE_CUDA
#if USE_FMM3D 
            !------------- if use_fmm then copy sparse nbr_corr tensor else copy normal demag tensor --------------
            if ( gb_problem%use_fmm) then
              call displayGUIMessage( 'copying nbr_corr for FMM nbor correction')  
              call cudaInit_sparse( gb_problem%K_fmm_s )    
            !----------------------------------------------------------------------------------------------
            !------------ else copy the normal demag tensor ------------------------------------------------
            else 
              !Initialize the Cuda arrays and load the demag tensors into the GPU memory
              if ( ( gb_problem%demag_approximation .eq. DemagApproximationThreshold ) .or. ( gb_problem%demag_approximation .eq. DemagApproximationThresholdFraction ) ) then
                  !If the matrices are sparse
                call displayGUIMessage( 'Initializing sparse matrices on CUDA' )
                  call cudaInit_sparse( gb_problem%K_s )            
              else
                  !if the matrices are dense 
                call displayGUIMessage( 'Initializing dense matrices on CUDA' )
                  call cudaInit_s( gb_problem%Kxx, gb_problem%Kxy, gb_problem%Kxz, gb_problem%Kyy, gb_problem%Kyz, gb_problem%Kzz )
              endif
            end if
            !----------------------------------------------------------------------------------------------
#else
            !Initialize the Cuda arrays and load the demag tensors into the GPU memory
            if ( ( gb_problem%demag_approximation .eq. DemagApproximationThreshold ) .or. ( gb_problem%demag_approximation .eq. DemagApproximationThresholdFraction ) ) then
                !If the matrices are sparse
                call displayGUIMessage( 'Initializing sparse matrices on CUDA' )
                call cudaInit_sparse( gb_problem%K_s )        
             
            else
                !if the matrices are dense 
                call displayGUIMessage( 'Initializing dense matrices on CUDA' )

                call cudaInit_s( gb_problem%Kxx, gb_problem%Kxy, gb_problem%Kxz, gb_problem%Kyy, gb_problem%Kyz, gb_problem%Kzz )

            !========== optionally dump CUDA version of demag tensor =????? ==================
            !     used for comparing GPU vs. CPU tensor
            !call cudaDumpDemagDense("CUDA_dense_K_matrices.bin")
            !=================================================================================
        endif
#endif
#else
            call displayGUIMessage( 'MagTense not compiled with CUDA - exiting!' )
            stop
#endif
        endif
    endif
   
    
    call displayGUIMessage( 'Initializing solution' )
    !Initialize the solution, i.e. allocate various arrays
    call initializeSolution( gb_problem, gb_solution )
         
    !Set the initial values for m (remember that M is organized such that mx = m(1:ntot), my = m(ntot+1:2*ntot), mz = m(2*ntot+1:3*ntot)
    allocate(gb_solution%t_out(size(gb_problem%t)))
        
    
    call displayGUIMessage( 'Running solution' )
    
    !Do the solution
    fct => dmdt_fct
    fct_thermal => updateThermal_wrapper
    cb_fct => displayGUIProgressMessage
    
    gb_solution%HextInd = 1
    
    !Number of timesteps (for each field value, right?)
    nt = size( gb_problem%t )

    !Whether to include the thermal field
    !If the temperature is nonzero anywhere, include stochastic field and normalise magnetisation every timestep
    if (any(gb_problem%temperature > 1.0e-15)) then
        gb_problem%includeThermal = .true.
        call displayGUIMessage( 'Including thermal noise' )
    else
        gb_problem%includeThermal = .false.
    end if
        
    if ( gb_problem%solver .eq. MicroMagSolverExplicit ) then
        !Go through several different applied fields and find the equilibrium solution for each of them
        nt_Hext = size(gb_problem%Hext, 1)          !The no. of applied fields to consider
        if (gb_problem%adaptiveHext) then
            ! Allocate one extra slot for the initial starting field and state
            nt_Hext = gb_problem%maxHextSteps + 1
        endif
    else if ( gb_problem%solver .eq. MicroMagSolverDynamic ) then
        !Simply do a time evolution as specified in the problem
        nt_Hext = 1

        if (gb_problem%includeThermal) then
            ! Calculate the prefactor for the thermal magnetic field
            t_step = gb_problem%t(nt) / nt    ! Timestep [s] (assumed constant)
            gamma = gb_problem%gamma          ! Gyromagnetic factor [m/(A*s)]
            alpha0 = gb_problem%alpha0        ! Damping constant [m/(A*s)]
            kB = 1.38064851 * 1e-23           ! The Boltzmann constant [J/K]
            ! Get cell volumes
            allocate( volCells(ntot) )
            if (gb_problem%grid%gridType /= gridTypeUniform) then
                volCells = gb_problem%grid%abc(:,1) * gb_problem%grid%abc(:,2) * gb_problem%grid%abc(:,3)
            else
                volCells = gb_problem%grid%dx * gb_problem%grid%dy * gb_problem%grid%dz
            end if
            ! Combined computation
            alphaGilbert = gamma/(2*alpha0) + 1/2 * sqrt((gamma/alpha0)**2 - 4)   ! Dimensionless Gilbert damping (second order polynomial equation)
            allocate( gb_problem%Tfact(ntot) )
            gb_problem%Tfact = sqrt(2*kB * gb_problem%temperature * alphaGilbert / (mu0 * gamma * volCells * gb_problem%Ms * t_step))
            !print *, "Finished calculating thermal prefactor"
            !print *, 'Tfact', gb_problem%Tfact(1)   ! Just for debugging
        else
            gb_problem%Tfact = 0
        end if
    endif
    
    allocate(M_out(3*ntot,nt,nt_Hext))   
    allocate(gb_solution%M_out(size(gb_problem%t),ntot,nt_Hext,3))
    !Allocate the arrays for the different fields
    !Only if these are to be returned are they saved at every time step
    if (gb_problem%useReturnHall .eq. useReturnHallTrue) then
        call displayGUIMessage( 'useReturnHall is True ' )
        allocate(gb_solution%H_exc(size(gb_problem%t),ntot,nt_Hext,3))
        allocate(gb_solution%H_ext(size(gb_problem%t),ntot,nt_Hext,3))
        allocate(gb_solution%H_dem(size(gb_problem%t),ntot,nt_Hext,3))
        allocate(gb_solution%H_ani(size(gb_problem%t),ntot,nt_Hext,3))  
    else
        call displayGUIMessage( 'useReturnHall is False ' )
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
    if ( gb_problem%dummy_run .eq. 1) then
       call displayGUIMessage( 'Performing dummy run - skipping time evolution ' )
    else 

        !CALL SYSTEM_CLOCK(c1)
    !$omp parallel default(shared)
    !$omp master
        
      if (gb_problem%adaptiveHext) then
          n_reject_max = max(100, 10 * gb_problem%maxHextSteps)
          call SolveAdaptiveHextLoop(fct, fct_thermal, cb_fct, M_out, ntot, nt, n_reject_max)
      else
          call SolveFixedHextLoop(fct, fct_thermal, cb_fct, M_out, ntot, nt, nt_Hext)
      endif

      !$omp end master
      !$omp end parallel

      if (gb_problem%adaptiveHext) then
          nt_Hext = gb_problem%nHextAccepted
      endif

      
        !CALL SYSTEM_CLOCK(c2)
        !call displayGUIMessage( 'Time simulation time integration:' )
        !write (prog_str,'(f10.3)') (c2 - c1)/rate
        !call displayGUIMessage( prog_str )


      !clean up
      deallocate(crossX,crossY,crossZ,HeffX,HeffY,HeffZ,HeffX2,HeffY2,HeffZ2, M_out)
      
      !clean up
      if (gb_problem%CV > 0) then
          deallocate(gb_solution%u1, gb_solution%u2, gb_solution%u3, gb_solution%u4, gb_solution%u5, gb_solution%u6)
      endif
      
      !clean-up
      stat = DftiFreeDescriptor(gb_problem%desc_hndl_FFT_M_H)
  end if

    #if USE_CUDA
      if ( gb_problem%useCuda .eqv. useCudaTrue ) then
          call cudaDestroy()
      endif
    #endif
    !Return the correct state
    sol = gb_solution
    prob = gb_problem    


    call trace%end( "SolveLandauLifshitzEquation", itimer=itimer, verbose=1 )
    end subroutine SolveLandauLifshitzEquation


    subroutine SolveFixedHextLoop(fct, fct_thermal, cb_fct, M_out, ntot, nt, nt_Hext)
    procedure(dydt_fct), pointer :: fct
    procedure(no_argument_fct), pointer :: fct_thermal
    procedure(callback_fct),pointer :: cb_fct
    real(DP),dimension(:,:,:),intent(inout) :: M_out
    integer,intent(in) :: ntot, nt, nt_Hext
    integer :: i
    character*(100) :: prog_str

      do i=1,nt_Hext
          !Applied field
          gb_solution%HextInd = i

          if (gb_problem%solver .eq. MicroMagSolverExplicit) then
              write(prog_str,'(A20, I5, A8, I5, A6, F6.2, A7)') 'External Field nr.: ', i, ' out of ', nt_Hext, ' i.e. ', real(i)/real(nt_Hext)*100,'% done'
              call displayGUIMessage( trim(prog_str) )
          endif

          ! This is where the LLG is actually integrated (short description by F. Durhuus)
          ! fct is a pointer to the dmdt_fct function which returns time derivatives of the normalised magnetic moments
          ! With m0 as initial condition, integrate from t(1) to t(nt) with effective field updated at each time coordinate t, then store result in t_out and M_out.
          ! t_conv is time coordinates where convergence is tested when computing equilibrium structure (explicit solver rather than dynamic). Also included in t_out.
          ! When thermal noise is included (includeThermal = .true.) the magnetisation is normalised each timestep to prevent thermal drift
          call MagTense_ODE( fct, gb_problem%t, gb_problem%m0, gb_solution%t_out, M_out(:,:,i), gb_problem%includeThermal, fct_thermal, cb_fct, &
                  gb_problem%setTimeDisplay, gb_problem%tol, gb_problem%thres_value, gb_problem%useCVODE, gb_problem%t_conv, gb_problem%conv_tol )

          !The initial state of the next solution is the previous solution result
          gb_problem%m0 = M_out(:,nt,i)

          !Store the solution
          gb_solution%M_out(:,:,i,1) =  transpose( M_out(1:ntot,:,i) )
          gb_solution%M_out(:,:,i,2) =  transpose( M_out((ntot+1):2*ntot,:,i) )
          gb_solution%M_out(:,:,i,3) =  transpose( M_out((2*ntot+1):3*ntot,:,i)  )

          call StoreHeffComponents ( gb_problem, gb_solution )
      enddo

      gb_problem%nHextAccepted = nt_Hext
    end subroutine SolveFixedHextLoop

    subroutine SolveAdaptiveHextLoop(fct, fct_thermal, cb_fct, M_out, ntot, nt, n_reject_max)
    procedure(dydt_fct), pointer :: fct
    procedure(no_argument_fct), pointer :: fct_thermal
    procedure(callback_fct),pointer :: cb_fct
    real(DP),dimension(:,:,:),intent(inout) :: M_out
    integer,intent(in) :: ntot, nt, n_reject_max
    integer :: i_acc, i_trial, n_reject, n_reject_total, ti
    real(DP) :: H_delta(3), H_dir(3), H_current(3), H_trial(3), H_remaining(3)
    real(DP) :: H_distance, remaining_distance, dH, dH_initial, dH_step, dH_initial_T, dH_T, dH_step_T, dM, m_parallel_before, m_parallel_trial
    real(DP),dimension(:),allocatable :: m_before, m_accepted, m_trial
    logical :: reject_step, reached_end
    character*(256) :: prog_str

      if (gb_problem%maxHextSteps <= 0) then
          call displayGUIMessage('Adaptive hysteresis requires maxHextSteps > 0')
          error stop 'Adaptive hysteresis requires maxHextSteps > 0'
      endif
      if (gb_problem%dH_initial <= 0.0_DP .or. gb_problem%dH_min <= 0.0_DP .or. gb_problem%dH_max <= 0.0_DP) then
          call displayGUIMessage('Adaptive hysteresis requires positive dH_initial, dH_min and dH_max')
          error stop 'Adaptive hysteresis requires positive dH_initial, dH_min and dH_max'
      endif
      if (gb_problem%dH_min > gb_problem%dH_max) then
          call displayGUIMessage('Adaptive hysteresis requires dH_min <= dH_max')
          error stop 'Adaptive hysteresis requires dH_min <= dH_max'
      endif

      H_delta = gb_problem%H_end - gb_problem%H_start
      H_distance = sqrt(sum(H_delta**2))
      if (H_distance <= tiny(1.0_DP)) then
          call displayGUIMessage('Adaptive hysteresis requires H_start and H_end to differ')
          error stop 'Adaptive hysteresis requires H_start and H_end to differ'
      endif

      H_dir = H_delta / H_distance
      H_current = gb_problem%H_start
      dH_initial = gb_problem%dH_initial
      dH = dH_initial
      allocate(m_before(3*ntot), m_accepted(3*ntot), m_trial(3*ntot))
      m_accepted = gb_problem%m0

    ! Save the starting field and magnetisation state (store at first time index)
    i_acc = 1
    gb_solution%HextInd = i_acc
    gb_problem%Hext(i_acc,1) = 0.0_DP
    gb_problem%Hext(i_acc,2:4) = H_current
    ! Store the initial condition for all time indices so the returned
    ! `gb_solution%M_out` contains a valid time-series for the initial field.
    do ti = 1, size(gb_problem%t)
        gb_solution%M_out(ti,:,i_acc,1) = m_accepted(1:ntot)
        gb_solution%M_out(ti,:,i_acc,2) = m_accepted(ntot+1:2*ntot)
        gb_solution%M_out(ti,:,i_acc,3) = m_accepted(2*ntot+1:3*ntot)
    end do
    call StoreHeffComponents ( gb_problem, gb_solution )

    n_reject = 0
      n_reject_total = 0
      reached_end = .false.

      do while (.not. reached_end .and. i_acc < gb_problem%maxHextSteps + 1)
          remaining_distance = dot_product(gb_problem%H_end - H_current, H_dir)
          if (remaining_distance <= max(1.0e-12_DP * H_distance, tiny(1.0_DP))) exit

          i_trial = i_acc + 1
          H_trial = H_current + min(dH, remaining_distance) * H_dir
          if (remaining_distance <= dH) H_trial = gb_problem%H_end
          dH_step = sqrt(sum((H_trial - H_current)**2))
          dH_initial_T = mu0 * dH_initial
          dH_T = mu0 * dH
          dH_step_T = mu0 * dH_step

          m_before = m_accepted
          gb_problem%m0 = m_before
          gb_solution%HextInd = i_trial
          gb_problem%Hext(i_trial,1) = 0.0_DP
          gb_problem%Hext(i_trial,2:4) = H_trial

          if (n_reject == 0) then
              write(prog_str,'(A31,I5,A10,F8.2,A7)') 'Adaptive External Field nr.: ', i_trial, '   progress ', real(100.0_DP*dot_product(H_trial - gb_problem%H_start, H_dir)/H_distance, DP), '% done'
              call displayGUIMessage( trim(prog_str) )
              write(prog_str,'(A22,3F10.6,A2)') '      mu0 H [T] = ', mu0*H_trial(1), mu0*H_trial(2), mu0*H_trial(3), ' '
              call displayGUIMessage( trim(prog_str) )
              write(prog_str,'(A27,F10.6,A9,F10.6,A13,F10.6)') '      dH_initial [T] = ', dH_initial_T, ' dH [T] = ', dH_T, ' actual [T] = ', dH_step_T
              call displayGUIMessage( trim(prog_str) )
          endif

          call MagTense_ODE( fct, gb_problem%t, gb_problem%m0, gb_solution%t_out, M_out(:,:,i_trial), gb_problem%includeThermal, fct_thermal, cb_fct, &
                  gb_problem%setTimeDisplay, gb_problem%tol, gb_problem%thres_value, gb_problem%useCVODE, gb_problem%t_conv, gb_problem%conv_tol )

          m_trial = M_out(:,nt,i_trial)
          dM = adaptiveStepMetric(m_before, m_trial, ntot)

          m_parallel_before = meanMagnetisationAlongField(m_before, ntot, H_dir)
          m_parallel_trial = meanMagnetisationAlongField(m_trial, ntot, H_dir)
          reject_step = .false.
          if (dM > gb_problem%dM_reject .and. dH > gb_problem%dH_min) reject_step = .true.
          if (gb_problem%use_switch_refine) then
              if (m_parallel_before * m_parallel_trial < 0.0_DP .and. dH > gb_problem%switch_refine_dH .and. dH > gb_problem%dH_min) reject_step = .true.
          endif

          if (reject_step) then
              if (n_reject > 0 .and. dH <= gb_problem%dH_min) then
                  reject_step = .false.
              else
                  dH = max(dH * gb_problem%dH_shrink, gb_problem%dH_min)
                  dH_T = mu0 * dH
                  write(prog_str,'(A27,F10.6,A4)') '   Retrying with lower dH = ', dH_T, ' T'
                  call displayGUIMessage( trim(prog_str) )
              endif
          endif

          if (reject_step) then
              gb_problem%m0 = m_before
              n_reject = n_reject + 1
              n_reject_total = n_reject_total + 1
              if (n_reject_total > n_reject_max) then
                  call displayGUIMessage('Adaptive hysteresis stopped after too many rejected field steps')
                  error stop 'Adaptive hysteresis stopped after too many rejected field steps'
              endif
              cycle
          endif

          if (dM > gb_problem%dM_reject .and. dH <= gb_problem%dH_min) then
              call displayGUIMessage('Adaptive hysteresis accepting a large magnetisation change at dH_min')
          endif

          i_acc = i_trial
          n_reject = 0
          H_current = H_trial
          m_accepted = m_trial
          gb_problem%m0 = m_accepted
          gb_solution%HextInd = i_acc
          gb_problem%Hext(i_acc,1) = 0.0_DP
          gb_problem%Hext(i_acc,2:4) = H_current

          gb_solution%M_out(:,:,i_acc,1) =  transpose( M_out(1:ntot,:,i_trial) )
          gb_solution%M_out(:,:,i_acc,2) =  transpose( M_out((ntot+1):2*ntot,:,i_trial) )
          gb_solution%M_out(:,:,i_acc,3) =  transpose( M_out((2*ntot+1):3*ntot,:,i_trial)  )

          call StoreHeffComponents ( gb_problem, gb_solution )

          H_remaining = gb_problem%H_end - H_current
          reached_end = sqrt(sum(H_remaining**2)) <= max(1.0e-10_DP * H_distance, 1.0e-9_DP)

          if (dM < gb_problem%dM_min) then
              dH = min(dH * gb_problem%dH_grow, gb_problem%dH_max)
          else if (dM > gb_problem%dM_target) then
              dH = max(dH * gb_problem%dH_shrink, gb_problem%dH_min)
          endif

          
          
      enddo

      gb_problem%nHextAccepted = i_acc
      if (.not. reached_end) then
          call displayGUIMessage('Adaptive hysteresis reached maxHextSteps before H_end')
      endif
      deallocate(m_before, m_accepted, m_trial)
    end subroutine SolveAdaptiveHextLoop

    real(DP) function adaptiveStepMetric(m_before, m_trial, ntot) result(dM_rms)
    real(DP),dimension(:),intent(in) :: m_before, m_trial
    integer,intent(in) :: ntot
    real(DP) :: M_sum_before(3), M_sum_trial(3)
      M_sum_before(1) = sum(m_before(1:ntot)) / real(ntot, DP)
      M_sum_before(2) = sum(m_before(ntot+1:2*ntot)) / real(ntot, DP)
      M_sum_before(3) = sum(m_before(2*ntot+1:3*ntot)) / real(ntot, DP)
      M_sum_trial(1)  = sum(m_trial(1:ntot)) / real(ntot, DP)
      M_sum_trial(2)  = sum(m_trial(ntot+1:2*ntot)) / real(ntot, DP)
      M_sum_trial(3)  = sum(m_trial(2*ntot+1:3*ntot)) / real(ntot, DP)
      dM_rms = sqrt(sum((M_sum_trial - M_sum_before)**2))
    end function adaptiveStepMetric

    real(DP) function meanMagnetisationAlongField(m, ntot, H_dir) result(M_parallel)
    real(DP),dimension(:),intent(in) :: m
    integer,intent(in) :: ntot
    real(DP),dimension(3),intent(in) :: H_dir
    real(DP) :: M_mean(3)
      M_mean(1) = sum(m(1:ntot)) / real(ntot, DP)
      M_mean(2) = sum(m(ntot+1:2*ntot)) / real(ntot, DP)
      M_mean(3) = sum(m(2*ntot+1:3*ntot)) / real(ntot, DP)
      M_parallel = dot_product(M_mean, H_dir)
    end function meanMagnetisationAlongField
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
        real(DP) :: mx_mean, my_mean, mz_mean, volume_total
        integer :: i
        integer, save :: itimer = 0, itimer_wait = 0
        call trace%begin( "dmdt_fct", itimer=itimer, verbose=1 )
        !------------------------------------------
        ntot = gb_problem%grid%nx * gb_problem%grid%ny * gb_problem%grid%nz

        if ( .not. allocated(crossX) ) then
            allocate( crossX(ntot), crossY(ntot), crossZ(ntot) )
            allocate( HeffX(ntot), HeffY(ntot), HeffZ(ntot) )
            allocate( HeffX2(ntot), HeffY2(ntot), HeffZ2(ntot) )
        endif

        !------------- Update magnetisation (load from m)-------------
        gb_solution%Mx = m(1:ntot)
        gb_solution%My = m(ntot+1:2*ntot)
        gb_solution%Mz = m(2*ntot+1:3*ntot)
        !-------------------------------------------------------------

        !------------ add exchange term -----------------------------
        !$omp task untied default(shared)
        call updateExchangeTerms( gb_problem, gb_solution )
        !$omp end task
        !-------------------------------------------------------------
        !-------------- add external field ---------------------------
        !$omp task untied default(shared)
        call updateExternalField( gb_problem, gb_solution, t )
        !$omp end task
        !-------------------------------------------------------------
        !-------------- add anisotropy term --------------------------
        !$omp task untied default(shared)
        call updateAnisotropy(  gb_problem, gb_solution )
        !$omp end task
        !-------------------------------------------------------------
        !-------------- add demagnetization field --------------------
                ! NOTE - we don't run updateDemagfield as a task as this is the most expensive part
                ! instead, we run this on the main thread and then.
#if USE_FMM3D
        !-------------------- if use_fmm then use FMM otherwise shortcircuit to normal demag field --------------
        if ( gb_problem%use_fmm)  then
          call updateDemagfieldFMM( gb_problem, gb_solution )
        else
           call updateDemagfield( gb_problem, gb_solution )
        end if
        !--------------------------------------------------------------------------------------------------------------
#else
         call updateDemagfield( gb_problem, gb_solution )
#endif
        !-------------------------------------------------------------
        !------------------ wait for all task to finish before calculating the effective field and the dm/dt ----------------
        call trace%begin( "dmdt_fct_taskwait", itimer=itimer_wait, verbose=1 )
        !NOTE:      Not needed when running FMM, since it already contains a taskwait
        !$omp taskwait
        call trace%end( "dmdt_fct_taskwait", itimer=itimer_wait, verbose=1 )
        !--------------------------------------------------------------------------------------------------------------------
        !--------------- combine to get effective field, Heff -------------
        HeffX = gb_solution%HhX + gb_solution%HjX + gb_solution%HmX + gb_solution%HkX + gb_solution%HtX
        HeffY = gb_solution%HhY + gb_solution%HjY + gb_solution%HmY + gb_solution%HkY + gb_solution%HtX
        HeffZ = gb_solution%HhZ + gb_solution%HjZ + gb_solution%HmZ + gb_solution%HkZ + gb_solution%HtX
        !-------------------------------------------------------------
        !--------------- calculate precession term: m x Heff -------------
        crossX = 1.0_DP * ( gb_solution%My * HeffZ - gb_solution%Mz * HeffY )
        crossY = 1.0_DP * ( gb_solution%Mz * HeffX - gb_solution%Mx * HeffZ )
        crossZ = 1.0_DP * ( gb_solution%Mx * HeffY - gb_solution%My * HeffX )
        !-------------------------------------------------------------
        !--------------- calculate damping term: m x (m x Heff) -------------
        HeffX2 = gb_solution%My * crossZ - gb_solution%Mz * crossY
        HeffY2 = gb_solution%Mz * crossX - gb_solution%Mx * crossZ
        HeffZ2 = gb_solution%Mx * crossY - gb_solution%My * crossX
        !-------------------------------------------------------------
        !--------------- assemble dm/dt -------------------------------------
        dmdt(1:ntot)               = -gb_problem%gamma * crossX - alpha(t,gb_problem) * HeffX2
        dmdt(ntot+1:2*ntot)        = -gb_problem%gamma * crossY - alpha(t,gb_problem) * HeffY2
        dmdt(2*ntot+1:3*ntot)      = -gb_problem%gamma * crossZ - alpha(t,gb_problem) * HeffZ2
        !---------------------------------------------------------------------
        call trace%end("dmdt_fct", itimer=itimer, verbose=1 )
    end subroutine dmdt_fct

    
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
    integer, save :: itimer = 0
    call trace%begin( "StoreHeffComponents", itimer=itimer, verbose=1 )
    
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
        call updateExternalField( gb_problem, gb_solution, gb_solution%t_out(j) )
        !Anisotropy term
        call updateAnisotropy(  gb_problem, gb_solution )
        !Demag. field
#if USE_FMM3D
        if ( gb_problem%use_fmm)  then
          call updateDemagfieldFMM( gb_problem, gb_solution )
        else
           call updateDemagfield( gb_problem, gb_solution )
        end if
#else
        call updateDemagfield( gb_problem, gb_solution )
#endif

    
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
    end do
    
    call trace%end( "StoreHeffComponents", itimer=itimer, verbose=1 )
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
    integer, save :: itimer = 0
    call trace%begin( "initializeSolution", itimer=itimer, verbose=1 )

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

        !Thermal field
        allocate( solution%HtX(ntot), solution%HtY(ntot), solution%HtZ(ntot) )
        solution%HtX(:) = 0.
        solution%HtY(:) = 0.
        solution%HtZ(:) = 0.
        
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
    
    call trace%end( "initializeSolution", itimer=itimer, verbose=1 )
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
    integer, save :: itimer = 0 

    call trace%begin( "updateExchangeTerms", itimer=itimer, verbose=1 )
    
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    descr%mode = SPARSE_FILL_MODE_FULL
    descr%diag = SPARSE_DIAG_NON_UNIT
    
    alpha = 2.! * solution%Jfact
    beta = 0.
    
    ntot = problem%grid%nx * problem%grid%ny * problem%grid%nz
    allocate(temp(ntot))
    !Effective field in the X-direction. Note that the scalar alpha is multiplied on from the left, such that
    !y = alpha * (A_exch * Mx )
    stat = mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%A_exch, descr, solution%Mx, beta, temp )
    solution%HjX = temp * problem%Jfact
    
    !Effective field in the Y-direction
    stat = mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%A_exch, descr, solution%My, beta, temp )
    solution%HjY = temp * problem%Jfact
    
    !Effective field in the Z-direction
    stat = mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%A_exch, descr, solution%Mz, beta, temp )
    solution%HjZ = temp * problem%Jfact
    
    deallocate(temp)
    

    call trace%end( "updateExchangeTerms", itimer=itimer, verbose=1 )
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
    integer, save :: itimer = 0 

    call trace%begin( "updateExternalField", itimer=itimer, verbose=1 )
    
    if ( problem%solver .eq. MicroMagSolverExplicit ) then
         !Assume the field to be constant in time (we are finding the equilibrium solution at a given applied field)
        solution%HhX = problem%Hext(solution%HextInd,2)
        solution%HhY = problem%Hext(solution%HextInd,3)
        solution%HhZ = problem%Hext(solution%HextInd,4)

    elseif ( problem%solver .eq. MicroMagSolverDynamic ) then
        
        !Interpolate to get the applied field at time t
        call interp1_MagTense( problem%Hext(:,1), problem%Hext(:,2), t, size(problem%Hext(:,1)), HextX )
        call interp1_MagTense( problem%Hext(:,1), problem%Hext(:,3), t, size(problem%Hext(:,1)), HextY )
        call interp1_MagTense( problem%Hext(:,1), problem%Hext(:,4), t, size(problem%Hext(:,1)), HextZ )

        solution%HhX = HextX
        solution%HhY = HextY
        solution%HhZ = HextZ
        
    elseif ( problem%solver .eq. MicroMagSolverImplicit ) then
        !not implemented yet
    endif
    
    call trace%end( "updateExternalField", itimer=itimer, verbose=1 )

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
    integer, save :: itimer = 0, itimer1 = 0, itimer2 = 0

    call trace%begin( "updateAnisotropy", itimer=itimer, verbose=1 )
    
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
    
    if (any(problem%K0 .ne. 0)) then
        solution%Hkx = 2. * problem%Kfact * ( problem%Axx * solution%Mx + problem%Axy * solution%My + problem%Axz * solution%Mz )
        solution%Hky = 2. * problem%Kfact * ( problem%Axy * solution%Mx + problem%Ayy * solution%My + problem%Ayz * solution%Mz )
        solution%Hkz = 2. * problem%Kfact * ( problem%Axz * solution%Mx + problem%Ayz * solution%My + problem%Azz * solution%Mz )
    else 
        call trace%begin( "updateAnisotropy_alloc", itimer=itimer1, verbose=2 )
        allocate(Mx_rot(size(solution%Mx)), My_rot(size(solution%My)), Mz_rot(size(solution%Mz)))
        allocate(Hkx_rot(size(Mx_rot)), Hky_rot(size(My_rot)), Hkz_rot(size(Mz_rot)))
        call trace%end( "updateAnisotropy_alloc", itimer=itimer1, verbose=2 )
        call trace%begin( "updateAnisotropy_calc", itimer=itimer2, verbose=2 )
        Mx_rot = problem%CrystalAxis(:,1,1)*solution%Mx(:) + problem%CrystalAxis(:,1,2)*solution%My(:) + problem%CrystalAxis(:,1,3)*solution%Mz(:)
        My_rot = problem%CrystalAxis(:,2,1)*solution%Mx(:) + problem%CrystalAxis(:,2,2)*solution%My(:) + problem%CrystalAxis(:,2,3)*solution%Mz(:)
        Mz_rot = problem%CrystalAxis(:,3,1)*solution%Mx(:) + problem%CrystalAxis(:,3,2)*solution%My(:) + problem%CrystalAxis(:,3,3)*solution%Mz(:)
        
        Hkx_rot = -2*Mx_rot(:)*(problem%Kfact_arr(:,1,1) + 2*problem%Kfact_arr(:,2,1)*Mx_rot(:)**2 + problem%Kfact_arr(:,3,3)*My_rot(:)**2 + problem%Kfact_arr(:,3,2)*Mz_rot(:)**2 + 3*problem%Kfact_arr(:,4,1)*Mx_rot(:)**4 + problem%Kfact_arr(:,5,2)*My_rot(:)**4 + problem%Kfact_arr(:,5,3)*Mz_rot(:)**4 + 2*problem%Kfact_arr(:,5,1)*Mx_rot(:)**2*My_rot(:)**2 + 2*problem%Kfact_arr(:,5,1)*Mx_rot(:)**2*Mz_rot(:)**2 + problem%Kfact_arr(:,6,1)*My_rot(:)**2*Mz_rot(:)**2 )
        Hky_rot = -2*My_rot(:)*(problem%Kfact_arr(:,1,2) + 2*problem%Kfact_arr(:,2,2)*My_rot(:)**2 + problem%Kfact_arr(:,3,3)*Mx_rot(:)**2 + problem%Kfact_arr(:,3,1)*Mz_rot(:)**2 + 3*problem%Kfact_arr(:,4,2)*My_rot(:)**4 + problem%Kfact_arr(:,5,1)*Mx_rot(:)**4 + problem%Kfact_arr(:,5,3)*Mz_rot(:)**4 + 2*problem%Kfact_arr(:,5,2)*Mx_rot(:)**2*My_rot(:)**2 + 2*problem%Kfact_arr(:,5,2)*My_rot(:)**2*Mz_rot(:)**2 + problem%Kfact_arr(:,6,1)*Mx_rot(:)**2*Mz_rot(:)**2 )
        Hkz_rot = -2*Mz_rot(:)*(problem%Kfact_arr(:,1,3) + 2*problem%Kfact_arr(:,2,3)*Mz_rot(:)**2 + problem%Kfact_arr(:,3,2)*Mx_rot(:)**2 + problem%Kfact_arr(:,3,1)*My_rot(:)**2 + 3*problem%Kfact_arr(:,4,3)*Mz_rot(:)**4 + problem%Kfact_arr(:,5,1)*Mx_rot(:)**4 + problem%Kfact_arr(:,5,2)*My_rot(:)**4 + 2*problem%Kfact_arr(:,5,3)*Mx_rot(:)**2*Mz_rot(:)**2 + 2*problem%Kfact_arr(:,5,3)*My_rot(:)**2*Mz_rot(:)**2 + problem%Kfact_arr(:,6,1)*Mx_rot(:)**2*My_rot(:)**2 )
        call trace%end( "updateAnisotropy_calc", itimer=itimer2, verbose=2 )

        solution%Hkx(:) = - problem%CrystalAxis(:,1,1)*Hkx_rot(:) - problem%CrystalAxis(:,2,1)*Hky_rot(:) - problem%CrystalAxis(:,3,1)*Hkz_rot(:)
        solution%Hky(:) = - problem%CrystalAxis(:,1,2)*Hkx_rot(:) - problem%CrystalAxis(:,2,2)*Hky_rot(:) - problem%CrystalAxis(:,3,2)*Hkz_rot(:)
        solution%Hkz(:) = - problem%CrystalAxis(:,1,3)*Hkx_rot(:) - problem%CrystalAxis(:,2,3)*Hky_rot(:) - problem%CrystalAxis(:,3,3)*Hkz_rot(:)
        
        deallocate(Mx_rot, My_rot, Mz_rot, Hkx_rot, Hky_rot, Hkz_rot)
    end if 
    call trace%end( "updateAnisotropy", itimer=itimer, verbose=1 )

    end subroutine updateAnisotropy    

    !>-----------------------------------------
    !> @author Frederik L. Durhuus, fladu@dtu.dk, DTU, 2026
    !> Lightly edited function from ChatGPT
    !> Fills an array, z, with Gaussian distributed random numbers of mean 0 and standard deviation 1
    !> Uses the Marsaglia polar method to map uniformly distributed random numbers onto a Gaussian distribution
    !> TODO : Put this in some utilities module rather than the core script
    !>-----------------------------------------
    subroutine GenerateGaussianArray(z)
    real, intent(out) :: z(:)

    integer :: n, filled, batch, i
    real, allocatable :: u(:), v(:), s(:), factor(:)
    logical, allocatable :: mask(:)

    n = size(z)
    filled = 0

    ! Batch size: typically similar to output size
    ! Lower limit to avoid efficiency issues and upper limit to avoid excessively large arrays
    batch = min(max(64, n), 8192)

    allocate(u(batch), v(batch), s(batch), factor(batch), mask(batch))

    do while (filled < n)

        ! Generate uniform random numbers
        call random_number(u)
        call random_number(v)

        ! Transform to [-1,1]
        u = 2.0*u - 1.0
        v = 2.0*v - 1.0

        ! Compute radius
        s = u*u + v*v

        ! Acceptance mask
        mask = (s < 1.0) .and. (s > 0.0)

        where(mask)
            factor = sqrt(-2.0*log(s)/s)
        end where

        ! Store accepted numbers
        do i = 1, batch
            if (mask(i)) then
                if (filled < n) then
                    filled = filled + 1
                    z(filled) = u(i)*factor(i)
                end if

                if (filled < n) then
                    filled = filled + 1
                    z(filled) = v(i)*factor(i)
                end if
            end if

            if (filled >= n) exit
        end do

    end do

    deallocate(u,v,s,factor,mask)

    end subroutine GenerateGaussianArray

    !>-----------------------------------------
    !> @author Frederik L. Durhuus, fladu@dtu.dk, DTU, 2026
    !> Calculates and returns the thermal field
    !>-----------------------------------------
    subroutine updateThermalField( problem, solution )
    type(MicroMagProblem),intent(in) :: problem         !> Problem data structure
    type(MicroMagSolution),intent(inout) :: solution    !> Solution data structure

    integer :: ntot
    real,dimension(:),allocatable :: GaussianArray

    ntot = problem%grid%nx * problem%grid%ny * problem%grid%nz
    allocate(GaussianArray(3*ntot))
    call GenerateGaussianArray(GaussianArray)

    solution%HtX = GaussianArray(1:ntot) * problem%Tfact
    solution%HtY = GaussianArray(ntot+1:2*ntot) * problem%Tfact
    solution%HtZ = GaussianArray(2*ntot+1:3*ntot) * problem%Tfact

    end subroutine updateThermalField

    !> Wrapper function so we can construct a pointer for updateThermalField
    !>  and future thermal update functions that does not require any arguments
    subroutine updateThermal_wrapper()
        call updateThermalField( gb_problem, gb_solution)
    end subroutine updateThermal_wrapper


subroutine updateDemagfieldFMM(problem, solution)
  implicit none
  !------------------ Arguments ------------------
  type(MicroMagProblem),  intent(inout) :: problem
  type(MicroMagSolution), intent(inout) :: solution

#if USE_FMM3D
  !------------------ Timing ---------------------
  !------------------ Locals ---------------------
  integer :: ntot, i
  real(DP) :: fourpi
  integer :: nd, ier
  double precision , contiguous, pointer :: source(:,:), dipvec(:,:,:)
  double precision , contiguous, pointer :: grad(:,:,:)
  double precision :: vol_i
  double precision :: mx, my, mz
  double precision :: eps_fmm

  class(FMM3DTree), pointer :: fmm_tree
  logical :: built_tree
  integer, save :: itimer = 0, itimer_wait = 0 
  !------------------------------------------------
    
  if ( problem%useDemag .eq. useDemagFalse ) then
        solution%HmX(:) = 0.
        solution%HmY(:) = 0.
        solution%HmZ(:) = 0.
        return
  endif
  
  call trace%begin( "updateDemagfieldFMM", itimer=itimer, verbose=1 )

  fourpi  = 12.566370614359172D0
  ntot = size(problem%grid%pts, dim=1)
  ier = 0

  !----------- cast Mx,My,Mz to single precision -----------
  solution%Mx_s = real(solution%Mx, SP)
  solution%My_s = real(solution%My, SP)
  solution%Mz_s = real(solution%Mz, SP)
  !--------------------------------------------------------------


  !-------------- reset field ---------------------------------------------
solution%HmX = 0.0_SP
solution%HmY = 0.0_SP
solution%HmZ = 0.0_SP
  !-------------------------------------------------------------------------

    !------------- add correction from neighbouring tiles ---------------------
    !$omp task untied default(shared)
    call add_near_field(problem, solution)
    !$omp end task
    !-------------------------------------------------------------------------


  !------------------ Allocate FMM work arrays ------------------
    built_tree = .false.
    if (.not. associated(solution%fmm_tree)) then
        allocate(solution%fmm_tree)
        built_tree = .true.
        allocate(source(3, ntot))
    end if
    fmm_tree => solution%fmm_tree

  allocate(dipvec(1, 3, ntot))
  allocate(grad(1, 3, ntot))
  grad(:,:,:) = 0.0_DP 
  !------------------------------------------------------------
  !------------------ Pack sources and dipoles ------------------
    if (built_tree) then
        do i = 1, ntot
            !------------------ positions -> FMM (3, N) ------------------
            source(1,i) = real(problem%grid%pts(i,1), DP)
            source(2,i) = real(problem%grid%pts(i,2), DP)
            source(3,i) = real(problem%grid%pts(i,3), DP)
            !------------------------------------------------------------
        end do
    end if
  do i = 1, ntot
    !------------- get cell volume ------------------------------
    ! TODO - implement non-uniform and non-cubic grids here
    if (problem%grid%gridType /= gridTypeUniform) then
      vol_i = real(problem%grid%abc(i,1), DP) * &
              real(problem%grid%abc(i,2), DP) * &
              real(problem%grid%abc(i,3), DP)
    else
    vol_i = real(problem%grid%dx, DP) * &
            real(problem%grid%dy, DP) * &
            real(problem%grid%dz, DP)
    end if
    !------------------------------------------------------------
    !--------- conert to dipole moment m = M * ΔV --------------
    mx = solution%Mx(i) * vol_i * problem%Ms(i)
    my = solution%My(i) * vol_i * problem%Ms(i)
    mz = solution%Mz(i) * vol_i * problem%Ms(i)
    !------------------------------------------------------------
    !-------------- pack dipole vector ---------------------------
    dipvec(1,1,i) = mx 
    dipvec(1,2,i) = my
    dipvec(1,3,i) = mz
    !------------------------------------------------------------ 
  end do
    !--------------------------------------------------------------
   !------------------ Call FMM (sources->sources) ------------------
   nd = 1
   fmm_tree%nterms_in = problem%fmm_nterms
   call fmm_tree%build_tree( source, problem%fmm_eps, problem%fmm_cells_per_node , ier, problem%ifunif, problem%nlmin, problem%nlmax)
   !------- only run if number of boxes > 9 -----------
   !--- NOTE - if nboxes <= 9 then we have all-to-all and the is no need for FMM ---
   if (fmm_tree%nboxes > 9) then
       call fmm_tree%make_and_eval(dipvec, grad)
   end if
   !--------------------------------------------------

  if (ier /= 0) then
    write(*,*) 'FMM3D returned error code in updateDemagfieldFMM: ier =', ier
    stop " FMM3D error"
  end if
  !-----------------------------------------------------------------
  !------------------ Map grad -> H (and to single) ----------------
  ! include factor 4pi to match Magtense units
   if (fmm_tree%nboxes > 9) then 
      !$omp critical (solution_update)
        solution%HmX = solution%HmX - real( grad(1,1,:) / fourpi, SP )
        solution%HmY = solution%HmY - real( grad(1,2,:) / fourpi, SP )
        solution%HmZ = solution%HmZ - real( grad(1,3,:) / fourpi, SP ) 
      !$omp end critical (solution_update)
  end if
  !-----------------------------------------------------------------

  !------------------- Wait for near-field task to complete ------------------
  call trace%begin( "updateDemagfieldFMM_taskwait", itimer=itimer_wait, verbose=1 )
  !$omp taskwait
  call trace%end( "updateDemagfieldFMM_taskwait", itimer=itimer_wait, verbose=1 )
  !---------------------------------------------------------------------------


    !------------------ Cleanup -------------------------------------
  if (.not. fmm_tree%keep_tree) then
      call fmm_tree%dealloc()
      deallocate( fmm_tree )
      nullify(solution%fmm_tree)
  end if
  if (built_tree) then
      deallocate( source )
  end if
  deallocate(dipvec, grad)
    !--------------------------------------------------------------
  !------------------ Optional field noise (CV) --------------------
  if (problem%CV > 0) then
      solution%HmX = solution%HmX + solution%HmX*problem%CV*sqrt(-2d0*log(solution%u1))*cos(2d0*pi*solution%u2)
      solution%HmY = solution%HmY + solution%HmY*problem%CV*sqrt(-2d0*log(solution%u3))*cos(2d0*pi*solution%u4)
      solution%HmZ = solution%HmZ + solution%HmZ*problem%CV*sqrt(-2d0*log(solution%u5))*cos(2d0*pi*solution%u6)
  end if
    call trace%end( "updateDemagfieldFMM", itimer=itimer, verbose=1 )
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
    type(MicroMagProblem),intent(in) :: problem         !> Problem data structure    
    type(MicroMagSolution),intent(inout) :: solution    !> Solution data structure
    integer :: stat,ntot,i
    type(matrix_descr) :: descr
    real(SP) :: pref,alpha,beta
    complex(kind=4) :: alpha_c, beta_c
    real(SP), dimension(:), allocatable :: temp, mx, my, mz
    character*(100) :: prog_str 
    integer, save :: itimer = 0
    real(DP),dimension(3) :: Mavg           !> Average magnetisation
    
    if ( problem%useDemag .eq. useDemagFalse ) then
        solution%HmX(:) = 0.
        solution%HmY(:) = 0.
        solution%HmZ(:) = 0.
        return
    endif
    
    call trace%begin( "updateDemagfield", itimer=itimer, verbose=1 )
    
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    descr%mode = SPARSE_FILL_MODE_FULL
    descr%diag = SPARSE_DIAG_NON_UNIT
    
    ntot = problem%grid%nx * problem%grid%ny * problem%grid%nz
    
    ! Convert the magnetization to single precision before the demag calculation
    solution%Mx_s = real(solution%Mx, SP)
    solution%My_s = real(solution%My, SP)
    solution%Mz_s = real(solution%Mz, SP)

    ! Get average magnetisation
    ! Average across micromagnetic cells, then multiply by the volume fraction of the cells so
    ! non-magnetic regions are included in the volume average. Assumes all cells have the same volume.
    Mavg(1) = sum(solution%Mx)/ntot * problem%VfracOcc
    Mavg(2) = sum(solution%My)/ntot * problem%VfracOcc
    Mavg(3) = sum(solution%Mz)/ntot * problem%VfracOcc
    
    if ( ( problem%demag_approximation .eq. DemagApproximationThreshold ) .or. ( problem%demag_approximation .eq. DemagApproximationThresholdFraction ) ) then
        if ( problem%useCuda .eq. useCudaFalse ) then
            !Do the matrix multiplications using sparse matrices
            alpha = 1.0
            beta = 0.
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(1)%A, descr, solution%Mx_s, beta, solution%HmX )
            beta = 1.0
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(2)%A, descr, solution%My_s, beta, solution%HmX )
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(3)%A, descr, solution%Mz_s, beta, solution%HmX )
        
           
            beta = 0.
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(2)%A, descr, solution%Mx_s, beta, solution%HmY )
            beta = 1.0
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(4)%A, descr, solution%My_s, beta, solution%HmY )
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(5)%A, descr, solution%Mz_s, beta, solution%HmY )
          
            beta = 0.
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(3)%A, descr, solution%Mx_s, beta, solution%HmZ )
            beta = 1.0
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(5)%A, descr, solution%My_s, beta, solution%HmZ )
            stat = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_s(6)%A, descr, solution%Mz_s, beta, solution%HmZ )
        
        else
#if USE_CUDA
            !Do the sparse matrix multiplication using CUDA
            !call displayGUIMessage( ' runing cuda demag sparse matrix multiplication' )
            pref = sngl(1 )                              
            call cudaMatrVecMult_sparse( solution%Mx_s, solution%My_s, solution%Mz_s, solution%HmX, solution%HmY, solution%HmZ, pref )
#endif
        endif
        
    elseif ( ( problem%demag_approximation .eq. DemagApproximationFFTThreshold ) .or. ( problem%demag_approximation .eq. DemagApproximationFFTThresholdFraction ) ) then

        if ( problem%useCuda .eq. useCudaFalse ) then
            !fourier transform Mx, My and Mz
            ntot = problem%grid%nx * problem%grid%ny * problem%grid%nz
                

            !======= Needs proper Ms correction 
            print *, "WARNING: not properly corrected for Ms"
            !==========================================

            !Convert to complex format
            do i=1,ntot
                solution%Mx_FT(i) = cmplx( solution%Mx_s(i), 0. )
                solution%My_FT(i) = cmplx( solution%My_s(i), 0. )
                solution%Mz_FT(i) = cmplx( solution%Mz_s(i), 0. )
            end do
        
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
            !solution%HmX = -problem%Mfact * real(solution%HmX_c)
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
            !solution%HmY = -problem%Mfact * real(solution%HmY_c)
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
           ! solution%HmZ = -problem%Mfact * real(solution%HmZ_c)
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
            !======== Hmx ===========
            alpha = 1.
            beta = 0.0
            call gemv( problem%Kxx, solution%Mx_s, solution%HmX, alpha, beta )        
            beta = 1.0
            call gemv( problem%Kxy, solution%My_s, solution%HmX, alpha, beta )
            call gemv( problem%Kxz, solution%Mz_s, solution%HmX, alpha, beta )
            !======= Hmy ==========
            beta = 0.0
            call gemv( problem%Kxy, solution%Mx_s, solution%HmY, alpha, beta )
            beta = 1.0
            call gemv( problem%Kyy, solution%My_s, solution%HmY, alpha, beta )
            call gemv( problem%Kyz, solution%Mz_s, solution%HmY, alpha, beta )
            !======= Hmz =========
            beta = 0.0
            call gemv( problem%Kxz, solution%Mx_s, solution%HmZ, alpha, beta )
            beta = 1.0
            call gemv( problem%Kyz, solution%My_s, solution%HmZ, alpha, beta )
            call gemv( problem%Kzz, solution%Mz_s, solution%HmZ, alpha, beta )
        else
            pref = sngl(1)
#if USE_CUDA
            call cudaMatrVecMult( solution%Mx_s, solution%My_s, solution%Mz_s, solution%HmX, solution%HmY, solution%HmZ, pref )
#endif
        endif 
    endif
        
    !Apply shape correction
    solution%HmX = solution%HmX + problem%Kxx_shape * Mavg(1) * problem%Mfact + problem%Kxy_shape * Mavg(2) * problem%Mfact + problem%Kxz_shape * Mavg(3) * problem%Mfact
    solution%HmY = solution%HmY + problem%Kxy_shape * Mavg(1) * problem%Mfact + problem%Kyy_shape * Mavg(2) * problem%Mfact + problem%Kyz_shape * Mavg(3) * problem%Mfact
    solution%HmZ = solution%HmZ + problem%Kxz_shape * Mavg(1) * problem%Mfact + problem%Kyz_shape * Mavg(2) * problem%Mfact + problem%Kzz_shape * Mavg(3) * problem%Mfact

    if (problem%CV > 0) then
        solution%HmX = solution%HmX + solution%HmX*problem%CV*sqrt(-2d0*log(solution%u1))*cos(2d0*pi*solution%u2)
        solution%HmY = solution%HmY + solution%HmY*problem%CV*sqrt(-2d0*log(solution%u3))*cos(2d0*pi*solution%u4)
        solution%HmZ = solution%HmZ + solution%HmZ*problem%CV*sqrt(-2d0*log(solution%u5))*cos(2d0*pi*solution%u6)
    endif
    
    call trace%end( "updateDemagfield", itimer=itimer, verbose=1 )

    end subroutine updateDemagfield
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Initializes the interaction matrices
    !> @param[inout] problem struct containing the problem information    
    !---------------------------------------------------------------------------   
    subroutine initializeInteractionMatrices( problem, solution )
    type(MicroMagProblem), intent(inout) :: problem         !> Struct containing the grid information
    type(MicroMagSolution),intent(inout) :: solution        !> Solution data structure

    integer, save :: itimer = 0
    call trace%begin( "initializeInteractionMatrices", itimer=itimer, verbose=1 )
    
    if ( problem%useDemag .eq. useDemagTrue ) then
        !Demagnetization tensor matrix
#if USE_FMM3D
        !------------- build neighbour demag tensor -------------------------------------------------------------------
        if (problem%use_fmm) then
            call BuildNeighbourDemagTensor( problem )
        end if
        
        !---------- if BuildNeighbourDemagTensor sets use_fmm to false then compute the demag tensor normally -----
        if (.not. problem%use_fmm) then
           call ComputeDemagfieldTensor( problem )
        end if
            !-----------------------------------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------------------------------------
#else
        call ComputeDemagfieldTensor( problem )
#endif
    endif
    

    !Shape correction tensor
    call ComputeShapeCorrectionTensor( problem )

    !Anisotropy matrix
    call ComputeAnisotropyTerm3D( problem )
    
    !Exhange term matrix
    call ComputeExchangeTerm3D( problem%grid, problem%A_exch, problem, solution )
    
    
    call trace%end( "initializeInteractionMatrices", itimer=itimer, verbose=1 )
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
    integer, save :: itimer = 0
    call trace%begin( "setupGrid", itimer=itimer, verbose=1 )
    
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
            end do
            
        else
            grid%x(:,:,:) = 0.
            grid%dx = grid%Lx
        endif
    
        if ( grid%ny .gt. 1 ) then
            grid%dy = grid%Ly / grid%ny
            do i=1,grid%ny
                grid%y(:,i,:) = -grid%Ly/2 + (i-1) * grid%dy + grid%dy/2
            end do
        else
            grid%y(:,:,:) = 0.
            grid%dy = grid%Ly
        endif

        if ( grid%nz .gt. 1 ) then
            grid%dz = grid%Lz / grid%nz
            do i=1,grid%nz
                grid%z(:,:,i) = -grid%Lz/2 + (i-1) * grid%dz + grid%dz/2
            end do
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
                end do
            end do
        end do
    endif
    call trace%end( "setupGrid", itimer=itimer, verbose=1 )
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
    real(DP), dimension(:,:),allocatable :: obs_size_arr          !> Array containing the size of the observation tile for the average prism demag tensor
    real(DP),dimension(:,:,:,:),allocatable :: Nout,Noutave       !> Temporary storage for the demag tensor            
    integer,dimension(4) :: indx_ele
    real :: rate
    !integer :: c1,c2,cr,cm
    character(10) :: prog_str
    integer, save :: itimer = 0
    integer,dimension(3) :: n_macro
    real(DP),dimension(3) :: shiftVec
    
    call trace%begin( "ComputeDemagfieldTensor", itimer=itimer, verbose=1 )
    
    ! First initialize the system_clock
    !call system_clock(count_rate=cr)
    !call system_clock(count_max=cm)
    !rate = REAL(cr)
    
    nx = problem%grid%nx
    ny = problem%grid%ny
    nz = problem%grid%nz
    ntot = nx * ny * nz
    
    !The number of tiles to average the receiving tile over
    nx_ave = problem%N_ave(1)
    ny_ave = problem%N_ave(2)
    nz_ave = problem%N_ave(3)

    !Macrogeometry info (added by F. Durhuus)
    n_macro = problem%macrogrid%n_macro
    shiftVec = problem%macrogrid%shiftVec

    !Demag tensor components
    allocate( problem%Kxx(ntot,ntot), problem%Kxy(ntot,ntot), problem%Kxz(ntot,ntot) )
    allocate( problem%Kyy(ntot,ntot), problem%Kyz(ntot,ntot) )
    allocate( problem%Kzz(ntot,ntot) )
        
        
    if ( problem%demagTensorLoadState .gt. DemagTensorReturnMemory ) then
        call loadDemagTensorFromDisk( problem )
    else

        !CALL SYSTEM_CLOCK(c1)
 
        !call mkl_set_num_threads(problem%nThreadsMatlab)
        !call omp_set_num_threads(problem%nThreadsMatlab)
        !call omp_set_num_threads(1)               
        if ( problem%grid%gridType .eq. gridTypeUniform ) then
            
            allocate(obs_size_arr(ntot,3))
            obs_size_arr(:,1) = problem%grid%dx
            obs_size_arr(:,2) = problem%grid%dy
            obs_size_arr(:,3) = problem%grid%dz
            
            !for each element find the tensor for all evaluation points (i.e. all elements)    
            !======== NOTE ===============
            ! this parallelization seem to give issues when compiled with matlab in debug mode
            !$OMP PARALLEL DO collapse(3) SHARED(problem, nx, ny, nz, ntot, obs_size_arr, n_macro, shiftVec) PRIVATE(ind, tile, H, Nout, pts_arr) default(none)
            do k=1,nz
                do j=1,ny                
                    do i=1,nx
                        if (problem%useAvgN .eq. useAvgNTrue) then
                            !Setup average N template tile
                            tile(1)%tileType = 8 !(for avgPrism)
                        else
                            !Setup template tile
                            tile(1)%tileType = 2 !(for prism)
                        endif
                        !Dimensions of the tile
                        tile(1)%a = problem%grid%dx
                        tile(1)%b = problem%grid%dy
                        tile(1)%c = problem%grid%dz
                        tile(1)%exploitSymmetry = 0 !0 for no and this is important
                        tile(1)%rotAngles(:) = 0.   !ensure that these are indeed zero
                        tile(1)%M(:) = 0.
                        
                        !Set the center of the tile to be the current point
                        tile(1)%offset(1) = problem%grid%x(i,j,k)
                        tile(1)%offset(2) = problem%grid%y(i,j,k)
                        tile(1)%offset(3) = problem%grid%z(i,j,k)
                        
                        allocate(Nout(1,ntot,3,3))
                        allocate(H(ntot,3))

                        allocate(pts_arr(ntot,3))
                        pts_arr(:,1) =  problem%grid%pts(:,1)
                        pts_arr(:,2) =  problem%grid%pts(:,2)
                        pts_arr(:,3) =  problem%grid%pts(:,3)

                        if (all(n_macro == 0.0)) then
                            call getFieldFromTiles( tile, H, problem%grid%pts, 1, ntot, Nout, .false. )
                        else
                            call getFieldFromTiles_PBC( tile, H, problem%grid%pts, 1, ntot, n_macro, &
                            shiftVec, Nout, .false. )
                        end if
                        
                        if (problem%useAvgN .eq. useAvgNTrue) then
                            !Use the average prism tensor
                            !call getFieldFromTiles( tile, H, problem%grid%pts, 1, ntot, Nout, .false., Obs_size=obs_size_arr)
                            call getFieldFromTiles( tile, H, pts_arr, 1, ntot, Nout, .false., Obs_size=obs_size_arr)
                        else
                            !Use the point prism tensor
                            !call getFieldFromTiles( tile, H, problem%grid%pts, 1, ntot, Nout, .false. )
                            call getFieldFromTiles( tile, H, pts_arr, 1, ntot, Nout, .false. )
                        endif
                        
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
                        deallocate(pts_arr)
                        deallocate(Nout)
                        deallocate(H)
                    enddo
                enddo
            enddo
            !$OMP END PARALLEL DO

            
            deallocate(obs_size_arr)
        elseif ( problem%grid%gridType .eq. gridTypeTetrahedron ) then
        
            if (nx_ave*ny_ave*nz_ave > 1) then
                call displayGUIMessage( 'Averaging the N_tensor not supported for this tile type' )
            endif
    
            !$OMP PARALLEL DO SHARED(problem) PRIVATE(ind, indx_ele, tile, H, Nout, pts_arr, i)
                        
            !for each tile find the tensor for all evaluation points (i.e. all tiles)
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
                
                !allocate(pts_arr(ntot,3))
                !pts_arr(:,1) =  problem%grid%pts(:,1)
                !pts_arr(:,2) =  problem%grid%pts(:,2)
                !pts_arr(:,3) =  problem%grid%pts(:,3)
                
                call getFieldFromTiles( tile, H, problem%grid%pts, 1, ntot, Nout, .false. )
                !call getFieldFromTiles( tile, H, pts_arr, 1, ntot, Nout, .false. )
                    
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
                !deallocate(pts_arr)
                deallocate(Nout)
                deallocate(H)
            end do
            
            !$OMP END PARALLEL DO
            
        elseif ( problem%grid%gridType .eq. gridTypeUnstructuredPrisms ) then
            !call displayGUIMessage( 'Constructing the Tensormap' )
            !call ConstructDemagTensorMap( problem )
            !call displayGUIMessage( 'Done constructing the Tensormap' )
            !$OMP PARALLEL DO SHARED(problem) PRIVATE(ind, tile, H, Nout,Noutave,dx,dy,dz,pts_arr)
            
            !for each tile find the tensor for all evaluation points (i.e. all tiles)
            do i=1,ntot
                if (problem%useAvgN .eq. useAvgNTrue) then
                    !Setup average N template tile
                    tile(1)%tileType = 8 !(for avgPrism)
                else
                    !Setup template tile
                    tile(1)%tileType = 2 !(for prism)
                endif
                tile(1)%exploitSymmetry = 0 !0 for no and this is important
                tile(1)%rotAngles(:) = 0.   !ensure that these are indeed zero
                tile(1)%M(:) = 0.
            
                !Dimensions of the tile
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
                !Note that abc is the full side length of the tile - it is divided with 2 in the demag tensor calculation to make it compatible to the expression in Smith_2010
                dx = problem%grid%abc(:,1)/(nx_ave+1)
                dy = problem%grid%abc(:,2)/(ny_ave+1)
                dz = problem%grid%abc(:,3)/(nz_ave+1)

                do k_a=1,nz_ave
                    do j_a=1,ny_ave                
                        do i_a=1,nx_ave
                            pts_arr(:,1) =  (problem%grid%pts(:,1)-problem%grid%abc(:,1)/2)+dx(:)*i_a
                            pts_arr(:,2) =  (problem%grid%pts(:,2)-problem%grid%abc(:,2)/2)+dy(:)*j_a
                            pts_arr(:,3) =  (problem%grid%pts(:,3)-problem%grid%abc(:,3)/2)+dz(:)*k_a
                            
                            if (problem%useAvgN .eq. useAvgNTrue) then
                                !Use the average prism tensor
                                call getFieldFromTiles( tile, H, pts_arr, 1, ntot, Nout, .false., obs_size = problem%grid%abc)
                            else
                                !Use the point prism tensor
                                !call getFieldFromTiles( tile, H, problem%grid%pts, 1, ntot, Nout, .false. )
                                call getFieldFromTiles( tile, H, pts_arr, 1, ntot, Nout, .false. )
                            endif

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
            end do
        
            !$OMP END PARALLEL DO
            
        endif
        
        !CALL SYSTEM_CLOCK(c2)        

        !---------- scale sources with the saturation magnetization -----------------
        do i=1,ntot
            problem%Kxx(i,:) = problem%Kxx(i,:) * problem%Ms(:)
            problem%Kxy(i,:) = problem%Kxy(i,:) * problem%Ms(:)
            problem%Kxz(i,:) = problem%Kxz(i,:) * problem%Ms(:)
            problem%Kyy(i,:) = problem%Kyy(i,:) * problem%Ms(:)
            problem%Kyz(i,:) = problem%Kyz(i,:) * problem%Ms(:)
            problem%Kzz(i,:) = problem%Kzz(i,:) * problem%Ms(:)
        end do
        !---------------------------------------------------------------------------


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
    
    call trace%end( "ComputeDemagfieldTensor", itimer=itimer, verbose=1 )
    
    end subroutine ComputeDemagfieldTensor

    !>-----------------------------------------
    !> @author Frederik L. Durhuus, fladu@dtu.dk, DTU, 2025
    !> Calculates and returns the shape correction tensor,
    !> which ensures that the correct sample shape is used
    !> for computing the demagnetisation field
    !> @param[inout] problem, the struct containing the problem
    !>-----------------------------------------
    subroutine ComputeShapeCorrectionTensor( problem )
    type(MicroMagProblem),intent(inout) :: problem       !> Grid data structure

    integer :: nx,ny,nz,ntot                             !> Internal counters and index variables
    real(DP),dimension(:,:),allocatable :: pts           !> Evaluation points
    real(DP),dimension(:,:,:),allocatable :: Nshape      !> Temporary storage for the demag tensor
    real(DP),allocatable :: aMacro,bMacro,cMacro         !> Macrogeometry shape parameters
    real(DP),allocatable :: aSample,bSample,cSample      !> Sample shape parameters
    integer :: ntotMacro                                 !> Number of domain copies in macrogeometry
    real(DP) :: Vcell, Vdomain                           !> Volume of single cell and simulated domain
    character*(100) :: prog_str 

    ! Get number of elements
    nx = problem%grid%nx
    ny = problem%grid%ny
    nz = problem%grid%nz
    ntot = nx * ny * nz

    ! Get macrogeometry shape parameters (approximated by rectangular prism)
    aMacro = problem%macrogrid%macroShape(1)
    bMacro = problem%macrogrid%macroShape(2)
    cMacro = problem%macrogrid%macroShape(3)
    ! Get sample shape parameters (approximated by rectangular prism)
    aSample = problem%macrogrid%sampleShape(1)
    bSample = problem%macrogrid%sampleShape(2)
    cSample = problem%macrogrid%sampleShape(3)
    ! Get array of points to evaluate tensor at
    pts = problem%grid%pts

    allocate(Nshape(ntot,3,3))
    call getShapeTensor(pts, ntot, aMacro, bMacro, cMacro, aSample, bSample, cSample, Nshape)

    ! Tensor components for shape correction tensor
    allocate( problem%Kxx_shape(ntot), problem%Kxy_shape(ntot), problem%Kxz_shape(ntot) )
    allocate( problem%Kyy_shape(ntot), problem%Kyz_shape(ntot) )
    allocate( problem%Kzz_shape(ntot) )

    ! Copy Nshape into the proper structure used by the micromag model
    problem%Kxx_shape = Nshape(:,1,1)
    problem%Kxy_shape = Nshape(:,1,2)
    problem%Kxz_shape = Nshape(:,1,3)

    ! By symmetry Kxy = Kyx
    problem%Kyy_shape = Nshape(:,2,2)
    problem%Kyz_shape = Nshape(:,2,3)

    ! By symmetry Kxz = Kzx and Kyz = Kzy
    problem%Kzz_shape = Nshape(:,3,3)

    ! Get volume fraction occupied
    ntotMacro = (1 + 2*problem%macrogrid%n_macro(1)) * (1 + 2*problem%macrogrid%n_macro(2)) * (1 + 2*problem%macrogrid%n_macro(3))
    Vdomain = aMacro * bMacro * cMacro / ntotMacro
    Vcell = problem%grid%dx * problem%grid%dy * problem%grid%dz
    allocate( problem%VfracOcc )
    problem%VfracOcc = ntot * Vcell / Vdomain
    call displayGUIMessage( 'Volume fraction occupied by magnetic material :' )
    write(prog_str,'(F20.10)') problem%VfracOcc
    call displayGUIMessage( trim(prog_str) )

    ! Clean up
    deallocate(Nshape)

    end subroutine ComputeShapeCorrectionTensor
    
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
        

    if ( problem%passExch .eq. passExchTrue) return ! Skip this if the exchange has already been passed from outside

    allocate( A0_normalized(size(problem%A0)) )  
    A0_normalized = problem%A0 / ( maxval(problem%A0) )   ! Normalized by the largest exchange factor
        
    if ( grid%gridType .eq. gridTypeUniform ) then
        call ComputeExchangeTerm3D_Uniform( grid, A, problem, solution )
    elseif (( grid%gridType .eq. gridTypeTetrahedron ) .or. (grid%gridType .eq. gridTypeUnstructuredPrisms)) then
        call computeDifferentialOperatorsFromMesh_DirectLap(solution%gridinfo, problem%exch_interpn, problem%exch_weight, problem%exch_method, A0_normalized, problem%A_exch)
    endif

    
    end subroutine ComputeExchangeTerm3D
    
        
    !>-----------------------------------------
    !> @author Rasmus Bjørk, rabj@dtu.dk, DTU, 2026
    !> @brief
    !> Calculates the exhange term matrix on a uniform grid
    !> which means it produces the differential operator d^2/dx^2 + d^2/dy^2 + d^2/dz^2 and returns this in the sparse matrix A
    !> Includes position dependent exchange stiffness A0 through the modified stencil described in 
    !> Heistracher et al. "Proposal for a micromagnetic standard problem: Domain wall pinning at phase boundaries", DOI: 10.1016/j.jmmm.2021.168875
    !---------------------------------------------------------------------------   
    subroutine ComputeExchangeTerm3D_Uniform( grid, A, problem, solution )
    type(MicroMagGrid),intent(in) :: grid              !> Struct containing the grid information    
    type(sparse_matrix_t),intent(inout) :: A           !> The returned matrix from the sparse matrix creator
    type(MicroMagProblem),intent(inout) :: problem     !> Problem data structure
    type(MicroMagSolution),intent(inout) :: solution   !> Solution data structure
        
    integer :: stat                                    !> Status value for the various sparse matrix operations        
    type(MagTenseSparse_d) :: d2dx2, d2dy2, d2dz2      !> Sparse matrices for the double derivatives with respect to x, y and z, respectively.
    type(sparse_matrix_t) :: tmp                       !> Temporary sparse matrices used for internal calculations
    integer :: ind,ntot,nNonZero,colInd,rowInd         !> Internal counter for indexing, the total no. of elements in the current sparse matrix being manipulated
    integer :: i,j,k,nx,ny,nz                          !> For-loop counters
    type(matrix_descr) :: descr                        !> Describes a sparse matrix operation
    real(DP) :: const
    character*(100) :: prog_str 
    integer, save :: itimer = 0
    logical :: PBCx, PBCy, PBCz                        !> Whether to have periodic boundary conditions along x, y and z respectively
    
    call trace%begin( "ComputeExchangeTerm3D_Uniform", itimer=itimer, verbose=1 )
    
    !Find the three sparse matrices for the individual directions. Then add them to get the total matrix
    !It is assumed that the magnetization vector to operate on is in fact a single column of Mx, My or Mz respectively.

    !The sparse matrices are encoded using an older version of Intels Compressed Sparse Row (CSR) format
    !cols(i) and values(i) give the column index and value of the i'th non-zero element
    !rows_start(i) and rows_end(i)-1 give the first and last index corresponding to the i'th row
    !If rows_start(i) = rows_end(i), the row is empty
    !For instance the matrix [[1, 0, 2], [0, -1, 0], [0, 0, 0], [-2, 5, 4]] has cols = (1, 3, 2, 1, 2, 3),
    !values = (1, 2, -1, -2, 5, 4), rows_start = (1, 3, 4, 4) and rows_end = (3, 4, 4, 7)
    !Note how size(rows_start) and size(rows_end) equal the number of rows, while size(cols) and size(values) equal the number of non-zero elements
    
    nx = grid%nx
    ny = grid%ny
    nz = grid%nz
    ntot = nx*ny*nz     ! Total number of micromagnetic tiles
    PBCx = problem%macrogrid%exchPBC(1)
    PBCy = problem%macrogrid%exchPBC(2)
    PBCz = problem%macrogrid%exchPBC(3)
    
    !----------------------------------d^2dx^2 begins -----------------------------!
    if ( nx .gt. 1 ) then
        ! Make the d^2/dx^2 matrix
        if ( PBCx ) then
            nNonZero = 3*ntot             ! For each tile, there are 3 non-zero elements corresponding to itself and the two neighbours along x
        else
            nNonZero = 3*ntot - 2*ny*nz   ! With a vacuum boundary, the 2*ny*nz tiles at the two end surfaces only have one neighbour each
        end if
        allocate(d2dx2%values(nNonZero),d2dx2%cols(nNonZero),d2dx2%rows_start(ntot),d2dx2%rows_end(ntot))
    
        ind = 1
        rowInd = 1
        colInd = 1
    
        do k=1,nz
            do j=1,ny
                ! The left boundary
                if ( PBCx ) then
                    ! Coupling to the left (loops around and connects with rightmost tile)
                    d2dx2%values(ind) = 2*problem%A0_map(nx,j,k)/(problem%A0_map(nx,j,k)+problem%A0_map(1,j,k))
                    d2dx2%cols(ind) = colInd + nx - 1
                    d2dx2%rows_start(rowInd) = ind
                    ind = ind + 1

                    ! Current tile
                    d2dx2%values(ind) = 2*(-( problem%A0_map(nx,j,k)/(problem%A0_map(nx,j,k)+problem%A0_map(1,j,k)) &
                            + problem%A0_map(2,j,k)/(problem%A0_map(2,j,k)+problem%A0_map(1,j,k)) ))
                    d2dx2%cols(ind) = colInd
                    ind = ind + 1

                    ! Coupling to the right
                    d2dx2%values(ind) = 2*problem%A0_map(2,j,k)/(problem%A0_map(2,j,k)+problem%A0_map(1,j,k))
                    d2dx2%cols(ind) = colInd + 1
                    d2dx2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                    ind = ind + 1
                else
                    !The left boundary
                    const = 2 *  problem%A0_map(2,j,k) / (problem%A0_map(2,j,k) + problem%A0_map(1,j,k))

                    d2dx2%values(ind) = -const
                    d2dx2%cols(ind) = colInd
                    d2dx2%rows_start(rowInd) = ind            
                    ind = ind + 1
                
                    d2dx2%values(ind) = const
                    d2dx2%cols(ind) = colInd + 1
                    d2dx2%rows_end(rowInd) = ind+1
                    rowInd = rowInd + 1
                    ind = ind + 1
                end if
                        
            
                !Go through one row at a time
                do i=2,nx-1
                            
                    !Left-most point
                    d2dx2%values( ind ) = 2 * problem%A0_map(i-1,j,k)/(problem%A0_map(i-1,j,k)+problem%A0_map(i,j,k))
                    
                    d2dx2%cols(ind) = colInd
                    !update where the row starts
                    d2dx2%rows_start(rowInd) = ind                           
                    ind = ind + 1
                
                    !Center point
                    d2dx2%values( ind ) = 2 * ( -(problem%A0_map(i-1,j,k)/(problem%A0_map(i-1,j,k)+problem%A0_map(i,j,k))+problem%A0_map(i+1,j,k)/(problem%A0_map(i+1,j,k)+problem%A0_map(i,j,k))))
                    d2dx2%cols(ind) = colInd + 1
                    ind = ind + 1
                
                    !Right-most point
                    d2dx2%values( ind ) = 2 * problem%A0_map(i+1,j,k)/(problem%A0_map(i+1,j,k)+problem%A0_map(i,j,k))
                    d2dx2%cols(ind) = colInd + 2
                    d2dx2%rows_end(rowInd) = ind+1
                    rowInd = rowInd + 1
                    ind = ind + 1

                    colInd = colInd + 1
                enddo            

                ! The right boundary
                if ( PBCx ) then

                    ! Coupling to the left
                    d2dx2%values(ind) = 2*problem%A0_map(nx-1,j,k)/(problem%A0_map(nx-1,j,k)+problem%A0_map(nx,j,k))
                    d2dx2%cols(ind) = colInd
                    ! Update where the row starts
                    d2dx2%rows_start(rowInd) = ind
                    ind = ind + 1

                    ! Current tile
                    d2dx2%values(ind) = 2*(-( problem%A0_map(nx-1,j,k)/(problem%A0_map(nx-1,j,k)+problem%A0_map(nx,j,k)) &
                            + problem%A0_map(1,j,k)/(problem%A0_map(1,j,k)+problem%A0_map(nx,j,k)) ))
                    d2dx2%cols(ind) = colInd + 1
                    ind = ind + 1

                    ! Coupling to the right (loops around and connects with leftmost tile)
                    d2dx2%values(ind) = 2*problem%A0_map(1,j,k)/(problem%A0_map(1,j,k)+problem%A0_map(nx,j,k))
                    d2dx2%cols(ind) = colInd - nx + 2
                    d2dx2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                    ind = ind + 1

                    colInd = colInd + 2
                else
                    const = 2 * problem%A0_map(nx-1,j,k) / (problem%A0_map(nx-1,j,k) + problem%A0_map(nx,j,k))
                    !The right boundary
                    d2dx2%values(ind) = const 
                    d2dx2%cols(ind) = colInd
                    !update where the row starts
                    d2dx2%rows_start(rowInd) = ind            
                    ind = ind + 1
            
                    d2dx2%values(ind) = -const
                    d2dx2%cols(ind) = colInd+1
                    d2dx2%rows_end(rowInd) = ind+1
                    rowInd = rowInd + 1
                    ind = ind + 1     
            
                    colInd = colInd + 2
                endif
            enddo
        enddo
    
        !Multiply by the discretization
        d2dx2%values = d2dx2%values * 1.0/grid%dx**2
        
        ! Create the sparse matrix for the d^2dx^2
        stat = mkl_sparse_d_create_csr ( d2dx2%A, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, d2dx2%rows_start, d2dx2%rows_end, d2dx2%cols, d2dx2%values)
    
    endif
    
    !----------------------------------d^2dx^2 ends ---------------------------- -!
    
    
    !----------------------------------d^2dy^2 begins ----------------------------!
    if ( ny .gt. 1 ) then
        ! Make the d^2/dy^2 matrix
        if ( PBCy ) then
            nNonZero = 3*ntot             ! For each tile, there are 3 non-zero elements corresponding to itself and the two neighbours along y
        else
            nNonZero = 3*ntot - 2*nx*nz   ! With a vacuum boundary, the 2*nx*nz tiles at the two end surfaces only have one neighbour each
        end if
        allocate(d2dy2%values(nNonZero),d2dy2%cols(nNonZero),d2dy2%rows_start(ntot),d2dy2%rows_end(ntot))
    
        ind = 1
        rowInd = 1
        colInd = 1
        do k=1,nz
            !The bottom boundary
             ! The backwards boundary (Imagine standing on the origin, eyes pointed along the positive y-axis)
            if ( PBCy ) then
                do i=1,nx
                    ! Backwards coupling (loops around and couples to forwardmost tile)
                    d2dy2%values(ind) = 2*problem%A0_map(i,ny,k)/(problem%A0_map(i,ny,k)+problem%A0_map(i,1,k))
                    d2dy2%cols(ind) = colInd + nx*(ny-1)
                    d2dy2%rows_start(rowInd) = ind
                    ind = ind + 1  ! Increment to next element

                    ! Current tile
                    d2dy2%values(ind) = 2*(-(problem%A0_map(i,ny,k)/(problem%A0_map(i,ny,k)+problem%A0_map(i,1,k)) &
                            + problem%A0_map(i,2,k)/(problem%A0_map(i,2,k)+problem%A0_map(i,1,k))))
                    d2dy2%cols(ind) = colInd
                    ind = ind + 1  ! Increment to next element

                    ! Forwards coupling
                    d2dy2%values(ind) = 2*problem%A0_map(i,2,k)/(problem%A0_map(i,2,k)+problem%A0_map(i,1,k))
                    d2dy2%cols(ind) = colInd + nx
                    d2dy2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                    ind = ind + 1  ! Increment to next element
                    colInd = colInd + 1
                end do
            else
            
                do i=1,nx
                    const = 2 * problem%A0_map(i,2,k) / (problem%A0_map(i,2,k) + problem%A0_map(i,1,k))

                    d2dy2%values(ind) = -const
                    d2dy2%cols(ind) = colInd
                    d2dy2%rows_start(rowInd) = ind
            
                    !increment to next element
                    ind = ind + 1
            
                    d2dy2%values(ind) = const
                    d2dy2%cols(ind) = colInd + nx
                    d2dy2%rows_end(rowInd) = ind+1
                    rowInd = rowInd + 1
            
                    !increment to next element
                    ind = ind + 1
            
                    colInd = colInd + 1
                enddo
            endif
            
            ! Everything in between
            do j=2,ny-1
                do i=1,nx
                
                    !lower value
                    d2dy2%values(ind) = 2 * problem%A0_map(i,j-1,k)/(problem%A0_map(i,j-1,k)+problem%A0_map(i,j,k))
                    d2dy2%cols(ind) = colInd-nx
                    d2dy2%rows_start(rowInd) = ind
                    ind = ind + 1  ! Increment to next element
                
                    !central value
                    d2dy2%values(ind) = 2*(-(problem%A0_map(i,j-1,k)/(problem%A0_map(i,j-1,k)+problem%A0_map(i,j,k))+problem%A0_map(i,j+1,k)/(problem%A0_map(i,j+1,k)+problem%A0_map(i,j,k))))
                    d2dy2%cols(ind) = colInd
                    ind = ind + 1  ! Increment to next element
            
            
                    !upper value
                    d2dy2%values(ind) = 2 * problem%A0_map(i,j+1,k)/(problem%A0_map(i,j+1,k)+problem%A0_map(i,j,k))
                    d2dy2%cols(ind) = colInd + nx
                    d2dy2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                    ind = ind + 1  ! Increment to next element
                    colInd = colInd + 1
                end do
            end do
        
            !The forwards boundary
            if ( PBCy ) then
                do i=1,nx
                    ! Backwards coupling
                    d2dy2%values(ind) = problem%A0_map(i,ny-1,k)/(problem%A0_map(i,ny-1,k) + problem%A0_map(i,ny,k))
                    d2dy2%cols(ind) = colInd - nx
                    d2dy2%rows_start(rowInd) = ind
                    ind = ind + 1  ! Increment to next element

                    ! Current tile
                    d2dy2%values(ind) = -(problem%A0_map(i,ny-1,k)/(problem%A0_map(i,ny-1,k) + problem%A0_map(i,ny,k)) &
                            + problem%A0_map(i,1,k)/(problem%A0_map(i,1,k) + problem%A0_map(i,ny,k)))
                    d2dy2%cols(ind) = colInd
                    ind = ind + 1  ! Increment to next element

                    ! Forwards coupling (loops around and couples to backmost tile)
                    d2dy2%values(ind) = problem%A0_map(i,1,k)/(problem%A0_map(i,1,k) + problem%A0_map(i,ny,k))
                    d2dy2%cols(ind) = colInd - nx*(ny-1)
                    d2dy2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                    ind = ind + 1  ! Increment to next element
                    colInd = colInd + 1
                enddo
            else
                do i=1,nx
                   const = 2 * problem%A0_map(i,ny-1,k) / (problem%A0_map(i,ny-1,k) + problem%A0_map(i,ny,k))
                    !lower element
                    d2dy2%values(ind) = const
                    d2dy2%cols(ind) = colInd - nx
                    d2dy2%rows_start(rowInd) = ind
            
                    !increment to next element
                    ind = ind + 1            
            
                    !central element
                    d2dy2%values(ind) = -const
                    d2dy2%cols(ind) = colInd
                    d2dy2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
            
                    ind = ind + 1
                    colInd = colInd + 1
                enddo
            endif
        enddo
                    
        !Multiply by the discretization
        d2dy2%values = d2dy2%values * 1./grid%dy**2
        
        ! Create the sparse matrix for the d^2dy^2
        stat = mkl_sparse_d_create_csr ( d2dy2%A, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, d2dy2%rows_start, d2dy2%rows_end, d2dy2%cols, d2dy2%values)
    endif
    
    !----------------------------------d^2dy^2 ends ----------------------------!
    
    
    !----------------------------------d^2dz^2 begins ----------------------------!
    if ( nz .gt. 1 ) then
        ! Make the d^2/dz^2 matrix
        if ( PBCz ) then
            nNonZero = 3*ntot             ! For each tile, there are 3 non-zero elements corresponding to itself and the two neighbours along z
        else
            nNonZero = 3*ntot - 2*nx*ny   ! With a vacuum boundary, the 2*nx*ny tiles at the two end surfaces only have one neighbour each
        end if
        allocate(d2dz2%values(nNonZero),d2dz2%cols(nNonZero),d2dz2%rows_start(ntot),d2dz2%rows_end(ntot))
    
        ind = 1
        rowInd = 1
        colInd = 1
        ! The bottom boundary
        if ( PBCz ) then
            do j=1,ny
                do i=1,nx
                    ! Downwards coupling (loops around and couples to topmost tile)
                    d2dz2%values(ind) = 2*problem%A0_map(i,j,nz)/(problem%A0_map(i,j,nz)+problem%A0_map(i,j,1))
                    d2dz2%cols(ind) = colInd + ny*nx*(nz-1)
                    d2dz2%rows_start(rowInd) = ind
                    ind = ind + 1   ! Increment to next element

                    ! Current tile
                    d2dz2%values(ind) = 2*(-(problem%A0_map(i,j,nz)/(problem%A0_map(i,j,nz) + problem%A0_map(i,j,1)) &
                            + problem%A0_map(i,j,2)/(problem%A0_map(i,j,2) + problem%A0_map(i,j,1))))
                    d2dz2%cols(ind) = colInd
                    d2dz2%rows_start(rowInd) = ind
                    ind = ind + 1   ! Increment to next element

                    ! Upwards coupling
                    d2dz2%values(ind) = 2*problem%A0_map(i,j,2)/(problem%A0_map(i,j,2)+problem%A0_map(i,j,1))
                    d2dz2%cols(ind) = nx*ny + colInd
                    d2dz2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                    ind = ind + 1   ! Increment to next element
                    colInd = colInd + 1
                end do
            end do
        else
            do j=1,ny
                do i=1,nx
                    const = 2 * problem%A0_map(i,j,2) / (problem%A0_map(i,j,2) + problem%A0_map(i,j,1))
                    !central value
                    d2dz2%values(ind) = -const
                    d2dz2%cols(ind) = colInd
                    d2dz2%rows_start(rowInd) = ind
            
                    !increment position
                    ind = ind + 1
            
                    !right-most value
                    d2dz2%values(ind) = const
                    d2dz2%cols(ind) = nx * ny + colInd
                    d2dz2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
            
                    !increment position
                    ind = ind + 1
                    colInd = colInd + 1
                end do
            end do
        endif

        ! Everything in between
        do k=2,nz-1            
            do j=1,ny
                do i=1,nx
                    !left-most value
                    d2dz2%values(ind) = 2 * problem%A0_map(i,j,k-1)/(problem%A0_map(i,j,k-1)+problem%A0_map(i,j,k))
                    d2dz2%cols(ind) = colInd - nx * ny
                    d2dz2%rows_start(rowInd) = ind
                    ind = ind + 1   ! Increment to next element
                
                    !central value
                    d2dz2%values(ind) = 2 * (-(problem%A0_map(i,j,k-1)/(problem%A0_map(i,j,k-1)+problem%A0_map(i,j,k))+problem%A0_map(i,j,k+1)/(problem%A0_map(i,j,k+1)+problem%A0_map(i,j,k))))
                    d2dz2%cols(ind) = colInd
                    ind = ind + 1   ! Increment to next element
                
                    !right-most value
                    d2dz2%values(ind) = 2 * problem%A0_map(i,j,k+1)/(problem%A0_map(i,j,k+1)+problem%A0_map(i,j,k))
                    d2dz2%cols(ind) = colInd + nx * ny
                    d2dz2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                    ind = ind + 1   ! Increment to next element
                    colInd = colInd + 1
                end do
            end do
        end do

        ! The top boundary
        if ( PBCz ) then
            do j=1,ny
                do i=1,nx
                    ! Downwards coupling
                    d2dz2%values(ind) = problem%A0_map(i,j,nz-1)/(problem%A0_map(i,j,nz-1) + problem%A0_map(i,j,nz))
                    d2dz2%cols(ind) = colInd - nx*ny
                    d2dz2%rows_start(rowInd) = ind
                    ind = ind + 1   ! Increment to next element

                    ! Current tile
                    d2dz2%values(ind) = -(problem%A0_map(i,j,nz-1)/(problem%A0_map(i,j,nz-1) + problem%A0_map(i,j,nz)) &
                            + problem%A0_map(i,j,1)/(problem%A0_map(i,j,1) + problem%A0_map(i,j,nz)))
                    d2dz2%cols(ind) = colInd
                    ind = ind + 1   ! Increment to next element

                    ! Upwards coupling (loops around and couples to bottommost tile)
                    d2dz2%values(ind) = problem%A0_map(i,j,1)/(problem%A0_map(i,j,1) + problem%A0_map(i,j,nz))
                    d2dz2%cols(ind) = colInd - nx*ny*(nz - 1)
                    d2dz2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                    ind = ind + 1   ! Increment to next element
                    colInd = colInd + 1
                end do
            end do
        else
            do j=1,ny
                do i=1,nx
                    const = 2 * problem%A0_map(i,j,nz-1) / (problem%A0_map(i,j,nz-1) + problem%A0_map(i,j,nz))

                    !left-most value
                    d2dz2%values(ind) = const
                    d2dz2%cols(ind) = colInd - nx * ny
                    d2dz2%rows_start(rowInd) = ind
                
                    !increment position
                    ind = ind + 1
                
                    !central value
                    d2dz2%values(ind) = -const
                    d2dz2%cols(ind) = colInd
                    d2dz2%rows_end(rowInd) = ind + 1
                    rowInd = rowInd + 1
                
                    !increment position
                    ind = ind + 1
                    colInd = colInd + 1
                end do
            end do
        endif
        
        !Multiply by the discretization
        d2dz2%values = d2dz2%values * 1./grid%dz**2
        
        !Create the sparse matrix for the d^2dz^2
        stat = mkl_sparse_d_create_csr ( d2dz2%A, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, d2dz2%rows_start, d2dz2%rows_end, d2dz2%cols, d2dz2%values)
    endif
    
    !----------------------------------d^2dz^2 ends ----------------------------!
    
            
    !Finally, add up the three diagonals and store in the output sparse matrix, A
    !Store the results temporarily in tmp4
    
    !call writeSparseMatrixToDisk( d2dz2%A, nx*ny*nz, 'd2dz2.dat' )
    !call writeSparseMatrixToDisk( d2dy2%A, nx*ny*nz, 'd2dy2.dat' )
    
    const = 1.
    
    if ( nx .gt. 1) then
        if ( ny .gt. 1) then
            stat = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, d2dx2%A, const, d2dy2%A, tmp)  
            deallocate(d2dy2%values,d2dy2%cols,d2dy2%rows_start,d2dy2%rows_end)
            stat = mkl_sparse_destroy (d2dy2%A)
            
            if ( nz .gt. 1) then
                !nx+ny+nz
                call displayGUIMessage( 'Exchange in x,y,z' )
                stat = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, d2dz2%A, const, tmp, A)       
                deallocate(d2dz2%values,d2dz2%cols,d2dz2%rows_start,d2dz2%rows_end)
                stat = mkl_sparse_destroy (d2dz2%A)
            else
                !nx+ny
                call displayGUIMessage( 'Exchange in x,y' )
                descr%type = SPARSE_MATRIX_TYPE_GENERAL
                descr%mode = SPARSE_FILL_MODE_FULL
                descr%diag = SPARSE_DIAG_NON_UNIT
                stat = mkl_sparse_copy ( tmp, descr, A )
                stat = mkl_sparse_destroy (tmp)
            endif
        else
            if ( nz .gt. 1) then
                !nx+nz
                call displayGUIMessage( 'Exchange in x,z' )
                stat = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, d2dx2%A, const, d2dz2%A, A)  
                deallocate(d2dz2%values,d2dz2%cols,d2dz2%rows_start,d2dz2%rows_end)
                stat = mkl_sparse_destroy (d2dz2%A)
            else
                !nx
                call displayGUIMessage( 'Exchange in x' )
                descr%type = SPARSE_MATRIX_TYPE_GENERAL
                descr%mode = SPARSE_FILL_MODE_FULL
                descr%diag = SPARSE_DIAG_NON_UNIT
                stat = mkl_sparse_copy ( d2dx2%A, descr, A )
            endif
        endif
        deallocate(d2dx2%values,d2dx2%cols,d2dx2%rows_start,d2dx2%rows_end)
        stat = mkl_sparse_destroy (d2dx2%A)
    else
        if ( ny .gt. 1) then
            if ( nz .gt. 1) then
                !ny+nz
                call displayGUIMessage( 'Exchange in y,z' )
                stat = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, d2dy2%A, const, d2dz2%A, A) 
                deallocate(d2dz2%values,d2dz2%cols,d2dz2%rows_start,d2dz2%rows_end)
                stat = mkl_sparse_destroy (d2dz2%A)
            else
                !ny
                call displayGUIMessage( 'Exchange in y' )
                descr%type = SPARSE_MATRIX_TYPE_GENERAL
                descr%mode = SPARSE_FILL_MODE_FULL
                descr%diag = SPARSE_DIAG_NON_UNIT
                stat = mkl_sparse_copy ( d2dy2%A, descr, A )
            endif
            deallocate(d2dy2%values,d2dy2%cols,d2dy2%rows_start,d2dy2%rows_end)
            stat = mkl_sparse_destroy (d2dy2%A)
        else
            if ( nz .gt. 1) then
                !nz
                call displayGUIMessage( 'Exchange in z' )
                descr%type = SPARSE_MATRIX_TYPE_GENERAL
                descr%mode = SPARSE_FILL_MODE_FULL
                descr%diag = SPARSE_DIAG_NON_UNIT
                stat = mkl_sparse_copy ( d2dz2%A, descr, A )
            else
                !no exchange terms
            endif
            deallocate(d2dz2%values,d2dz2%cols,d2dz2%rows_start,d2dz2%rows_end)
            stat = mkl_sparse_destroy (d2dz2%A)
        endif
    endif
    
    !call writeSparseMatrixToDisk( tmp, nx*ny*nz, 'A_exch.dat' )
    
    !call writeSparseMatrixToDisk( A, nx*ny*nz, 'A_exch.dat' )
    
    call create_COO_values_from_CSR(A,solution%gridinfo)
    

    call trace%end( "ComputeExchangeTerm3D_Uniform", itimer=itimer, verbose=1 )
    end subroutine ComputeExchangeTerm3D_Uniform
       
    !>-----------------------------------------
    !> @author Rasmus Bjørk, rabj@dtu.dk, DTU, 2020
    !> @brief
    !> Calculates the anisotropy term matrix on any grid  
    !> @param[inout] problem the data structure containing the problem
    !---------------------------------------------------------------------------   
    subroutine ComputeAnisotropyTerm3D( problem )
    type(MicroMagProblem),intent(inout) :: problem             !> Struct containing the problem
    
    integer :: nx,ny,nz,ntot, i
    integer, save :: itimer = 0
    
    ntot = size(problem%Ms)

    call trace%begin( "ComputeAnisotropyTerm3D", itimer=itimer, verbose=1 )
    !--- We use the general formulation introduced in updateAnisotropy
    !--- If the user has specified a value for the uniaxial anisotropy or the cubic anisotropy, we use transform
    !--- those to the general matrix formulation
    allocate(problem%Kfact_arr(size(problem%Ms),6,3))
    problem%Kfact_arr(:,:,:) = 0.0
    
    if (any(problem%K0 .ne. 0)) then !Uniaxial anisotropy
        call displayGUIMessage( 'Assuming uniaxial anisotropy' )
        do i = 1,size(problem%Ms)
            problem%Kfact_arr(i,1,3) = problem%K0(i) / ( mu0 * problem%Ms(i) )
        end do

        !Set the crystal axis as the user only specified u_ea
        problem%CrystalAxis(:,3,1) =  problem%u_ea(:,1)
        problem%CrystalAxis(:,3,2) =  problem%u_ea(:,2)
        problem%CrystalAxis(:,3,3) =  problem%u_ea(:,3)

        allocate(problem%Axx(ntot),problem%Axy(ntot),problem%Axz(ntot), &
                 problem%Ayy(ntot),problem%Ayz(ntot),problem%Azz(ntot))
        problem%Axx = problem%u_ea(:,1) * problem%u_ea(:,1)
        problem%Axy = problem%u_ea(:,1) * problem%u_ea(:,2)
        problem%Axz = problem%u_ea(:,1) * problem%u_ea(:,3)
        problem%Ayy = problem%u_ea(:,2) * problem%u_ea(:,2)
        problem%Ayz = problem%u_ea(:,2) * problem%u_ea(:,3)
        problem%Azz = problem%u_ea(:,3) * problem%u_ea(:,3)

        
    elseif (any(problem%K1 .ne. 0)) then !Cubic anisotropy
        call displayGUIMessage( 'Assuming cubic anisotropy' )
        do i = 1,size(problem%Ms)
            problem%Kfact_arr(i,3,1) = problem%K1(i) / ( mu0 * problem%Ms(i) )
            problem%Kfact_arr(i,3,2) = problem%K1(i) / ( mu0 * problem%Ms(i) )
            problem%Kfact_arr(i,3,3) = problem%K1(i) / ( mu0 * problem%Ms(i) )
            problem%Kfact_arr(i,6,1) = problem%K2(i) / ( mu0 * problem%Ms(i) )
        end do
    else !General anisotropy
        do i = 1,size(problem%Ms)
            problem%Kfact_arr(i,:,:) = problem%K0_arr(i,:,:) / ( mu0 * problem%Ms(i) )
        end do
    endif

    call trace%end( "ComputeAnisotropyTerm3D", itimer=itimer, verbose=1 )
    
    end subroutine ComputeAnisotropyTerm3D

    


    !================== Helper functions for FMM correction =========================================
!>======================================================================
!> Add demagnitization field from nearfield 
!> uses K_nbrcorr, which depending on the setup contains either
!>      1. The correction tensors (Nnbr - Kdip*V) stored in sparse format
!>      2. The demagnetization tensor Nnbr stored in sparse format
!> adds computed field to the existing solution%H
!>======================================================================
subroutine add_near_field(problem, solution)
  implicit none
  type(MicroMagProblem),  intent(in)    :: problem
  type(MicroMagSolution), intent(inout) :: solution

  !-----------------------------------------------------------------------------
  real(SP), contiguous, pointer :: hx_tmp(:), hy_tmp(:), hz_tmp(:), temp(:)
  integer :: ntot
  real(SP) :: pref
  integer :: stat
  real(SP) :: alpha, beta
  integer, save :: itimer = 0
  call trace%begin( "add_near_field", itimer=itimer, verbose=2 )


  ntot = size(problem%grid%pts, dim=1)


#if USE_CUDA
    allocate( hx_tmp(ntot), hy_tmp(ntot), hz_tmp(ntot) )
    hx_tmp = 0.0_SP
    hy_tmp = 0.0_SP
    hz_tmp = 0.0_SP

    pref = sngl(1.0)
    call cudaMatrVecMult_sparse( solution%Mx_s , solution%My_s , solution%Mz_s , hx_tmp, hy_tmp, hz_tmp, pref )

    !$omp critical (solution_update)
    solution%HmX = solution%HmX + hx_tmp 
    solution%HmY = solution%HmY + hy_tmp 
    solution%HmZ = solution%HmZ + hz_tmp 
    !$omp end critical (solution_update)

    deallocate(hx_tmp, hy_tmp, hz_tmp)

#else
    allocate(temp(ntot))

    alpha = 1.0_SP
    ! ---------------- Hx correction = xx*Mx + xy*My + xz*Mz ----------------
    beta = 0.0_SP
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(1)%A, problem%K_fmm_descr_s, solution%Mx_s , beta, temp)
    beta = 1.0_SP
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(2)%A, problem%K_fmm_descr_s, solution%My_s, beta, temp)
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(3)%A, problem%K_fmm_descr_s, solution%Mz_s, beta, temp)

    !$omp critical (solution_update)
    solution%HmX = solution%HmX - temp
    !$omp end critical (solution_update)
    ! ---------------- Hy correction = xy*Mx + yy*My + yz*Mz ----------------
    beta = 0.0_SP
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(2)%A, problem%K_fmm_descr_s, solution%Mx_s , beta, temp)
    beta = 1.0_SP
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(4)%A, problem%K_fmm_descr_s, solution%My_s, beta, temp)
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(5)%A, problem%K_fmm_descr_s, solution%Mz_s, beta, temp)

    !$omp critical (solution_update)
    solution%HmY = solution%HmY - temp
    !$omp end critical (solution_update)
    ! ---------------- Hz correction = xz*Mx + yz*My + zz*Mz ----------------
    beta = 0.0_SP
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(3)%A, problem%K_fmm_descr_s, solution%Mx_s , beta, temp)
    beta = 1.0_SP
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(5)%A, problem%K_fmm_descr_s, solution%My_s, beta, temp)
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(6)%A, problem%K_fmm_descr_s, solution%Mz_s, beta, temp)

    !$omp critical (solution_update)
    solution%HmZ = solution%HmZ - temp
    !$omp end critical (solution_update)
    deallocate(temp)
#endif

    call trace%end( "add_near_field", itimer=itimer, verbose=2 )
end subroutine add_near_field

        !------------------------------------------------------------------------------------------------
    !============================================================================================================




    
end module LandauLifshitzSolution
    