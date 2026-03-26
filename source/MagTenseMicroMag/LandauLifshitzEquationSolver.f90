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
    type(MicroMagProblem),intent(inout) :: prob     !> The problem data structure
    type(MicroMagSolution),intent(inout) :: sol     !> The solution data structure    
    type(MicroMagGridInfo) :: gridinfo              !> The grid information structure
    integer :: ntot,i,j,k,ind,nt,nt_Hext,stat       !> total no. of tiles
    procedure(dydt_fct), pointer :: fct             !> Input function pointer for the function to be integrated
    procedure(callback_fct),pointer :: cb_fct       !> Callback function for displaying progress
    real(DP),dimension(:,:,:),allocatable :: M_out  !> Internal buffer for the solution (M) on the form (3*ntot,nt)
    real(DP) :: A0_max
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
        enddo
    endif   
    
    if ( gb_problem%grid%gridType .eq. gridTypeUniform ) then
        !Setup the grid
        call setupGrid( gb_problem%grid )
        
        !Set up the mesh for an uniform grid
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
                    enddo
                enddo
            enddo
        endif
    endif
      
    call displayGUIMessage( 'Initializing matrices' )
    !Calculate the interaction matrices
    call initializeInteractionMatrices( gb_problem, gb_solution )
    
    !write(prog_str,'(A37, F8.4, A5)') 'Demagnetization tensor memory usage: ', 6*storage_size(gb_problem%Kxx)*ntot/(8*2**30), ' gigabytes'
    !call displayGUIMessage( trim(prog_str) )    
    
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

            !TODO - make this as option depending on flag/use input
            !call cudaDumpDemagDense("CUDA_dense_K_matrices.bin")

        endif
#endif
#else
        call displayGUIMessage( 'MagTense not compiled with CUDA - exiting!' )
        stop
#endif
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
    if ( gb_problem%dummy_run .eq. 1) then
       call displayGUIMessage( 'Performing dummy run - skipping time evolution ' )
    else 

        CALL SYSTEM_CLOCK(c1)

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



            ! open (11, file="sparse_CUDA_H.bin",  &
            !         status='unknown', form='unformatted', &
            !         access='direct', recl=1*ntot*ntot)
            ! write(11,rec=1) gb_solution%HmX 
            ! write(11,rec=2) gb_solution%HmY
            ! write(11,rec=3) gb_solution%HmZ
            ! close(11)
            ! error stop " test stop after cuda sparse"
              
      enddo

      
        CALL SYSTEM_CLOCK(c2)
        call displayGUIMessage( 'Time simulation time integration:' )
        write (prog_str,'(f10.3)') (c2 - c1)/rate
        call displayGUIMessage( prog_str )


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
        integer, save :: itimer = 0
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
        call updateExchangeTerms( gb_problem, gb_solution )
        !-------------------------------------------------------------
        !-------------- add external field ---------------------------
        call updateExternalField( gb_problem, gb_solution, t )
        !-------------------------------------------------------------
        !-------------- add anisotropy term --------------------------
        call updateAnisotropy(  gb_problem, gb_solution )
        !-------------------------------------------------------------
        !-------------- add demagnetization field --------------------
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
        !--------------- combine to get effective field, Heff -------------
        HeffX = gb_solution%HhX + gb_solution%HjX + gb_solution%HmX + gb_solution%HkX
        HeffY = gb_solution%HhY + gb_solution%HjY + gb_solution%HmY + gb_solution%HkY
        HeffZ = gb_solution%HhZ + gb_solution%HjZ + gb_solution%HmZ + gb_solution%HkZ
        !-------------------------------------------------------------
        !--------------- calculate precession term: m x Heff -------------
        crossX = -1.0_DP * ( gb_solution%My * HeffZ - gb_solution%Mz * HeffY )
        crossY = -1.0_DP * ( gb_solution%Mz * HeffX - gb_solution%Mx * HeffZ )
        crossZ = -1.0_DP * ( gb_solution%Mx * HeffY - gb_solution%My * HeffX )
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
        call updateExternalField( gb_problem, gb_solution, gb_problem%t(j) )
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
    enddo
    
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
    
    alpha = -2.! * solution%Jfact
    !alpha = -2. * problem%A0 / ( mu0 * 8.0e5 )
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
        solution%Hkx = -2. * problem%Kfact * ( problem%Axx * solution%Mx + problem%Axy * solution%My + problem%Axz * solution%Mz )
        solution%Hky = -2. * problem%Kfact * ( problem%Axy * solution%Mx + problem%Ayy * solution%My + problem%Ayz * solution%Mz )
        solution%Hkz = -2. * problem%Kfact * ( problem%Axz * solution%Mx + problem%Ayz * solution%My + problem%Azz * solution%Mz )
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

        solution%Hkx(:) = problem%CrystalAxis(:,1,1)*Hkx_rot(:) + problem%CrystalAxis(:,2,1)*Hky_rot(:) + problem%CrystalAxis(:,3,1)*Hkz_rot(:)
        solution%Hky(:) = problem%CrystalAxis(:,1,2)*Hkx_rot(:) + problem%CrystalAxis(:,2,2)*Hky_rot(:) + problem%CrystalAxis(:,3,2)*Hkz_rot(:)
        solution%Hkz(:) = problem%CrystalAxis(:,1,3)*Hkx_rot(:) + problem%CrystalAxis(:,2,3)*Hky_rot(:) + problem%CrystalAxis(:,3,3)*Hkz_rot(:)
        
        deallocate(Mx_rot, My_rot, Mz_rot, Hkx_rot, Hky_rot, Hkz_rot)
    end if 
    call trace%end( "updateAnisotropy", itimer=itimer, verbose=1 )

    end subroutine updateAnisotropy    


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
  double precision , contiguous, pointer :: pot(:,:), grad(:,:,:)
!   double precision , pointer :: source(:,:), dipvec(:,:)
!   double precision , pointer :: pot(:,:), grad(:,:)
  double precision :: vol_i
  double precision :: mx, my, mz
  double precision :: eps_fmm

  class(FMM3DTree), pointer :: fmm_tree
  logical :: built_tree
  integer, save :: itimer = 0
  !------------------------------------------------

    call trace%begin( "updateDemagfieldFMM", itimer=itimer, verbose=1 )

  fourpi  = 12.566370614359172D0
  ntot = size(problem%grid%pts, dim=1)
  ier = 0

  !----------- cast Mx,My,Mz to single precision -----------
  solution%Mx_s = real(solution%Mx, SP)
  solution%My_s = real(solution%My, SP)
  solution%Mz_s = real(solution%Mz, SP)
  !--------------------------------------------------------------

  !------------------ Allocate FMM work arrays ------------------
    built_tree = .false.
    if (.not. associated(solution%fmm_tree)) then
        allocate(solution%fmm_tree)
        built_tree = .true.
        allocate(source(3, ntot))
    end if
    fmm_tree => solution%fmm_tree

  !allocate(pot(1, ntot))
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
    ! TODO - is Ms the correct scaling here?
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
    !$omp parallel do default(shared)
    do i = 1, ntot
      solution%HmX(i) = real( grad(1,1,i) / fourpi, SP )
      solution%HmY(i) = real( grad(1,2,i) / fourpi, SP )
      solution%HmZ(i) = real( grad(1,3,i) / fourpi, SP )
      ! solution%HmX(i) = real( grad(1,i) / fourpi, SP )
      ! solution%HmY(i) = real( grad(2,i) / fourpi, SP )
      ! solution%HmZ(i) = real( grad(3,i) / fourpi, SP )
    end do
    !$omp end parallel do
  else
    solution%HmX = 0.0_SP
    solution%HmY = 0.0_SP
    solution%HmZ = 0.0_SP
  end if
  !-----------------------------------------------------------------
  !------------- add correction from neighbouring tiles ---------------------
call add_near_field(problem, solution)
  !--------------------------------------------------------------------------
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
    type(MicroMagProblem),intent(inout) :: problem         !> Problem data structure    
    type(MicroMagSolution),intent(inout) :: solution    !> Solution data structure
    integer :: stat,ntot,i
    type(matrix_descr) :: descr
    real(SP) :: pref,alpha,beta
    complex(kind=4) :: alpha_c, beta_c
    real(SP), dimension(:), allocatable :: temp
    character*(100) :: prog_str 
    integer, save :: itimer = 0
    call trace%begin( "updateDemagfield", itimer=itimer, verbose=1 )
    
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
            !call displayGUIMessage( ' runing cuda demag sparse matrix multiplication' )
            pref = sngl(-1 )!* problem%Mfact)                                
            call cudaMatrVecMult_sparse( solution%Mx_s, solution%My_s, solution%Mz_s, solution%HmX, solution%HmY, solution%HmZ, pref )
            temp = solution%HmX * problem%Mfact
            solution%HmX = temp
            temp = solution%HmY * problem%Mfact
            solution%HmY = temp
            temp = solution%HmZ * problem%Mfact
            solution%HmZ = temp

            ! open (11, file="sparse_CUDA_H.bin",  &
            !         status='unknown', form='unformatted', &
            !         access='direct', recl=1*ntot)
            ! write(11,rec=1) solution%HmX 
            ! write(11,rec=2) solution%HmY
            ! write(11,rec=3) solution%HmZ
            ! close(11)
            ! error stop " test stop after cuda sparse"
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
    
    !Demagnetization tensor matrix
#if USE_FMM3D
    !------------- build neighbour demag tensor -------------------------------------------------------------------
    call BuildNeighbourDemagTensor( problem)
        !---------- if BuildNeighbourDemagTensor sets use_fmm to false then compute the demag tensor normally -----
    if (.not. problem%use_fmm) then
       call ComputeDemagfieldTensor( problem )
    end if
        !-----------------------------------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------------------------------------
#else
    call ComputeDemagfieldTensor( problem )
#endif

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
    real(DP),dimension(:,:,:,:),allocatable :: Nout,Noutave             !> Temporary storage for the demag tensor            
    integer,dimension(4) :: indx_ele
    real :: rate
    integer :: c1,c2,cr,cm
    character(10) :: prog_str
    integer, save :: itimer = 0
    call trace%begin( "ComputeDemagfieldTensor", itimer=itimer, verbose=1 )
    
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
        !call omp_set_num_threads(1)               
        if ( problem%grid%gridType .eq. gridTypeUniform ) then
            
            if (nx_ave*ny_ave*nz_ave > 1) then
                call displayGUIMessage( 'Averaging the N_tensor not supported for this tile type' )
            endif
        
            !for each element find the tensor for all evaluation points (i.e. all elements)    
            !======== NOTE ===============
            ! this parallelization seem to give issues when compiled with matlab in debug mode
            !$OMP PARALLEL DO collapse(3) SHARED(problem, nx, ny, nz, ntot) PRIVATE(ind, tile, H, Nout, pts_arr) default(none)
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
                        
                        allocate(pts_arr(ntot,3))
                        pts_arr(:,1) =  problem%grid%pts(:,1)
                        pts_arr(:,2) =  problem%grid%pts(:,2)
                        pts_arr(:,3) =  problem%grid%pts(:,3)
                        
                        call getFieldFromTiles( tile, H, pts_arr, 1, ntot, Nout, .false. )
                        
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

            
            open(21,file='Kxx.txt',status='unknown',form='formatted',action='write')
            do i=1,size(problem%Kxx,1)
                do j=1,size(problem%Kxx,2)
                    write(21,*)  problem%Kxx(i,j)
                enddo
            enddo
            close(21)
            
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



    !-------------- for debug write the dense matrices to binary files --------------
    ! TODO - add option to control this and maybe an "auxiliary" module to do this
    !open (11, file="dense_std_ref_nocuda.bin",  &
    !       status='unknown', form='unformatted', &
    !       access='direct', recl=1*ntot*ntot)
    ! write(11,rec=1) problem%Kxx
    ! write(11,rec=2) problem%Kxy
    ! write(11,rec=3) problem%Kxz
    ! write(11,rec=4) problem%Kyy
    ! write(11,rec=5) problem%Kyz
    ! write(11,rec=6) problem%Kzz
    ! close(11)
    !error stop " test stop after writing dense ref"
    !-------------- end debug write the dense matrices to binary files --------------

    call trace%end( "ComputeDemagfieldTensor", itimer=itimer, verbose=1 )
    
    end subroutine ComputeDemagfieldTensor
    
    
    
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
        

    if ( gb_problem%passExch .eq. passExchTrue) return ! Skip this if the exchange has already been passed from outside

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
    type(MicroMagProblem),intent(inout) :: problem      !> Problem data structure 
    type(MicroMagSolution),intent(inout) :: solution    !> Solution data structure
        
    integer :: stat                                    !> Status value for the various sparse matrix operations        
    type(MagTenseSparse_d) :: d2dx2, d2dy2, d2dz2      !> Sparse matrices for the double derivatives with respect to x, y and z, respectively.
    type(sparse_matrix_t) :: tmp                       !> Temporary sparse matrices used for internal calculations
    integer :: ind, ntot,colInd,rowInd                 !> Internal counter for indexing, the total no. of elements in the current sparse matrix being manipulated    
    integer :: i,j,k,nx,ny,nz                          !> For-loop counters
    type(matrix_descr) :: descr                        !> Describes a sparse matrix operation
    real(DP) :: const
    character*(100) :: prog_str 
    integer, save :: itimer = 0
    call trace%begin( "ComputeExchangeTerm3D_Uniform", itimer=itimer, verbose=1 )
    
    !Find the three sparse matrices for the the individual directions, respectively. Then add them to get the total matrix
    !It is assumed that the magnetization vector to operate on is in fact a single column of Mx, My and Mz respectively.
    
    nx = grid%nx
    ny = grid%ny
    nz = grid%nz
    
    !----------------------------------d^2dx^2 begins -----------------------------!
    if ( nx .gt. 1 ) then
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
                    d2dx2%values( ind ) = problem%A0_map(i-1,j,k)/(problem%A0_map(i-1,j,k)+problem%A0_map(i,j,k))
                    d2dx2%cols(ind) = colInd
                    !update where the row starts
                    d2dx2%rows_start(rowInd) = ind                           
                    ind = ind + 1
                
                    !Center point
                    d2dx2%values( ind ) = -(problem%A0_map(i-1,j,k)/(problem%A0_map(i-1,j,k)+problem%A0_map(i,j,k))+problem%A0_map(i+1,j,k)/(problem%A0_map(i+1,j,k)+problem%A0_map(i,j,k)))
                    d2dx2%cols(ind) = colInd + 1
                    ind = ind + 1
                
                    !Right-most point
                    d2dx2%values( ind ) = problem%A0_map(i+1,j,k)/(problem%A0_map(i+1,j,k)+problem%A0_map(i,j,k))
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
        d2dx2%values = d2dx2%values * 2./grid%dx**2
        
        !Create the sparse matrix for the d^2dx^2
        stat = mkl_sparse_d_create_csr ( d2dx2%A, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, d2dx2%rows_start, d2dx2%rows_end, d2dx2%cols, d2dx2%values)
    
    endif
    
    !----------------------------------d^2dx^2 ends ---------------------------- -!
    
    
    !----------------------------------d^2dy^2 begins ----------------------------!
    if ( ny .gt. 1 ) then
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
                    d2dy2%values(ind) = problem%A0_map(i,j-1,k)/(problem%A0_map(i,j-1,k)+problem%A0_map(i,j,k))
                    d2dy2%cols(ind) = colInd-nx
                    d2dy2%rows_start(rowInd) = ind
                    !increment to next element
                    ind = ind + 1  
                
                    !central value
                    d2dy2%values(ind) = -(problem%A0_map(i,j-1,k)/(problem%A0_map(i,j-1,k)+problem%A0_map(i,j,k))+problem%A0_map(i,j+1,k)/(problem%A0_map(i,j+1,k)+problem%A0_map(i,j,k)))
                    d2dy2%cols(ind) = colInd
            
                    !increment to next element
                    ind = ind + 1
            
                    !upper value
                    d2dy2%values(ind) = problem%A0_map(i,j+1,k)/(problem%A0_map(i,j+1,k)+problem%A0_map(i,j,k))
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
        d2dy2%values = d2dy2%values * 2./grid%dy**2
        
        !Create the sparse matrix for the d^2dy^2
        stat = mkl_sparse_d_create_csr ( d2dy2%A, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, d2dy2%rows_start, d2dy2%rows_end, d2dy2%cols, d2dy2%values)
    endif
    
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
                    d2dz2%values(ind) = problem%A0_map(i,j,k-1)/(problem%A0_map(i,j,k-1)+problem%A0_map(i,j,k))
                    d2dz2%cols(ind) = colInd - nx * ny
                    d2dz2%rows_start(rowInd) = ind
            
                    !increment position
                    ind = ind + 1
                
                    !central value
                    d2dz2%values(ind) = -(problem%A0_map(i,j,k-1)/(problem%A0_map(i,j,k-1)+problem%A0_map(i,j,k))+problem%A0_map(i,j,k+1)/(problem%A0_map(i,j,k+1)+problem%A0_map(i,j,k)))
                    d2dz2%cols(ind) = colInd
                            
                    !increment position
                    ind = ind + 1
                
                    !right-most value
                    d2dz2%values(ind) = problem%A0_map(i,j,k+1)/(problem%A0_map(i,j,k+1)+problem%A0_map(i,j,k))
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
        d2dz2%values = d2dz2%values * 2./grid%dz**2
        
        !Create the sparse matrix for the d^2dz^2
        stat = mkl_sparse_d_create_csr ( d2dz2%A, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, d2dz2%rows_start, d2dz2%rows_end, d2dz2%cols, d2dz2%values)
    endif
    
    !----------------------------------d^2dz^2 ends ----------------------------!
    
            
    !Finally, add up the three diagonals and store in the output sparse matrix, A
    !store the results temporarily in tmp4
    
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
        enddo
        
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
        enddo
    else !General anisotropy
        do i = 1,size(problem%Ms)
            problem%Kfact_arr(i,:,:) = problem%K0_arr(i,:,:) / ( mu0 * problem%Ms(i) )
        enddo
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
  real(SP), contiguous, pointer :: mxm(:), mym(:), mzm(:)
  real(SP), contiguous, pointer :: hx_tmp(:), hy_tmp(:), hz_tmp(:), temp(:)
  integer :: ntot
  real(SP) :: pref
  integer :: stat
  real(SP) :: alpha, beta
  integer, save :: itimer = 0
  call trace%begin( "add_near_field", itimer=itimer, verbose=1 )


  ntot = size(problem%grid%pts, dim=1)


#if USE_CUDA
    allocate( hx_tmp(ntot), hy_tmp(ntot), hz_tmp(ntot) )
    hx_tmp = 0.0_SP
    hy_tmp = 0.0_SP
    hz_tmp = 0.0_SP

    pref = sngl(-1)
    call cudaMatrVecMult_sparse( solution%Mx_s , solution%My_s , solution%Mz_s , hx_tmp, hy_tmp, hz_tmp, pref )

    solution%HmX = solution%HmX - hx_tmp  * problem%Mfact
    solution%HmY = solution%HmY - hy_tmp  * problem%Mfact
    solution%HmZ = solution%HmZ - hz_tmp  * problem%Mfact

    deallocate(hx_tmp, hy_tmp, hz_tmp)

#else
    allocate(mxm(ntot), mym(ntot), mzm(ntot))
    allocate(temp(ntot))

    ! Pre-scale magnetisation by Ms (as you do for CUDA)
    mxm = solution%Mx_s * real(problem%Ms, SP)
    mym = solution%My_s * real(problem%Ms, SP)
    mzm = solution%Mz_s * real(problem%Ms, SP)

    alpha = 1.0_SP
    ! ---------------- Hx correction = xx*Mx + xy*My + xz*Mz ----------------
    beta = 0.0_SP
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(1)%A, problem%K_fmm_descr_s, mxm, beta, temp)
    beta = 1.0_SP
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(2)%A, problem%K_fmm_descr_s, mym, beta, temp)
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(3)%A, problem%K_fmm_descr_s, mzm, beta, temp)

    solution%HmX = solution%HmX - temp
    ! ---------------- Hy correction = xy*Mx + yy*My + yz*Mz ----------------
    beta = 0.0_SP
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(2)%A, problem%K_fmm_descr_s, mxm, beta, temp)
    beta = 1.0_SP
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(4)%A, problem%K_fmm_descr_s, mym, beta, temp)
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(5)%A, problem%K_fmm_descr_s, mzm, beta, temp)

    solution%HmY = solution%HmY - temp
    ! ---------------- Hz correction = xz*Mx + yz*My + zz*Mz ----------------
    beta = 0.0_SP
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(3)%A, problem%K_fmm_descr_s, mxm, beta, temp)
    beta = 1.0_SP
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(5)%A, problem%K_fmm_descr_s, mym, beta, temp)
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%K_fmm_s(6)%A, problem%K_fmm_descr_s, mzm, beta, temp)

    solution%HmZ = solution%HmZ - temp
    deallocate(mxm, mym, mzm, temp)
#endif

    call trace%end( "add_near_field", itimer=itimer, verbose=1 )
end subroutine add_near_field

        !------------------------------------------------------------------------------------------------
    !============================================================================================================




    
end module LandauLifshitzSolution
    