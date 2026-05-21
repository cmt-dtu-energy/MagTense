!>Overall module for interacting with various Fortran implementations of different ODE solvers
module ODE_Solvers
    
    !======= Inclusions ===========
    use rksuite_90
    use integrationDataTypes
    use SPECIALFUNCTIONS
    use IO_GENERAL

    !======= Declarations =========
    implicit none

    !integer,parameter :: ODE_Solver_RKSUITE=1, ODE_Solver_CVODE=2
    procedure(dydt_fct), pointer :: MTdmdt  !>Input function pointer for the function to be integrated
    real, allocatable, dimension(:) :: MTy_out, MTf_vec
    private MTdmdt, MTy_out, MTf_vec

    contains
    
    !---------------------------------------------------------------------------    
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> General entry point for getting MagTense to solve a set of ODE's
    !> @param[in] fct procedure pointer to the function that gives dy_i/dt for i=1...n
    !> where n is the no. of equations
    !> @param[in] t array of size m that holds the times at which the solution is requested
    !> @param[in] y0 array of size n holding the initial values of y_i 
    !> @param[inout] t_out array of size m (pre-allocated) that holds the times where y_i are actually found
    !> @param[inout] y_out array of size [n,m] (pre-allocated) holding the n y_i values at the m times
    !> @param[in] callback function pointer for progress updates to Matlab
    !> @param[in] useCVODE optional flag for choosing solvers. 
    !> more parameters to come as we progress in the build-up of this function (error, options such as tolerances etc)
    !---------------------------------------------------------------------------
    subroutine MagTense_ODE( fct, t, y0, t_out, y_out, includeThermal, fct_thermal, callback, callback_display, tol, thres_value, useCVODE, t_conv, conv_tol, rksuite_verbose )
    
        !======= Inclusions ===========
        use, intrinsic :: iso_c_binding
        
        !======= Declarations =========
        procedure(dydt_fct), pointer :: fct                     !>Input function pointer for the function to be integrated
        procedure(no_argument_fct), pointer :: fct_thermal      !>Function pointer for the function that updates the thermal field
        procedure(callback_fct), pointer :: callback            !>Callback function
        real,dimension(:),intent(in) :: t, y0                   !>Requested time (size m) and initial values ofy (size n)
        real,dimension(:),intent(inout) :: t_out                !>Actual time values at which the y_i are found, size m
        real,dimension(:,:),intent(inout) :: y_out              !>Function values at the times t_out, size [n,m]
        integer,intent(in) :: callback_display                  !>Sets at what time index values Fortran displays the results in Matlab
        real,intent(in) :: tol                                  !>Relative tolerance
        real,intent(in) :: thres_value                          !>When a solution component Y(L) is less in magnitude than thres_value its set to zero
        integer,intent(in),optional :: useCVODE                 !>Flag that determines if the CVODE solver is to be used or not
        real,dimension(:),intent(in) :: t_conv                  !>Array for the time values where the solution will be checked for convergence
        real,intent(in) :: conv_tol                             !>Converge criteria on difference between magnetization at different timesteps
        logical,intent(in),optional :: rksuite_verbose                    !>If .true., RKSuite will display diagnostic messages (default .false.)
        integer :: solver_flag
        logical :: includeThermal
        
        integer :: neq, nt, nt_conv
        real, allocatable, dimension(:,:) :: yderiv_out         !>The derivative of y_i wrt t at each time step
    
        !find the no. of equations and the no. of requested timesteps
        neq = size(y0)
        nt = size(t)
        nt_conv = size(t_conv)
        
        !useCVODE defaults to true. Inelegant but functional.
        if ( present(useCVODE) ) then
            solver_flag = useCVODE
        else
            solver_flag = useCVODETrue
        endif
        
        !Call the solver
        if ( solver_flag .eq. useCVODEFalse ) then

            !Allocate the derivative output array
            allocate(yderiv_out(neq, nt))
            yderiv_out(:,:) = 0

            !Call the solver
            call MagTense_ODE_RKSuite( fct, neq, t, nt, y0, t_out, y_out, yderiv_out, includeThermal, fct_thermal, callback, callback_display, tol, thres_value, nt_conv, t_conv, conv_tol, rksuite_verbose )
            
            !clean-up
            deallocate(yderiv_out)

        else if ( solver_flag .eq. useCVODETrue ) then
#if USE_CVODE
            !Do the magic for CVODE
            MTdmdt => fct
            !internal temporary arrays
            allocate( MTy_out(neq), MTf_vec(neq) )

            !call solver...
            call MagTense_CVODEsuite( int(neq, kind=c_long), neq, real(t, kind=c_double), int(nt, kind=c_long), nt, real(y0, kind=c_double), t_out, y_out, real(tol, kind=c_double), callback, int(callback_display, kind=c_long) )

            !clean-up
            deallocate(MTy_out, MTf_vec)
#else
            call callback( 'MagTense not compiled with CVODE - exiting!', 0 )
            stop
#endif    
        endif

    end subroutine MagTense_ODE
    
    
    !---------------------------------------------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Specific implementation for solving ODE's using the RK suite.
    !> @param[in] fct the function that gives the neq ODE derivatives
    !> @param[in] neq no. of equations
    !> @param[in] t, time array of size nt
    !> @param[in] nt, no. of elements in t
    !> @param[in] ystart, the initial condition of the y-vector (size neq)
    !> @param[inout] t_out output array with the times at which y_i are found
    !> @param[inout] y_out output array with the y_i values
    !> @param[inut] yderiv_out output array with dy_i/dt at each time
    !> @param[in] callback procedure pointer to callback to Matlab for progress updates
    !---------------------------------------------------------------------------
    subroutine MagTense_ODE_RKSuite( fct, neq, t, nt, ystart,  t_out, y_out, yderiv_out, includeThermal, fct_thermal, callback, callback_display, tol, thres_value, nt_conv, t_conv, conv_tol, rksuite_verbose )
        
        !======= Declarations =========
        procedure(dydt_fct), pointer :: fct                  !>Input function pointer for the function to be integrated
        procedure(no_argument_fct), pointer :: fct_thermal   !>Function pointer for the function that updates the thermal field
        integer,intent(in) :: neq, nt, nt_conv               !>Input no. of equations, no. of time steps and no. of time steps in the check for convergence array
        real,dimension(nt),intent(in) :: t                   !>Input time array, size nt
        real,dimension(neq),intent(in) :: ystart             !>Input initial conditions (y at t=0), size neq
        real,dimension(nt),intent(inout) :: t_out            !>Array for the actual time values where the solution was found
        real,dimension(neq,nt),intent(inout) :: y_out        !>Array returning the solution at the times in t_out
        real,dimension(neq,nt),intent(inout) :: yderiv_out   !>Array returning dy/dt at the times in t_out
        procedure(callback_fct), pointer :: callback         !>Callback function
        integer,intent(in) :: callback_display               !>Sets at what time index values Fortran displays the results in Matlab
        real,intent(in) :: tol                               !>Relative tolerance
        real,intent(in) :: thres_value                       !>When a solution component Y(L) is less in magnitude than thres_value its set to zero
        real,dimension(nt_conv),intent(in) :: t_conv         !>Array for the time values where the solution will be checked for convergence
        real,intent(in) :: conv_tol                          !>Converge criteria on difference between magnetization at different timesteps
        logical,intent(in),optional :: rksuite_verbose                 !>If .true., RKSuite will display diagnostic messages (default .false.)
        
        real,dimension(:),allocatable :: thres          !>arrays used by the initiater     
        real,dimension(neq) :: y_last                   !>Array containing the solution in the last returned convergence timestep
        real,dimension(neq) :: y_step                   !>Array containing the solution in the current timestep
        real,dimension(neq) :: yderiv_step              !>Array containing dy/dt in the current timestep
        real,dimension(:),allocatable :: y_norm         !>Array of vector magnitudes used for normalising y_step. neq/3 = ntot
        real,dimension(nt+nt_conv) :: t_comb            !>The concatenated time array of the output times and the convergence times
        real,dimension(nt+nt_conv) :: t_comb_out        !>The concatenated time array of the output times and the convergence times
        real,dimension(:),allocatable :: t_comb_unique  !>The concatenated time array of the output times and the convergence times, only unique values
        integer,dimension(:),allocatable :: ind         !>The indices for sort
        
        character(len=1) :: task,method             !>Which version of the solver to use. = 'u' or 'U' for normal and 'C' or 'c' for complicated, Which RK method to use. 1 = RK23, 2 = RK45 and 3 = RK78
        logical :: errass,message                   !>whether to assess the true error or not, give message on errors
        real :: hstart                              !>Whether the code should choose the size of the first step. Set to 0.0d if so (recommended)    
        real :: t_step                              !>The current time at the end of a time step
        real :: conv_error                          !>The maximum error in the current time step
        integer :: ntot                             !>Total number of micromagnetic cells (ntot = neq/3)
        type(rk_comm_real_1d) :: setup_comm         !>Stores all the stuff used by setup
        integer :: flag                             !>Flag indicating how the integration went
        integer :: message_flag                     !>Flag indicating if the error message has already been displayed to the user
        integer :: i, k                             !>Counter variable
        character*(100) :: prog_str                 !>Variable holding the output string
        logical :: includeThermal                   !>Whether to normalise result and update thermal field after each timestep
        !integer,parameter :: n_write=100
        !Perform allocations. 
        allocate(thres(neq))
        
        !set thres to the low value
        thres(:) = thres_value
        
        !Set the method to RK45 as default (will be parameterized later as the code evolves)
        !L or l for 23, M or m for 45 and H or h for 67
        method = 'M'
        !Set default solver to normal (R or r) or complex (S or s)
        task = 'R'
        !Set the assessing of the true error to not be done
        errass = .false.
        !Set the flag so that the code selects the starting step size
        hstart = 0.0
        !Set output message based on rksuite_verbose flag (default .false. to suppress RKSuite diagnostic messages)
        if ( present(rksuite_verbose) ) then
            message = rksuite_verbose
        else
            message = .false.
        endif
        
        !Call the setup function in order to initiate the solver    
        call setup(  setup_comm, t(1), ystart, t(nt), tol, thres, method, task, errass, hstart, message)

        !Calculate the magnetization to use for the first convergence calculation
        y_last = ystart
        
        !Concertinate the t_out and t_conv arrays
        t_comb(1:size(t)) = t
        t_comb((size(t)+1):(nt+nt_conv)) = t_conv
        
        !Remove the duplicate elements, if any from the t_comb array
        k = 1
        call simple_unique( t_comb, t_comb_out, k)
        
        !Sort the array
        allocate(t_comb_unique(k))
        allocate(ind(k))
        call simple_sort( t_comb_out(1:k), t_comb_unique, ind )
                
        !First time is the same as the input
        t_out(1) = t(1)
        y_out(:,1) = ystart
        y_step = ystart         !Initialize y_step with the starting magnetization


        !Allocate norm array
        ntot = neq/3
        if ( includeThermal ) then
            allocate( y_norm(ntot) )
            call fct_thermal()  ! Generate the stochastic thermal field
        end if
        
        k = 2
        !Call the integrator
        do i=2,size(t_comb_unique)                
            call range_integrate( setup_comm, fct, t_comb_unique(i), t_step, y_step, yderiv_step, flag )

            if ( includeThermal ) then
                y_norm = sqrt(y_step(1:ntot)**2 + y_step(ntot+1:2*ntot)**2 + y_step(2*ntot+1:3*ntot)**2)
                y_step(1:ntot) = y_step(1:ntot)/y_norm                     ! Normalise x-component
                y_step(ntot+1:2*ntot) = y_step(ntot+1:2*ntot)/y_norm       ! Normalise y-component
                y_step(2*ntot+1:3*ntot) = y_step(2*ntot+1:3*ntot)/y_norm   ! Normalise z-component
                call fct_thermal()  ! Update the stochastic thermal field
            end if
            
            message_flag = 0
            !Keep integrating until we reach the target time
            do while (t_step .lt. t_comb_unique(i))
                call range_integrate( setup_comm, fct, t_comb_unique(i), t_step, y_step, yderiv_step, flag )

                !Check flag and handle warnings/errors
                if (flag .eq. 1) then
                    !Success - continue
                    continue
                else if (flag .eq. 2) then
                    !Inefficiency warning - continue but log
                    if ( mod(i, callback_display) .eq. 0 ) then
                        write(prog_str,'(A)') 'Warning: Integration is inefficient. Consider using METHOD = M with interpolation.'
                        call callback( prog_str, -1 )
                    endif
                else if (flag .eq. 3) then
                    !Too much work warning - continue but log
                    if ( mod(i, callback_display) .eq. 0 ) then
                        write(prog_str,'(A)') 'Warning: Approximately 15000 function evaluations used. Continuing...'
                        call callback( prog_str, -1 )
                    endif
                else if (flag .eq. 4) then
                    !Stiffness detected - RKSuite does not have a stiff solver, so we log warning and continue
                    if ( mod(i, callback_display) .eq. 0 ) then
                        if (message_flag .eq. 0) then
                            write(prog_str,'(A)') 'Warning: Stiff problem detected. Consider using the CVODE solver.'
                            call callback( prog_str, -1 )
                            message_flag = 1
                        endif
                    endif
                else if (flag .eq. 5 .or. flag .eq. 6) then
                    !Fatal error - cannot continue
                    write(prog_str,'(A,I1,A)') 'Error: Integration failed with flag = ', flag, '. Stopping integration.'
                    call callback( prog_str, -1 )
                    exit
                endif

                !Check if we've reached the target time (within tolerance)
                if (abs(t_step - t_comb_unique(i)) .lt. 1.0e-12) then
                    exit
                endif
            enddo

            if ( mod(i, callback_display) .eq. 0 ) then
                write(prog_str,'(A6, F8.2, A15, I4.1, A1, I4.1)') 'Time: ', t_step*1e9, ' ns, i.e. step ', i, '/', size(t_comb_unique)
                call callback( prog_str, -1 )
            endif
            
            !Check if the time which the solution is returned is part of the array that checks for converge or if it part of the times where the simulation is to be saved
            !ind = findloc(t_conv,t_comb_unique(i))
            !if (maxval(ind) .gt. 0) then !If the time is part of the converge array, we check for converge
            !    conv_error = maxval(abs(y_step-y_last))
            !    if (conv_error < conv_tol) then
            !        !Save the current state before exiting
            !        y_out(:,k) = y_step
            !        yderiv_out(:,k) = yderiv_step
            !        t_out(k) = t_step
            !        exit
            !    endif
            !    !Save the new magnetization to use for the next converge calculation
            !    y_last = y_step
            !endif
            
            !Check if we should save the result in the array to be returned
            !ind = findloc(t,t_comb_unique(i))
            !if (maxval(ind) .gt. 0) then !If the time is part of the save array
                y_out(:,k) = y_step
                yderiv_out(:,k) = yderiv_step
                t_out(k) = t_step
                k = k+1
            !endif
        enddo
        !Clean up
        deallocate(thres)
        deallocate(t_comb_unique, ind)
        call collect_garbage(setup_comm)
    
    end subroutine MagTense_ODE_RKSuite
    
#if USE_CVODE
    !---------------------------------------------------------------------------
    !> @author Emil B. Poulsen, emib@dtu.dk, DTU, 2019
    !> @brief
    !> Specific implementation for solving ODE's using the CVODE suite.
    !> @param[in] neq no. of equations
    !> @param[in] t, time array of size nt
    !> @param[in] nt, no. of elements in t
    !> @param[in] ystart, the initial condition of the y-vector (size neq)
    !> @param[inout] t_out output array with the times at which y_i are found
    !> @param[inout] y_out output array with the y_i values
    !---------------------------------------------------------------------------
    subroutine MagTense_CVODEsuite( neq, neq_f, t, nt, nt_f, ystart, t_inout, y_out, rtol, callback, callback_display )
	use, intrinsic :: iso_c_binding
    
    use fsundials_core_mod           ! Fortran interface to data types and constants
    use fcvode_mod                   ! Fortran interface to CVODE
    use fnvector_serial_mod          ! Fortran interface to serial N_Vector
    use fsunmatrix_dense_mod         ! Fortran interface to dense SUNMatrix
    use fsunlinsol_spgmr_mod         ! Fortran interface to Generalised minimum residual solver

    ! local variables
    integer(c_long), intent(in) :: nt                     ! number of timesteps
    integer(c_long), intent(in) :: neq                    ! number of equations
    integer, intent(in) :: nt_f                           ! number of timesteps
    integer, intent(in) :: neq_f                          ! number of equations
    real(c_double), dimension(nt), intent(in) :: t        ! initial time
    real, dimension(nt_f), intent(inout) :: t_inout
    real(c_double), dimension(nt) :: t_out                ! output time
    real(c_double)  :: tcur(1)                            ! current time
    procedure(callback_fct), pointer :: callback          ! Callback function
    real(c_double), intent(in) :: rtol                    ! relative tolerance
    real(c_double) :: atol                                ! absolute tolerance
    real(c_double) :: hlast(1)
    real(c_double) :: dum
    real(c_double) :: dum_t(1)
    character(100) :: err_str                             ! error string for diagnostics

    integer(c_int) :: ierr                                ! error flags from C functions
    integer(c_long) :: nsteps(1)                          ! debugging variables
    integer(c_long) :: nlinconvfails(1)
    integer(c_long) :: linconvfails(1)

    integer(c_long) :: outstep                            ! output loop counter
    integer(c_long), intent(in) :: callback_display       ! Sets at what time index values Fortran displays the results in Matlab
    
    type(N_Vector), pointer :: sunvec_y                   ! sundials vector
    type(SUNMatrix), pointer :: sunmat_A                  ! sundials matrix
    type(SUNLinearSolver), pointer :: sunls               ! sundials linear solver
    type(c_ptr) :: cvode_mem                              ! CVODE memory
    type(c_ptr) :: user_data                              ! Data for the RhsFn function
    type(c_ptr) :: ctx                                    ! SUNContext type introduced in SUNDIALS 6.0.0

    ! solution vector, neq is set in the ode_functions module
    real, dimension(neq_f, nt_f), intent(inout) :: y_out
    real(c_double), dimension(neq), intent(in) :: ystart
    real(c_double), dimension(neq) :: y_cur, y_norm
    real(c_double) :: max_norm_dev

    !======= Internals ============
    ! create the SUNDIALS context
    ierr = FSUNContext_Create(SUN_COMM_NULL, ctx)
    ! set relative and absolute tolerances
    atol = rtol !Used to be 1.0d-10

    ! initialize solution vector
    y_cur = ystart
    
    ! create a serial vector
    sunvec_y => FN_VMake_Serial(neq, y_cur, ctx)
    if (.not. associated(sunvec_y)) then
        print *, 'ERROR: sunvec = NULL'
        stop 1
    end if
    
    ! create a dense matrix
    sunmat_A => FSUNDenseMatrix(neq, neq, ctx)
    ! sunmat_A = FSUNDenseMatrix(neq, neq, ctx)
    if (.not. associated(sunmat_A)) then
        print *, 'ERROR: sunmat = NULL'
        stop 1
    end if

    ! create a dense linear solver. 0 -- no preconditioning, lmin set to 300 according to
    ! D. Suess -- Time resolved micromagnetics using a preconditioned time integration method
    ! https://doi.org/10.1016/S0304-8853(02)00341-4
    sunls => FSUNLinSol_SPGMR(sunvec_y, 0, 300, ctx)
    if (.not. associated(sunls)) then
        print *, 'ERROR: sunls = NULL'
        stop 1
    end if
    
    ! create CVode memory. CV_BDF recommended for stiff problems (along with default newton iteration)
    ! CV_ADAMS recommended for nonstiff problems (along with fixed point solver)
    cvode_mem = FCVodeCreate(CV_BDF, ctx)
    if (.not. c_associated(cvode_mem)) then
        print *, 'ERROR: cvode_mem = NULL'
        stop 1
    end if

    ! initiate user data. Just neq for now
    user_data = transfer(neq, user_data)
    ierr = FCVodeSetUserData(cvode_mem, user_data);
    if (ierr /= 0) then
		call CVODE_error('Error in FCVodeSetUserData, ierr = ', ierr, callback )
	    stop
    end if
  
    ! initialize CVode
    ierr = FCVodeInit(cvode_mem, c_funloc(RhsFn), t(1), sunvec_y)
    if (ierr /= 0) then
	    call CVODE_error('Error in FCVodeInit, ierr = ', ierr, callback )
	    stop
    end if

    ierr = FCVodeSStolerances(cvode_mem, rtol, atol)
    if (ierr /= 0) then
	    call CVODE_error('Error in FCVodeSStolerances, ierr = ', ierr, callback ) 
	    stop
    end if
    
    ! set maximum amount of convergence fails (default: 10)
    ierr = FCVodeSetMaxConvFails(cvode_mem, 100)
    if (ierr /= 0) then
	    call CVODE_error('Error in FCVodeSetMaxConvFails, ierr = ', ierr, callback ) 
	    stop
    end if
    
    ! set maximum number of steps (default: 500)
    ierr = FCVodeSetMaxNumSteps(cvode_mem, 15000)
    if (ierr /= 0) then
	    call CVODE_error('Error in FCVodeSetMaxNumSteps, ierr = ', ierr, callback ) 
		stop
    end if
    
    ! set maximum order of BDF method (default: 5)
    ierr = FCVodeSetMaxOrd(cvode_mem, 5)    !Used to be 2
    if (ierr /= 0) then
	    call CVODE_error('Error in FCVodeSetMaxOrd, ierr = ', ierr, callback ) 
		stop
    end if
    
    ! Input because otherwise we get error that in- and output times are too close
    dum_t(1) = (t(2) - t(1)) / 2
    ierr = FCVodeSetInitStep(cvode_mem, dum_t(1))
    if (ierr /= 0) then
        call CVODE_error('Error in FCVodeSetInitStep, ierr = ', ierr, callback ) 
		stop
    end if

    ! attach linear solver
    ierr = FCVodeSetLinearSolver(cvode_mem, sunls, sunmat_A);
    if (ierr /= 0) then
	    call CVODE_error('Error in FCVodeSetLinearSolver, ierr = ', ierr, callback ) 
	    stop
    end if
    
    ! ensure FCVode doesn't overstep on the last evolution
    dum_t(1) = t(nt)
    ierr = FCVodeSetStopTime(cvode_mem, dum_t(1))
    if (ierr /= 0) then
        call CVODE_error( 'Error in FCVodeSetStopTime, ierr = ', ierr, callback  )
        write(err_str,'(A30, G10.3)') 'Attempted stop time: ', dum_t
        call CVODE_error(trim(err_str), int(nt, c_int), callback ) 
        stop
    end if
    call callback('Finished initialization, starting time steps', 0 )
    
    t_out = real(t_inout, kind=c_double)
    t_out(1) = t(1)
    t_inout(1) = real(t(1))
    y_out(:, 1) = real(y_cur)

    do outstep = 2, nt
	    ! call CVode
        if ( mod(outstep, callback_display) .eq. 0 ) then
            ! Update user
            call callback( 'Time', int(outstep, c_int) )
        endif
	    ierr = FCVode(cvode_mem, t(outstep), sunvec_y, t_out(outstep), CV_NORMAL)
        t_inout(outstep) = real(t(outstep))

        ! http://sundials.wikidot.com/return-time
	    if (ierr .lt. 0) then
            call CVODE_error('Error in FCVODE, ierr = ', ierr, callback ) 
            ierr = FCVodeGetLastStep(cvode_mem, hlast)
            if (ierr .eq. 0) then
                write(err_str,'(A30, G10.3)') 'Last stepsize, hlast = ', hlast
                call CVODE_error(err_str, 0, callback ) 
            else
                call CVODE_error('Error in FCVodeGetLastStep, ierr = ', ierr, callback )
            endif
            ierr = FCVodeGetCurrentTime(cvode_mem, dum_t(1));
            if (ierr .eq. 0) then
                write(err_str,'(A30, G10.3)') 'Last time, t = ', dum_t
                call CVODE_error(err_str, 0, callback ) 
            else
                call CVODE_error('Error in FCVodeGetCurrentTime, ierr = ', ierr, callback )
            endif
            ierr = FCVodeGetNumNonlinSolvConvFails(cvode_mem, nlinconvfails)
            if (ierr .eq. 0) then
                call CVODE_error('Nonlinear solver fails: ', int(nlinconvfails(1)), callback )
            else
                call CVODE_error('Error in FCVodeGetNumNonlinSolvConvFails, ierr = ', ierr, callback )
            endif
            ierr = FCVodeGetNumLinConvFails(cvode_mem, linconvfails)
            if (ierr .eq. 0) then
              call CVODE_error('Linear solver fails: ', int(linconvfails(1)), callback )
            else
                call CVODE_error('Error in FCVodeGetNumLinConvFails, ierr = ', ierr, callback )
            endif
            ierr = FCVodeGetNumSteps(cvode_mem, nsteps)
            if (ierr .eq. 0) then
                call CVODE_error('Number of steps taken: ', int(nsteps(1)), callback ) 
            else
                call CVODE_error('Error in FCVodeGetNumSteps, ierr = ', ierr, callback )
            endif

		    exit
        endif
        
        ! Check norm (under the assumption that y_cur is a magnetization vector containing x-, y- and z- components)
        call Norm_error(neq, y_cur, y_norm, max_norm_dev)
        ! If norm deviates from 1 by more than relative tolerance, rescale and reinitiate CVODE
        if (max_norm_dev .gt. rtol) then
            write(err_str,'(A26, G10.3, A14)') 'Norm - 1 > rtol (max dev = ', max_norm_dev, '). Rescale @ t'
            call CVODE_error(err_str, int(outstep, c_int), callback ) 
            
            ! and reinitialize solution vector ...
            y_cur = y_cur / y_norm
            ierr = FCVodeReInit(cvode_mem, t(outstep), sunvec_y)
            if (ierr /= 0) then
	            call CVODE_error('Error in FCVodeReInit, ierr = ', ierr, callback )
	            stop
            endif
        endif
        y_out(:, outstep) = real(y_cur)
    enddo

    ! diagnostics output
    ! call CVodeStats(cvode_mem)

    ! clean up
    call FCVodeFree(cvode_mem)
    call FSUNMatDestroy(sunmat_A)
    call FN_VDestroy(sunvec_y)
    ierr = FSUNLinSolFree(sunls)
    ierr = FSUNContext_Free(ctx)
	
    end subroutine MagTense_CVODEsuite
    
    
    !Wrapper for CVODE to call fct
    integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
        result(ierr) bind(C, name='RhsFn')

        !======= Inclusions ===========
        use, intrinsic :: iso_c_binding
        use fsundials_core_mod
        use fnvector_serial_mod

        !======= Declarations =========
        implicit none

        ! calling variables
        real(c_double), value :: tn        ! current time
        type(N_Vector)        :: sunvec_y  ! solution N_Vector
        type(N_Vector)        :: sunvec_f  ! rhs N_Vector
        type(c_ptr), value    :: user_data ! user-defined data. So far only neq.

        ! pointers to data in SUNDAILS vectors
        integer(c_long) :: i, neq
        real(c_double), pointer :: yvec(:)
        real(c_double), pointer :: fvec(:)

        neq = transfer(user_data, neq)

        ! get data arrays from SUNDIALS vectors
        yvec => FN_VGetArrayPointer(sunvec_y)
        fvec => FN_VGetArrayPointer(sunvec_f)

        !copy the yvec
        MTy_out(1:neq) = yvec(1:neq)
        
        ! fill RHS vector
        ! fvec(1) = lamda * yvec(1) + 1.0 / (1.0 + tn * tn) - lamda * atan(tn);
        call MTdmdt( real(tn), MTy_out, MTf_vec )  

        !copy the result back
        fvec(1:neq) = MTf_vec(1:neq)
        
        ! return success
        ierr = 0
        return

    end function RhsFn
	
	subroutine CVODE_error( msg, ierr, callback )
        use, intrinsic :: iso_c_binding
		integer(c_int), intent(in) :: ierr
		character(*), intent(in) :: msg
        logical :: exist
        procedure(callback_fct), pointer :: callback         !> Callback function
		
        !FILE *fp;
        
        call callback(msg, ierr)
        !inquire(file="error.txt", exist=exist)
        !if(exist) then
        !    open(12,file='error.txt',status='old',access='append',form='formatted',action='write')    
        !else
        !    open(12,file='error.txt',status='new',form='formatted',action='write')
        !end if
        !write(12,*) msg
        !write(12,*) ierr
        !close(12)
    end subroutine CVODE_error

    subroutine Norm_error(neq, y_cur, y_norm, max_norm_dev)
        use, intrinsic :: iso_c_binding
        integer(c_long), intent(in) :: neq
        real(c_double), dimension(:), intent(in) :: y_cur
        real(c_double), intent(out) :: max_norm_dev
        real(c_double) :: tmp_norm
        real(c_double),dimension(:), intent(out) :: y_norm
        integer :: k, nvar
        
        nvar = neq/3
        max_norm_dev=0.0
        do k=1,nvar
            tmp_norm=sqrt(y_cur(k)*y_cur(k)+y_cur(nvar+k)*y_cur(nvar+k)+y_cur(2*nvar+k)*y_cur(2*nvar+k))
            max_norm_dev=max(max_norm_dev,abs(tmp_norm-1))
            y_norm(k)=tmp_norm
            y_norm(nvar+k)=tmp_norm
            y_norm(2*nvar+k)=tmp_norm
        enddo
    end subroutine Norm_error

#endif

end module ODE_Solvers