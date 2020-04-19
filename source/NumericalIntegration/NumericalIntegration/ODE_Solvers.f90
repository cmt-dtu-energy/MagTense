!>Overall module for interacting with various Fortran implementations of different ODE solvers
module ODE_Solvers
    
    use rksuite_90
    use integrationDataTypes
implicit none
    
    
!integer,parameter :: ODE_Solver_RKSUITE=1,ODE_Solver_CVODE=2

procedure(dydt_fct), pointer :: MTdmdt                     !>Input function pointer for the function to be integrated
real,allocatable,dimension(:) :: MTy_out,MTf_vec
private MTdmdt, MTy_out,MTf_vec
    contains
    
  
    
    !---------------------------------------------------------------------------    
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
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
    subroutine MagTense_ODE( fct, t, y0, t_out, y_out, callback, callback_display, useCVODE )
	use, intrinsic :: iso_c_binding
    procedure(dydt_fct), pointer :: fct                     !>Input function pointer for the function to be integrated
    procedure(callback_fct), pointer :: callback            !> Callback function
    real,dimension(:),intent(in) :: t,y0                    !>requested time (size m) and initial values ofy (size n)
    real,dimension(:),intent(inout) :: t_out                !>actual time values at which the y_i are found, size m
    real,dimension(:,:),intent(inout) :: y_out              !>Function values at the times t_out, size [n,m]
    integer,intent(in) :: callback_display                  !>Sets at what time index values Fortran displays the results in Matlab
    integer,intent(in),optional :: useCVODE
	integer :: solver_flag
    
    integer :: neq, nt    
    real,allocatable,dimension(:,:) :: yderiv_out           !>The derivative of y_i wrt t at each time step
    real(c_double),dimension(size(t)) :: t_out_double,t_double! = t_out
    real(c_double),dimension(size(y0),size(t)) :: y_out_double! = y_out
    real(c_double) :: t1b,t2b
    real :: t1a,t2a
    
    
    !find the no. of equations and the no. of requested timesteps
    neq = size(y0)
    nt = size(t)
    
    !Allocate the derivative output array
    allocate(yderiv_out(neq,nt))
    yderiv_out(:,:) = 0
    
    !useCVODE defaults to true. Inelegant but functional.
	if ( present(useCVODE) ) then
		solver_flag = useCVODE
	else
		solver_flag = useCVODETrue
	endif
	
	!Call the solver
    if ( solver_flag .eq. useCVODEFalse ) then
    
        call MagTense_ODE_RKSuite( fct, neq, t, nt, y0, t_out, y_out, yderiv_out, callback, callback_display )
    else if ( solver_flag .eq. useCVODETrue ) then
        !Do the magic for CVODE
        MTdmdt => fct
        !internal temporary arrays
        allocate(MTy_out(neq),MTf_vec(neq) )
        
        !call solver...
        !call MagTense_CVODEsuite( integer(neq,c_int), real(t,c_double), integer(nt,c_int), real(y0,c_double), real(t_out,c_double), real(y_out,c_double) )
        !call MagTense_CVODEsuite( int(neq,c_int), real(t,c_double), int(nt,c_int), real(y0,c_double), t_out, y_out )
        t_out_double = real(t_out,c_double)
        y_out_double = real(y_out,c_double)
        t_double = real(t,c_double)
        t1a = t(1)
        t2a = t(2)
        t1b = t_double(1)
        t2b = t_double(2)
        call MagTense_CVODEsuite( int(neq,c_int), real(t,c_double), int(nt,c_int), real(y0,c_double), t_out_double, y_out_double )
        t_out = real(t_out_double)
        y_out = real(y_out_double)
        
        !clean-up
        deallocate(MTy_out,MTf_vec)
    endif
    
    !clean-up
    deallocate(yderiv_out)
    
    
    
    end subroutine MagTense_ODE
    
    
        
    !---------------------------------------------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
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
    subroutine MagTense_ODE_RKSuite( fct, neq, t, nt, ystart,  t_out, y_out, yderiv_out, callback, callback_display )
    procedure(dydt_fct), pointer :: fct                  !>Input function pointer for the function to be integrated
    integer,intent(in) :: neq,nt                         !>Input no. of equations and no. of time steps
    real,dimension(nt),intent(in) :: t                   !>Input time array, size nt
    real,dimension(neq),intent(in) :: ystart             !>Input initial conditions (y at t=0), size neq
    real,dimension(nt),intent(inout) :: t_out            !>Array for the actual time values where the solution was found
    real,dimension(neq,nt),intent(inout) :: y_out        !>Array returning the solution at the times in t_out
    real,dimension(neq,nt),intent(inout) :: yderiv_out   !>Array returning dy/dt at the times in t_out
    procedure(callback_fct), pointer :: callback         !>Callback function
    integer,intent(in) :: callback_display               !>Sets at what time index values Fortran displays the results in Matlab
    
    real :: tol                                 !>Relative tolerance. Will be parameterized when the code is running properly
    real,dimension(:),allocatable :: thres      !>arrays used by the initiater     
    
    character(len=1) :: task,method             !>Which version of the solver to use. = 'u' or 'U' for normal and 'C' or 'c' for complicated, Which RK method to use. 1 = RK23, 2 = RK45 and 3 = RK78
    logical :: errass,message                   !>whether to assess the true error or not, give message on errors
    real :: hstart                              !>Whether the code should choose the size of the first step. Set to 0.0d if so (recommended)    
    type(rk_comm_real_1d) :: setup_comm         !>Stores all the stuff used by setup
    integer :: flag                             !>Flag indicating how the integration went
    integer :: i                                !>Counter variable
    !integer,parameter :: n_write=100
    !Perform allocations. 
    allocate(thres(neq))
    
    !tolerance
    tol = 1e-6
    
    !set thres to the low value
    thres(:) = 1e-10
    
    !Set the method to RK45 as default (will be parameterized later as the code evolves)
    !L or l for 23, M or m for 45 and H o h for 67
    method = 'M'
    !Set default solver to normal (R or r) or complex (S or s)
    task = 'R'
    !Set the assessing of the true error to not be done
    errass = .false.
    !Set the flag so that the code selects the starting step size
    hstart = 0.0
    !Set output message to true
    message = .true.
    
    !Call the setup function in order to initiate the solver    
    call setup(  setup_comm, t(1), ystart, t(nt), tol, thres, method,task,errass, hstart,message)

    !First time is the same as the input
    t_out(1) = t(1)
    y_out(:,1) = ystart
    !Call the integrator
    do i=2,nt                
        call range_integrate( setup_comm, fct, t(i), t_out(i), y_out(:,i), yderiv_out(:,i), flag )
        
        if ( mod(i,callback_display) .eq. 0 ) then
            call callback( 'Time', i )
        endif
    enddo
    
    !Clean up
    deallocate(thres)
    
    end subroutine MagTense_ODE_RKSuite
    
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
    subroutine MagTense_CVODEsuite( neq, t, nt, ystart,  t_out, y_out )
	use, intrinsic :: iso_c_binding

    use fcvode_mod             ! Fortran interface to CVODE
    use fnvector_serial_mod    ! Fortran interface to serial N_Vector
    use fsunmatrix_dense_mod      ! Fortran interface to dense SUNMatrix
    use fsunlinsol_dense_mod   ! Fortran interface to dense SUNLinearSolver
    use fsunnonlinsol_fixedpoint_mod ! Fortran interface to fixed-point nonlinear solver
    !use ode_mod               ! ODE functions

    ! local variables
    integer(c_int),intent(in) :: neq,nt       ! number of eq. and timesteps
    real(c_double),dimension(nt),intent(in) :: t     ! initial time
    real(c_double),dimension(nt),intent(inout) :: t_out       ! output time
    real(c_double) :: rtol, atol,hlast,dum1,dum2 ! relative and absolute tolerance

      integer(c_int) :: ierr, ierr2,ierr3,ierr4,ierr5       ! error flags from C functions
    integer(c_int) :: nlinconvfails,linconvfails,nsteps     ! debugging variables

      integer :: outstep           ! output loop counter
    
    type(c_ptr) :: sunvec_y      ! sundials vector
    type(c_ptr) :: sunmat_A      ! sundials matrix
    type(c_ptr) :: sunlinsol_LS  ! sundials linear solver
    type(c_ptr) :: cvode_mem     ! CVODE memory
    type(c_ptr) :: NLS  	     ! nonlinear solver
    type(c_ptr) :: user_data     ! Data for the RhsFn function

    ! solution vector, neq is set in the ode_functions module
    real(c_double),dimension(neq,nt),intent(inout) :: y_out
    real(c_double),dimension(neq),intent(in) :: ystart
    real(c_double),dimension(neq) :: y_cur

    !======= Internals ============
    dum1=t(1)
    dum2=t(2)
    ! initialize solution vector
    y_cur = ystart

    ! create a serial vector
    sunvec_y = FN_VMake_Serial(neq, y_cur)
    if (.not. c_associated(sunvec_y)) print *,'ERROR: sunvec = NULL'

    ! create a dense matrix
    sunmat_A = FSUNDenseMatrix(neq, neq);
    if (.not. c_associated(sunmat_A)) print *,'ERROR: sunmat = NULL'

    ! create a dense linear solver
    sunlinsol_LS = FSUNDenseLinearSolver(sunvec_y, sunmat_A);
    if (.not. c_associated(sunlinsol_LS)) print *,'ERROR: sunlinsol = NULL'
    
    ! create CVode memory. CV_BDF recommended for stiff problems (along with default newton iteration)
    cvode_mem = FCVodeCreate(CV_ADAMS) ! CV_ADAMS recommended for nonstiff problems (along with fixed point solver)
    if (.not. c_associated(cvode_mem)) print *,'ERROR: cvode_mem = NULL'

    ! initiate user data. Just neq for now
    user_data = transfer(neq,user_data)
    ierr = FCVodeSetUserData(cvode_mem, user_data);
    if (ierr /= 0) then
        write(*,*) 'Error in FCVodeSetUserData, ierr = ', ierr, '; halting'
	    stop
    end if
  
    ! initialize CVode
    ierr = FCVodeInit(cvode_mem, c_funloc(RhsFn), t(1), sunvec_y)
    if (ierr /= 0) then
	    write(*,*) 'Error in FCVodeInit, ierr = ', ierr, '; halting'
	    stop
    end if
    
    ! create fixed point nonlinear solver object
    NLS = FSUNNonlinSol_FixedPoint(sunvec_y, 0)
    if (.not. c_associated(NLS)) print *,'ERROR: sunnls = NULL'

    ! attache nonlinear solver object to CVode
    ierr = FCVodeSetNonlinearSolver(cvode_mem, NLS)
    if (ierr /= 0) then
        write(*,*) 'Error in FCVodeSetNonlinearSolver, ierr = ', ierr, '; halting'
        stop
    end if

    ! set relative and absolute tolerances
    rtol = 1.0d-6
    atol = 1.0d-10

    ierr = FCVodeSStolerances(cvode_mem, rtol, atol)
    if (ierr /= 0) then
	    write(*,*) 'Error in FCVodeSStolerances, ierr = ', ierr, '; halting'
	    stop
    end if
    
    ! set maximum amount of convergence fails (default: 10)
    ierr = FCVodeSetMaxConvFails(cvode_mem, 100)
    if (ierr /= 0) then
	    write(*,*) 'Error in FCVodeSetMaxConvFails, ierr = ', ierr, '; halting'
	    stop
    end if
    
    ierr = FCVodeSetMaxNumSteps(cvode_mem, 10000)
    if (ierr /= 0) then
	    write(*,*) 'Error in FCVodeSetMaxNumSteps, ierr = ', ierr, '; halting'
	    stop
    end if
    
    ierr = FCVodeSetInitStep(cvode_mem, (t(2)-t(1))/2)
    if (ierr /= 0) then
        write(*,*) 'Error in FCVodeSetInitStep, ierr = ', ierr, '; halting'
        stop
    end if

    ! attach linear solver
    ierr = FCVodeSetLinearSolver(cvode_mem, sunlinsol_LS, sunmat_A);
    if (ierr /= 0) then
	    write(*,*) 'Error in FCVodeSetLinearSolver, ierr = ', ierr, '; halting'
	    stop
    end if

    ! start time stepping
    print *, '   '
    print *, 'Finished initialization, starting time steps'

    do outstep = 1,(nt-1+1)
	    ! call CVode
        if (outstep == nt) then
            write(*,*) 'break'
            ! ensure FCVode doesn't overstep on the last evolution
            ierr = FCVodeSetStopTime(cvode_mem, t(nt))
            if (ierr /= 0) then
                write(*,*) 'Error in FCVodeSetStopTime, ierr = ', ierr, '; halting'
                stop
            end if

        endif
        
	    ierr = FCVode(cvode_mem, t(outstep), sunvec_y, t_out(outstep), CV_NORMAL)
        ! http://sundials.wikidot.com/return-time
        !ierr = FCVode(cvode_mem, t(outstep), sunvec_y, t_out(outstep), CV_NORMAL_TSTOP)
	    if (ierr /= 0) then
            rtol=t_out(1)
            atol=t_out(2)
            dum1=t(1)
            dum2=t(2)
        ierr2 = FCVodeGetLastStep(cvode_mem, hlast)
        ierr3 = FCVodeGetNumNonlinSolvConvFails(cvode_mem, nlinconvfails)
        ierr4 = FCVodeGetNumLinConvFails(cvode_mem, linconvfails)
        ierr5 = FCVodeGetNumSteps(cvode_mem, nsteps)
		write(*,*) 'Error in FCVODE, ierr = ', ierr, '; halting'
		stop
	    endif
	    y_out(:,outstep) = y_cur    ! Todo: Find a way to write directly to y_out
    enddo

    ! diagnostics output (custom function. Not written)
    !call CVodeStats(cvode_mem)

    ! clean up
    call FCVodeFree(cvode_mem)
    call FSUNLinSolFree_Dense(sunlinsol_LS)
    call FSUNMatDestroy_Dense(sunmat_A)
    call FN_VDestroy_Serial(sunvec_y)
	
	end subroutine MagTense_CVODEsuite
    
    !Wrapper for CVODE to call fct
    integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
       result(ierr) bind(C,name='RhsFn')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fnvector_serial_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn        ! current time
    type(c_ptr), value    :: sunvec_y  ! solution N_Vector
    type(c_ptr), value    :: sunvec_f  ! rhs N_Vector
    type(c_ptr), value    :: user_data ! user-defined data. So far only neq.

    ! pointers to data in SUNDAILS vectors
    real(c_double), pointer :: yvec(:)
    real(c_double), pointer :: fvec(:)

    integer(c_int) :: i, neq
    neq = transfer(user_data,neq)
    
    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    call FN_VGetData_Serial(sunvec_y, yvec)
    call FN_VGetData_Serial(sunvec_f, fvec)

    !copy the yvec
    do i=1,neq
        MTy_out(i) = yvec(i)
    enddo
    
    
    ! fill RHS vector
    !fvec(1) = lamda*yvec(1) + 1.0/(1.0+tn*tn) - lamda*atan(tn);
    call MTdmdt( tn, MTy_out, MTf_vec )  

    !copy the result back
    do i=1,neq
        fvec(i) = MTf_vec(i)
    enddo
    
    
    ! return success
    ierr = 0
    return

    end function RhsFn
end module ODE_Solvers