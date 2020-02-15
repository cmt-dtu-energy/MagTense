!>Overall module for interacting with various Fortran implementations of different ODE solvers
module ODE_Solvers
    
    use rksuite_90
    use integrationDataTypes
implicit none
    
    

  
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
    !> more parameters to come as we progress in the build-up of this function (error, options such as tolerances etc)
    !---------------------------------------------------------------------------
    subroutine MagTense_ODE( fct, t, y0, t_out, y_out, callback, callback_display )
    procedure(dydt_fct), pointer :: fct                     !>Input function pointer for the function to be integrated
    procedure(callback_fct), pointer :: callback            !> Callback function
    real,dimension(:),intent(in) :: t,y0                    !>requested time (size m) and initial values ofy (size n)
    real,dimension(:),intent(inout) :: t_out                !>actual time values at which the y_i are found, size m
    real,dimension(:,:),intent(inout) :: y_out              !>Function values at the times t_out, size [n,m]
    integer,intent(in) :: callback_display                  !>Sets at what time index values Fortran displays the results in Matlab
    
    integer :: neq, nt    
    real,allocatable,dimension(:,:) :: yderiv_out           !>The derivative of y_i wrt t at each time step
    
    
    
    
    !find the no. of equations and the no. of requested timesteps
    neq = size(y0)
    nt = size(t)
    
    !Allocate the derivative output array
    allocate(yderiv_out(neq,nt))
    yderiv_out(:,:) = 0
    
    
    !Call the solver
    call MagTense_ODE_RKSuite( fct, neq, t, nt, y0, t_out, y_out, yderiv_out, callback, callback_display )
    
    deallocate(yderiv_out)
    
    
    
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
    
    
end module ODE_Solvers