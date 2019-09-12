
    !>Overall module for interacting with various Fortran implementations of different ODE solvers
module ODE_Solvers
    
    use rksuite_90

implicit none
    
    contains

    !> General entry point for getting MagTense to solve a set of ODE's
    subroutine MagTense_ODE()
    integer :: neq, nt
    real,allocatable,dimension(:) :: t,t_out,ystart
    real,allocatable,dimension(:,:) :: y_out, yderiv_out
    
    integer :: i
    
    ystart = 0.
    
    neq = 1
    nt = 10
    
    allocate(t(nt),t_out(nt),y_out(neq,nt),yderiv_out(neq,nt),ystart(neq))
    t(1) = 0
    ystart(:) = 0
    do i=2,nt
        t(i) = t(i-1) + 0.1
    enddo
    
    
    !Run a test case
    call MagTense_ODE_RKSuite( neq, t, nt, ystart, t_out, y_out, yderiv_out )
    
    deallocate(t_out,y_out,yderiv_out,t,ystart)
    
    end subroutine MagTense_ODE
    
    function const_acc( t, y )
    real, intent(in) :: t
    real,dimension(:), intent(in) :: y
    real,dimension(size(y)) :: const_acc
    
    const_acc(:) = t
    
    end function const_acc
    
        
    !---------------------------------------------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Specific implementation for solving ODE's using the RK suite.
    !> @param[in] neq no. of equations
    !> @param[in] t, time array of size nt
    !> @param[in] nt, no. of elements in t
    !> @param[in] ystart, the initial condition of the y-vector (size neq)
    !---------------------------------------------------------------------------
    subroutine MagTense_ODE_RKSuite( neq, t, nt, ystart,  t_out, y_out, yderiv_out )
    
    integer,intent(in) :: neq,nt                    !>Input no. of equations and no. of time steps
    real,dimension(nt),intent(in) :: t              !>Input time array, size nt
    real,dimension(neq),intent(in) :: ystart        !>Input initial conditions (y at t=0), size neq
    real,dimension(nt),intent(inout) :: t_out       !>Array for the actual time values where the solution was found
    real,dimension(neq,nt),intent(inout) :: y_out       !>Array returning the solution at the times in t_out
    real,dimension(neq,nt),intent(inout) :: yderiv_out  !>Array returning dy/dt at the times in t_out
    !real,dimension(neq),intent(in) :: f             !>The function that gived the dy/dt array at a specific t and y value
    
    real,parameter :: tol = 0.0001              !>Parameter for the relative tolerance. Will be parameterized when the code is running properly
    real,dimension(:),allocatable :: thres      !>arrays used by the initiater     
    
    character(len=1) :: task,method             !>Which version of the solver to use. = 'u' or 'U' for normal and 'C' or 'c' for complicated, Which RK method to use. 1 = RK23, 2 = RK45 and 3 = RK78
    logical :: errass,message                   !>whether to assess the true error or not, give message on errors
    real :: hstart                              !>Whether the code should choose the size of the first step. Set to 0.0d if so (recommended)    
    type(rk_comm_real_1d) :: setup_comm         !>Stores all the stuff used by setup
    integer :: flag                             !>Flag indicating how the integration went
    integer :: i                                !>Counter variable
    
    !Perform allocations. 
    allocate(thres(neq))
    
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
    
    !Call the integrator
    do i=2,nt
        call range_integrate( setup_comm, const_acc, t(i), t_out(i), y_out(:,i), yderiv_out(:,i), flag )
    enddo
    
    !Clean up
    deallocate(thres)
    
    end subroutine MagTense_ODE_RKSuite
    
    
end module ODE_Solvers