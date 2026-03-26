
    module integrationDataTypes
    !!@todo This file needs commenting
    implicit none

    integer,parameter :: coord_sys_carth=1,coord_sys_cyl=2,coord_sys_cone=3
    !!@todo Do NOT have useCVODETrue/-False variables both here and in MicroMagParameters.
    integer,parameter :: useCVODETrue=1,useCVODEFalse=0

    !---------------------------------------------------------------------------    
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Defines the interface of the function to be integrated.
    !> Should return n values
    !> @param[in] t the time at which the derivative is requested
    !> @param[in] y array size n holding the y_i values corresponding to the time t
    !> @param[inout] dydt array size n for the derivatives at the time t
    !---------------------------------------------------------------------------    
    abstract interface
          subroutine dydt_fct ( t, y, dydt )  
             real,intent(in) :: t
             real,dimension(:),intent(in) :: y
             real,dimension(:),intent(inout) :: dydt
         
          end subroutine dydt_fct
    end interface
    
     !---------------------------------------------------------------------------    
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Defines an interface for progress callback to Matlab
    !> Should return n values
    !> @param[in] mess message to display
    !> @param[in] prog current progress
    !---------------------------------------------------------------------------    
    abstract interface
          subroutine callback_fct ( mess, prog )  
             character(*),intent(in) :: mess
             integer,intent(in) :: prog
         
          end subroutine callback_fct
    end interface
    
    
     !---------------------------------------------------------------------------    
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2020
    !> @brief
    !> Defines an interface for message callback to Matlab
    !> Should return n values
    !> @param[in] mess message to display
    !---------------------------------------------------------------------------    
    abstract interface
          subroutine callback_mess_fct ( mess )  
             character(*),intent(in) :: mess
         
          end subroutine callback_mess_fct
    end interface
    
    
    
    type dataCollectionBase
        real :: r1, r2, theta1, theta2, z1, z2, rs, thetas, zs
        real :: x1,x2,y1,y2,x,y,z,epsabs, epsrel,abserr_x,abserr_y,b   
        integer(kind=8) :: neval_x, neval_y, ier_x, ier_y
        integer:: progCallbackCnt
        real :: z0,cone_angle,cone_z0
        real,dimension(1,3) :: n_vec
        procedure (func), pointer, nopass :: f_ptr => null ()
        procedure (func_vec), pointer, nopass :: f_ptr_vec => null ()
        procedure (func), pointer, nopass :: progCallback => null ()
        real,dimension(2) :: abserr_tot
        integer,dimension(2) :: ier,neval       
        class( dataCollectionModelBase ), pointer :: model => null ()
        integer :: coord_sys
    end type dataCollectionBase

    type dataCollectionModelBase
        real :: tmp    
    end type dataCollectionModelBase

    abstract interface
          function func ( dat )
             import dataCollectionBase
             real :: func
             class(dataCollectionBase), intent(inout), target :: dat
          end function func
    end interface

    abstract interface
          subroutine func_vec ( yy, dat, n, res)
             import dataCollectionBase
             integer,intent(in) :: n
             real,dimension(n),intent(in) :: yy
             class(dataCollectionBase), intent(inout), target :: dat
             real,dimension(n) :: res
         
          end subroutine func_vec
    end interface
    
    abstract interface
          function f_int_dat ( x, dat )
             import dataCollectionBase
             real :: f_int_dat
             real, intent (in) :: x
             class(dataCollectionBase), intent(inout), target :: dat
          end function f_int_dat
    end interface

    !! Vectorized version of the integrated function
    interface
          function f_int_dat_vec ( x, dat, n )
             import dataCollectionBase
             integer,intent(in) :: n
             real,intent (in),dimension(n) :: x
             real,dimension(n) :: f_int_dat_vec        
             class(dataCollectionBase), intent(inout), target :: dat
          end function f_int_dat_vec
    end interface

    abstract interface
        subroutine error_handler( dat, abserr )
            import dataCollectionBase
            class( dataCollectionBase ), intent(inout), target :: dat
            real, intent(inout),dimension(2) :: abserr
        end subroutine error_handler
    end interface

    type func_ptr
        procedure ( func ), pointer, nopass :: f_ptr => null ()    
    end type func_ptr

    type surf_carth
        real,dimension(2) :: x,y,z,theta
        real,dimension(3) :: r  
        integer(kind=4) :: coord,n_surfaces
        real :: cone_angle,z0
        integer,dimension(3) :: retVec
    end type surf_carth

    !! A custom type encapsulating a dataCollectionBase pointer in order to make an array
    type dat_ptr
        class(dataCollectionBase), pointer :: dat
    end type dat_ptr

    
    end module integrationDataTypes
