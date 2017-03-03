module integrationDataTypes
implicit none



type dataCollectionBase
    real :: x1,x2,y1,y2,x,y,epsabs, epsrel,abserr_x,abserr_y
    integer :: neval_x, neval_y, ier_x, ier_y
    real*8 :: z0
    real*8,dimension(3) :: n_vec
    procedure (func), pointer, nopass :: f_ptr => null ()
    real*8,dimension(2) :: abserr_tot
    integer*8,dimension(2) :: ier,neval       
    class( dataCollectionModelBase ), pointer :: model => null ()
end type dataCollectionBase

type dataCollectionModelBase
    real :: tmp    
end type dataCollectionModelBase


abstract interface
      function func ( dat )
         import dataCollectionBase
         real :: func
         class( dataCollectionBase ) :: dat
      end function func
end interface

abstract interface
      function f_int_dat ( x, dat )
         import dataCollectionBase
         real :: f_int_dat
         real, intent (in) :: x
         class( dataCollectionBase ) :: dat
      end function f_int_dat
end interface

abstract interface
    subroutine error_handler( dat, abserr )
        import dataCollectionBase
        class( dataCollectionBase ), intent(in), target :: dat
        real, intent(inout),dimension(2) :: abserr
    end subroutine error_handler
end interface

type func_ptr
    procedure ( func ), pointer, nopass :: f_ptr => null ()    
end type func_ptr

type surf_carth
    real,dimension(2) :: x,y,z    
end type surf_carth

!::A custom type encapsulating a dataCollectionBase pointer in order to make an array
type dat_ptr
    class( dataCollectionBase ), pointer :: dat  
end type dat_ptr



    
end module integrationDataTypes