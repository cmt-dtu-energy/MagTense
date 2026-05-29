module TilePlanarCoilTensor
    use QuadPack
    use SpecialFunctions
    use TileTensorHelperFunctions
    
    implicit none
    
    contains
    
    
    !::Based on the solution from Bartberger, J. Appl. Phys. 21, 1108 (1950)
    subroutine getTensorFromSingleLoop( dat, R_loop, r, z, Nout )
    real,intent(in) :: R_loop,r,z
    real,intent(inout),dimension(2) :: Nout
    class(dataCollectionBase), intent(inout), target :: dat
    
    real :: H0,s_sq,I1,I2
    procedure (f_int_dat), pointer :: f_ptr => null ()
    
        !::Tensor at the center of the current loop
	    H0 = 1. / ( 2. * R_loop )
    
        !::Parameter for        
        s_sq = ( R_loop**2 + z**2 + r**2 ) / R_loop**2

        dat%b = 2 * R_loop * r / ( R_loop**2 + z**2 + r**2 )

        
        f_ptr => int_fctI1
        call qags_x ( f_ptr, dat, 0., pi, dat%epsabs, dat%epsrel, I1, dat%abserr_x, dat%neval_x, dat%ier_x )

        I1 = 1./pi * I1

        f_ptr => int_fctI2
        call qags_x ( f_ptr, dat, 0., pi, dat%epsabs, dat%epsrel, I2, dat%abserr_x, dat%neval_x, dat%ier_x )
        I2 = 1./pi * I2
            
        !::z-component
        Nout(2) = H0 / s_sq**(3./2.) * ( I1 - r / R_loop * I2 )
        !::r-component
        Nout(1) = H0 * z / ( s_sq**(3./2.) * R_loop ) * I2
            
        
    end subroutine getTensorFromSingleLoop

    
    function int_fctI1( theta, dat )
    real,intent(in) :: theta
    class(dataCollectionBase), intent(inout), target :: dat
    real :: int_fctI1
    
    int_fctI1 = 1. / ( 1 - dat%b * cos(theta) )**(3./2.)
    
    end function int_fctI1
    
    
    function int_fctI2( theta, dat )
    real,intent(in) :: theta
    class(dataCollectionBase), intent(inout), target :: dat
    real :: int_fctI2
    
    int_fctI2 = cos(theta) / ( 1 - dat%b * cos(theta) )**(3./2.)
    
    end function int_fctI2
    
end module
    