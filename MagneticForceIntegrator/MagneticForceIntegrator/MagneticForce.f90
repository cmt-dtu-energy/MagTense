module MagneticForce
    
    use MagForceMaxwellStressTensor
    implicit none

    
    
    contains
    
    
    
    !::
    !::Calculates the magnetic force vector, F, given the tiles (in datModel) by integrating Maxwell's stress tensor over the surface surf.
    !::
    subroutine getForce( datModel, surf, F )    
    type( surf_carth ), intent(in) :: surf
    real,dimension(3),intent(out) :: F
    
    real,parameter :: eps_abs = 1e-4,eps_rel = 1e-3
    class( magStatModel ), target :: datModel
    type( dat_ptr ), dimension(18) :: dat_arr
    integer,dimension(2) :: ier,neval
    integer,dimension(3) :: retV
    integer :: i,j,ind
    
    !::Always get all three force components
    retV(:) = 0
            
      !::Initialize each member of the dat_arr
      do i=1,3
          do j=1,surf%n_surfaces
              ind = (i-1) * surf%n_surfaces + j
              
              allocate( dat_arr(ind)%dat )
              
              dat_arr(ind)%dat%model => datModel
              
              dat_arr(ind)%dat%epsabs = eps_abs
              dat_arr(ind)%dat%epsrel = eps_rel
              dat_arr(ind)%dat%cone_angle = surf%cone_angle
              dat_arr(ind)%dat%cone_z0 = surf%z0
              
              dat_arr(ind)%dat%coord_sys = surf%coord              
              
              !::progress callback function
              dat_arr(ind)%dat%progCallback => MagForceProgCallback
          enddo
       enddo
       select case ( surf%coord )
        case ( coord_sys_cone )
        !x-component
        dat_arr(1)%dat%f_ptr => F_int_11 !theta-z surface, x part of diff. area
        dat_arr(2)%dat%f_ptr => F_int_12 !theta-z surface, y part of diff. area
        dat_arr(3)%dat%f_ptr => F_int_13 !theta-z surface, z part of diff. area
        dat_arr(4)%dat%f_ptr => F_int_13 !r-theta surface, bottom
        dat_arr(5)%dat%f_ptr => F_int_13 !r-theta surface, top
        !y-component
        dat_arr(6)%dat%f_ptr => F_int_21 !theta-z surface, x part of diff. area
        dat_arr(7)%dat%f_ptr => F_int_22 !theta-z surface, y part of diff. area
        dat_arr(8)%dat%f_ptr => F_int_23 !theta-z surface, z part of diff. area
        dat_arr(9)%dat%f_ptr => F_int_23 !r-theta surface, bottom
        dat_arr(10)%dat%f_ptr => F_int_23 !r-theta surface, top
        !z-component
        dat_arr(11)%dat%f_ptr => F_int_31 !theta-z surface, x part of diff. area
        dat_arr(12)%dat%f_ptr => F_int_32 !theta-z surface, y part of diff. area
        dat_arr(13)%dat%f_ptr => F_int_33 !theta-z surface, z part of diff. area
        dat_arr(14)%dat%f_ptr => F_int_33 !r-theta surface, bottom
        dat_arr(15)%dat%f_ptr => F_int_33!r-theta surface, top
        call surface_integral_cone( surf, dat_arr(1:15), retV, handleError, F, ier, neval )
       case ( coord_sys_cyl )
           !x-component
          dat_arr(1)%dat%f_ptr => F_int_11
          dat_arr(2)%dat%f_ptr => F_int_12
          dat_arr(3)%dat%f_ptr => F_int_13      
          dat_arr(4)%dat%f_ptr => F_int_13
          !y-component
          dat_arr(5)%dat%f_ptr => F_int_21
          dat_arr(6)%dat%f_ptr => F_int_22
          dat_arr(7)%dat%f_ptr => F_int_23
          dat_arr(8)%dat%f_ptr => F_int_23
          !z-component
          dat_arr(9)%dat%f_ptr => F_int_31
          dat_arr(10)%dat%f_ptr => F_int_32
          dat_arr(11)%dat%f_ptr => F_int_33
          dat_arr(12)%dat%f_ptr => F_int_33
          
          call surface_integral_cyl( surf, dat_arr, retV, handleError, F, ier, neval )
       case ( coord_sys_carth )
          dat_arr(1)%dat%f_ptr => F_int_11
          dat_arr(2)%dat%f_ptr => F_int_11
          dat_arr(3)%dat%f_ptr => F_int_12
          dat_arr(4)%dat%f_ptr => F_int_12
          dat_arr(5)%dat%f_ptr => F_int_13
          dat_arr(6)%dat%f_ptr => F_int_13
      
          dat_arr(7)%dat%f_ptr => F_int_21
          dat_arr(8)%dat%f_ptr => F_int_21
          dat_arr(9)%dat%f_ptr => F_int_22
          dat_arr(10)%dat%f_ptr => F_int_22
          dat_arr(11)%dat%f_ptr => F_int_23
          dat_arr(12)%dat%f_ptr => F_int_23
      
          dat_arr(13)%dat%f_ptr => F_int_31
          dat_arr(14)%dat%f_ptr => F_int_31
          dat_arr(15)%dat%f_ptr => F_int_32
          dat_arr(16)%dat%f_ptr => F_int_32
          dat_arr(17)%dat%f_ptr => F_int_33
          dat_arr(18)%dat%f_ptr => F_int_33
          
          call surface_integral_carth( surf, dat_arr, retV, handleError, F, ier, neval )
       end select
    
    
    end subroutine getForce
    
    
    
end module MagneticForce
    