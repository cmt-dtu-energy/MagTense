module MagForceMaxwellStressTensor
    use IntegrationDataTypes
    use NumInt
    use MagStat2GetSolution
    use MagForceIO
    implicit none
    
    contains

    !::handles how to store the numerical error from the integrator
    subroutine handleError( dat, abserr_tot )      
      class( dataCollectionBase ), intent(in), target :: dat
      real,intent(inout),dimension(2) :: abserr_tot
      
      if ( dat%abserr_x .lt. abserr_tot(1) ) then
        abserr_tot(1) = dat%abserr_x
      endif      
      if ( dat%abserr_y .lt. abserr_tot(2) ) then
        abserr_tot(2) = dat%abserr_y
      endif      
      if ( dat%ier_x .gt. 0 .or. dat%ier_y .gt. 0 ) then
        !stop
      endif
      end subroutine handleError
    
      !::The nine tensor components of the Maxwell stress tensor  
      function F_int_12( dat )
      real :: F_int_12
      class( dataCollectionBase ), intent(in), target :: dat
      real :: jacobi
      real,dimension(3) :: B 
          
      !::Get the field from the model solution
      call getB_field( dat, B, jacobi )
                                 
      !::Calculate the off diagonal component. 
       F_int_12 = 1. /mu0 * jacobi * B(1) * B(2)
    
       return

      end function F_int_12
      
      function F_int_13( dat )
      real :: F_int_13
      class( dataCollectionBase ), intent(in), target :: dat
      
      real,dimension(3) :: B
      real :: jacobi
          
      !::Get the field from the model solution
      call getB_field( dat, B, jacobi )
                                 
      !::Calculate the off diagonal component
       F_int_13 = 1. /mu0 * jacobi * B(1) * B(3)
    
       return
      end function F_int_13
      
      function F_int_23( dat )
      real :: F_int_23
      class( dataCollectionBase ), intent(in), target :: dat
      
      real,dimension(3) :: B
      real :: jacobi
          
      !::Get the field from the model solution
      call getB_field( dat, B, jacobi  )
                                 
      !::Calculate the off diagonal component
       F_int_23 = 1. /mu0 * jacobi * B(2) * B(3)
    
       return

      end function F_int_23
      
      function F_int_21( dat )
      real :: F_int_21
      class( dataCollectionBase ), intent(in), target :: dat
      
      F_int_21 = F_int_12( dat )
      
      end function F_int_21
      
      function F_int_31( dat )
      real :: F_int_31
      class( dataCollectionBase ), intent(in), target :: dat
      
      F_int_31 = F_int_13( dat )
      
      end function F_int_31
      
      function F_int_32( dat )
      real :: F_int_32
      class( dataCollectionBase ), intent(in), target :: dat
      
      F_int_32 = F_int_23( dat )
      
      end function F_int_32
                  
      function F_int_11( dat )
      real :: F_int_11
      class( dataCollectionBase ), intent(in), target :: dat
      real,dimension(3) :: B
      real :: Bnorm,jacobi
      
      call getB_field( dat, B, jacobi )
      
      Bnorm = sqrt( B(1)**2 + B(2)**2 + B(3)**2 )
                                 
      F_int_11 = 1./mu0 * jacobi * ( B(1)**2 - 0.5 * Bnorm**2 )
      
      
      end function F_int_11
      
      function F_int_22( dat )
      real :: F_int_22
      class( dataCollectionBase ), intent(in), target :: dat
      real,dimension(3) :: B
      real :: Bnorm, jacobi
      
      call getB_field( dat, B, jacobi )
      
      Bnorm = sqrt( B(1)**2 + B(2)**2 + B(3)**2 )
                                 
      F_int_22 = 1./mu0 * jacobi * ( B(2)**2 - 0.5 * Bnorm**2 )
      
      end function F_int_22
      
      function F_int_33( dat )
      real :: F_int_33
      class( dataCollectionBase ), intent(in), target :: dat
      real,dimension(3) :: B
      real :: Bnorm, jacobi
      
      call getB_field( dat, B, jacobi )
      
      Bnorm = sqrt( B(1)**2 + B(2)**2 + B(3)**2 )
                                 
      F_int_33 = 1./mu0 * jacobi * ( B(3)**2 - 0.5 * Bnorm**2 )
      
      end function F_int_33
      
       !::Returns the B-field given the specific dat model and jacobi-multiplier
        subroutine getB_field( dat, B, jacobi )
        class( dataCollectionBase ), intent(in), target :: dat
        real,dimension(3),intent(inout) :: B
        real,intent(inout) :: jacobi
        real,dimension(3) :: Hext
        
        if ( dat%coord_sys .eq.  coord_sys_carth ) then
            call getB_field_cart( dat, B )
            jacobi = 1
        elseif (dat%coord_sys .eq. coord_sys_cyl ) then
            call getB_field_cyl( dat, B, jacobi )     
        elseif ( dat%coord_sys .eq. coord_sys_cone ) then
            call getB_field_cone( dat, B, jacobi )
        endif
                
        end subroutine getB_field
        
        !::Get the field from the model solution in Cartesian coordinates
        subroutine getB_field_cart( dat, B )
        class( dataCollectionBase ), intent(in), target :: dat        
        real,dimension(3),intent(inout) :: B
        
        class( magStatModel ), pointer :: model
        real,dimension(3) :: Hsol, solPts,hext
        integer,dimension(3) :: inds
        
        !::Make the down cast to get the actual data model set
        call getDataModel( dat%model, model )
        
        !::Make sure the coordinate indices are mapped according to the normal vector
        if ( abs(dat%n_vec(1,1)) .eq. 1 ) then
            !yz-plane
            inds(1) = 2
            inds(2) = 3
            inds(3) = 1
        else if ( abs(dat%n_vec(1,2)) .eq. 1 ) then
            !xz-plane
            inds(1) = 1
            inds(2) = 3
            inds(3) = 2
        else if ( abs(dat%n_vec(1,3)) .eq. 1 ) then
            !xy-plane
            inds(1) = 1
            inds(2) = 2
            inds(3) = 3
        endif
        
        
        !::Reset the solution array
        Hsol(:) = 0
        
        solPts(inds(1)) = dat%x
        solPts(inds(2)) = dat%y
        solPts(inds(3)) = dat%z0;
        
        Hext(:) = 0.
        !::Get the solution from the MagStatVersion2 library
        call getFieldFromTiles( model%tiles, Hsol, solPts, model%n_tiles, 1 )                
	
        B = mu0 * Hsol
        end subroutine getB_field_cart
        
        
        !::Get the field from the model solution in cylindrical coordinates
        subroutine getB_field_cyl( dat, B, jacobi )
        class( dataCollectionBase ), intent(in), target :: dat
        class( magStatModel ), pointer :: model
        real,dimension(3),intent(inout) :: B
        real,intent(inout) :: jacobi
        
        real,dimension(3) :: Hsol, solPts 
        real :: x,y,z,r,theta,normMult
        
        !::Make the down cast to get the actual data model set
        call getDataModel( dat%model, model )
        
        !::Make sure the coordinate indices are mapped according to the normal vector
        if ( abs(dat%n_vec(1,1)) .eq. 1 ) then
            !theta-z plane
            !::dat%x represents theta, dat%y represents z and dat%z0 represents R, the constant radius of the
            !::cylinder surface
            r = dat%z0
            theta = dat%x            
            z = dat%y         
            !::part of the normal vector in the x-direction
            normMult = cos(theta)
        else if ( abs(dat%n_vec(1,2)) .eq. 1 ) then
            !theta-z plane, y-component
            !::dat%x represents theta, dat%y represents z and dat%z0 represents R, the constant radius of the
            !::cylinder surface
            r = dat%z0
            theta = dat%x            
            z = dat%y         
            !::part of the normal vector in the y-direction
            normMult = sin(theta)
        else if ( abs(dat%n_vec(1,3)) .eq. 1 ) then
            !r-theta plane
            !::dat%x represents r, dat%y represents theta and dat%z0 represents z
            r = dat%x
            theta = dat%y
            
            z = dat%z0
            normMult = 1
        endif
        
        jacobi = r * normMult
        
        x = r * cos( theta )
        y = r * sin( theta )
        
        !::Reset the solution array
        Hsol(:) = 0        
        
        solPts(1) = x
        solPts(2) = y
        solPts(3) = z
        
        !::Get the solution from the MagStatVersion2 library
        call getFieldFromTiles( model%tiles, Hsol, solPts, model%n_tiles, 1 )
        
        B = mu0 * Hsol
       
        end subroutine getB_field_cyl
       
        
        !::Get the field from the model solution on the surface of a cone
        subroutine getB_field_cone( dat, B, jacobi )
        class( dataCollectionBase ), intent(in), target :: dat
        class( magStatModel ), pointer :: model
        real,dimension(3),intent(inout) :: B
        real,intent(inout) :: jacobi        
        real,dimension(3) :: Hsol, solPts
        real :: x,y,z,r,theta,normMult
        
        !::Make the down cast to get the actual data model set
        call getDataModel( dat%model, model )
        
        !::Make sure the coordinate indices are mapped according to the normal vector
        if ( abs(dat%n_vec(1,1)) .eq. 1 ) then
            !theta-z plane
            !::dat%x represents theta, dat%y represents z and dat%z0 represents R, the constant radius of the
            !::cylinder surface            
            theta = dat%x            
            z = dat%y      
            
            r = (z - dat%cone_z0) * tan( dat%cone_angle )
            !::part of the normal vector in the x-direction
            normMult = cos(theta) * cos( dat%cone_angle )
        elseif ( abs(dat%n_vec(1,2)) .eq. 1 ) then
            !theta-z plane, y-component
            !::dat%x represents theta, dat%y represents z and dat%z0 represents R, the constant radius of the
            !::cylinder surface            
            theta = dat%x            
            z = dat%y         
            r = (z - dat%cone_z0) * tan( dat%cone_angle )
            !::part of the normal vector in the y-direction
            normMult = sin(theta) * cos( dat%cone_angle )
        elseif ( abs(dat%n_vec(1,3)) .eq. 1 ) then
            !r-theta plane
            !::dat%x represents r, dat%y represents theta and dat%z0 represents z
            r = dat%x
            theta = dat%y            
            z = dat%z0
            normMult = 1
        endif
        
        jacobi = r * normMult
        
        x = r * cos( theta )
        y = r * sin( theta )
        
        !::Reset the solution array
        Hsol(:) = 0
        
        solPts(1) = x
        solPts(2) = y
        solPts(3) = z
    
        !::Get the solution from the MagStatVersion2 library
        call getFieldFromTiles( model%tiles, Hsol, solPts, model%n_tiles, 1 )
	
        B = mu0 * Hsol
     
        end subroutine getB_field_cone
        
        
        !::Make the down cast
        subroutine getDataModel( modelIn, modelOut )
        class( dataCollectionModelBase ), intent(in), target :: modelIn
        class( magStatModel ), intent(out), pointer :: modelOut
        
            select type( modelIn )
            class is (magStatModel)
                modelOut => modelIn
                class default
                    write(*,*) 'Problem with the type cast!?'
            end select
        end subroutine getDataModel
    
end module MagForceMaxwellStressTensor
    