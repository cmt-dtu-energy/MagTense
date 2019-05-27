
    module NumInt
    
    use QUADPACK
    use integrationDataTypes

    implicit none
    
    contains

    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Surface_integral_cone calculates the surface integral on a cone
    !! @param surf Defines the min and max vals of the x, y and z integration limits
    !! @params dat_arr An array of 18 elements each being a pointer to a data structure used by the integration routine
    !! @params handleError A pointer to a specific subroutine that deals with any errors from the integration
    !! @params vOut A 3-element vector that returns the solution of the three vector components
    !! @params ier A 2-element array with error-codes from the x- and y-integrations
    !! @params neval A 2-element arrant with the maximum number of function evaluations from the x- and y-integrations
    subroutine surface_integral_cone( surf, dat_arr, handleError, vOut, ier, neval )
        type( surf_carth ), intent(in) :: surf
        type( dat_ptr ), intent(inout), dimension(15) :: dat_arr
        procedure (error_handler), intent(in), pointer :: handleError
    
        real,dimension(3),intent(inout) :: vOut
    
        !! Two surfaces on top and bottom, respectively, each with one integral for each force component
        !! plus one surface, the cylindrical surface, each with two integrals for each force component
        integer,parameter :: n = 5
    
        integer,intent(inout),dimension(2) :: ier,neval
        real,dimension(n,3) :: n_vec
        real,dimension(n,2) :: xp,yp    !< defines the plane of integration
        real,dimension(n) :: zp,vecSign !< defines the position of the plane on the third axis
        integer :: i,j,ind
        real :: res
    
        vOut(:) = 0.    !< Reset input array
    
        neval = 0
      
        !! Setup the three planes to be integrated over
        !! theta-z plane, first integral (x-part of the normal vector)
        xp(1,:) = surf%theta
        yp(1,:) = surf%z
        zp(1) = surf%r(2)
        vecSign(1) = 1
        n_vec(1,:) = (/ 1, 0, 0 /)
    
        !! theta-z plane, second integral (y-part of the normal vector)
        xp(2,:) = surf%theta
        yp(2,:) = surf%z
        zp(2) = surf%r(2)
        vecSign(2) = 1
        n_vec(2,:) = (/ 0, 1, 0 /)
    
        !! theta-z plane, third integral (z-part of the normal vector)
        xp(3,:) = surf%theta
        yp(3,:) = surf%z
        zp(3) = surf%r(2)
        vecSign(3) = -1
        n_vec(3,:) = (/ 0, 0, -1 /)
      
        !! r-theta- plane, the smaller of the two conical circles
        xp(4,1) = surf%r(1)
        xp(4,2) = surf%r(3)
        yp(4,:) = surf%theta
        zp(4) = surf%z(1)
        vecSign(4) = -1
        n_vec(4,:) = (/ 0, 0, -1 /)
      
        !! r-theta+ plane, the bigger of the two conical circles
        xp(5,1) = surf%r(1)
        xp(5,2) = surf%r(2)
        yp(5,:) = surf%theta
        zp(5) = surf%z(2)
        vecSign(5) = 1
        n_vec(5,:) = (/ 0, 0, 1 /)
      
        call surf_int( dat_arr, xp, yp, zp, handleError, surf%retVec, n, vecSign, n_vec, vOut, ier, neval )
    
    end subroutine surface_integral_cone
    
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Surface_integral_cyl calculates the surface integral for a cylinder
    !! @param surf Defines the min and max vals of the x, y and z integration limits
    !! @param dat_arr An array of 18 elements each being a pointer to a data structure used by the integration routine    
    !! @param handleError A pointer to a specific subroutine that deals with any errors from the integration
    !! @param vOut A 3-element vector that returns the solution of the three vector components
    !! @param ier A 2-element array with error-codes from the x- and y-integrations
    !! @param neval A 2-element arrant with the maximum number of function evaluations from the x- and y-integrations
    subroutine surface_integral_cyl( surf, dat_arr, handleError, vOut, ier, neval )
        type( surf_carth ), intent(in) :: surf
        type( dat_ptr ), intent(inout), dimension(12) :: dat_arr
        procedure (error_handler), intent(in), pointer :: handleError
    
        real,dimension(3),intent(inout) :: vOut
    
        !! two surfaces on top and bottom, respectively, each with one integral for each force component
        !! plus one surface, the cylindrical surface, each with two integrals for each force component
        integer,parameter :: n = 4
    
        integer,intent(inout),dimension(2) :: ier,neval
        real,dimension(4,3) :: n_vec
        real,dimension(4,2) :: xp,yp    !! defines the plane of integration
        real,dimension(4) :: zp,vecSign !! defines the position of the plane on the third axis
        integer :: i,j,ind
        real :: res
    
        vOut(:) = 0.    !< Reset input array
    
        neval = 0
      
        !! Setup the three planes to be integrated over
        !! theta-z plane, first integral (x-part of the normal vector)
        xp(1,:) = surf%theta
        yp(1,:) = surf%z
        zp(1) = surf%r(2)
        vecSign(1) = 1
        n_vec(1,:) = (/ 1, 0, 0 /)
    
        !! theta-z plane, second integral (y-part of the normal vector)
        xp(2,:) = surf%theta
        yp(2,:) = surf%z
        zp(2) = surf%r(2)
        vecSign(2) = 1
        n_vec(2,:) = (/ 0, 1, 0 /)
      
        !! r-theta- plane
        xp(3,:) = surf%r(1:2)
        yp(3,:) = surf%theta
        zp(3) = surf%z(1)
        vecSign(3) = -1
        n_vec(3,:) = (/ 0, 0, -1 /)
      
        !! r-theta+ plane
        xp(4,:) = surf%r(1:2)
        yp(4,:) = surf%theta
        zp(4) = surf%z(2)
        vecSign(4) = 1
        n_vec(4,:) = (/ 0, 0, 1 /)
      
        call surf_int( dat_arr, xp, yp, zp, handleError, surf%retVec, n, vecSign, n_vec, vOut, ier, neval )
    
    end subroutine surface_integral_cyl
    
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Subroutine that does the integral not depending on the coordinate base
    !! @param dat_arr has n components (3 vector components and then for each of these a number of surfaces integrated over)
    !! @param xp,yp and zp contains the limits to the surfaces integrated over
    !! @param handleError A subroutine pointer to the specific errorhandler
    !! @param retVector Determines which components to be calculated (0 for yes, otherwise don't)
    !! @param n The total number of surface integrals to be performed (n_vec * n_surf)
    !! @param vOut The output vector, i.e. the result of the calculation
    !! @param ier Contains any error values
    !! @param neval Contains the number of function evaluations
    subroutine surf_int( dat_arr, xp, yp, zp, handleError, retVector, n, vecSign, n_vec, vOut, ier, neval )
        type( dat_ptr ), intent(inout), dimension(n*3) :: dat_arr
        procedure (error_handler), intent(in), pointer :: handleError
        integer,intent(in),dimension(3) :: retVector
        integer,intent(in) :: n
        real,dimension(n),intent(in) :: vecSign
        real,dimension(n,3),intent(in) :: n_vec
        real,dimension(3),intent(inout) :: vOut
        integer,intent(inout),dimension(2) :: ier,neval
        real,dimension(n,2) :: xp,yp
        real,dimension(n) :: zp
        integer :: i, j, ind
        real :: res
      
        vOut(:) = 0.
        
        !! j loops over the three vector components     
        do j=1,3                        
            do i=1,n
                if ( retVector(j) .eq. 0 ) then
                    ind = n*(j-1) + i
                
                    dat_arr(ind)%dat%x1 = xp(i,1)
                    dat_arr(ind)%dat%x2 = xp(i,2)
                    dat_arr(ind)%dat%y1 = yp(i,1)
                    dat_arr(ind)%dat%y2 = yp(i,2)
                
                    dat_arr(ind)%dat%z0 = zp(i)
                
                    dat_arr(ind)%dat%n_vec(1,:) = n_vec(i,:)

                    call integral2( dat_arr(ind)%dat, res )
                
                    call handleError(dat_arr(ind)%dat, dat_arr(ind)%dat%abserr_tot) 
    
                    vOut(j) = vOut(j) + vecSign(i) * res           
    
                endif                
            enddo                          
        enddo   
    
        neval = 0
        ier = 0
        do i=1,n*3
            neval(1) = neval(1) + dat_arr(i)%dat%neval_x
            neval(2) = neval(2) + dat_arr(i)%dat%neval_y
        
            if ( dat_arr(i)%dat%ier_x .ne. 0 ) then
                ier(1) = dat_arr(i)%dat%ier_x
            endif
        
            if ( dat_arr(i)%dat%ier_y .ne. 0 ) then
                ier(2) = dat_arr(i)%dat%ier_y
            endif
      
        enddo

    end subroutine surf_int
    
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! @param surf defines the min and max vals of the x, y and z integration limits
    !! @param dat_arr is an array of 18 elements each being a pointer to a data structure used by the integration routine
    !! @param handleError is a pointer to a specific subroutine that deals with any errors from the integration
    !! @param vOut is a 3-element vector that returns the solution of the three vector components
    !! @param ier is a 2-element array with error-codes from the x- and y-integrations
    !! @param neval is a 2-element arrant with the maximum number of function evaluations from the x- and y-integrations
    subroutine surface_integral_carth( surf, dat_arr, handleError, vOut, ier, neval )
        type( surf_carth ), intent(in) :: surf
        type( dat_ptr ), intent(inout), dimension(18) :: dat_arr
        procedure (error_handler), intent(in), pointer :: handleError
    
        real,dimension(3),intent(inout) :: vOut
    
        integer,intent(inout),dimension(2) :: ier,neval
        real,dimension(6,3) :: n_vec
        real,dimension(6,2) :: xp,yp    !< defines the plane of integration
        real,dimension(6) :: zp,vecSign !< defines the position of the plane on the third axis
        integer :: i,j,ind,n
        real :: res
        
        n = 6   !< number of surfaces
    
        vOut(:) = 0.    !< Reset input array
    
        neval = 0
      
        !! Setup the six planes to be integrated over
        !! yz- plane
        xp(1,:) = surf%y
        yp(1,:) = surf%z
        zp(1) = surf%x(1)
        vecSign(1) = -1
        n_vec(1,:) = (/ -1, 0, 0 /)
      
        !! yz+ plane
        xp(2,:) = surf%y
        yp(2,:) = surf%z
        zp(2) = surf%x(2)
        vecSign(2) = 1
        n_vec(2,:) = (/ 1, 0, 0 /)
      
        !! xz- plane
        xp(3,:) = surf%x
        yp(3,:) = surf%z
        zp(3) = surf%y(1)
        vecSign(3) = -1
        n_vec(3,:) = (/ 0, -1, 0 /)
      
        !! xz+ plane
        xp(4,:) = surf%x      
        yp(4,:) = surf%z
        zp(4) = surf%y(2)
        vecSign(4) = 1
        n_vec(4,:) = (/ 0, 1, 0 /)
      
        !! xy- plane
        xp(5,:) = surf%x
        yp(5,:) = surf%y
        zp(5) = surf%z(1)
        vecSign(5) = -1
        n_vec(5,:) = (/ 0, 0, -1 /)
      
        !! xy+ plane
        xp(6,:) = surf%x
        yp(6,:) = surf%y
        zp(6) = surf%z(2)
        vecSign(6) = 1
        n_vec(6,:) = (/ 0, 0, 1 /)
    
        call surf_int( dat_arr, xp, yp, zp, handleError, surf%retVec, n, vecSign, n_vec, vOut, ier, neval )
    
    end subroutine surface_integral_carth
    
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Numerical double integral
    !! Performs the numerical double integral of user supplied function fun on the rectangular interval [xl,xh], [yl,yh]
    !! Returns the result in res (kind = real)
    !! The user supplied function fun should take two arguments, x and y: fun=fun(x,y)
    subroutine integral2( dat, res )
        real,intent(inout) :: res
        class(dataCollectionBase),intent(inout), target :: dat
        procedure (f_int_dat), pointer :: f_ptr => null ()
        !procedure (f_int_dat_vec), pointer :: f_ptr => null ()

        res = 0
        !f_ptr => h_vec
        f_ptr => h

        dat%progCallbackCnt = 0
        call qags_x ( f_ptr, dat, dat%x1, dat%x2, dat%epsabs, dat%epsrel, res, dat%abserr_x, dat%neval_x, dat%ier_x )
        !res = qromb_mod( f_ptr, dat, dat%x1, dat%x2 )

    end subroutine integral2

    
    !!@todo These function from NR should be in the NR project?

    !::Adopted from NR
    function h(xx, dat)
    real :: h
    real,intent(in) :: xx
    real :: ss, res
    class(dataCollectionBase), target :: dat
    procedure (f_int_dat_vec), pointer :: f_ptr => null ()
    real :: tmp

        f_ptr => f_vec
        dat%x = xx

    
        call qags_y_vec( f_ptr, dat, dat%y1, dat%y2, dat%epsabs, dat%epsrel, ss, dat%abserr_y, dat%neval_y, dat%ier_y )    
        h = ss
        
        if ( mod( dat%progCallbackCnt, 100 ) .eq. 0 ) then
            !::Progress callback
            tmp = dat%progCallback( dat )
        endif
    
    
        dat%progCallbackCnt = dat%progCallbackCnt + 1
    
    return
    end function h

    !function h_vec(xx, dat, n)
    !integer,intent(in) :: n
    !real,dimension(n) :: h_vec
    !real,intent(in),dimension(n) :: xx
    !real :: ss, res
    !class(dataCollectionBase), target :: dat
    !procedure (f_int_dat_vec), pointer :: f_ptr => null ()
    !real :: tmp
    !integer :: i
    !
    !    do i=1,n
    !        f_ptr => f_vec
    !        dat%x = xx(i)
    !
    !        h_vec(i) = qromb_mod( f_ptr, dat, dat%y1, dat%y2 )
    !    
    !        if ( mod( dat%progCallbackCnt, 100 ) .eq. 0 ) then
    !            !::Progress callback
    !            tmp = dat%progCallback( dat )
    !        endif
    !        
    !        dat%progCallbackCnt = dat%progCallbackCnt + 1
    !    enddo
    !    
    !    
    !return
    !end function h_vec


    !::adopted from NR
    function f( yy, dat )
    real :: f
    real,intent(in) :: yy
    class(dataCollectionBase),intent(inout), target :: dat


    dat%y = yy

    f = dat%f_ptr(dat)

    return
    end function f

    function f_vec( yy, dat, n )
    integer,intent(in) :: n
    real,intent(in),dimension(n) :: yy
    real,dimension(n) :: f_vec
    class(dataCollectionBase),intent(inout), target :: dat

    call dat%f_ptr_vec(yy, dat, n, f_vec )


    return
    end function f_vec


end module NumInt
