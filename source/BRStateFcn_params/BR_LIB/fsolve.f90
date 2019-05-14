module BR_fzero
    implicit none
! This module is for finding a root of function(s). This contains 
! following items:
! 
! FUNCTIONS:
!   brent (fun, x0, tol, LowerBound, UpperBound, detail)
!   Find root of a function using Brent's method. This function tries to 
!   find x1 and x2 such that sign of f(x1) and f(x2) are not equal. Then it 
!   calls brent0(). 
! 
!   brent0 (fun, x1, x2, tol, detail)
!   Fund root of function fun given TWO intial values.
! 
!   fun1(x) ... fun3(x)
!   Test functions for maximization problem
! 
! SUBROUTINES:
!   test_functions()
!   It tests functions defined in this module. 




contains
        function brent0 (fun, x1, x2, tol, detail)
        ! Find zero in a fucntion fun using the Brent's method
        !
        ! Description of algorithm: 
        ! Find a root of function f(x) given intial bracketing interval [a,b]
        ! where f(a) and f(b) must have opposite signs. At a typical step we have
        ! three points a, b, and c such that f(b)f(c)<0, and a may coincide with
        ! c. The points a, b, and c change during the algorithm, and the root
        ! always lies in either [b,c] or [c, b]. The value b is the best
        ! approximation to the root and a is the previous value of b. 
        ! 
        ! The iteration uses following selection of algorithms 
        ! when bracket shrinks reasonablly fast,
        !  - Linear interporation if a == b
        !  - Quadratic interporation if a != b and the point is in the bracket.
        ! othrwise 
        ! - Bisection.
        ! 
        ! Inputs:
        !   fun: function to be solved
        !   x0, x1: Upper bound and lower bound for a function
        !   tol: error tolerance. 
        !   detail (optional): output result of iteration if detail is some value. 
        !  
        ! Date: Jan 30, 09
        ! Based on zeroin.f in netlib

        
        integer, parameter :: d = selected_real_kind(p=13,r=200)
        real :: brent0
        real, intent(IN) :: x1, x2, tol
        real, external :: fun
        integer, intent(IN), optional :: detail
        integer :: i, exitflag, disp
        real :: a, b, c, diff,e, fa, fb, fc, p, q, r, s, tol1,xm,tmp 
        real, parameter :: EPS = epsilon(a)
        integer, parameter :: imax = 100  ! maximum number of iteration
        ! values from Numerical Recipe
        
        exitflag = 0
        if (present(detail) .and. detail /= 0) then
                disp = 1
        else
                disp = 0
        end if
                
        ! intialize values
        a = x1
        b = x2
        c = x2
        fa = fun(a)
        fb = fun(b)
        fc=fb
                
        
        
        ! check sign
        if ( (fa .gt. 0. .and. fb .gt. 0. )  .or.  (fa .lt. 0. .and. fb .lt. 0. )) then            
                write(*,*)  'Error (brent.f90): Root must be bracked by two imputs'
                write(*,*) fa,fb
                write(*,*) x1,x2
                write(*,*) 'press any key to halt the program'
                read(*,*)
                stop
        end if
        
        if (disp == 1 ) then 
                write(*,*) 'Brents method to find a root of f(x)'
                write(*,*) ' '
                write(*,*) '  i           x          bracketsize            f(x)'
        end if
        
        ! main iteration
        do i = 1, imax
                ! rename c and adjust bounding interval if both a(=b) and c are same sign
                if ((fb > 0.  .and. fc > 0) .or. (fb <0. .and. fc < 0. ) ) then 
                        c = a
                        fc = fa
                        e = b-a
                        diff = e
                end if
                
                ! if c is better guess than b, use it. 
                if (abs(fc) < abs(fb) ) then 
                        a=b
                        b=c
                        c=a
                        fa=fb
                        fb=fc
                        fc= fa
                end if
                
                ! convergence check
                tol1=2.0* EPS * abs(b) + 0.5*tol
                xm = 0.5 * (c - b)
                if (abs(xm) < tol1 .or. fb == 0.0 )  then
                        exitflag = 1
                        exit
                end if
                
                if (disp == 1) then 
                        tmp = c-b
                        write(*,"('  ', 1I2, 3F16.6)") i, b, abs(b-c), fb
                end if
           
                ! try inverse quadratic interpolation
                if (abs(e) >= tol1 .and. abs(fa) > abs(fb) ) then 
                        s = fb/fa
                                if (abs(a - c) < EPS) then 
                                p = 2.0 *xm * s
                                q = 1.0  - s
                        else
                                q = fa/fc
                                r = fb/fc
                                p = s * (2.0 * xm * q * (q -r ) - (b - a) * (r - 1.0))
                                q = (q - 1.0 ) * (r - 1.0) * (s - 1.0) 
                        end if
                        
                        ! accept if q is not too small to stay in bound
                        if (p > 0.0) q = -q
                        p = abs(p)                
                        if (2.0 * p < min(3.0 * xm * q - abs(tol1* q), abs(e *q))) then 
                                e = d
                                diff = p / q
                        else   ! interpolation failed. use bisection
                                diff= xm 
                                e = d
                        end if
                else  ! quadratic interpolation bounds moves too slowly, use bisection
                        diff = xm
                        e = d
                end if
                
                ! update last bound
                a = b
                fa = fb
                
                ! move the best guess
                if (abs(d) > tol1) then 
                        b = b + diff
                else
                        b = b + sign(tol1, xm)
                end if
                
                ! evaluate new trial root
                fb = fun(b)        
        end do 
        
        ! case for non convergence
        if (exitflag /= 1 ) then 
                write(*,*) 'Error (brent.f90) :  convergence was not attained'
                write(*,*) 'Initial value:'
                write(*,"(4F10.5)" )   x1, x2, fun(x1), fun(x2)
                write(*,*) ' '
                write(*,*) 'final value:'
                !write(*,"('x = '  ,1F6.4, ':                f(x1) = ' ,  1F6.4  )" )  b,  fb  
        else if( disp == 1) then
                write(*,*) 'Brents method was converged.'
                write(*,*) ''
        end if
        brent0 = b
        return
        
        end function brent0
        
        function brent (fun, x0, tol, LowerBound, UpperBound, detail)
        ! Root finding using brent(). 
        ! find an initial guess of bracket and call brent()
        ! 
        ! Inputs
        ! fun: function to evaluate
        ! x0: Initial guess
        !
        ! Optional Inputs: 
        ! LowerBound, UpperBound : Lower and upper bound of the function
        ! detail : Output result of iteration if detail is there. 
        ! 
        ! Date: Jan 30, 09
        
        
        !integer, parameter :: d = selected_real_kind(p=13,r=200)
        real :: brent
        real, intent(IN) :: x0, tol
        real, external :: fun
        real, intent(IN), optional :: LowerBound, UpperBound
        integer, intent(IN), optional :: detail
        
        real :: a , b , olda, oldb, fa, fb
        real, parameter :: sqrt2 = sqrt(2.0)! change in dx
        integer, parameter :: maxiter = 40
        real :: dx  ! change in bracket
        integer :: iter, exitflag, disp
        real :: sgn
        
        a  = x0  ! lower bracket
        b =  x0 ! upper bracket
        olda = a  
        oldb = b 
        exitflag = 0  ! flag to see we found the bracket
        sgn  =fun(x0) ! sign of initial guess
        
        ! set disp variable
        if (present(detail) .and. detail /= 0) then
                disp = 1
        else
                disp = 0
        end if
                
        
        ! set initial change dx
        if (abs(x0)<0.00000002) then 
                dx = 1.0/50.0
        else
                dx = 1.0/50.0 * x0
        end if
        
        if (disp == 1) then 
                write(*,*) 'Search for initial guess for Brents method'
                write(*,*) 'find two points whose sign for f(x) is different '
                write(*,*) 'x1 searches downwards, x2 searches upwards with increasing increment'
                write(*,*) ' '
                write(*,*) '   i         x1         x2         f(x1)          f(x2)'
        end if
        
        
        ! main loop to extend a and b
        do iter = 1, maxiter
                fa = fun(a)
                fb = fun(b)
                
                if (disp == 1) write(*,"(1I4,4F14.7)") iter, a, b, fa, fb
                
                ! check if sign of functions changed or not
                if ( (sgn >= 0 ) .and.  (fa <= 0) ) then  ! sign of a changed 
                        ! use a and olda as bracket
                        b = olda
                        exitflag = 1
                        exit
                else if  ( (sgn <= 0 ) .and.  (fa >= 0  ) ) then ! sign of b changed
                        b = olda
                        exitflag = 1
                        exit
                else if  ( (sgn >= 0 ) .and.  (fb <= 0  ) ) then ! sign of a changed
                        a = oldb
                        exitflag = 1
                        exit
                else if  ( (sgn <= 0 ) .and.  (fb >= 0  ) ) then ! sign of a changed
                        a = oldb
                        exitflag = 1
                        exit
                end if
                
                ! update boundary
                olda = a 
                oldb = b
                a = a - dx
                b = b+ dx
                dx = dx * sqrt2
                
                ! boundary check
                if (present(LowerBound)) then
                        if (a < LowerBound ) a = LowerBound + tol
                end if
                if (present(UpperBound) ) then 
                        if (b > UpperBound ) b = UpperBound  -  tol
                end if
        end do
        
        
        if (exitflag /=  1 ) then 
                write(*,*) ' Error (brent2) : Proper initial value for Brents method could not be found'
                write(*,*) ' Change initial guess and try again. '
                write(*,*) ' You might want to try disp = 1 option too'
                write(*,*) 'i              x1            x2          fx1              fx2'
                write(*,"(1I4,4F12.7)") iter, a, b, fa, fb
                write(*,*) '  press any key to abort the program'
                read(*,*) 
                stop
        else if (disp == 1) then
                write(*,*) '  Initial guess was found.'
                write(*,*) ''
        end if
        
        if (present (detail)) then 
                brent = brent0(fun,a,b,tol,detail)
        else
                brent = brent0(fun,a,b,tol)
        end if
        
        end function brent
        
        
        
        
        subroutine TestFunctions()
        ! test various functions defined in this module
        implicit none
        !integer, parameter :: d = selected_real_kind(p=13,r=200)
        
        real :: out 
        
        ! test function 1 with a few more options
!         write(*,*) 'Function 1: answer =  0.09534 '
!         out = brent(fun =fun1, x0 = 4.0, tol = 0.0000001, & 
!                 LowerBound= 0.0, detail = 1)

        ! test function 2
!         write(*,*) 'Function 2: answer =  -3 '
!         out = brent(fun2, 3.0,0.00001, detail = 1)

        write(*,*) 'Function 3: answer =  8.10451, 10, 11.8955' 
        out = brent(fun3, 0.0, 0.000001, detail = 34)


        write(*,*) ' '
        
        write(*,*) 'press any key to continue'
        read(*,*) 
        
        end subroutine testfunctions
        
        function fun1  (x)
        ! test function for root finding
        ! f(x) = exp(x)  -  1 / (10 * x) ** 2
        ! answer = 0.095344
        implicit none
        !integer, parameter :: d = selected_real_kind(p=13,r=200)
        real :: fun1
        real , intent(IN) :: x
        
        fun1 =  exp(x)  -  1 / (10 * x) ** 2
        end function fun1

        function fun2  (x)
        ! test function for root finding
        ! f(x) = (x+3) *  (x - 1) ** 2
        ! answer:  -3, 1
        implicit none
        !integer, parameter :: d = selected_real_kind(p=13,r=200)
        real :: fun2
        real , intent(IN) :: x
        
        fun2 =  (x+3) *  (x - 1) ** 2
        end function fun2

        function fun3  (x)
        ! test function for root finding
        ! f(x) = sin(x-10)  - 0.5 * (x - 10) 
        ! answer = 10
        implicit none
        !integer, parameter :: d = selected_real_kind(p=13,r=200)
        real :: fun3
        real , intent(IN) :: x
        
        fun3 =  sin(x-10)  - 0.5 * (x - 10)
        end function fun3

        
end module BR_fzero