module BR_SMITH_CALL
    use BR_fzero
implicit none
    !::Particular BR parameters that are only available within this scope
    !::(which are constant parameters)
    real,private,parameter :: pi=3.14159265358979,muB=9.274e-24,kB=1.38e-23,mu0=4*pi*1e-7
    !::(which are constant)
    real,private :: kappa,eta,J,g,Ns,rho
    !::(which are variable and not thread safe when finding M)
    real,private :: beta, lambda0, T_cst, H_cst, p_cst, A,B,gamma,NV
    integer,private :: signBeta
    
    
contains
    
!::Fill out the particular parameters relevant for the BR model
subroutine initializeBRParameters( kappa_, eta_, J_, g_, Ns_, rho_, signBeta_  )
real,intent(in) :: kappa_,eta_,J_,g_,Ns_,rho_
integer,intent(in) :: signBeta_

      

!    open (11, file='debug.txt',status='old', access='sequential',form='formatted',action='write',position='append' )
    kappa = kappa_
    eta = eta_
    J = J_
    g = g_
    Ns = Ns_
    rho = rho_
    signBeta = signBeta_

    !::Number of spins per volume
    NV = Ns * rho
    gamma = g * muB
    !::Constants for the Br function
    A = ( 2 * J + 1 ) / ( 2 * J )
    B = 1 / ( 2 * J )

!    write(11,*) 'BR parameters initialized'
!    write(11,*) kappa,eta,J,g,Ns,rho,NV,gamma,A,B,signBeta
!    close(11)
end subroutine initializeBRParameters
    
!::Bean-Rodbell modified MFT that includes a distribution in Tc
!::This version assumes H, T and Tc_arr to be of the same size and thus
!::representing a batch of individual particles that are not directly interacting
subroutine getM_BR_batch( T, H, p, Tc_arr, M, S, deltaF, n )
integer,intent(in) :: n

integer,parameter :: maxCnt=2
real,intent(in),dimension(n) :: T,H,Tc_arr
real,intent(in) :: p

real,intent(inout),dimension(n,maxCnt) :: M, S, deltaF
integer :: ii,uu

real,dimension(maxCnt) :: Mval
real :: startInt,deltaM, nextPoint, endInt,tol,Fmin,Fcurr,time_start,time_end
real :: testLow,testHigh,Mtmp,Fmax
integer :: cnt



call cpu_time( time_start )
deltaM = 1e4
endInt = 1e7
tol = 1e-4

p_cst = p

!::Reset M
M(:,:) = 0

do ii=1,n
        !::Main loop for doing the BR solution            
        !::The beta-parameter from BR                
        beta = signBeta * sqrt( 1./40. * eta * 1. / ( NV * kB * Tc_arr(ii) * kappa ) * ( ( 2 * J + 1 )**4 - 1 )/( J * ( J + 1 ) )**2 )
        
        
        !::The lambda0 parameter    
        lambda0 = Tc_arr(ii) * 3. * kB / ( J * ( J + 1 ) * gamma**2 ) * 1/NV
        
        !::Temperature and field
        
        T_cst = T(ii)
        H_cst = H(ii)        
        !::Reset the result array
        Mval = 0
        !::Starting value
        !::Has to be greater than 0 since for H = 0 and T>Tc
        !::M will be zero and that will be found below
        !::Otherwise, when H = 0 and T<Tc then dFdM starts negative and we won't find the zero point
        startInt = 100
        nextPoint = startInt + deltaM
        cnt = 1
            
        !Make sure to include the first M as zero to ensure that possible minimum
        if ( dFdM( deltaM ) .gt. 0 .AND. H_cst .eq. 0 ) then
            !Then M is zero but S_mag is not            
            cnt = cnt + 1
            
        endif
            
        do                            
            if ( nextPoint .lt. endInt .AND. cnt .le. maxCnt ) then
                testLow = dFdM(startInt)
                testHigh = dFdM(nextPoint)
                    
                if ( testLow .gt. 0 .AND. testHigh .lt. 0 .or. testLow .lt. 0 .AND. testHigh .gt. 0 ) then                    
                    Mtmp = brent0 ( dFdM, startInt, nextPoint, tol )
        
                    if ( testLow .lt. 0 ) then
                        Mval(cnt) = Mtmp
                            
                        cnt = cnt + 1
                    else
                        !::The global maximum
                        Fmax = F_M( Mtmp )
                    endif
                        
                        
                    startInt = Mtmp*1.01
                    do
                        if ( dFdM(startInt) .ne. 0 ) then
                            exit
                        else
                            startInt = startInt * 1.01
                        endif                            
                    enddo
                        
                    nextPoint = startInt + deltaM
                        
                else
                    nextPoint = nextPoint + deltaM
                endif
            else
                exit
            endif
        enddo
        
        !Numerical problem patch
        if ( Mval(2) .eq. 0 .AND. H_cst .gt. 0 ) then
            Mval(2) = Mval(1)
        elseif ( Mval(1) .eq. 0 .AND. H_cst .gt. 0 ) then
            Mval(1) = Mval(2)                
        endif
         
        !::The heating curve should be the first element in the array, i.e.
        !::the maximum value of M represents this
        if ( Mval(1) .lt. Mval(2) ) then
            M(ii,1) = M(ii,1) + Mval(2)
            S(ii,1) = S(ii,1) + Smag_( Mval(2) )
            
            M(ii,2) = M(ii,2) + Mval(1)
            S(ii,2) = S(ii,2) + Smag_( Mval(1) )
            
        else
            M(ii,1) = M(ii,1) + Mval(1)
            S(ii,1) = S(ii,1) + Smag_( Mval(1) )
            
            M(ii,2) = M(ii,2) + Mval(2)
            S(ii,2) = S(ii,2) + Smag_( Mval(2) )
            
        endif
                                                
        deltaF(ii,1) = abs( F_M(Mval(1)) - Fmax ) / ( T_cst * kB * NV )
        deltaF(ii,2) = abs( F_M(Mval(2)) - Fmax ) / ( T_cst * kB * NV )
            
        !::The second entry into the array should be the global minimum, 
        !::which is essentially the cooling curve
        !if ( deltaF(ii,1) .ge. deltaF(ii,2) ) then
        !    M(ii,2) = M(ii,2) + Mval(1)
        !    S(ii,2) = S(ii,2) + Smag_( Mval(1) )
        !else
        !    M(ii,2) = M(ii,2) + Mval(2)
        !    S(ii,2) = S(ii,2) + Smag_( Mval(2) )
        !endif
enddo


call cpu_time( time_end )

!write(*,*) 'time for BR',time_end-time_start
end subroutine getM_BR_batch


!::Magnetic part of the total entropy
function Smag_( M )
real,intent(in) :: M
real :: Smag_
real :: x,Bri
    x = getX( T_cst, H_cst, p_cst, M, lambda0, beta )
    if ( M .gt. 0 ) then            
        Bri = getBr( x, A, B )
        
        Smag_ = kB * NV / rho * ( log( sinh( A * x ) / sinh( B * x ) )  - x * Bri );
    else
        Smag_ = kB * NV / rho * log( 2*J+1 )        
    endif
    
end function Smag_

!::The free energy
function F_M( M )
real,intent(in) :: M
real :: F_M
real :: x
    if ( M .eq. 0 .and. H_cst .eq. 0 ) then
        F_M = -kB * T_cst * NV * log( 2 * J + 1 )
    else
        
        x = getX( T_cst, H_cst, p_cst, M, lambda0, beta )
        
        !17-8-2016, kaki: updated the F_M term to include pressure
        !F_M = -kB * T_cst * NV * log( sinh( A * x ) / sinh( B * x ) ) + 0.5 * lambda0 * M**2 + 3./8. * lambda0**2 * beta**2 * kappa * M**4
        
        F_M = -kB * T_cst * NV * log( sinh( A * x ) / sinh( B * x ) ) + 0.5 * lambda0 * M**2 *( 1 - beta * p_cst * kappa ) + 3./8. * lambda0**2 * beta**2 * kappa * M**4 - 0.5 * kappa * p_cst**2
    endif
    
end function F_M
!::The derivative of the free energy with respect to magnetization
function dFdM( M )
real,intent(in) :: M
real :: dFdM
real :: Bri,x
!    write(*,*) H_cst
    x = getX( T_cst, H_cst, p_cst, M, lambda0, beta )
    
    Bri = getBr( x, A, B )
    
    !dFdM = -NV * gamma * J * lambda0 * ( 1 + 3./2. * lambda0 * beta**2 * kappa * M**2 ) * Bri + lambda0 * M + 3./2. * lambda0**2 * beta**2 * kappa * M**3
    
    !17-8-2016, kaki: Updated the dFdM term to include pressure
    dFdM = -NV * gamma * J * lambda0 * ( 1 + 3./2. * lambda0 * beta**2 * kappa * M**2 - beta * p_cst * kappa ) * Bri + 3./2. * lambda0**2 * beta**2 * kappa * M**3 + lambda0 * M * (1 - p_cst * kappa * beta )
    
    
end function dFdM        

!::Brillouin function
function getBr( x, A, B )
real :: getBR
real,intent(in) :: x, A, B
    
    getBr = A * 1./tanh( A * x ) - B * 1./tanh( B * x )

end function getBr


!::The parameter x from BR-Smith formulation
function getX( T, H, p, M, lambda0, beta )
real,intent(in) :: T, H, p, M, lambda0, beta
real :: getX

    !17-8-2016, kaki: updated the x-term to include pressure
    gamma = g * muB    
    !getX = gamma * J / ( kB * T ) * ( H + lambda0 * M + 0.5 * lambda0**2 * beta**2 * kappa * M**3 )
    
    getX = gamma * J / ( kB * T ) * ( H + lambda0 * M + 0.5 * lambda0**2 * beta**2 * kappa * M**3 - lambda0 * M * beta * p * kappa )
    
end function getX
end module BR_SMITH_CALL
