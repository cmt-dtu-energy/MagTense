module BR_parameters_call

    implicit none
    
    real,parameter :: recl=2,pi=3.14159265358979
    
    real :: Tmin,Tmax,dT
    
    real :: Hmin,Hmax,dH
    
    real :: eta
          
    real :: kappa_,eta_,J_,g_,Ns_,rho_
    real :: Tdeb_,Na_c_,mol_,gammaSommerfeld_
    
    real :: T0,T0_std
    
    real,dimension(:),allocatable :: T0_arr
    
    integer :: nT0
    
    integer :: calcS_isoT_isoH
    
    
    character(len=1000) :: file_input,file_out,file_Tc_distr
    
    real,dimension(:),allocatable :: T,H
    real,dimension(:,:,:),allocatable :: S,M,deltaF
    real,dimension(:,:),allocatable :: dTad,deltaS
    
    integer :: nT,nH
    
contains
    
    subroutine setupParameters()
    integer :: i
    real,dimension(:),allocatable :: Z0
    
    nT = ceiling( (Tmax-Tmin)/dT )+1
    
    nH = ceiling( (Hmax-Hmin)/dH )+1
    write(*,*) 'nH',nH
    
    allocate( S(nT,nH,2),M(nT,nH,2),dTad(nT,nH),deltaS(nT,nH), deltaF(nT,nH,2) )
    
    allocate( T(nT), H(nH) )
    allocate( T0_arr(nT0) )
    
    S(:,:,:) = 0
    M(:,:,:) = 0
    deltaF(:,:,:) = 0
    
    
    dTad(:,:) = 0
    deltaS(:,:) = 0
    
        
    T = 0
    H = 0    
    T0_arr = 0
    if ( nT0 .eq. 1 ) then
        T0_arr(1) = T0
    else            
        
        call getNormalRandomNumbers( nT0, Z0 )
    
        T0_arr = Z0 * T0_std + T0
    
    endif
    !::Setup T and H arrays
    T(1) = Tmin
    do i=2,nT
        T(i) = T(i-1) + dT
    enddo
    write(*,*) T
    H(1) = Hmin
    do i=2,nH
        H(i) = H(i-1) + dH
    enddo
    write(*,*) H
    
    write(*,*) nT,nH
    
    end subroutine setupParameters
    
    subroutine getNormalRandomNumbers( nT0, Z0 )
    integer,intent(in) :: nT0
    real,intent(inout),dimension(:),allocatable :: Z0
    
    real,dimension(:),allocatable :: U1, U2
    
    allocate( U1(nT0),U2(nT0),Z0(nT0) )
        call random_seed()
        call RANDOM_NUMBER(U1)
        call RANDOM_NUMBER(U2)
    
        Z0 = sqrt( -2*log(U1) ) * cos( 2*pi * U2)
    
    deallocate( U1, U2 )
    end subroutine getNormalRandomNumbers
    
end module BR_parameters_call

    