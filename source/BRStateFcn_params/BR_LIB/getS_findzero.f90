module BR_GETS_FINDZERO_CALL
    use BR_SMITH_CALL
    use BR_fzero
    use BR_DEBYE_CALL
    implicit none
    
    !::Private members for this module
    real,dimension(1) :: Hcurr, TcVal
    
    integer :: T_ind,H_ind
    real,dimension(:,:),allocatable :: Mout,Sout,deltaFout
    
contains

subroutine getS_fromZero( T, H, nT, nH, Tc_in, M, S, deltaF, dTad )
real,intent(in),dimension(nT) :: T
real,intent(in),dimension(nH) :: H
real,intent(inout),dimension(nT,nH) :: M, S, deltaF, dTad
real,intent(in),dimension(1) :: Tc_in
integer,intent(in) :: nT,nH
real :: tol
integer :: i,j

real,dimension(10) :: Mtmp,Stmp

TcVal = Tc_in

tol = 0.0001;

allocate( Mout(nT,nH),Sout(nT,nH) )
Mout(:,:) = 0
Sout(:,:) = 0

!::First get S as a function of T at H=H(1) (typically = 0)
call getM_BR( H(1), Mout(:,1), Sout(:,1), deltaF(:,1), T, Tc_in, nT, 1, 1 )
!::Get the background specific heat
call getDebye_S( T, Sout(:,1), nT, 1 )

!::Second, for each field find dTad, i.e. so that
!::S(T,H=0) == S(T+dTad,H=H)
do i=1,nT
    T_ind = i
    do j=2,nH
        Hcurr = H(j)
        H_ind = j
        dTad(i,j) = brent0 ( optS_dTad, T(i), T(i)+200, tol ) - T(i)
    enddo    
enddo

M = Mout
S = Sout

end subroutine getS_fromZero


function optS_dTad( T_in )
real,intent(in) :: T_in
real :: optS_dTad
real,dimension(1) :: T_inp

T_inp = T_in
!::Remember to reset S and M before calling!
Mout(T_ind,H_ind) = 0
Sout(T_ind,H_ind) = 0
deltaFout(T_ind,H_ind) = 0
!::Get the background contribution
call getDebye_S( T_inp, Sout(T_ind,H_ind), 1, 1 )

!::Find S at the given T and H
call getM_BR( Hcurr, Mout(T_ind,H_ind), Sout(T_ind,H_ind), deltaFout(T_ind,H_ind), T_inp, TcVal, 1, 1, 1 )


!::Subtract the zero-field S from this value, which is then what
!::the zero point is found on
optS_dTad = Sout(T_ind,H_ind) - Sout(T_ind,1)




end function optS_dTad
    
end module BR_GETS_FINDZERO_CALL
