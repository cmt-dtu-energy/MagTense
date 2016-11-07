module BR_UTIL_CALL

    implicit none
    
contains

subroutine get_dTad_deltaS( T, S, dTad, deltaS, nT, nH )
real,intent(in),dimension(nT) :: T
real,intent(in),dimension(nT,nH) :: S
real,intent(inout),dimension(nT,nH) :: dTad,deltaS
integer, intent(in) :: nT, nH
integer :: i,j,ind
real :: lin

do i=1,nT
    do j=2,nH
        call locate( S(:,j), nT, S(i,1), ind )
        
        lin = ( S(i,1) - S(ind,j) ) / ( S(ind+1,j) - S(ind,j) )
        
        dTad(i,j) = T(ind+1) * lin + (1-lin) * T(ind) - T(i)
        
    enddo
enddo

do j=2,nH
    deltaS(:,j) = S(:,j) - S(:,1)
enddo


end subroutine get_dTad_deltaS
    
end module BR_UTIL_CALL
