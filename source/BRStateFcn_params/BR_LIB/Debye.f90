module BR_DEBYE_CALL
    
implicit none    
    
    real,private :: Tdeb, Na_c, mol,gammaSommerfeld
    
    real,private,parameter :: kB=1.38e-23,Na=6.02e23        
contains
 
    
subroutine getDebye_S( T, S, nT,nH )
real,intent(in),dimension(nT) :: T
real,intent(inout),dimension(nT,nH) :: S
integer,intent(in) :: nT,nH
real :: integS
integer :: i

    do i=1,nT
        call qromb( debyeS, 0.0001,Tdeb / T(i), integS )
	    
        S(i,:) = S(i,:) + kB * Na * Na_c /mol * ( -3 * log( 1 - exp(-Tdeb/T(i) )) + 12 * (T(i) / Tdeb )**3 * integS ) + T(i) * gammaSommerfeld
        
    enddo
    
end subroutine getDebye_S


function debyeS( x )
real,intent(in) :: x
real :: debyeS

    debyeS = x**3 / ( exp(x) - 1 )

end function debyeS

subroutine initializeDebye( Tdeb_, Na_c_, mol_, gammaSommerfeld_)
real,intent(in) :: Tdeb_, Na_c_, mol_,gammaSommerfeld_

Tdeb = Tdeb_
Na_c = Na_c_
mol = mol_
gammaSommerfeld = gammaSommerfeld_


end subroutine initializeDebye

end module BR_DEBYE_CALL
