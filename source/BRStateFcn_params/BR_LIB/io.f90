module BR_IO_CALL
use BR_parameters_call
    implicit none
    
contains

subroutine writeData()

open (11, file=file_out,	&
			   status='unknown', form='unformatted',	&
			   access='stream')

    write(11) nT
    write(11) nH
    write(11) nT0
    write(11) T
    write(11) H
    write(11) M
    write(11) S
    write(11) dTad
    write(11) deltaS
    write(11) deltaF

close(11)


end subroutine writeData
    
subroutine readInput( file )
character(len=1000),intent(in) :: file

    
    open(11,file=file,status='old',access='sequential',form='formatted',action='read')
    
    read(11,*) file_input
    read(11,*) file_out    
    
    close(11)
    
    open(11,file=file_input,status='old',access='sequential',form='formatted',action='read')
    
    read(11,*) kappa_, J_, g_
    write(*,*) 'kappa_, J_, g_'
    write(*,*) kappa_, J_, g_
    
    read(11,*) Ns_, rho_,Tdeb_
    write(*,*) 'Ns_, rho_,Tdeb_'
    write(*,*) Ns_, rho_,Tdeb_
    
    read(11,*) Na_c_, mol_,gammaSommerfeld_
    write(*,*) 'Na_c_, mol_,gammaSommerfeld_'
    write(*,*) Na_c_, mol_,gammaSommerfeld_
    
    read(11,*) Tmin,Tmax,dT
    write(*,*) 'Tmin,Tmax,dT'
    write(*,*) Tmin,Tmax,dT
    
    read(11,*) Hmin,Hmax,dH
    write(*,*) 'Hmin,Hmax,dH'
    write(*,*) Hmin,Hmax,dH
    
    read(11,*) eta
    write(*,*) 'eta'
    write(*,*) eta
        
    read(11,*) nT0,T0,T0_std
    write(*,*) 'nT0,T0,T0std'
    write(*,*) nT0,T0,T0_std
    
    read(11,*) calcS_isoT_isoH
    write(*,*) 'calcS_isoT_isoH (0 = yes, 1 = use findzero to get dTad)'
    write(*,*) calcS_isoT_isoH
    
    close(11)
    

    
    
end subroutine readInput

end module BR_IO_CALL
