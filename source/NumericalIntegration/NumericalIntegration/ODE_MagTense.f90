module ODE_MagTense_CALL

    implicit none
    
    contains
    
    !> Entry point for the integration of any set of ODE's using MagTense
    subroutine ODE_MagTense()
    
    end subroutine ODE_MagTense
    
    !> Specific implementation of solving a set of ODE's using the RK Suite
    subroutine ODE_RK_Suite()
    
    !>First call the setup
    call SETUP(NEQ,TSTART,YSTART,TEND,TOL,THRES,METHOD,TASK,ERRASS,
                       HSTART,WORK,LENWRK,MESAGE)

    !> Then loop over the desired time points and call the integrator
    call UT(F,TWANT,TGOT,YGOT,YPGOT,YMAX,WORK,UFLAG)
    
    end subroutine ODE_RK_Suite
    
end module ODE_MagTense_CALL
    