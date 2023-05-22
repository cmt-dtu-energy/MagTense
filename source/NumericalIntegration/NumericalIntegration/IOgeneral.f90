#if USE_MATLAB
#include "fintrf.h"
#endif

module IO_GENERAL

    implicit none
    contains
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !---------------------------------------------------------------- Message back to the GUI ---------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    !>----------------------------------------
    !> Rasmus Bjørk, rabj@dtu.dk, August 2022
    !> Callback to display message to the GUI launching the simulation
    !> @params[in] mess The message to display to the GUI
    !>----------------------------------------    
    
    subroutine displayGUIMessage( mess )
        character(*),intent(in) :: mess
        logical :: ex

        
        !! test if we are inside Matlab or called from a stand-alone. The simple test
        !! is whether io.txt exists in the current path (false for Matlab, true for stand-alone)
        !! nothing bad should happen if in fact we are called from ML but the file somehow
        !! exists - the written output will just not be shown to the user

        inquire( file='io.txt', EXIST=ex )

        if ( ex .eq. .true. ) then
            write(*,*) mess
        endif
        
#if USE_MATLAB    
        call displayMatlabMessage( mess )
#endif 
    
    end subroutine displayGUIMessage

    
    subroutine displayGUIProgressMessage( mess, prog )
        character(*),intent(in) :: mess
        integer,intent(in) :: prog
        character*(100) :: prog_str
        
        if (prog .lt. 0) then
            call displayGUIMessage( trim(mess) )
        else
            write(prog_str,'(A50, I4.0)') mess, prog           
            call displayGUIMessage( trim(prog_str) )
        endif

    end subroutine displayGUIProgressMessage
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !---------------------------------------------------------------- Message Matlab interface ---------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#if USE_MATLAB  

    !! Routine for displaying progress within Matlab
    subroutine displayMatlabMessage( mess )
        character(*),intent(in) :: mess

        integer :: mexCallMATLAB, nlhs_cb, nrhs_cb, tmp
        mwPointer plhs_cb(1), prhs_cb(1),mxCreateString
        character*(4) functionName_cb
        
        functionName_cb = "disp"
        nlhs_cb = 0
        nrhs_cb = 1
        
        prhs_cb(1) = mxCreateString(mess)
        
        tmp = mexCallMATLAB(nlhs_cb, plhs_cb, nrhs_cb, prhs_cb, "disp")
        

    end subroutine displayMatlabMessage
   
#endif
    
end module IO_GENERAL