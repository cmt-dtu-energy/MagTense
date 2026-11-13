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

        !! test if we are inside Matlab or called from python or a stand-alone
        !! nothing bad should happen if in fact we are called from ML but the file somehow
        !! exists - the written output will just not be shown to the user

#if USE_MATLAB
        call displayMatlabMessage( mess )
#else
        write(*,*) mess
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


    !>----------------------------------------
    !> Rasmus Bjørk, rabj@dtu.dk, October 2026
    !> Stops the calculation because of an error in the input or an unsupported option.
    !> Inside a MEX file an error stop terminates the whole MATLAB process, and the message
    !> printed just before it is lost with the command window. There the message is raised as
    !> a MATLAB error instead, which unwinds the MEX call and returns to the MATLAB prompt.
    !> Elsewhere (python, stand-alone) it is an error stop.
    !> Must not be called from inside an OpenMP parallel region, as the MATLAB error can only
    !> be raised from the thread that MATLAB called the MEX file on.
    !> @params[in] mess The error message
    !>----------------------------------------
    subroutine abortMagTense( mess )
        character(*),intent(in) :: mess

#if USE_MATLAB
        call mexErrMsgIdAndTxt( 'MagTense:abort', mess )
#endif
        error stop mess

    end subroutine abortMagTense

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

        !! mexCallMATLAB does not take ownership of prhs. MATLAB reclaims the
        !! string when the MEX call returns, but a solve is one long call that
        !! emits a message per field step, so release it here instead.
        call mxDestroyArray( prhs_cb(1) )

    end subroutine displayMatlabMessage

#endif

end module IO_GENERAL
