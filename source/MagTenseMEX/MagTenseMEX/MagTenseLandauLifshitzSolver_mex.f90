#include "fintrf.h"  
    
    
    !-----------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Entry point for getting MagTense so solve the Landau-Lifshitz equation
    !> @param[in] nlhs no. of output arguments. Should be 1 (the output solution)
    !> @param[in] plhs pointer array to the output arguments
    !> @param[in] nrhs no. of input arguments. Should be 2 (the input problem and starting solution)
    !> @param[in] prhs pointer array to the input arguments. 
    !> Prhs(1) is the problem struct while prhs(2) is the initial solution struct
    !> Plhs(1) is the output solution struct
    !-----------------------------------
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use LandauLifshitzSolution
      use integrationDataTypes
      use MagTenseMicroMagIO
      
      implicit none
      mwPointer :: plhs(*), prhs(*)             !>Defines the input/output pointers
      mwPointer :: mxGetPr, mxCreateNumericArray!>Defines functions for interacting with MEX
      integer*4 :: nlhs, nrhs                   !>Defines the no. of I/O arguments       
                  
      mwSize :: sx                              !>MEX integer for keeping the size correct
      mwSize,dimension(3) :: dims               !>Used for copying back to ML
      integer*4 ComplexFlag,classid             !>Flags for communicating with ML
      integer*4 mxIsDouble, mxIsInt32, mxIsStruct, mxClassIDFromClassName  !>Various MEX functions
      type(MicroMagProblem) :: problem          !> The problem structure
      type(MicroMagSolution) :: solution        !> The solution structure
      
      
    
      !Check the input parameters
      if ( nrhs .ne. 2 ) then
        call mexErrMsgIdAndTxt ('MATLAB:MagTensePDE:nInput',
     +                           'Two inputs are required.')
      elseif ( nlhs .ne. 1 ) then
          call mexErrMsgIdAndTxt ('MATLAB:MagTensePDE:nOutput',
     +                           'One output is required.')
      endif
      
      
      if ( .NOT. mxIsStruct(prhs(1)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input one should be a struct')
      else if ( .NOT. mxIsStruct(prhs(2)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input two should be a struct')
      endif
      
      !Load the problem from Matlab into Fortran
      call loadMicroMagProblem( prhs(1), problem )
            
      !Call the ODE solver
      call SolveLandauLifshitzEquation( problem, solution )    
    
      call returnMicroMagSolution( solution, plhs(1) )
      
      
      
      end subroutine 