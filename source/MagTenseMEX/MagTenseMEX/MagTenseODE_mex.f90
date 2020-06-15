#include "fintrf.h"  
    
    
    !-----------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Entry point for the MagTense ODE solvers
    !> @param[in] nlhs no. of output arguments. Should be 2
    !> @param[in] plhs pointer array to the output arguments
    !> @param[in] nrhs no. of input arguments. Should be 4
    !> @param[in] prhs pointer array to the input arguments. 
    !> Prhs(1) is the time array while prhs(2) is the initial values of y
    !> Prhs(3) is the no. of time points and prhs(4) is the no. of equations
    !-----------------------------------
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use ODE_Solvers
      use integrationDataTypes
    
      implicit none
      mwPointer :: plhs(*), prhs(*)             !>Defines the input/output pointers
      mwPointer :: mxGetPr, mxCreateNumericArray!>Defines functions for interacting with MEX
      integer*4 :: nlhs, nrhs                   !>Defines the no. of I/O arguments       
      procedure(dydt_fct), pointer :: fct       !>Procedural pointer used for evaluating dydt
      integer*4 :: neq,nt                       !>no. of equations, no. of output times
      real,dimension(:),allocatable :: t,y0     !>Time and y0 arrays
      real,dimension(:),allocatable :: t_out    !>Time out array
      real,dimension(:,:),allocatable :: y_out  !>y output array
      mwSize :: sx, sizevars                    !>MEX integer for keeping the size correct
      mwSize,dimension(3) :: dims               !>Used for copying back to ML
      integer*4 ComplexFlag,classid             !>Flags for communicating with ML
      integer*4 mxIsDouble, mxIsInt32, mxIsStruct, mxClassIDFromClassName  !>Various MEX functions
      
      !Set the function to be integrated and currently only one choice is avaiable
      fct => dydt_ML
    
      !Check the input parameters
      if ( nrhs .ne. 4 ) then
        call mexErrMsgIdAndTxt ('MATLAB:MagTensePDE:nInput','Four inputs are required.')
      elseif ( nlhs .ne. 2 ) then
          call mexErrMsgIdAndTxt ('MATLAB:MagTensePDE:nOutput','Two outputs are required.')
      endif
      
      
      if ( .NOT. mxIsInt32(prhs(3)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input three should be an integer')
      else if ( .NOT. mxIsInt32(prhs(4)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input four should be an integer')
      endif
      
      
      !Load the neq and nt
      sizevars = 1
      call mxCopyPtrToInteger4(mxGetPr(prhs(3)), nt, sizevars )
      call mxCopyPtrToInteger4(mxGetPr(prhs(4)), neq, sizevars )
    
      !Allocate the t and y0 arrays
      allocate(t(nt),y0(neq))
      !Allocate the output time and y arrays
      allocate(t_out(nt),y_out(neq,nt))
      t_out(:) = 0
      y_out(:,:) = 0

      !Copy the time array
      sx = nt
      call mxCopyPtrToReal8(mxGetPr(prhs(1)), t, sx )
      !Copy the y0 array
      sx = neq
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), y0, sx )
      
      
      !Call the ODE solver
      call MagTense_ODE( fct, t, y0, t_out, y_out )
    
    
      !Copy the solution back to Matlab
      ComplexFlag = 0
      classid = mxClassIDFromClassName('double')
      
      !The time output
      dims(1) = nt
      plhs(1) = mxCreateNumericArray( 1, dims, classid, ComplexFlag )
      sx = nt      
      call mxCopyReal8ToPtr( t_out, mxGetPr( plhs(1) ), sx )
      
      !The y output
      dims(1) = neq
      dims(2) = nt
      plhs(2) = mxCreateNumericArray( 2, dims, classid, ComplexFlag )
      sx = neq * nt      
      call mxCopyReal8ToPtr( y_out, mxGetPr( plhs(2) ), sx )
      
      
      !Clean-up
      deallocate(t,y0,t_out,y_out)
      
      end subroutine 