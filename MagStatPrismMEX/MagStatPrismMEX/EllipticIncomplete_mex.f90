#include "fintrf.h"  
    !::
    !::Four inputs are required: 
    !::prhs(1) is the k parameter (with the relationship m=k^2 and Matlab conventionally uses the m-parameter while Maple uses the k-parameter)
    !::prhs(2) the phi parameter (in radians)
    !::prhs(3) the n parameter (used only by the incomplete elliptical integral of the third kind
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use SPECIALFUNCTIONS

      implicit none
      
!     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer*4 nlhs, nrhs

!     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      mwPointer mxCreateNumericArray
      integer*4 mxIsNumeric
      integer*4 mxIsDouble
      integer*4 mxIsInt32
      integer*4 mxIsStruct
      integer*4 mxIsCell
      integer*4 mxClassIDFromClassName
    
      mwPointer mxGetField
      mwPointer mxGetFieldByNumber
      
      mwPointer mxGetM, mxGetN
      mwSize mxGetNumberOfDimensions,dimIteOut(1), sx
      mwPointer mxGetDimensions
      
      integer*4 ComplexFlag,classid

!     Pointers to input/output mxArrays:
      mwPointer :: stPr
      
      real*8,dimension(:),allocatable :: HK,phi,C
      real*8,dimension(:),allocatable :: fe,ee,el3
      integer*4 :: i,n_ele
      real*8,parameter :: pi=3.141592653589793
      
!     Check for proper number of arguments. 
      if( nrhs .ne. 4) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nInput',
     +                           'four inputs are required.')
      elseif(nlhs .gt. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nOutput',
     +                           'Too many output arguments.')
      endif

!Check the type of the inputs      
      if ( .NOT. mxIsDouble(prhs(1)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input one should be a double')
      else if ( .NOT. mxIsDouble(prhs(2)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input two should be a double')
      else if ( .NOT. mxIsDouble(prhs(3)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input three should be a double')
      else if ( .NOT. mxIsInt32(prhs(4)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input four should be an integer')
      endif                  
      
      !::Copy the input parameters
      call mxCopyPtrToInteger4(mxGetPr(prhs(4)), n_ele,1)
          allocate(HK(n_ele),phi(n_ele),C(n_ele),fe(n_ele),ee(n_ele),el3(n_ele))
      sx = n_ele
      call mxCopyPtrToReal8(mxGetPr(prhs(1)), HK, sx )
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), phi, sx )
      call mxCopyPtrToReal8(mxGetPr(prhs(3)), C, sx )
      
      
      !::calculate the elliptical integrals
      do i=1,sx
          !call elit ( HK(i), phi(i), fe(i), ee(i) )
          !call elit3 ( phi(i), HK(i), C(i), el3(i) )
          fe(i) = ellf( phi(i), HK(i) )
          ee(i) = elle( phi(i), HK(i) )
          el3(i) = ellpi( phi(i), C(i), HK(i) )
      enddo
      
     
      !::Load the result back to Matlab
      
      ComplexFlag = 0
      
      plhs(1) = mxCreateDoubleMatrix( sx, 1, ComplexFlag )
      plhs(2) = mxCreateDoubleMatrix( sx, 1, ComplexFlag )
      plhs(3) = mxCreateDoubleMatrix( sx, 1, ComplexFlag )
           
      call mxCopyReal8ToPtr( fe, mxGetPr( plhs(1) ), sx )
      call mxCopyReal8ToPtr( ee, mxGetPr( plhs(2) ), sx )
      call mxCopyReal8ToPtr( el3, mxGetPr( plhs(3) ), sx )
            
      end subroutine
      
      
      
            
      
      
      