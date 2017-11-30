#include "fintrf.h"  
    
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
      integer*4 :: i
      real*8,parameter :: pi=3.141592654
      
!     Check for proper number of arguments. 
      if( nrhs .ne. 4) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nInput',
     +                           'three inputs are required.')
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
      call mxCopyPtrToInteger4(mxGetPr(prhs(4)), sx,1)
      allocate(HK(sx),phi(sx),C(sx),fe(sx),ee(sx),el3(sx))
      !sx = 1
      call mxCopyPtrToReal8(mxGetPr(prhs(1)), HK, sx )
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), phi, sx )
      call mxCopyPtrToReal8(mxGetPr(prhs(3)), C, sx )
      
      !::convert from radians to degrees
      phi = 180 / pi * phi
      !::calculate the elliptical integrals
      do i=1,sx
          call elit ( HK(i), phi(i), fe(i), ee(i) )
          call elit3 ( phi(i), HK(i), C(i), el3(i) )
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
      
      
      
            
      
      
      