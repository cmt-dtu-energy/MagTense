#include "fintrf.h"  
      

    
    !::
    !::prhs(1) is a struct containing the tile information
    !::prhs(2) is an [n,3] array containing the pts at which the solution is required
    !::prhs(3) is an integer representing the number of points, n
    !::plhs(1) is a [3,3,n] double array containing the output tensor field
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use MagStat2GetSolution
      use MagTileIO
      implicit none
      
!     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer*4 nlhs, nrhs

!     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateNumericArray
      integer*4 mxIsNumeric, mxIsDouble, mxIsInt32, mxIsStruct, mxClassIDFromClassName
          
      integer*4 ComplexFlag,classid

!     Pointers to input/output mxArrays:
      real*8,dimension(:,:),allocatable :: pts,H
      real*8,dimension(:,:,:),allocatable :: N
      type(MagTile),dimension(1):: cylTile
      integer*4 :: n_ele
      mwSize,dimension(3) :: dims      
      
!     Check for proper number of arguments. 
      if( nrhs .ne. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nInput',
     +                           'three inputs are required.')
      elseif(nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nOutput',
     +                           'Too many output arguments.')
      endif

!Check the type of the inputs      
      if ( .NOT. mxIsStruct(prhs(1)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input one should be a struct')
      else if ( .NOT. mxIsDouble(prhs(2)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input two should be a double')      
      else if ( .NOT. mxIsInt32(prhs(3)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input three should be a integer')      
      endif                  
            
      !::Copy the input parameters
      !::Copy the input parameters      
      call loadMagTile( prhs(1), cylTile, 1 )
      
      !::copy the number of points where the tensor field is required
      call mxCopyPtrToInteger4(mxGetPr(prhs(3)), n_ele,1)
      allocate(pts(n_ele,3),N(n_ele,3,3),H(n_ele,3))      
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), pts, n_ele*3)
      
      !::do the calculation
      N(:,:,:) = 0
      H(:,:) = 0      
      if ( cylTile(1)%tileType .eq. tileTypeCylPiece ) then
            call getFieldFromCylTile( cylTile(1), H, pts, n_ele, N, .false. )
      else if (cylTile(1)%tileType .eq. tileTypeCircPiece ) then          
          call getFieldFromCircPieceTile( cylTile(1), H, pts, n_ele, N, .false. )
      endif
      
      
      !::Load the result back to Matlab
      
      ComplexFlag = 0
      
      dims(1) = 3
      dims(2) = 3
      dims(3) = n_ele
      classid = mxClassIDFromClassName('double')
      
      !Return the N field
      plhs(1) = mxCreateNumericArray( 3, dims, classid, ComplexFlag )
      call mxCopyReal8ToPtr( N, mxGetPr( plhs(1) ), 3*3*n_ele )
      
      deallocate(pts,N)
            
      end subroutine
      
      
      
            
      
      
      