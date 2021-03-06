#include "fintrf.h"  
      

    
    !::
    !::prhs(1) is a struct array of size n_tiles containing the tile information
    !::prhs(2) is an [n,3] array containing the pts at which the solution is required
    !::prhs(3) is an integer representing the number of tiles, n_tiles
    !::prhs(4) is an integer representing the number of points, n
    !::plhs(1) is a [n,3] double matrix containing the output magnetic field
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use MagTileIO
      use DemagFieldGetSolution
      implicit none

!     mexFunction arguments:
      mwPointer :: plhs(*), prhs(*)
      integer*4 :: nlhs, nrhs

!     Function declarations:
      mwPointer :: mxGetPr, mxCreateNumericArray
      integer*4 mxIsDouble, mxIsInt32, mxIsStruct, mxClassIDFromClassName    
      mwSize sx            
      integer*4 ComplexFlag,classid

!     Pointers to input/output mxArrays:
      real*8,dimension(:,:),allocatable :: pts,H            
      type(MagTile),allocatable,dimension(:) :: tile      
      integer*4 :: n_ele,n_tiles
      mwSize,dimension(3) :: dims
      mwIndex :: i
      mwSize sizevars
      
!     Check for proper number of arguments. 
      if( nrhs .ne. 4) then
         call mexErrMsgIdAndTxt ('MATLAB:MagTense_mex:nInput','Four inputs are required.')
      elseif(nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:MagTense_mex:nOutput','Too many output arguments.')
      endif


!Check the type of the inputs      
      if ( .NOT. mxIsStruct(prhs(1)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input one should be a struct')
      else if ( .NOT. mxIsDouble(prhs(2)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input two should be a double')      
      else if ( .NOT. mxIsInt32(prhs(3)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input three should be a integer')
      else if ( .NOT. mxIsInt32(prhs(4)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input four should be a integer')
      endif                  
      
      !::Copy the number of tiles
      sizevars = 1
      call mxCopyPtrToInteger4(mxGetPr(prhs(3)), n_tiles,sizevars )
      
      !::allocate the tiles
      allocate( tile(n_tiles) )
      !::Copy the input parameters      
      call loadMagTile( prhs(1), tile, n_tiles )
      
      !::Copy the points at which the solution is required
      sizevars = 1
      call mxCopyPtrToInteger4(mxGetPr(prhs(4)), n_ele, sizevars )
      allocate(pts(n_ele,3),H(n_ele,3))      
      sx = n_ele * 3
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), pts, sx )
      
      !::do the calculation
      H(:,:) = 0
      call getFieldFromTiles( tile, H, pts, n_tiles, n_ele )
      
      !::Load the result back to Matlab      
      ComplexFlag = 0
      
      dims(1) = n_ele
      dims(2) = 3      
      classid = mxClassIDFromClassName('double')

      !Return the H field
      sx = 3 * n_ele
      sizevars = 2
      plhs(1) = mxCreateNumericArray( sizevars , dims, classid, ComplexFlag )
      call mxCopyReal8ToPtr( H, mxGetPr( plhs(1) ), sx )
      
      deallocate(tile,pts,H)
            
      end subroutine
      
      
      
            
      
      
      