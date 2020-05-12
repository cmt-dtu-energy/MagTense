#include "fintrf.h"  
      

    
    !::
    !::prhs(1) is a struct containing the tile information
    !::prhs(2) is an [n,3] array containing the pts at which the solution is required
    !::prhs(3) is an integer representing the number of points, n
    !::plhs(1) is a [3,3,n] double array containing the output tensor field
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use DemagFieldGetSolution
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
      real*8,dimension(:,:,:,:),allocatable :: N
      type(MagTile),dimension(1):: tile
      integer*4 :: n_ele
      mwSize,dimension(3) :: dims      
      
!     Check for proper number of arguments. 
      if( nrhs .ne. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:MagTense_mex:nInput','Three inputs are required.')
      elseif(nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:MagTense_mex:nOutput','Too many output arguments.')
      endif

!Check the type of the inputs      
      if ( .NOT. mxIsStruct(prhs(1)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input one should be a struct')
      else if ( .NOT. mxIsDouble(prhs(2)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input two should be a double')      
      else if ( .NOT. mxIsInt32(prhs(3)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input three should be an integer')            
      endif                  
            
      !::Copy the input parameters
      !::Copy the input parameters      
      call loadMagTile( prhs(1), tile, 1 )
      
      !::copy the number of points where the tensor field is required
      call mxCopyPtrToInteger4(mxGetPr(prhs(3)), n_ele,1)
      allocate(pts(n_ele,3),N(1,n_ele,3,3),H(n_ele,3))      
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), pts, n_ele*3)
      
      
      
      
      !::do the calculation
      N(:,:,:,:) = 0
      H(:,:) = 0      
      if ( tile(1)%tileType .eq. tileTypeCylPiece ) then
          call getFieldFromCylTile( tile(1), H, pts, n_ele, N, .false. )
      else if (tile(1)%tileType .eq. tileTypeCircPiece ) then          
          call getFieldFromCircPieceTile( tile(1), H, pts, n_ele, N, .false. )
      else if (tile(1)%tileType .eq. tileTypeCircPieceInverted ) then          
          call getFieldFromCircPieceInvertedTile( tile(1), H, pts, n_ele, N, .false. )
      else if (tile(1)%tileType .eq. tileTypePrism ) then          
          call getFieldFromRectangularPrismTile( tile(1), H, pts, n_ele, N, .false. )
      else if (tile(1)%tileType .eq. tileTypeTetrahedron ) then          
          call getFieldFromTetrahedronTile( tile(1), H, pts, n_ele, N, .false. )
      else if (tile(1)%tileType .eq. tileTypeSphere ) then          
          call getFieldFromSphereTile( tile(1), H, pts, n_ele, N, .false. )
      else if (tile(1)%tileType .eq. tileTypeSpheroid ) then          
          call getFieldFromSpheroidTile( tile(1), H, pts, n_ele, N, .false. )
      endif
      
      
      !::Load the result back to Matlab
      
      ComplexFlag = 0
      
      dims(1) = n_ele
      dims(2) = 3
      dims(3) = 3
      classid = mxClassIDFromClassName('double')
      
      !Return the N field
      plhs(1) = mxCreateNumericArray( 3, dims, classid, ComplexFlag )
      call mxCopyReal8ToPtr( N, mxGetPr( plhs(1) ), 3*3*n_ele )
      
      deallocate(pts,N)
            
      end subroutine
      
      
      
            
      
      
      