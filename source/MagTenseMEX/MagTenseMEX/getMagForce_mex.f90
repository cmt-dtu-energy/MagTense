#include "fintrf.h"  
      

    
    !::
    !::prhs(1) is a struct array of size n_tiles containing the tile information    
    !::prhs(2) is an integer representing the number of tiles, n_tiles
    !::prhs(3) is a struct representing the surface over which to integrate
    !::plhs(1) is a [n,3] double matrix containing the output magnetic field
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use MagTileIO
      use MagForceIO
      use MagneticForce
      implicit none
      
!     mexFunction arguments:
      mwPointer :: plhs(*), prhs(*)
      integer*4 :: nlhs, nrhs

!     Function declarations:
      mwPointer :: mxGetPr, mxCreateNumericArray
      integer*4 mxIsDouble, mxIsInt32, mxIsStruct, mxClassIDFromClassName    
      mwSize sx            
      integer*4 ComplexFlag,classid

      integer*4,dimension(18,2) :: ier, neval
      
!     Pointers to input/output mxArrays:        
      type( magStatModel ),pointer :: magModel
      type(surf_carth) :: surf
      real*8,dimension(3) :: Fout
      mwSize,dimension(3) :: dims
      mwIndex :: i
      
!     Check for proper number of arguments. 
      if( nrhs .ne. 3 ) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nInput',
     +                           'three inputs are required.')
      elseif(nlhs .gt. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nOutput',
     +                           'Too many output arguments.')
      endif

!Check the type of the inputs      
      if ( .NOT. mxIsStruct(prhs(1)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input one should be a struct')
      else if ( .NOT. mxIsInt32(prhs(2)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input two should be an integer')      
      else if ( .NOT. mxIsStruct(prhs(3)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input three should be a struct')      
      endif                  
      
      !::Load the integration surface
      call loadMagForceSurf( prhs(3), surf )
      
      !::allocate the magstatmodel
      allocate(magModel)
      !::Copy the number of tiles
      call mxCopyPtrToInteger4(mxGetPr(prhs(2)), magModel%n_tiles,1)
      
      !::allocate the tiles
      allocate( magModel%tiles(magModel%n_tiles) )
      !::Copy the input parameters      
      call loadMagTile( prhs(1), magModel%tiles, magModel%n_tiles )
      Fout(:) = 0
      call getForce( magModel, surf, Fout, ier, neval )
      
      !::Load the result back to Matlab      
      ComplexFlag = 0
      
      dims(1) = 1
      dims(2) = 3      
      classid = mxClassIDFromClassName('double')
      
      !Return the force vector
      sx = 3
      plhs(1) = mxCreateNumericArray( 2, dims, classid, ComplexFlag )
      call mxCopyReal8ToPtr( Fout, mxGetPr( plhs(1) ), sx )
      
      if ( nrhs .ge. 2 ) then
          !Return the integer errors
          dims(1) = 18
          dims(2) = 2
          classid = mxClassIDFromClassName('int32')
          plhs(2) = mxCreateNumericArray( 2, dims, classid, ComplexFlag )
          sx = 36
          call mxCopyReal8ToPtr( ier, mxGetPr( plhs(2) ), sx )
      endif
      
      if ( nrhs .ge. 3 ) then
          !Return the no. of evaluations
          dims(1) = 18
          dims(2) = 2
          classid = mxClassIDFromClassName('int32')
          plhs(3) = mxCreateNumericArray( 2, dims, classid, ComplexFlag )
          sx = 36
          call mxCopyReal8ToPtr( neval, mxGetPr( plhs(3) ), sx )
      endif
      
      deallocate(magModel)
            
      end subroutine
      
      
      
            
      
      
      