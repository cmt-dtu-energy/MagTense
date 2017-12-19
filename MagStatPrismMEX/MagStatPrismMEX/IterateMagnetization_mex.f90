#include "fintrf.h"  
      

    
    !::
    !::prhs(1) is an array of structs containing the tile information, size n_tiles
    !::prhs(2) is an integer representing the number of tiles, n_tiles
    !::prhs(3) is a double defining the error tolerance for the iteration
    !::plhs(1) is an array of structs containing the updated tiles, size n_tiles
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use IterateMagnetSolution
      use MagTileIO
      implicit none
      
!     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer*4 nlhs, nrhs

!     Function declarations:                        
      integer*4 mxIsDouble
      integer*4 mxIsInt32
      integer*4 mxIsStruct      
    
      mwPointer mxGetPr
    
!     Pointers to input/output mxArrays:      
      type(MagTile),allocatable,dimension(:) :: cylTile      
      
      integer*4 :: n_tiles      
      real*8 :: err_max
      
!     Check for proper number of arguments. 
      if( nrhs .ne. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nInput',
     +                           'two inputs are required.')
      elseif(nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nOutput',
     +                           'Too many output arguments.')
      endif

!Check the type of the inputs      
      if ( .NOT. mxIsStruct(prhs(1)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input one should be a struct')      
      else if ( .NOT. mxIsInt32(prhs(2)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input two should be an integer')      
      else if ( .NOT. mxIsDouble(prhs(3)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input three should be a double')      
      endif                  
      
      
      !::Copy the number of tiles
      call mxCopyPtrToInteger4(mxGetPr(prhs(2)), n_tiles,1)
      
      !::Copy the max error
      call mxCopyPtrToReal8(mxGetPr(prhs(3)), err_max, 1 )
      
      !::allocate the tiles
      allocate( cylTile(n_tiles) )
      !::Copy the input parameters      
      call loadMagTile( prhs(1), cylTile, n_tiles )
            
      !::do the calculation
      call iterateMagnetization( cylTile, n_tiles, err_max )
      
      !::Return the updated struct array to matlab
      call returnMagTile( cylTile, n_tiles, plhs(1) )
      
      deallocate(cylTile)
            
      end subroutine
      
      
      
            
      
      
      