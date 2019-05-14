#include "fintrf.h"  
      

    
    !::
    !::prhs(1) is an array of structs containing the tile information, size n_tiles
    !::prhs(2) is an integer representing the number of tiles, n_tiles
    !::prhs(3) is an array of structs containing the state function information, size n_statefunctions
    !::prhs(4) is an integer representing the number of state functions, n_statefunctions
    !::prhs(5) is a double defining the error tolerance for the iteration
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
      type(MagStatStateFunction),allocatable,dimension(:) :: stateFunctions
      
      integer*4 :: n_tiles,n_statefunctions,max_ite
      real*8 :: err_max,T,resumeIteration
      
!     Check for proper number of arguments. 
      if( nrhs .lt. 6) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nInput',
     +                           'at least six inputs are required.')
      elseif(nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nOutput',
     +                           'Too many output arguments.')
      endif

      max_ite = 100
      !! default value is zero, i.e. do not resume iteration
      resumeIteration = 0.
      
!Check the type of the inputs      
      if ( .NOT. mxIsStruct(prhs(1)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input one should be a struct')      
      else if ( .NOT. mxIsInt32(prhs(2)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input two should be an integer')      
      else if ( .NOT. mxIsStruct(prhs(3)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input three should be a struct')      
      else if ( .NOT. mxIsInt32(prhs(4)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input four should be an integer')      
      else if ( .NOT. mxIsDouble(prhs(5)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input five should be a double')      
      else if ( .NOT. mxIsDouble(prhs(6)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input six should be a double')      
      else if ( nrhs .gt. 6 ) then
           if ( .NOT. mxIsInt32(prhs(7)) ) then
                call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input seven (optional) should be an integer')                      
           else
               call mxCopyPtrToInteger4(mxGetPr(prhs(7)), max_ite,1)
           endif            
           if ( nrhs .gt. 7 ) then          
               if ( .NOT. mxIsDouble(prhs(8)) ) then
                     call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input eight (optional) should be a real')                      
                 else
                     call mxCopyPtrToReal8(mxGetPr(prhs(8)), resumeIteration,1)
                 endif            
           endif          
      endif                  
      
      
      !::Copy the number of tiles
      call mxCopyPtrToInteger4(mxGetPr(prhs(2)), n_tiles,1)
      
      !::Copy the number of statefunctions
      call mxCopyPtrToInteger4(mxGetPr(prhs(4)), n_statefunctions,1)

      !::Copy the temperature
      call mxCopyPtrToReal8(mxGetPr(prhs(5)), T, 1 )
      
      !::Copy the max error
      call mxCopyPtrToReal8(mxGetPr(prhs(6)), err_max, 1 )
      
      !::allocate the tiles
      allocate( cylTile(n_tiles) )
      !::allocate the state functions
      allocate( stateFunctions(n_stateFunctions) )
      !::Copy the input parameters      
      call loadMagTile( prhs(1), cylTile, n_tiles )
      !::copy the state functions
      call loadMagStateFunction( prhs(3), stateFunctions, n_statefunctions )
            
      !::do the calculation
      call iterateMagnetization( cylTile, n_tiles,stateFunctions, n_statefunctions, T, err_max, max_ite, 
     +     displayIteration_Matlab, resumeIteration )
      
      !::Return the updated struct array to matlab
      call returnMagTile( cylTile, n_tiles, plhs(1) )
      
      deallocate(cylTile,stateFunctions)
            
      end subroutine
      
      
      
            
      
      
      