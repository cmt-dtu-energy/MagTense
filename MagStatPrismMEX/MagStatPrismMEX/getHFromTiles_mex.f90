#include "fintrf.h"  
      

    
    !::
    !::prhs(1) is a struct containing the tile information
    !::prhs(2) is an [n,3] array containing the pts at which the solution is required
    !::prhs(3) is an integer representing the number of points, n
    !::plhs(1) is a [n,3] double matrix containing the output magnetic field
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use MagStat2GetSolution
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
      mwSize mxGetNumberOfDimensions, sx
      mwPointer mxGetDimensions
      
      integer*4 ComplexFlag,classid

!     Pointers to input/output mxArrays:
      mwPointer :: stPr
      real*8,dimension(:,:),allocatable :: pts,H      
      type(MagTile),allocatable,dimension(:) :: cylTile
      mwPointer :: r0Ptr,theta0Ptr,z0Ptr,drPtr,dthetaPtr,dzPtr,MPtr,u_eaPtr,u_oa1Ptr,u_oa2Ptr,mur_eaPtr,mur_oaPtr,MremPtr
      integer*4 :: n_ele,n_tiles
      mwSize,dimension(3) :: dims
      mwIndex :: i
      
!     Check for proper number of arguments. 
      if( nrhs .ne. 4) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nInput',
     +                           'four inputs are required.')
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
      else if ( .NOT. mxIsInt32(prhs(4)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType', 'Input four should be a integer')
      endif                  
      
      !::Copy the number of tiles
      call mxCopyPtrToInteger4(mxGetPr(prhs(3)), n_tiles,1)
      
      !::allocate the tiles
      allocate( cylTile(n_tiles) )
      !::Copy the input parameters
      sx = 1            
      do i=1,n_tiles
          r0Ptr = mxGetField(prhs(1),i,'r0')
          theta0Ptr = mxGetField(prhs(1),i,'theta0')
          z0Ptr = mxGetField(prhs(1),i,'z0')
          drPtr = mxGetField(prhs(1),i,'dr')
          dthetaPtr = mxGetField(prhs(1),i,'dtheta')
          dzPtr = mxGetField(prhs(1),i,'dz')
          MPtr = mxGetField(prhs(1),i,'M')
          u_eaPtr = mxGetField(prhs(1),i,'u_ea')
          u_oa1Ptr = mxGetField(prhs(1),i,'u_oa1')
          u_oa2Ptr = mxGetField(prhs(1),i,'u_oa2')
          mur_eaPtr = mxGetField(prhs(1),i,'mu_r_ea')
          mur_oaPtr = mxGetField(prhs(1),i,'mu_r_oa')
          MremPtr = mxGetField(prhs(1),i,'Mrem')
      
          call mxCopyPtrToReal8(mxGetPr(r0Ptr), cylTile(i)%r0, 1 )
          call mxCopyPtrToReal8(mxGetPr(theta0Ptr), cylTile(i)%theta0, 1 )
          call mxCopyPtrToReal8(mxGetPr(z0Ptr), cylTile(i)%z0, 1 )
          call mxCopyPtrToReal8(mxGetPr(drPtr), cylTile(i)%dr, 1 )
          call mxCopyPtrToReal8(mxGetPr(dthetaPtr), cylTile(i)%dtheta, 1 )
          call mxCopyPtrToReal8(mxGetPr(dzPtr), cylTile(i)%dz, 1 )
          call mxCopyPtrToReal8(mxGetPr(MPtr), cylTile(i)%M, 3 )
          call mxCopyPtrToReal8(mxGetPr(u_eaPtr), cylTile(i)%u_ea, 3 )
          call mxCopyPtrToReal8(mxGetPr(u_oa1Ptr), cylTile(i)%u_oa1, 3 )
          call mxCopyPtrToReal8(mxGetPr(u_oa2Ptr), cylTile(i)%u_oa2, 3 )
          call mxCopyPtrToReal8(mxGetPr(mur_eaPtr), cylTile(i)%mu_r_ea, 1 )
          call mxCopyPtrToReal8(mxGetPr(mur_oaPtr), cylTile(i)%mu_r_oa, 1 )
          call mxCopyPtrToReal8(mxGetPr(MremPtr), cylTile(i)%Mrem, 1 )
      enddo
      
      
      call mxCopyPtrToInteger4(mxGetPr(prhs(4)), n_ele,1)
      allocate(pts(n_ele,3),H(n_ele,3))      
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), pts, n_ele*3)
      
      !::do the calculation
      H(:,:) = 0
      call getFieldFromTiles( cylTile, H, pts, n_tiles, n_ele )
      !::Load the result back to Matlab
      
      ComplexFlag = 0
      
      dims(1) = n_ele
      dims(2) = 3      
      classid = mxClassIDFromClassName('double')
      
      !Return the H field
      plhs(1) = mxCreateNumericArray( 2, dims, classid, ComplexFlag )
      call mxCopyReal8ToPtr( H, mxGetPr( plhs(1) ), 3*n_ele )
      
      deallocate(cylTile,pts,H)
            
      end subroutine
      
      
      
            
      
      
      