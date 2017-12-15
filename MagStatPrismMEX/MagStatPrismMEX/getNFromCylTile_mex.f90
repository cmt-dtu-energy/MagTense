#include "fintrf.h"  
      

    
    !::
    !::prhs(1) is a struct containing the tile information
    !::prhs(2) is an [n,3] array containing the pts at which the solution is required
    !::prhs(3) is an integer representing the number of points, n
    !::plhs(1) is a [3,3,n] double array containing the output tensor field
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
      real*8,dimension(:,:,:),allocatable :: N
      real*8,dimension(3) :: rectDims
      type(MagTile):: cylTile
      mwPointer :: r0Ptr,theta0Ptr,z0Ptr,drPtr,dthetaPtr,dzPtr,MPtr,u_eaPtr,u_oa1Ptr,u_oa2Ptr,mur_eaPtr,mur_oaPtr,MremPtr
      mwPointer :: tileTypePtr,offsetPtr,rotAnglesPtr,rectDimsPtr
      integer*4 :: n_ele
      mwSize,dimension(3) :: dims
      mwIndex :: i
      
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
      sx = 1            
      
        r0Ptr = mxGetField(prhs(1),sx,'r0')
        theta0Ptr = mxGetField(prhs(1),sx,'theta0')
        z0Ptr = mxGetField(prhs(1),sx,'z0')
        drPtr = mxGetField(prhs(1),sx,'dr')
        dthetaPtr = mxGetField(prhs(1),sx,'dtheta')
        dzPtr = mxGetField(prhs(1),sx,'dz')
        MPtr = mxGetField(prhs(1),sx,'M')
        u_eaPtr = mxGetField(prhs(1),sx,'u_ea')
        u_oa1Ptr = mxGetField(prhs(1),sx,'u_oa1')
        u_oa2Ptr = mxGetField(prhs(1),sx,'u_oa2')
        mur_eaPtr = mxGetField(prhs(1),sx,'mu_r_ea')
        mur_oaPtr = mxGetField(prhs(1),sx,'mu_r_oa')
        MremPtr = mxGetField(prhs(1),sx,'Mrem')
        tileTypePtr = mxGetField(prhs(1),sx,'tileType')
        offsetPtr = mxGetField(prhs(1),sx,'offset')
        rotAnglesPtr = mxGetField(prhs(1),sx,'rotAngles')
        rectDimsPtr =  mxGetField(prhs(1),sx,'abc')
      
        call mxCopyPtrToReal8(mxGetPr(r0Ptr), cylTile%r0, 1 )
        call mxCopyPtrToReal8(mxGetPr(theta0Ptr), cylTile%theta0, 1 )
        call mxCopyPtrToReal8(mxGetPr(z0Ptr), cylTile%z0, 1 )
        call mxCopyPtrToReal8(mxGetPr(drPtr), cylTile%dr, 1 )
        call mxCopyPtrToReal8(mxGetPr(dthetaPtr), cylTile%dtheta, 1 )
        call mxCopyPtrToReal8(mxGetPr(dzPtr), cylTile%dz, 1 )
        call mxCopyPtrToReal8(mxGetPr(MPtr), cylTile%M, 3 )
        call mxCopyPtrToReal8(mxGetPr(u_eaPtr), cylTile%u_ea, 3 )
        call mxCopyPtrToReal8(mxGetPr(u_oa1Ptr), cylTile%u_oa1, 3 )
        call mxCopyPtrToReal8(mxGetPr(u_oa2Ptr), cylTile%u_oa2, 3 )
        call mxCopyPtrToReal8(mxGetPr(mur_eaPtr), cylTile%mu_r_ea, 1 )
        call mxCopyPtrToReal8(mxGetPr(mur_oaPtr), cylTile%mu_r_oa, 1 )
        call mxCopyPtrToReal8(mxGetPr(MremPtr), cylTile%Mrem, 1 )
        call mxCopyPtrToInteger4(mxGetPr(tileTypePtr), cylTile%tileType, 1 )
        call mxCopyPtrToReal8(mxGetPr(offsetPtr), cylTile%offset, 3 )
        call mxCopyPtrToReal8(mxGetPr(rotAnglesPtr), cylTile%rotAngles, 3 )
        call mxCopyPtrToReal8(mxGetPr(rectDimsPtr), rectDims, 3 )
        cylTile%a = rectDims(1)
        cylTile%b = rectDims(2)
        cylTile%c = rectDims(3)
      
      !::copy the number of points where the tensor field is required
      call mxCopyPtrToInteger4(mxGetPr(prhs(3)), n_ele,1)
      allocate(pts(n_ele,3),N(3,3,n_ele),H(n_ele,3))      
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), pts, n_ele*3)
      
      !::do the calculation
      N(:,:,:) = 0
      H(:,:) = 0      
      call getFieldFromCylTile( cylTile, H, pts, n_ele, N )
      
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
      
      
      
            
      
      
      