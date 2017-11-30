#include "fintrf.h"  
      module model_helper
      
      implicit none                  
       
        contains
        subroutine getRotZ( phi, Rz )
        real*8,intent(in) :: phi
        real*8,dimension(3,3),intent(inout) :: Rz
        Rz(:,:) = 0
    
        Rz(1,1) = cos(phi)
        Rz(1,2) = -sin(phi)
        Rz(2,1) = sin(phi)
        Rz(2,2) = cos(phi)
    
        end subroutine getRotZ
        end module model_helper

    
    !::
    !::prhs(1) is a struct containing the tile information
    !::prhs(2) is the x-coordinate at which the field is required
    !::plhs(1) is a 3x3 double matrix containing the output N-tensor
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use TileNComponents
      use model_helper
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
      real*8,dimension(:),allocatable :: x,phi,r
      real*8,dimension(3,3) :: Rz
      real*8,dimension(3) :: M_orig,M_tmp
      real*8 :: phi_orig,z_orig
      type(CylPiece) :: cylTile
      mwPointer :: r0Ptr,theta0Ptr,z0Ptr,drPtr,dthetaPtr,dzPtr,MPtr
      integer*4 :: n_ele,i
      mwSize,dimension(3) :: dims
      
!     Check for proper number of arguments. 
      if( nrhs .ne. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nInput',
     +                           'three inputs are required.')
      elseif(nlhs .gt. 2) then
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
      
      call mxCopyPtrToReal8(mxGetPr(r0Ptr), cylTile%r0, 1 )
      call mxCopyPtrToReal8(mxGetPr(theta0Ptr), cylTile%theta0, 1 )
      call mxCopyPtrToReal8(mxGetPr(z0Ptr), cylTile%z0, 1 )
      call mxCopyPtrToReal8(mxGetPr(drPtr), cylTile%dr, 1 )
      call mxCopyPtrToReal8(mxGetPr(dthetaPtr), cylTile%dtheta, 1 )
      call mxCopyPtrToReal8(mxGetPr(dzPtr), cylTile%dz, 1 )
      call mxCopyPtrToReal8(mxGetPr(MPtr), cylTile%M, 3 )
      
      call mxCopyPtrToInteger4(mxGetPr(prhs(3)), n_ele,1)
      allocate(pts(n_ele,3),N(3,3,n_ele))      
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), pts, n_ele*3)
      N(:,:,:) = 0
      
      
      !::Run the calculation
      allocate( r(n_ele), x(n_ele), phi(n_ele), H(n_ele,3) )
      x(:) = 0.
      phi(:) = 0.
      !the length of the radius vector
      r = sqrt( pts(:,1)**2 + pts(:,2)**2 )
      !Find the rotation angle. It is assumed that r != 0 in all points since the tensor-field diverges here      
      phi = acos( pts(:,1) / r )
      !Change the sign if y is negative
      where ( pts(:,2) .lt. 0 )
          phi = -phi
      endwhere
      
      !save the original orientation of the tile
      phi_orig = cylTile%theta0
      z_orig = cylTile%z0      
      M_orig = cylTile%M
      do i=1,n_ele
          !Offset the angle
          cylTile%theta0 = phi_orig - phi(i)
          !Offset the z-coordinate
          cylTile%z0 = z_orig - pts(i,3)
          call getN( cylTile, r(i), N(:,:,i) )
          !Get the rotation vector
          call getRotZ( -phi(i), Rz )
          !Rotate the magnetization vector
          M_tmp = matmul( Rz, M_orig )
          !find the field at the rotated point
          H(i,:) = 1./(4.*pi) * matmul( N(:,:,i), M_tmp )
          !get the negative rotation
          call getRotZ( phi(i), Rz )
          !rotate the field back to the original orientation
          H(i,:) = matmul( Rz, H(i,:) )
      enddo
                 
      !::Load the result back to Matlab
      
      ComplexFlag = 0
      
      dims(1) = n_ele
      dims(2) = 3      
      classid = mxClassIDFromClassName('double')
      
      !Return the H field
      plhs(1) = mxCreateNumericArray( 2, dims, classid, ComplexFlag )
      call mxCopyReal8ToPtr( H, mxGetPr( plhs(1) ), 3*n_ele )
      
      !Return the N vector
      dims(1) = 3
      dims(2) = 3
      dims(3) = n_ele
      plhs(2) = mxCreateNumericArray( 3, dims, classid, ComplexFlag )
      
      call mxCopyReal8ToPtr( N, mxGetPr( plhs(2) ), 3*3*n_ele )
            
      end subroutine
      
      
      
            
      
      
      