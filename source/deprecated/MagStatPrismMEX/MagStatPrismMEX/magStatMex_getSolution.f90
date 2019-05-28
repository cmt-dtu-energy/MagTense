
      
#include "fintrf.h"
      
      

!     Gateway routine
!     The prhs is a pointer array with the following content
!     prhs(1)  solPts (l,3), solution points
!     prhs(2)  rotAngles (n,3), rotation angles for each prism
!     prhs(3)  M (n,3) magnetization for each prism
!     prhs(4)  dims (n,3) dimensions of each prism
!     prhs(5)  pos(n,3) positions of each prism
!     prhs(6)  n, integer
!     prhs(7)  l, integer
!     prhs(8)  spacedim, integer
!     prhs(9)  Hext (l,3), external field
!     prhs(10) formType, the type of object
      
!     plhs(1) Hout(l,3) output field

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      USE MagPrism_Call
      USE Parameters_CALL
!     Declarations
      implicit none
      
!     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

!     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      mwPointer mxCreateNumericArray
      integer mxIsNumeric
      integer mxIsDouble
      integer mxIsInt32
      integer mxIsStruct
      integer mxIsCell
      integer mxClassIDFromClassName
      integer writeToConsole
      mwPointer mxGetField
      mwPointer mxGetFieldByNumber
      
      mwPointer mxGetM, mxGetN
      mwSize mxGetNumberOfDimensions,dimIteOut(1)
      mwPointer mxGetDimensions
      
      integer*4 ComplexFlag,d_dims

!     Pointers to input/output mxArrays:
      mwPointer :: stPr
      integer*4 :: nn, mm, ll, spacedim
      real*8,dimension(:,:),allocatable :: pos,dims,solPts,Hout,Hext
      real*8,dimension(:,:),allocatable :: M,rotAngles
      real*8,dimension(:,:,:),allocatable :: rotMat,rotMatInv
      integer*4,dimension(:),allocatable :: formType
      
      real*8,dimension(2) :: dt
      real*8 :: t1,t2,t3
      
	  
	  !-----------------------------------------------------------------------	  	       
	  !   open(11,file='version_magStatMEX_getSolution.txt',
   !  +   status='unknown',access='sequential',form='formatted',
   !  +                       position='rewind',action='write')     
	  !write(11,*) "magStatMEX_getSolution compiled on:"
   !   write(11,*) __DATE__
   !   write(11,*) __TIME__
   !   close (11)
      !-----------------------------------------------------------------------
      
      !-----------------------------------------------------------------------
!     Check for proper number of arguments. 
      if(nrhs .ne. 10) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nInput',
     +                           '8 inputs are required.')
      elseif(nlhs .gt. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nOutput',
     +                           'Too many output arguments.')
      endif

!Check the type of the inputs      
      if ( .NOT. mxIsDouble(prhs(1)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input one should be double')
      elseif ( .NOT. mxIsDouble(prhs(2)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input two should be double')      
      elseif ( .NOT. mxIsDouble(prhs(3)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input three should be double')
      elseif ( .NOT. mxIsDouble(prhs(4)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input four should be Double')
      elseif ( .NOT. mxIsDouble(prhs(5)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input five should be double')
      elseif ( .NOT. mxIsInt32(prhs(6)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input six should be an integer')
      elseif ( .NOT. mxIsInt32(prhs(7)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input seven should be an integer')
      elseif ( .NOT. mxIsInt32(prhs(8)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input eight should be an integer')      
      elseif ( .NOT. mxIsDouble(prhs(9)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input nine should be an double')     
      elseif ( .NOT. mxIsInt32(prhs(10)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input ten should be an integer') 
      endif
      
      call mxCopyPtrToInteger4(mxGetPr(prhs(6)) , nn, 1)
      call mxCopyPtrToInteger4(mxGetPr(prhs(7)), ll, 1)
      call mxCopyPtrToInteger4(mxGetPr(prhs(8)), spacedim, 1)
      
      
      allocate( solPts(ll,3), Hout(ll,3), Hext(ll,3) )
      allocate( M(nn,3), dims(nn,3), pos(nn,3), rotAngles(nn,3) )
      allocate( formType(nn) )
       
      Hout(:,:) = 0
            
      
      call mxCopyPtrToReal8(mxGetPr(prhs(1)), solPts, ll * 3 )
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), rotAngles, nn * 3 )
      call mxCopyPtrToReal8(mxGetPr(prhs(3)), M, nn * 3 )
      call mxCopyPtrToReal8(mxGetPr(prhs(4)), dims, nn * 3 )
      call mxCopyPtrToReal8(mxGetPr(prhs(5)), pos, nn * 3 )
      call mxCopyPtrToReal8(mxGetPr(prhs(9)), Hext, ll * 3)
      call mxCopyPtrToInteger4(mxGetPr(prhs(10)) , formType, nn)
      !CALL CPU_TIME(t1)
      !::Call the subroutine to find the rotations matrices
      call findRotationMatrices( rotAngles, nn, rotMat, rotMatInv)
      !CALL CPU_TIME(t2)
      !::call the subroutine to find the solution at the specified points
      call getSolution( solPts, ll, M, dims, pos, nn, spacedim, Hout, 
     +                  Hext, rotMat, rotMatInv, formType )
      !CALL CPU_TIME(t3)
      
      !dt(1) = t2-t1;
      !dt(2) = t3-t2;
      
      !::Load the result back to Matlab
            
      d_dims = 3;
      ComplexFlag = 0;
      
      plhs(1) = mxCreateDoubleMatrix( ll, d_dims, ComplexFlag )
      
      call mxCopyReal8ToPtr( Hout ,mxGetPr( plhs(1) ), ll * 3 )
      
      
      if(nlhs .eq. 2) then
        plhs(2) = mxCreateDoubleMatrix( 2, 1, ComplexFlag )   
        call mxCopyReal8ToPtr( dt, mxGetPr( plhs(2) ), 2 )
      endif
      
      deallocate( solpts, Hout, M, dims, pos, rotAngles )
      
      end subroutine
      
      
      

      
      