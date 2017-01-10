#include "fintrf.h"  
      module model_glob
      USE MagPrism_Call
      USE Parameters_CALL

      implicit none
      
      type modelStruct
      integer,dimension(3) :: inds
      real*8 :: z0
      integer :: n,spacedim
      real*8,dimension(:,:),allocatable :: dims,pos,M,rotAngles
      real*8,dimension(:,:,:),allocatable :: rotMat,rotMatInv
      endtype
      
      type(modelStruct) :: model
      
      real*8,parameter :: mu0=4*pi*1e-7
      
      contains
      !::Returns the integrand for the off-diagonal part of the Maxwell stress
      !::tensor dotted with the surface area
      function F_int_off_diag( x, y )
      real*8 :: F_int_off_diag
      real*8,intent(in) :: x,y
      
      real*8,dimension(3) :: solPts,B,Hsol

      !::setup the points where we wish to get the solution
      solPts(model%inds(1)) = x
      solPts(model%inds(2)) = y
      solPts(model%inds(3)) = model%z0;

      !::call the subroutine to find the solution at the specified points
      call getSolution( solPts, 1, model%M, model%dims, 
     +                  model%pos, model%n, model%spacedim,
     +                  Hsol, model%rotMat, model%rotMatInv )
      
                                 
       B = mu0 * Hsol
    
       F_int_off_diag = 
     + 1. /mu0 * B(model%inds(2)) * B(model%inds(3))
    
       return

      end function F_int_off_diag
      
        !::Returns the integrand of the Maxwell stress tensor along the diagonal over
        !::the given surface area. 
        !::
        !::integration variables while .zInd gives the index of the plane
        !::1 = x, 2 = y and 3 = z
        function F_int_diag( x, y )
        real*8 :: F_int_diag
        real*8,intent(in) : x,y
        
        real*8,dimension(3) :: solPts,B,Hsol
        real*8 :: Bnorm
            
        
        solPts(model%inds(1)) = x
        solPts(model%inds(2)) = y
        solPts(model%inds(3)) = model%z0;

    
        call getSolution( solPts, 1, model%M, model%dims, 
     +                  model%pos, model%n, model%spacedim,
     +                  Hsol, model%rotMat, model%rotMatInv )
	
        B = mu0 * Hsol
        Bnorm = sqrt( B(1)**2 + B(2)**2 + B(3)**2 )
                                 
        F_int_off_diag = 1./mu0 * ( B(model%inds(3))**2 - 0.5 * Bnorm**2 )
    
        
        end function F_int_off_diag
      
      end module model_glob


      
      
!     Gateway routine to do a box integral using data points from the MagStatPrism library call
!     The prhs is a pointer array with the following content
!     prhs(1) rotAngles (n,3), rotation angles for each prism
!     prhs(2) M (n,3) magnetization for each prism
!     prhs(3) dims (n,3) dimensions of each prism
!     prhs(4) pos(n,3) positions of each prism
!     prhs(5) n, integer
!     prhs(6) spacedim, integer
!     prhs(7) xlim (2), x bounds
!     prhs(8) ylim (2), y bounds
!     prhs(9) zlim (2), z bounds

!     plhs(1) Fout(3) output force vector [Fx,Fy,Fz]

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      !USE box_integrator
      use model_glob
!      USE MagPrism_Call
!      USE Parameters_CALL

      implicit none
!      type modelStruct
!      integer,dimension(3) :: inds
!      real*8 :: z0
!      integer :: n,spacedim
!      real*8,dimension(:,:),allocatable :: dims,pos,M,rotAngles
!      real*8,dimension(:,:,:),allocatable :: rotMat,rotMatInv
!      endtype
     
      
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
      integer*4 :: nn, mm, spacedim
      real*8,dimension(2) :: xlim,ylim,zlim
      
      real*8,dimension(:,:),allocatable :: Fout !
      
	  
	  !-----------------------------------------------------------------------	  	       
	     open(11,file='version_integral2_mex.txt',
     +   status='unknown',access='sequential',form='formatted',
     +                       position='rewind',action='write')     
	  write(11,*) "integral2_mex compiled on:"
      write(11,*) __DATE__
      write(11,*) __TIME__
      !-----------------------------------------------------------------------
      
      !-----------------------------------------------------------------------
!     Check for proper number of arguments. 
      if(nrhs .ne. 9) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nInput',
     +                           '9 inputs are required.')
      elseif(nlhs .gt. 1) then
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
      elseif ( .NOT. mxIsInt32(prhs(5)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input five should be an integer')
      elseif ( .NOT. mxIsInt32(prhs(6)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input six should be an integer')
      elseif ( .NOT. mxIsDouble(prhs(7)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input seven should be a double')
      elseif ( .NOT. mxIsDouble(prhs(8)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input eight should be a double')      
      elseif ( .NOT. mxIsDouble(prhs(9)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input nine should be a double')      
      endif
      
      !::Allocate input data
      call mxCopyPtrToInteger4(mxGetPr(prhs(5)), model%n, 1)
      call mxCopyPtrToInteger4(mxGetPr(prhs(6)), model%spacedim, 1)
      
      allocate( model%M(model%n,3), model%dims(model%n,3), 
     +                model%pos(model%n,3), model%rotAngles(model%n,3) )
      !::Copy input data
      call mxCopyPtrToReal8(mxGetPr(prhs(1)), model%rotAngles, 
     +                      model%n * 3 )
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), model%M, model%n * 3 )
      call mxCopyPtrToReal8(mxGetPr(prhs(3)), model%dims, model%n * 3 )
      call mxCopyPtrToReal8(mxGetPr(prhs(4)), model%pos, model%n * 3 )
      
      call mxCopyPtrToReal8(mxGetPr(prhs(7)), xlim, 2 )
      call mxCopyPtrToReal8(mxGetPr(prhs(8)), ylim, 2 )
      call mxCopyPtrToReal8(mxGetPr(prhs(9)), zlim, 2 )
      
      !::Call the subroutine to find the rotations matrices
      !call findRotationMatrices( rotAngles, nn, rotMat, rotMatInv)
      
  !    call integral2( fun, xl, xh, yl, yh, res, abserr, neval, ier )
      
      !::call the subroutine to find the solution at the specified points
  !    call getSolution( solPts, ll, M, dims, pos, nn, spacedim, Hout, 
  !   +                  rotMat, rotMatInv )
      
      !::Load the result back to Matlab
      !call boxIntegral( xlim, ylim, zlim, model )
      d_dims = 3;
      ComplexFlag = 0;
      
      !plhs(1) = mxCreateDoubleMatrix( ll, d_dims, ComplexFlag )
           
      !call mxCopyReal8ToPtr( Hout ,mxGetPr( plhs(1) ), ll * 3 )
      
      !deallocate(  M, dims, pos, rotAngles )
         
      end subroutine
      
      !::
      !::Performs the box integral around the region specified 
      !::and over the Maxwell stress tensor using integral2 via quadpack
      !::xPos(2) defines min and max of the x range
      !::yPos(2) defines min and max of the y range
      !::zPos(2) defines min and max of the z range
      !::mdStr is a struct with data from the model run and information to be shared
      !::between integration function calls
      subroutine boxIntegral( xPos, yPos, zPos )
      use model_glob
      USE MagPrism_Call
      USE Parameters_CALL
      real*8,dimension(2),intent(in) :: xPos,yPos,zPos
      
      !::find the rotation matrices once and for all
      call findRotationMatrices( model%rotAngles, model%n, 
     +                            model%rotMat, model%rotMatInv)
      
      end subroutine boxIntegral
      
            
      
      
      