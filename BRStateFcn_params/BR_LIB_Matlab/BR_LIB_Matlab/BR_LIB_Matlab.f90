  
#include "fintrf.h"
      
      

!     Gateway routine
!     The prhs is a pointer array with the following content
!     prhs(1) BRStruct, struct containing the BR parameters
!     prhs(2) n, integer specifying the number of elements in the input arrays
!     prhs(3) T(n) input temperature array
!     prhs(4) H(n) input field array
!     prhs(5) p input pressure
!     prhs(6) Tc(n) input T0 array

!     plhs(1) M(n,2) output magnetization array
!     plhs(2) S(n,2) output entropy array

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      USE BR_SMITH_Call
      
!     Declarations
      implicit none
      
!     mexFunction arguments:
      mwPointer plhs(*), prhs(*), xx
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
      mwPointer kappaPtr,etaPtr,JPtr,gPtr,NsPtr,rhoPtr,signBetaPtr
      
      mwPointer mxGetM, mxGetN
      mwSize mxGetNumberOfDimensions,dimIteOut(1)
      mwPointer mxGetDimensions


      integer :: n,signBeta
      real :: kappa,eta,J,g,Ns,rho,p
      real,dimension(:),allocatable :: T,H,Tc
      real,dimension(:,:),allocatable :: M,S,deltaF
      
      
!      open (11, file='debug.txt',status='replace', access='sequential',form='formatted',action='write' )
      
      !-----------------------------------------------------------------------	  	       
	  open(11,file='version_BL_LIB_Matlab.txt',status='unknown',access='sequential',
     +              form='formatted',position='rewind',action='write')     
	  write(11,*) "BL_LIB_Matlab compiled on:"
      write(11,*) __DATE__
      write(11,*) __TIME__
      close(11)
      !-----------------------------------------------------------------------
      
      !-----------------------------------------------------------------------
!     Check for proper number of arguments. 
      if(nrhs .ne. 6) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nInput',
     +                           '6 inputs required.')
      elseif(nlhs .gt. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nOutput',
     +                           'Too many output arguments.')
      endif

!Check the type of the inputs      
      if ( .NOT. mxIsStruct(prhs(1)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input one should be a struct')
      elseif ( .NOT. mxIsInt32(prhs(2)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input two should be an integer')      
      elseif ( .NOT. mxIsDouble(prhs(3)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input three should be double')
      elseif ( .NOT. mxIsDouble(prhs(4)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input four should be double')
      elseif ( .NOT. mxIsDouble(prhs(5)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input five should be double')
      elseif ( .NOT. mxIsDouble(prhs(6)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input six should be double')      
      endif
     
     
      !::Get the BR parameters from the input struct
      kappaPtr = mxGetField( prhs(1), 1, 'kappa' )
      etaPtr = mxGetField( prhs(1), 1, 'eta' )
      JPtr = mxGetField( prhs(1), 1, 'J' )
      gPtr = mxGetField( prhs(1), 1, 'g' )
      NsPtr = mxGetField( prhs(1), 1, 'Ns' )
      rhoPtr = mxGetField( prhs(1), 1, 'rho' )
      signBetaPtr = mxGetField( prhs(1), 1, 'signBeta' )
      
      call mxCopyPtrToReal8( mxGetPr(kappaPtr), kappa, 1 )
      call mxCopyPtrToReal8( mxGetPr(etaPtr), eta, 1 )
      call mxCopyPtrToReal8( mxGetPr(JPtr), J, 1 )
      call mxCopyPtrToReal8( mxGetPr(gPtr), g, 1 )
      call mxCopyPtrToReal8( mxGetPr(NsPtr), Ns, 1 )
      call mxCopyPtrToReal8( mxGetPr(rhoPtr), rho, 1 )
      call mxCopyPtrToInteger4( mxGetPr(signBetaPtr), signBeta, 1 )
      !write(11,*) 'sign',signBeta
      !::Get the array dimension
      call mxCopyPtrToInteger4( mxGetPr(prhs(2)), n, 1 )
      allocate( T(n), H(n), Tc(n), M(n,2), S(n,2), deltaF(n,2) )
      T = 0
      H = 0
      Tc = 0
      M = 0
      S = 0
      deltaF = 0
      !::Copy the input arrays
      call mxCopyPtrToReal8( mxGetPr(prhs(3)), T, n )
      call mxCopyPtrToReal8( mxGetPr(prhs(4)), H, n )
      call mxCopyPtrToReal8( mxGetPr(prhs(5)), p, 1 )
      call mxCopyPtrToReal8( mxGetPr(prhs(6)), Tc, n )
      
      !write(11,*) minval(T),maxval(T)
      !write(11,*) minval(H),maxval(H)
      !write(11,*) p
      !write(11,*) minval(Tc),maxval(Tc)
      !write(11,*) n
      !close(11)
      
      !::Initializes the BR parameters
      call initializeBRParameters( kappa, eta, J, g, Ns, rho, signBeta )
      
      !::Input: T (n), H(n), p, Tc(n) and n
      !::Output M(n), S(n), deltaF(n) are output
      !::
      
      call getM_BR_batch( T, H, p, Tc, M, S, deltaF, n )
     
      !call mexWarnMsgIdAndTxt('MATLAB:Matlab_single_mex:Debug', 'Calculation done. Returning.')
      
      !open (11, file='debug.txt',status='old', access='sequential',form='formatted',action='write',position='append' )
      !write(11,*) minval(M),maxval(M)
      !write(11,*) minval(S),maxval(S)
      !write(11,*) minval(deltaF),maxval(deltaF)
      !close(11)
      !::Load the result back to Matlab
      
      plhs(1) = mxCreateDoubleMatrix(n,2,0)
      
      plhs(2) = mxCreateDoubleMatrix(n,2,0)
      
      plhs(3) = mxCreateDoubleMatrix(n,2,0)
      
      call mxCopyReal8ToPtr( M, mxGetPr(plhs(1)), n*2 )
      
      call mxCopyReal8ToPtr( S, mxGetPr(plhs(2)), n*2 )
      
      call mxCopyReal8ToPtr( deltaF, mxGetPr(plhs(3)), n*2 )
      
      
      
      end subroutine
      
      
      
      