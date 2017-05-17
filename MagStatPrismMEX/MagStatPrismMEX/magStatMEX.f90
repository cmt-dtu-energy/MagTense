  
#include "fintrf.h"
      
      

!     Gateway routine
!     The prhs is a pointer array with the following content
!     prhs(1) pos(n,3), real
!     prhs(2) dims(n,3), real
!     prhs(3) dir(n,3), real
!     prhs(4) magType(n,3), integer
!     prhs(5) Hext(n,3), real
!     prhs(6) T(n,1), real
!     prhs(7) n, integer
!     prhs(8) stateFcn(m),array of struct with elements T, H and M(length(T),length(H))
!     prhs(9) stateFcnIndices(n), integer
!     prhs(10) m, integer
!     prhs(11) solPts(l,3), real
!     prhs(12) l, integer
!     prhs(13) spacedim, integer
!--------return values--------
!     plhs(1) Hout(n+l,3), real
!     plhs(2) Mout(n,3),real
!( pos, dims, dir, magType, Hext,T, n, stateFcn, stateFcnIndices, m, solPts, l, Hout, Mout )
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
      mwPointer stateFcnPtr,TPtr,HPtr,MCoolPtr,MHeatPtr,nTPtr,nHPtr
      mwPointer :: nCritCoolPtr,TCritCoolPtr,HCritCoolPtr
      mwPointer :: nCritHeatPtr,TCritHeatPtr,HCritHeatPtr
      
      mwPointer mxGetM, mxGetN
      mwSize mxGetNumberOfDimensions,dimIteOut(1)
      mwPointer mxGetDimensions
      
      mwSize nnn, d_dims
      integer*4 ComplexFlag

!     Pointers to input/output mxArrays:
      mwPointer :: stPr
      integer*4 :: nn, mm, ll, spacedim, i, dummy,nIte
      real*8,dimension(:,:),allocatable :: pos,dims,dir,Hext,solPts,Hout
      real*8,dimension(:,:),allocatable :: Mout,rotAngles,H_init
      real*8,dimension(:),allocatable :: T
      integer*4,dimension(:),allocatable :: magType,stateFcnIndices
      integer*4,dimension(:),allocatable :: hyst_map_init, hyst_map
      type(stateFunction),dimension(:),allocatable :: stateFcn
      
      integer*4, external :: debugMLCallback
      
!     !      open (11, file='debug.txt',status='replace',
!     !+  access='sequential',form='formatted',action='write' )
!     !      write(11,*) 'START'
!     !       close(11)
      
      !-----------------------------------------------------------------------	  	       
	     open(11,file='version_magStatMEX.txt',status='unknown',
     +   access='sequential',form='formatted',position='rewind',
     +                       action='write')     
	  write(11,*) "magStatMEX compiled on:"
      write(11,*) __DATE__
      write(11,*) __TIME__
      close(11)
      !-----------------------------------------------------------------------

!     Check for proper number of arguments. 
      if(nrhs .lt. 15) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nInput',
     +                           '15 inputs or more are required.')
      elseif(nlhs .gt. 4) then
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
      elseif ( .NOT. mxIsInt32(prhs(4)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input four should be integer')
      elseif ( .NOT. mxIsDouble(prhs(5)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input five should be double')
      elseif ( .NOT. mxIsDouble(prhs(6)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input six should be double')
      elseif ( .NOT. mxIsInt32(prhs(7)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input seven should be integer')
      elseif ( .NOT. mxIsStruct(prhs(8)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input eight should be a struct')
      elseif ( .NOT. mxIsInt32(prhs(9)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input nine should be integer')
      elseif ( .NOT. mxIsInt32(prhs(10)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input Ten should be integer')
      elseif ( .NOT. mxIsDouble(prhs(11)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input eleven should be double')
      elseif ( .NOT. mxIsInt32(prhs(12)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input twelve should be integer')
      elseif ( .NOT. mxIsInt32(prhs(13)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input thirteen should be integer')     
      elseif ( .NOT. mxIsDouble(prhs(14)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input fourteen should be double')           
      elseif ( .NOT. mxIsInt32(prhs(15)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input fifteen should be integer')     
      endif
      if ( nrhs .eq. 16 ) then
        if ( .NOT. mxIsDouble(prhs(16)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input sixteen should be double')  
      endif
      endif
	 
       !write(11,*) 'HUND'
      call mxCopyPtrToInteger4(mxGetPr(prhs(7)) , nn, 1)
      call mxCopyPtrToInteger4(mxGetPr(prhs(10)), mm, 1)
      call mxCopyPtrToInteger4(mxGetPr(prhs(12)), ll, 1)
      call mxCopyPtrToInteger4(mxGetPr(prhs(13)), spacedim, 1)

      allocate(pos(nn,3),dims(nn,3),dir(nn,3),magType(nn))
      allocate(hext(nn,3),T(nn),stateFcnIndices(nn),solPts(ll,3))
      allocate(rotAngles(nn,3) )
      allocate(stateFcn(2*mm))
      allocate(hyst_map_init(nn),hyst_map(nn))
      
      call mxCopyPtrToInteger4(mxGetPr(prhs(4)) , magType, nn)
      call mxCopyPtrToInteger4(mxGetPr(prhs(9)) , stateFcnIndices,nn)
      call mxCopyPtrToInteger4(mxGetPr(prhs(15)) , hyst_map_init, nn)
      
      call mxCopyPtrToReal8(mxGetPr(prhs(1)) ,pos,nn*3)
      call mxCopyPtrToReal8(mxGetPr(prhs(2)) ,dims,nn*3)
      call mxCopyPtrToReal8(mxGetPr(prhs(3)) ,dir,nn*3)
      call mxCopyPtrToReal8(mxGetPr(prhs(5)) ,Hext,nn*3)
      call mxCopyPtrToReal8(mxGetPr(prhs(6)) ,T,nn)
      call mxCopyPtrToReal8(mxGetPr(prhs(11)),solPts,ll*3)
      call mxCopyPtrToReal8(mxGetPr(prhs(14)),rotAngles,nn*3)
      
      
      
      do i=1,mm          
          TPtr = mxGetField(prhs(8),i,'T')          
          HPtr = mxGetField(prhs(8),i,'H')
          MCoolPtr = mxGetField(prhs(8),i,'MCool')
          MHeatPtr = mxGetField(prhs(8),i,'MHeat')
          nTPtr = mxGetField(prhs(8),i,'nT')
          nHPtr = mxGetField(prhs(8),i,'nH')
          
          nCritCoolPtr = mxGetField(prhs(8),i,'nCritCool')          
          TCritCoolPtr = mxGetField(prhs(8),i,'TCritCool')
          HCritCoolPtr = mxGetField(prhs(8),i,'HCritCool')
          
          nCritHeatPtr = mxGetField(prhs(8),i,'nCritHeat')
          TCritHeatPtr = mxGetField(prhs(8),i,'TCritHeat')
          HCritHeatPtr = mxGetField(prhs(8),i,'HCritHeat')
          
          call mxCopyPtrToInteger4(mxGetPr(nTPtr), stateFcn(i)%nT, 1)                             
          call mxCopyPtrToInteger4(mxGetPr(nHPtr), stateFcn(i)%nH, 1)
          
          call mxCopyPtrToInteger4(mxGetPr(nCritCoolPtr), 
     + stateFcn(i)%nCritCool, 1)
          call mxCopyPtrToInteger4(mxGetPr(nCritHeatPtr), 
     + stateFcn(i)%nCritHeat, 1)
          
          call mxCopyPtrToInteger4(mxGetPr(nTPtr), 
     + stateFcn(i+mm)%nT, 1)
          call mxCopyPtrToInteger4(mxGetPr(nHPtr), 
     + stateFcn(i+mm)%nH, 1)
          
          allocate( stateFcn(i)%T(stateFcn(i)%nT) )
          allocate( stateFcn(i)%H(stateFcn(i)%nH) )
          allocate( stateFcn(i)%M(stateFcn(i)%nT,stateFcn(i)%nH) )
          allocate( stateFcn(i)%Tcrit_cool(stateFcn(i)%nCritCool) )
          allocate( stateFcn(i)%Hcrit_cool(stateFcn(i)%nCritCool) )
          allocate( stateFcn(i)%Tcrit_heat(stateFcn(i)%nCritHeat) )
          allocate( stateFcn(i)%Hcrit_heat(stateFcn(i)%nCritHeat) )
          
          
          allocate( stateFcn(i+mm)%T(stateFcn(i+mm)%nT) )
          allocate( stateFcn(i+mm)%H(stateFcn(i+mm)%nH) )
          allocate( stateFcn(i+mm)%M(stateFcn(i+mm)%nT,
     + stateFcn(i+mm)%nH) )
          
          stateFcn(i)%T = 0
          stateFcn(i)%H = 0
          stateFcn(i)%M = 0
          
          stateFcn(i+mm)%T = 0
          stateFcn(i+mm)%H = 0
          stateFcn(i+mm)%M = 0
          
          
          call mxCopyPtrToReal8(mxGetPr(TPtr),stateFcn(i)%T,
     + stateFcn(i)%nT)
          call mxCopyPtrToReal8(mxGetPr(HPtr),
     + stateFcn(i)%H,stateFcn(i)%nH)
          call mxCopyPtrToReal8(mxGetPr(MCoolPtr),
     + stateFcn(i)%M,stateFcn(i)%nT*stateFcn(i)%nH)
     
     
          call mxCopyPtrToReal8(mxGetPr(TPtr),stateFcn(i+mm)%T,
     + stateFcn(i+mm)%nT)
          call mxCopyPtrToReal8(mxGetPr(HPtr),
     + stateFcn(i+mm)%H,stateFcn(i+mm)%nH)
          call mxCopyPtrToReal8(mxGetPr(MHeatPtr),
     + stateFcn(i+mm)%M,stateFcn(i+mm)%nT*stateFcn(i+mm)%nH)
     
     
          call mxCopyPtrToReal8(mxGetPr(TCritCoolPtr),
     + stateFcn(i)%Tcrit_cool,stateFcn(i)%nCritCool )
     
          call mxCopyPtrToReal8(mxGetPr(TCritHeatPtr),
     + stateFcn(i)%Tcrit_heat,stateFcn(i)%nCritHeat )
     
          call mxCopyPtrToReal8(mxGetPr(HCritCoolPtr),
     + stateFcn(i)%Hcrit_cool,stateFcn(i)%nCritCool )
     
          call mxCopyPtrToReal8(mxGetPr(HCritHeatPtr),
     + stateFcn(i)%Hcrit_heat,stateFcn(i)%nCritHeat )
          
      enddo
      
      if ( nrhs .eq. 16 ) then
          allocate( H_init(nn,3) )
          call mxCopyPtrToReal8(mxGetPr(prhs(16)) ,H_init,nn*3)
          !::call mexWarnMsgIdAndTxt('MATLAB:Matlab_single_mex:Debug', 'Doing calculation')
      call calcH( pos, dims, dir, magType, Hext,T, nn, stateFcn,
     +            stateFcnIndices, mm, solPts, ll, 
     +            spacedim, Hout, Mout, nIte, rotAngles,
     +            hyst_map_init, hyst_map, H_init )    
      else
      !::call mexWarnMsgIdAndTxt('MATLAB:Matlab_single_mex:Debug', 'Doing calculation')
      call calcH( pos, dims, dir, magType, Hext,T, nn, stateFcn,
     +            stateFcnIndices, mm, solPts, ll, 
     +            spacedim, Hout, Mout, nIte, rotAngles,
     +            hyst_map_init, hyst_map )    
      endif
      
      
      
      !call mexWarnMsgIdAndTxt('MATLAB:Matlab_single_mex:Debug', 'Calculation done. Returning.')
      
      !::Load the result back to Matlab
      !:: Ensure that the variables used to call mxCreateDoubleMatrix are the right type
      nnn = nn+ll;
      
      !!call calcH( pos, dims, dir, magType, Hext,T, nn, stateFcn, stateFcnIndices, mm, solPts, ll, 
    !!+            spacedim, Hout, Mout, nIte, rotAngles )
      
      d_dims = 3;
      ComplexFlag = 0;
      plhs(1) = mxCreateDoubleMatrix(nnn,d_dims,ComplexFlag)
      plhs(2) = mxCreateDoubleMatrix(nn,d_dims,ComplexFlag)
      dimIteOut(1) = 1
      plhs(3) = mxCreateNumericArray(1,
     + dimIteOut,mxClassIDFromClassName('int32'),0)
       
       dimIteOut(1) = nn
       plhs(4) = mxCreateNumericArray(1,dimIteOut,
     + mxClassIDFromClassName('int32'),ComplexFlag)
     
     
      call mxCopyReal8ToPtr(Hout,mxGetPr(plhs(1)),(nn+ll)*3)
      call mxCopyReal8ToPtr(Mout,mxGetPr(plhs(2)),nn*3)
      call mxCopyInteger8ToPtr(nIte,mxGetPr(plhs(3)),1)
      call mxCopyInteger4ToPtr(hyst_map,mxGetPr(plhs(4)),nn)
         
      end subroutine
      
      function debugMLCallback( message )
      integer*4, external :: mexPrintf
      character(len=80) :: message
      integer*4 :: k
      integer*4 :: debugMLCallback
      
       k = mexPrintf(message)
       
       debugMLCallback = 0
      
      end function debugMLCallback
      
      !function 
      !integer*4 :: mexCallMATLAB(nlhs, plhs, nrhs, prhs, functionName)
      !integer*4 :: nlhs, nrhs
      !mwPointer :: plhs(1), prhs(1)
      !character*(*) :: functionName
      !character(*) :: message
      !mwPointer,parameter :: NULL = 0
      
      
      !mexCallMATLAB( 0, NULL, 1, prhs, 'mexDebugMessag' )
      

      
      