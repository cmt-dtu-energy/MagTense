  
#include "fintrf.h"
      
      

!     Gateway routine
!     The prhs is a pointer array with the following content
!     prhs(1) pos(n,3), real
!     prhs(2) dims(n,3), real
!     prhs(3) dir(n,3), real
!     prhs(prhs_in(4)) magType(n,3), integer
!     prhs(prhs_in(6)) Hext(n,3), real
!     prhs(6) T(n,1), real
!     prhs(prhs_in(8)) n, integer
!     prhs(8) stateFcn(m),array of struct with elements T, H and M(length(T),length(H))
!     prhs(prhs_in(10)) stateFcnIndices(n), integer
!     prhs(prhs_in(11)) m, integer
!     prhs(prhs_in(12)) solPts(l,3), real
!     prhs(prhs_in(13)) l, integer
!     prhs(prhs_in(14)) spacedim, integer
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
      integer(kind=4) nlhs, nrhs

!     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      mwPointer mxCreateNumericArray
      integer(kind=4) mxIsNumeric
      integer(kind=4) mxIsDouble
      integer(kind=4) mxIsInt32
      integer(kind=4) mxIsStruct
      integer(kind=4) mxIsCell
      integer(kind=4) mxClassIDFromClassName
      integer(kind=4) writeToConsole
      mwPointer mxGetField
      mwPointer mxGetFieldByNumber
      mwPointer mxGetFieldNumber
      mwPointer :: stateFcnPtr,TPtr,HPtr,MCoolPtr,MHeatPtr,nTPtr,nHPtr
      mwPointer :: nCritCoolPtr,TCritCoolPtr,HCritCoolPtr
      mwPointer :: nCritHeatPtr,TCritHeatPtr,HCritHeatPtr
      
      mwPointer mxGetM, mxGetN
      mwSize mxGetNumberOfDimensions,dimIteOut(1),sx
      mwPointer mxGetDimensions
      
      mwSize nnn, d_dims
      integer*4 ComplexFlag

!     Pointers to input/output mxArrays:
      mwPointer :: stPr
      mwIndex :: i
      integer*4 :: nn, mm, ll, spacedim, dummy,nIte
      integer*4 :: fieldnum, bool
      real*8,dimension(:,:),allocatable :: pos,dims,dir,Hext,solPts,Hout
      real*8,dimension(:,:),allocatable :: Mout,rotAngles,H_init
      real*8,dimension(:),allocatable :: T
      integer*4,dimension(:),allocatable :: magType,stateFcnIndices,formType
      integer*4,dimension(:),allocatable :: hyst_map_init, hyst_map
      type(stateFunction),dimension(:),allocatable :: stateFcn
      
      integer*4, external :: debugMLCallback
      
      integer*4,dimension(16) :: prhs_in
      
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

      
      prhs_in(1)  = 1        !--- pos
      prhs_in(2)  = 2        !--- dims
      prhs_in(3)  = 3        !--- dir
      prhs_in(4)  = 4        !--- magtype
      prhs_in(5)  = 5        !--- formType
      prhs_in(6)  = 6        !--- Hext
      prhs_in(7)  = 7        !--- T
      prhs_in(8)  = 8        !--- n
      prhs_in(9)  = 9        !--- stateFcn
      prhs_in(10) = 10       !--- stateFcnIndices
      prhs_in(11) = 11       !--- m
      prhs_in(12) = 12       !--- solPts
      prhs_in(13) = 13       !--- l
      prhs_in(14) = 14       !--- spacedim
      prhs_in(15) = 15       !--- rotAngles
      prhs_in(16) = 16       !--- hyst_map_init
      prhs_in(17) = 17       !--- H_init
      
      
!     Check for proper number of arguments. 
      if(nrhs .lt. 16) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nInput',
     +                           '15 inputs or more are required.')
      elseif(nlhs .gt. 4) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nOutput',
     +                           'Too many output arguments.')
      endif

!Check the type of the inputs      
      if ( .NOT. mxIsDouble(prhs(prhs_in(1))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input one should be double')
      elseif ( .NOT. mxIsDouble(prhs(prhs_in(2))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input two should be double')      
      elseif ( .NOT. mxIsDouble(prhs(prhs_in(3))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input three should be double')
      elseif ( .NOT. mxIsInt32(prhs(prhs_in(4))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input four should be integer')
      elseif ( .NOT. mxIsInt32(prhs(prhs_in(5))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input five should be integer')    
      elseif ( .NOT. mxIsDouble(prhs(prhs_in(6))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input six should be double')
      elseif ( .NOT. mxIsDouble(prhs(prhs_in(7))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input seven should be double')
      elseif ( .NOT. mxIsInt32(prhs(prhs_in(8))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input eight should be integer')
      elseif ( .NOT. mxIsStruct(prhs(prhs_in(9))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input nine should be a struct')
      elseif ( .NOT. mxIsInt32(prhs(prhs_in(10))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input ten should be integer')
      elseif ( .NOT. mxIsInt32(prhs(prhs_in(11))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input eleven should be integer')
      elseif ( .NOT. mxIsDouble(prhs(prhs_in(12))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input twelve should be double')
      elseif ( .NOT. mxIsInt32(prhs(prhs_in(13))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input thirteen should be integer')
      elseif ( .NOT. mxIsInt32(prhs(prhs_in(14))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input fourteen should be integer')     
      elseif ( .NOT. mxIsDouble(prhs(prhs_in(15))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input fifteen should be double')           
      elseif ( .NOT. mxIsInt32(prhs(prhs_in(16))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input sixteen should be integer')      
      endif
      if ( nrhs .eq. 16 ) then
        if ( .NOT. mxIsDouble(prhs(prhs_in(17))) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input seventen should be double')  
      endif
      endif
	        
      call mxCopyPtrToInteger4(mxGetPr(prhs(prhs_in(8))) , nn, 1)
      call mxCopyPtrToInteger4(mxGetPr(prhs(prhs_in(11))), mm, 1)
      call mxCopyPtrToInteger4(mxGetPr(prhs(prhs_in(13))), ll, 1)
      call mxCopyPtrToInteger4(mxGetPr(prhs(prhs_in(14))), spacedim, 1)

      allocate(pos(nn,3),solPts(ll,3))
      allocate(dims(nn,3),dir(nn,3),magType(nn),formType(nn))
      allocate(Hext(nn+ll,3),T(nn),stateFcnIndices(nn))
      allocate(rotAngles(nn,3) )
      allocate(stateFcn(2*mm))
      allocate(hyst_map_init(nn),hyst_map(nn))
      
      sx = nn * 3
      call mxCopyPtrToReal8(mxGetPr(prhs(prhs_in(1))) ,pos,sx)
      call mxCopyPtrToReal8(mxGetPr(prhs(prhs_in(2))) ,dims,sx)
      call mxCopyPtrToReal8(mxGetPr(prhs(prhs_in(3))) ,dir,sx)
      sx = (nn + ll) * 3
      call mxCopyPtrToReal8(mxGetPr(prhs(prhs_in(6))) ,Hext,sx)
      sx = nn
      call mxCopyPtrToReal8(mxGetPr(prhs(prhs_in(7))) ,T,sx)
      sx = ll * 3
      call mxCopyPtrToReal8(mxGetPr(prhs(prhs_in(12))),solPts,sx)
      sx = nn * 3
      call mxCopyPtrToReal8(mxGetPr(prhs(prhs_in(15))),rotAngles,sx)
      sx = nn
      call mxCopyPtrToInteger4(mxGetPr(prhs(prhs_in(4))) , magType, sx)
      call mxCopyPtrToInteger4(mxGetPr(prhs(prhs_in(5))) , formType, sx)
      call mxCopyPtrToInteger4(mxGetPr(prhs(prhs_in(10))) , stateFcnIndices,sx)
      call mxCopyPtrToInteger4(mxGetPr(prhs(prhs_in(16))) , hyst_map_init, sx)
      
      
      
      
      do i=1,mm  
          
                  
          nTPtr = mxGetField(prhs(prhs_in(9)), i, 'nT')
          nHPtr = mxGetField(prhs(prhs_in(9)),i,'nH')
          MCoolPtr = mxGetField(prhs(prhs_in(9)),i,'MCool')
          MHeatPtr = mxGetField(prhs(prhs_in(9)),i,'MHeat')
          TPtr = mxGetField(prhs(prhs_in(9)),i,'T')          
          HPtr = mxGetField(prhs(prhs_in(9)),i,'H')
          
          nCritCoolPtr = mxGetField(prhs(prhs_in(9)),i,'nCritCool')          
          TCritCoolPtr = mxGetField(prhs(prhs_in(9)),i,'TCritCool')
          HCritCoolPtr = mxGetField(prhs(prhs_in(9)),i,'HCritCool')
          
          nCritHeatPtr = mxGetField(prhs(prhs_in(9)),i,'nCritHeat')
          TCritHeatPtr = mxGetField(prhs(prhs_in(9)),i,'TCritHeat')
          HCritHeatPtr = mxGetField(prhs(prhs_in(9)),i,'HCritHeat')
          
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
          
          sx = stateFcn(i)%nT
          call mxCopyPtrToReal8(mxGetPr(TPtr),stateFcn(i)%T,sx)
          sx = stateFcn(i)%nH
          call mxCopyPtrToReal8(mxGetPr(HPtr),stateFcn(i)%H,sx)
          sx = stateFcn(i)%nT*stateFcn(i)%nH
          call mxCopyPtrToReal8(mxGetPr(MCoolPtr),stateFcn(i)%M,sx)
     
          sx = stateFcn(i+mm)%nT
          call mxCopyPtrToReal8(mxGetPr(TPtr),stateFcn(i+mm)%T,sx)
          sx = stateFcn(i+mm)%nH
          call mxCopyPtrToReal8(mxGetPr(HPtr),stateFcn(i+mm)%H,sx)
          sx = stateFcn(i+mm)%nT*stateFcn(i+mm)%nH
          call mxCopyPtrToReal8(mxGetPr(MHeatPtr),stateFcn(i+mm)%M,sx)
     
          sx = stateFcn(i)%nCritCool
          call mxCopyPtrToReal8(mxGetPr(TCritCoolPtr),
     + stateFcn(i)%Tcrit_cool,sx )
          sx = stateFcn(i)%nCritHeat
          call mxCopyPtrToReal8(mxGetPr(TCritHeatPtr),
     + stateFcn(i)%Tcrit_heat,sx )
          sx = stateFcn(i)%nCritCool
          call mxCopyPtrToReal8(mxGetPr(HCritCoolPtr),
     + stateFcn(i)%Hcrit_cool,sx )
          sx = stateFcn(i)%nCritHeat
          call mxCopyPtrToReal8(mxGetPr(HCritHeatPtr),
     + stateFcn(i)%Hcrit_heat,sx )
          
      enddo
      
      if ( nrhs .eq. 16 ) then
          allocate( H_init(nn,3) )
          sx = nn * 3
          call mxCopyPtrToReal8(mxGetPr(prhs(prhs_in(17))) ,H_init,sx)
          !::call mexWarnMsgIdAndTxt('MATLAB:Matlab_single_mex:Debug', 'Doing calculation')
      call calcH( pos, dims, dir, magType, formType, Hext, T, nn, 
     +            stateFcn, stateFcnIndices, mm, solPts, ll, 
     +            spacedim, Hout, Mout, nIte, rotAngles,
     +            hyst_map_init, hyst_map, H_init )    
      else
      !::call mexWarnMsgIdAndTxt('MATLAB:Matlab_single_mex:Debug', 'Doing calculation')
      call calcH( pos, dims, dir, magType, formType, Hext,T, nn, 
     +            stateFcn, stateFcnIndices, mm, solPts, ll, 
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
      sx = nnn
      plhs(1) = mxCreateDoubleMatrix(sx,d_dims,ComplexFlag)
      sx = nn
      plhs(2) = mxCreateDoubleMatrix(sx,d_dims,ComplexFlag)
      dimIteOut(1) = 1
      plhs(3) = mxCreateNumericArray(1,
     + dimIteOut,mxClassIDFromClassName('int32'),0)
       
       dimIteOut(1) = nn
       plhs(4) = mxCreateNumericArray(1,dimIteOut,
     + mxClassIDFromClassName('int32'),ComplexFlag)
     
      sx = (nn+ll)*3
      call mxCopyReal8ToPtr(Hout,mxGetPr(plhs(1)),sx)
      sx = nn * 3
      call mxCopyReal8ToPtr(Mout,mxGetPr(plhs(2)),sx)
      call mxCopyInteger8ToPtr(nIte,mxGetPr(plhs(3)),1)
      sx = nn
      call mxCopyInteger4ToPtr(hyst_map,mxGetPr(plhs(4)),sx)
         
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
      !mwPointer :: plhs(1), prhs(prhs_in(1))
      !character*(*) :: functionName
      !character(*) :: message
      !mwPointer,parameter :: NULL = 0
      
      
      !mexCallMATLAB( 0, NULL, 1, prhs, 'mexDebugMessag' )
      

      
      