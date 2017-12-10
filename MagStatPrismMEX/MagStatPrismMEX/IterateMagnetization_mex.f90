#include "fintrf.h"  
      

    
    !::
    !::prhs(1) is an array of structs containing the tile information, size n_tiles
    !::prhs(2) is an integer representing the number of tiles, n_tiles
    !::prhs(3) is a double defining the error tolerance for the iteration
    !::plhs(1) is an array of structs containing the updated tiles, size n_tiles
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use IterateMagnetSolution
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
      mwPointer mxGetFieldByNumber,mxCreateStructArray
      
      mwPointer mxGetM, mxGetN
      mwSize mxGetNumberOfDimensions, sx
      mwPointer mxGetDimensions
      
      integer*4 ComplexFlag,classid

!     Pointers to input/output mxArrays:
      mwPointer :: stPr      
      type(MagTile),allocatable,dimension(:) :: cylTile
      mwPointer :: r0Ptr,theta0Ptr,z0Ptr,drPtr,dthetaPtr,dzPtr,MPtr,u_eaPtr,u_oa1Ptr,u_oa2Ptr,mur_eaPtr,mur_oaPtr,MremPtr
      mwPointer,dimension(:,:),allocatable :: pvalue
      integer*4 :: n_tiles
      mwSize,dimension(2) :: dims
      mwSize :: s1,s2
      mwIndex :: i
      real*8 :: err_max
      
      integer*4,parameter :: nfields=13
      character(len=10),dimension(nfields) :: fieldnames
      
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
      
      !::Setup the names of the members of the input struct
      fieldnames(1) = 'r0'
      fieldnames(2) = 'theta0'
      fieldnames(3) = 'z0'
      fieldnames(4) = 'dr'
      fieldnames(5) = 'dtheta'
      fieldnames(6) = 'dz'
      fieldnames(7) = 'M'
      fieldnames(8) = 'u_ea'
      fieldnames(9) = 'u_oa1'
      fieldnames(10) = 'u_oa2'
      fieldnames(11) = 'mu_r_ea'
      fieldnames(12) = 'mu_r_oa'
      fieldnames(13) = 'Mrem'
      
      !::Copy the number of tiles
      call mxCopyPtrToInteger4(mxGetPr(prhs(2)), n_tiles,1)
      
      !::Copy the max error
      call mxCopyPtrToReal8(mxGetPr(prhs(3)), err_max, 1 )
      
      !::allocate the tiles
      allocate( cylTile(n_tiles) )
      !::Copy the input parameters
      sx = 1            
      do i=1,n_tiles
          r0Ptr = mxGetField(prhs(1),i,fieldnames(1))
          theta0Ptr = mxGetField(prhs(1),i,fieldnames(2))
          z0Ptr = mxGetField(prhs(1),i,fieldnames(3))
          drPtr = mxGetField(prhs(1),i,fieldnames(4))
          dthetaPtr = mxGetField(prhs(1),i,fieldnames(5))
          dzPtr = mxGetField(prhs(1),i,fieldnames(6))
          MPtr = mxGetField(prhs(1),i,fieldnames(7))
          u_eaPtr = mxGetField(prhs(1),i,fieldnames(8))
          u_oa1Ptr = mxGetField(prhs(1),i,fieldnames(9))
          u_oa2Ptr = mxGetField(prhs(1),i,fieldnames(10))
          mur_eaPtr = mxGetField(prhs(1),i,fieldnames(11))
          mur_oaPtr = mxGetField(prhs(1),i,fieldnames(12))
          MremPtr = mxGetField(prhs(1),i,fieldnames(13))
      
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
      
      !::do the calculation
      call iterateMagnetization( cylTile, n_tiles, err_max )
      
      !::Load the result back to Matlab
      
      ComplexFlag = 0
      
      dims(1) = 1
      dims(2) = n_tiles
      sx = 2
      !::create the return array of structs      
      plhs(1) = mxCreateStructArray( sx, dims, nfields, fieldnames)
      
      !::pointer for each data entry
      allocate( pvalue(n_tiles,nfields) )
      
      !::Fill in the data
      do i=1,n_tiles
          s1 = 1
          s2 = 1
          pvalue(i,1) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
          call mxCopyReal8ToPtr( cylTile(i)%r0, mxGetPr( pvalue(i,1) ), 1 )
          call mxSetField( plhs(1), i, fieldnames(1), pvalue(i,1) )
          
          pvalue(i,2) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%theta0, mxGetPr( pvalue(i,2) ), 1 )
          call mxSetField( plhs(1), i, fieldnames(2), pvalue(i,2) )
          
          pvalue(i,3) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
          call mxCopyReal8ToPtr( cylTile(i)%z0, mxGetPr( pvalue(i,3) ), 1 )
          call mxSetField( plhs(1), i, fieldnames(3), pvalue(i,3) )
          
          pvalue(i,4) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%dr, mxGetPr( pvalue(i,4) ), 1 )
          call mxSetField( plhs(1), i, fieldnames(4), pvalue(i,4))
          
          pvalue(i,5) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
          call mxCopyReal8ToPtr( cylTile(i)%dtheta, mxGetPr( pvalue(i,5) ), 1 )
          call mxSetField( plhs(1), i, fieldnames(5), pvalue(i,5) )
          
          pvalue(i,6) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
          call mxCopyReal8ToPtr( cylTile(i)%dz, mxGetPr( pvalue(i,6) ), 1 )
          call mxSetField( plhs(1), i, fieldnames(6), pvalue(i,6) )
          
          s2 = 3
          pvalue(i,7) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%M, mxGetPr( pvalue(i,7) ), 3 )
          call mxSetField( plhs(1), i, fieldnames(7), pvalue(i,7))
          
          pvalue(i,8) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%u_ea, mxGetPr( pvalue(i,8) ), 3 )
          call mxSetField( plhs(1), i, fieldnames(8), pvalue(i,8) )
          
          pvalue(i,9) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%u_oa1, mxGetPr( pvalue(i,9) ), 3 )
          call mxSetField( plhs(1), i, fieldnames(9), pvalue(i,9) )
          
          pvalue(i,10) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%u_oa2, mxGetPr( pvalue(i,10) ), 3 )
          call mxSetField( plhs(1), i, fieldnames(10), pvalue(i,10) )

          s2 = 1
          pvalue(i,11) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%mu_r_ea, mxGetPr( pvalue(i,11) ), 1 )
          call mxSetField( plhs(1), i, fieldnames(11), pvalue(i,11) )
          
          pvalue(i,12) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%mu_r_oa, mxGetPr( pvalue(i,12) ), 1 )
          call mxSetField( plhs(1), i, fieldnames(12), pvalue(i,12) )
          
          pvalue(i,13) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%Mrem, mxGetPr( pvalue(i,13) ), 1 )
          call mxSetField( plhs(1), i, fieldnames(13), pvalue(i,13) )
          
      enddo
      
      
      deallocate(pvalue,cylTile)
            
      end subroutine
      
      
      
            
      
      
      