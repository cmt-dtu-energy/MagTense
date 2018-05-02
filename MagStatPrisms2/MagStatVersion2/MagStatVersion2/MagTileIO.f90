#include "fintrf.h"  
module MagTileIO

use TileNComponents
use MagStatParameters
implicit none
    
    contains
    
    subroutine displayProgress( prog )
    integer,intent(in) :: prog
    
    integer :: mexCallMATLAB, nlhs_cb, nrhs_cb, tmp
    mwPointer plhs_cb(1), prhs_cb(1),mxCreateString
    character*(4) functionName_cb,prog_str   
    character(len=8) :: fmt

    fmt = '(I4.2)'
    write (prog_str,fmt) prog
    
    
    functionName_cb = "disp"
    nlhs_cb = 0
    nrhs_cb = 1
      
    prhs_cb(1) = mxCreateString(prog_str)
      
    tmp = mexCallMATLAB(nlhs_cb, plhs_cb, nrhs_cb, prhs_cb, "disp")
    
    end subroutine displayProgress
    
    subroutine displayIteration( err, err_max )
    real,intent(in) :: err,err_max
    
    integer :: mexCallMATLAB, nlhs_cb, nrhs_cb, tmp
    mwPointer plhs_cb(1), prhs_cb(1),mxCreateString
    character*(4) :: functionName_cb
    character*(100) :: prog_str   
    character(len=8) :: fmt

    fmt = '(f15.7)'
    write (prog_str,*) "err: ",err,"Max. err: ",err_max
    
    
    functionName_cb = "disp"
    nlhs_cb = 0
    nrhs_cb = 1
      
    prhs_cb(1) = mxCreateString(prog_str)
      
    tmp = mexCallMATLAB(nlhs_cb, plhs_cb, nrhs_cb, prhs_cb, "disp")
    
    end subroutine displayIteration
    
    !::
    !::Loads the data from Matlab into a fortran struct
    !::
    subroutine loadMagTile( prhs, cylTile, n_tiles )
    mwPointer, intent(in) :: prhs
    character(len=10),dimension(:),allocatable :: fieldnames    
    integer,intent(in) :: n_tiles    
    type(MagTile),intent(inout),dimension(n_tiles) :: cylTile
    
    mwIndex :: i,j,k,l
    mwSize :: sx
    integer :: nfields,cnt,n
    real*8,dimension(3) :: rectDims
    mwPointer :: r0Ptr,theta0Ptr,z0Ptr,drPtr,dthetaPtr,dzPtr,MPtr,u_eaPtr,u_oa1Ptr,u_oa2Ptr,mur_eaPtr,mur_oaPtr,MremPtr
    mwPointer :: tileTypePtr,offsetPtr,rotAnglesPtr,rectDimsPtr,magnetTypePtr,stateFunctionIndexPtr
    mwPointer :: mxGetField, mxGetPr
    real,dimension(:),allocatable :: r, theta, z
    real :: dr,dtheta,dz
    
    
        call getTileFieldnames( fieldnames, nfields )
        do i=1,n_tiles
            r0Ptr = mxGetField(prhs,i,fieldnames(1))
            theta0Ptr = mxGetField(prhs,i,fieldnames(2))
            z0Ptr = mxGetField(prhs,i,fieldnames(3))
            drPtr = mxGetField(prhs,i,fieldnames(4))
            dthetaPtr = mxGetField(prhs,i,fieldnames(5))
            dzPtr = mxGetField(prhs,i,fieldnames(6))
            MPtr = mxGetField(prhs,i,fieldnames(7))
            u_eaPtr = mxGetField(prhs,i,fieldnames(8))
            u_oa1Ptr = mxGetField(prhs,i,fieldnames(9))
            u_oa2Ptr = mxGetField(prhs,i,fieldnames(10))
            mur_eaPtr = mxGetField(prhs,i,fieldnames(11))
            mur_oaPtr = mxGetField(prhs,i,fieldnames(12))
            MremPtr = mxGetField(prhs,i,fieldnames(13))
            tileTypePtr = mxGetField(prhs,i,fieldnames(14))
            offsetPtr = mxGetField(prhs,i,fieldnames(15))
            rotAnglesPtr = mxGetField(prhs,i,fieldnames(16))
            rectDimsPtr =  mxGetField(prhs,i,fieldnames(17))
            magnetTypePtr =  mxGetField(prhs,i,fieldnames(18))
            stateFunctionIndexPtr =  mxGetField(prhs,i,fieldnames(19))
      
            sx = 1
            call mxCopyPtrToReal8(mxGetPr(r0Ptr), cylTile(i)%r0, sx )
            call mxCopyPtrToReal8(mxGetPr(theta0Ptr), cylTile(i)%theta0, sx )
            call mxCopyPtrToReal8(mxGetPr(z0Ptr), cylTile(i)%z0, sx )
            call mxCopyPtrToReal8(mxGetPr(drPtr), cylTile(i)%dr, sx )
            call mxCopyPtrToReal8(mxGetPr(dthetaPtr), cylTile(i)%dtheta, sx )
            call mxCopyPtrToReal8(mxGetPr(dzPtr), cylTile(i)%dz, sx )
            sx = 3
            call mxCopyPtrToReal8(mxGetPr(MPtr), cylTile(i)%M, sx )
            call mxCopyPtrToReal8(mxGetPr(u_eaPtr), cylTile(i)%u_ea, sx )
            call mxCopyPtrToReal8(mxGetPr(u_oa1Ptr), cylTile(i)%u_oa1, sx )
            call mxCopyPtrToReal8(mxGetPr(u_oa2Ptr), cylTile(i)%u_oa2, sx )
            sx = 1
            call mxCopyPtrToReal8(mxGetPr(mur_eaPtr), cylTile(i)%mu_r_ea, sx )
            call mxCopyPtrToReal8(mxGetPr(mur_oaPtr), cylTile(i)%mu_r_oa, sx )
            call mxCopyPtrToReal8(mxGetPr(MremPtr), cylTile(i)%Mrem, sx )
            call mxCopyPtrToInteger4(mxGetPr(tileTypePtr), cylTile(i)%tileType, sx )
            sx = 3
            call mxCopyPtrToReal8(mxGetPr(offsetPtr), cylTile(i)%offset, sx )
            call mxCopyPtrToReal8(mxGetPr(rotAnglesPtr), cylTile(i)%rotAngles, sx )
            call mxCopyPtrToReal8(mxGetPr(rectDimsPtr), rectDims, sx )
            cylTile(i)%a = rectDims(1)
            cylTile(i)%b = rectDims(2)
            cylTile(i)%c = rectDims(3)
            sx = 1
            call mxCopyPtrToInteger4(mxGetPr(magnetTypePtr), cylTile(i)%magnetType, sx )            
            call mxCopyPtrToInteger4(mxGetPr(stateFunctionIndexPtr), cylTile(i)%stateFunctionIndex, sx )
            
            !cylTile(i)%fieldEvaluation = fieldEvaluationCentre
            if ( cyltile(i)%tiletype == tiletypecylpiece ) then
                cyltile(i)%fieldevaluation  = fieldevaluationaverage
                cyltile(i)%n_ave(1) = 5
                cyltile(i)%n_ave(2) = 5
                cyltile(i)%n_ave(3) = 5
                n = cyltile(i)%n_ave(1)*cyltile(i)%n_ave(2)*cyltile(i)%n_ave(3)
                allocate(cyltile(i)%H_ave_pts(n,3),cyltile(i)%H_ave(n,3),r(n),theta(n),z(n))                
                dr = cyltile(i)%dr / cyltile(i)%n_ave(1)
                dtheta = cyltile(i)%dtheta / cyltile(i)%n_ave(2)
                dz = cyltile(i)%dz / cyltile(i)%n_ave(3)
                
                !::set as default value so other users of the code don't have to update.
                cylTile(i)%isIterating = .false.
                
                cnt = 1
                do j=1,cyltile(i)%n_ave(1)
                    do k=1,cyltile(i)%n_ave(2)
                        do l=1,cyltile(i)%n_ave(3)
                            r(cnt) = dr/2 + (j-1)*dr + cylTile(i)%r0 - cylTile(i)%dr / 2
                            theta(cnt) = dtheta/2 + (k-1)*dtheta + cylTile(i)%theta0 - cylTile(i)%dtheta / 2
                            z(cnt) = dz/2 + (l-1)*dz + cylTile(i)%z0 - cylTile(i)%dz / 2
                            cnt = cnt + 1
                        enddo
                    enddo
                enddo
                cyltile(i)%h_ave_pts(:,1) = r * cos( theta ) + cyltile(i)%offset(1)
                cyltile(i)%h_ave_pts(:,2) = r * sin( theta ) + cyltile(i)%offset(2)
                cyltile(i)%h_ave_pts(:,3) = z + cyltile(i)%offset(3)
                cylTile(i)%H_ave(:,:) = 0.
                deallocate(r,theta,z)
            endif
            
            
            
        enddo        
        deallocate(fieldnames)
    end subroutine loadMagTile
    
    !::
    !::Returns the data to matlab
    !::
    subroutine returnMagTile( cylTile, n_tiles, plhs )
    type(MagTile),intent(in),dimension(n_tiles) :: cylTile
    integer, intent(in) :: n_tiles
    mwPointer,intent(inout) :: plhs
    
    integer :: ComplexFlag,classid,mxClassIDFromClassName
    mwSize,dimension(2) :: dims
    mwSize :: s1,s2,sx
    mwPointer,dimension(:,:),allocatable :: pvalue
    mwPointer :: mxCreateStructArray, mxCreateDoubleMatrix,mxGetPr,mxCreateNumericMatrix    
    mwIndex :: i
    
    real,dimension(3) :: rectDims
    character(len=10),dimension(:),allocatable :: fieldnames    
    integer :: nfields 
    
    
    
      call getTileFieldnames( fieldnames, nfields )
    
      !::Load the result back to Matlab      
      ComplexFlag = 0
      
      dims(1) = 1
      dims(2) = n_tiles
      sx = 2
      !::create the return array of structs      
      plhs = mxCreateStructArray( sx, dims, nfields, fieldnames)
      
      !::pointer for each data entry
      allocate( pvalue(n_tiles,nfields) )
      
      !::Fill in the data
      do i=1,n_tiles
          s1 = 1
          s2 = 1
          sx = 1
          pvalue(i,1) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
          call mxCopyReal8ToPtr( cylTile(i)%r0, mxGetPr( pvalue(i,1) ), sx )
          call mxSetField( plhs, i, fieldnames(1), pvalue(i,1) )
          
          pvalue(i,2) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%theta0, mxGetPr( pvalue(i,2) ), sx )
          call mxSetField( plhs, i, fieldnames(2), pvalue(i,2) )
          
          pvalue(i,3) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
          call mxCopyReal8ToPtr( cylTile(i)%z0, mxGetPr( pvalue(i,3) ), sx )
          call mxSetField( plhs, i, fieldnames(3), pvalue(i,3) )
          
          pvalue(i,4) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%dr, mxGetPr( pvalue(i,4) ), sx )
          call mxSetField( plhs, i, fieldnames(4), pvalue(i,4))
          
          pvalue(i,5) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
          call mxCopyReal8ToPtr( cylTile(i)%dtheta, mxGetPr( pvalue(i,5) ), sx )
          call mxSetField( plhs, i, fieldnames(5), pvalue(i,5) )
          
          pvalue(i,6) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
          call mxCopyReal8ToPtr( cylTile(i)%dz, mxGetPr( pvalue(i,6) ), sx )
          call mxSetField( plhs, i, fieldnames(6), pvalue(i,6) )
          
          s2 = 3
          sx = 3
          pvalue(i,7) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%M, mxGetPr( pvalue(i,7) ), sx )
          call mxSetField( plhs, i, fieldnames(7), pvalue(i,7))
          
          pvalue(i,8) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%u_ea, mxGetPr( pvalue(i,8) ), sx )
          call mxSetField( plhs, i, fieldnames(8), pvalue(i,8) )
          
          pvalue(i,9) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%u_oa1, mxGetPr( pvalue(i,9) ), sx )
          call mxSetField( plhs, i, fieldnames(9), pvalue(i,9) )
          
          pvalue(i,10) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%u_oa2, mxGetPr( pvalue(i,10) ), sx )
          call mxSetField( plhs, i, fieldnames(10), pvalue(i,10) )

          s2 = 1
          sx = 1
          pvalue(i,11) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%mu_r_ea, mxGetPr( pvalue(i,11) ), sx )
          call mxSetField( plhs, i, fieldnames(11), pvalue(i,11) )
          
          pvalue(i,12) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%mu_r_oa, mxGetPr( pvalue(i,12) ), sx )
          call mxSetField( plhs, i, fieldnames(12), pvalue(i,12) )
          
          pvalue(i,13) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%Mrem, mxGetPr( pvalue(i,13) ), sx )
          call mxSetField( plhs, i, fieldnames(13), pvalue(i,13) )
          
          classid = mxClassIDFromClassName('int32')
          pvalue(i,14) = mxCreateNumericMatrix(s1,s2,classid,ComplexFlag)
          call mxCopyInteger4ToPtr( cylTile(i)%tileType, mxGetPr( pvalue(i,14) ), sx )
          call mxSetField( plhs, i, fieldnames(14), pvalue(i,14) )
          
          s2 = 3
          sx = 3
          pvalue(i,15) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%offset, mxGetPr( pvalue(i,15) ), sx )
          call mxSetField( plhs, i, fieldnames(15), pvalue(i,15) )
          
          pvalue(i,16) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          call mxCopyReal8ToPtr( cylTile(i)%rotAngles, mxGetPr( pvalue(i,16) ), sx )
          call mxSetField( plhs, i, fieldnames(16), pvalue(i,16) )
          
          pvalue(i,17) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
          rectDims(1) = cylTile(i)%a
          rectDims(2) = cylTile(i)%b
          rectDims(3) = cylTile(i)%c
          call mxCopyReal8ToPtr( rectDims, mxGetPr( pvalue(i,17) ), sx )
          call mxSetField( plhs, i, fieldnames(17), pvalue(i,17) )
          
          s2 = 1
          sx = 1
          classid = mxClassIDFromClassName('int32')
          pvalue(i,18) = mxCreateNumericMatrix(s1,s2,classid,ComplexFlag)
          call mxCopyInteger4ToPtr( cylTile(i)%magnetType, mxGetPr( pvalue(i,18) ), sx )
          call mxSetField( plhs, i, fieldnames(18), pvalue(i,18) )
          
          classid = mxClassIDFromClassName('int32')
          pvalue(i,19) = mxCreateNumericMatrix(s1,s2,classid,ComplexFlag)
          call mxCopyInteger4ToPtr( cylTile(i)%stateFunctionIndex, mxGetPr( pvalue(i,19) ), sx )
          call mxSetField( plhs, i, fieldnames(19), pvalue(i,19) )
          
      enddo
      
      deallocate(pvalue,fieldnames)
    
    end subroutine returnMagTile
    
    !::
    !::Returns an array with the names of the fields expected in the tile struct
    subroutine getTileFieldnames( fieldnames, nfields)
    integer,intent(out) :: nfields
    integer,parameter :: nf=19
    character(len=10),dimension(:),intent(out),allocatable :: fieldnames
            
        nfields = nf
        allocate(fieldnames(nfields))
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
        fieldnames(14) = 'tileType'
        fieldnames(15) = 'offset'
        fieldnames(16) = 'rotAngles'
        fieldnames(17) = 'abc'
        fieldnames(18) = 'magnetType'
        fieldnames(19) = 'stfcnIndex'
        
    
        
    end subroutine getTileFieldnames
    
    !::
    !::Loads the statefunction struct from Matlab to Fortran
    !::
    subroutine loadMagStateFunction( prhs, stateFunction, n )
    mwPointer, intent(in) :: prhs
    type(MagStatStateFunction),dimension(n),intent(inout) :: stateFunction
    integer,intent(in) :: n
    
    mwPointer :: mxGetField, mxGetPr
    
    mwIndex :: i
    mwSize :: sx
    mwPointer :: TPtr,HPtr,MPtr,nTPtr,nHPtr
    character(len=10),dimension(:),allocatable :: fieldnames    
    integer :: nfields
    
    call getStateFunctionFieldnames( fieldnames, nfields )
    
    
    do i=1,n    
        TPtr = mxGetField(prhs,i,fieldnames(1))
        HPtr = mxGetField(prhs,i,fieldnames(2))
        MPtr = mxGetField(prhs,i,fieldnames(3))
        nTPtr = mxGetField(prhs,i,fieldnames(4))
        nHPtr = mxGetField(prhs,i,fieldnames(5))
              
        sx = 1
        call mxCopyPtrToInteger4(mxGetPr(nTPtr), stateFunction(i)%nT, sx )                    
        call mxCopyPtrToInteger4(mxGetPr(nHPtr), stateFunction(i)%nH, sx )            
        
        allocate( stateFunction(i)%T(stateFunction(i)%nT) )
        allocate( stateFunction(i)%H(stateFunction(i)%nH) )
        allocate( stateFunction(i)%M(stateFunction(i)%nT,stateFunction(i)%nH) )
        
        sx = stateFunction(i)%nT
        call mxCopyPtrToReal8(mxGetPr(TPtr), stateFunction(i)%T, sx )
        sx = stateFunction(i)%nH
        call mxCopyPtrToReal8(mxGetPr(HPtr), stateFunction(i)%H, sx )
        
        sx = stateFunction(i)%nT*stateFunction(i)%nH
        call mxCopyPtrToReal8(mxGetPr(MPtr), stateFunction(i)%M, sx )
        
        
    enddo
    
    deallocate( fieldnames )
    
    end subroutine loadMagStateFunction
    
    !::Returns an array with the field names expected in the state function struct
    subroutine getStateFunctionFieldnames( fieldnames, nfields )
    integer,intent(out) :: nfields
    integer,parameter :: nf=5
    character(len=10),dimension(:),intent(out),allocatable :: fieldnames
    
    nfields = nf
    
    allocate( fieldnames(nfields) )
    
    fieldnames(1) = 'T'
    fieldnames(2) = 'H'
    fieldnames(3) = 'M'
    fieldnames(4) = 'nT'
    fieldnames(5) = 'nH'
    
    end subroutine getStateFunctionFieldnames
    
end module MagTileIO
    