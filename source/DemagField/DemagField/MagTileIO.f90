#include "fintrf.h"
    module MagTileIO

    use TileNComponents
    use MagParameters
    implicit none
    
    contains
    
    !!@todo This file could in general use more comments
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Routine for displaying progress within Matlab
    subroutine displayProgress( prog )
        integer,intent(in) :: prog
    
        integer :: mexCallMATLAB, nlhs_cb, nrhs_cb, tmp
        mwPointer plhs_cb(1), prhs_cb(1),mxCreateString
        character*(4) functionName_cb,prog_str   
        character(len=8) :: fmt
        logical :: ex
    
        fmt = '(I4.2)'
        write (prog_str,fmt) prog
    
        !! test if we are inside Matlab or called from a stand-alone. The simple test
        !! is whether io.txt exists in the current path (false for Matlab, true for stand-alone)
        !! nothing bad should happen if in fact we are called from ML but the file somehow
        !! exists - the written output will just not be shown to the user
    
        inquire( file='io.txt', EXIST=ex )
    
        if ( ex .eq. .true. ) then
            write(*,*) prog_str
        else            
            functionName_cb = "disp"
            nlhs_cb = 0
            nrhs_cb = 1
      
            prhs_cb(1) = mxCreateString(prog_str)
      
            tmp = mexCallMATLAB(nlhs_cb, plhs_cb, nrhs_cb, prhs_cb, "disp")
        endif
    
    end subroutine displayProgress
    
   
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Display the iteration within Matlab
    !!
    function displayIteration_Matlab( err, err_max )
        real,intent(in) :: err,err_max
        integer :: displayIteration_Matlab
        integer :: mexCallMATLAB, nlhs_cb, nrhs_cb, tmp
        mwPointer plhs_cb(1), prhs_cb(1),mxCreateString
        character*(4) :: functionName_cb
        character*(100) :: prog_str   
        character(len=8) :: fmt
        logical :: ex

        fmt = '(f15.7)'
        write (prog_str,*) "err: ",err,"Max. err: ",err_max
    
        functionName_cb = "disp"
        nlhs_cb = 0
        nrhs_cb = 1
      
        prhs_cb(1) = mxCreateString(prog_str)
      
        tmp = mexCallMATLAB(nlhs_cb, plhs_cb, nrhs_cb, prhs_cb, "disp")
    
        displayIteration_Matlab = 1
    
    end function displayIteration_Matlab
    
    
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Loads the data from Matlab into a fortran struct
    !!
    subroutine loadMagTile( prhs, cylTile, n_tiles )
        mwPointer, intent(in) :: prhs
        character(len=10),dimension(:),allocatable :: fieldnames    
        integer,intent(in) :: n_tiles    
        type(MagTile),intent(inout),dimension(n_tiles) :: cylTile
    
        mwIndex :: i
        mwSize :: sx
        integer :: nfields
        real*8,dimension(3) :: rectDims
        mwPointer :: r0Ptr,theta0Ptr,z0Ptr,drPtr,dthetaPtr,dzPtr,MPtr,u_eaPtr,u_oa1Ptr,u_oa2Ptr,mur_eaPtr,mur_oaPtr,MremPtr
        mwPointer :: tileTypePtr,offsetPtr,rotAnglesPtr,rectDimsPtr,magnetTypePtr,stateFunctionIndexPtr,includeInIterationPtr
        mwPointer :: mxGetField, mxGetPr,colorPtr,symmetryPtr,symmetryOpsPtr,MrelPtr,vertPtr
    
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
            includeInIterationPtr =  mxGetField(prhs,i,fieldnames(20))
            colorPtr =  mxGetField(prhs,i,fieldnames(21))
            symmetryPtr =  mxGetField(prhs,i,fieldnames(22))
            symmetryOpsPtr =  mxGetField(prhs,i,fieldnames(23))
            MrelPtr =  mxGetField(prhs,i,fieldnames(24))
            VertPtr =  mxGetField(prhs,i,fieldnames(25))
      
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
            call mxCopyPtrToInteger4(mxGetPr(includeInIterationPtr), cylTile(i)%includeInIteration, sx )
            call mxCopyPtrToInteger4(mxGetPr(symmetryPtr), cylTile(i)%exploitSymmetry, sx )
            
            sx = 3
            call mxCopyPtrToReal8(mxGetPr(offsetPtr), cylTile(i)%offset, sx )
            call mxCopyPtrToReal8(mxGetPr(rotAnglesPtr), cylTile(i)%rotAngles, sx )
            call mxCopyPtrToReal8(mxGetPr(rectDimsPtr), rectDims, sx )
            cylTile(i)%a = rectDims(1)
            cylTile(i)%b = rectDims(2)
            cylTile(i)%c = rectDims(3)
            
            call mxCopyPtrToReal8(mxGetPr(colorPtr), cylTile(i)%color, sx )
            
            sx = 1
            call mxCopyPtrToInteger4(mxGetPr(magnetTypePtr), cylTile(i)%magnetType, sx )            
            call mxCopyPtrToInteger4(mxGetPr(stateFunctionIndexPtr), cylTile(i)%stateFunctionIndex, sx )
            
            sx = 3
            call mxCopyPtrToReal8(mxGetPr(symmetryOpsPtr), cylTile(i)%symmetryOps, sx )
            
            sx = 1
            call mxCopyPtrToReal8(mxGetPr(MrelPtr), cylTile(i)%Mrel, sx )
            
            sx = 12
            call mxCopyPtrToReal8(mxGetPr(VertPtr), cylTile(i)%vert, sx )
            
            !cylTile(i)%fieldEvaluation = fieldEvaluationCentre
            if ( cyltile(i)%tiletype == tiletypecylpiece ) then
                call setupEvaluationPoints( cylTile(i) )
            endif
            
        enddo        
        deallocate(fieldnames)
    end subroutine loadMagTile
    
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Returns the data to matlab
    !!
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
    
        !! Load the result back to Matlab      
        ComplexFlag = 0
      
        dims(1) = 1
        dims(2) = n_tiles
        sx = 2
        
        !! Create the return array of structs      
        plhs = mxCreateStructArray( sx, dims, nfields, fieldnames)
      
        !! Pointer for each data entry
        allocate( pvalue(n_tiles,nfields) )
      
        !! Fill in the data for each tile
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
          
          
            classid = mxClassIDFromClassName('int32')
            pvalue(i,20) = mxCreateNumericMatrix(s1,s2,classid,ComplexFlag)
            call mxCopyInteger4ToPtr( cylTile(i)%includeInIteration, mxGetPr( pvalue(i,20) ), sx )
            call mxSetField( plhs, i, fieldnames(20), pvalue(i,20) )
          
            !! load the color assumed to be an RGB triplet, i.e. an array with 3 elements
            s1 = 1
            s2 = 3
            sx = 3
            pvalue(i,21) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
            call mxCopyReal8ToPtr( cylTile(i)%color, mxGetPr( pvalue(i,21) ), sx )
            call mxSetField( plhs, i, fieldnames(21), pvalue(i,21) )
          
            !! load the symmetry flag (if false don't assume symmetry, if true (1) then assume symmetry
            s1 = 1
            s2 = 1
            sx = 1
            classid = mxClassIDFromClassName('int32')
            pvalue(i,22) = mxCreateNumericMatrix(s1,s2,classid,ComplexFlag)
            call mxCopyInteger4ToPtr( cylTile(i)%exploitSymmetry, mxGetPr( pvalue(i,22) ), sx )
            call mxSetField( plhs, i, fieldnames(22), pvalue(i,22) )
          
            s1 = 1
            s2 = 3
            sx = 3
            pvalue(i,23) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
            call mxCopyReal8ToPtr( cylTile(i)%symmetryOps, mxGetPr( pvalue(i,23) ), sx )
            call mxSetField( plhs, i, fieldnames(23), pvalue(i,23) )
          
            s1 = 1
            s2 = 1
            sx = 1
            pvalue(i,24) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
            call mxCopyReal8ToPtr( cylTile(i)%Mrel, mxGetPr( pvalue(i,24) ), sx )
            call mxSetField( plhs, i, fieldnames(24), pvalue(i,24) )
          
            
            s1 = 3
            s2 = 4
            sx = 12
            pvalue(i,25) = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
            call mxCopyReal8ToPtr( cylTile(i)%vert, mxGetPr( pvalue(i,25) ), sx )
            call mxSetField( plhs, i, fieldnames(25), pvalue(i,25) )
            
        enddo
      
        deallocate(pvalue,fieldnames)
    
    end subroutine returnMagTile
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Returns an array with the names of the fields expected in the tile struct
    !!
    subroutine getTileFieldnames( fieldnames, nfields)
        integer,intent(out) :: nfields
        integer,parameter :: nf=25
        character(len=10),dimension(:),intent(out),allocatable :: fieldnames
            
        nfields = nf
        allocate(fieldnames(nfields))
        
        !! Setup the names of the members of the input struct
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
        fieldnames(20) = 'inclIter'
        fieldnames(21) = 'color'
        fieldnames(22) = 'useSymm'
        fieldnames(23) = 'symmOps'
        fieldnames(24) = 'Mrel'
        fieldnames(25) = 'vertices'
        
    end subroutine getTileFieldnames
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Loads the statefunction struct from Matlab to Fortran
    subroutine loadMagStateFunction( prhs, stateFunction, n )
        mwPointer, intent(in) :: prhs
        type(MagStateFunction),dimension(n),intent(inout) :: stateFunction
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
            allocate( stateFunction(i)%y2a(stateFunction(i)%nH) )
        
            sx = stateFunction(i)%nT
            call mxCopyPtrToReal8(mxGetPr(TPtr), stateFunction(i)%T, sx )
            sx = stateFunction(i)%nH
            call mxCopyPtrToReal8(mxGetPr(HPtr), stateFunction(i)%H, sx )
        
            sx = stateFunction(i)%nT*stateFunction(i)%nH
            call mxCopyPtrToReal8(mxGetPr(MPtr), stateFunction(i)%M, sx )
        
            !! make the spline derivatives for later interpolation
            !call splie2( sngl(stateFunction(i)%T), sngl(stateFunction(i)%H), sngl(stateFunction(i)%M), sngl(stateFunction(i)%nT), sngl(stateFunction(i)%nH), sngl(stateFunction(i)%y2a) )
            !!@todo If this code is deprecated, the spline.f90 code can be removed as it is only called here.
            !    call spline( stateFunction(i)%H, stateFunction(i)%M(1,:), 1e30, 1e30, stateFunction(i)%y2a, stateFunction(i)%nH )
        
        enddo
    
        deallocate( fieldnames )
    
    end subroutine loadMagStateFunction
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !!Returns an array with the field names expected in the state function struct
    !!
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
    