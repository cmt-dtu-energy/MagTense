#include "fintrf.h"
    
    module MagTenseMicroMagIO
    use MicroMagParameters
        
    implicit none
    
    contains
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> Loads the data struct problem from Matlab into a Fortran struct
    !> @param[in] prhs pointer to the Matlab data struct
    !> @param[in] problem struct for the internal Fortran represantation of the problem
    !>-----------------------------------------
    subroutine loadMicroMagProblem( prhs, problem )
        mwPointer, intent(in) :: prhs
        type(MicroMagProblem),intent(inout) :: problem
        
        character(len=10),dimension(:),allocatable :: problemFields
        mwIndex :: i
        mwSize :: sx
        integer :: nFieldsProblem, ntot, nt, nt_Hext, useCuda, status, nt_alpha, useCVODE, nt_conv, nnodes, nvalues, nrows, usePrecision, useReturnHall
        mwPointer :: nGridPtr, LGridPtr, dGridPtr, typeGridPtr, ueaProblemPtr, modeProblemPtr, solverProblemPtr
        mwPointer :: A0ProblemPtr, MsProblemPtr, K0ProblemPtr, gammaProblemPtr, alpha0ProblemPtr, MaxT0ProblemPtr
        mwPointer :: ntProblemPtr, m0ProblemPtr, HextProblemPtr, alphaProblemPtr, tProblemPtr, useCudaPtr, useCVODEPtr, nThreadPtr
        mwPointer :: mxGetField, mxGetPr, mxGetM, mxGetN, mxGetNzmax, mxGetIr, mxGetJc
        mwPointer :: ntHextProblemPtr, demThresProblemPtr, demApproxPtr, setTimeDisplayProblemPtr, CVThresProblemPtr
        mwPointer :: NFileReturnPtr, NReturnPtr, NLoadPtr, mxGetString, NFileLoadPtr
        mwPointer :: tolProblemPtr, thres_valueProblemPtr
        mwPointer :: exch_matProblemPtr, irPtr, jcPtr
        mwPointer :: genericProblemPtr
        mwPointer :: ptsGridPtr, nodesGridPtr, elementsGridPtr, nnodesGridPtr
        mwPointer :: valuesPtr, rows_startPtr, rows_endPtr,  colsPtr, nValuesSparsePtr, nRowsSparsePtr
        mwPointer :: usePrecisionPtr, N_aveProblemPtr, useReturnHallProblemPtr
        mwPointer :: demag_ignore_stepsProblemPtr
        integer,dimension(3) :: int_arr
        real(DP),dimension(3) :: real_arr
        real(DP) :: demag_fac, CV
    
        !Get the expected names of the fields
        call getProblemFieldnames( problemFields, nFieldsProblem)
                   
        sx = 3
        i = 1
        nGridPtr = mxGetField( prhs, i, problemFields(1) )
        call mxCopyPtrToInteger4( mxGetPr(nGridPtr), int_arr, sx )
        problem%grid%nx = int_arr(1)
        problem%grid%ny = int_arr(2)
        problem%grid%nz = int_arr(3)
        ntot = product(int_arr)
        
        LGridPtr = mxGetField( prhs, i, problemFields(2) )
        call mxCopyPtrToReal8( mxGetPr(LGridPtr), real_arr, sx )
        problem%grid%Lx = real_arr(1)
        problem%grid%Ly = real_arr(2)
        problem%grid%Lz = real_arr(3)
        
        
        problem%grid%dx = problem%grid%Lx / problem%grid%nx
        problem%grid%dy = problem%grid%Ly / problem%grid%ny
        problem%grid%dz = problem%grid%Lz / problem%grid%nz
        
        
        sx = 1
        typeGridPtr = mxGetField( prhs, i, problemFields(3) )
        call mxCopyPtrToInteger4(mxGetPr(typeGridPtr), problem%grid%gridType, sx )
        
        !Load additional things for a tetrahedron grid
        if ( problem%grid%gridType .eq. gridTypeTetrahedron ) then
            !The center points of all the tetrahedron elements           
            allocate( problem%grid%pts(ntot,3) )
            sx = ntot * 3
            ptsGridPtr = mxGetField( prhs, i, problemFields(35) )
            call mxCopyPtrToReal8(mxGetPr(ptsGridPtr), problem%grid%pts, sx )
            
            !The elements of all the tetrahedron elements
            allocate( problem%grid%elements(4,ntot) )
            sx = ntot * 4
            nodesGridPtr = mxGetField( prhs, i, problemFields(36) )
            call mxCopyPtrToInteger4(mxGetPr(nodesGridPtr), problem%grid%elements, sx )
            
            !The number of nodes in the tetrahedron mesh
            sx = 1
            nnodesGridPtr = mxGetField( prhs, i, problemFields(38) )
            call mxCopyPtrToInteger4(mxGetPr(nnodesGridPtr), problem%grid%nnodes, sx )
            
            !The nodes of all the tetrahedron elements
            nnodes = problem%grid%nnodes
            allocate( problem%grid%nodes(3,nnodes) )
            sx = nnodes * 3
            nodesGridPtr = mxGetField( prhs, i, problemFields(37) )
            call mxCopyPtrToReal8(mxGetPr(nodesGridPtr), problem%grid%nodes, sx )
            
            !the number of nodes in the tetrahedron mesh
            sx = 1
            nnodesGridPtr = mxGetField( prhs, i, problemFields(38) )
            call mxCopyPtrToInteger4(mxGetPr(nnodesGridPtr), problem%grid%nnodes, sx )
        endif
        
        !Load additional things for a grid of unstructured prisms
        if ( problem%grid%gridType .eq. gridTypeUnstructuredPrisms ) then
            !The center points of all the prisms elements           
            allocate( problem%grid%pts(ntot,3) )
            sx = ntot * 3
            ptsGridPtr = mxGetField( prhs, i, problemFields(35) )
            call mxCopyPtrToReal8(mxGetPr(ptsGridPtr), problem%grid%pts, sx )
            
            
            !The side lengths of all the prisms
            allocate( problem%grid%abc(ntot,3) )
            sx = ntot * 3
            nodesGridPtr = mxGetField( prhs, i, problemFields(45) )
            call mxCopyPtrToReal8(mxGetPr(nodesGridPtr), problem%grid%abc, sx )
            
        endif
                
        !Finished loading the grid------------------------------------------
        
        !Start loading the problem
        !Allocate memory for the easy axis vectors
        allocate( problem%u_ea(ntot,3) )
        ueaProblemPtr = mxGetField(prhs,i,problemFields(4))
        sx = ntot * 3
        call mxCopyPtrToReal8(mxGetPr(ueaProblemPtr), problem%u_ea, sx )
        
        sx = 1
        modeProblemPtr = mxGetField( prhs, i, problemFields(5) )
        call mxCopyPtrToInteger4(mxGetPr(modeProblemPtr), problem%ProblemMode, sx )
        
        sx = 1
        solverProblemPtr = mxGetField( prhs, i, problemFields(6) )
        call mxCopyPtrToInteger4(mxGetPr(solverProblemPtr), problem%solver, sx )
        
        sx = 1
        A0ProblemPtr = mxGetField( prhs, i, problemFields(7) )
        call mxCopyPtrToReal8(mxGetPr(A0ProblemPtr), problem%A0, sx )
        
        allocate( problem%Ms(ntot) )
        sx = ntot
        MsProblemPtr = mxGetField( prhs, i, problemFields(8) )
        call mxCopyPtrToReal8(mxGetPr(MsProblemPtr), problem%Ms, sx )
        
        allocate( problem%K0(ntot) )
        sx = ntot
        K0ProblemPtr = mxGetField( prhs, i, problemFields(9) )
        call mxCopyPtrToReal8(mxGetPr(K0ProblemPtr), problem%K0, sx )
        
        sx = 1
        gammaProblemPtr = mxGetField( prhs, i, problemFields(10) )
        call mxCopyPtrToReal8(mxGetPr(gammaProblemPtr), problem%gamma, sx )
        
        sx = 1
        alpha0ProblemPtr = mxGetField( prhs, i, problemFields(11) )
        call mxCopyPtrToReal8(mxGetPr(alpha0ProblemPtr), problem%alpha0, sx )
        
        sx = 1
        MaxT0ProblemPtr = mxGetField( prhs, i, problemFields(12) )
        call mxCopyPtrToReal8(mxGetPr(MaxT0ProblemPtr), problem%MaxT0, sx )                
        
        !load the no. of time steps in the applied field
        sx = 1
        ntHextProblemPtr = mxGetField( prhs, i, problemFields(13) )
        call mxCopyPtrToInteger4(mxGetPr(ntHextProblemPtr), nt_Hext, sx )
        
        !Applied field as a function of time evaluated at the timesteps specified in nt_Hext
        !problem%Hext(:,1) is the time grid while problem%Hext(:,2:4) are the x-,y- and z-components of the applied field
        sx = nt_Hext * 4
        allocate( problem%Hext(nt_Hext,4) )
        HextProblemPtr = mxGetField( prhs, i, problemFields(14) )
        call mxCopyPtrToReal8(mxGetPr(HextProblemPtr), problem%Hext, sx )
                
        !Load the no. of time steps required
        sx = 1
        ntProblemPtr = mxGetField( prhs, i, problemFields(15) )
        call mxCopyPtrToInteger4(mxGetPr(ntProblemPtr), nt, sx )
        
        allocate( problem%t(nt) )
        tProblemPtr = mxGetField(prhs,i,problemFields(16) )
        sx = nt
        call mxCopyPtrToReal8(mxGetPr(tProblemPtr), problem%t, sx )
        
        !Initial magnetization
        allocate( problem%m0(3*ntot) )
        m0ProblemPtr = mxGetField(prhs,i,problemFields(17))
        sx = ntot * 3
        call mxCopyPtrToReal8(mxGetPr(m0ProblemPtr), problem%m0, sx )
        
        !Demagnetization threshold value        
        demThresProblemPtr = mxGetField(prhs,i,problemFields(18))
        sx = 1
        call mxCopyPtrToReal8(mxGetPr(demThresProblemPtr), demag_fac, sx )
            
        problem%demag_threshold = sngl(demag_fac)
        
        sx = 1
        useCudaPtr = mxGetField( prhs, i, problemFields(19) )
        call mxCopyPtrToInteger4(mxGetPr(useCudaPtr), useCuda, sx )
        if ( useCuda .eq. 1 ) then
            problem%useCuda = useCudaTrue
        else
            problem%useCuda = useCudaFalse
        endif
               
        sx = 1
        demApproxPtr = mxGetField( prhs, i, problemFields(20) )
        call mxCopyPtrToInteger4(mxGetPr(demApproxPtr), problem%demag_approximation, sx )
        
        !flag whether the demag tensor should be returned and if so how
        sx = 1
        NReturnPtr = mxGetField( prhs, i, problemFields(21) )
        call mxCopyPtrToInteger4(mxGetPr(NReturnPtr), problem%demagTensorReturnState, sx )
        
        !File for returning the demag tensor to a file on disk (has to have length>2)
        if ( problem%demagTensorReturnState .gt. 2 ) then
            !Length of the file name
            sx = problem%demagTensorReturnState
            NFileReturnPtr = mxGetField( prhs, i, problemFields(22) )            
            status = mxGetString( NFileReturnPtr, problem%demagTensorFileOut, sx )
        endif
        
        !flag whether the demag tensor should be loaded
        sx = 1
        NLoadPtr = mxGetField( prhs, i, problemFields(23) )
        call mxCopyPtrToInteger4(mxGetPr(NLoadPtr), problem%demagTensorLoadState, sx )
        
        !File for loading the demag tensor to a file on disk (has to have length>2)
        if ( problem%demagTensorLoadState .gt. 2 ) then
            !Length of the file name
            sx = problem%demagTensorLoadState
            NFileLoadPtr = mxGetField( prhs, i, problemFields(24) )            
            status = mxGetString( NFileLoadPtr, problem%demagTensorFileIn, sx )
        endif
        
        
        problem%setTimeDisplay = 100
        
        !Set how often to display the timestep in Matlab
        sx = 1
        setTimeDisplayProblemPtr = mxGetField( prhs, i, problemFields(25) )
        call mxCopyPtrToInteger4(mxGetPr(setTimeDisplayProblemPtr), problem%setTimeDisplay, sx )
        
        !Load the no. of times in the alpha function
        sx = 1
        ntProblemPtr = mxGetField( prhs, i, problemFields(26) )
        call mxCopyPtrToInteger4(mxGetPr(ntProblemPtr), nt_alpha, sx )
        
        !alpha as a function of time evaluated at the timesteps
        !problem%alpha(:,1) is the time grid while problem%alpha(:,2) are the alpha values
        sx = nt_alpha * 2
        allocate( problem%alpha(nt_alpha,2) )
        alphaProblemPtr = mxGetField( prhs, i, problemFields(27) )
        call mxCopyPtrToReal8(mxGetPr(alphaProblemPtr), problem%alpha, sx )
        
        sx = 1
        tolProblemPtr = mxGetField( prhs, i, problemFields(28) )
        call mxCopyPtrToReal8(mxGetPr(tolProblemPtr), problem%tol, sx )
        
        sx = 1
        thres_valueProblemPtr = mxGetField( prhs, i, problemFields(29) )
        call mxCopyPtrToReal8(mxGetPr(thres_valueProblemPtr), problem%thres_value, sx )

        sx = 1
        useCVODEPtr = mxGetField( prhs, i, problemFields(30) )
        call mxCopyPtrToInteger4(mxGetPr(useCVODEPtr), useCVODE, sx )
        if ( useCVODE .eq. 1 ) then
            problem%useCVODE = useCVODETrue
        else
            problem%useCVODE = useCVODEFalse
        endif
        
        !File for loading the sparse exchange tensor from Matlab (for non-uniform grids)
        if (( problem%grid%gridType .eq. gridTypeTetrahedron ) .or. (problem%grid%gridType .eq. gridTypeUnstructuredPrisms)) then
            ! Load the CSR sparse information from Matlab
            sx = 1
            nValuesSparsePtr = mxGetField( prhs, i, problemFields(39) )
            call mxCopyPtrToInteger4(mxGetPr(nValuesSparsePtr), problem%grid%A_exch_load%nvalues, sx )
       
            sx = 1
            nRowsSparsePtr = mxGetField( prhs, i, problemFields(40) )
            call mxCopyPtrToInteger4(mxGetPr(nRowsSparsePtr), problem%grid%A_exch_load%nrows, sx )
            
            nvalues = problem%grid%A_exch_load%nvalues
            nrows = problem%grid%A_exch_load%nrows
            allocate( problem%grid%A_exch_load%values(nvalues), problem%grid%A_exch_load%rows_start(nrows) , problem%grid%A_exch_load%rows_end(nrows) , problem%grid%A_exch_load%cols(nvalues) )
             
            sx = nvalues
            valuesPtr = mxGetField( prhs, i, problemFields(41) )
            call mxCopyPtrToReal8(mxGetPr(valuesPtr), problem%grid%A_exch_load%values, sx )
            
            sx = nrows
            rows_startPtr = mxGetField( prhs, i, problemFields(42) )
            call mxCopyPtrToInteger4(mxGetPr(rows_startPtr), problem%grid%A_exch_load%rows_start, sx )
        
            sx = nrows
            rows_endPtr = mxGetField( prhs, i, problemFields(43) )
            call mxCopyPtrToInteger4(mxGetPr(rows_endPtr), problem%grid%A_exch_load%rows_end, sx )
        
            sx = nvalues
            colsPtr = mxGetField( prhs, i, problemFields(44) )
            call mxCopyPtrToInteger4(mxGetPr(colsPtr), problem%grid%A_exch_load%cols, sx )
        endif
          
        !Load the no. of time steps in the time convergence array
        sx = 1
        ntProblemPtr = mxGetField( prhs, i, problemFields(32) )
        call mxCopyPtrToInteger4(mxGetPr(ntProblemPtr), nt_conv, sx )
        
        allocate( problem%t_conv(nt_conv) )
        genericProblemPtr = mxGetField(prhs,i,problemFields(33) )
        sx = nt_conv
        call mxCopyPtrToReal8(mxGetPr(genericProblemPtr), problem%t_conv, sx )
        
        sx = 1
        genericProblemPtr = mxGetField( prhs, i, problemFields(34) )
        call mxCopyPtrToReal8(mxGetPr(genericProblemPtr), problem%conv_tol, sx )

        sx = 1
        usePrecisionPtr = mxGetField( prhs, i, problemFields(46) )
        call mxCopyPtrToInteger4(mxGetPr(usePrecisionPtr), usePrecision, sx )
        if ( usePrecision .eq. 1 ) then
            problem%usePrecision = usePrecisionTrue
        else
            problem%usePrecision = usePrecisionFalse
        endif
        
        sx = 1
        nThreadPtr = mxGetField( prhs, i, problemFields(47) )
        call mxCopyPtrToInteger4(mxGetPr(nThreadPtr), problem%nThreadsMatlab, sx )
        
        sx = 3
        N_aveProblemPtr = mxGetField( prhs, i, problemFields(48) )
        call mxCopyPtrToInteger4(mxGetPr(N_aveProblemPtr), problem%N_ave, sx )
        
        !Coefficient of variation value
        problem%CV = 0
        
        sx = 1
        CVThresProblemPtr = mxGetField(prhs,i,problemFields(49))
        call mxCopyPtrToReal8(mxGetPr(CVThresProblemPtr), CV, sx )
            
        problem%CV = sngl(CV)
        
        !Parameter to determine if the specific H_fields are returned (exchange, demag, etc.)
        sx = 1
        useReturnHallProblemPtr = mxGetField(prhs,i,problemFields(50))
        
        call mxCopyPtrToInteger4(mxGetPr(useReturnHallProblemPtr), useReturnHall, sx )
        if ( useReturnHall .eq. 1 ) then
            problem%useReturnHall = useReturnHallTrue
        else
            problem%useReturnHall = useReturnHallFalse
        endif
        
        sx = 1
        demag_ignore_stepsProblemPtr = mxGetField( prhs, i, problemFields(51) )
        call mxCopyPtrToInteger4(mxGetPr(demag_ignore_stepsProblemPtr), problem%demag_ignore_steps, sx )
        
        !Clean-up 
        deallocate(problemFields)
    end subroutine loadMicroMagProblem
    
    
   !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> Returns the solution data struct from Fortran to Matlab
    !> @param[in] solution struct for the internal Fortran represantation of the solution
    !> @param[in] plhs pointer to the Matlab data struct    
    !>-----------------------------------------
    subroutine returnMicroMagSolution( solution, plhs )
        type(MicroMagSolution),intent(in) :: solution           !> Solution to be copied to Matlab        
        mwPointer,intent(inout) :: plhs
    
        integer :: ComplexFlag,classid,mxClassIDFromClassName
        mwSize,dimension(1) :: dims
        mwSize :: s1,s2,sx,ndim
        mwSize,dimension(4) :: dims_4
        mwPointer :: pt,pm,pp,pdem,pext,pexc,pani
        mwPointer :: mxCreateStructArray, mxCreateDoubleMatrix,mxGetPr,mxCreateNumericMatrix,mxCreateNumericArray
        mwIndex :: ind
        character(len=10),dimension(:),allocatable :: fieldnames    
        integer :: nfields,ntot ,nt
    
        call getSolutionFieldnames( fieldnames, nfields)
    
        nt = size(solution%t_out)
        ntot = size(solution%M_out(1,:,1,1))
        
        ! Load the result back to Matlab      
        ComplexFlag = 0
      
        dims(1) = 1        
        sx = 1        
        ! Create the return array of structs      
        plhs = mxCreateStructArray( sx, dims, nfields, fieldnames)
      
        ind = 1
        
        s1 = nt
        s2 = 1
        pt = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
        sx = s1 * s2
        call mxCopyReal8ToPtr( solution%t_out, mxGetPr( pt ), sx )
        call mxSetField( plhs, ind, fieldnames(1), pt )
          
        
        ndim = 4
        dims_4(1) = nt
        dims_4(2) = ntot
        dims_4(3) = size( solution%M_out(1,1,:,1) )
        dims_4(4) = 3
        classid = mxClassIDFromClassName( 'double' )
        !pm = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
        pm = mxCreateNumericArray( ndim, dims_4, classid, ComplexFlag)
        sx = dims_4(1) * dims_4(2) * dims_4(3) * dims_4(4)
        call mxCopyReal8ToPtr( solution%M_out, mxGetPr( pm ), sx )
        call mxSetField( plhs, ind, fieldnames(2), pm )
      
        
        s1 = ntot
        s2 = 3
        pp = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
        sx = s1 * s2
        call mxCopyReal8ToPtr( solution%pts, mxGetPr( pp ), sx )
        call mxSetField( plhs, ind, fieldnames(3), pp )
        
        
        
        ndim = 4
        dims_4(1) = size( solution%H_exc(:,1,1,1) )
        dims_4(2) = size( solution%H_exc(1,:,1,1) )
        dims_4(3) = size( solution%H_exc(1,1,:,1) )
        dims_4(4) = 3
        classid = mxClassIDFromClassName( 'double' )
        pexc = mxCreateNumericArray( ndim, dims_4, classid, ComplexFlag)
        sx = dims_4(1) * dims_4(2) * dims_4(3) * dims_4(4)
        call mxCopyReal8ToPtr( solution%H_exc, mxGetPr( pexc ), sx )
        call mxSetField( plhs, ind, fieldnames(4), pexc )
        
        
        
        ndim = 4
        dims_4(1) = size( solution%H_ext(:,1,1,1) )
        dims_4(2) = size( solution%H_ext(1,:,1,1) )
        dims_4(3) = size( solution%H_ext(1,1,:,1) )
        dims_4(4) = 3
        classid = mxClassIDFromClassName( 'double' )
        pext = mxCreateNumericArray( ndim, dims_4, classid, ComplexFlag)
        sx = dims_4(1) * dims_4(2) * dims_4(3) * dims_4(4)
        call mxCopyReal8ToPtr( solution%H_ext, mxGetPr( pext ), sx )
        call mxSetField( plhs, ind, fieldnames(5), pext )
        
        
        
        ndim = 4
        dims_4(1) = size( solution%H_dem(:,1,1,1) )
        dims_4(2) = size( solution%H_dem(1,:,1,1) )
        dims_4(3) = size( solution%H_dem(1,1,:,1) )
        dims_4(4) = 3
        classid = mxClassIDFromClassName( 'double' )
        pdem = mxCreateNumericArray( ndim, dims_4, classid, ComplexFlag)
        sx = dims_4(1) * dims_4(2) * dims_4(3) * dims_4(4)
        call mxCopyReal8ToPtr( solution%H_dem, mxGetPr( pdem ), sx )
        call mxSetField( plhs, ind, fieldnames(6), pdem )
        
        
        
        ndim = 4
        dims_4(1) = size( solution%H_ani(:,1,1,1) )
        dims_4(2) = size( solution%H_ani(1,:,1,1) )
        dims_4(3) = size( solution%H_ani(1,1,:,1) )
        dims_4(4) = 3
        classid = mxClassIDFromClassName( 'double' )
        pani = mxCreateNumericArray( ndim, dims_4, classid, ComplexFlag)
        sx = dims_4(1) * dims_4(2) * dims_4(3) * dims_4(4)
        call mxCopyReal8ToPtr( solution%H_ani, mxGetPr( pani ), sx )
        call mxSetField( plhs, ind, fieldnames(7), pani )
        
       
        !Clean up
        deallocate(fieldnames)
    
    end subroutine returnMicroMagSolution
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> Returns an array with the names of the fields expected in the MicroMagProblem struct
    !> @param[inout] fieldnames, array of the names of the fields
    !> @param[inout] nfields the no. of elements in fieldnames
    !>-----------------------------------------
    subroutine getProblemFieldnames( fieldnames, nfields)
        integer,intent(out) :: nfields        
        integer,parameter :: nf=51
        character(len=10),dimension(:),intent(out),allocatable :: fieldnames
            
        nfields = nf
        allocate(fieldnames(nfields))
        
        !! Setup the names of the members of the input struct
        fieldnames(1) = 'grid_n'
        fieldnames(2) = 'grid_L'
        fieldnames(3) = 'grid_type'
        fieldnames(4) = 'u_ea'
        fieldnames(5) = 'ProblemMode'
        fieldnames(6) = 'solver'
        fieldnames(7) = 'A0'
        fieldnames(8) = 'Ms'
        fieldnames(9) = 'K0'
        fieldnames(10) = 'gamma'
        fieldnames(11) = 'alpha'
        fieldnames(12) = 'MaxT0'
        fieldnames(13) = 'nt_Hext'
        fieldnames(14) = 'Hext'    
        fieldnames(15) = 'nt'
        fieldnames(16) = 't'
        fieldnames(17) = 'm0'
        fieldnames(18) = 'dem_thres'
        fieldnames(19) = 'useCuda'
        fieldnames(20) = 'dem_appr'        
        fieldnames(21) = 'N_ret'
        fieldnames(22) = 'N_file_out'
        fieldnames(23) = 'N_load'
        fieldnames(24) = 'N_file_in'
        fieldnames(25) = 'setTimeDis'
        fieldnames(26) = 'nt_alpha'
        fieldnames(27) = 'alphat'
        fieldnames(28) = 'tol'
        fieldnames(29) = 'thres'
        fieldnames(30) = 'useCVODE'
        fieldnames(31) = 'exch_mat'
        fieldnames(32) = 'nt_conv'
        fieldnames(33) = 't_conv'
        fieldnames(34) = 'conv_tol'
        fieldnames(35) = 'grid_pts'
        fieldnames(36) = 'grid_ele'
        fieldnames(37) = 'grid_nod'
        fieldnames(38) = 'grid_nnod'
        fieldnames(39) = 'exch_nval'
        fieldnames(40) = 'exch_nrow'
        fieldnames(41) = 'exch_val'
        fieldnames(42) = 'exch_rows'
        fieldnames(43) = 'exch_rowe'
        fieldnames(44) = 'exch_col'
        fieldnames(45) = 'grid_abc'
        fieldnames(46) = 'usePres'
        fieldnames(47) = 'nThreads'
        fieldnames(48) = 'N_ave'
        fieldnames(49) = 'CV'
        fieldnames(50) = 'ReturnHall'
        fieldnames(51) = 'demigstp'
        
    end subroutine getProblemFieldnames
    
   
    
     !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> Returns an array with the names of the fields expected in the MicroMagSolution struct
    !> @param[inout] fieldnames, array of the names of the fields
    !> @param[inout] nfields the no. of elements in fieldnames
    !>-----------------------------------------
    subroutine getSolutionFieldnames( fieldnames, nfields)
        integer,intent(out) :: nfields
        integer,parameter :: nf=7
        character(len=10),dimension(:),intent(out),allocatable :: fieldnames
            
        nfields = nf
        allocate(fieldnames(nfields))
        
        !! Setup the names of the members of the input struct
        fieldnames(1) = 't'
        fieldnames(2) = 'M'
        fieldnames(3) = 'pts'
        fieldnames(4) = 'H_exc'
        fieldnames(5) = 'H_ext'
        fieldnames(6) = 'H_dem'
        fieldnames(7) = 'H_ani'
        
        
        
    end subroutine getSolutionFieldnames
    
    !>----------------------------------------
    !> Kaspar K. Nielsen, kasparkn@gmail.com, January 2020
    !> Writes the demag tensors to disk given a filename in problem
    !> @params[in] problem the struct containing the entire problem
    !>----------------------------------------    
    subroutine writeDemagTensorToDisk( problem )
    type(MicroMagProblem), intent(in) :: problem
    
    integer :: n            !> No. of elements in the grid
    
    
    n = problem%grid%nx * problem%grid%ny * problem%grid%nz
        
        open (11, file=problem%demagTensorFileOut,	&
                status='unknown', form='unformatted',	&
                access='direct', recl=1*n*n)

        write(11,rec=1) problem%Kxx
        write(11,rec=2) problem%Kxy
        write(11,rec=3) problem%Kxz
        write(11,rec=4) problem%Kyy
        write(11,rec=5) problem%Kyz
        write(11,rec=6) problem%Kzz

        close(11)
        

    end subroutine writeDemagTensorToDisk


    !>----------------------------------------
    !> Kaspar K. Nielsen, kasparkn@gmail.com, January 2020
    !> Loads the demag tensors from disk given a file in problem
    !> @params[inout] problem the struct containing the entire problem
    !>----------------------------------------
    subroutine loadDemagTensorFromDisk( problem )
    type( MicroMagProblem ), intent(inout) :: problem
    integer :: n

            n = problem%grid%nx * problem%grid%ny * problem%grid%nz
            
            
                open (11, file=problem%demagTensorFileIn,	&
                        status='unknown', form='unformatted',	&
                        access='direct', recl=1*n*n)

            read(11,rec=1) problem%Kxx
            read(11,rec=2) problem%Kxy
            read(11,rec=3) problem%Kxz
            read(11,rec=4) problem%Kyy
            read(11,rec=5) problem%Kyz
            read(11,rec=6) problem%Kzz

            close(11)



    end subroutine loadDemagTensorFromDisk

end module MagTenseMicroMagIO
    