#include "fintrf.h"
    
    module MagTenseMicroMagIO
    use MicroMagParameters
    use IO_GENERAL
    
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
        
        character(len=12),dimension(:),allocatable :: problemFields
        mwIndex :: i
        mwSize :: sx
        integer :: nFieldsProblem, ntot, nt, nt_Hext, useCuda, status, nt_alpha, useCVODE, nt_conv, nnodes, nvalues, nrows, useReturnHall, passExch, UseFMM, useDemag, useAvgN
        mwPointer :: nGridPtr, LGridPtr, dGridPtr, typeGridPtr, ueaProblemPtr, modeProblemPtr, solverProblemPtr
        mwPointer :: exch_weightProblemPtr, exch_methodProblemPtr, exch_interpnProblemPtr
        mwPointer :: A0ProblemPtr, MsProblemPtr, K0ProblemPtr, K1ProblemPtr, K2ProblemPtr, gammaProblemPtr, alpha0ProblemPtr
        mwPointer :: nPhaseProblemPtr, phaseIdProblemPtr, AIntProblemPtr
        real(DP),dimension(:),allocatable :: A_int_flat
        real(DP),dimension(:),allocatable :: phase_id_real
        integer :: ip, jp
        mwPointer :: ntProblemPtr, m0ProblemPtr, HextProblemPtr, alphaProblemPtr, tProblemPtr, useCudaPtr, useCVODEPtr, usePassExchPtr, useAvgNProblemPtr
        mwPointer :: mxGetField, mxGetPr, mxGetM, mxGetN, mxGetNzmax, mxGetIr, mxGetJc
        mwPointer :: ntHextProblemPtr, demThresProblemPtr, demApproxPtr, setTimeDisplayProblemPtr, CVThresProblemPtr
        mwPointer :: NFileReturnPtr, NReturnPtr, NLoadPtr, mxGetString, NFileLoadPtr
        mwPointer :: temperatureProblemPtr, n_macroVecProblemPtr, shiftVecProblemPtr, macroShapeProblemPtr, sampleShapeProblemPtr, exchPBCProblemPtr, dummy_runProblemPtr, fmm_ntermsProblemPtr
        mwPointer :: tolProblemPtr, thres_valueProblemPtr
        mwPointer :: exch_matProblemPtr, irPtr, jcPtr
        mwPointer :: genericProblemPtr
        mwPointer :: ptsGridPtr, nodesGridPtr, elementsGridPtr, nnodesGridPtr
        mwPointer :: valuesPtr, rowsPtr, colsPtr, nValuesSparsePtr, nRowsSparsePtr, nColsSparsePtr
        mwPointer :: N_aveProblemPtr, useReturnHallProblemPtr
        mwPointer :: CrystalAxisProblemPtr, K0_arrProblemPtr
        mwPointer :: fmm_cellsProblemPtr,fmm_epsProblemPtr,ifunifProblemPtr,nlminProblemPtr,nlmaxProblemPtr,use_fmmlProblemPtr,fmm_shortProblemPtr,fmm_min_nProblemPtr
        mwPointer :: useDemagPtr
        mwPointer :: window_enaProblemPtr, window_intProblemPtr, trace_enaProblemPtr, flush_eachProblemPtr, trace_verbProblemPtr
        mwPointer :: N_log_dirPtr, log_dirPtr, N_timer_logPtr, timer_logPtr, N_trace_logPtr, trace_logPtr
        integer :: N_timer_log, N_trace_log, N_log_dir
        mwPointer :: adaptiveHextPtr, maxHextStepsPtr, H_startPtr, H_endPtr
        mwPointer :: dH_initialPtr, dH_minPtr, dH_maxPtr, dH_growPtr, dH_shrinkPtr
        mwPointer :: dM_minPtr, dM_targetPtr, dM_rejectPtr, switch_refine_dHPtr, use_switch_refinePtr
        mwPointer :: rng_seedPtr
        mwPointer :: minimizerPtr
        integer :: use_switch_refine
        integer,dimension(3) :: int_arr
        real(DP),dimension(3) :: real_arr
        real(DP) :: demag_fac, CV, pi, mu0
        character*(40) :: prog_str
            
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
            ptsGridPtr = mxGetField( prhs, i, problemFields(34) )
            call mxCopyPtrToReal8(mxGetPr(ptsGridPtr), problem%grid%pts, sx )
            
            !The elements of all the tetrahedron elements
            allocate( problem%grid%elements(4,ntot) )
            sx = ntot * 4
            nodesGridPtr = mxGetField( prhs, i, problemFields(35) )
            call mxCopyPtrToInteger4(mxGetPr(nodesGridPtr), problem%grid%elements, sx )
            
            !The number of nodes in the tetrahedron mesh
            sx = 1
            nnodesGridPtr = mxGetField( prhs, i, problemFields(37) )
            call mxCopyPtrToInteger4(mxGetPr(nnodesGridPtr), problem%grid%nnodes, sx )
            
            !The nodes of all the tetrahedron elements
            nnodes = problem%grid%nnodes
            allocate( problem%grid%nodes(3,nnodes) )
            sx = nnodes * 3
            nodesGridPtr = mxGetField( prhs, i, problemFields(36) )
            call mxCopyPtrToReal8(mxGetPr(nodesGridPtr), problem%grid%nodes, sx )
            
            !the number of nodes in the tetrahedron mesh
            sx = 1
            nnodesGridPtr = mxGetField( prhs, i, problemFields(37) )
            call mxCopyPtrToInteger4(mxGetPr(nnodesGridPtr), problem%grid%nnodes, sx )
        endif
        
        !Load additional things for a grid of unstructured prisms
        if ( problem%grid%gridType .eq. gridTypeUnstructuredPrisms ) then
            !The center points of all the prisms elements           
            allocate( problem%grid%pts(ntot,3) )
            sx = ntot * 3
            ptsGridPtr = mxGetField( prhs, i, problemFields(34) )
            call mxCopyPtrToReal8(mxGetPr(ptsGridPtr), problem%grid%pts, sx )
            
            
            !The side lengths of all the prisms
            allocate( problem%grid%abc(ntot,3) )
            sx = ntot * 3
            nodesGridPtr = mxGetField( prhs, i, problemFields(43) )
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
        
        allocate( problem%A0(ntot) )
        sx = ntot
        A0ProblemPtr = mxGetField( prhs, i, problemFields(7) )
        call mxCopyPtrToReal8(mxGetPr(A0ProblemPtr), problem%A0, sx )
        
        !----------------- Interface exchange between two materials ------------------------
        !Every field above is required, but these three are not: a problem struct saved
        !before the feature existed, or built by hand, will not have them and mxGetField
        !returns a null pointer. Reading through that would be a hard crash, so the absence
        !is treated as 'one material', which is exactly the previous behaviour.
        problem%n_phase = 1
        nPhaseProblemPtr = mxGetField( prhs, i, problemFields(100) )
        phaseIdProblemPtr = mxGetField( prhs, i, problemFields(101) )
        AIntProblemPtr = mxGetField( prhs, i, problemFields(102) )
        if ( nPhaseProblemPtr .ne. 0 .and. phaseIdProblemPtr .ne. 0 .and. AIntProblemPtr .ne. 0 ) then
            sx = 1
            call mxCopyPtrToInteger4(mxGetPr(nPhaseProblemPtr), problem%n_phase, sx )

            if ( problem%n_phase .gt. 1 ) then
                !MATLAB hands over doubles, so the phase indices arrive as reals and are
                !rounded rather than truncated - 2.9999999 has to become 3, not 2.
                allocate( phase_id_real(ntot) )
                sx = ntot
                call mxCopyPtrToReal8(mxGetPr(phaseIdProblemPtr), phase_id_real, sx )
                allocate( problem%phase_id(ntot) )
                problem%phase_id = nint(phase_id_real)
                deallocate( phase_id_real )

                if ( minval(problem%phase_id) .lt. 1 .or. maxval(problem%phase_id) .gt. problem%n_phase ) then
                    call displayGUIMessage( 'MagTense: phase_id must be between 1 and n_phase' )
                    error stop 'loadMicroMagProblem: phase_id out of range'
                endif

                !A MATLAB matrix arrives in column-major order, which is also how the
                !n_phase x n_phase table is stored here, so the flat copy reshapes directly.
                allocate( A_int_flat(problem%n_phase*problem%n_phase) )
                sx = problem%n_phase*problem%n_phase
                call mxCopyPtrToReal8(mxGetPr(AIntProblemPtr), A_int_flat, sx )
                allocate( problem%A_int(problem%n_phase,problem%n_phase) )
                do jp = 1, problem%n_phase
                    do ip = 1, problem%n_phase
                        problem%A_int(ip,jp) = A_int_flat(ip + (jp-1)*problem%n_phase)
                    end do
                end do
                deallocate( A_int_flat )

                !An asymmetric table would make the exchange across a face depend on which of
                !the two cells is asked, which is not a physical operator.
                if ( maxval(abs(problem%A_int - transpose(problem%A_int))) .gt. 0.0_DP ) then
                    call displayGUIMessage( 'MagTense: the interface exchange table must be symmetric' )
                    error stop 'loadMicroMagProblem: A_int is not symmetric'
                endif
            endif
        endif
        !-----------------------------------------------------------------------------------

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
        
        !load the no. of time steps in the applied field
        sx = 1
        ntHextProblemPtr = mxGetField( prhs, i, problemFields(12) )
        call mxCopyPtrToInteger4(mxGetPr(ntHextProblemPtr), nt_Hext, sx )
        
        !Applied field as a function of time evaluated at the timesteps specified in nt_Hext
        !problem%Hext(:,1) is the time grid while problem%Hext(:,2:4) are the x-,y- and z-components of the applied field
        sx = nt_Hext * 4
        allocate( problem%Hext(nt_Hext,4) )
        HextProblemPtr = mxGetField( prhs, i, problemFields(13) )
        call mxCopyPtrToReal8(mxGetPr(HextProblemPtr), problem%Hext, sx )
                
        !Load the no. of time steps required
        sx = 1
        ntProblemPtr = mxGetField( prhs, i, problemFields(14) )
        call mxCopyPtrToInteger4(mxGetPr(ntProblemPtr), nt, sx )
        
        allocate( problem%t(nt) )
        tProblemPtr = mxGetField(prhs,i,problemFields(15) )
        sx = nt
        call mxCopyPtrToReal8(mxGetPr(tProblemPtr), problem%t, sx )
        
        !Initial magnetization
        allocate( problem%m0(3*ntot) )
        m0ProblemPtr = mxGetField(prhs,i,problemFields(16))
        sx = ntot * 3
        call mxCopyPtrToReal8(mxGetPr(m0ProblemPtr), problem%m0, sx )
        
        !Demagnetization threshold value        
        demThresProblemPtr = mxGetField(prhs,i,problemFields(17))
        sx = 1
        call mxCopyPtrToReal8(mxGetPr(demThresProblemPtr), demag_fac, sx )
            
        problem%demag_threshold = sngl(demag_fac)
        
        sx = 1
        useCudaPtr = mxGetField( prhs, i, problemFields(18) )
        call mxCopyPtrToInteger4(mxGetPr(useCudaPtr), useCuda, sx )
        if ( useCuda .eq. 1 ) then
            problem%useCuda = useCudaTrue
        else
            problem%useCuda = useCudaFalse
        endif
               
        sx = 1
        demApproxPtr = mxGetField( prhs, i, problemFields(19) )
        call mxCopyPtrToInteger4(mxGetPr(demApproxPtr), problem%demag_approximation, sx )
        
        !flag whether the demag tensor should be returned and if so how
        sx = 1
        NReturnPtr = mxGetField( prhs, i, problemFields(20) )
        call mxCopyPtrToInteger4(mxGetPr(NReturnPtr), problem%demagTensorReturnState, sx )
        
        !File for returning the demag tensor to a file on disk (has to have length>2)
        if ( problem%demagTensorReturnState .gt. 2 ) then
            !Length of the file name
            sx = problem%demagTensorReturnState
            NFileReturnPtr = mxGetField( prhs, i, problemFields(21) )            
            status = mxGetString( NFileReturnPtr, problem%demagTensorFileOut, sx )
        endif
        
        !flag whether the demag tensor should be loaded
        sx = 1
        NLoadPtr = mxGetField( prhs, i, problemFields(22) )
        call mxCopyPtrToInteger4(mxGetPr(NLoadPtr), problem%demagTensorLoadState, sx )
        
        !File for loading the demag tensor to a file on disk (has to have length>2)
        if ( problem%demagTensorLoadState .gt. 2 ) then
            !Length of the file name
            sx = problem%demagTensorLoadState
            NFileLoadPtr = mxGetField( prhs, i, problemFields(23) )            
            status = mxGetString( NFileLoadPtr, problem%demagTensorFileIn, sx )
        endif
        
        
        problem%setTimeDisplay = 100
        
        !Set how often to display the timestep in Matlab
        sx = 1
        setTimeDisplayProblemPtr = mxGetField( prhs, i, problemFields(24) )
        call mxCopyPtrToInteger4(mxGetPr(setTimeDisplayProblemPtr), problem%setTimeDisplay, sx )
        
        !Load the no. of times in the alpha function
        sx = 1
        ntProblemPtr = mxGetField( prhs, i, problemFields(25) )
        call mxCopyPtrToInteger4(mxGetPr(ntProblemPtr), nt_alpha, sx )
        
        !alpha as a function of time evaluated at the timesteps
        !problem%alpha(:,1) is the time grid while problem%alpha(:,2) are the alpha values
        sx = nt_alpha * 2
        allocate( problem%alpha(nt_alpha,2) )
        alphaProblemPtr = mxGetField( prhs, i, problemFields(26) )
        call mxCopyPtrToReal8(mxGetPr(alphaProblemPtr), problem%alpha, sx )
        
        sx = 1
        tolProblemPtr = mxGetField( prhs, i, problemFields(27) )
        call mxCopyPtrToReal8(mxGetPr(tolProblemPtr), problem%tol, sx )
        
        sx = 1
        thres_valueProblemPtr = mxGetField( prhs, i, problemFields(28) )
        call mxCopyPtrToReal8(mxGetPr(thres_valueProblemPtr), problem%thres_value, sx )

        sx = 1
        useCVODEPtr = mxGetField( prhs, i, problemFields(29) )
        call mxCopyPtrToInteger4(mxGetPr(useCVODEPtr), useCVODE, sx )
        if ( useCVODE .eq. 1 ) then
            problem%useCVODE = useCVODETrue
        else
            problem%useCVODE = useCVODEFalse
        endif
        
        sx = 1
        usePassExchPtr = mxGetField( prhs, i, problemFields(50) )
        call mxCopyPtrToInteger4(mxGetPr(usePassExchPtr), passExch, sx )
        if ( passExch .eq. 1 ) then
            problem%passExch = passExchTrue
        else
            problem%passExch = passExchFalse
        endif
        
        
        !File for loading the sparse exchange tensor from Matlab (for non-uniform grids)
        if (( problem%grid%gridType .eq. gridTypeTetrahedron ) .or. (problem%grid%gridType .eq. gridTypeUnstructuredPrisms)) then
            if (problem%passExch .eq. passExchTrue ) then
                sx = 1
                nRowsSparsePtr = mxGetField( prhs, i, problemFields(39) )
                call mxCopyPtrToInteger4(mxGetPr(nRowsSparsePtr), problem%grid%A_exch_load%nrows, sx )
                
                sx = 1
                nColsSparsePtr = mxGetField( prhs, i, problemFields(51) )
                call mxCopyPtrToInteger4(mxGetPr(nColsSparsePtr), problem%grid%A_exch_load%ncols, sx )
                
                sx = 1
                nValuesSparsePtr = mxGetField( prhs, i, problemFields(38) )
                call mxCopyPtrToInteger4(mxGetPr(nValuesSparsePtr), problem%grid%A_exch_load%nvalues, sx )
                
                nvalues = problem%grid%A_exch_load%nvalues
                allocate( problem%grid%A_exch_load%values(nvalues), problem%grid%A_exch_load%rows(nvalues) , problem%grid%A_exch_load%cols(nvalues) )
               
                sx = nvalues
                valuesPtr = mxGetField( prhs, i, problemFields(40) )
                call mxCopyPtrToReal8(mxGetPr(valuesPtr), problem%grid%A_exch_load%values, sx )
            
                sx = nvalues
                rowsPtr = mxGetField( prhs, i, problemFields(41) )
                call mxCopyPtrToInteger4(mxGetPr(rowsPtr), problem%grid%A_exch_load%rows, sx )
        
                sx = nvalues
                colsPtr = mxGetField( prhs, i, problemFields(42) )
                call mxCopyPtrToInteger4(mxGetPr(colsPtr), problem%grid%A_exch_load%cols, sx )       
            endif
        endif
          
        !Load the no. of time steps in the time convergence array
        sx = 1
        ntProblemPtr = mxGetField( prhs, i, problemFields(31) )
        call mxCopyPtrToInteger4(mxGetPr(ntProblemPtr), nt_conv, sx )
        
        allocate( problem%t_conv(nt_conv) )
        genericProblemPtr = mxGetField(prhs,i,problemFields(32) )
        sx = nt_conv
        call mxCopyPtrToReal8(mxGetPr(genericProblemPtr), problem%t_conv, sx )
        
        sx = 1
        genericProblemPtr = mxGetField( prhs, i, problemFields(33) )
        call mxCopyPtrToReal8(mxGetPr(genericProblemPtr), problem%conv_tol, sx )

        sx = 3
        N_aveProblemPtr = mxGetField( prhs, i, problemFields(44) )
        call mxCopyPtrToInteger4(mxGetPr(N_aveProblemPtr), problem%N_ave, sx )
        
        !Coefficient of variation value       
        sx = 1
        CVThresProblemPtr = mxGetField(prhs,i,problemFields(45))
        call mxCopyPtrToReal8(mxGetPr(CVThresProblemPtr), CV, sx )
            
        problem%CV = sngl(CV)
        
        !Parameter to determine if the specific H_fields are returned (exchange, demag, etc.)
        sx = 1
        useReturnHallProblemPtr = mxGetField(prhs,i,problemFields(46))
        call mxCopyPtrToInteger4(mxGetPr(useReturnHallProblemPtr), useReturnHall, sx )
        if ( useReturnHall .eq. 1 ) then
            problem%useReturnHall = useReturnHallTrue
        else
            problem%useReturnHall = useReturnHallFalse
        endif
              
        ! 3x3 matrix specifying the local coordinate system
        !                             [v1_x v2_x v3_x]
        !problem%CrystalAxis(i,:,:) = [v1_y v2_y v3_y]
        !                             [v1_z v2_z v3_z]
        sx = ntot * 3 * 3
        allocate( problem%CrystalAxis(ntot,3,3) )
        CrystalAxisProblemPtr = mxGetField( prhs, i, problemFields(52) )
        call mxCopyPtrToReal8(mxGetPr(CrystalAxisProblemPtr), problem%CrystalAxis, sx )
        
        sx = ntot * 6 * 3
        allocate( problem%K0_arr(ntot,6,3) )
        K0_arrProblemPtr = mxGetField( prhs, i, problemFields(53) )
        call mxCopyPtrToReal8(mxGetPr(K0_arrProblemPtr), problem%K0_arr, sx )
        
        allocate( problem%K1(ntot) )
        sx = ntot
        K1ProblemPtr = mxGetField( prhs, i, problemFields(54) )
        call mxCopyPtrToReal8(mxGetPr(K1ProblemPtr), problem%K1, sx )
        
        allocate( problem%K2(ntot) )
        sx = ntot
        K2ProblemPtr = mxGetField( prhs, i, problemFields(55) )
        call mxCopyPtrToReal8(mxGetPr(K2ProblemPtr), problem%K2, sx )
        
        sx = 1
        exch_weightProblemPtr = mxGetField( prhs, i, problemFields(47) )
        call mxCopyPtrToReal8(mxGetPr(exch_weightProblemPtr), problem%exch_weight, sx )
        
        sx = 1
        exch_methodProblemPtr = mxGetField( prhs, i, problemFields(48) )
        call mxCopyPtrToInteger4(mxGetPr(exch_methodProblemPtr), problem%exch_method, sx )
        
        sx = 1
        exch_interpnProblemPtr = mxGetField( prhs, i, problemFields(49) )
        call mxCopyPtrToInteger4(mxGetPr(exch_interpnProblemPtr), problem%exch_interpn, sx )
        
        !Parameter to determine if the average tensor is used for the prisms
        sx = 1
        useAvgNProblemPtr = mxGetField(prhs,i,problemFields(65))
        call mxCopyPtrToInteger4(mxGetPr(useAvgNProblemPtr), useAvgN, sx )
        if ( useAvgN .eq. 1 ) then
            problem%useAvgN = useAvgNTrue
        else
            problem%useAvgN = useAvgNFalse
        endif
        
        !>-----------------------------------------
        !Calculate the local scaled coefficients for the LLG equation
        !"J" : exchange term
        pi = 3.141592653589793
        mu0 = 4*pi*1e-7
        problem%Jfact = problem%A0 / ( mu0 * problem%Ms )
        !"M" : demagnetization term
        problem%Mfact = problem%Ms
        !"K" : anisotropy term
        problem%Kfact = problem%K0 / ( mu0 * problem%Ms )
        
        !FMM parameters
        sx = 1
        fmm_cellsProblemPtr = mxGetField( prhs, i, problemFields(56) )
        call mxCopyPtrToInteger4(mxGetPr(fmm_cellsProblemPtr), problem%fmm_cells_per_node, sx )
        
        sx = 1
        fmm_epsProblemPtr = mxGetField( prhs, i, problemFields(57) )
        call mxCopyPtrToReal8(mxGetPr(fmm_epsProblemPtr), problem%fmm_eps, sx )
        
        sx = 1
        ifunifProblemPtr = mxGetField( prhs, i, problemFields(58) )
        call mxCopyPtrToInteger4(mxGetPr(ifunifProblemPtr), problem%ifunif, sx )
        
        sx = 1
        nlminProblemPtr = mxGetField( prhs, i, problemFields(59) )
        call mxCopyPtrToInteger4(mxGetPr(nlminProblemPtr), problem%nlmin, sx )
        
        sx = 1
        nlmaxProblemPtr = mxGetField( prhs, i, problemFields(60) )
        call mxCopyPtrToInteger4(mxGetPr(nlmaxProblemPtr), problem%nlmax, sx )
        
        sx = 1
        use_fmmlProblemPtr = mxGetField(prhs,i,problemFields(61))
        call mxCopyPtrToInteger4(mxGetPr(use_fmmlProblemPtr), UseFMM, sx )
        if ( UseFMM .eq. 1 ) then
            problem%use_fmm = useFMMTrue
        else
            problem%use_fmm = useFMMFalse
        endif
        
        sx = 1
        fmm_shortProblemPtr = mxGetField( prhs, i, problemFields(62) )
        call mxCopyPtrToInteger4(mxGetPr(fmm_shortProblemPtr), problem%allow_fmm_short_circuit, sx )
        
        sx = 1
        fmm_min_nProblemPtr = mxGetField( prhs, i, problemFields(63) )
        call mxCopyPtrToInteger4(mxGetPr(fmm_min_nProblemPtr), problem%fmm_min_n, sx )
        
        sx = 1
        useDemagPtr = mxGetField( prhs, i, problemFields(64) )
        call mxCopyPtrToInteger4(mxGetPr(useDemagPtr), useDemag, sx )
        if ( useDemag .eq. 1 ) then
            problem%useDemag = useDemagTrue
        else
            problem%useDemag = useDemagFalse
            call displayGUIMessage( 'NOT using demag field in calculations' )
        endif
        
        allocate( problem%temperature(ntot) )
        sx = ntot
        temperatureProblemPtr = mxGetField( prhs, i, problemFields(66) )
        call mxCopyPtrToReal8(mxGetPr(temperatureProblemPtr), problem%temperature, sx )
       
        sx = 3
        n_macroVecProblemPtr = mxGetField( prhs, i, problemFields(67) )
        call mxCopyPtrToInteger4(mxGetPr(n_macroVecProblemPtr), problem%macrogrid%n_macro, sx )
        
        sx = 3
        shiftVecProblemPtr = mxGetField( prhs, i, problemFields(68) )
        call mxCopyPtrToReal8(mxGetPr(shiftVecProblemPtr), problem%macrogrid%shiftVec, sx )
        
        sx = 3
        macroShapeProblemPtr = mxGetField( prhs, i, problemFields(69) )
        call mxCopyPtrToReal8(mxGetPr(macroShapeProblemPtr), problem%macrogrid%macroShape, sx )
        
        sx = 3
        sampleShapeProblemPtr = mxGetField( prhs, i, problemFields(70) )
        call mxCopyPtrToReal8(mxGetPr(sampleShapeProblemPtr), problem%macrogrid%sampleShape, sx )
        
        ! exchPBC has one entry per direction, like n_macro and shiftVec above. Copying only a
        ! single element left the y and z entries at zero, so from MATLAB the periodic exchange
        ! could only ever be switched on along x.
        sx = 3
        exchPBCProblemPtr = mxGetField( prhs, i, problemFields(71) )
        call mxCopyPtrToInteger4(mxGetPr(exchPBCProblemPtr), problem%macrogrid%exchPBC, sx )
        
        sx = 1
        dummy_runProblemPtr = mxGetField( prhs, i, problemFields(72) )
        call mxCopyPtrToInteger4(mxGetPr(dummy_runProblemPtr), problem%dummy_run, sx )
        
        sx = 1
        fmm_ntermsProblemPtr = mxGetField( prhs, i, problemFields(73) )
        call mxCopyPtrToInteger4(mxGetPr(fmm_ntermsProblemPtr), problem%fmm_nterms, sx )
        
        sx = 1
        window_enaProblemPtr = mxGetField( prhs, i, problemFields(77) )
        call mxCopyPtrToInteger4(mxGetPr(window_enaProblemPtr), problem%window_ena, sx )
        
        sx = 1
        window_intProblemPtr = mxGetField( prhs, i, problemFields(78) )
        call mxCopyPtrToReal8(mxGetPr(window_intProblemPtr), problem%window_int, sx )
        
        sx = 1
        trace_enaProblemPtr = mxGetField( prhs, i, problemFields(79) )
        call mxCopyPtrToInteger4(mxGetPr(trace_enaProblemPtr), problem%trace_ena, sx )
        
        sx = 1
        flush_eachProblemPtr = mxGetField( prhs, i, problemFields(80) )
        call mxCopyPtrToInteger4(mxGetPr(flush_eachProblemPtr), problem%flush_each, sx )
        
        sx = 1
        trace_verbProblemPtr = mxGetField( prhs, i, problemFields(81) )
        call mxCopyPtrToInteger4(mxGetPr(trace_verbProblemPtr), problem%trace_verb, sx )
        
        !flag whether the demag tensor should be loaded
        sx = 1
        N_log_dirPtr = mxGetField( prhs, i, problemFields(82) )
        call mxCopyPtrToInteger4(mxGetPr(N_log_dirPtr), N_log_dir, sx )
        !Length of the file name
        sx = N_log_dir
        log_dirPtr = mxGetField( prhs, i, problemFields(74) )            
        status = mxGetString( log_dirPtr, problem%log_dir, sx )
        
        !flag whether the demag tensor should be loaded
        sx = 1
        N_timer_logPtr = mxGetField( prhs, i, problemFields(83) )
        call mxCopyPtrToInteger4(mxGetPr(N_timer_logPtr), N_timer_log, sx )
        !Length of the file name
        sx = N_timer_log
        timer_logPtr = mxGetField( prhs, i, problemFields(75) )            
        status = mxGetString( timer_logPtr, problem%timer_log, sx )
        
        !flag whether the demag tensor should be loaded
        sx = 1
        N_trace_logPtr = mxGetField( prhs, i, problemFields(84) )
        call mxCopyPtrToInteger4(mxGetPr(N_trace_logPtr), N_trace_log, sx )
        !Length of the file name
        sx = N_trace_log
        trace_logPtr = mxGetField( prhs, i, problemFields(76) )            
        status = mxGetString( trace_logPtr, problem%trace_log, sx )
        
        !Load adaptive hysteresis parameters
        sx = 1
        adaptiveHextPtr = mxGetField( prhs, i, problemFields(85) )
        call mxCopyPtrToInteger4(mxGetPr(adaptiveHextPtr), use_switch_refine, sx )
        problem%adaptiveHext = (use_switch_refine .ne. 0)

        sx = 1
        maxHextStepsPtr = mxGetField( prhs, i, problemFields(86) )
        call mxCopyPtrToInteger4(mxGetPr(maxHextStepsPtr), problem%maxHextSteps, sx )

        sx = 3
        H_startPtr = mxGetField( prhs, i, problemFields(87) )
        call mxCopyPtrToReal8(mxGetPr(H_startPtr), problem%H_start, sx )

        sx = 3
        H_endPtr = mxGetField( prhs, i, problemFields(88) )
        call mxCopyPtrToReal8(mxGetPr(H_endPtr), problem%H_end, sx )

        sx = 1
        dH_initialPtr = mxGetField( prhs, i, problemFields(89) )
        call mxCopyPtrToReal8(mxGetPr(dH_initialPtr), problem%dH_initial, sx )

        sx = 1
        dH_minPtr = mxGetField( prhs, i, problemFields(90) )
        call mxCopyPtrToReal8(mxGetPr(dH_minPtr), problem%dH_min, sx )

        sx = 1
        dH_maxPtr = mxGetField( prhs, i, problemFields(91) )
        call mxCopyPtrToReal8(mxGetPr(dH_maxPtr), problem%dH_max, sx )

        sx = 1
        dH_growPtr = mxGetField( prhs, i, problemFields(92) )
        call mxCopyPtrToReal8(mxGetPr(dH_growPtr), problem%dH_grow, sx )

        sx = 1
        dH_shrinkPtr = mxGetField( prhs, i, problemFields(93) )
        call mxCopyPtrToReal8(mxGetPr(dH_shrinkPtr), problem%dH_shrink, sx )

        sx = 1
        dM_minPtr = mxGetField( prhs, i, problemFields(94) )
        call mxCopyPtrToReal8(mxGetPr(dM_minPtr), problem%dM_min, sx )

        sx = 1
        dM_targetPtr = mxGetField( prhs, i, problemFields(95) )
        call mxCopyPtrToReal8(mxGetPr(dM_targetPtr), problem%dM_target, sx )

        sx = 1
        dM_rejectPtr = mxGetField( prhs, i, problemFields(96) )
        call mxCopyPtrToReal8(mxGetPr(dM_rejectPtr), problem%dM_reject, sx )

        sx = 1
        switch_refine_dHPtr = mxGetField( prhs, i, problemFields(97) )
        call mxCopyPtrToReal8(mxGetPr(switch_refine_dHPtr), problem%switch_refine_dH, sx )

        sx = 1
        use_switch_refinePtr = mxGetField( prhs, i, problemFields(98) )
        call mxCopyPtrToInteger4(mxGetPr(use_switch_refinePtr), use_switch_refine, sx )
        problem%use_switch_refine = (use_switch_refine .ne. 0)

        !Seed for the stochastic thermal field
        sx = 1
        rng_seedPtr = mxGetField( prhs, i, problemFields(99) )
        call mxCopyPtrToInteger4(mxGetPr(rng_seedPtr), problem%rng_seed, sx )

        !Energy minimizer settings. Optional in the struct, so that a problem struct from an older
        !DefaultMicroMagProblem keeps working with the defaults of the Fortran type.
        sx = 1
        minimizerPtr = mxGetField( prhs, i, problemFields(103) )
        if ( minimizerPtr .ne. 0 ) call mxCopyPtrToReal8(mxGetPr(minimizerPtr), problem%min_tol, sx )
        minimizerPtr = mxGetField( prhs, i, problemFields(104) )
        if ( minimizerPtr .ne. 0 ) call mxCopyPtrToInteger4(mxGetPr(minimizerPtr), problem%min_maxiter, sx )
        minimizerPtr = mxGetField( prhs, i, problemFields(105) )
        if ( minimizerPtr .ne. 0 ) call mxCopyPtrToReal8(mxGetPr(minimizerPtr), problem%min_maxrot, sx )
        minimizerPtr = mxGetField( prhs, i, problemFields(106) )
        if ( minimizerPtr .ne. 0 ) call mxCopyPtrToInteger4(mxGetPr(minimizerPtr), problem%min_fallback, sx )
        minimizerPtr = mxGetField( prhs, i, problemFields(107) )
        if ( minimizerPtr .ne. 0 ) call mxCopyPtrToInteger4(mxGetPr(minimizerPtr), problem%min_saddle_check, sx )
        minimizerPtr = mxGetField( prhs, i, problemFields(109) )
        if ( minimizerPtr .ne. 0 ) call mxCopyPtrToInteger4(mxGetPr(minimizerPtr), problem%min_predictor, sx )

        !Timing log switch (optional in the struct; the Fortran default is off)
        minimizerPtr = mxGetField( prhs, i, problemFields(108) )
        if ( minimizerPtr .ne. 0 ) call mxCopyPtrToInteger4(mxGetPr(minimizerPtr), problem%timer_ena, sx )

        !Clean-up
        deallocate(problemFields)
    end subroutine loadMicroMagProblem
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> Returns an array with the names of the fields expected in the MicroMagProblem struct
    !> @param[inout] fieldnames, array of the names of the fields
    !> @param[inout] nfields the no. of elements in fieldnames
    !>-----------------------------------------
    subroutine getProblemFieldnames( fieldnames, nfields)
        integer,intent(out) :: nfields        
        integer,parameter :: nf=109
        character(len=12),dimension(:),intent(out),allocatable :: fieldnames
            
        nfields = nf
        allocate(fieldnames(nfields))
        
        !! Setup the names of the members of the input struct
        fieldnames(1) = 'grid_n'
        fieldnames(2) = 'grid_L'
        fieldnames(3) = 'grid_type'
        fieldnames(4) = 'u_ea'
        fieldnames(5) = 'ProblemMod'
        fieldnames(6) = 'solver'
        fieldnames(7) = 'A0'
        fieldnames(8) = 'Ms'
        fieldnames(9) = 'K0'
        fieldnames(10) = 'gamma'
        fieldnames(11) = 'alpha'
        fieldnames(12) = 'nt_Hext'
        fieldnames(13) = 'Hext'    
        fieldnames(14) = 'nt'
        fieldnames(15) = 't'
        fieldnames(16) = 'm0'
        fieldnames(17) = 'dem_thres'
        fieldnames(18) = 'useCuda'
        fieldnames(19) = 'dem_appr'        
        fieldnames(20) = 'N_ret'
        fieldnames(21) = 'N_file_out'
        fieldnames(22) = 'N_load'
        fieldnames(23) = 'N_file_in'
        fieldnames(24) = 'setTimeDis'
        fieldnames(25) = 'nt_alpha'
        fieldnames(26) = 'alphat'
        fieldnames(27) = 'tol'
        fieldnames(28) = 'thres'
        fieldnames(29) = 'useCVODE'
        fieldnames(30) = 'exch_mat'
        fieldnames(31) = 'nt_conv'
        fieldnames(32) = 't_conv'
        fieldnames(33) = 'conv_tol'
        fieldnames(34) = 'grid_pts'
        fieldnames(35) = 'grid_ele'
        fieldnames(36) = 'grid_nod'
        fieldnames(37) = 'grid_nnod'
        fieldnames(38) = 'exch_nval'
        fieldnames(39) = 'exch_nrow'
        fieldnames(40) = 'exch_val'
        fieldnames(41) = 'exch_rows'
        fieldnames(42) = 'exch_cols'
        fieldnames(43) = 'grid_abc'
        fieldnames(44) = 'N_ave'
        fieldnames(45) = 'CV'
        fieldnames(46) = 'ReturnHall'
        fieldnames(47) = 'exch_weigh'
        fieldnames(48) = 'exch_meth'
        fieldnames(49) = 'exch_intpn'
        fieldnames(50) = 'passExch'
        fieldnames(51) = 'exch_ncol'
        fieldnames(52) = 'CrysAxis'
        fieldnames(53) = 'K0_arr'
        fieldnames(54) = 'K1'
        fieldnames(55) = 'K2'
        fieldnames(56) = 'fmm_cells'
        fieldnames(57) = 'fmm_eps'
        fieldnames(58) = 'ifunif'
        fieldnames(59) = 'nlmin'
        fieldnames(60) = 'nlmax'
        fieldnames(61) = 'use_fmm'
        fieldnames(62) = 'fmm_short'
        fieldnames(63) = 'fmm_min_n'
        fieldnames(64) = 'useDemag'
        fieldnames(65) = 'useAvgN'
        fieldnames(66) = 'temperature'
        fieldnames(67) = 'n_macro'
        fieldnames(68) = 'shiftVec'
        fieldnames(69) = 'macroShape'
        fieldnames(70) = 'sampleShape'
        fieldnames(71) = 'exchPBC'
        fieldnames(72) = 'dummy_run'
        fieldnames(73) = 'fmm_nterms'
        fieldnames(74) = 'log_dir'
        fieldnames(75) = 'timer_log'
        fieldnames(76) = 'trace_log'
        fieldnames(77) = 'window_ena'
        fieldnames(78) = 'window_int'
        fieldnames(79) = 'trace_ena'
        fieldnames(80) = 'flush_each'
        fieldnames(81) = 'trace_verb'
        fieldnames(82) = 'N_log_dir'
        fieldnames(83) = 'N_timer_log'
        fieldnames(84) = 'N_trace_log'
        fieldnames(85) = 'adaptiveHext'
        fieldnames(86) = 'maxHextSteps'
        fieldnames(87) = 'H_start'
        fieldnames(88) = 'H_end'
        fieldnames(89) = 'dH_initial'
        fieldnames(90) = 'dH_min'
        fieldnames(91) = 'dH_max'
        fieldnames(92) = 'dH_grow'
        fieldnames(93) = 'dH_shrink'
        fieldnames(94) = 'dM_min'
        fieldnames(95) = 'dM_target'
        fieldnames(96) = 'dM_reject'
        fieldnames(97) = 'switch_refdH'
        fieldnames(98) = 'use_sw_ref'
        fieldnames(99) = 'rng_seed'

        !Optional exchange at the interface between two materials. A problem struct that
        !predates the feature simply does not have these, which the loader handles.
        fieldnames(100) = 'n_phase'
        fieldnames(101) = 'phase_id'
        fieldnames(102) = 'A_int'

        !Energy minimizer settings, optional as well
        fieldnames(103) = 'min_tol'
        fieldnames(104) = 'min_maxiter'
        fieldnames(105) = 'min_maxrot'
        fieldnames(106) = 'min_fallback'
        fieldnames(107) = 'min_saddle'
        fieldnames(109) = 'min_pred'

        !Timing log switch, optional so that older problem structs keep working
        fieldnames(108) = 'timer_ena'

    end subroutine getProblemFieldnames
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> Returns the solution data struct from Fortran to Matlab
    !> @param[in] solution struct for the internal Fortran represantation of the solution
    !> @param[in] plhs pointer to the Matlab data struct    
    !>-----------------------------------------
    subroutine returnMicroMagSolution( solution, plhs, problem )
        type(MicroMagSolution),intent(in) :: solution           !> Solution to be copied to Matlab        
        type(MicroMagProblem),intent(in) :: problem             !> Problem struct needed for n_Hext_accepted
        mwPointer,intent(inout) :: plhs
    
        integer :: ComplexFlag,classid,mxClassIDFromClassName
        mwSize,dimension(1) :: dims
        mwSize :: s1,s2,sx,ndim
        mwSize,dimension(4) :: dims_4
        mwSize,dimension(3) :: dims_3
        mwPointer :: pt,pm,pp,pdem,pext,pexc,pani,pnHext
        mwPointer :: pE,pnfe,pmit,pmtq,pmst
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

        s1 = 1
        s2 = 1
        pnHext = mxCreateNumericMatrix(s1, s2, mxClassIDFromClassName('int32'), ComplexFlag)
        sx = s1 * s2
        call mxCopyInteger4ToPtr( problem%nHextAccepted, mxGetPr( pnHext ), sx )
        call mxSetField( plhs, ind, fieldnames(8), pnHext )

        !Energies (nt, nt_Hext, 4): exchange, external, demag, anisotropy [J]
        ndim = 3
        dims_3(1) = size( solution%E_out, 1 )
        dims_3(2) = size( solution%E_out, 2 )
        dims_3(3) = size( solution%E_out, 3 )
        classid = mxClassIDFromClassName( 'double' )
        pE = mxCreateNumericArray( ndim, dims_3, classid, ComplexFlag)
        sx = dims_3(1) * dims_3(2) * dims_3(3)
        call mxCopyReal8ToPtr( solution%E_out, mxGetPr( pE ), sx )
        call mxSetField( plhs, ind, fieldnames(9), pE )

        !Relaxation diagnostics, one entry per applied field
        s1 = size( solution%n_feval )
        s2 = 1
        sx = s1 * s2
        pnfe = mxCreateNumericMatrix(s1, s2, mxClassIDFromClassName('int32'), ComplexFlag)
        call mxCopyInteger4ToPtr( solution%n_feval, mxGetPr( pnfe ), sx )
        call mxSetField( plhs, ind, fieldnames(10), pnfe )

        pmit = mxCreateNumericMatrix(s1, s2, mxClassIDFromClassName('int32'), ComplexFlag)
        call mxCopyInteger4ToPtr( solution%min_iter, mxGetPr( pmit ), sx )
        call mxSetField( plhs, ind, fieldnames(11), pmit )

        pmtq = mxCreateDoubleMatrix(s1, s2, ComplexFlag)
        call mxCopyReal8ToPtr( solution%min_torque, mxGetPr( pmtq ), sx )
        call mxSetField( plhs, ind, fieldnames(12), pmtq )

        pmst = mxCreateNumericMatrix(s1, s2, mxClassIDFromClassName('int32'), ComplexFlag)
        call mxCopyInteger4ToPtr( solution%min_status, mxGetPr( pmst ), sx )
        call mxSetField( plhs, ind, fieldnames(13), pmst )

        !Clean up
        deallocate(fieldnames)

    end subroutine returnMicroMagSolution
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> Returns an array with the names of the fields expected in the MicroMagSolution struct
    !> @param[inout] fieldnames, array of the names of the fields
    !> @param[inout] nfields the no. of elements in fieldnames
    !>-----------------------------------------
    subroutine getSolutionFieldnames( fieldnames, nfields)
        integer,intent(out) :: nfields
        integer,parameter :: nf=13
        character(len=10),dimension(:),intent(out),allocatable :: fieldnames

        nfields = nf
        allocate(fieldnames(nfields))

        !! Setup the names of the members of the output struct
        fieldnames(1) = 't'
        fieldnames(2) = 'M'
        fieldnames(3) = 'pts'
        fieldnames(4) = 'H_exc'
        fieldnames(5) = 'H_ext'
        fieldnames(6) = 'H_dem'
        fieldnames(7) = 'H_ani'
        fieldnames(8) = 'n_Hext_acc'
        !Energies and relaxation diagnostics
        fieldnames(9) = 'E'
        fieldnames(10) = 'n_feval'
        fieldnames(11) = 'min_iter'
        fieldnames(12) = 'min_torque'
        fieldnames(13) = 'min_status'

    end subroutine getSolutionFieldnames
    
    
    
    !>-----------------------------------------
    !> @author Rasmus Bj�rk, rabj@dtu.dk, DTU, 2025
    !> Returns the Micromagnetic grid, e.g. GridInfo, from Fortran to Matlab
    !> @param[in] solution struct for the internal Fortran represantation of the solution
    !> @param[in] plhs pointer to the Matlab data struct    
    !>-----------------------------------------
    subroutine returnMicroMagGrid( gridinfo, plhs )
        type(MicroMagGridInfo),intent(in) :: gridinfo           !> The GridInfo is saved in Problem, and is to be copied to Matlab        
        mwPointer,intent(inout) :: plhs
    
        integer :: ComplexFlag,classid,mxClassIDFromClassName
        mwSize,dimension(1) :: dims
        mwSize :: s1,s2,sx,ndim
        mwSize,dimension(4) :: dims_4
        mwPointer :: ptfNormX,ptfNormY,ptfNormZ,ptAreaFaces,ptVolumes,ptXel,ptYel,ptZel
        mwPointer :: ptXf,ptYf,ptZf,ptDimsF,ptTheTs,ptTheDs,ptTheSigns
        mwPointer :: mxCreateStructArray, mxCreateDoubleMatrix,mxGetPr,mxCreateNumericMatrix,mxCreateNumericArray
        mwPointer :: ptExch_mat_r, ptExch_mat_c, ptExch_mat_v, ptExch_mat_nr, ptExch_mat_nc
        mwIndex :: ind
        character(len=10),dimension(:),allocatable :: fieldnames    
        integer :: nfields, nel, nfaces
    
        call getGridInfoFieldnames( fieldnames, nfields)
    
        nfaces = size(gridinfo%fNormX)
        nel = size(gridinfo%Xel)       
                        
        ! Load the result back to Matlab      
        ComplexFlag = 0
      
        dims(1) = 1        
        sx = 1        
        ! Create the return array of structs      
        plhs = mxCreateStructArray( sx, dims, nfields, fieldnames)
      
        ind = 1
        
        s1 = nfaces
        s2 = 1
        ptfNormX = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
        sx = s1 * s2
        call mxCopyReal8ToPtr( gridinfo%fNormX, mxGetPr( ptfNormX ), sx )
        call mxSetField( plhs, ind, fieldnames(1), ptfNormX )
        
        ptfNormY = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
        call mxCopyReal8ToPtr( gridinfo%fNormY, mxGetPr( ptfNormY ), sx )
        call mxSetField( plhs, ind, fieldnames(2), ptfNormY )
        
        ptfNormZ = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
        call mxCopyReal8ToPtr( gridinfo%fNormZ, mxGetPr( ptfNormZ ), sx )
        call mxSetField( plhs, ind, fieldnames(3), ptfNormZ )
        
        ptAreaFaces = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
        call mxCopyReal8ToPtr( gridinfo%AreaFaces, mxGetPr( ptAreaFaces ), sx )
        call mxSetField( plhs, ind, fieldnames(4), ptAreaFaces )
        
        s1 = nel
        s2 = 1
        ptVolumes = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
        sx = s1 * s2
        call mxCopyReal8ToPtr( gridinfo%Volumes, mxGetPr( ptVolumes ), sx )
        call mxSetField( plhs, ind, fieldnames(5), ptVolumes )
        
        ptXel = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
        call mxCopyReal8ToPtr( gridinfo%Xel, mxGetPr( ptXel ), sx )
        call mxSetField( plhs, ind, fieldnames(6), ptXel )
        
        ptYel = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
        call mxCopyReal8ToPtr( gridinfo%Yel, mxGetPr( ptYel ), sx )
        call mxSetField( plhs, ind, fieldnames(7), ptYel )
        
        ptZel = mxCreateDoubleMatrix(s1,s2,ComplexFlag)
        call mxCopyReal8ToPtr( gridinfo%Zel, mxGetPr( ptZel ), sx )
        call mxSetField( plhs, ind, fieldnames(8), ptZel )
                
        s1 = nfaces
        s2 = 1
        ptXf = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
        sx = s1 * s2
        call mxCopyReal8ToPtr( gridinfo%Xf, mxGetPr( ptXf ), sx )
        call mxSetField( plhs, ind, fieldnames(9), ptXf )
        
        ptYf = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
        call mxCopyReal8ToPtr( gridinfo%Yf, mxGetPr( ptYf ), sx )
        call mxSetField( plhs, ind, fieldnames(10), ptYf )
        
        ptZf = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
        call mxCopyReal8ToPtr( gridinfo%Zf, mxGetPr( ptZf ), sx )
        call mxSetField( plhs, ind, fieldnames(11), ptZf )
                
        s1 = nfaces
        s2 = 3
        ptDimsF = mxCreateDoubleMatrix(s1,s2,ComplexFlag)    
        sx = s1 * s2
        call mxCopyReal8ToPtr( gridinfo%DimsF, mxGetPr( ptDimsF ), sx )
        call mxSetField( plhs, ind, fieldnames(12), ptDimsF )
        
        s1 = size(gridinfo%TheTs,1)
        s2 = size(gridinfo%TheTs,2)
        ptTheTs = mxCreateNumericMatrix(s1, s2, mxClassIDFromClassName('int32'), ComplexFlag)
        sx = s1 * s2
        call mxCopyInteger4ToPtr( gridinfo%TheTs, mxGetPr( ptTheTs ), sx )
        call mxSetField( plhs, ind, fieldnames(13), ptTheTs )
        
        s1 = size(gridinfo%TheDs,1)
        s2 = size(gridinfo%TheDs,2)
        ptTheDs = mxCreateNumericMatrix(s1, s2, mxClassIDFromClassName('int32'), ComplexFlag)
        sx = s1 * s2
        call mxCopyInteger4ToPtr( gridinfo%TheDs, mxGetPr( ptTheDs ), sx )
        call mxSetField( plhs, ind, fieldnames(14), ptTheDs )
        
        s1 = size(gridinfo%TheSigns,1)
        s2 = size(gridinfo%TheSigns,2)
        ptTheSigns = mxCreateNumericMatrix(s1, s2, mxClassIDFromClassName('int32'), ComplexFlag)
        sx = s1 * s2
        call mxCopyInteger4ToPtr( gridinfo%TheSigns, mxGetPr( ptTheSigns ), sx )
        call mxSetField( plhs, ind, fieldnames(15), ptTheSigns )
        
        !The COO exchange matrix only exists when the solver built the exchange operator itself;
        !with passExch it came from Matlab and these arrays are unallocated, so return them empty.
        if ( allocated(gridinfo%Exch_mat_r) ) then
            s1 = size(gridinfo%Exch_mat_r,1)
        else
            s1 = 0
        endif
        s2 = 1
        sx = s1 * s2
        ptExch_mat_r = mxCreateNumericMatrix(s1, s2, mxClassIDFromClassName('int32'), ComplexFlag)    
        ptExch_mat_c = mxCreateNumericMatrix(s1, s2, mxClassIDFromClassName('int32'), ComplexFlag)    
        ptExch_mat_v = mxCreateDoubleMatrix(s1, s2, ComplexFlag)    
        if ( s1 .gt. 0 ) then
            call mxCopyInteger4ToPtr( gridinfo%Exch_mat_r, mxGetPr( ptExch_mat_r ), sx )
            call mxCopyInteger4ToPtr( gridinfo%Exch_mat_c, mxGetPr( ptExch_mat_c ), sx )
            call mxCopyReal8ToPtr( gridinfo%Exch_mat_v, mxGetPr( ptExch_mat_v ), sx )
        endif
        call mxSetField( plhs, ind, fieldnames(16), ptExch_mat_r )
        call mxSetField( plhs, ind, fieldnames(17), ptExch_mat_c )
        call mxSetField( plhs, ind, fieldnames(18), ptExch_mat_v )
        
        s1 = 1
        s2 = 1
        sx = s1 * s2
        ptExch_mat_nr = mxCreateNumericMatrix(s1, s2, mxClassIDFromClassName('int32'), ComplexFlag)   
        call mxCopyInteger4ToPtr( gridinfo%Exch_mat_nr, mxGetPr( ptExch_mat_nr ), sx )
        call mxSetField( plhs, ind, fieldnames(19), ptExch_mat_nr )
        
        ptExch_mat_nc = mxCreateNumericMatrix(s1, s2, mxClassIDFromClassName('int32'), ComplexFlag)    
        call mxCopyInteger4ToPtr( gridinfo%Exch_mat_nc, mxGetPr( ptExch_mat_nc ), sx )
        call mxSetField( plhs, ind, fieldnames(20), ptExch_mat_nc )
        
        !Clean up
        deallocate(fieldnames)
    
    end subroutine returnMicroMagGrid
    
    
    !>-----------------------------------------
    !> @author Rasmus Bj�rk, rabj@dtu.dk, DTU, 2025
    !> Returns an array with the names of the fields expected in the GridInfo struct
    !> @param[inout] fieldnames, array of the names of the fields
    !> @param[inout] nfields the no. of elements in fieldnames
    !>-----------------------------------------
    subroutine getGridInfoFieldnames( fieldnames, nfields)
        integer,intent(out) :: nfields
        integer,parameter :: nf=20
        character(len=10),dimension(:),intent(out),allocatable :: fieldnames
            
        nfields = nf
        allocate(fieldnames(nfields))
        
        !! Setup the names of the members of the output struct
        fieldnames(1)  = 'fNormX'
        fieldnames(2)  = 'fNormY'
        fieldnames(3)  = 'fNormZ'
        fieldnames(4)  = 'AreaFaces'
        fieldnames(5)  = 'Volumes'
        fieldnames(6)  = 'Xel'
        fieldnames(7)  = 'Yel'
        fieldnames(8)  = 'Zel'
        fieldnames(9)  = 'Xf'
        fieldnames(10) = 'Yf'
        fieldnames(11) = 'Zf'
        fieldnames(12) = 'DimsF'
        fieldnames(13) = 'TheTs'
        fieldnames(14) = 'TheDs'
        fieldnames(15) = 'TheSigns'
        fieldnames(16) = 'ExchMat_r'
        fieldnames(17) = 'ExchMat_c'
        fieldnames(18) = 'ExchMat_v'
        fieldnames(19) = 'ExchMat_nr'
        fieldnames(20) = 'ExchMat_nc'
        
    end subroutine getGridInfoFieldnames
    
    
    
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
    