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
        integer :: nFieldsProblem, ntot, nt, nt_Hext, useCuda, status, nt_alpha, useCVODE, nt_conv
        mwPointer :: nGridPtr, LGridPtr, dGridPtr, typeGridPtr, ueaProblemPtr, modeProblemPtr, solverProblemPtr
        mwPointer :: A0ProblemPtr, MsProblemPtr, K0ProblemPtr, gammaProblemPtr, alpha0ProblemPtr, MaxT0ProblemPtr
        mwPointer :: ntProblemPtr, m0ProblemPtr, HextProblemPtr, tProblemPtr, useCudaPtr, useCVODEPtr
        mwPointer :: mxGetField, mxGetPr, mxGetM, mxGetN, mxGetNzmax, mxGetIr, mxGetJc
        mwPointer :: ntHextProblemPtr, demThresProblemPtr, demApproxPtr, setTimeDisplayProblemPtr
        mwPointer :: NFileReturnPtr, NReturnPtr, NLoadPtr, mxGetString, NFileLoadPtr
        mwPointer :: tolProblemPtr, thres_valueProblemPtr
        mwPointer :: exch_matProblemPtr, irPtr, jcPtr
        mwPointer :: genericProblemPtr
        integer,dimension(3) :: int_arr
        real*8,dimension(3) :: real_arr
        real*8 :: demag_fac
        
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
        
        sx = 1
        MsProblemPtr = mxGetField( prhs, i, problemFields(8) )
        call mxCopyPtrToReal8(mxGetPr(MsProblemPtr), problem%Ms, sx )
        
        sx = 1
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
        HextProblemPtr = mxGetField( prhs, i, problemFields(27) )
        call mxCopyPtrToReal8(mxGetPr(HextProblemPtr), problem%alpha, sx )
        
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
        if ( problem%grid%gridType .eq. gridTypeNonUniform ) then
            ! Number of non-zero entries in sparse matrix
            exch_matProblemPtr = mxGetField( prhs, i, problemFields(31) )
            problem%grid%nu_exch_mat%length = mxGetNzmax(exch_matProblemPtr)
            
            problem%grid%nu_exch_mat%cols = mxGetN(exch_matProblemPtr)
            problem%grid%nu_exch_mat%rows = mxGetM(exch_matProblemPtr)
            
            ! Row indices for elements
            sx = problem%grid%nu_exch_mat%length
            allocate( problem%grid%nu_exch_mat%ir( problem%grid%nu_exch_mat%length ) )
            irPtr = mxGetIr(exch_matProblemPtr)
            call mxCopyPtrToInteger4(mxGetPr(irPtr), problem%grid%nu_exch_mat%ir, sx)
            
            ! Column index information
            sx = problem%grid%nu_exch_mat%cols + 1
            allocate( problem%grid%nu_exch_mat%jc( sx ) )
            jcPtr = mxGetJc(exch_matProblemPtr)
            call mxCopyPtrToInteger4(mxGetPr(jcPtr), problem%grid%nu_exch_mat%jc, sx)
            
            sx = problem%grid%nu_exch_mat%length
            allocate( problem%grid%nu_exch_mat%values( problem%grid%nu_exch_mat%length ) )
            call mxCopyPtrToReal8(mxGetPr(exch_matProblemPtr), problem%grid%nu_exch_mat%values, sx )
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
        dims_4(1) = nt
        dims_4(2) = ntot
        dims_4(3) = size( solution%H_exc(1,1,:,1) )
        dims_4(4) = 3
        classid = mxClassIDFromClassName( 'double' )
        pexc = mxCreateNumericArray( ndim, dims_4, classid, ComplexFlag)
        sx = dims_4(1) * dims_4(2) * dims_4(3) * dims_4(4)
        call mxCopyReal8ToPtr( solution%H_exc, mxGetPr( pexc ), sx )
        call mxSetField( plhs, ind, fieldnames(4), pexc )
        
        
        
        ndim = 4
        dims_4(1) = nt
        dims_4(2) = ntot
        dims_4(3) = size( solution%H_ext(1,1,:,1) )
        dims_4(4) = 3
        classid = mxClassIDFromClassName( 'double' )
        pext = mxCreateNumericArray( ndim, dims_4, classid, ComplexFlag)
        sx = dims_4(1) * dims_4(2) * dims_4(3) * dims_4(4)
        call mxCopyReal8ToPtr( solution%H_ext, mxGetPr( pext ), sx )
        call mxSetField( plhs, ind, fieldnames(5), pext )
        
        
        
        ndim = 4
        dims_4(1) = nt
        dims_4(2) = ntot
        dims_4(3) = size( solution%H_dem(1,1,:,1) )
        dims_4(4) = 3
        classid = mxClassIDFromClassName( 'double' )
        pdem = mxCreateNumericArray( ndim, dims_4, classid, ComplexFlag)
        sx = dims_4(1) * dims_4(2) * dims_4(3) * dims_4(4)
        call mxCopyReal8ToPtr( solution%H_dem, mxGetPr( pdem ), sx )
        call mxSetField( plhs, ind, fieldnames(6), pdem )
        
        
        
        ndim = 4
        dims_4(1) = nt
        dims_4(2) = ntot
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
        integer,parameter :: nf=34
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
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !! Routine for displaying progress within Matlab
    subroutine displayMatlabMessage( mess )
        character(*),intent(in) :: mess
    
        integer :: mexCallMATLAB, nlhs_cb, nrhs_cb, tmp
        mwPointer plhs_cb(1), prhs_cb(1),mxCreateString
        character*(4) functionName_cb
        logical :: ex
    
        
        !! test if we are inside Matlab or called from a stand-alone. The simple test
        !! is whether io.txt exists in the current path (false for Matlab, true for stand-alone)
        !! nothing bad should happen if in fact we are called from ML but the file somehow
        !! exists - the written output will just not be shown to the user
    
        inquire( file='io.txt', EXIST=ex )
    
        if ( ex .eq. .true. ) then
            write(*,*) mess
        else            
            functionName_cb = "disp"
            nlhs_cb = 0
            nrhs_cb = 1
      
            prhs_cb(1) = mxCreateString(mess)
      
            tmp = mexCallMATLAB(nlhs_cb, plhs_cb, nrhs_cb, prhs_cb, "disp")
        endif
    
    end subroutine displayMatlabMessage
    
    subroutine displayMatlabProgressMessage( mess, prog )
        character(*),intent(in) :: mess
        integer,intent(in) :: prog
        character*(4) :: prog_str   
        character(len=8) :: fmt
        
    
        fmt = '(I4.2)'
        write (prog_str,fmt) prog
        
        call displayMatlabMessage( mess )
        call displayMatlabMessage( prog_str )
        
    end subroutine displayMatlabProgressMessage
    
    subroutine displayMatlabProgessTime( mess, time  )
        character(*),intent(in) :: mess
        real,intent(in) :: time
        character*(4) :: prog_str   
                
            
        write (prog_str,'(F4.2)') time
        
        call displayMatlabMessage( mess )
        call displayMatlabMessage( prog_str )
        
    end subroutine displayMatlabProgessTime
    
    end module MagTenseMicroMagIO
    