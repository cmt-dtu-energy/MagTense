#include "fintrf.h"
    
    module MagTenseMicroMagIO
    use MicroMagParameters
        
    implicit none
    
    contains
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
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
        integer :: nFieldsProblem,ntot,nt, nt_Hext,useCuda
        mwPointer :: nGridPtr,LGridPtr,dGridPtr,typeGridPtr, ueaProblemPtr, modeProblemPtr,solverProblemPtr
        mwPointer :: A0ProblemPtr,MsProblemPtr,K0ProblemPtr,gammaProblemPtr,alpha0ProblemPtr,MaxT0ProblemPtr
        mwPointer :: ntProblemPtr, m0ProblemPtr,HextProblemPtr,tProblemPtr,useCudaPtr
        mwPointer :: mxGetField, mxGetPr, ntHextProblemPtr,demThresProblemPtr,demApproxPtr
        integer,dimension(3) :: int_arr
        real,dimension(3) :: real_arr
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
        
        
        !Clean-up 
        deallocate(problemFields)
    end subroutine loadMicroMagProblem
    
    
   !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
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
        mwPointer :: pt,pm,pp
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
        
        !Clean up
        deallocate(fieldnames)
    
    end subroutine returnMicroMagSolution
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> Returns an array with the names of the fields expected in the MicroMagProblem struct
    !> @param[inout] fieldnames, array of the names of the fields
    !> @param[inout] nfields the no. of elements in fieldnames
    !>-----------------------------------------
    subroutine getProblemFieldnames( fieldnames, nfields)
        integer,intent(out) :: nfields
        integer,parameter :: nf=20
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
        
        
    end subroutine getProblemFieldnames
    
   
    
     !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> Returns an array with the names of the fields expected in the MicroMagSolution struct
    !> @param[inout] fieldnames, array of the names of the fields
    !> @param[inout] nfields the no. of elements in fieldnames
    !>-----------------------------------------
    subroutine getSolutionFieldnames( fieldnames, nfields)
        integer,intent(out) :: nfields
        integer,parameter :: nf=3
        character(len=10),dimension(:),intent(out),allocatable :: fieldnames
            
        nfields = nf
        allocate(fieldnames(nfields))
        
        !! Setup the names of the members of the input struct
        fieldnames(1) = 't'
        fieldnames(2) = 'M'
        fieldnames(3) = 'pts'
        
        
        
    end subroutine getSolutionFieldnames
    
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
    
    subroutine displayMatlabProgessMessage( mess, prog )
        character(*),intent(in) :: mess
        integer,intent(in) :: prog
        character*(4) :: prog_str   
        character(len=8) :: fmt
        
    
        fmt = '(I4.2)'
        write (prog_str,fmt) prog
        
        call displayMatlabMessage( mess )
        call displayMatlabMessage( prog_str )
        
    end subroutine displayMatlabProgessMessage
    
    subroutine displayMatlabProgessTime( mess, time  )
        character(*),intent(in) :: mess
        real,intent(in) :: time
        character*(4) :: prog_str   
                
            
        write (prog_str,'(F4.2)') time
        
        call displayMatlabMessage( mess )
        call displayMatlabMessage( prog_str )
        
    end subroutine displayMatlabProgessTime
    
    end module MagTenseMicroMagIO
    