module MagTenseMicroMagPyIO
use MicroMagParameters
    
implicit none

contains


subroutine loadMicroMagProblem( ntot, grid_n, grid_L, grid_type, u_ea, ProblemMode, solver, A0, Ms, K0, &
    gamma, alpha, MaxT0, nt_Hext, Hext, nt, t, m0, dem_thres, useCuda, dem_appr, N_ret, N_file_out, &
    N_load, N_file_in, setTimeDis, nt_alpha, alphat, tol, thres, useCVODE, nt_conv, t_conv, &
    conv_tol, grid_pts, grid_ele, grid_nod, grid_nnod, exch_nval, exch_nrow, exch_val, exch_rows, &
    exch_rowe, exch_col, grid_abc, usePrecision, nThreadsMatlab, N_ave, problem )
    !DEC$ ATTRIBUTES ALIAS:"loadmicromagproblem_" :: loadMicroMagProblem
    integer(4), intent(in) :: ntot, nt_conv, grid_type, nt_Hext, nt_alpha, nt, grid_nnod, exch_nval, exch_nrow
    integer(4),dimension(3),intent(in) :: grid_n
    real(8),dimension(3),intent(in) :: grid_L
    real(8),dimension(ntot, 3),intent(in) :: grid_pts
    integer(4),dimension(4, ntot),intent(in) :: grid_ele
    real(8),dimension(grid_nnod, 3),intent(in) :: grid_nod
    real(8),dimension(ntot, 3),intent(in) :: grid_abc
    real(8),dimension(ntot, 3),intent(in) :: u_ea
    real(8),dimension(nt_Hext, 4),intent(in) :: Hext
    real(8),dimension(nt),intent(in) :: t
    real(8),dimension(3*ntot),intent(in) :: m0
    real(8),dimension(nt_alpha,2),intent(in) :: alphat
    integer(4),dimension(exch_nval),intent(in) :: exch_val
    integer(4),dimension(exch_nrow),intent(in) :: exch_rows, exch_rowe
    integer(4),dimension(exch_nval),intent(in) :: exch_col
    real(8),dimension(nt_conv),intent(in) :: t_conv
    integer(4),intent(in) :: ProblemMode, solver, useCuda, dem_appr, usePrecision, nThreadsMatlab
    integer(4),intent(in) :: N_ret, N_load, setTimeDis, useCVODE
    real(8),intent(in) :: A0, Ms, K0, gamma, alpha, MaxT0, tol, thres, conv_tol, dem_thres
    integer(4), dimension(3) :: N_ave
    real(8) :: demag_fac
    character*256,intent(in) :: N_file_in, N_file_out

    type(MicroMagProblem),intent(inout) :: problem

    problem%grid%nx = grid_n(1)
    problem%grid%ny = grid_n(2)
    problem%grid%nz = grid_n(3)

    problem%grid%Lx = grid_L(1)
    problem%grid%Ly = grid_L(2)
    problem%grid%Lz = grid_L(3)

    problem%grid%dx = problem%grid%Lx / problem%grid%nx
    problem%grid%dy = problem%grid%Ly / problem%grid%ny
    problem%grid%dz = problem%grid%Lz / problem%grid%nz

    problem%grid%gridType = grid_type

    !Load additional things for a tetrahedron grid
    if ( problem%grid%gridType .eq. gridTypeTetrahedron ) then
        !The center points of all the tetrahedron elements           
        allocate( problem%grid%pts(ntot,3) )
        problem%grid%pts = grid_pts
        
        !The elements of all the tetrahedron elements
        allocate( problem%grid%elements(4,ntot) )
        problem%grid%elements = grid_ele
        
        !The number of nodes in the tetrahedron mesh
        problem%grid%nnodes = grid_nnod
        
        !The nodes of all the tetrahedron elements
        allocate( problem%grid%nodes(3,grid_nnod) )
        problem%grid%nodes = grid_nod
        
        !the number of nodes in the tetrahedron mesh
        problem%grid%nnodes = grid_nnod
    endif
    
    !Load additional things for a grid of unstructured prisms
    if ( problem%grid%gridType .eq. gridTypeUnstructuredPrisms ) then
        !The center points of all the prisms elements           
        allocate( problem%grid%pts(ntot,3) )
        problem%grid%pts = grid_pts
        
        !The side lengths of all the prisms
        allocate( problem%grid%abc(ntot,3) )
        problem%grid%abc = grid_abc
        
    endif

    !Finished loading the grid------------------------------------------
    
    !Start loading the problem
    !Allocate memory for the easy axis vectors
    allocate( problem%u_ea(ntot,3) )
    problem%u_ea = u_ea
    
    problem%ProblemMode = ProblemMode
    problem%solver = solver
    problem%A0 = A0
    allocate( problem%Ms(ntot) )
    problem%Ms = Ms
    allocate( problem%K0(ntot) )
    problem%K0 = K0
    problem%gamma = gamma
    problem%alpha0 = alpha
    problem%MaxT0 = MaxT0
    
    !Applied field as a function of time evaluated at the timesteps specified in nt_Hext
    !problem%Hext(:,1) is the time grid while problem%Hext(:,2:4) are the x-,y- and z-components of the applied field
    problem%Hext = Hext
    
    allocate( problem%t(nt) )
    problem%t = t
    
    !Initial magnetization
    allocate( problem%m0(3*ntot) )
    problem%m0 = m0
    
    !Demagnetization threshold value
    demag_fac = dem_thres
    problem%demag_threshold = sngl(demag_fac)
    
    if ( useCuda .eq. 1 ) then
        problem%useCuda = useCudaTrue
    else
        problem%useCuda = useCudaFalse
    endif
            
    problem%demag_approximation = dem_appr
    
    !flag whether the demag tensor should be returned and if so how
    problem%demagTensorReturnState = N_ret
    
    !File for returning the demag tensor to a file on disk (has to have length>2)
    if ( problem%demagTensorReturnState .gt. 2 ) then
        !Length of the file name
        problem%demagTensorFileOut = N_file_out
    endif
    
    !flag whether the demag tensor should be loaded
    problem%demagTensorLoadState = N_load
    
    !File for loading the demag tensor to a file on disk (has to have length>2)
    if ( problem%demagTensorLoadState .gt. 2 ) then
        !Length of the file name
        problem%demagTensorFileIn = N_file_in
    endif
    
    problem%setTimeDisplay = 100
    
    !Set how often to display the timestep in Matlab
    problem%setTimeDisplay = setTimeDis
    
    !alpha as a function of time evaluated at the timesteps
    !problem%alpha(:,1) is the time grid while problem%alpha(:,2) are the alpha values
    problem%alpha = alphat
    
    problem%tol = tol
    problem%thres_value = thres

    if ( useCVODE .eq. 1 ) then
        problem%useCVODE = useCVODETrue
    else
        problem%useCVODE = useCVODEFalse
    endif

    !File for loading the sparse exchange tensor from Matlab (for non-uniform grids)
    if (( problem%grid%gridType .eq. gridTypeTetrahedron ) .or. (problem%grid%gridType .eq. gridTypeUnstructuredPrisms)) then
        ! Load the CSR sparse information from Matlab
        problem%grid%A_exch_load%nvalues = exch_nval
        problem%grid%A_exch_load%nrows = exch_nrow
        
        allocate( problem%grid%A_exch_load%values(exch_nval) )
        allocate( problem%grid%A_exch_load%rows_start(exch_nrow) )
        allocate( problem%grid%A_exch_load%rows_end(exch_nrow) )
        allocate( problem%grid%A_exch_load%cols(exch_nval) )
            
        problem%grid%A_exch_load%values = exch_val
        problem%grid%A_exch_load%rows_start = exch_rows
        problem%grid%A_exch_load%rows_end = exch_rowe
        problem%grid%A_exch_load%cols = exch_col
    endif
        
    !Load the no. of time steps in the time convergence array        
    allocate( problem%t_conv(nt_conv) )
    problem%t_conv = t_conv
    problem%conv_tol = conv_tol

    if ( usePrecision .eq. 1 ) then
        problem%usePrecision = usePrecisionTrue
    else
        problem%usePrecision = usePrecisionFalse
    endif
    
    problem%nThreadsMatlab = nThreadsMatlab
    problem%N_ave = N_ave

end subroutine loadMicroMagProblem


subroutine returnMicroMagSolutionPy( solution, nt, ntot, ndim, t, M, pts, H_exc, H_ext, H_dem, H_ani )
    type(MicroMagSolution),intent(in) :: solution
    real(8),dimension(nt),intent(out) :: t
    real(8),dimension(nt,ntot,ndim,3),intent(out) :: M, H_exc, H_ext, H_dem, H_ani
    real(8),dimension(ntot,nt),intent(out) :: pts
    integer,intent(in) :: ntot, nt, ndim

    t = solution%t_out
    M = solution%M_out
    pts = solution%pts
    H_exc = solution%H_exc
    H_ext = solution%H_ext
    H_dem = solution%H_dem
    H_ani = solution%H_ani

end subroutine returnMicroMagSolutionPy


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

end module MagTenseMicroMagPyIO
