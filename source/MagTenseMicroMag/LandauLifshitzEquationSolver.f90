!include 'mkl_spblas.f90'
    module LandauLifshitzSolution
    use ODE_Solvers
    use integrationDataTypes
    use MKL_SPBLAS
    use MicroMagParameters
    use MagTenseMicroMagIO
    !use spline
    use DemagFieldGetSolution
    implicit none
    
   
    
    !>Module variables
    type(MicroMagSolution) :: gb_solution
    type(MicroMagProblem) :: gb_problem
    
    real,dimension(:),allocatable :: crossX,crossY,crossZ   !>Cross product terms
    real,dimension(:),allocatable :: HeffX,HeffY,HeffZ      !>Effective fields
    real,dimension(:),allocatable :: HeffX2,HeffY2,HeffZ2      !>Effective fields
    
    private :: gb_solution,gb_problem,crossX,crossY,crossZ,HeffX,HeffY,HeffZ,HeffX2,HeffY2,HeffZ2
    
    contains
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> @param[inout] prob data structure containing the problem to be solved
    !> @param[inout] sol data structure containing the solution    
    !>-----------------------------------------
    subroutine SolveLandauLifshitzEquation( prob, sol )    
    type(MicroMagProblem),intent(inout) :: prob     !> The problem data structure
    type(MicroMagSolution),intent(inout) :: sol     !> The solution data structure    
    integer :: ntot                                 !> total no. of tiles
    procedure(dydt_fct), pointer :: fct             !> Input function pointer for the function to be integrated
    procedure(callback_fct),pointer :: cb_fct       !> Callback function for displaying progress
    !Save internal representation of the problem and the solution
    gb_solution = sol
    gb_problem = prob
    
    call displayMatlabMessage( 'Initializing matrices' )
    !Calculate the interaction matrices
    call initializeInteractionMatrices( gb_problem )
    
    call displayMatlabMessage( 'Initializing solution' )
    !Initialize the solution, i.e. allocate various arrays
    call initializeSolution( gb_problem, gb_solution )
        
    ntot = gb_problem%grid%nx * gb_problem%grid%ny * gb_problem%grid%nz
    !Set the initial values for m (remember that M is organized such that mx = m(1:ntot), my = m(ntot+1:2*ntot), mz = m(2*ntot+1:3*ntot)
    allocate(gb_solution%t_out(size(gb_problem%t)),gb_solution%M_out(3*ntot,size(gb_problem%t)))
    
    
    call displayMatlabMessage( 'Running solution' )
    !Do the solution
    fct => dmdt_fct
    cb_fct => displayMatlabProgessMessage
    call MagTense_ODE( fct, gb_problem%t, gb_problem%m0, gb_solution%t_out, gb_solution%M_out, cb_fct )
    
    
    !Clean up
    deallocate(crossX,crossY,crossZ,HeffX,HeffY,HeffZ,HeffX2,HeffY2,HeffZ2)
    
    !Make sure to return the correct state
    sol = gb_solution
    prob = gb_problem
    
    end subroutine SolveLandauLifshitzEquation

    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Defines the function that gives the derivative dmdt that is to be integrated
    !> in the Landau-Lifshitz equation. This function only works on a uniform grid
    !> @param[in] t the time at which the derivative is requested
    !> @param[in] m array size n holding the m_i values corresponding to the time t
    !> @param[inout] dmdt array size n for the derivatives at the time t
    !---------------------------------------------------------------------------    
    subroutine dmdt_fct ( t, m, dmdt )  
    real,intent(in) :: t
    real,dimension(:),intent(in) :: m
    real,dimension(:),intent(inout) :: dmdt
    integer :: ntot
    

    
    ntot = gb_problem%grid%nx * gb_problem%grid%ny * gb_problem%grid%nz
    if ( .not. allocated(crossX) ) then
        allocate( crossX(ntot), crossY(ntot), crossZ(ntot) )
        allocate( HeffX(ntot), HeffY(ntot), HeffZ(ntot) )
        allocate( HeffX2(ntot), HeffY2(ntot), HeffZ2(ntot) )
    endif
    
    
    !Update the magnetization
    gb_solution%Mx = m(1:ntot)
    gb_solution%My = m(ntot+1:2*ntot)
    gb_solution%Mz = m(2*ntot+1:3*ntot)
             
    !Exchange term    
    call updateExchangeTerms( gb_problem, gb_solution )
    
    !External field
    call updateExternalField( gb_problem, gb_solution, t )
    
    !Demag. field
    call updateDemagfield( gb_problem, gb_solution )
    
    !Anisotropy term
    call updateAnisotropy(  gb_problem, gb_solution )
    
    !Effective field
    HeffX = gb_solution%HhX + gb_solution%HjX + gb_solution%HmX + gb_solution%HkX
    HeffY = gb_solution%HhY + gb_solution%HjY + gb_solution%HmY + gb_solution%HkY
    HeffZ = gb_solution%HhZ + gb_solution%HjZ + gb_solution%HmZ + gb_solution%HkZ
    
    !Compute m x heff (Precession term)    
    crossX = -1. * ( gb_solution%My * HeffZ - gb_solution%Mz * HeffY )
    crossY = -1. * ( gb_solution%Mz * Heffx - gb_solution%Mx * HeffZ )
    crossZ = -1. * ( gb_solution%Mx * HeffY - gb_solution%My * HeffX )
    
    !Compute m x m x heff (Damping term)
    HeffX2 = gb_solution%My * crossZ - gb_solution%Mz * crossY
    HeffY2 = gb_solution%Mz * crossX - gb_solution%Mx * crossZ
    HeffZ2 = gb_solution%Mx * crossY - gb_solution%My * crossX
    
    !Compute the time derivative of m
    !dMxdt
    dmdt(1:ntot) = alpha(t,gb_problem) * HeffX2 + gb_problem%gamma * crossX
    !dMydt
    dmdt(ntot+1:2*ntot) = alpha(t,gb_problem) * HeffY2 + gb_problem%gamma * crossY
    !dMzdt
    dmdt(2*ntot+1:3*ntot) = alpha(t,gb_problem) * HeffZ2 + gb_problem%gamma * crossZ
    
             
    

    end subroutine dmdt_fct
    
    !>--------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Calulates the alpha-coefficient for the Landau-Lifshitz equation (can be time dependent)
    !> @param[in] t the time at which to evaluate alpha
    !> @param[in] problem the problem on which the solution is solved
    function alpha( t, problem )
    real :: alpha
    real,intent(in) :: t
    type(MicroMagProblem),intent(in) :: problem
    
    alpha = problem%alpha0 * 10**( 7 * min(t,problem%MaxT0)/problem%MaxT0 )
    !MySim.alpha = @(t) -65104e-17*(10.^(7*min(t,MaxT0)/MaxT0)); 
    
    end function alpha
    
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Initializes the solution arrays
    !> @param[in] problem the problem from which the solution is found
    !> @param[inout] solution the solution data structure
    subroutine initializeSolution( problem, solution )
    type(MicroMagProblem),intent(in) :: problem
    type(MicroMagSolution),intent(inout) :: solution
    
    integer :: ntot
    
    if ( problem%problemMode .eq. ProblemModeNew ) then
        !No. of grid points
        ntot = problem%grid%nx * problem%grid%ny * problem%grid%nz
        !Magnetization
        allocate( solution%Mx(ntot), solution%My(ntot), solution%Mz(ntot) )
        solution%Mx(:) = 0.
        solution%My(:) = 0.
        solution%Mz(:) = 0.
        !Exchange effective field
        allocate( solution%HjX(ntot), solution%HjY(ntot), solution%HjZ(ntot) )
        solution%HjX(:) = 0.
        solution%HjY(:) = 0.
        solution%HjZ(:) = 0.
        
        !External effective field
        allocate( solution%HhX(ntot), solution%HhY(ntot), solution%HhZ(ntot) )
        solution%HhX(:) = 0.
        solution%HhY(:) = 0.
        solution%HhZ(:) = 0.
        
        !Anisotropy
        allocate( solution%HkX(ntot), solution%HkY(ntot), solution%HkZ(ntot) )
        solution%HkX(:) = 0.
        solution%HkY(:) = 0.
        solution%HkZ(:) = 0.
        
        !Demag field
        allocate( solution%HmX(ntot), solution%HmY(ntot), solution%HmZ(ntot) )
        solution%HmX(:) = 0.
        solution%HmY(:) = 0.
        solution%HmZ(:) = 0.
        
    endif
    
    !"J" : exchange term
    solution%Jfact = problem%A0 / ( mu0 * problem%Ms )
    !"H" : external field term (b.c. user input is in Tesla)
    solution%Hfact = 1./mu0
    !"M" : demagnetization term
    solution%Mfact = problem%Ms
    !"K" : anisotropy term
    solution%Kfact = problem%K0 / ( mu0 * problem%Ms )
    
    end subroutine initializeSolution
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Calculates and returns the exchange terms
    !> @param[in] problem, the struct containing the current problem
    !> @param[inout] solution, struct containing the current solution    
    !>-----------------------------------------
    subroutine updateExchangeTerms( problem, solution )
    type(MicroMagProblem),intent(in) :: problem
    type(MicroMagSolution),intent(inout) :: solution
    
    integer :: stat
    type(MATRIX_DESCR) :: descr
    real :: alpha
    
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    descr%mode = SPARSE_FILL_MODE_FULL
    descr%diag = SPARSE_DIAG_NON_UNIT
    
    
    
    alpha = 2 * solution%Jfact
    
    !Effective field in the X-direction. Note that the scalar alpha is multiplied on from the left, such that
    !y = alpha * (A_exch * Mx )
    stat = mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%A_exch, descr, solution%Mx, 0., solution%HjX )
    
    !Effective field in the Y-direction
    stat = mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%A_exch, descr, solution%My, 0., solution%HjY )
    
    !Effective field in the Z-direction
    stat = mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, alpha, problem%A_exch, descr, solution%Mz, 0., solution%HjZ )
    
    end subroutine updateExchangeTerms
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Calculates and returns the external field
    !> @param[in] problem, the struct containing the current problem
    !> @param[inout] solution, struct containing the current solution    
    !> @param[in] t the current time
    !>-----------------------------------------
    subroutine updateExternalField( problem, solution, t )
    type(MicroMagProblem),intent(in) :: problem         !> Problem data structure    
    type(MicroMagSolution),intent(inout) :: solution    !> Solution data structure
    real,intent(in) :: t                                !> The time
    
    
    if ( problem%solver .eq. MicroMagSolverExplicit .OR. problem%solver .eq. MicroMagSolverDynamic ) then
        
        
        
        solution%HhX = solution%Hfact * problem%HextX
        solution%HhY = solution%Hfact * problem%HextY
        solution%HhZ = solution%Hfact * problem%HextZ
        
        
    elseif ( problem%solver .eq. MicroMagSolverImplicit ) then
    !not implemented yet
    endif
    
    
    end subroutine updateExternalField
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Calculates and returns the effective field from the anisotropy
    !> @param[in] problem, the struct containing the current problem
    !> @param[inout] solution, struct containing the current solution        
    !>-----------------------------------------
    subroutine updateAnisotropy( problem, solution)
    type(MicroMagProblem),intent(in) :: problem         !> Problem data structure    
    type(MicroMagSolution),intent(inout) :: solution    !> Solution data structure
    
    real :: alpha                                       !> Multiplicative scalar factor
    type(MATRIX_DESCR) :: descr                         !>descriptor for the sparse matrix-vector multiplication
    
    
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    descr%mode = SPARSE_FILL_MODE_FULL
    descr%diag = SPARSE_DIAG_NON_UNIT
    
    
    alpha = -2.*solution%Kfact
    
    !Notice that the anisotropy matrix is symmetric and so Axy = Ayx etc.
    solution%Hkx = alpha * ( problem%Axx * solution%Mx + problem%Axy * solution%My + problem%Axz * solution%Mz )
    solution%Hky = alpha * ( problem%Axy * solution%Mx + problem%Ayy * solution%My + problem%Ayz * solution%Mz )
    solution%Hkz = alpha * ( problem%Axz * solution%Mx + problem%Ayz * solution%My + problem%Azz * solution%Mz )
    

    
    end subroutine updateAnisotropy
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Calculates and returns the effective demag field
    !> @param[in] problem, the struct containing the current problem
    !> @param[inout] solution, struct containing the current solution        
    !>-----------------------------------------
    subroutine updateDemagfield( problem, solution)
    type(MicroMagProblem),intent(in) :: problem         !> Problem data structure    
    type(MicroMagSolution),intent(inout) :: solution    !> Solution data structure
    
    if ( problem%grid%gridType .eq. GridTypeUniform ) then
        call updateDemagfield_uniform( problem, solution)
    endif
    
    
    end subroutine updateDemagfield
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Calculates and returns the effective demag field on a uniform grid
    !> @param[in] problem, the struct containing the current problem
    !> @param[inout] solution, struct containing the current solution        
    !>-----------------------------------------
    subroutine updateDemagfield_uniform( problem, solution)
    type(MicroMagProblem),intent(in) :: problem         !> Problem data structure    
    type(MicroMagSolution),intent(inout) :: solution    !> Solution data structure
    
    !Needs to be checked for proper matrix calculation (Kxx is an n x n matrix while Mx should be n x 1 column vector and the result an n x 1 column vector)
    !Note that the demag tensor is symmetric such that Kxy = Kyx and we only store what is needed.
    solution%HmX = - solution%Mfact * ( matmul( problem%Kxx, solution%Mx ) + matmul( problem%Kxy, solution%My ) + matmul( problem%Kxz, solution%Mz ) )
    solution%HmY = - solution%Mfact * ( matmul( problem%Kxy, solution%Mx ) + matmul( problem%Kyy, solution%My ) + matmul( problem%Kyz, solution%Mz ) )
    solution%HmZ = - solution%Mfact * ( matmul( problem%Kxz, solution%Mx ) + matmul( problem%Kxy, solution%My ) + matmul( problem%Kzz, solution%Mz ) )
    
    
    
    end subroutine updateDemagfield_uniform
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Initializes the interaction matrices
    !> @param[inout] problem struct containing the problem information    
    !---------------------------------------------------------------------------   
    subroutine initializeInteractionMatrices( problem )
    type(MicroMagProblem), intent(inout) :: problem         !> Struct containing the grid information
    
    
    !Setup the grid
    call setupGrid( problem%grid )
    
    !Demagnetization tensor matrix
    call ComputeDemagfieldTensor( problem )
    
    !Anisotropy matrix
    call ComputeAnisotropyTerm3D( problem )
    
    !Exhange term matrix
    call ComputeExchangeTerm3D( problem%grid, problem%A_exch )
    
    
    end subroutine initializeInteractionMatrices
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Sets up the grid for the model. Only uniform grid is supported currently, but this will change in due time
    !> @param[inout] grid the struct that contains the information about the current grid
    !>-----------------------------------------
    
    subroutine setupGrid( grid )
    type(MicroMagGrid),intent(inout) :: grid            !> Grid information to be generated
    integer :: i
    
    !Setup the grid depending on which type of grid it is
    if ( grid%gridType .eq. gridTypeUniform ) then
        
        !Allocate the grid
        allocate( grid%x(grid%nx,grid%ny,grid%nz),grid%y(grid%nx,grid%ny,grid%nz),grid%z(grid%nx,grid%ny,grid%nz) )
        
        if ( grid%nx .gt. 1 ) then
            grid%dx = grid%Lx / (grid%nx-1)                        
            do i=1,grid%nx
                grid%x(i,:,:) = -grid%Lx/2 + (i-1) * grid%dx
            enddo
            
        else
            grid%x(:,:,:) = 0.
            grid%dx = grid%Lx
        endif
    
        if ( grid%ny .gt. 1 ) then
            grid%dy = grid%Ly / (grid%ny-1)                        
            do i=1,grid%ny
                grid%y(:,i,:) = -grid%Ly/2 + (i-1) * grid%dy
            enddo            
        else
            grid%y(:,:,:) = 0.
            grid%dy = grid%Ly
        endif

        if ( grid%nz .gt. 1 ) then
            grid%dz = grid%Lz / (grid%nz-1)                        
            do i=1,grid%nz
                grid%z(:,:,i) = -grid%Lz/2 + (i-1) * grid%dz
            enddo            
        else
            grid%z(:,:,:) = 0.
            grid%dz = grid%Lz
        endif
    endif
    
    
    end subroutine setupGrid
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Calculates and returns the demag field tensor
    !> @param[inout] problem, the struct containing the problem
    
    !>-----------------------------------------
    subroutine ComputeDemagfieldTensor( problem )
    type(MicroMagProblem),intent(inout) :: problem      !> Grid data structure    
    
    
    type(MagTile),dimension(1) :: tile                  !> Tile representing the current tile under consideration
    real,dimension(:,:),allocatable :: H,pts            !> The field and the corresponding evaluation point arrays
    integer :: i,j,k,nx,ny,nz,ntot,ind                  !> Internal counters and index variables
    real,dimension(:,:,:,:),allocatable :: Nout         !> Temporary storage for the demag tensor    
    integer,dimension(1) :: shp
    
    if ( problem%grid%gridType .eq. gridTypeUniform ) then
        nx = problem%grid%nx
        ny = problem%grid%ny
        nz = problem%grid%nz
        ntot = nx * ny * nz
        
        !Demag tensor components
        allocate( problem%Kxx(ntot,ntot), problem%Kxy(ntot,ntot), problem%Kxz(ntot,ntot) )
        allocate( problem%Kyy(ntot,ntot), problem%Kyz(ntot,ntot) )
        allocate( problem%Kzz(ntot,ntot) )
        
        
        allocate(H(ntot,3),pts(ntot,3))
        
        !Setup template tile
        tile(1)%tileType = 2 !(for prism)
        !dimensions of the tile
        tile(1)%a = problem%grid%dx
        tile(1)%b = problem%grid%dy
        tile(1)%c = problem%grid%dz
        tile(1)%exploitSymmetry = 0 !0 for no and this is important
        tile(1)%rotAngles(:) = 0. !ensure that these are indeed zero
        tile(1)%M(:) = 0.
        
        !Set up the points at which the field is to be evaluated (the centers of all the tiles)
        shp(1) = ntot
        pts(:,1) = reshape( problem%grid%x, shp )
        pts(:,2) = reshape( problem%grid%y, shp )
        pts(:,3) = reshape( problem%grid%z, shp )
        
        !for each element find the tensor for all evaluation points (i.e. all elements)
        do k=1,nz
            do j=1,ny                
                do i=1,nx
                    !Set the center of the tile to be the current point
                    tile(1)%offset(1) = problem%grid%x(i,j,k)
                    tile(1)%offset(2) = problem%grid%y(i,j,k)
                    tile(1)%offset(3) = problem%grid%z(i,j,k)
                    !Nout will be allocated by the subroutine. Should be de-allocated afterwards for consistency
                    call getFieldFromTiles( tile, H, pts, 1, ntot, Nout )
                    
                    !Copy Nout into the proper structure used by the micro mag model
                    ind = (i-1) * ny * nz + (j-1) * nz + k
                    
                    problem%Kxx(ind,:) = Nout(1,:,1,1)
                    problem%Kxy(ind,:) = Nout(1,:,1,2)
                    problem%Kxz(ind,:) = Nout(1,:,1,3)
                    
                    !Not stored due to symmetry  (Kxy = Kyx)
                    !Kyx(ind,:) = Nout(1,:,2,1)
                    problem%Kyy(ind,:) = Nout(1,:,2,2)
                    problem%Kyz(ind,:) = Nout(1,:,2,3)
                    
                    !Not stored due to symmetry (Kxz = Kzx)
                    !Kzx(ind,:) = Nout(1,:,3,1)
                    !Not stored due to symmetry (Kyz = Kzy)
                    !Kzy(ind,:) = Nout(1,:,3,2)
                    problem%Kzz(ind,:) = Nout(1,:,3,3)
                    
                    deallocate(Nout)
                enddo
            enddo
        enddo
        
        !Clean up
        deallocate(H,pts)
    endif
    
    
    end subroutine ComputeDemagfieldTensor
    
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Calculates the exhange term matrix
    !> which means it produces the differential operator d^2/dx^2 + d^2/dy^2 + d^2/dz^2 and returns this in the sparse matrix A
    !---------------------------------------------------------------------------   
    subroutine ComputeExchangeTerm3D( grid, A )
    type(MicroMagGrid),intent(in) :: grid             !> Struct containing the grid information    
    type(sparse_matrix_t),intent(inout) :: A          !> The returned matrix from the sparse matrix creator
    
    if ( grid%gridType .eq. gridTypeUniform ) then
        call ComputeExchangeTerm3D_Uniform( grid, A )
    endif
    
    
    end subroutine ComputeExchangeTerm3D
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Calculates the exhange term matrix on a uniform grid
    !> which means it produces the differential operator d^2/dx^2 + d^2/dy^2 + d^2/dz^2 and returns this in the sparse matrix A
    !---------------------------------------------------------------------------   
    subroutine ComputeExchangeTerm3D_Uniform( grid, A )
    type(MicroMagGrid),intent(in) :: grid             !> Struct containing the grid information    
    type(sparse_matrix_t),intent(inout) :: A          !> The returned matrix from the sparse matrix creator
    
    integer :: stat                                   !> Status value for the various sparse matrix operations    
    real,dimension(:),allocatable :: values           !> Array containing the on-zero values of the sparse matrix
    type(sparse_matrix_t) :: d2dx2, d2dy2, d2dz2      !> Sparse matrices for the double derivatives with respect to x, y and z, respectively.
    type(sparse_matrix_t) :: tmp                      !> Temporary sparse matrices used for internal calculations
    integer :: ind, ntot,colInd,rowInd                !> Internal counter for indexing, the total no. of elements in the current sparse matrix being manipulated
    integer,dimension(:),allocatable :: cols          !> Columns array keeping the index of the first element in values of the i'th column of a sparse matrix
    integer,dimension(:),allocatable :: rows_start,rows_end! > rows_start and rows_end are relevant for the sparse matrix definition
    integer :: i,j,k,nx,ny,nz                         !> For-loop counters
    type(matrix_descr) :: descr                         !> Describes a sparse matrix operation
    
    !Find the three sparse matrices for the the individual directions, respectively. Then add them to get the total matrix
    !It is assumed that the magnetization vector to operate on is in fact a single column of Mx, My and Mz respectively.
    
    nx = grid%nx
    ny = grid%ny
    nz = grid%nz
    
    !----------------------------------d^2dx^2 begins -----------------------------!
    !Make the d^2/dx^2 matrix. The no. of non-zero elements is 3 * nx*ny*nz - 2 * ny * nz
    ntot = 3 * nx*ny*nz - 2*ny*nz
    allocate(values(ntot),cols(ntot),rows_start(nx*ny*nz),rows_end(nx*ny*nz))
    
    ind = 1
    rowInd = 1
    do k=1,nz
        do j=1,ny
            
            colInd = (k-1)*ny + (j-1)*nx + 1
            
            !The left boundary
            values(ind) = -1.
            cols(ind) = colInd
            rows_start(rowInd) = ind            
            ind = ind + 1
            
            values(ind) = 1.
            cols(ind) = colInd + 1
            rows_end(rowInd) = ind
            rowInd = rowInd + 1
            ind = ind + 1
            
            !Go through one row at a time
            do i=2,nx-1
                            
                !Left-most point
                values( ind ) = 1.
                cols(ind) = colInd
                !update where the row starts
                rows_start(rowInd) = ind           
                
                ind = ind + 1
                
                !Center point
                values( ind + 1 ) = -2.
                cols(ind) = colInd + 1
                ind = ind + 1
                
                !Right-most point
                values( ind + 2 ) = 1.
                cols(ind) = colInd + 2
                rows_end(rowInd) = ind
                rowInd = rowInd + 1
                ind = ind + 1
                
                colInd = colInd + 1
            enddo            
            !The right boundary
            values(ind) = 1.
            cols(ind) = colInd
            !update where the row starts
            rows_start(rowInd) = ind            
            ind = ind + 1
            values(ind) = -1.
            cols(ind) = colInd+1
            rows_end(rowInd) = ind
            rowInd = rowInd + 1
            ind = ind + 1            
        enddo
    enddo
    !Multiply by the discretization
    values = values * 1./grid%dx**2
        
    !Create the sparse matrix for the d^2dx^2
    stat = mkl_sparse_d_create_csr ( d2dx2, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, rows_start, rows_end, cols, values)
       
    
    !----------------------------------d^2dx^2 ends -----------------------------!
    
    
    !----------------------------------d^2dy^2 begins ----------------------------!
    
    !Make the d^2/dy^2 matrix. The no. of non-zero elements is 3 * nx*ny*nz - 2 * ny * nz just as for d^2dx^2
    values(:) = 0.
    cols(:) = 0
    rows_start(:) = 0
    rows_end(:) = 0
    ind = 1
    rowInd = 1
    do k=1,nz
        !The bottom boundary
        do i=1,nx
            values(ind) = -1.
            cols(ind) = (k-1) * nx * ny + i
            rows_start(rowInd) = ind
            
            !increment to next element
            ind = ind + 1
            
            values(ind) = 1.
            cols(ind) = (k-1) * nx * ny + nx + i
            rows_end(rowInd) = ind
            rowInd = rowInd + 1
            
            !increment to next element
            ind = ind + 1
        enddo
        
        !Everything in between
        do j=2,ny-1
                    
            
            do i=1,nx
                
                values(ind) = 1.
                cols(ind) = (k-1) * nx * ny + (j-2) * nx + i
                rows_start(rowInd) = ind
                !increment to next element
                ind = ind + 1  
                
                values(ind) = -2.
                cols(ind) = (k-1) * nx * ny + (j-1) * nx + i
            
                !increment to next element
                ind = ind + 1
            
                values(ind) = 1.
                cols(ind) = (k-1) * nx * ny + j * nx + i
                rows_end(rowInd) = ind
                rowInd = rowInd + 1
                
                !increment to next element
                ind = ind + 1
                
            enddo
                        
        
        enddo
        
        !The top boundary    
        do i=1,nx
            
            values(ind) = 1.
            cols(ind) = k * nx * ny - 2*nx + i
            rows_start(rowInd) = ind
            
            !increment to next element
            ind = ind + 1            
            
            values(ind) = -1.
            cols(ind) = k * nx * ny - nx + i
            rows_end(rowInd) = ind
            rowInd = rowInd + 1
            
            ind = ind + 1
            
            
        enddo
    enddo
    
    !Multiply by the discretization
    values = values * 1./grid%dy**2
        
    !Create the sparse matrix for the d^2dy^2
    stat = mkl_sparse_d_create_csr ( d2dy2, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, rows_start, rows_end, cols, values)
    
    !----------------------------------d^2dy^2 ends ----------------------------!
    
    !----------------------------------d^2dz^2 begins ----------------------------!
    !Make the d^2/dz^2 matrix. The no. of non-zero elements is 3 * nx*ny*nz - 2 * ny * nz just as for d^2dx^2 and d^2dy^2
    values(:) = 0.
    cols(:) = 0
    rows_start(:) = 0
    rows_end(:) = 0
    ind = 1
    rowInd = 1
    if ( nz .gt. 1 ) then
        !The z=1 face
        do i=1,nx
            do j=1,ny
                !central value
                values(ind) = -1.
                cols(ind) = (j-1) * nx + i
                rows_start(rowInd) = ind
            
                !increment position
                ind = ind + 1
            
                !right-most value
                values(ind) = 1.
                cols(ind) = nx * ny + (j-1)*nx + i
                rows_end(rowInd) = ind
                rowInd = rowInd + 1
            
                !increment position
                ind = ind + 1
            enddo
        enddo
        !Everything in between
        do k=2,nz-1
            do i=1,nx
                do j=1,ny
                    !left-most value
                    values(ind) = 1.
                    cols(ind) = nx * ny * k + (j-1) * nx + i
                    rows_start(rowInd) = ind
            
                    !increment position
                    ind = ind + 1
                
                    !central value
                    values(ind) = -2.
                    cols(ind) = nx * ny * (k-1) + (j-1) * nx + i
                
            
                    !increment position
                    ind = ind + 1
                
                
                    values(ind) = -1.
                    cols(ind) = nx * ny * (k-2) + (j-1) * nx + i
                    rows_end(rowInd) = ind
                    rowInd = rowInd + 1
                
                    !increment position
                    ind = ind + 1
                
                enddo
            enddo
        
        enddo
    
    
        !The z=nz face
        do i=1,nx
            do j=1,ny
            
                !left-most value
                values(ind) = 1.
                cols(ind) = nx * ny * (nz-2) + (j-1) * nx + i
                rows_start(rowInd) = ind
            
                !increment position
                ind = ind + 1
            
                !central value
                values(ind) = -1.
                cols(ind) = nx * ny * (nz-1) + (j-1) * nx + i
                rows_end(rowInd) = ind
                rowInd = rowInd + 1
            
                !increment position
                ind = ind + 1
            
            
            enddo
        enddo
    
    
        !Multiply by the discretization
        values = values * 1./grid%dz**2
        
        !Create the sparse matrix for the d^2dz^2
        stat = mkl_sparse_d_create_csr ( d2dz2, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, rows_start, rows_end, cols, values)
    endif
    
    !----------------------------------d^2dz^2 ends ----------------------------!
    
            
    !Finally, add up the three diagonals and store in the output sparse matrix, A
    !store the results temporarily in tmp4
    stat = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, d2dx2, 1., d2dy2, tmp)    
    if ( nz .gt. 1 ) then    
        stat = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, d2dz2, 1., tmp, A)
        !clean up
        stat = mkl_sparse_destroy (d2dz2)    
    else
        descr%type = SPARSE_MATRIX_TYPE_GENERAL
        descr%mode = SPARSE_FILL_MODE_FULL
        descr%diag = SPARSE_DIAG_NON_UNIT
        stat = mkl_sparse_copy ( tmp, descr, A )
    endif
    
    !clean up
    deallocate(values,cols,rows_start,rows_end)
    
    stat = mkl_sparse_destroy (d2dx2)
    stat = mkl_sparse_destroy (d2dy2)
    
    
    end subroutine ComputeExchangeTerm3D_Uniform
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Calculates the anisotropy term sparse matrix assuming the effective field anisotropy is linear in m    
    !> @param[inout] problem the data structure containing the problem
    !---------------------------------------------------------------------------   
    subroutine ComputeAnisotropyTerm3D( problem )
    type(MicroMagProblem),intent(inout) :: problem       !> Struct containing the problem
    
    
    if ( problem%grid%gridType .eq. gridTypeUniform ) then
        call ComputeAnisotropyTerm3D_Uniform( problem )
    endif
    
    
    end subroutine ComputeAnisotropyTerm3D
    
    !>-----------------------------------------
    !> @author Kaspar K. Nielsen, kaki@dtu.dk, DTU, 2019
    !> @brief
    !> Calculates the anisotropy term matrix on a uniform grid  
    !> @param[inout] problem the data structure containing the problem
    !---------------------------------------------------------------------------   
    subroutine ComputeAnisotropyTerm3D_Uniform( problem )
    type(MicroMagProblem),intent(inout) :: problem             !> Struct containing the problem
    
    
    integer :: nx,ny,nz,ntot
    
        nx = problem%grid%nx
        ny = problem%grid%ny
        nz = problem%grid%nz
        ntot = nx * ny * nz
        
        !Allocate the anisotropy vectors and note that the operation is symmetric so we do not need to store three of the nine components
        allocate( problem%Axx(ntot),problem%Axy(ntot),problem%Axz(ntot),problem%Ayy(ntot),problem%Ayz(ntot),problem%Azz(ntot) )
        
        problem%Axx = problem%u_ea(:,1) * problem%u_ea(:,1)
        problem%Axy = problem%u_ea(:,1) * problem%u_ea(:,2)
        problem%Axz = problem%u_ea(:,1) * problem%u_ea(:,3)
        problem%Ayy = problem%u_ea(:,2) * problem%u_ea(:,2)
        problem%Ayz = problem%u_ea(:,2) * problem%u_ea(:,3)
        problem%Azz = problem%u_ea(:,3) * problem%u_ea(:,3)
    
    end subroutine ComputeAnisotropyTerm3D_Uniform
    
end module LandauLifshitzSolution
    