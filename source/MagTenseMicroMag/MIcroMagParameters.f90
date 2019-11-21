include 'mkl_spblas.f90'
    
    module MicroMagParameters
    use MKL_SPBLAS
    
    !>------------------
    !> Custom types
    !>------------------
    
    type MicroMagGrid
        integer :: nx, ny, nz
        real :: Lx,Ly,Lz
        real :: dx,dy,dz
        real,dimension(:,:,:),allocatable :: x,y,z
        real,dimension(:), allocatable :: dV
        
        real,dimension(:,:),allocatable :: pts  !> Array with the x,y,z points on list form, i.e. pts(i,:) is the x,y,z components of the i'th point
        integer :: gridType
    end type MicroMagGrid
    
    !> Stores a table in one variable
    type MicroMagTable1D
        real,dimension(:),allocatable :: x,y        
    end type MicroMagTable1D
    
    !>---------------
    !> Wrapper type for sparse matrices. It appears that Intel's MKL
    !> SPBLAS does not copy the input arrays with the values etc when making a
    !> sparse matrix and we therefore need to store these values close to the sparse
    !> matrix handle    
    type MagTenseSparse
        type(sparse_matrix_t) :: A                                      !> Sparse matrix handle to MKL
        real,dimension(:),allocatable :: values                         !> the non-zero values
        integer,dimension(:),allocatable :: rows_start                  !> array of length no. of rows containing the index into values of the first non-zero value in that row
        integer,dimension(:),allocatable :: rows_end                    !> array of length no of rows containing the index into values of the last non-zero value in that row plus one, i.e. the starting value of the next row
        integer,dimension(:),allocatable :: cols                        !> Array of same length as values containing the column no. of the i'th value
    end type MagTenseSparse
    
    
    
    !>-----------------
    !> Overall data structure for a micro magnetism problem.
    !> The design intention is such that a problem may be restarted given the information stored in this struct
    !>-----------------
    type MicroMagProblem
        !Below is stuff that needs to be provided by the "user":
        type(MicroMagGrid) :: grid                  !> Grid of the problem
        
        real,dimension(:,:),allocatable :: u_ea     !> Easy axis vectors that should have the dimensions (n,3) where n is the no. of grid points and thus u_ea(i,3) is the i'th point's z-component
        
        integer :: ProblemMode                      !> Defines the problem mode (new or continued from previous solution)
        
        integer :: solver                           !> Determines what type of solver to use
        
        real :: A0,Ms,K0,gamma,alpha0,MaxT0         !> User defined coefficients determining part of the problem.
        
        real,dimension(:,:),allocatable :: Hext     !> Applied field as a function of time. Size (nt,3) with the latter dimension specifying the spatial dimensions.
        
        real,dimension(:),allocatable :: t          !> Time array for the desired output times
        real,dimension(:),allocatable :: m0         !>Initial value of the magnetization
        
        real :: demag_threshold                     !> Used for specifying whether the demag tensors should be converted to sparse matrices by defining values below this value to be zero
        
        integer :: setTimeDisplay                   !> Determines how often the timestep is shown in Matlab
        
        !Below is stuff that is computed when the solver initializes
        
        type(sparse_matrix_t) :: A_exch         !> Exchange term matrix
        
        type(MagTenseSparse),dimension(6) :: K_s !> Sparse matrices (used if the threshold is >0 )
        
        real,dimension(:,:),allocatable :: Kxx,Kxy,Kxz  !> Demag field tensor split out into the nine symmetric components
        real,dimension(:,:),allocatable :: Kyy,Kyz      !> Demag field tensor split out into the nine symmetric components
        real,dimension(:,:),allocatable :: Kzz          !> Demag field tensor split out into the nine symmetric components
        
        real,dimension(:),allocatable :: Axx,Axy,Axz,Ayy,Ayz,Azz    !> Anisotropy vectors assuming local anisotropy only, i.e. no interaction between grains
        
        
    end type MicroMagProblem
    
    !>-----------------
    !> Data structure for a micro magnetism solution.
    !> The design intention is such that a problem may be restarted given the information stored in this struct
    !>-----------------
    type MicroMagSolution
        real,dimension(:),allocatable :: HjX,HjY,HjZ    !> Effective fields for the exchange term (X,Y and Z-directions, respectively)
        real,dimension(:),allocatable :: HhX,HhY,HhZ    !> Effective fields for the external field (X,Y and Z-directions, respectively)
        real,dimension(:),allocatable :: HkX,HkY,HkZ    !> Effective fields for the anisotropy energy term (X,Y and Z-directions, respectively)        
        real,dimension(:),allocatable :: HmX,HmY,HmZ    !> Effective fields for the demag energy term (X,Y and Z-directions, respectively)        
        real,dimension(:),allocatable :: Mx,My,Mz       !> The magnetization components used internally as the solution progresses
                        
        real,dimension(:),allocatable :: t_out          !> Output times at which the solution was computed
        real,dimension(:,:,:,:),allocatable :: M_out        !> The magnetization at each of these times (n,3,nt,nt)
        
        real,dimension(:,:),allocatable :: pts          !> n,3 array with the points (x,y,z) of the centers of the tiles
        
        real :: Jfact,Mfact,Kfact                 !> Constant factors used for the determination of the effective fields
        
        integer :: HextInd                              !> Index specifying which external field in the input array we have reached in the explicit method
    end type MicroMagSolution
    
    
    !>------------
    !> Parameters
    !>------------
    
    integer,parameter :: gridTypeUniform=1
    integer,parameter :: ProblemModeNew=1,ProblemModeContinued=2
    integer,parameter :: MicroMagSolverExplicit=1,MicroMagSolverDynamic=2,MicroMagSolverImplicit=3
    
    
end module MicroMagParameters    