include 'mkl_spblas.f90'
include "mkl_dfti.f90"    

    module MicroMagParameters
    use MKL_SPBLAS
    Use MKL_DFTI   
    INTEGER, PARAMETER :: SP = KIND(1.0E0)
    INTEGER, PARAMETER :: DP = KIND(1.0D0)
    
    !>------------------
    !> Custom types
    !>------------------
    !type SparseMatlabMat
    !    integer :: length !> Number of non-zero values
    !    integer :: rows, cols !> Number of rows, cols
    !    integer,dimension(:),allocatable :: ir, jc !> Row- and column index information
    !    real(DP),dimension(:),allocatable :: values
    !end type
    
   !> Stores a table in one variable
    type MicroMagTable1D
        real(DP),dimension(:),allocatable :: x,y        
    end type MicroMagTable1D
    
    !>---------------
    !> Wrapper type for sparse matrices. It appears that Intel's MKL
    !> SPBLAS does not copy the input arrays with the values etc when making a
    !> sparse matrix and we therefore need to store these values close to the sparse
    !> matrix handle    
    type MagTenseSparse
        type(sparse_matrix_t) :: A                                      !> Sparse matrix handle to MKL
        real*4,dimension(:),allocatable :: values                       !> the non-zero values
        integer,dimension(:),allocatable :: rows_start                  !> array of length no. of rows containing the index into values of the first non-zero value in that row
        integer,dimension(:),allocatable :: rows_end                    !> array of length no of rows containing the index into values of the last non-zero value in that row plus one, i.e. the starting value of the next row
        integer,dimension(:),allocatable :: cols                        !> Array of same length as values containing the column no. of the i'th value
        integer :: nvalues                                              !> the number of elements in values
        integer :: nrows                                                !> the number of elements in the row arrays
    end type MagTenseSparse
    
    !Complex version
    type MagTenseSparse_c
        type(sparse_matrix_t) :: A                                      !> Sparse matrix handle to MKL
        complex(kind=4),dimension(:),allocatable :: values              !> the non-zero values
        integer,dimension(:),allocatable :: rows_start                  !> array of length no. of rows containing the index into values of the first non-zero value in that row
        integer,dimension(:),allocatable :: rows_end                    !> array of length no of rows containing the index into values of the last non-zero value in that row plus one, i.e. the starting value of the next row
        integer,dimension(:),allocatable :: cols                        !> Array of same length as values containing the column no. of the i'th value
        integer :: nvalues                                              !> the number of elements in values
        integer :: nrows                                                !> the number of elements in the row arrays
    end type MagTenseSparse_c
    
    !The grid information
    type MicroMagGrid
        integer :: nx, ny, nz
        real(DP) :: Lx,Ly,Lz
        real(DP) :: dx,dy,dz
        real(DP),dimension(:,:,:),allocatable :: x,y,z
        real(DP),dimension(:), allocatable :: dV
        real(DP),dimension(:,:), allocatable :: nodes       !> Arrys with the nodes for a tetrahedron grid
        integer,dimension(:,:), allocatable :: elements     !> Arrys with the elements for a tetrahedron grid, i.e. which nodes belong to which element
        real(DP),dimension(:,:),allocatable :: pts          !> Array with the x,y,z points on list form, i.e. pts(i,:) is the x,y,z components of the i'th point
        integer :: gridType
        integer :: nnodes                                   !> The number of nodes in a tetrahedral grid
        !type(SparseMatlabMat) :: nu_exch_mat                !> Sparse exchange matrix for non-uniform grids (generated in Matlab). Consider moving to problem.
        type(MagTenseSparse) :: A_exch_load                 !> The exchange matrix as read from Matlab
     end type MicroMagGrid
     
    !>-----------------
    !> Overall data structure for a micro magnetism problem.
    !> The design intention is such that a problem may be restarted given the information stored in this struct
    !>-----------------
    type MicroMagProblem
        !Below is stuff that needs to be provided by the "user":
        type(MicroMagGrid) :: grid                  !> Grid of the problem
        
        real(DP),dimension(:,:),allocatable :: u_ea     !> Easy axis vectors that should have the dimensions (n,3) where n is the no. of grid points and thus u_ea(i,3) is the i'th point's z-component
        
        integer :: ProblemMode                      !> Defines the problem mode (new or continued from previous solution)
        
        integer :: solver                           !> Determines what type of solver to use
        
        real(DP) :: A0,Ms,K0,gamma,alpha0,MaxT0         !> User defined coefficients determining part of the problem.
        real(DP) :: tol,thres_value                     !> User defined coefficients for the ODE solver
        
        real(DP),dimension(:,:),allocatable :: Hext     !> Applied field as a function of time. Size (nt,3) with the latter dimension specifying the spatial dimensions.
        real(DP),dimension(:,:),allocatable :: alpha    !> A time dependent dampning parameter, i.e. as a function of time. Size (nt,1).
        
        real(DP),dimension(:),allocatable :: t          !> Time array for the desired output times
        real(DP),dimension(:),allocatable :: m0         !> Initial value of the magnetization
        
        real(DP),dimension(:),allocatable :: t_conv     !> Time array with the time values where the solution will be checked for convergence compared to the last timestep
        real(DP) :: conv_tol                            !> Converge criteria on difference between magnetization at different timesteps
        
        real*4 :: demag_threshold                     !> Used for specifying whether the demag tensors should be converted to sparse matrices by defining values below this value to be zero
        
        integer :: setTimeDisplay                               !> Determines how often the timestep is shown in Matlab
        integer :: useCuda                                      !> Defines whether to attempt using CUDA or not
        integer :: useCVODE                                     !> Defines whether to attempt using CVODE or not
        integer :: demag_approximation                          !> Flag for how to approximate the demagnetization tensor as specified in the parameters below
        integer :: demagTensorReturnState                       !> Flag describing how or if the demag tensor should be returned
        integer :: demagTensorLoadState                         !> Flag describing how or if to load the demag tensor (from disk e.g.)
        character*256 :: demagTensorFileOut, demagTensorFileIn  !> Filename (including path) for output (input) of demag tensor if it is to be returned as a file (demagTensorReturnState >2 and the value is equal to the length of the file including path)
        
        
        !Below is stuff that is computed when the solver initializes
        
        type(sparse_matrix_t) :: A_exch         !> Exchange term matrix
        
        type(MagTenseSparse),dimension(6) :: K_s           !> Sparse matrices (used if the threshold is >0 )
        type(MagTenseSparse_c),dimension(6) :: K_s_c       !> Sparse matrices (used if the threshold is >0 ), complex version
        
        real(SP),dimension(:,:),allocatable :: Kxx,Kxy,Kxz  !> Demag field tensor split out into the nine symmetric components
        real(SP),dimension(:,:),allocatable :: Kyy,Kyz      !> Demag field tensor split out into the nine symmetric components
        real(SP),dimension(:,:),allocatable :: Kzz          !> Demag field tensor split out into the nine symmetric components
        
        real(DP),dimension(:),allocatable :: Axx,Axy,Axz,Ayy,Ayz,Azz    !> Anisotropy vectors assuming local anisotropy only, i.e. no interaction between grains
        
        
        
        type(DFTI_DESCRIPTOR), POINTER :: desc_hndl_FFT_M_H       !> Handle for the FFT MKL stuff
        
    end type MicroMagProblem
    
    !>-----------------
    !> Data structure for a micro magnetism solution.
    !> The design intention is such that a problem may be restarted given the information stored in this struct
    !>-----------------
    type MicroMagSolution
        real*4,dimension(:),allocatable :: HjX,HjY,HjZ                     !> Effective fields for the exchange term (X,Y and Z-directions, respectively)
        real(DP),dimension(:),allocatable :: HhX,HhY,HhZ                   !> Effective fields for the external field (X,Y and Z-directions, respectively)
        real(DP),dimension(:),allocatable :: HkX,HkY,HkZ                   !> Effective fields for the anisotropy energy term (X,Y and Z-directions, respectively)        
        real*4,dimension(:),allocatable :: HmX,HmY,HmZ                     !> Effective fields for the demag energy term (X,Y and Z-directions, respectively)        
        real*4,dimension(:),allocatable :: Mx,My,Mz                        !> The magnetization components used internally as the solution progresses
        complex(kind=4),dimension(:),allocatable :: Mx_FT, My_FT, Mz_FT    !> Fourier transform of Mx, My and Mz (complex)
        complex(kind=4),dimension(:),allocatable :: HmX_c,HmY_c,HmZ_c      !> Complex version of the demag field, used for the Fourier cut-off approach
        
        real(DP),dimension(:),allocatable :: t_out          !> Output times at which the solution was computed
        real(DP),dimension(:,:,:,:),allocatable :: M_out    !> The magnetization at each of these times (nt,ntot,nt_Hext,3)
        real(DP),dimension(:,:,:,:),allocatable :: H_exc    !> The exchange field at each of these times (nt,ntot,nt_Hext,3)
        real(DP),dimension(:,:,:,:),allocatable :: H_ext    !> The external field at each of these times (nt,ntot,nt_Hext,3)
        real(DP),dimension(:,:,:,:),allocatable :: H_dem    !> The demagnetization field at each of these times (nt,ntot,nt_Hext,3)
        real(DP),dimension(:,:,:,:),allocatable :: H_ani    !> The anisotropy field at each of these times (nt,ntot,nt_Hext,3)
        
        real(DP),dimension(:,:),allocatable :: pts          !> n,3 array with the points (x,y,z) of the centers of the tiles
        
        real(DP) :: Jfact,Mfact,Kfact
        
        integer :: HextInd                              !> Index specifying which external field in the input array we have reached in the explicit method
    end type MicroMagSolution
    
    
    !>------------
    !> Parameters
    !>------------
    
    integer,parameter :: gridTypeUniform=1,gridTypeTetrahedron=2
    integer,parameter :: ProblemModeNew=1,ProblemModeContinued=2
    integer,parameter :: MicroMagSolverExplicit=1,MicroMagSolverDynamic=2,MicroMagSolverImplicit=3
    integer,parameter :: useCudaTrue=1,useCudaFalse=0
    integer,parameter :: DemagApproximationNothing=1,DemagApproximationThreshold=2,DemagApproximationFFTThreshold=3,DemagApproximationThresholdFraction=4,DemagApproximationFFTThresholdFraction=5
    integer,parameter :: DemagTensorReturnNot=1,DemagTensorReturnMemory=2
    !!@todo Do NOT have useCVODETrue/-False variables both here and in IntegrationDataTypes.
    integer,parameter :: useCVODETrue=1,useCVODEFalse=0
    
end module MicroMagParameters    