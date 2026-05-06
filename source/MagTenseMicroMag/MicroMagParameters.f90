include 'mkl_spblas.f90'
include "mkl_dfti.f90"

    module MicroMagParameters
    use MKL_SPBLAS
    Use MKL_DFTI
#if USE_FMM3D
    use fmm3d_tree_mod
#endif
    INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(6, 37)
    INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(15, 307)
    
    !>------------------
    !> Custom types
    !>------------------
       
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
        real(SP),dimension(:),allocatable :: values                     !> the non-zero values
        integer,dimension(:),allocatable :: rows_start                  !> array of length no. of rows containing the index into values of the first non-zero value in that row
        integer,dimension(:),allocatable :: rows_end                    !> array of length no of rows containing the index into values of the last non-zero value in that row plus one, i.e. the starting value of the next row
        integer,dimension(:),allocatable :: cols                        !> Array of same length as values containing the column no. of the i'th value
        integer :: nvalues                                              !> the number of elements in values
        integer :: nrows                                                !> the number of elements in the row arrays
        integer :: ncols                                                !> the number of columns
    end type MagTenseSparse
    
    !Double version
    type MagTenseSparse_d
        type(sparse_matrix_t) :: A                                      !> Sparse matrix handle to MKL
        real(DP),dimension(:),allocatable :: values                     !> the non-zero values
        integer,dimension(:),allocatable :: rows_start                  !> array of length no. of rows containing the index into values of the first non-zero value in that row
        integer,dimension(:),allocatable :: rows_end                    !> array of length no of rows containing the index into values of the last non-zero value in that row plus one, i.e. the starting value of the next row
        integer,dimension(:),allocatable :: cols                        !> Array of same length as values containing the column no. of the i'th value
        integer :: nvalues                                              !> the number of elements in values
        integer :: nrows                                                !> the number of elements in the row arrays
        integer :: ncols                                                !> the number of columns
    end type MagTenseSparse_d
    
    !Complex version
    type MagTenseSparse_c
        type(sparse_matrix_t) :: A                                      !> Sparse matrix handle to MKL
        complex(kind=4),dimension(:),allocatable :: values              !> the non-zero values
        integer,dimension(:),allocatable :: rows_start                  !> array of length no. of rows containing the index into values of the first non-zero value in that row
        integer,dimension(:),allocatable :: rows_end                    !> array of length no of rows containing the index into values of the last non-zero value in that row plus one, i.e. the starting value of the next row
        integer,dimension(:),allocatable :: cols                        !> Array of same length as values containing the column no. of the i'th value
        integer :: nvalues                                              !> the number of elements in values
        integer :: nrows                                                !> the number of elements in the row arrays
        integer :: ncols                                                !> the number of columns
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
        real(DP),dimension(:,:),allocatable :: abc          !> Array with the side lengths a,b,c in list form
        integer :: gridType
        integer :: nnodes                                   !> The number of nodes in a tetrahedral grid
        type(MagTenseSparse_d) :: A_exch_load                 !> The exchange matrix as read from Matlab
    end type MicroMagGrid
    
    !Additional information about the grid
    type MicroMagGridInfo
        real(dp), allocatable :: fNormX(:)
        real(dp), allocatable :: fNormY(:)
        real(dp), allocatable :: fNormZ(:)
        real(dp), allocatable :: AreaFaces(:)
        real(dp), allocatable :: Volumes(:)
        integer, allocatable  :: TheTs(:,:)                                      !> Indices that can be supplied to create sparse matrix handle to MKL
        integer, allocatable  :: TheDs(:,:)                                      !> Indices that can be supplied to create sparse matrix handle to MKL
        integer, allocatable  :: TheSigns(:,:)                                   !> Indices that can be supplied to create sparse matrix handle to MKL
        real(dp), allocatable :: Xel(:)
        real(dp), allocatable :: Yel(:)
        real(dp), allocatable :: Zel(:)
        real(dp), allocatable :: Xf(:)
        real(dp), allocatable :: Yf(:)
        real(dp), allocatable :: Zf(:)
        real(dp), allocatable :: DimsF(:,:)
        integer :: Exch_mat_nr                   !> Number of rows in the exchange coupling matrix
        integer :: Exch_mat_nc                   !> Number of columns in the exchange coupling matrix
        integer :: Exch_mat_ntot                 !> Number of elements in the exchange coupling matrix
        integer, allocatable  :: Exch_mat_r(:)   !> Row indices for the exchange coupling matrix
        integer, allocatable  :: Exch_mat_c(:)   !> Column indices for the exchange coupling matrix
        real(dp), allocatable :: Exch_mat_v(:)   !> Values for the exchange coupling matrix
     end type MicroMagGridInfo

    !Additional information about the macrogeometry and sample shape
    ! used for computing the demagnetisation field
    !Note that current version approximates the overall shape of sample and macrogeometry by rectangular prisms
    type MicroMagMacrogeometry
        !> Number of copies on both sides of the intial domain which together form the macrogeometry
        integer,allocatable :: n_macro(:)
        real(dp),allocatable :: shiftVec(:)    !> How far to shift neighbouring domain copies along x, y, z
        real(dp),allocatable :: macroShape(:)  !> Sidelengths of macrogeometry prism
        real(dp),allocatable :: sampleShape(:) !> Sidelengths of sample prism
        ! Settings for exchange coupling between copies of the simulated domain
        logical,allocatable :: exchPBC(:)      !> Periodic boundary conditions along x, y and z for the exchange coupling
    end type MicroMagMacrogeometry

    !>-----------------
    !> Overall data structure for a micro magnetism problem.
    !> The design intention is such that a problem may be restarted given the information stored in this struct
    !>-----------------
    type MicroMagProblem
        !Below is stuff that needs to be provided by the "user":
        type(MicroMagGrid) :: grid                        !> Grid of the problem
        type(MicroMagMacrogeometry) :: macrogrid          !> Info on macrogeometry grid and sample shape
        
        real(DP),dimension(:,:),allocatable :: u_ea      !> Easy axis vectors that should have the dimensions (n,3) where n is the no. of grid points and thus u_ea(i,3) is the i'th point's z-component
        
        integer :: ProblemMode                            !> Defines the problem mode (new or continued from previous solution)
        
        integer :: solver                                 !> Determines what type of solver to use
        
        real(DP) :: exch_weight                           !> Sets the weight for the exchange operator calculation
        integer  :: exch_method                           !> Determines what type of exchange operator method to use
        integer  :: exch_interpn                          !> Determines what type of exchange interpolation method to use
    
        real(DP) :: gamma,alpha0,MaxT0                    !> User defined coefficients determining part of the problem.
        real(DP) :: tol,thres_value                       !> User defined coefficients for the ODE solver
        real(DP),dimension(:),allocatable :: Jfact,Kfact
        real(SP),dimension(:),allocatable :: Mfact
        real(DP),dimension(:),allocatable :: temperature  !> User defined system temperature
        real(DP),dimension(:),allocatable :: Tfact        !> Prefactor of the termal magnetic field
        logical :: includeThermal                         !> Whether thermal noise is included in the simulations
        
        real(DP),dimension(:,:),allocatable :: Hext       !> Applied field as a function of time. Size (nt,3) with the latter dimension specifying the spatial dimensions.
        real(DP),dimension(:,:),allocatable :: alpha      !> A time dependent damping parameter, i.e. as a function of time. Size (nt,1).
        
        real(DP),dimension(:),allocatable :: t              !> Time array for the desired output times
        real(DP),dimension(:),allocatable :: m0             !> Initial value of the magnetization
        real(DP),dimension(:),allocatable :: Ms,K0,K1,K2,A0 !> User defined coefficients determining part of the problem.
        real(DP),dimension(:,:,:),allocatable :: Kfact_arr  !> n,6,3 array specifying the local anisotropy constants. See updateAnisotropy for details
        
        real(DP),dimension(:),allocatable :: t_conv        !> Time array with the time values where the solution will be checked for convergence compared to the last timestep
        real(DP) :: conv_tol                               !> Converge criteria on difference between magnetization at different timesteps
        
        real(SP) :: demag_threshold                        !> Used for specifying whether the demag tensors should be converted to sparse matrices by defining values below this value to be zero
        real(SP) :: CV                                   !> The coefficient of variation (CV), i.e. the ratio of the standard deviation to the mean, which can be used to add an error to the demag field
        integer :: demag_ignore_steps                    !> Only compute the demag tensor every demag_ignore_steps'th-step in a calculation using the hysteresis-model. Otherwise the parameter is ignore (i.e. in the dynamic solver)
        
        integer :: setTimeDisplay                               !> Determines how often the timestep is shown in Matlab
        integer :: useCuda                                      !> Defines whether to attempt using CUDA or not
        integer :: useCVODE                                     !> Defines whether to attempt using CVODE or not
        integer :: usePrecision                                 !> Defines whether to use single (false) or double precision (true)
        integer :: useReturnHall                                !> Defines whether to return all the specific H-fields (exchange, demag) �(true) or not (false)
        integer :: passExch                                     !> Defines whether the exchange matrix is passed from Matlab/Python (true) or calculated localled (false).
        integer :: demag_approximation                          !> Flag for how to approximate the demagnetization tensor as specified in the parameters below
        integer :: demagTensorReturnState                       !> Flag describing how or if the demag tensor should be returned
        integer :: demagTensorLoadState                         !> Flag describing how or if to load the demag tensor (from disk e.g.)
        integer :: nThreadsMatlab                               !> Number of threads to use in the OpenMP demag tensor allocation
        integer,dimension(3) :: N_ave                           !> Number of points to average the demag tensor in in the recieving tile, N_ave(1) = N_x etc
        character*256 :: demagTensorFileOut, demagTensorFileIn  !> Filename (including path) for output (input) of demag tensor if it is to be returned as a file (demagTensorReturnState >2 and the value is equal to the length of the file including path)
        
        
        !Below is stuff that is computed when the solver initializes
        
        type(sparse_matrix_t) :: A_exch         !> Exchange term matrix
        
        type(MagTenseSparse),dimension(6) :: K_s           !> Sparse matrices (used if the threshold is >0 )
        type(MagTenseSparse_c),dimension(6) :: K_s_c       !> Sparse matrices (used if the threshold is >0 ), complex version
        
        real(SP),dimension(:,:),allocatable :: Kxx,Kxy,Kxz  !> Demag field tensor split out into the nine symmetric components
        real(SP),dimension(:,:),allocatable :: Kyy,Kyz      !> Demag field tensor split out into the nine symmetric components
        real(SP),dimension(:,:),allocatable :: Kzz          !> Demag field tensor split out into the nine symmetric components

        real(SP),dimension(:),allocatable :: Kxx_shape,Kxy_shape,Kxz_shape  !> Shape correction tensor components
        real(SP),dimension(:),allocatable :: Kyy_shape,Kyz_shape            !> Shape correction tensor components
        real(SP),dimension(:),allocatable :: Kzz_shape                      !> Shape correction tensor components

        integer,dimension(:,:),allocatable :: tensorMap     !> A map of the unique entries in the demagnetization tensor
        logical,dimension(:,:),allocatable :: tensorMapX, tensorMapY, tensorMapZ   !> The sign of the different components in the demagnetization tensor map
        
        real(DP),dimension(:),allocatable :: Axx,Axy,Axz,Ayy,Ayz,Azz    !> Anisotropy vectors assuming local anisotropy only, i.e. no interaction between grains
        real(DP),dimension(:,:,:),allocatable :: CrystalAxis, K0_arr !> The local crystal coordinates and the local anisotropy constants. See updateAnisotropy for details        
        real(DP),dimension(:,:,:),allocatable :: A0_map     !> The specified exchange constant, mapped out on a grid, for computing the finite difference exchange operator for uniform grids        
        
        type(DFTI_DESCRIPTOR), POINTER :: desc_hndl_FFT_M_H       !> Handle for the FFT MKL stuff





        !--------------- demag tensor neighbour stuff ---------------
        integer :: fmm_cells_per_node
        real(DP) :: fmm_eps 
        integer :: ifunif = 1 
        integer :: nlmin = 0
        integer :: nlmax = 5
        logical :: use_fmm = .true.
        integer :: allow_fmm_short_circuit = 0
        integer :: fmm_min_n = 20000    !> Minimum number of cells to use FMM

        integer, dimension(:,:), pointer :: nbr_idx(:,:) => null()    !> Neighbour indices for each target cell
        integer, dimension(:), pointer :: n_nbors(:) => null()
        real(SP),dimension(:,:,:,:), pointer:: Nnbr(:,:,:,:) => null()
        real(SP),dimension(:,:,:,:), pointer:: diffTens(:,:,:,:) => null()

        type(MagTenseSparse), dimension(6) :: K_fmm_s
        logical :: K_fmm_built = .false.
        type(matrix_descr) :: K_fmm_descr_s

        real(SP),dimension(:,:), contiguous, pointer :: Kxx_fmm,Kxy_fmm,Kxz_fmm  !> Demag field tensor split out into the nine symmetric components
        real(SP),dimension(:,:), contiguous, pointer :: Kyy_fmm,Kyz_fmm      !> Demag field tensor split out into the nine symmetric components
        real(SP),dimension(:,:), contiguous, pointer :: Kzz_fmm          !> Demag field tensor split out into the nine symmetric components
        integer :: dummy_run = 0    !> Flag to indicate whether the demag tensor neighbour test setup has been done

        
    end type MicroMagProblem
    
    !>-----------------
    !> Data structure for a micro magnetism solution.
    !> The design intention is such that a problem may be restarted given the information stored in this struct
    !>-----------------
    type MicroMagSolution
        type(MicroMagGridInfo) :: gridinfo                                 !> GridInfo of the problem
        
        real(DP),dimension(:),allocatable :: HjX,HjY,HjZ                   !> Effective fields for the exchange term (X,Y and Z-directions, respectively)
        real(DP),dimension(:),allocatable :: HhX,HhY,HhZ                   !> Effective fields for the external field (X,Y and Z-directions, respectively)
        real(DP),dimension(:),allocatable :: HkX,HkY,HkZ                   !> Effective fields for the anisotropy energy term (X,Y and Z-directions, respectively)        
        real(SP),dimension(:),allocatable :: HmX,HmY,HmZ                   !> Effective fields for the demag energy term (X,Y and Z-directions, respectively)        
        real(DP),dimension(:),allocatable :: Mx,My,Mz                      !> The magnetization components used internally as the solution progresses
        real(SP),dimension(:),allocatable :: Mx_s,My_s,Mz_s                !> The magnetization components used internally as the solution progresses in single precision
        complex(kind=4),dimension(:),allocatable :: Mx_FT, My_FT, Mz_FT    !> Fourier transform of Mx, My and Mz (complex)
        complex(kind=4),dimension(:),allocatable :: HmX_c,HmY_c,HmZ_c      !> Complex version of the demag field, used for the Fourier cut-off approach
        
        real(DP),dimension(:),allocatable :: t_out          !> Output times at which the solution was computed
        real(DP),dimension(:,:,:,:),allocatable :: M_out    !> The magnetization at each of these times (nt,ntot,nt_Hext,3)
        real(DP),dimension(:,:,:,:),allocatable :: H_exc    !> The exchange field at each of these times (nt,ntot,nt_Hext,3)
        real(DP),dimension(:,:,:,:),allocatable :: H_ext    !> The external field at each of these times (nt,ntot,nt_Hext,3)
        real(DP),dimension(:,:,:,:),allocatable :: H_dem    !> The demagnetization field at each of these times (nt,ntot,nt_Hext,3)
        real(DP),dimension(:,:,:,:),allocatable :: H_ani    !> The anisotropy field at each of these times (nt,ntot,nt_Hext,3)
        
        real(DP),dimension(:,:),allocatable :: pts          !> n,3 array with the points (x,y,z) of the centers of the tiles
        
        real(SP),dimension(:),allocatable :: u1,u2,u3,u4,u5,u6  !> Random vectors to add noise to the demagnetization field
        
        integer :: HextInd                              !> Index specifying which external field in the input array we have reached in the explicit method

#if USE_FMM3D
        class(FMM3DTree), pointer :: fmm_tree => null()    !> FMM tree for computing the demag field using FMM
#endif
    end type MicroMagSolution
    
    
    !>------------
    !> Parameters
    !>------------
    
    integer,parameter :: gridTypeUniform=1,gridTypeTetrahedron=2,gridTypeUnstructuredPrisms=3
    integer,parameter :: ProblemModeNew=1,ProblemModeContinued=2
    integer,parameter :: MicroMagSolverExplicit=1,MicroMagSolverDynamic=2,MicroMagSolverImplicit=3
    integer,parameter :: MicroMagExchMethodDirectLaplacianNeumann=1,MicroMagExchMethodGGNeumann=2
    integer,parameter :: MicroMagExchInterpnExtended=1,MicroMagExchInterpnCompact=2
    integer,parameter :: useCudaTrue=1,useCudaFalse=0
    integer,parameter :: DemagApproximationNothing=1,DemagApproximationThreshold=2,DemagApproximationFFTThreshold=3,DemagApproximationThresholdFraction=4,DemagApproximationFFTThresholdFraction=5
    integer,parameter :: DemagTensorReturnNot=1,DemagTensorReturnMemory=2
    !!@todo Do NOT have useCVODETrue/-False variables both here and in IntegrationDataTypes.
    integer,parameter :: useCVODETrue=1,useCVODEFalse=0
    integer,parameter :: passExchTrue=1,passExchFalse=0
    integer,parameter :: usePrecisionTrue=1,usePrecisionFalse=0
    integer,parameter :: useReturnHallTrue=1,useReturnHallFalse=0
    integer,parameter :: useFMMTrue=1,useFMMFalse=0
    
end module MicroMagParameters    