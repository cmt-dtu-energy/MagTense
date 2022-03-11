import numpy as np


def get_micromag_solver(nm=None):
    solver_dict = {
        None: -1,
        'explicit': 1,
        'dynamic': 2,
        'implicit': 3,
    }
    return solver_dict[nm]


def get_micromag_problem_mode(nm=None):
    mode_dict = {
        None: -1,
        'new': 1,
        'old': 2,
    }
    return mode_dict[nm]


def get_micromag_gridtype(nm=None):
    grid_dict = {
        None: -1,
        'uniform': 1,
        'tetrahedron': 2,
        'unstructuredPrisms': 3
    }
    return grid_dict[nm]


def get_micromag_demag_approx(mode):
    '''
    None: No approximation applied to the demag tensor
    Threshold: cut-off all values below threshold and make the matrix sparse
    FFT threshold: apply the threshold in fourier-space through:
        NP = FT * N * IFT
        NP(NP<threshold) = 0
        NP = sparse(NP)
        H = IFT ( NP * FT * M )
        with FT = fft( eye(n,n) ) and IFT = ifft( eye(n,n) )
    Threshold fraction: Cut-off all values below a certain fraction specified
        by threshold and make the matrix sparse
    FFT threshold fraction: Cut-off all values below a certain fraction specified
        by threshold and make the matrix sparse
    '''
    approx_dict = {
        None: 1,
        'threshold': 2,
        'fft_thres': 3,
        'threshold_fraction': 4,
        'fft_threshold_fraction': 5
    }
    return approx_dict[mode]


class DefaultMicroMagProblem():
    '''
    Defines a class with the default settings for running a MicroMag
    problem using the Fortran implementation in MagTense.
    Note that the naming convention is dictated in the subroutine 
    getProblemFieldNames in the module MagTenseMicroMagIO defined in
    the file MagTenseMicroMagIO.f90 in the sub-project MagTenseMicroMag
    to the main project MagTense.
    '''
    def __init__(
        self,
        res=[3,3,1],
        hext_fct=lambda t: np.atleast_2d(t).T * [0, 0, 0],
        alpha_fct=lambda t: np.atleast_2d(t).T * 0
    ):
        self.grid_n = res
        self.ntot = np.prod(self.grid_n)
        # Size of the (rectangular) domain in each direction [m]
        self.grid_L = [500e-9, 125e-9, 3e-9]

        # These are set to zero for a non-tetrahedron grid
        self.grid_pts = np.zeros(shape=(self.ntot,3), dtype=np.float64, order='F')
        self.grid_ele = np.zeros(shape=(4, self.ntot), dtype=np.float64, order='F')
        self.grid_nnod = 0
        self.grid_nod = np.zeros(shape=(self.grid_nnod,3), dtype=np.float64, order='F')
        self.grid_abc = np.zeros(shape=(self.ntot,3), dtype=np.float64, order='F')
        
        # Defines the grid type which currently only supports "uniform"
        self.grid_type = get_micromag_gridtype('uniform')

        # easy axes of each cell
        self.u_ea = np.zeros(shape=(np.prod(res), 3), dtype=np.float64, order='F')

        # new or old problem
        self.prob_mode = get_micromag_problem_mode('new')

        # solver type ('explicit', 'implicit' or 'dynamic')
        self.solver = get_micromag_solver('dynamic')

        # Exchange term constant
        self.A0 = 1.3e-11
        # Demag magnetization constant [A/m]
        self.Ms = 8e5
        # Anisotropy constant
        self.K0 = 0 

        # Precession constant [m/A*s]
        self.gamma = 0
        self.alpha = 0.02

        # if set to zero then the alpha parameter remains constant.
        # if MaxT0 > 0 then alpha = alpha0 * 10^( 7 * min(t,MaxT0)/MaxT0 )
        # thus scaling with the solution time. This is used in the explicit
        # solver for tuning into the correct time scale of the problem
        self.MaxT0 = 2
        
        # Solution times
        self.nt = 1000
        self.t = np.linspace(0, 1, self.nt)
        
        self.nt_explicit = 1
        self.t_explicit  = 0
                
        self.nt_conv = 1
        self.t_conv = 0
        
        self.setTimeDis = 10
        
        # Initial magnetization (mx = m0(1:n), my = m(n+1:2n), mz = m(2n+1:3n)
        # with n = no. of elements )
        self.m0 = np.zeros(shape=(self.ntot, 3), dtype=np.float64, order='F')
        self.theta = np.pi * np.random.rand(self.ntot)
        self.phi = 2 * np.pi * np.random.rand(self.ntot)
        self.m0[:,0] = np.sin(self.theta) * np.cos(self.phi)
        self.m0[:,1] = np.sin(self.theta) * np.sin(self.phi)
        self.m0[:,2] = np.cos(self.theta)

        # Initial value of the demag threshold is zero, i.e. it is not used
        self.dem_thres = 0
        # Set use cuda to default not
        self.use_CUDA = 0
		# Set use CVODE to default
        self.use_CVODE = 0
        # Set the demag approximation to the default, i.e. use no approximation
        self.dem_appr = get_micromag_demag_approx(None)
        
        #  Set alpha as function of time
        self.set_alpha(alpha_fct, np.zeros(shape=(1)))
        
        self.save = 0
        self.show = 1

        self.dir_name = ''
        self.simulationName = ''
        self.fileName = ''
        self.recompute = 0
        self.ext_mesh = 0

        # Type of external mesh (Voronoi/Tetra)
        self.mesh_type = ''
        self.ext_mesh_filename = ''
        self.N_filename = 0

        # Set the default return N behavior
        self.set_return_N_filename('t')
        self.set_load_N_filename('t')

        # Relative tolerance for the Fortran ODE solver
        self.tol = 1e-4
        
        # Fortran ODE solver, when a solution component Y(L) is less in
        # magnitude than thres_value its set to zero
        self.thres = 1e-6
        
        # The convergence tolerence, which is the maximum change in 
        # magnetization between two timesteps
        self.conv_tol = 1e-4
    
        # The exchange operator matrix
        self.exch_nval = 0
        self.exch_nrow = 0
        self.exch_val = 0.
        self.exch_rows = 0
        self.exch_rowe = 0
        self.exch_col = 0
    
        # Defines if and how to return the N tensor (1 do not return, 2 return
        # in memory and >2 return as a file with file length = N_ret
        # N_ret {mustBeGreaterThan(N_ret,0),mustBeInteger(N_ret)}= 1
        
        # Defines whether the N tensor should be loaded rather than calculated
        # 2 = load from memory, >2 load from file with filename length N_load
        # N_load {mustBeGreaterThan(N_load,0),mustBeInteger(N_load)}= 1
        
        # Defines certain constant parameters for easy use and guarantee of being
        # consistent with the Fortran counter parts
        # Defines the constant for not returning N in any way (default)
        # returnNModeNot = 1            
        # Defines the constant for attempting to return the N tensor directly in memory
        # returnNModeMemory = 2            
        # Return as a file and the length of the filename is equal to mode (converted to an int32)
        # returnNModeNFile = 3

    def set_Hext(self, fct, t_Hext):
        '''
        Calculates the applied field as a function of time on the time grid
        already specified in this instance (self) of the class given the
        function handle fct
        '''
        self.nt_Hext = len(t_Hext)
        self.Hext = np.zeros(shape=(self.nt_Hext, 4))
        self.Hext[:,0] = t_Hext
        self.Hext[:,1:4] = fct(t_Hext)
    
    def set_ddHext(self, fct, t_ddHext):
        self.nt_ddHext = len(t_ddHext)        
        self.ddHext = np.zeros(shape=(self.nt_ddHext, 4))
        self.ddHext[:,0] = t_ddHext
        self.ddHext[:,1:4] = fct(t_ddHext)
        
    def set_Hext_time(self, nt):
        self.nt_Hext = nt
        self.t_Hext = np.linspace(self.t[0], self.t[-1], self.nt_Hext)
    
    def set_time(self, t):
        self.t = t
        self.nt = len(t)

    def set_alpha(self, fct, t_alpha):
        self.nt_alpha = len(t_alpha)
        self.alphat = np.zeros(shape=(self.nt_alpha, 2))
        self.alphat[:,0] = t_alpha
        self.alphat[:,1] = fct(t_alpha)
     
    def set_load_N_filename(self, filename):
        self.N_load = len(filename)
        self.N_file_in = filename
    
    def set_return_N_filename(self, filename):
        self.N_ret = len(filename)
        self.N_file_out = filename

    def set_use_CUDA(self, enabled):
       self.use_CUDA = 1 if enabled else 0

    def set_use_CVODE(self, enabled):
       self.use_CVODE = 1 if enabled else 0

    def set_save(self, enabled):
       self.save = 1 if enabled else 0

    def set_show(self, enabled):
       self.save = 1 if enabled else 0
    
    def set_convergence_check_time(self, t_conv):
        self.t_conv  = t_conv
        self.nt_conv = len(t_conv)
    
    def set_time_explicit(self, t_explicit):
        self.t_explicit  = t_explicit
        self.nt_explicit = len(t_explicit)
