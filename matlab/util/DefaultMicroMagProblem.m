classdef DefaultMicroMagProblem
    
    %% Defines a class with the default settings for running a MicroMag
    %% problem using the Fortran implementation in MagTense.
    %% Note that the naming convention is dictated in the subroutine 
    %% getProblemFieldNames in the module MagTenseMicroMagIO defined in
    %% the file MagTenseMicroMagIO.f90 in the sub-project MagTenseMicroMag
    %% to the main project MagTense.

properties
    %grid resolution
    grid_n
    %grid_n = [int32(3), int32(3), int32(1)];
    ntot
    %size of the (rectangular) domain in each direction
    grid_L

    %the pts for a tetrahedron grid, i.e. the center of the elements
    grid_pts
    %the elements for a tetrahedron grid
    grid_ele
    %the nodes for a tetrahedron grid
    grid_nod
    %the number of nodes in the tetrahedron grid
    grid_nnod
    
    %the dimensions (abc) on the unstructured prisms
    grid_abc
    
    %the exchange operator matrix
    exch_mat
    
    exch_nval
    exch_nrow
    exch_val
    exch_rows
    exch_rowe
    exch_col
    
    %defines the grid type which currently only supports "uniform"
    grid_type

    %easy axes of each cell
    u_ea
    %new or old problem
    ProblemMod
    %solver type ('Explicit', 'Implicit' or 'Dynamic')
    solver

    %Exchange term constant
    A0
    %demag magnetization constant
    Ms
    %Anisotropy constant
    K0

    %precession constant
    gamma

    %
    alpha

    MaxT0

    %Sets how often timestep is displayed from Fortran
    setTimeDis
    
    %initial magnetization (mx = m0(1:n), my = m(n+1:2n), mz = m(2n+1:3n) with
    %n = no. of elements )
    m0
    %x-direction
    theta
    phi    
    
    %threshold for whether to attempt making the demag tensors sparse
    % if dem_thres>0 then all values in the demag tensors that are
    % absolutely below dem_thres are set to exactly 0 and subsequently the
    % tensors are made into sparse matrices, i.e. if abs(K) < dem_thres
    % then K = 0
    dem_thres
    
    %defines which approximation (if any) to use for the demag tensor 
    dem_appr
    
    %defines if and how to return the N tensor (1 do not return, 2 return
    %in memory and >2 return as a file with file length = N_ret
    N_ret {mustBeGreaterThan(N_ret,0),mustBeInteger(N_ret)}=int32(1);
    
    %defines whether the N tensor should be loaded rather than calculated
    %2 = load from memory, >2 load from file with filename length N_load
    N_load {mustBeGreaterThan(N_load,0),mustBeInteger(N_load)}=int32(1);
    
    %filename to which N is written to
    N_file_out char = '';
    
    %filename from which N is loaded from
    N_file_in char = '';
    
    %filename for the Matlab save file
    DirectoryFilename = '';
    
    %simulation for the Matlab save file
    SimulationName = '';
    
    %filename for the Matlab save file
    FileName = '';
    
    %relative tolerance for the Fortran ODE solver
    tol = 1e-4;
    
    %Fortran ODE solver, when a solution component Y(L) is less in
    %magnitude than thres_value its set to zero
    thres = 1e-6;
    
    %the convergence tolerence, which is the maximum change in 
    %magnetization between two timesteps
    conv_tol = 1e-4;
    
    %defines whether to recompute the Interaction Matrices or not
    RecomputeInteractionMatrices = 0 ;
    
    %defines whether to use an External Mesh or not
    ExternalMesh = 0 ; 
    
    %type of external mesh (Voronoi/Tetra)
    MeshType = '' ;
    
    % filename of external mesh
    ExternalMeshFileName = '' ;
    
    % filename of Interaction matrices (demag tensor)
    DemagTensorFileName = 0;
   
    % function handle for external field
    HextFct = [] ;

    %The number of threads used by OpenMP for building the demag tensor
    nThreads = int32(1);
end

properties (SetAccess=private,GetAccess=public)
    %applied field, should be (nt_Hext,3)
    Hext
    %no. of time steps in the applied field    
    nt_Hext
    
    %derivative of the applied field (nt_Hext,3)
    ddHext
    %no. of time steps in the derivative of the applied field    
    nt_ddHext
    
    %solution times
    nt
    %should have size (nt,1)
    t
    
    %check for convergence times
    nt_conv
    %should have size (nt_conv,1)
    t_conv
    
    %alpha as function of time, should be (nt_alpha,2)
    alphat
    %size of the alpha array
    nt_alpha
    
    %solution times for the explicit solver
    t_explicit;
    %should have size (nt,1)
    nt_explicit 
    
    %defines whether to attempt using CUDA (will crash if no appropriate
    %NVIDIA driver is present or if insufficient memory is available. 0 for
    %do not use cuda, 1 for do use (int32)
    useCuda
    MagTenseLandauLifshitzSolver_mex
    
    %defines whether to attempt using CVODE (it is currently uncertain if
    %the code even compiles in this version without having CVODE installed,
    %regardless of whether or not CVODE is used). 0 for do not use CVODE, 
    %1 for do use (int32). Currently implemented alternative is RK_SUITE,
    %active if CVODE is not used.
    useCVODE

    %defines what precision is used for the demag tensor. Right now only
    %single is supported. All other varibales are double.
    usePres
    
    %defines whether to save the result or not
    SaveTheResult
    
    %defines whether to show the result or not
    ShowTheResult
    
    %defines which solver to use in Matlab
    SolverType
    
    
end

%defines certain constant parameters for easy use and guarantee of being
%consistent with the Fortran counter parts
properties (Constant)
    %defines the constant for not returning N in any way (default)
    returnNModeNot = int32(1);            
    %defines the constant for attempting to return the N tensor directly in
    %memory
    returnNModeMemory = int32(2);            
    %return as a file and the length of the filename is equal to mode (converted to an int32)
    returnNModeNFile = int32(3);            
    
end

methods
    function obj = DefaultMicroMagProblem(nx,ny,nz,HextFct)
        nx = int32(nx);
        ny = int32(ny);
        nz = int32(nz);
        %grid resolution
        obj.grid_n = [nx,ny,nz];
        %grid_n = [int32(3), int32(3), int32(1)];
        obj.ntot = prod(obj.grid_n);
        %size of the (rectangular) domain in each direction
        obj.grid_L = [500e-9,125e-9,3e-9];%m

        %these are set to zero for a non-tetrahedron grid
        obj.grid_pts = 0;
        obj.grid_ele = int32(0);
        obj.grid_nod = 0;
        obj.grid_nnod = int32(0);
        obj.grid_abc = 0;
        
        
        %defines the grid type which currently only supports "uniform"
        obj.grid_type = getMicroMagGridType('uniform');


        %easy axes of each cell
        obj.u_ea = zeros( obj.ntot, 3 );
        %new or old problem
        obj.ProblemMod = getMicroMagProblemMode( 'new' );
        %solver type ('Explicit', 'Implicit' or 'Dynamic')
        obj.solver = getMicroMagSolver( 'Dynamic' );

        % Exchange term constant
        obj.A0 = 1.3e-11;
        % demag magnetization constant
        obj.Ms = 8e5*ones(obj.ntot,1); %A/m
        %Anisotropy constant
        obj.K0 = zeros(obj.ntot,1); 

        %precession constant
        obj.gamma = 0; %m/A*s

        %
        obj.alpha = 0.02;

        %if set to zero then the alpha parameter remains constant.
        %if MaxT0 > 0 then alpha = alpha0 * 10^( 7 * min(t,MaxT0)/MaxT0 )
        %thus scaling with the solution time. This is used in the explicit
        %solver for tuning into the correct time scale of the problem
        obj.MaxT0 = 2;
        
        %solution times
        obj.nt = int32(1000);
        obj.t = linspace(0,1,obj.nt);
        
        obj.nt_explicit = int32(1);
        obj.t_explicit  = 0;
                
        obj.nt_conv = int32(1);
        obj.t_conv = 0;
        
        obj.setTimeDis = int32(10);

        if ~exist('HextFct')
            HextFct = @(t) (t>=0)' * [0,0,0];
        end
        obj.nt_Hext = 10;
        obj = obj.setHext(HextFct, obj.nt_Hext);
        
        %initial magnetization (mx = m0(1:n), my = m(n+1:2n), mz = m(2n+1:3n) with
        %n = no. of elements )
        obj.m0 = zeros(obj.ntot, 3 );
        %x-direction
        obj.theta = pi * rand(obj.ntot,1);
        obj.phi = 2*pi*rand(obj.ntot,1);
        obj.m0(1:obj.ntot) = sin(obj.theta) .* cos(obj.phi);
        obj.m0(obj.ntot+1:2*obj.ntot) = sin(obj.theta) .* sin(obj.phi);
        obj.m0(2*obj.ntot+1:3*obj.ntot) = cos(obj.theta);
        %initial value of the demag threshold is zero, i.e. it is not used
        obj.dem_thres = 0;
        %set use cuda to default not
        obj.useCuda = int32(0);
		%set use CVODE to default
        obj.useCVODE = int32(0);
        %set use CVODE to default
        obj.usePres = int32(0);
        %set the demag approximation to the default, i.e. use no
        %approximation
        obj.dem_appr = getMicroMagDemagApproximation('none');
        

        %--- Set alpha as function of time
        if ~exist('AlphaFct')
            AlphaFct = @(t) (t>=0)' * 0;
        end
        obj = obj.setAlpha(AlphaFct, 0);
        
        
        %set the default return N behavior
        obj = obj.setReturnNFilename('t');
        obj = obj.setLoadNFilename('t');
        
        
        obj.SaveTheResult = int32(0);
        obj.ShowTheResult = int32(1);

        obj.DirectoryFilename = '';
        obj.RecomputeInteractionMatrices = int32(0) ;
        obj.ExternalMesh = int32(0) ;
        obj.MeshType = '' ;
        obj.ExternalMeshFileName = '' ;
        obj.DemagTensorFileName = 0 ;
    end
    
    %%Calculates the applied field as a function of time on the time grid
    %%already specified in this instance (obj) of the class given the
    %%function handle fct
    function obj = setHext( obj, fct, t_Hext )
        obj.nt_Hext = int32(length(t_Hext));
        
        obj.Hext = zeros( obj.nt_Hext, 4 );
        obj.Hext(:,1) = t_Hext;
        obj.Hext(:,2:4) = fct( t_Hext );
    end
    
    function obj = setddHext( obj, fct, t_ddHext )
        obj.nt_ddHext = int32(length(t_ddHext));
        
        obj.ddHext = zeros( obj.nt_Hext, 4 );
        obj.ddHext(:,1) = t_ddHext;
        obj.ddHext(:,2:4) = fct( t_ddHext );
    end
        
    function obj = setHextTime( obj, nt )
        obj.nt_Hext = int32( nt );
        obj.t_Hext  = linspace( obj.t(1), obj.t(end), obj.nt_Hext );
    end
    
    function obj = setTime( obj, t )
        obj.t  = t;
        obj.nt = int32(length(t));
    end
    
    function obj = setConvergenceCheckTime( obj, t_conv )
        obj.t_conv  = t_conv;
        obj.nt_conv = int32(length(t_conv));
    end
    
    function obj = setTimeExplicit( obj, t_explicit )
        obj.t_explicit  = t_explicit;
        obj.nt_explicit = int32( length(t_explicit) );
    end
    
    function obj = setAlpha( obj, fct, t_alpha )
        obj.nt_alpha = int32(length(t_alpha));
       
        obj.alphat = zeros( obj.nt_alpha, 2 );
        obj.alphat(:,1) = t_alpha;
        obj.alphat(:,2) = fct( t_alpha );
    end
    
    function obj = setUseCuda( obj, enabled )
       if enabled
           obj.useCuda = int32(1);
           obj.MagTenseLandauLifshitzSolver_mex = @MagTenseLandauLifshitzSolver_mex;
       else
           obj.useCuda = int32(0);
%            obj.MagTenseLandauLifshitzSolver_mex = @MagTenseLandauLifshitzSolverNoCUDA_mex;
           obj.MagTenseLandauLifshitzSolver_mex = @MagTenseLandauLifshitzSolver_mex;
       end
    end
    
    function obj = setUseCVODE( obj, enabled )
       if enabled
           obj.useCVODE = int32(1);
       else
           obj.useCVODE = int32(0);
       end
    end
    
    function obj = setLoadNFilename( obj, filename )
        obj.N_load = int32(length(filename));
        obj.N_file_in = filename;
    end
    
    function obj = setReturnNFilename( obj, filename )
        obj.N_ret = int32(length(filename));
        obj.N_file_out = filename;
    end
    
    function obj = setSaveTheResult( obj, enabled )
       if enabled
           obj.SaveTheResult = int32(1);
       else
           obj.SaveTheResult = int32(0);
       end
    end
    
    function obj = setShowTheResult( obj, enabled )
       if enabled
           obj.ShowTheResult = int32(1);
       else
           obj.ShowTheResult = int32(0);
       end
    end
    
    function obj = setSolverType( obj, type_var )
        switch type_var
            case 'UseImplicitSolver'
                obj.SolverType = 1;
            case 'UseImplicitStepsSolver'
                obj.SolverType = 2;
            case 'UseExplicitSolver'
                obj.SolverType = 3;
            case 'UseDynamicSolver'
                obj.SolverType = 4;
        end
    end

    function obj = setExchangeMatrixSparse( obj, ExchangeMatrix )
    % Convert the Exchange matrix to CSR and store it in the problem statement
        [v,c,rs,re]   = convertToCSR(ExchangeMatrix);
        obj.exch_nval = int32(numel(v));
        obj.exch_nrow = int32(numel(rs));
        obj.exch_val  = double(v);
        obj.exch_rows = int32(rs);
        obj.exch_rowe = int32(re);
        obj.exch_col  = int32(c);

        disp(['The demag tensor will require around ' num2str(((3*numel(rs)*(3*numel(rs) + 1)/2))*4/(10^9)) ' Gb'])
    end
    
    %Override struct function for a final check before handing to Fortran
    function obj2 = struct(obj)
        if length(obj.Ms)==1 % Check if Ms is vectorized
            obj.Ms=obj.Ms*ones(obj.ntot,1);
        end
        if length(obj.K0)==1 % Check if K0 is vectorized
            obj.K0=obj.K0*ones(obj.ntot,1);
        end
        mnorm=vecnorm(obj.m0,2,2); % Check if input array is normalized
        normcondfail=abs(mnorm-ones(obj.ntot,1)) >= obj.tol;
        if any(normcondfail)
            if all(mnorm(normcondfail)==0*mnorm(normcondfail))
                warning('Zero magnetization in initial array')
            else
                warning('Initial array not normalized -- Normalizing')
                obj.m0=obj.m0./mnorm;
            end
        end
        warning('off','MATLAB:structOnObject')
        obj2=builtin('struct',obj); % Actual struct conversion
        warning('on','MATLAB:structOnObject')
    end
    
end
end