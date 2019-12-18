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

    %defines the grid type which currently only supports "uniform"
    grid_type


    %easy axes of each cell
    u_ea
    %new or old problem
    ProblemMod
    %solver type ('Explicit', 'Implicit' or 'Dynamic')
    solver

    % Exchange term constant
    A0
    % demag magnetization constant
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
    
    %defines whether to attempt using CUDA (will crash if no appropriate
    %NVIDIA driver is present or if insufficient memory is available. 0 for
    %do not use cuda, 1 for do use (int32)
    useCuda
    
    %defines which approximation (if any) to use for the demag tensor 
    dem_appr
end

properties (SetAccess=private,GetAccess=public)
    %applied field, should be (nt_Hext,3)
    Hext
    %no. of time steps in the applied field    
    nt_Hext
    %time axis of the applied field
    t_Hext
    
    %solution times
    nt
    %should have size (nt,1)
    t
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
        obj.Ms = 8e5; %A/m
        %Anisotropy constant
        obj.K0 = 0; 

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
        
        obj.setTimeDis = int32(10);

        if ~exist('HextFct')
            HextFct = @(t) (t>=0)' * [0,0,0];
        end
        obj = obj.setHextTime( 10 * obj.nt );
        obj = obj.setHext(HextFct);
        
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
        %set the demag approximation to the default, i.e. use no
        %approximation
        obj.dem_appr = getMicroMagDemagApproximation('none');
    end
    
    %%Calculates the applied field as a function of time on the time grid
    %%already specified in this instance (obj) of the class given the
    %%function handle fct
    function obj = setHext( obj, fct )
        
        obj = obj.setHextTime( 10 * obj.nt );
        obj.Hext = zeros( obj.nt_Hext, 4 );
        
        obj.Hext(:,1) = obj.t_Hext;
        
        obj.Hext(:,2:4) = fct( obj.t_Hext );
    end
    
    function obj = setTime( obj, t )
        obj.t = t;
        obj.nt = int32(length(t));
    end
    
    function obj = setHextTime( obj, nt )
        obj.nt_Hext = int32( nt );
        obj.t_Hext = linspace( min(obj.t), max(obj.t), obj.nt_Hext );
    end
    
    function obj = setUseCuda( obj, enabled )
        
       if enabled
           obj.useCuda = int32(1);
       else
           obj.useCuda = int32(0);
       end
    end
    
end
end