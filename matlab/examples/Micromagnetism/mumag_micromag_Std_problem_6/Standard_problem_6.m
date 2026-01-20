function Switching_field = mumag_std_6_Rasmus(settings, x_steps, field_steps, cart_dir, options)
%STANDARD_PROBLEM_6 
%A function script to setup and simulate proposed standard problem 6
%The problem is described in "Proposal for a micromagnetic standard problem: domain wall pinning at phase boundaries" (2022), Heistracher et al. [https://arxiv.org/abs/2107.07855] 
%
%Syntax:
%------
%   Standard_problem_6()
%   [Switching_field] = Standard_problem_6( settings, xsteps, tsteps, cart_dir, options)
%
%Description of syntax:
%------
%   Standard_problem_6() 
%       Uses the default parameters to solve mumag problem 6 and displays the results on screen
%
%   Standard_problem_4( settings, xsteps, tsteps, cart_dir, options)
%       Takes 1 or 2 input argument which specifies the material parameters varied in the left side of the modelled domain, the time and space discretization, and the sample orientation. Additional options can also be specified
%
%   Switching_field = Standard_problem_6( settings, xsteps, tsteps, cart_dir, options)
%       As above but returns the switching field
%
%Input arguments:
%------
%   settings : A combination of 'a', 'j' and 'k'
%      settings detail which of the three material parameters change in the left side of the simulation domain. 
%      a: Exchange interaction strength
%      k: Anisotropy constant
%      j: Saturation magnetization
%      any combination may be used, such as 'ak', 'km', 'akj', '', etc.

%-------
%   The script setups up and runs the mumag standard problem 6
%
%Version: 1.1.0
%Author:  Rasmus Bjørk
%Date:    2026.01.19

arguments
    settings char                                      = 'akj';        %--- Parameter controlling the experiment
    x_steps (1,1) {mustBeNumeric}                      = 80;           %--- The spatial resolution
    field_steps (1,1) {mustBeNumeric}                  = 201;          %--- The field resolution
    cart_dir (1,1) string                              = 'x';          %--- The directions along which the geometry is oriented. Has no influence on the results, but can test the physics is correct in different directions
    options.use_CUDA {mustBeNumericOrLogical}          = true;         %--- Use CUDA for the calculations
    options.use_uniform_mesh {mustBeNumericOrLogical}  = true;         %--- Use a uniform or unstructured mesh
    options.use_CVODE {mustBeNumericOrLogical}         = true;         %--- Use CVODE for the numerical time evolution
    options.ShowTheResult {mustBeNumericOrLogical}     = true          %--- Show the result
end

mu0 = 4*pi*1e-7;

addpath('../../../../../MagTense/matlab/MEX_files');
addpath('../../../../../MagTense/matlab/util');

%% Known analytical solutions
HPs=containers.Map({'akj','ak','aj','a','kj','k','j',''},[1.568,1.089,1.206,0.838,1.005,0.565,0,0]); % Theoretical pinning fields for the different cases.
T_HP=HPs(settings); 
title_str=['"',settings,'" -- ', num2str(field_steps),' steps, theoretical pinning field: ',num2str(T_HP),' T'];

%% Physical Parameters
Ms = 1/mu0 ;    % 1 T
K0 = 1e6 ;      % [J/m3]
A0 = 1e-11 ;    % [J/m3]

Ms_soft = 0.25/mu0; % 0.25 T
A0_soft = 0.25e-11; % [J/m]
K0_soft = 1e5;      % [J/m3]

gamma0 = 2.2128e5;
alpha0 = 1 ;
gamma = gamma0/(1+alpha0^2);
alpha = alpha0*gamma0/(1+alpha0^2);

%% Mesh
thisGridL = [1e-9,1e-9,1e-9]; % Grid size
resolution = [1 1 1];

switch cart_dir
    case 'x'
        dim = 1;
    case 'y'
        dim = 2;
    case 'z'
        dim = 3;
end

if (options.use_uniform_mesh)
    resolution(dim) = x_steps;
else
    x_start = 1/2*80e-9/x_steps;
    x_end   = 80e-9-1/2*80e-9/x_steps;

    pos_out  = zeros(x_steps,3);
    dims_out = 1e-9*ones(x_steps,3);

    pos_out(:,dim)  = linspace(x_start,x_end,x_steps);
    dims_out(:,dim) = 80e-9/x_steps*ones(x_steps,1);

    resolution(1) = x_steps;
end

thisGridL(dim)  = 80e-9;


%% Problem structure creation
problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
problem = problem.setUseCuda( options.use_CUDA );
problem = problem.setMicroMagDemagApproximation('threshold_fraction');  % Turn off demag field
problem.dem_thres = 2;                                                  % Turn off demag field
if ~(options.use_uniform_mesh)
    problem = problem.setMicroMagGridType('unstructuredPrisms');
    problem.exch_weigh = 8;
    problem = problem.setMicroMagExchMethod( 'DirectLaplacianNeumann'  );
    problem = problem.setMicroMagExchInterpn( 'Extended'  ); 

    problem.grid_pts = pos_out ;
    problem.grid_abc = dims_out ;
end
problem = problem.setSolverType( 'UseDynamicSolver' );
problem = problem.setMicroMagSolver( 'Dynamic' );
problem.ReturnHall = int32(1);
problem.grid_L   = thisGridL ;

%% General problem parameters
problem.alpha = alpha;
problem.gamma = gamma;
problem.Ms = Ms*ones(prod(resolution),1);
problem.K0 = K0*ones(prod(resolution),1);
problem.A0 = A0*ones(prod(resolution),1);

% Set lower properties for left region
n_middle = x_steps/2;
if contains(settings,'a')
    problem.A0(1:n_middle) = A0_soft;
end
if contains(settings,'k')
    problem.K0(1:n_middle) = K0_soft;
end
if contains(settings,'j') 
    problem.Ms(1:n_middle) = Ms_soft;
end

%% Grain anisotropies
switch cart_dir
    case 'x'
        easyX = 1 ;
        easyY = 0 ;
        easyZ = 0 ;
    case 'y'
        easyX = 0 ;
        easyY = 1 ;
        easyZ = 0 ;
    case 'z'
        easyX = 0 ;
        easyY = 0 ;
        easyZ = 1 ;
end
problem.u_ea = repmat([easyX,easyY,easyZ],x_steps,1);

%% Applied Field
HystDir = normalize([easyX,easyY,easyZ],'norm') ;

% Time-dependent applied field
Hext0= 2e7*HystDir/mu0 .* 20e-9; % Initial applied field for ill-conditioned 'k' and 'kj' cases. Could be 0 in any other case.
HextFct = @(t) 2e7*HystDir/mu0 .* t' + Hext0;
problem = problem.setHext( HextFct, linspace(0,100e-9,field_steps) );

% Initial State
switch cart_dir
    case 'x'
        init_stat=[-1,0.3,0];
    case 'y'
        init_stat=[0,-1,0.3];
    case 'z'
        init_stat=[0,0.3,-1];
end
problem.m0(:,1) = init_stat(1)/norm(init_stat) ;
problem.m0(:,2) = init_stat(2)/norm(init_stat) ;
problem.m0(:,3) = init_stat(3)/norm(init_stat) ;

problem.m0(1:n_middle,dim) = -problem.m0(1:n_middle,dim);

% Time/field grid on which to solve the problems
problem = problem.setUseCVODE( options.use_CVODE );
problem = problem.setTime( linspace(0,100e-9,field_steps) );


%% Run simulation
solution = struct();
prob_struct = struct(problem);  %convert the class obj to a struct so it can be loaded into fortran
[solution, GridInfo] = MagTenseLandauLifshitzSolver_mex( prob_struct, solution );

[Mx,My,Mz,mx,my,mz] = computeMagneticMomentGeneralMesh(solution.M,GridInfo.Volumes) ;

switch cart_dir
case 'x'
    M=Mx;
case 'y'
    M=My;
case 'z'
    M=Mz;
end

%% Plot average magnetization
H = -mu0*solution.H_ext(:,1,1,dim);
Switching_field = min(H(M > 1-1e-3));

if (options.ShowTheResult)
    figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
    plot(fig1,H,M)
    xlabel('\mu_0H_{app} [T]')
    ylabel('\langle m \rangle')
    title(title_str)
end

end