function mumag_std_6_tetra(settings,useCVODE,tsteps)
%%% Calculates the "standard" micromagnetic problem number 6 on a 
%%% tetrahedral mesh. The problem is described in
%%% "Proposal for a micromagnetic standard problem: domain wall pinning at
%%% phase boundaries" (2022), Heistracher et al. [https://arxiv.org/abs/2107.07855] 
%%% 
%%% Variables:
%%% 
%%% settings detail which of the three parameters change.
%%%     a: Exchange interaction strength
%%%     k: Anisotropy constant
%%%     j (or m): Saturation magnetization
%%% any combination may be used, such as 'ak', 'km', 'akj', '', etc.
%%% NOTE that cases 'kj' and 'k' are ill posed.
%%% 
%%% useCVODE is a switch between the solvers. CVODE (true) and RK (false)
%%% 
%%% tsteps is the number of time steps used in the problem. CVODE is almost
%%% immune, but RK is heavily dependent on the time resolution asked for.

addpath('../../../../../MagTense/matlab/MEX_files');
addpath('../../../../../MagTense/matlab/util');
addpath('../../../../../MagTense/matlab/micromagnetism_matlab_only_implementation');
%addpath('../../micromagnetism');

fnameSave = 'mumag_std_6_tetra' ;
extra=''; % Extra text fragments for graph titles and filenames
runFortranPart = true;

% Theoretical domain wall pinning fields
HPs=containers.Map({'akj','ak','aj','a','kj','k','j',''},[1.568,1.089,1.206,0.838,1.005,0.565,0,0]);

%% Defaults
if ~exist('settings','var')
    settings='akj';
end
solver='CVODE';
if ~exist('useCVODE','var')
    useCVODE=true;
elseif ~useCVODE
    solver='RK';
end
if ~exist('tsteps','var')
    tsteps=401; % Overkill, but gives RK a fighting chance.
end
T_HP=HPs(settings); % Theoretical pinning field in our case.
title_str=['"',settings,'" -- ',solver, ', ', num2str(tsteps),' steps, theoretical pinning field: ',num2str(T_HP),' T'];
rng(9)
%% Physical Parameters
demag_approx = 'threshold_fraction';
dem_thres = 2;

mu0 = 4*pi*1e-7;

Ms = 1/mu0 ; % 1 T
K0 = 1e6 ;    % [J/m3]
A0 = 1e-11 ;  % [J/m3]

Ms_soft=0.25/mu0; % 0.25 T
A0_soft=0.25e-11; % [J/m]
K0_soft=1e5; % [J/m3]

gamma0 = 2.2128e5;
alpha0 = 1 ;
gamma = gamma0/(1+alpha0^2);
alpha = alpha0*gamma0/(1+alpha0^2);

%% Mesh generation. Generate Voronoi regions with intergrain region
thisGridL = [80e-9,1e-9,1e-9];
tetname = 'TetraModel_80_1_1_4GrainsRegular'; % basename for mesh file
load(tetname,'model','GridInfo');

%--- Setup the problem
resolution = [size(model.Mesh.Elements,2),1,1];
disp(['Tetra N_grid = ' num2str(prod(resolution))])

%% Problem structure creation
problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
problem = problem.setMicroMagGridType('tetrahedron');

problem = problem.setUseCuda( true );
problem = problem.setMicroMagDemagApproximation(demag_approx);  % turn off demag field
problem.dem_thres = dem_thres;                                  % turn off demag field
problem = problem.setSolverType( 'UseDynamicSolver' );
problem = problem.setMicroMagSolver( 'Dynamic' );

%--- Information on the grid
problem.grid_pts    = [GridInfo.Xel, GridInfo.Yel, GridInfo.Zel] ;
problem.grid_ele    = int32(model.Mesh.Elements(:,:)) ;
problem.grid_nod    = model.Mesh.Nodes(:,:) ;      
problem.grid_nnod   = int32(length(problem.grid_nod));
problem.grid_L   = thisGridL ;

%% General problem parameters
problem.alpha = alpha;
problem.gamma = gamma;
problem.Ms = Ms*ones(prod(resolution),1);
problem.K0 = K0*ones(prod(resolution),1);
problem.A0 = A0;

% region division
left_region=GridInfo.Xel<thisGridL(1)/2;
%% Grain anisotropies
easyX = 1 ;
easyY = 0 ;
easyZ = 0 ;

% set anisotropy direction
problem.u_ea=repmat([easyX,easyY,easyZ],numel(GridInfo.Xel),1);
if contains(settings,'k')% soft anisotropy in left phase
    problem.K0(left_region,:)=K0_soft;
end

%% Exchange matrix
InteractionMatrices.GridInfo = GridInfo;
InteractionMatrices.X = GridInfo.Xel ;
InteractionMatrices.Y = GridInfo.Yel ;
InteractionMatrices.Z = GridInfo.Zel ;

Aexch=ones(numel(GridInfo.Xel),1);
if contains(settings,'a') % soft exchange stiffness in left phase
    Aexch(left_region)= A0_soft/A0 ; % Hack so exchange strength varies
end

if numel(unique(GridInfo.Zel))-1
    [D2X,D2Y,D2Z] = computeDifferentialOperatorsFromMesh_DirectLap(GridInfo,'extended',8,"DirectLaplacianNeumann",Aexch);
    InteractionMatrices.A2 = D2X + D2Y + D2Z ;
elseif numel(unique(GridInfo.Yel))-1
    [D2X,D2Y] = computeDifferentialOperatorsFromMesh_DirectLap(GridInfo,'extended',8,"DirectLaplacianNeumann",Aexch);
    InteractionMatrices.A2 = D2X + D2Y ;
else
    [D2X] = computeDifferentialOperatorsFromMesh_DirectLap(GridInfo,'extended',8,"DirectLaplacianNeumann",Aexch);
    InteractionMatrices.A2 = D2X ;
end

%--- Convert the exchange matrix to sparse
problem = problem.setExchangeMatrixSparse( InteractionMatrices.A2 );

%% Applied Field

HystDir = normalize([easyX,easyY,easyZ],'norm') ;

%time-dependent applied field
Hext0=0*2e6*1e-7*HystDir/mu0;    % Initial applied field for ill-conditioned 'k' and 'kj' cases.
                                 % Could be 0 in any other case.
HextFct = @(t) 2e7*HystDir/mu0 .* t'+Hext0;
problem = problem.setHext( HextFct, linspace(0,100e-9,2) );

problem.HextFct{1} = @(t) 2e7*HystDir(1)/mu0 .* t';
problem.HextFct{2} = @(t) 2e7*HystDir(2)/mu0 .* t';
problem.HextFct{3} = @(t) 2e7*HystDir(3)/mu0 .* t';

% Initial State
init_stat=[-1,0.3,0];
problem.m0(:,1) = init_stat(1)/norm(init_stat) ;
problem.m0(:,2) = init_stat(2)/norm(init_stat) ;
problem.m0(:,3) = init_stat(3)/norm(init_stat) ;

if contains(settings,'m') || contains(settings,'j')% Setting lower Ms for left region
    problem.Ms(left_region) = Ms_soft;
end
% setting left half init state to point right.
problem.m0(left_region,1) = -problem.m0(left_region,1);

% time grid on which to solve the problems
if useCVODE
    problem = problem.setUseCVODE( true );
end
problem = problem.setTime( linspace(0,100e-9,tsteps) );
problem.setTimeDis = int32(10); % Progress update intervals. Every 10th time step.

% Extra parameters to be played with if problem doesn't converge
% problem = problem.setConvergenceCheckTime( linspace(0,40e-9,2) );
% problem.conv_tol = 1e-6;
% problem.tol = 5e-7;
problem.tol = 5e-4;
extra=[extra,'reltol_5e-4'];
% fnameSave=['_inc_tol_1e5_',fnameSave];
% solver=['_inc_tol_1e5_',solver];

extra=[extra,'_Tetra'];

%% RUN !!
if runFortranPart
    solution = struct();
    w = warning ('off','MATLAB:structOnObject');
    prob_struct = struct(problem);  %convert the class obj to a struct so it can be loaded into fortran
    warning(w)
    tic
    solution = MagTenseLandauLifshitzSolver_mex( prob_struct, solution );
    RunTime = toc ;

    save([fnameSave,solver,extra,'Solution.mat'],'solution')
    
    %% Get hysteresis loop
    MxMean2 = zeros(tsteps,1) ;
    MyMean2 = zeros(tsteps,1) ;
    MzMean2 = zeros(tsteps,1) ;

    for ij = 1:tsteps
        Mx_arr = (solution.M(ij,:,1,1)) ;
        My_arr = (solution.M(ij,:,1,2)) ;
        Mz_arr = (solution.M(ij,:,1,3)) ;

        MxMean2(ij) = sum(GridInfo.Volumes(:).*Mx_arr(:))/sum(GridInfo.Volumes(:)) ;
        MyMean2(ij) = sum(GridInfo.Volumes(:).*My_arr(:))/sum(GridInfo.Volumes(:)) ;
        MzMean2(ij) = sum(GridInfo.Volumes(:).*Mz_arr(:))/sum(GridInfo.Volumes(:)) ;
    end
    H=-mu0*solution.H_ext(:,1,1,1);
    M=MxMean2;
    figure
    plot(H,M)
    xlabel('\mu_0H_{app} [T]')
    ylabel('\langle m_x \rangle')
    title(title_str)
    Hp = min(H(MxMean2>1-1e-3));
    save(['Std_prob_6\',fnameSave,solver,extra,'SolutionSetting',settings,'_',num2str(tsteps),'_tsteps.mat'],'solution','Hp','H','M') ;
    savefig(['Std_prob_6\',solver,extra,'_',num2str(tsteps),'_',settings,'_no_precond'])
end

close all
end