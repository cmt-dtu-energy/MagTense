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
useFinDif = ~true; % Switch between Direct Laplacian method and finite difference method
runFortranPart = true;
runMatlabPart = false;

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
T_HP=HPs(settings);
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
X = [0,thisGridL(1)] ;
Y = [0,thisGridL(2)] ;
Z = [0,thisGridL(3)] ;

% Base filename
fname = 'temp\ThisVoronoiMeshMumag06a' ;

tetname = 'TetraModel_80_1_1_4GrainsRegular';
if ~true
    % Load Comsol created mesh
    warning('Mesh is loaded, not generated')
    % model is a structure passed to TetrahedralMeshAnalysis
    load(tetname,'model') ; % Voronoi & Tetra
    model.Mesh.Nodes = meshTetra.Vert.' ;  % Positions of the vertexes
    model.Mesh.Elements = meshTetra.A.' ;  % Vertices-elements connectivity matrix (starts from 1)
elseif true
    % Load given mesh
    nodes=importdata('donaumesh\paul.knt')*1e-9;
    nodes=nodes+40e-9*repmat([1 0 0],size(nodes,1),1); % shift positions to start at x=0
    elements=importdata('donaumesh\paul.ijk');
    model = createpde();
    geometryFromMesh(model,nodes',elements(:,1:4)');
elseif true
    % Create mesh in Matlab
    mesh_res = thisGridL(1)/80;
    % Fancy stuff to ensure exact boundary in the middle
    model = CreateTetraMeshStdProb6(thisGridL,mesh_res);
else
    % Create mesh in Comsol
    Res = thisGridL(1)/160;
    L = 1;
    oldFolder = cd('../../Library');
    %% Generate Tetra Mesh
    meshTetra = GetMeshFromComsol01(X,Y,Z,fname,Res,L);
    meshTetra = GetDomainsFromMeshHull(meshTetra,VoronoiStruct) ;
    cd(oldFolder)
    model.Mesh.Nodes = meshTetra.Vert.' ;  % Positions of the vertexes
    model.Mesh.Elements = meshTetra.A.' ;  % Vertices-elements connectivity matrix (starts from 1)
    save(tetname,'model')
end

% Plot the mesh
figure ; pdeplot3D(model,'FaceAlpha',0.1) ;

%--- Setup the problem
resolution = [size(model.Mesh.Elements,2),1,1];
disp(['Tetra N_grid = ' num2str(prod(resolution))])

%% Problem structure creation
problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
problem.grid_type = getMicroMagGridType('tetrahedron');

problem = problem.setUseCuda( true );
problem.dem_appr = getMicroMagDemagApproximation(demag_approx); % turn off demag field
problem.dem_thres = dem_thres;                                  % turn off demag field
problem = problem.setSolverType( 'UseDynamicSolver' );
problem.solver = getMicroMagSolver( 'Dynamic' );

%--- Information on the grid
GridInfo = TetrahedralMeshAnalysis(model) ;
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
if runMatlabPart
    %% Matlab
    problem = problem.setSolverType( 'UseDynamicSolver' );

    problem = problem.setShowTheResult(true);
    problem = problem.setSaveTheResult(false);
    problem.DirectoryFilename = ['Matlab_simulations/Matlab_resolution_' num2str(resolution(1)) '_' num2str(resolution(2)) '_' num2str(resolution(3))];
    problem.SimulationName = [fnameSave,'MATLAB_ode45',extra,'SolSetting',settings,'_',num2str(tsteps),'_tsteps'];

    oldPath = cd('../../../../../MagTense/matlab/micromagnetism_matlab_only_implementation');
    
    problem.DemagTensorFileName= [oldPath,'/',problem.DirectoryFilename,problem.SimulationName,'_demag.mat'];
    prob_struct=struct(problem);
    % Pre-calculate interaction matrices
    %if ~exist(problem.DemagTensorFileName,'file')
        [ProblemSetupStruct, ~] = SetupProblem(prob_struct);
        ProblemSetupStruct.thresholdFract=dem_thres; % Turn off demagnetization
        ProblemSetupStruct.FFTdims=[];
        %ProblemSetupStruct.use_single=true;
        %ProblemSetupStruct.use_sparse=false;
        InteractionMatrices = CalculateInteractionMatrices(ProblemSetupStruct);
        InteractionMatrices.A2=D2X + D2Y + D2Z ; % Supply our own exchange matrix
        save(ProblemSetupStruct.DemagTensorFileName,'InteractionMatrices') ;
    %end
    tic
    SigmaSol1 = ComputeTheSolution(prob_struct);
    toc
    cd(oldPath);

    % Get hysteresis loop
    for k=1:size(SigmaSol1,1) 
        Sigma = SigmaSol1(k,:).' ;
        NN = round(numel(Sigma)/3) ;

        SigmaX = Sigma(0*NN+[1:NN]) ;
        SigmaY = Sigma(1*NN+[1:NN]) ;
        SigmaZ = Sigma(2*NN+[1:NN]) ;
        SigmaN = sqrt(SigmaX.^2+SigmaY.^2+SigmaZ.^2) ;
        Mx(k) = mean(SigmaX./SigmaN) ;
        My(k) = mean(SigmaY./SigmaN) ;
        Mz(k) = mean(SigmaZ./SigmaN) ;
    end
    figure
    plot(mu0*HextFct(problem.t),Mx)
    xlabel('\mu_0H_{app} [T]')
    ylabel('\langle m_x \rangle')
    title_str=['"',settings,'" -- MATLAB Ode45, ', num2str(tsteps),' steps, theoretical pinning field: ',num2str(T_HP),' T'];
	title(title_str)
    fnameSave='mumag_std_6_tetraMATLAB_ode45';
    solver='MATLAB_ODE45';
    Hp = mu0*HextFct(min(problem.t(Mx>1-1e-3)));
    Hp = Hp(:,1);
    save(['Std_prob_6\',fnameSave,'SolutionSetting',settings,'_',num2str(tsteps),'_tsteps.mat'],'problem','Hp','SigmaSol1') ;
    savefig(['Std_prob_6\',solver,'_',num2str(tsteps),'_',settings,'_no_precond'])
    if (ShowTheResult)
        plot(fig1, problem.t,Mx,'ro')
        plot(fig1, problem.t,My,'go')
        plot(fig1, problem.t,Mz,'bo')
    end
end
close all
end