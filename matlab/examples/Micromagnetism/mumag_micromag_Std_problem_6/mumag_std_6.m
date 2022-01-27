function mumag_std_6(settings,useCVODE,tsteps)

addpath('../../../../../MagTense/matlab/MEX_files');
addpath('../../../../../MagTense/matlab/util');
%addpath('../../micromagnetism');

fnameSave = 'mumag_std_6_prismal' ;
extra='';
useFinDif = ~true;
runFortranPart = true;
runMatlabPart = false;
HPs=containers.Map({'akj','ak','aj','a','kj','k','j',''},[1.568,1.089,1.206,0.838,1.005,0.565,0,0]);

if ~exist('settings','var')
    settings='akj';%'akj';
end
solver='CVODE';
if ~exist('useCVODE','var')
    useCVODE=true;
elseif ~useCVODE
    solver='RK';
end
if ~exist('tsteps','var')
    tsteps=401;
end
T_HP=HPs(settings);
title_str=['"',settings,'" -- ',solver, ', ', num2str(tsteps),' steps, theoretical pinning field: ',num2str(T_HP),' T'];
rng(9)
%% Physical Parameters
demag_approx = 'threshold_fraction';
dem_thres = 2;

mu0 = 4*pi*1e-7;

Ms = 1/mu0 ; % 1 T
K0 = 1e6 ;    % J/m3
A0 = 1e-11 ;  % J/m3

Ms_soft=0.25/mu0; %0.25 T
A0_soft=0.25e-11; %J/m
K0_soft=1e5; % J/m3

gamma0 = 2.2128e5;%0;
alpha0 = 1 ; %1e15 ; 
gamma = gamma0/(1+alpha0^2);
alpha = alpha0*gamma0/(1+alpha0^2);

%% Generate Voronoi regions with intergrain region

thisGridL = [80e-9,1e-9,1e-9];
X = [0,thisGridL(1)] ;
Y = [0,thisGridL(2)] ;
Z = [0,thisGridL(3)] ;

nGrains = 4 ; %5 ; % number of grains
offSetD = 0 ; %2.*7e-9 ; % [m]
x_gen = rand(nGrains,1).*thisGridL(1) ;
y_gen = rand(nGrains,1).*thisGridL(2) ;
z_gen = rand(nGrains,1).*thisGridL(3) ;

[x_gen,y_gen,z_gen,~] = DoLloydIteration(X,Y,Z,x_gen,y_gen,z_gen,70);%3) ;
% x_gen = thisGridL(1)*[1/4;3/4];
% y_gen = thisGridL(2)*[1/2;1/2];
% z_gen = thisGridL(3)*[1/2;1/2];
x_gen = thisGridL(1)*[1/8;3/8;5/8;7/8];
% y_gen = thisGridL(2)*1/2*ones(size(x_gen));
% z_gen = thisGridL(3)*1/2*ones(size(x_gen));

fname = 'temp\ThisVoronoiMeshMumag06a' ;
VoronoiStruct = GenerateVoronoiCells(X,Y,Z,x_gen,y_gen,z_gen,offSetD,fname);

%NNs = [80,7,1] ; %--- NNs is a 3-elements array corresponding to the desired grid resolution (before refinment).
% NNs = [80,2,2] ; %--- NNs is a 3-elements array corresponding to the desired grid resolution (before refinment).
NNs = [80,1,1] ; %--- NNs is a 3-elements array corresponding to the desired grid resolution (before refinment).

ress=3 ;
mesh_cart = refineRectGrid_v6( [X(1),X(end)],[Y(1),Y(end)],[Z(1),Z(end)],[],[],[],NNs,ress,VoronoiStruct);
% mesh_cart = refineRectGrid_v6( [X(1),X(end)],[Y(1),Y(end)],[Z(1),Z(end)],x_gen,y_gen,z_gen,NNs,ress,VoronoiStruct);

%--- Analyze and plot Cartesian Unstructured Mesh
mesh_cart = GetDomainsFromMeshHull(mesh_cart,VoronoiStruct) ;

GridInfo = CartesianUnstructuredMeshAnalysis(mesh_cart.pos_out,mesh_cart.dims_out) ;

resolution = [length(mesh_cart.pos_out) 1 1];

%% Problem structure creation

problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
problem = problem.setUseCuda( true );
problem.dem_appr = getMicroMagDemagApproximation(demag_approx);
problem.dem_thres = dem_thres;
problem.grid_type = getMicroMagGridType('unstructuredPrisms');
problem = problem.setSolverType( 'UseDynamicSolver' );
problem.solver = getMicroMagSolver( 'Dynamic' );

problem.grid_pts = mesh_cart.pos_out ;
problem.grid_abc = mesh_cart.dims_out ;
problem.grid_L   = thisGridL ;

%% Save the parameters

problem.alpha = alpha;
problem.gamma = gamma;
problem.Ms = Ms*ones(prod(resolution),1);
problem.K0 = K0*ones(prod(resolution),1);
problem.A0 = A0;%*ones(prod(resolution),1);

%% Grain anisotropies
easyX = 1 ;
easyY = 0 ;
easyZ = 0 ;

for i=1:nGrains
    problem.u_ea(mesh_cart.iIn{i},:)=repmat([easyX,easyY,easyZ],numel(mesh_cart.iIn{i}),1);
end
if contains(settings,'k')% soft anisotropy in left phase
    for i=1:nGrains/2
    problem.K0(mesh_cart.iIn{i},:)=K0_soft;
    end
end

%% Exchange matrix

InteractionMatrices.GridInfo = GridInfo;
InteractionMatrices.X = GridInfo.Xel ;
InteractionMatrices.Y = GridInfo.Yel ;
InteractionMatrices.Z = GridInfo.Zel ;

GridInfo.TheTs=GridInfo.TheDs;

Aexch=ones(numel(GridInfo.Xel),1);
if contains(settings,'a')
    for i=1:nGrains/2
        Aexch(mesh_cart.iIn{i})= A0_soft/A0 ;
        %problem.A0(mesh_cart.iIn{i},:)=A0_soft;
    end
end
if ~useFinDif
if numel(unique(GridInfo.Zel))-1
    [D2X,D2Y,D2Z] = ComputeDifferentialOperatorsFromMesh_Neumann04_GGDirLap(GridInfo,'extended',8,"DirectLaplacianNeumann",Aexch);
    InteractionMatrices.A2 = D2X + D2Y + D2Z ;
elseif numel(unique(GridInfo.Yel))-1
    [D2X,D2Y] = ComputeDifferentialOperatorsFromMesh_Neumann04_GGDirLap(GridInfo,'extended',8,"DirectLaplacianNeumann",Aexch);
    InteractionMatrices.A2 = D2X + D2Y ;
else
    [D2X] = ComputeDifferentialOperatorsFromMesh_Neumann04_GGDirLap(GridInfo,'extended',8,"DirectLaplacianNeumann",Aexch);
    InteractionMatrices.A2 = D2X ;
end

%--- Convert the exchange matrix to sparse
problem = PrepareExchangeMatrix(InteractionMatrices.A2,problem) ;

else
% Finite difference stencil
    D2A = ComputeDifferentialOperatorsFromMesh_Neumann_FiniteDifference(GridInfo,Aexch);
    problem = PrepareExchangeMatrix(D2A,problem) ;
end

%% Applied Field

HystDir = normalize([easyX,easyY,easyZ],'norm') ;

%time-dependent applied field
Hext0= 2e7*HystDir/mu0 .* 20e-9;
HextFct = @(t) 2e7*HystDir/mu0 .* t' + Hext0;
%problem = problem.setHext( HextFct, linspace(0,100e-9,2001) );
problem = problem.setHext( HextFct, linspace(0,100e-9,2) );

problem.HextFct{1} = @(t) 2e7*HystDir(1)/mu0 .* t' + Hext0;
problem.HextFct{2} = @(t) 2e7*HystDir(2)/mu0 .* t' + Hext0;
problem.HextFct{3} = @(t) 2e7*HystDir(3)/mu0 .* t' + Hext0;

% Initial State
init_stat=[-1,0.3,0];
problem.m0(:,1) = init_stat(1)/norm(init_stat) ;
problem.m0(:,2) = init_stat(2)/norm(init_stat) ;
problem.m0(:,3) = init_stat(3)/norm(init_stat) ;

for i=1:nGrains/2 
    if contains(settings,'m') || contains(settings,'j')% Setting lower Ms for left region -- will this work?
        problem.Ms(mesh_cart.iIn{i}) = Ms_soft;
    end
    % setting left half init state to point right.
    problem.m0(mesh_cart.iIn{i},1) = -problem.m0(mesh_cart.iIn{i},1);
end

% time grid on which to solve the problems
if useCVODE
    problem = problem.setUseCVODE( true );
end
problem = problem.setTime( linspace(0,100e-9,tsteps) );
problem.setTimeDis = int32(10);
% problem = problem.setConvergenceCheckTime( linspace(0,40e-9,2) );

% problem.conv_tol = 1e-6;
% problem.tol = 5e-7;
% problem.tol = 1e-5;
% fnameSave=['_inc_tol_1e5_',fnameSave];
% solver=['_inc_tol_1e5_',solver];

if useFinDif
    extra=[extra,'FinDifAexch'];
end
    extra=[extra,'_1DTest'];

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
    %% Hysteresis Loop

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
    % plot(problem.Hext(1:10:end,1),MxMean2)
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
    %problem.m0(:) = SigmaInit(end,:);
    %problem.m0 = single(problem.m0);

    problem = problem.setShowTheResult(true);
    problem = problem.setSaveTheResult(false);
    problem.DirectoryFilename = ['Matlab_simulations/Matlab_resolution_' num2str(resolution(1)) '_' num2str(resolution(2)) '_' num2str(resolution(3))];
    problem.SimulationName = [fnameSave,'MATLAB_ode45',extra,'SolSetting',settings,'_',num2str(tsteps),'_tsteps'];

    oldPath = cd('../../../../../MagTense/matlab/micromagnetism_matlab_only_implementation');
    
    problem.DemagTensorFileName= [oldPath,'/',problem.DirectoryFilename,problem.SimulationName,'_demag.mat'];
    prob_struct=struct(problem);
    %if ~exist(problem.DemagTensorFileName,'file')
        [ProblemSetupStruct, ~] = SetupProblem(prob_struct);
        ProblemSetupStruct.thresholdFract=dem_thres;
        ProblemSetupStruct.FFTdims=[];
        %ProblemSetupStruct.use_single=true;
        %ProblemSetupStruct.use_sparse=false;
        InteractionMatrices = CalculateInteractionMatrices(ProblemSetupStruct);
        InteractionMatrices.A2=D2A;
        save(ProblemSetupStruct.DemagTensorFileName,'InteractionMatrices') ;
    %end
    tic
    SigmaSol1 = ComputeTheSolution(prob_struct);
    toc
    cd(oldPath);

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
    fnameSave='mumag_std_6_prismalMATLAB_ode45FinDifAexch';
    solver='MATLAB_ODE45';
    Hp = mu0*HextFct(min(problem.t(Mx>1-1e-3)));
    Hp = Hp(:,1);
    save(['Std_prob_6\',fnameSave,'SolutionSetting',settings,'_',num2str(tsteps),'_tsteps_large_Hext.mat'],'problem','Hp','SigmaSol1') ;
    savefig(['Std_prob_6\',solver,'_',num2str(tsteps),'_',settings,'_no_precond'])
    if (ShowTheResult)
        plot(fig1, problem.t,Mx,'ro')
        plot(fig1, problem.t,My,'go')
        plot(fig1, problem.t,Mz,'bo')
    end
end
close all
%{
pause(inf)
%% Plot
hF = figure('color','w') ;
 ppsz = [2*20,19] ;
    ppps = [0,0,2*20,19] ;
    set(hF,'PaperUnits','centimeters','PaperSize',ppsz,'PaperPosition',ppps) ;

hA1 = subplot(1,2,1) ;


CartesianUnstructuredMeshPlot(mesh_cart.pos_out,mesh_cart.dims_out,GridInfo,mesh_cart.iIn,fname,hA1)  ;
set(hA1,'visible','off') ;

hA2 = subplot(2,2,2) ;



hL1 = plot(H,M,'.-k','linewidth',1.5,'markersize',12,'parent',hA2) ;
hold on ; hL2 = plot(Hc,0,'hw','markerfacecolor','r','markersize',9,'parent',hA2) ;
text(Hc,0,'   H_C','horizontalalignment','left','verticalalignment','bottom','color','r','parent',hA2,'fontsize',13)

hL3 = plot(-H_N,0,'hw','markerfacecolor','b','markersize',9,'parent',hA2) ;
text(-H_N,0,'H_N   ','horizontalalignment','right','verticalalignment','bottom','color','b','parent',hA2,'fontsize',13) 
set(gca,'ylim',[-1,+1].*1.2,'xlim',[min(Hvalues),max(Hvalues)]) ;
legend([hL2,hL3],{'Coercive Field','Nucleation Field'},'fontsize',13,'location','northwest')
grid on
xlabel('H [T]','fontsize',13)
ylabel('M / M_s [-]','fontsize',13)
set(hA2,'linewidth',1.5,'fontsize',13)
set(gcf,'units','normalized','position',[0.05,0.05,.9,.9]) ;


theTextBoxA = annotation(hF,'textbox','position',[.6-.02,.35+.05,.3,.1],'string',[num2str(nGrains),' Grains'],'LineStyle','none','fontsize',13,'FontWeight','bold') ; 

theTextBox1 = annotation(hF,'textbox','position',[.6,.3+.05,.3,.1],'string',[],'LineStyle','none','fontsize',13) ; 
set(theTextBox1,'string',[...
    'M_s^{ } = ',DoTheUnitThing(Ms*mu0,'T'),newline,...
    'A_0^{ } = ',DoTheUnitThing(A0,'J/m^3'),newline,...
    'K_0^{ } = ',DoTheUnitThing(K0,'J/m^3'),newline,...
    '\sigma_{} = ',num2str(sigmaEasy),newline,...
    ]) ;
if offSetD > 0 
theTextBoxB = annotation(hF,'textbox','position',[.75-.02,.35+.05,.3,.1],'string','Inter-Grain','LineStyle','none','fontsize',13,'FontWeight','bold') ; 
theTextBox2 = annotation(hF,'textbox','position',[.75,.3+.05,.3,.1],'string',[],'LineStyle','none','fontsize',13) ; 
set(theTextBox2,'string',[...
    'M_s^* = ',DoTheUnitThing(Msig*mu0,'T'),newline,...
    'A_0^* = ',DoTheUnitThing(A0*A0igFactor,'J/m^3'),newline,...
    'K_0^* = ','0 J/m^3',newline,...
    '\delta_{} = ',DoTheUnitThing(2*offSetD,'m'),newline,...
     ]) ;
end
 theTextBoxC = annotation(hF,'textbox','position',[.6-.02,.12+.05,.3,.1],'string',['Unstructured Cartesian Mesh (',num2str(ress),' refinments)'],'LineStyle','none','fontsize',13,'FontWeight','bold') ; 
 
 theTextBox3 = annotation(hF,'textbox','position',[.6,.07+.05,.3,.1],'string',[],'LineStyle','none','fontsize',13) ; 
set(theTextBox3,'string',[...
    'L_x = ',DoTheUnitThing(thisGridL(1),'m'),' ',newline,...
    'L_y = ',DoTheUnitThing(thisGridL(2),'m'),' ',newline,...
    'L_z = ',DoTheUnitThing(thisGridL(3),'m'),' ',newline,...
    ]) ;

 theTextBox4 = annotation(hF,'textbox','position',[.75,.07+.05,.3,.1],'string',[],'LineStyle','none','fontsize',13) ; 
set(theTextBox4,'string',[...
    'N_x base = ',num2str(NNs(1)),' ',newline,...
    'N_y base = ',num2str(NNs(2)),' ',newline,...
    'N_z base = ',num2str(NNs(3)),' ',newline,...
    ]) ;

    set(gcf,'PaperUnits','centimeters','PaperSize',ppsz,'PaperPosition',ppps) ;
    set(gcf,'PaperUnits','inches','Renderer','painters') ;
    figure(gcf) ;
    eval(['print -dpdf ',fnameSave,'Result','.pdf']) ;

%}
end