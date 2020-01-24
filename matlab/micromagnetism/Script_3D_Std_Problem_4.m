function Script_3D_Std_Problem_4(fig1,MySim)
disp('Running Matlab model')
mu0 = 4*pi*1e-7;

if (isfield(MySim,'ShowTheResult'))
    if ~(MySim.ShowTheResult == 0)
        if (~exist('fig1','var'))
            figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
        end
        if (~isempty(fig1))
            figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
        end
    end
end
if (~exist('MySim','var'))   
    HystDir = 1/mu0*[1,1,1] ;
    MySim.nGrid = [36,9,1]';
    MySim.Field_dir = HystDir;
end

%% Calculates dynamic solution to mumag std prob 4
clearvars -except fig1 MySim

%% addpaths
addpath('../MEX_files');
addpath('../examples/util');


%% Create Preset
% a = memory
if exist('MySim','var')
    if ~isfield(MySim,'SaveTheResult')
        MySim.SaveTheResult = 0;
    end
end
MySim.K0 = 0 ;
MySim.Kz = 0 ;
MySim.A0 = 1.3e-11 ;
MySim.Ms = 8.0e5 ;
% MySim.MaxH = 0.1 ;
% MySim.MaxHx = MySim.MaxH*sqrt(3) ;  
% MySim.nHalfTimes = round([20;5;0]*1);
% MySim.nGrid = round([36;9;1]);
% MySim.nGrid = round([32;8;1]);

% MySim.DemagTensorFileName = 'DemagStdProb4_x109y26z0' ;
MySim.DemagTensorFileName = nan;
MySim.HeffLimMagn = 0.14; % The magic number that kinda sorta works

MySim.Lx = 500e-9 ;
MySim.Ly = 125e-9 ;
MySim.Lz = 3e-9 ; 

%% Get initial state
MySim.MaxT = 1e-6;
% MaxT0 = 2 ;
% MySim.alpha = @(t) -65104e-17*(10.^(7*min(t,MaxT0)/MaxT0)); 
MySim.alpha = 4.42e3;
MySim.tt = 0 ;

MySim.InitialState = 'OldSol' ;
OldSolution = load('Stdprob4Initx51y13z1'); 
MySim.Nold = OldSolution.N;
MySim.SigmaIN = OldSolution.SigmaInitial;
MySim.HystDir = [1,1,1]./sqrt(3) ;

MySim.HsX = @(t) MySim.HystDir(1).*t ;
MySim.HsY = @(t) MySim.HystDir(2).*t ;
MySim.HsZ = @(t) MySim.HystDir(3).*t ;

MySim.UseImplicitSolver = 0;
MySim.UseExplicitSolver = 1;
MySim.UseDynamicSolver = 0;

% FFT 
%     MySim.FFTdims = [1,0,0] ;
%     MySim.thresholdFract = 0.9 ; 
%     MySim.use_sparse = 1 ;

MySim.SimulationName = 'PhysParams_TestDynamicsStdProb4_init' ;

[SigmaInit,~,~,InteractionMatrices] = ComputeTheSolution(MySim);

for k=1:size(SigmaInit,1) 
    Sigma = SigmaInit(k,:).' ;
    NN = round(numel(Sigma)/3) ;

    SigmaX = Sigma(0*NN+[1:NN]) ;
    SigmaY = Sigma(1*NN+[1:NN]) ;
    SigmaZ = Sigma(2*NN+[1:NN]) ;
    SigmaN = sqrt(SigmaX.^2+SigmaY.^2+SigmaZ.^2) ;
    Mx(k) = mean(SigmaX./SigmaN) ;
    My(k) = mean(SigmaY./SigmaN) ;
    Mz(k) = mean(SigmaZ./SigmaN) ;
    Mk(k) = Mx(k)*MySim.HystDir(1) + My(k)*MySim.HystDir(2) + Mz(k)*MySim.HystDir(3) ;
end
% disp(num2str([Mx(end);My(end);Mz(end)],10)) ;
% figure ; plot(Mx) ; hold on ; plot(My) ; plot(Mz) ;
if (MySim.ShowTheResult)
    figure; quiver(InteractionMatrices.X(:),InteractionMatrices.Y(:),SigmaX,SigmaY); axis equal;  title('Starting state - Matlab')
end

%% Perform dynamic calculation with first field
MySim.alpha = 4.42e3 ;
MySim.gamma = 2.21e5 ;

%--- Field 1
%    MySim.Field_dir = 1/mu0*[-24.6,4.3,0]/1000 ;
%--- Field 2
%     MySim.Field_dir = 1/mu0*[-35.5,-6.3,0]/1000 ;

MySim.HsX = @(t) MySim.Field_dir(1) + 0.*t ;
MySim.HsY = @(t) MySim.Field_dir(2) + 0.*t ;
MySim.HsZ = @(t) MySim.Field_dir(3) + 0.*t ;
    
MySim.t = linspace(0,1e-9,200);    % 1 ns for now

MySim.InitialState = 'OldSol' ;
MySim.SigmaIN = SigmaInit(end,:) ;
MySim.Nold = MySim.nGrid ;

MySim.UseImplicitSolver = 0;
MySim.UseExplicitSolver = 0;
MySim.UseDynamicSolver = 1;

MySim.HeffLimMagn = -21; % The magic number that kinda sorta works

MySim.SimulationName = 'PhysParams_TestDynamicsStdProb4_1st_field' ;

tic
SigmaSol1 = ComputeTheSolution(MySim);
toc

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
    Mk(k) = Mx(k)*MySim.HystDir(1) + My(k)*MySim.HystDir(2) + Mz(k)*MySim.HystDir(3) ;
end

if (MySim.ShowTheResult)
    plot(fig1, MySim.t,Mx,'ro')
    plot(fig1, MySim.t,My,'go')
    plot(fig1, MySim.t,Mz,'bo')
end

end