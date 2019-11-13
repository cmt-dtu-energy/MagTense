function Script_3D_TestDynamicsStdProbl4(fig1)
disp('Running Matlab model')

if (~exist('fig1','var'))
    figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
end

%% Calculates dynamic solution to mumag std prob 4
clearvars -except fig1

%% Create Preset
a = memory
MySim.K0 = 0 ;
MySim.Kz = 0 ;
MySim.A0 = 1.3e-11 ;
MySim.Ms = 8.0e5 ;
MySim.MaxH = 0.1 ;
MySim.MaxHx = MySim.MaxH*sqrt(3) ;  
% MySim.nHalfTimes = round([20;5;0]*1);
MySim.nGrid = round([36;9;1]);

% MySim.DemagTensorFileName = 'DemagStdProb4_x109y26z0' ;
MySim.DemagTensorFileName = nan;
MySim.HeffLimMagn = 0.14; % The magic number that kinda sorta works

MySim.Lx = 500e-9 ;
MySim.Ly = 125e-9 ;
MySim.Lz = 3e-9 ; 

%% Get initial state
MySim.MaxT = 1 ; MySim.MaxT = 1e3 ; %% HERE !!!!!
MaxT0 = 2 ;
MySim.alpha = @(t) -65104e-17*(10.^(7*min(t,MaxT0)/MaxT0)); 
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
figure; quiver(InteractionMatrices.X(:),InteractionMatrices.Y(:),SigmaX,SigmaY); axis equal;  title('Starting state - Matlab')

%% Perform dynamic calculation with first field
MySim.alpha = -4.42e-6 ;
MySim.gamma = -2.21e-4 ;

MySim.Field_dir = -[-24.6,4.3,0]./1000 ;

MySim.HsX = @(t) MySim.Field_dir(1) + 0.*t ;
MySim.HsY = @(t) MySim.Field_dir(2) + 0.*t ;
MySim.HsZ = @(t) MySim.Field_dir(3) + 0.*t ;
    
MySim.t = linspace(0,1,200);    % 1 ns for now, ~3 ns for proper solution

MySim.InitialState = 'OldSol' ;
MySim.SigmaIN = SigmaInit(end,:) ;
MySim.Nold = MySim.nGrid ;

MySim.UseImplicitSolver = 0;
MySim.UseExplicitSolver = 0;
MySim.UseDynamicSolver = 1;

MySim.HeffLimMagn = -21; % The magic number that kinda sorta works

MySim.SimulationName = 'PhysParams_TestDynamicsStdProb4_1st_field' ;

SigmaSol1 = ComputeTheSolution(MySim);

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
 
plot(fig1, MySim.t,Mx,'ro')
plot(fig1, MySim.t,My,'go')
plot(fig1, MySim.t,Mz,'bo')
ylabel('<M_i>/M_s')
xlabel('Time [ns]')

%-- Compare with published solutions
data = load('FIELD_1_SM_DT25.txt');
% data2 = load('MicroMagus_TimeDep_AllMagnProj_Field2.txt');
plot(fig1,data(:,1),data(:,2),'r-');
plot(fig1,data(:,1),data(:,3),'g-');
plot(fig1,data(:,1),data(:,4),'b-');

xlim(fig1,[0 1])

end