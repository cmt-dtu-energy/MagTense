%% Calculates dynamic solution to mumag std prob 4

clearvars

addpath('..\matlab\MEX_files')

%% Create Preset
a = memory
MySim.K0 = 0 ;
MySim.Kz = 0 ;
MySim.A0 = 1.3e-11 ;
MySim.Ms = 8.0e5 ;
MySim.MaxH = 0.1 ;
MySim.MaxHx = MySim.MaxH*sqrt(3) ;  
%nHalfTimes = [1,2;1,2;0,1] ;
%nHalfTimes = [0,19;0,4;0,1] ;
% nHalfTimes = [6;3;0] ;
% nHalfTimes = [19;4;1] ; %OLD GRID
% MySim.nHalfTimes = [51;13;1] ;
% MySim.nHalfTimes = [40;9;1] ;
% MySim.nHalfTimes = [40;9;0] ;
% MySim.nHalfTimes = [30;6;1] ;
% MySim.nHalfTimes = [19;4;1] ;
% MySim.nHalfTimes = [24;5;1] ;
% MySim.nHalfTimes = [51;13;3] ;
% MySim.nHalfTimes = [86;21;0] ;
% MySim.nHalfTimes = [109;26;0] ;
MySim.nHalfTimes = round([40;10;0]*1);

% nHalfTimes = [6;2;0] ;
% MySim.DemagTensorFileName = 'DemagStdProb4_x109y26z0' ;
MySim.DemagTensorFileName = nan;
MySim.HeffLimMagn = 0.14; % The magic number that kinda sorta works

MySim.Lx = 500e-9 ;
MySim.Ly = 125e-9 ;
MySim.Lz = 3e-9 ; 
% nHalfTimes = [12] ;
% MySim.FreqH = pi/2 ;
%% Get initial state

MySim.MaxT = 1 ;
MaxT0 = 2 ;
MySim.alpha = @(t) -65104e-17*(10.^(7*min(t,MaxT0)/MaxT0)); 
% MySim.tt = [linspace(0,1),linspace(1,0)]; %[linspace(0,1),
MySim.tt = 0 ;
% MySim.tt = -[linspace(1,0,26)]; %[linspace(0,1),
% MySim.SigmaXfun = @(X,Y,Z) 1/sqrt(3)+0.*X ;
% MySim.SigmaYfun = @(X,Y,Z) 1/sqrt(3)+0.*Y ;
% MySim.SigmaZfun = @(X,Y,Z) 1/sqrt(3)+0.*Z ;
% MySim.SigmaXfun = @(X,Y,Z) 1+0.*X ;
% MySim.SigmaYfun = @(X,Y,Z) 0.25+0.*Y ;
% MySim.SigmaZfun = @(X,Y,Z) 0.1+0.*Z ;

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
MySim.GifFilename = 'ThisPhysTest.gif' ;

SigmaInit = ComputeTheSolution(MySim);

%% Perform dynamic calculation with first field
MySim.alpha = -4.42e-6 ;
MySim.gamma = -2.21e-4 ;

MySim.HystDir = -[-24.6,4.3,0]./1000 ;

MySim.HsX = @(t) MySim.HystDir(1) + 0.*t ;
MySim.HsY = @(t) MySim.HystDir(2) + 0.*t ;
MySim.HsZ = @(t) MySim.HystDir(3) + 0.*t ;
    
MySim.t = linspace(0,1,200);    % 1 ns for now, ~3 ns for proper solution

MySim.InitialState = 'OldSol' ;
MySim.SigmaIN = SigmaInit(end,:) ;
MySim.Nold = 2*MySim.nHalfTimes + 1 ;

% MySim.SigmaInitialFileName = 'Stdprob4Initx51y13z1' ;

MySim.UseImplicitSolver = 0;
MySim.UseExplicitSolver = 0;
MySim.UseDynamicSolver = 1;

MySim.HeffLimMagn = -21; % The magic number that kinda sorta works

MySim.SimulationName = 'PhysParams_TestDynamicsStdProb4_1st_field' ;

SigmaSol1 = ComputeTheSolution(MySim);

for k=1:size(SigmaSol1,1) 
Sigma = SigmaSol1(k,:).' ;
NN = round(numel(Sigma)/3) ;
% N = round(NN^(1/3)) ;
SigmaX = Sigma(0*NN+[1:NN]) ;
SigmaY = Sigma(1*NN+[1:NN]) ;
SigmaZ = Sigma(2*NN+[1:NN]) ;
SigmaN = sqrt(SigmaX.^2+SigmaY.^2+SigmaZ.^2) ;
Mx(k) = mean(SigmaX./SigmaN) ;
My(k) = mean(SigmaY./SigmaN) ;
Mz(k) = mean(SigmaZ./SigmaN) ;
Mk(k) = Mx(k)*MySim.HystDir(1) + My(k)*MySim.HystDir(2) + Mz(k)*MySim.HystDir(3) ;
end
% 
figure
plot(MySim.t,My)
legend('Dynamic method, H_x = -24.6 mT, H_y = 4.3 mT')
ylabel('<M_y>/M_s')
xlabel('Time [ns]')
'' ;

% %% Perform dynamic calculation with second field
% % MySim.HystDir = -[-35.5,-6.3,0]./1000 ;
% MySim.HystDir = -36*[cosd(190),sind(190),0]./1000 ;
% % HystDir = [1,0,0] ;
% 
% MySim.HsX = @(t) MySim.HystDir(1) + 0.*t ;
% MySim.HsY = @(t) MySim.HystDir(2) + 0.*t ;
% MySim.HsZ = @(t) MySim.HystDir(3) + 0.*t ;
% % ThisPreSetFileNameCOMP = 'PhysParams_TestDirectCompare05.mat' ;
% MySim.SimulationName = 'PhysParams_TestDynamicsStdProb4_2nd_field' ;
% 
% SigmaSol2 = ComputeTheSolution(MySim);
% 
% for k=1:size(SigmaSol2,1) 
% Sigma = SigmaSol2(k,:).' ;
% NN = round(numel(Sigma)/3) ;
% % N = round(NN^(1/3)) ;
% SigmaX = Sigma(0*NN+[1:NN]) ;
% SigmaY = Sigma(1*NN+[1:NN]) ;
% SigmaZ = Sigma(2*NN+[1:NN]) ;
% Mx2(k) = mean(SigmaX) ;
% My2(k) = mean(SigmaY) ;
% Mz2(k) = mean(SigmaZ) ;
% Mk2(k) = Mx2(k)*MySim.HystDir(1) + My2(k)*MySim.HystDir(2) + Mz2(k)*MySim.HystDir(3) ;
% end
% %% Plot the solutions. Minus signs added to correspond to regular hysteresis plots.
% figure
% plot(MySim.t,My2,'or')
% legend('Dynamic method, H_x = -35.5 mT, H_y = -6.3 mT')
% ylabel('<M_y>/M_s')
% xlabel('Time [ns]')