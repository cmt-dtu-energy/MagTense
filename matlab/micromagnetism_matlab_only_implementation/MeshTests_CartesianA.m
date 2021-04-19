function MeshTests_Cartesian
disp('Running Matlab model')

%% Calculates dynamic solution to mumag std prob 4

mu0 = 4*pi*1e-7;

%% addpaths
addpath('../MEX_files');
addpath('../util');

%% Create Preset

resolution = [1*9,1*7,5]; % NOT STANDARD
MySim_ini = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
%     MySim_ini.A0 = 0
MySim_ini.grid_L = [125e-9,125e-9/2,125e-9/4];% NOT STANDARD
MySim_ini.DemagTensorFileName = 'CompareInteractionMatrices.mat' ;
MySim_ini.RecomputeInteractionMatrices = 1 ;
MySim_ini.dem_appr = getMicroMagDemagApproximation('none');
MySim_ini.alpha = 4.42e3;
MySim_ini.gamma = 0;

%initial magnetization
MySim_ini.m0(:) = 1/sqrt(3);

%time grid on which to solve the problem
MySim_ini = MySim_ini.setTime( linspace(0,100e-9,200) );
MySim_ini.setTimeDis = int32(100);
HystDir = 1/mu0*[1,1,1] ;

%time-dependent applied field
HextFct = @(t) (1e-9-t)' .* HystDir .* (t<1e-9)';

MySim_ini = MySim_ini.setHext( HextFct, linspace(0,100e-9,3) );
%     for n=1:3 ; MySim_ini.HextFct{n} = @(t) (1e-9-t)' .* HystDir(n) .* (t<1e-9)' ; end




%% Run the simulation
MySim_ini = MySim_ini.setSolverType( 'UseExplicitSolver' );
% close all
[SigmaInit,~,InteractionMatrices] = ComputeTheSolution(MySim_ini);
dV = (InteractionMatrices.X(2,1,1)-InteractionMatrices.X(1,1,1))*(InteractionMatrices.Y(1,2,1)-InteractionMatrices.Y(1,1,1))*(InteractionMatrices.Z(1,1,2)-InteractionMatrices.Z(1,1,1)) ;
[Mx,My,Mz] = ComputeMagneticMomentGeneralMesh(SigmaInit,repmat(dV,size(SigmaInit,2)/3,1)) ;


figure ; hold on
title('CARTESIAN A') ;
plot([MySim_ini.m0(1),Mx],'ro-') ;
plot([MySim_ini.m0(2),My],'go-')
plot([MySim_ini.m0(3),Mz],'bo-')
set(gca,'ylim',[-.1,1.1]) ; grid on
end