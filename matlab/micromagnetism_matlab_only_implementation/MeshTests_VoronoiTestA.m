function MeshTests_VoronoiTestA
% close all ;
disp('Running Matlab model')

%% Calculates dynamic solution to mumag std prob 4

mu0 = 4*pi*1e-7;

%% addpaths
addpath('../MEX_files');
addpath('../util');

%% TetraMesh
thisGridL = [125e-9,125e-9/2,125e-9/4] ;
n_Tc = 7 ; % number of grains (i.e. number of Voronoi cells)
%     base_res = 1.*[4,4,4] ;
base_res = [1*9,1*7,5] ;
L = 1 ; % number of refinment iterations from base_res


[voronoi_map] = generateVoronoiRefinedDomain2( n_Tc,   base_res, L,thisGridL);
plot_grid(voronoi_map.pos, voronoi_map.dims, voronoi_map.ndim, voronoi_map.Tc_distr, voronoi_map.vor_c, voronoi_map.dTc)
VoronoiMeshFileName = 'TestVoronoiMesh01.mat' ;
save( VoronoiMeshFileName,'voronoi_map') ;


%% Create Preset

resolution = [size(voronoi_map.pos,1),1,1];
MySim_ini = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
%         MySim_ini.A0 = 0  ;
%         MySim_ini.Ms = 1e-3 ;
MySim_ini.grid_L = thisGridL;% NOT STANDARD
MySim_ini.RecomputeInteractionMatrices = 1 ;
MySim_ini.ExternalMesh = 1 ;
MySim_ini.MeshType = 'Voronoi' ;
MySim_ini.ExternalMeshFileName = VoronoiMeshFileName ;
MySim_ini.DemagTensorFileName = 'VoronoiInteractionMatrices.mat' ;
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

MySim_ini = MySim_ini.setHext( HextFct, linspace(0,100e-9,3) ); % was 2000 !
%     for n=1:3 ; MySim_ini.HextFct{n} = @(t) (1e-9-t)' .* HystDir(n) .* (t<1e-9)' ; end




%% Run the simulation
MySim_ini = MySim_ini.setSolverType( 'UseExplicitSolver' );

% close all
[SigmaInit,~,InteractionMatrices] = ComputeTheSolution(MySim_ini);
% figure; quiver3(InteractionMatrices.X(:),InteractionMatrices.Y(:),InteractionMatrices.Z(:),MySim_ini.m0(:,1),MySim_ini.m0(:,2),MySim_ini.m0(:,3)); axis equal;  title('Starting state - Matlab')

[Mx,My,Mz] = ComputeMagneticMomentGeneralMesh(SigmaInit,prod(voronoi_map.dims,2)) ;
%     figure; quiver3(InteractionMatrices.X(:),InteractionMatrices.Y(:),InteractionMatrices.Z(:),SigmaX(:),SigmaY(:),SigmaZ(:)); axis equal;  title('Starting state - Matlab')
figure ; hold on
title('VORONOI A') ;
plot([MySim_ini.m0(1),Mx],'ro-') ;
plot([MySim_ini.m0(2),My],'go-')
plot([MySim_ini.m0(3),Mz],'bo-')
set(gca,'ylim',[-.1,1.1]) ; grid on
end