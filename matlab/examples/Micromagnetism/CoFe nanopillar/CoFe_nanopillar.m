function [solution, M, Hext, Hc, Mr, prob_struct] = CoFe_nanopillar(dimensions, resolution, HystDir, field_steps, options)

arguments
    dimensions (1,3) {mustBeNumeric}                = [250e-9, 1500e-9, 30];      %--- [Radius, height, phi angle] of the hexagon
    resolution (1,3) {mustBeInteger}                = [10,10,5];                  %--- [nx,ny,nz] of the grid
    HystDir (1,3) {mustBeNumeric}                   = [0.01,0.01,1]     %--- Direction of the applied field
    field_steps (1,1) {mustBeNumeric}               = 0.1               %--- Size of the field steps
    options.use_CUDA {mustBeNumericOrLogical}       = true              %--- Use CUDA for the calculations
    options.use_CVODE {mustBeNumericOrLogical}      = false;            %--- Use CVODE for the numerical time evolution
end


%% DEFINE THE GEOMETRY AND GET THE MESHINFO
addpath('../../../../../MagTense/matlab/MEX_files');
addpath('../../../../../MagTense/matlab/util');
addpath('../../../../matlab/util');

fnameSave = ['CoFe_piller_' sprintf('%03.3i',dimensions(1)) '_' ...
            sprintf('%03.3i',dimensions(2)) '_' ...
            sprintf('%03.3i',dimensions(3)) '_' ...
            sprintf('%03.3i',resolution(1)) '_' ...
            sprintf('%03.3i',resolution(2)) '_' ...
            sprintf('%03.3i',resolution(3)) '_' ...
            sprintf('%03.3i',HystDir(1)*100) '_' ...
            sprintf('%03.3i',HystDir(2)*100) '_' ...
            sprintf('%03.3i',HystDir(3)*100) '_' ...
            sprintf('%03.3i',field_steps*100)];

info.resolution = resolution;
info.dimensions = dimensions;

easy_axis = [0,0,1];
easy_axisNorm = sqrt(easy_axis(1).^2+easy_axis(2).^2+easy_axis(3).^2) ;
easy_axis = easy_axis/easy_axisNorm;

HystDirNorm = sqrt(HystDir(1).^2+HystDir(2).^2+HystDir(3).^2) ;
HystDir = HystDir./HystDirNorm;
info.HystDir = HystDir;

%% Physical Parameters

gamma = 0;
alpha = 4000 ; 

mu0 = 4*pi*1e-7;

% ---- In the grains:
Ms = 1.44*10^6 ;     % A/m
A0 = 1.3*10^(-11) ; % J/m3

[mesh_cart,GridInfo,mesh_params] = CreateHexagonCartesianMesh_R_h_phi(dimensions, resolution);

%% Problem structure creation

problem = DefaultMicroMagProblem(mesh_params.resolution(1),mesh_params.resolution(2),mesh_params.resolution(3));
problem = problem.setUseCuda( options.use_CUDA );
problem = problem.setUseCVODE( options.use_CVODE );
problem = problem.setMicroMagDemagApproximation('none');
problem = problem.setMicroMagGridType('unstructuredPrisms');
problem = problem.setSolverType( 'UseExplicitSolver' );
problem = problem.setMicroMagSolver( 'Explicit' );

problem.grid_pts = mesh_cart.pos_out;
problem.grid_abc = mesh_cart.dims_out;

problem.nThreads = int32(8);

%% Save the physical parameters in the problem structure

problem.alpha = alpha;
problem.gamma = gamma;
problem.Ms = Ms;
problem.K0 = 0;
problem.K1 = 2.7*10^3; % J/m3
problem.K2 = 2.7*10^3; % J/m3
problem.A0 = A0;

%% Exchange matrix

InteractionMatrices.GridInfo = GridInfo;
InteractionMatrices.X = GridInfo.Xel ;
InteractionMatrices.Y = GridInfo.Yel ;
InteractionMatrices.Z = GridInfo.Zel ;

% Generate Laplace operator:
[D2X,D2Y,D2Z] = computeDifferentialOperatorsFromMesh_DirectLap(GridInfo, 'extended',6,"DirectLaplacianNeumann");
InteractionMatrices.A2 = D2X + D2Y + D2Z ;

%--- Convert the exchange matrix to sparse
problem = problem.setExchangeMatrixSparse( InteractionMatrices.A2 );

%% Applied field
problem.m0(:,1) = HystDir(1) ;
problem.m0(:,2) = HystDir(2) ;
problem.m0(:,3) = HystDir(3) ;

%time-dependent applied field
HextFct = @(t) HystDir/mu0 .* t';

% time grid on which to solve the problems
problem = problem.setTime( linspace(0,4e-9,2) );
problem = problem.setConvergenceCheckTime( linspace(0,4e-9,2) );

problem.conv_tol = 1e-6;

% time-dependent applied field
Hvalues = sort(unique([1:-0.05:-1]),'descend') ; % [T]
problem = problem.setHext( HextFct, Hvalues );

%% RUN !!
save([fnameSave,'_Problem.mat'],'problem','GridInfo','mesh_cart','mesh_params','info','-v7.3')

tic
solution = struct();
prob_struct = struct(problem);  %convert the class obj to a struct so it can be loaded into fortran

solution = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );

elapsedTime = toc

save([fnameSave,'_Solution.mat'],'solution','-v7.3')

%% Compute Hysteresis Loop

[Mx,My,Mz,mx,my,mz] = computeMagneticMomentGeneralMesh(solution.M,GridInfo.Volumes) ;

M = problem.Ms.*(HystDir(1).*Mx + HystDir(2).*My + HystDir(3).*Mz) ;
H = squeeze(sign(problem.Hext(:,1)).*sqrt(problem.Hext(:,2).^2+problem.Hext(:,3).^2+problem.Hext(:,4).^2)*mu0) ;

solution = rmfield(solution,'H_exc');
solution = rmfield(solution,'H_ext');
solution = rmfield(solution,'H_dem');
solution = rmfield(solution,'H_ani');

try  Hc = interp1(M,H,0) ; catch ; Hc = nan ; end
save([fnameSave,'_Solution.mat'],'problem','solution','Hc','H','M','-v7.3') ;

%% Plot
figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
quiver3(fig1,solution.pts(:,1),solution.pts(:,2),solution.pts(:,3),squeeze(solution.M(end,:,end,1))',squeeze(solution.M(end,:,end,2))',squeeze(solution.M(end,:,end,3))')

makePDFResults_Hexagon([fnameSave '_']); 

%% Get the field from the tip of the nanopillar and specifically at 100 nm away as function of applied field
z = dimensions(2)/2:1e-9:(dimensions(2)/2+200e-9);

pts = zeros( numel(z), 3 );
pts(:,3) = z;

%%Get a default tile from MagTense
tile = getDefaultMagTile();
    
%ensure the tile is a permanent magnet
tile.magnetType = getMagnetType('hard');

%set the geometry to be a rectangular prism
tile.tileType = getMagTileType('prism');

figure2= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig2 = axes('Parent',figure2,'Layer','top','FontSize',16); hold on; grid on; box on
figure3= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig3 = axes('Parent',figure3,'Layer','top','FontSize',16); hold on; grid on; box on
color_arr = turbo(length(solution.M(end,1,:,1)));

for j = 1:length(solution.M(end,1,:,1))
    clear tiles
    for i = 1:length(mesh_cart.pos_out(:,1))
        tiles(i) = tile;
        tiles(i).abc = mesh_cart.dims_out(i,:);
        tiles(i).offset = mesh_cart.pos_out(i,:);
        tiles(i).M = problem.Ms.*solution.M(end,i,j,:);
    end
    
    %get the field
    H = getHFromTiles_mex( tiles, pts, int32( length(tiles) ), int32( length(pts(:,1)) ) );
    plot(fig2,z-dimensions(2)/2,4*pi*1e-7*H(:,3)+Hvalues(j),'.','color',color_arr(j,:))
    
    H_100nm(j,1) = 4*pi*1e-7*H(101,3);
    H_100nm(j,2) = 4*pi*1e-7*H(101,3)+Hvalues(j);
end
xlabel(fig2,'Distance from tip of cylinder [m]')
ylabel(fig2,'\mu_0H_{z,CoFe} + \mu_0H_{z,applied} [T]')
title(fig2,'Color indicate value of applied field')
figure(figure2)
% print('-dpng',[fnameSave '_field_z.png'])

plot(fig3,Hvalues,H_100nm(:,1),'r.','MarkerSize',20)
plot(fig3,Hvalues,H_100nm(:,2),'b.','MarkerSize',20)
xlabel(fig3,'Applied field [T]')
ylabel(fig3,'\mu_0H_{z} [T]')
legend(fig3,{'H_{CoFe}','H_{CoFe} + H_{z,applied}'},'Location','NorthWest')
title(fig3,'Field at 100 nm from tip')
figure(figure3)
% print('-dpng',[fnameSave '_field_100nm.png'])
end