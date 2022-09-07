function [elapsedTime_part1,elapsedTime_part2,problem_ini,solution_ini,problem_dym,solution_dym,mesh,GridInfo] = Standard_problem_4_unstructured_cart( NIST_field, n, L, mesh_scale )

clearvars -except NIST_field resolution use_CUDA SaveTheResult ShowTheResult weight_direct n L mesh_scale
% close all

%--- Use either field 1 or field 2 from the NIST example
if ~exist('NIST_field','var')
    NIST_field = 1;
end
if ~exist('resolution','var')
    resolution = [1*36,1*9,1];
end
if ~exist('use_CUDA','var')
    use_CUDA = true;
end
if ~exist('ShowTheResult','var')
    ShowTheResult = 1;
end
if ~exist('SaveTheResult','var')
    SaveTheResult = 0;
end

if use_CUDA
    dir_m = 'GPU';
else
    dir_m = 'CPU';
end
use_CVODE = false;

if (ShowTheResult)
    figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
end
mu0 = 4*pi*1e-7;

addpath('../../../MEX_files');
addpath('../../../util');
addpath('../../../micromagnetism_matlab_only_implementation');

%% --------------------------------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------- MAGTENSE ---------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
load('Std_prob_4_unstructured_mesh_grains_6_res_80_20_ref_2.mat');
% load('Std_prob_4_unstructured_mesh_grains_6_res_100_25_ref_3.mat');
cartesianUnstructuredMeshPlot(mesh.pos_out,mesh.dims_out,GridInfo,mesh.iIn);

%--- Setup the problem
resolution = [length(mesh.pos_out) 1 1];
disp(['Prisms N_grid = ' num2str(prod(resolution))])
problem_ini = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
problem_ini = problem_ini.setMicroMagGridType('unstructuredPrisms');

%--- Information on the grid
problem_ini.grid_pts    = mesh.pos_out;
problem_ini.grid_abc    = mesh.dims_out;

%--- Calculate the exchange matrix
InteractionMatrices.GridInfo = GridInfo;
InteractionMatrices.X = GridInfo.Xel ;
InteractionMatrices.Y = GridInfo.Yel ;
InteractionMatrices.Z = GridInfo.Zel ;

[D2X,D2Y] = computeDifferentialOperatorsFromMesh_DirectLap(GridInfo,'extended',8,"DirectLaplacianNeumann");
InteractionMatrices.A2 = D2X + D2Y ;

%--- Convert the exchange matrix to sparse
problem_ini = problem_ini.setExchangeMatrixSparse( InteractionMatrices.A2 );

%% Setup the problem for the initial configuration
problem_ini = problem_ini.setMicroMagDemagApproximation('none');
problem_ini = problem_ini.setUseCuda( use_CUDA );
problem_ini = problem_ini.setUseCVODE( use_CVODE );

problem_ini.alpha = 4.42e3;
problem_ini.gamma = 0;

%initial magnetization
problem_ini.m0(:) = 1/sqrt(3);
    
%time grid on which to solve the problem
problem_ini = problem_ini.setTime( linspace(0,100e-9,200) );
problem_ini.setTimeDis = int32(100);
HystDir = 1/mu0*[1,1,1] ;

%time-dependent applied field
HextFct = @(t) (1e-9-t)' .* HystDir .* (t<1e-9)';
problem_ini = problem_ini.setHext( HextFct, linspace(0,100e-9,2000) );


%% Create Preset
solution_ini = struct();
%convert the class obj to a struct so it can be loaded into fortran
prob_struct = struct(problem_ini);

tic
solution_ini = problem_ini.MagTenseLandauLifshitzSolver_mex( prob_struct, solution_ini );
% [Mx,My,Mz,mx,my,mz] = computeMagneticMomentGeneralMesh(solution_ini.M,GridInfo.Volumes) ;
elapsedTime_part1 = toc
if (ShowTheResult)
    figure; M_end = squeeze(solution_ini.M(end,:,:)); quiver(solution_ini.pts(:,1),solution_ini.pts(:,2),M_end(:,1),M_end(:,2)); axis equal; title('Starting state - Fortran')
end

%% Setup problem for the time-dependent solver
problem_dym = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
problem_dym = problem_dym.setMicroMagGridType('unstructuredPrisms');

%--- Information on the grid
problem_dym.grid_pts    = problem_ini.grid_pts;
problem_dym.grid_abc    = problem_ini.grid_abc ;

%--- Convert the exchange matrix to sparse
problem_dym.exch_nval = problem_ini.exch_nval;
problem_dym.exch_nrow = problem_ini.exch_nrow;
problem_dym.exch_val = problem_ini.exch_val;
problem_dym.exch_rows = problem_ini.exch_rows;
problem_dym.exch_rowe = problem_ini.exch_rowe;
problem_dym.exch_col = problem_ini.exch_col;

problem_dym.alpha = 4.42e3 ;
problem_dym.gamma = 2.21e5 ;
problem_dym = problem_dym.setUseCuda( use_CUDA );
problem_dym = problem_dym.setUseCVODE( use_CVODE );
problem_dym = problem_dym.setMicroMagDemagApproximation('none');
problem_dym = problem_dym.setTime( linspace(0,1e-9,200) );
problem_dym.setTimeDis = int32(10);

if (NIST_field == 1)
    %field 1
    HystDir = 1/mu0*[-24.6,4.3,0]/1000 ;
end
if (NIST_field == 2)
    %field 2
    HystDir = 1/mu0*[-35.5,-6.3,0]/1000 ;
end

HextFct = @(t) (t>-1)' .*HystDir;
problem_dym = problem_dym.setHext( HextFct, linspace(0,1e-9,2000) );

problem_dym.m0(:) = solution_ini.M(end,:,:);

%convert the class obj to a struct so it can be loaded into fortran
solution_dym = struct();
prob_struct = struct(problem_dym);

tic
solution_dym = problem_dym.MagTenseLandauLifshitzSolver_mex( prob_struct, solution_dym );
elapsedTime_part2 = toc

[Mx,My,Mz,mx,my,mz] = computeMagneticMomentGeneralMesh(solution_dym.M,GridInfo.Volumes) ;
if (ShowTheResult)
    plot(fig1,solution_dym.t,Mx,'rx');
    plot(fig1,solution_dym.t,My,'gx');
    plot(fig1,solution_dym.t,Mz,'bx');
end

%% --------------------------------------------------------------------------------------------------------------------------------------
%% --------------------------------------------------------------------  mumag -----------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
%% Compare with published solutions available from mumag webpage
if (ShowTheResult)
    t=linspace(0,1,1000);
    data = load(['Published_solutions_field' num2str(NIST_field)]);
    mean_avg=data.mean_avg; err=data.err;
    colours = [[1 0 0];[0 1 0];[0 0 1]];
    weak_colours = colours + ~colours*0.75;
    fill_ts=[t,fliplr(t)];  
    for j=1:3
        std_errors{j}(1:2,:)=[mean_avg(1,:,j)+err(1,:,j);mean_avg(1,:,j)-err(1,:,j)];
        interval = [std_errors{j}(1,:),fliplr(std_errors{j}(2,:))];
        fill(fig1,1e-9*fill_ts,interval,weak_colours(j,:),'linestyle','none')
        plot(fig1,1e-9*t,mean_avg(1,:,j),'color',colours(j,:))
    end

    legend(fig1,'Fortran Mx','Fortran My','Fortran Mz','Matlab Mx','Matlab My','Matlab Mz','\mu{}mag \sigma{}(Mx)','\mu{}mag <Mx>','\mu{}mag \sigma{}(My)','\mu{}mag <My>','\mu{}mag \sigma{}(Mz)','\mu{}mag <Mz>');
    ylabel(fig1,'<M_i>/M_s')
    xlabel(fig1,'Time [ns]')

    xlim(fig1,[0 1e-9])
    figure(figure1)
end

end