function [elapsedTime_part1,elapsedTime_part2,problem_ini,solution_ini,problem_dym,solution_dym] = Standard_problem_4( NIST_field, resolution, use_CUDA, ShowTheResult, SaveTheResult )
clearvars -except NIST_field resolution use_CUDA SaveTheResult ShowTheResult
% close all
%--- Use either field 1 or field 2 from the mumag example
if ~exist('NIST_field','var')
    NIST_field = 1;
end
if ~exist('resolution','var')
    resolution = [1*36,1*36,1];
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
    figure3= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig3 = axes('Parent',figure3,'Layer','top','FontSize',16); hold on; grid on; box on
end
mu0 = 4*pi*1e-7;
addpath('../../../MEX_files');
addpath('../../../util');
%% --------------------------------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------- MAGTENSE ---------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
%% Setup problem for the time-dependent solver
problem_dym = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
problem_dym.grid_L = [500e-9,500e-9,3e-9];%m
problem_dym.alpha = 4.42e3 ;
problem_dym.gamma = 2.21e5 ;
problem_dym = problem_dym.setUseCuda( use_CUDA );
problem_dym = problem_dym.setUseCVODE( use_CVODE );
problem_dym.dem_appr = getMicroMagDemagApproximation('none');
problem_dym = problem_dym.setTime( linspace(0,2e-9,500) );
problem_dym.setTimeDis = int32(10);

HystDir = 1/mu0*[-25,5,0]/1000 ;
HextFct = @(t) (t>-1)' .*HystDir;
problem_dym = problem_dym.setHext( HextFct, linspace(0,2e-9,2000) );

rng('default');
for i = 1:length(problem_dym.m0(:,1))
    % temp_vec = [2*rand(2,1)-1; 0];
    temp_vec = 2*rand(3,1)-1;
    problem_dym.m0(i,:) = 1/norm(temp_vec)*temp_vec;
end

%convert the class obj to a struct so it can be loaded into fortran
solution_dym = struct();
prob_struct = struct(problem_dym);
tic
solution_dym = problem_dym.MagTenseLandauLifshitzSolver_mex( prob_struct, solution_dym );
elapsedTime_part2 = toc
if (ShowTheResult)
    plot(fig1,solution_dym.t,mean(solution_dym.M(:,:,1),2),'rx');
    plot(fig1,solution_dym.t,mean(solution_dym.M(:,:,2),2),'gx');
    plot(fig1,solution_dym.t,mean(solution_dym.M(:,:,3),2),'bx');
    M_1 = squeeze(solution_dym.M(1,:,:)); figure(figure3); subplot(2,1,1); quiver(solution_dym.pts(:,1),solution_dym.pts(:,2),M_1(:,1),M_1(:,2)); axis equal; title('Fortran starting magnetization')
    M_end = squeeze(solution_dym.M(end,:,:)); figure(figure3); subplot(2,1,2); quiver(solution_dym.pts(:,1),solution_dym.pts(:,2),M_end(:,1),M_end(:,2)); axis equal; title('Fortran ending magnetization')

    % figure; hold all; axis equal; for i=1:200; M_i = squeeze(solution_dym.M(i,:,:)); quiver(solution_dym.pts(:,1),solution_dym.pts(:,2),M_i(:,1),M_i(:,2)); pause(1); cla;
end

end