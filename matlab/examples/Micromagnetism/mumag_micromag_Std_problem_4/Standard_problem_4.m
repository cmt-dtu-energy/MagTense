function [elapsedTime_part1,elapsedTime_part2,problem_ini,solution_ini,problem_dym,solution_dym] = Standard_problem_4( NIST_field, resolution, use_CUDA, ShowTheResult )

clearvars -except NIST_field resolution use_CUDA ShowTheResult
% close all

%--- Use either field 1 or field 2 from the mumag example
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

use_CVODE = false;

if (ShowTheResult)
    figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
end
mu0 = 4*pi*1e-7;

addpath('../../../MEX_files');
addpath('../../../util');

%% --------------------------------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------- MAGTENSE ---------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
%% Setup the problem for the initial configuration
%takes the size of the grid as arguments (nx,ny,nz) and a function handle
%that produces the desired field (if not present zero applied field is inferred)
problem_ini = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
problem_ini.grid_L = [500e-9,125e-9,3e-9];%m
problem_ini.nThreads = int32(8);

problem_ini.dem_appr = getMicroMagDemagApproximation('none');
problem_ini = problem_ini.setUseCuda( use_CUDA );
problem_ini = problem_ini.setUseCVODE( use_CVODE );

% loadFile = 'N_out_test.dat';
% if exist(loadFile,'file')
%     problem_ini = problem_ini.setLoadNFilename( loadFile );
% end
% problem_ini = problem_ini.setReturnNFilename( loadFile );

problem_ini.Ms = 8e5*ones(prod(resolution),1);
problem_ini.K0 = 0*zeros(prod(resolution),1);
problem_ini.A0 = 1.3e-11;

problem_ini.alpha = 4.42e3;
problem_ini.gamma = 0;

%initial magnetization
problem_ini.m0(:,1) = 1/sqrt(3);
problem_ini.m0(:,2) = 1/sqrt(3);
problem_ini.m0(:,3) = 1/sqrt(3);

% problem_ini.u_ea(:,:) = 0;
    
%time grid on which to solve the problem
problem_ini = problem_ini.setTime( linspace(0,100e-9,200) );
problem_ini.setTimeDis = int32(100);
HystDir = 1/mu0*[1,1,1] ;

%time-dependent applied field
HextFct = @(t) (1e-9-t)' .* HystDir .* (t<1e-9)';
problem_ini = problem_ini.setHext( HextFct, linspace(0,100e-9,2000) );


problem_ini = problem_ini.setConvergenceCheckTime( linspace(0,40e-9,2) );
% problem_ini.conv_tol = 1e-6;


%% Create Preset
solution_ini = struct();
%convert the class obj to a struct so it can be loaded into fortran
prob_struct = struct(problem_ini);

tic
solution_ini = problem_ini.MagTenseLandauLifshitzSolver_mex( prob_struct, solution_ini );
elapsedTime_part1 = toc
if (ShowTheResult)
    figure; M_end = squeeze(solution_ini.M(end,:,:)); quiver(solution_ini.pts(:,1),solution_ini.pts(:,2),M_end(:,1),M_end(:,2)); axis equal; title('Starting state - Fortran')
end

%% Setup problem for the time-dependent solver
problem_dym = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
problem_dym.alpha = 4.42e3 ;
problem_dym.gamma = 2.21e5 ;
problem_dym = problem_dym.setUseCuda( use_CUDA );
problem_dym = problem_dym.setUseCVODE( use_CVODE );
problem_dym.dem_appr = getMicroMagDemagApproximation('none');
problem_dym = problem_dym.setTime( linspace(0,1e-9,200) );
problem_dym.setTimeDis = int32(10);
problem_dym.Ms = 8e5*ones(prod(resolution),1);
problem_dym.K0 = 0*zeros(prod(resolution),1);
problem_dym.A0 = 1.3e-11;
problem_dym.nThreads = problem_ini.nThreads;

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

if (ShowTheResult)
    plot(fig1,solution_dym.t,mean(solution_dym.M(:,:,1),2),'rx'); 
    plot(fig1,solution_dym.t,mean(solution_dym.M(:,:,2),2),'gx'); 
    plot(fig1,solution_dym.t,mean(solution_dym.M(:,:,3),2),'bx'); 
    % figure; hold all; for i=2:4; plot(problem.Hext(:,1),problem.Hext(:,i),'.'); end;
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