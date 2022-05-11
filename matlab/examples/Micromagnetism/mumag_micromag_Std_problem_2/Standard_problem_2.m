function [elapsedTime,problem,solution,results] = Standard_problem_2(resolution, use_CUDA, ShowTheResult, SaveTheResult, run_single_curve, run_MrHc )

clearvars -except resolution use_CUDA SaveTheResult ShowTheResult run_single_curve run_MrHc
% close all

if ~exist('resolution','var')
    resolution = [20,4,1];
end
if ~exist('use_CUDA','var')
    use_CUDA = false;
end
if ~exist('ShowTheResult','var')
    ShowTheResult = 1;
end
if ~exist('SaveTheResult','var')
    SaveTheResult = 0;
end
if ~exist('run_single_curve','var')
    run_single_curve = 1;
end
if ~exist('run_MrHc','var')
    run_MrHc = 0;
end

if use_CUDA
    dir_m = 'GPU';
else
    dir_m = 'CPU';
end

if (ShowTheResult)
    figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
end
mu0 = 4*pi*1e-7;

addpath('../../../MEX_files');
addpath('../../../util');

%% --------------------------------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------- MAGTENSE ---------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
%% Setup the problem
problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));

problem.dem_appr = getMicroMagDemagApproximation('none');
problem = problem.setUseCuda( use_CUDA );
problem.alpha = 1e3;
% problem.gamma = 2.21e5;

MaxH = 0.1;

%time-dependent applied field
HystDir = 1/mu0*[1,1,1]/sqrt(3) ;
HextFct = @(t) HystDir .* t';

%initial magnetization
problem.m0(:) = 1/sqrt(3);

problem = problem.setSolverType( 'UseExplicitSolver' );
problem.solver = getMicroMagSolver( 'Explicit' );

problem = problem.setHext( HextFct, linspace(MaxH,-MaxH,40) );
problem = problem.setTime( linspace(0,40e-9,2) );
problem = problem.setConvergenceCheckTime( linspace(0,40e-9,2) );
problem.conv_tol = 1e-6;
    
problem.K0 = 0 ;
problem.Ms = 1000e3 ;
problem.A0 = 1.74532925199e-10;

%time-dependent alpha parameter, to ensure faster convergence
AlphaFct = @(t) problem.alpha * 10.^( 5 * min(t,2e-9)/2e-9 );
problem = problem.setAlpha( AlphaFct, linspace(0,10e-9,100) );
problem.alpha = 0;

if run_single_curve
    %% Run a single hysterhesis curve for standard problem 2
    d_loop = 0.5;
end
if run_MrHc
    d_loop = linspace(0.05,0.5,10);
end

tic
for i = 1:length(d_loop)

    problem.grid_L = [5e-6,1e-6,1e-7]*d_loop(i);%m
    results.dlex(i) = problem.grid_L(2)/sqrt(problem.A0/(1/2*mu0*problem.Ms^2));
    
    solution = struct();
    prob_struct = struct(problem);  %convert the class obj to a struct so it can be loaded into fortran

    disp(['Running d/l_ex = ' num2str(results.dlex(i)) ', i.e. ' num2str(i) '/' num2str(length(d_loop))])
    
    solution = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );

    for j = 1:problem.nt_Hext 
        Mx_arr = solution.M(end,:,j,1) ;
        My_arr = solution.M(end,:,j,2) ;
        Mz_arr = solution.M(end,:,j,3) ;
        MN = sqrt(Mx_arr.^2+My_arr.^2+Mz_arr.^2) ;
        Mx(j) = mean(Mx_arr./MN) ;
        My(j) = mean(My_arr./MN) ;
        Mz(j) = mean(Mz_arr./MN) ;
        M(j) = Mx(j) + My(j) + Mz(j) ;
    end
    if (ShowTheResult)
        plot(fig1,sign(problem.Hext(:,1)).*sqrt(problem.Hext(:,2).^2+problem.Hext(:,3).^2+problem.Hext(:,4).^2)/problem.Ms,M,'rp') %Minus signs added to correspond to regular hysteresis plots.
    end
    
    results.Mxr(i) = interp1(sign(problem.Hext(:,1)).*sqrt(problem.Hext(:,2).^2+problem.Hext(:,3).^2+problem.Hext(:,4).^2)/problem.Ms,Mx,0);
    results.Myr(i) = interp1(sign(problem.Hext(:,1)).*sqrt(problem.Hext(:,2).^2+problem.Hext(:,3).^2+problem.Hext(:,4).^2)/problem.Ms,My,0);
    results.Hc(i)  = interp1(M,sign(problem.Hext(:,1)).*sqrt(problem.Hext(:,2).^2+problem.Hext(:,3).^2+problem.Hext(:,4).^2)/problem.Ms,0);
end
elapsedTime = toc

if run_MrHc
    if (ShowTheResult)
        figure; plot(results.dlex,results.Mxr,'.'); xlabel('d/l_{ex}'); ylabel('M_{xr}/M_s');
        figure; plot(results.dlex,results.Myr,'.'); xlabel('d/l_{ex}'); ylabel('M_{yr}/M_s');
        figure; plot(results.dlex,results.abs(Hc),'.'); xlabel('d/l_{ex}'); ylabel('|H_c|/M_s');
    end
end

if run_single_curve
%% --------------------------------------------------------------------------------------------------------------------------------------
%% --------------------------------------------------------------------  mumag -----------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
%% Compare with published solutions available from mumag webpage for single curve for d/l_ex = 30
    load('OOMMF_Hysteresis2D_dlex30.mat');
    plot(fig1,mu0*H,M,'k>');
    load('OOMMF_HysteresisQuasi3D_dlex30.mat');
    plot(fig1,mu0*H,M,'k<');
    load('OOMMF_Hysteresis3D_dlex30.mat');
    plot(fig1,mu0*H,M,'k^');

    % legend(fig1,'Fortran Mx','Fortran My','Fortran Mz','Matlab Mx','Matlab My','Matlab Mz','\mu{}mag \sigma{}(Mx)','\mu{}mag <Mx>','\mu{}mag \sigma{}(My)','\mu{}mag <My>','\mu{}mag \sigma{}(Mz)','\mu{}mag <Mz>');
    legend(fig1,'"Fortran Explicit method"','OOMMF 2D','OOMMF Quasi3D','OOMMF 3D','Location','SouthEast');
    ylabel(fig1,'<M_i>/M_s')
    xlabel(fig1,'\mu_{0}H_{applied} [T]')
    xlim(fig1,[-0.1 0.1])
    figure(figure1)
end

end