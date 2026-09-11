function [elapsedTime,problem,solution,results] = Standard_problem_2(resolution, d_loop, options)

arguments
    resolution (1,3) {mustBeInteger}                = [20,4,1];     %--- [nx,ny,nz] of the grid
    d_loop (1,:) {mustBeNumeric}                    = linspace(0.05,0.5,10); %--- The values of d to run the model for
    options.use_CUDA {mustBeNumericOrLogical}       = true          %--- Use CUDA for the calculations
    options.use_CVODE {mustBeNumericOrLogical}      = false;        %--- Use CVODE for the numerical time evolution
    options.ShowTheResult {mustBeNumericOrLogical}  = true          %--- Show the result
end

if (isscalar(d_loop) && d_loop == 0.5)
    run_single_curve = 1;
else
    run_single_curve = 0;
end

mu0 = 4*pi*1e-7;

addpath('../../../MEX_files');
addpath('../../../util');

%% --------------------------------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------- MAGTENSE ---------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
%% Setup the problem
problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));

problem = problem.setMicroMagDemagApproximation('none');
problem = problem.setUseCuda( options.use_CUDA );
problem = problem.setUseCVODE( options.use_CVODE );
problem.alpha = 1e3;
% problem.gamma = 2.21e5;

MaxH = 0.1;

%time-dependent applied field
HystDir = 1/mu0*[1,1,1]/sqrt(3) ;
HextFct = @(t) HystDir .* t';

%initial magnetization
problem.m0(:) = 1/sqrt(3);

problem = problem.setSolverType( 'UseExplicitSolver' );
problem = problem.setMicroMagSolver( 'Explicit' );

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
        M(j) = Mx(j)*HystDir(1) + My(j)*HystDir(2) + Mz(j)*HystDir(3) ;
    end
    
    results.Mxr(i) = interp1(sign(problem.Hext(:,1)).*sqrt(problem.Hext(:,2).^2+problem.Hext(:,3).^2+problem.Hext(:,4).^2)/problem.Ms,Mx,0);
    results.Myr(i) = interp1(sign(problem.Hext(:,1)).*sqrt(problem.Hext(:,2).^2+problem.Hext(:,3).^2+problem.Hext(:,4).^2)/problem.Ms,My,0);
    results.Hc(i)  = interp1(M,sign(problem.Hext(:,1)).*sqrt(problem.Hext(:,2).^2+problem.Hext(:,3).^2+problem.Hext(:,4).^2)/problem.Ms,0);
end
elapsedTime = toc


if (options.ShowTheResult)
    figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
    figure2= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig2 = axes('Parent',figure2,'Layer','top','FontSize',16); hold on; grid on; box on
    figure3= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig3 = axes('Parent',figure3,'Layer','top','FontSize',16); hold on; grid on; box on
    
    plot(fig1,results.dlex,results.Mxr,'k.','Markersize',20); xlabel(fig1,'d/l_{ex}'); ylabel(fig1,'M_{xr}/M_s');
    plot(fig2,results.dlex,results.Myr,'k.','Markersize',20); xlabel(fig2,'d/l_{ex}'); ylabel(fig2,'M_{yr}/M_s');
    plot(fig3,results.dlex,abs(results.Hc),'k.','Markersize',20); xlabel(fig3,'d/l_{ex}'); ylabel(fig3,'|H_c|/M_s');
    
    %% --------------------------------------------------------------------------------------------------------------------------------------
    %% --------------------------------------------------------------------  mumag -----------------------------------------------------------
    %% --------------------------------------------------------------------------------------------------------------------------------------
    mumag_data_names = {'Streibl','McMichael','Lopez-Diaz','Donahue'};
    for i = 1:length(mumag_data_names)
        mumag_data = load(['../../../../documentation/examples_mumag_validation/Validation_standard_problem_2/' mumag_data_names{i} '.txt']);
        plot(fig1,mumag_data(:,1),mumag_data(:,2),'d');
        plot(fig2,mumag_data(:,1),mumag_data(:,3),'d');
        plot(fig3,mumag_data(:,1),mumag_data(:,4),'d');
    end
end

if run_single_curve
%% --------------------------------------------------------------------------------------------------------------------------------------
%% --------------------------------------------------------------------  mumag -----------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
%% Compare with published solutions available from mumag webpage for single curve for d/l_ex = 30
    figure4= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig4 = axes('Parent',figure4,'Layer','top','FontSize',16); hold on; grid on; box on
    plot(fig1,sign(problem.Hext(:,1)).*sqrt(problem.Hext(:,2).^2+problem.Hext(:,3).^2+problem.Hext(:,4).^2)/problem.Ms,mu0*M,'rp') %Minus signs added to correspond to regular hysteresis plots.

    load('OOMMF_Hysteresis2D_dlex30.mat');
    plot(fig4,mu0*H,M,'k>');
    load('OOMMF_HysteresisQuasi3D_dlex30.mat');
    plot(fig4,mu0*H,M,'k<');
    load('OOMMF_Hysteresis3D_dlex30.mat');
    plot(fig4,mu0*H,M,'k^');

    legend(fig4,'"MagTense"','OOMMF 2D','OOMMF Quasi3D','OOMMF 3D','Location','SouthEast');
    ylabel(fig4,'<M_i>/M_s')
    xlabel(fig4,'\mu_{0}H_{applied} [T]')
    xlim(fig4,[-0.1 0.1])
end

end