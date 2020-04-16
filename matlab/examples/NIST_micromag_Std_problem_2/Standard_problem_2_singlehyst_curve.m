clearvars
% close all

run_single_curve = 1;
run_MrHc = 0;

%--- Use either field 1 or field 2 from the NIST example
use_CUDA = false;

i_res = 1;
resolution = [i_res*20,i_res*4,1];

figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
mu0 = 4*pi*1e-7;

addpath('../../MEX_files');
addpath('../../util');

problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));

problem.dem_appr = getMicroMagDemagApproximation('none');
problem = problem.setUseCuda( use_CUDA );
problem.alpha = 1e3;
problem.MaxT0 = 2e-9;

MaxH = 0.1;

%time-dependent applied field
HystDir = 1/mu0*[1,1,1]/sqrt(3) ;
HextFct = @(t) HystDir .* t';

%initial magnetization
problem.m0(:) = 1/sqrt(3);

problem = problem.setSolverType( 'UseExplicitSolver' );

problem = problem.setTime( linspace(MaxH,-MaxH,26) );
problem = problem.setHext( HextFct, numel(problem.t) );
problem = problem.setTimeExplicit( linspace(0,1,numel(problem.t)) );

t_explicit = 5e-9 ;
problem = problem.setTime( linspace(0,t_explicit,2) );

problem.solver = getMicroMagSolver( 'Explicit' );
    
problem.K0 = 0 ;
problem.Ms = 1000e3 ;
problem.A0 = 1.74532925199e-10;


AlphaFct = @(t) problem.alpha * 10.^( 5 * min(t,problem.MaxT0)/problem.MaxT0 );
problem = problem.setAlpha( AlphaFct, linspace(0,t_explicit,100) );
problem.alpha = 0;

if run_single_curve
    %% Run a single hysterhesis curve for standard problem 2
    d_loop = 1;
end
if run_MrHc
    d_loop = linspace(0.1,1,20);
end
    
for i = 1:length(d_loop)

    problem.grid_L = [5e-6,1e-6,1e-7]*d_loop(i);%m
    
    dlex(i) = problem.grid_L(2)/sqrt(problem.A0/(1/2*mu0*problem.Ms^2));
    
    solution = struct();
    prob_struct = struct(problem);  %convert the class obj to a struct so it can be loaded into fortran

    tic
    solution = MagTenseLandauLifshitzSolver_mex( prob_struct, solution );
    toc

    for j = 1:problem.nt_Hext 
        Mx_arr = solution.M(end,:,j,1) ;
        My_arr = solution.M(end,:,j,2) ;
        Mz_arr = solution.M(end,:,j,3) ;
        MN = sqrt(Mx_arr.^2+My_arr.^2+Mz_arr.^2) ;
        Mx(j) = mean(Mx_arr./MN) ;
        My(j) = mean(My_arr./MN) ;
        Mz(j) = mean(Mz_arr./MN) ;
        M(j) = Mx(j) + My(j) + Mz(j) ;
        Mk(j) = Mx(j)*HystDir(1) + My(j)*HystDir(2) + Mz(j)*HystDir(3) ;
    end
    plot(fig1,problem.Hext(:,1),mu0*Mk,'rp') %Minus signs added to correspond to regular hysteresis plots.

    Mxr(i) = interp1(problem.Hext(:,1),Mx,0);
    Myr(i) = interp1(problem.Hext(:,1),My,0);
    Hc(i)  = interp1(M,problem.Hext(:,1),0);
end

if run_MrHc
    figure; plot(dlex,Mxr,'.'); xlabel('d/l_{ex}'); ylabel('M_{xr}/M_s');
    figure; plot(dlex,Myr,'.'); xlabel('d/l_{ex}'); ylabel('M_{yr}/M_s');
    figure; plot(dlex,abs(Hc),'.'); xlabel('d/l_{ex}'); ylabel('|H_c|/M_s');
end

if run_single_curve
    %% Run the Matlab version of the micromagnetism code
    addpath('..\..\micromagnetism')
    tic
    Matlab_model_params.nGrid = resolution';
    % Matlab_model_params.Field_dir = HystDir;
    Script_3D_Std_Problem_2(fig1,Matlab_model_params);
    toc

    %% Compare with published solutions available from NIST webpage
    load('OOMMF_Hysteresis2D_dlex30.mat');
    plot(fig1,mu0*H,M,'k>');
    load('OOMMF_HysteresisQuasi3D_dlex30.mat');
    plot(fig1,mu0*H,M,'k<');
    load('OOMMF_Hysteresis3D_dlex30.mat');
    plot(fig1,mu0*H,M,'k^');

    % legend(fig1,'Fortran Mx','Fortran My','Fortran Mz','Matlab Mx','Matlab My','Matlab Mz','NIST \sigma{}(Mx)','NIST <Mx>','NIST \sigma{}(My)','NIST <My>','NIST \sigma{}(Mz)','NIST <Mz>');
    legend(fig1,'"Fortran Explicit method"','"Matlab Implicit method"', 'Matlab Explicit method','OOMMF 2D','OOMMF Quasi3D','OOMMF 3D','Location','SouthEast');
    ylabel(fig1,'<M_i>/M_s')
    xlabel(fig1,'\mu_{0}H_{applied} [T]')
    xlim(fig1,[-0.1 0.1])
    figure(figure1)
end

