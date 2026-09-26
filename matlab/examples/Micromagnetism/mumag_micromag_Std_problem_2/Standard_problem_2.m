function [elapsedTime,problem,solution,results] = Standard_problem_2(resolution, d_loop, options)

arguments
    resolution (1,3) {mustBeInteger}                = [5*20,5*4,1];     %--- [nx,ny,nz] of the grid
    d_loop (1,:) {mustBeNumeric}                    = linspace(0.05,0.5,10); %--- The values of d to run the model for
    options.use_CUDA {mustBeNumericOrLogical}       = true          %--- Use CUDA for the calculations
    options.use_CVODE {mustBeNumericOrLogical}      = false;        %--- Use CVODE for the numerical time evolution
    options.ShowTheResult {mustBeNumericOrLogical}  = true          %--- Show the result
    options.use_minimizer {mustBeNumericOrLogical}  = true;         %--- Relax with the energy minimizer instead of integrating the LL equation in time
    options.use_adaptive {mustBeNumericOrLogical}   = true;        %--- Sweep the field with the adaptive stepping instead of the fixed table of 40 fields
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
%--- The field table below is a list of constant fields either way; the choice is how the
%--- equilibrium at each of them is found: by integrating the LL equation over the time window set
%--- with setTime ('ExplicitLL'), or by the Barzilai-Borwein energy minimizer ('Explicit'), which
%--- ignores the time window and stops when the largest torque is below problem.min_tol.
if options.use_minimizer
    problem = problem.setMicroMagSolver( 'Explicit' );
else
    problem = problem.setMicroMagSolver( 'ExplicitLL' );
end

problem = problem.setHext( HextFct, linspace(MaxH,-MaxH,40) );
if options.use_adaptive
    %--- The adaptive sweep walks from H_start to H_end and picks the step itself: it starts at the
    %--- spacing of the fixed table, grows to dH_max where the magnetization hardly changes, and is
    %--- refined to dH_min across the coercive field, where the mean magnetization along the field
    %--- changes sign. The accepted fields come back in solution.H_ext, which needs ReturnHall.
    problem.adaptiveHext = int32(1);
    problem.maxHextSteps = int32(400);
    problem.H_start      =  MaxH*HystDir;
    problem.H_end        = -MaxH*HystDir;
    problem.dH_initial   = 0.005/mu0;
    problem.dH_min       = 0.0005/mu0;
    problem.dH_max       = 0.02/mu0;
    problem.switch_refdH = 0.0005/mu0;
    problem.use_sw_ref   = int32(1);
    %--- The step-control thresholds on the change of the mean magnetization. The defaults suit a
    %--- hard grain that barely moves between switching events; this soft bar changes by a few
    %--- percent per step everywhere, so looser thresholds keep the step at the table spacing in
    %--- the smooth parts and leave the refinement to the switch test.
    problem.dM_min       = 0.01;
    problem.dM_target    = 0.05;
    problem.dM_reject    = 0.15;
    problem.ReturnHall   = int32(1);
end
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

    %--- The cost of the relaxation at each field is measured in effective-field evaluations,
    %--- which is what scales with the problem size, so it is the number to compare between the
    %--- two relaxation methods. The minimizer also reports its iterations and whether it converged.
    results.n_feval(i) = sum(solution.n_feval);
    if options.use_minimizer
        disp(['   Minimizer: ' num2str(sum(solution.n_feval)) ' field evaluations, ' num2str(sum(solution.min_iter)) ...
              ' iterations, ' num2str(sum(solution.min_status == 2)) ' fields not converged'])
    else
        disp(['   LL relaxation: ' num2str(sum(solution.n_feval)) ' field evaluations'])
    end

    %--- The fields that were visited, signed along the sweep direction, in units of Ms
    if options.use_adaptive
        n_fields = double(solution.n_Hext_acc);
        H_acc = squeeze(solution.H_ext(end,1,1:n_fields,:));            % n_fields x 3 [A/m]
        Hn = (H_acc * HystDir') * mu0 / problem.Ms;                     % HystDir has the magnitude 1/mu0
        disp(['   ' num2str(n_fields) ' adaptive field steps'])
    else
        n_fields = problem.nt_Hext;
        Hn = sign(problem.Hext(:,1)).*sqrt(problem.Hext(:,2).^2+problem.Hext(:,3).^2+problem.Hext(:,4).^2)/problem.Ms;
    end
    Mx = zeros(n_fields,1); My = Mx; Mz = Mx; M = Mx;
    for j = 1:n_fields
        Mx_arr = solution.M(end,:,j,1) ;
        My_arr = solution.M(end,:,j,2) ;
        Mz_arr = solution.M(end,:,j,3) ;
        MN = sqrt(Mx_arr.^2+My_arr.^2+Mz_arr.^2) ;
        Mx(j) = mean(Mx_arr./MN) ;
        My(j) = mean(My_arr./MN) ;
        Mz(j) = mean(Mz_arr./MN) ;
        M(j) = Mx(j)*HystDir(1) + My(j)*HystDir(2) + Mz(j)*HystDir(3) ;
    end
    
    results.Mxr(i) = interp1(Hn,Mx,0);
    results.Myr(i) = interp1(Hn,My,0);
    results.Hc(i)  = interp1(M,Hn,0);
    results.n_fields(i) = n_fields;
    results.Hn = Hn;
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

end