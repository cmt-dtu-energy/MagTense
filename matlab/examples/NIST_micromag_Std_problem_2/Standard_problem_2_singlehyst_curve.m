clearvars
close all

%--- Use either field 1 or field 2 from the NIST example
resolution = [20,4,1];

figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
mu0 = 4*pi*1e-7;

addpath('../../MEX_files');
addpath('../../util');

tic
%% Run a single hysterhesis curve for standard problem 2
problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));

%--- Solve the first step
solution = struct();
prob_struct = struct(problem);  %convert the class obj to a struct so it can be loaded into fortran

solution = MagTenseLandauLifshitzSolver_mex( prob_struct, solution );


for i = 1:nT
    problem.m0(:) = solution.M(end,:,:);

    %convert the class obj to a struct so it can be loaded into fortran
    prob_struct = struct(problem);

    solution_t = MagTenseLandauLifshitzSolver_mex( prob_struct, solution );
end
% figure; M_end = squeeze(solution.M(end,:,:)); quiver(solution.pts(:,1),solution.pts(:,2),M_end(:,1),M_end(:,2)); axis equal; title('Starting state - Fortran')


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
legend(fig1,'"Matlab Implicit method"', 'Matlab Explicit method','OOMMF 2D','OOMMF Quasi3D','OOMMF 3D','Location','SouthEast');
ylabel(fig1,'<M_i>/M_s')
xlabel(fig1,'\mu_{0}H_{applied} [T]')
xlim(fig1,[-0.1 0.1])
figure(figure1)


