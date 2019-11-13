
clear all
close all

thres = [0.1,0.01,0.001,0.0001,0.00001,0.000001,0.0000001,0];
M_out = cell(length(thres),1);
elap_time = zeros(length(thres),1);
addpath('../../Mex_files');
addpath('../util');
for i=1:length(thres)
    tic
    %test the Fortran implementation of the LL-ODE solver
    %get the default problem
    %takes the size of the grid as arguments (nx,ny,nz) and a function handle
    %that produces the desired field (if not present zero applied field is
    %inferred)
    problem = DefaultMicroMagProblem(39,9,3);


    %setup specific problem corresponding to the test provided in
    %Script_3D_TestPhysParamsDynTestAppl01.m
    problem.gamma = 1e20; % precession term
    problem.alpha = -3e19;
    problem.K0 = 5e5; % anisotropy term
    %specify the anisotropy to be along the negative z-direction
    problem.u_ea(:,3) = -1;

    %initial magnetization
    problem.m0(:,:) = 1/sqrt(3);
    %time grid on which to solve the problem
    problem = problem.setTime( linspace(0,2,200) );
    HystDir = [0,1,1] ;
    HystDir = HystDir./norm(HystDir);

    HextFct = @(t) (t>1)' * HystDir;

    problem = problem.setHext( HextFct );
    problem.dem_thres = thres(i);
    solution = struct();

    %call the solver
    %convert the class obj to a struct so it can be loaded into fortran
    prob_struct = struct(problem);
    solution = MagTenseLandauLifshitzSolver_mex( prob_struct, solution );
    elap_time(i) = toc

    %ml_sol = load('SigmaRes.mat');
    %nx=81;ny=21;n=nx*ny;nt=1000;Mx=reshape(ml_sol.SigmaSol(:,1:n),[nt,nx,ny]);My=reshape(ml_sol.SigmaSol(:,n+1:2*n),[nt,nx,ny]);Mz=reshape(ml_sol.SigmaSol(:,2*n+1:3*n),[nt,nx,ny]);
    %close all;
    ind=1;
    %figure;hold on;
%    plot(solution.t,solution.M(:,ind,1));
    M_out{i} = solution.M;
end
getFigure();
i=length(thres)-1;
plot(M_out{i}(:)-M_out{end}(:),'.','displayname',['\delta = ' num2str(thres(i)) ', \Delta t = ' num2str(elap_time(i)) ' s']);

legend('show','location','S');