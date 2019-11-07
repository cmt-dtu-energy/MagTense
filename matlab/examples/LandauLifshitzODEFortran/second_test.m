clear all
close all

addpath('../../Mex_files');
addpath('../util');
tic
%test the Fortran implementation of the LL-ODE solver
%get the default problem
%takes the size of the grid as arguments (nx,ny,nz)
problem = DefaultMicroMagProblem(39,9,3);


%setup specific problem corresponding to the test provided in
%Script_3D_TestPhysParamsDynTestAppl01.m
problem.gamma = 1e20; % precession term
problem.alpha = -3e19;
problem.K0 = 5e5; % anisotropy term
%specify the anisotropy to be along the negative z-direction
problem.u_ea(:,3) = -1;

problem.m0(:) = 1/sqrt(3);
problem.t = linspace(0,1,1000);
problem.Hext = [0,0,0];
solution = struct();
%call the solver
problem.nt = int32(length(problem.t));
%convert the class obj to a struct so it can be loaded into fortran
prob_struct = struct(problem);
solution = MagTenseLandauLifshitzSolver_mex( prob_struct, solution );
toc

%ml_sol = load('SigmaRes.mat');
%nx=81;ny=21;n=nx*ny;nt=1000;Mx=reshape(ml_sol.SigmaSol(:,1:n),[nt,nx,ny]);My=reshape(ml_sol.SigmaSol(:,n+1:2*n),[nt,nx,ny]);Mz=reshape(ml_sol.SigmaSol(:,2*n+1:3*n),[nt,nx,ny]);
close all;ind=1;
figure;hold on;plot(solution.t,solution.M(:,ind,1));plot(ml_sol.t,My(:,ind),'--');