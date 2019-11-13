clear all
close all

addpath('../../Mex_files');
addpath('../util');
tic
%test the Fortran implementation of the LL-ODE solver
%get the default problem
problem = DefaultMicroMagProblem(81,21,1);

Sigma = load('sigma_test1.mat');
problem.m0 = reshape( Sigma.Sigma, [81*21,3] );
problem.setTime( linspace(0,1,1000) );


%problem.alpha = -1e-10;
solution = struct();
%call the solver
prob = struct(problem);
solution = MagTenseLandauLifshitzSolver_mex( prob, solution );
toc

ml_sol = load('SigmaRes.mat');
nx=81;ny=21;n=nx*ny;nt=1000;Mx=reshape(ml_sol.SigmaSol(:,1:n),[nt,nx,ny]);My=reshape(ml_sol.SigmaSol(:,n+1:2*n),[nt,nx,ny]);Mz=reshape(ml_sol.SigmaSol(:,2*n+1:3*n),[nt,nx,ny]);
close all;ind=1;
figure;hold on;plot(solution.t,solution.M(:,ind,2));plot(ml_sol.t,My(:,ind),'--');