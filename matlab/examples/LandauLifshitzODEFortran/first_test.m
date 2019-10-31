clear all
close all

addpath('../../Mex_files');
addpath('../util');
tic
%test the Fortran implementation of the LL-ODE solver
%get the default problem
problem = getDefaultMicroMagProblem();

Sigma = load('sigma_test1.mat');
problem.m0 = Sigma.Sigma;
problem.t = linspace(0,1,1000);
problem.Hext = [0,0,0];
solution = struct();
%call the solver
problem.nt = int32(length(problem.t));
solution = MagTenseLandauLifshitzSolver_mex( problem, solution );
toc

ml_sol = load('SigmaRes.mat');