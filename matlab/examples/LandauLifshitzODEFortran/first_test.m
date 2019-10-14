clear all
close all

addpath('../../Mex_files');
addpath('../util');

%test the Fortran implementation of the LL-ODE solver
%get the default problem
problem = getDefaultMicroMagProblem();
problem.nt = 5;
problem.t = linspace(0,1,5);
solution = struct();
%call the solver
solution = MagTenseLandauLifshitzSolver_mex( problem, solution );