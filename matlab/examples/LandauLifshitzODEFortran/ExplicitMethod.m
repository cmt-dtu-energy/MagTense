clearvars
close all

addpath('../../MEX_files');
addpath('../../util');
tic
%test the Fortran implementation of the LL-ODE solver
%get the default problem
%takes the size of the grid as arguments (nx,ny,nz) and a function handle
%that produces the desired field (if not present zero applied field is
%inferred)
problem = DefaultMicroMagProblem(36,9,1);

problem.alpha = -4.42e-6;
problem.gamma = 0;

%initial magnetization
problem.m0(:) = 1/sqrt(3);


%time grid on which to solve the problem
problem = problem.setTime( linspace(0,50,20) );
HystDir = -[1,1,1] ;
problem.solver = getMicroMagSolver( 'Explicit' );

t_H_final = 100;
%time-dependent applied field
HextFct = @(t) (t_H_final-t)' .* HystDir .* (t<t_H_final)'./t_H_final;

problem = problem.setHext( HextFct );

solution = struct();
%convert the class obj to a struct so it can be loaded into fortran
prob_struct = struct(problem);


solution = MagTenseLandauLifshitzSolver_mex( prob_struct, solution );

%find the mean value of the last nt time-evolution points
nt = 10;
M_mean = mean(solution.M(end-nt:end,:,:,:),1);

%find the rms value of the last nt values

rms = sqrt( sum(( M_mean(1,:,:,:)-solution.M(end-nt:end,:,:,:) ).^2,1) );
figure
hold on
col = lines(3);
plot(squeeze(mean(rms(1,:,:,1))),'o','color',col(1,:));
plot(squeeze(mean(rms(1,:,:,2))),'o','color',col(2,:));
plot(squeeze(mean(rms(1,:,:,3))),'o','color',col(3,:));

