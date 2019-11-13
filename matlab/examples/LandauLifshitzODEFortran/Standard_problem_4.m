clearvars
% close all

figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on

addpath('MEX_files');
addpath('util');
tic
%test the Fortran implementation of the LL-ODE solver
%get the default problem
%takes the size of the grid as arguments (nx,ny,nz) and a function handle
%that produces the desired field (if not present zero applied field is
%inferred)
problem = DefaultMicroMagProblem(36,9,1);

problem.alpha = -4.42e-6;

%initial magnetization
problem.m0(:) = 1/sqrt(3);


%time grid on which to solve the problem
problem = problem.setTime( linspace(0,100,200) );
HystDir = -[1,1,1] ;


%time-dependent applied field
HextFct = @(t) (1-t)' .* HystDir .* (t<1)';

problem = problem.setHext( HextFct );

solution = struct();
%convert the class obj to a struct so it can be loaded into fortran
prob_struct = struct(problem);


solution = MagTenseLandauLifshitzSolver_mex( prob_struct, solution );

%setup problem for the time-dependent solver

problem = DefaultMicroMagProblem(36,9,1);
% problem.gamma = -2.21e5;
% problem.alpha = -4.42e3/problem.Ms;
problem.alpha = -4.42e-6 ;
problem.gamma = -2.21e-4 ;

% problem = problem.setTime( linspace(0,1e15,100) ); % 1e16
problem = problem.setTime( linspace(0,1,200) ); %
HystDir = -[-24.6,4.3,0]/1000 ;
HextFct = @(t) (t>-1)' .*HystDir;
problem = problem.setHext( HextFct );

problem.m0(:) = solution.M(end,:,:);

solution_t = struct();
%convert the class obj to a struct so it can be loaded into fortran
prob_struct = struct(problem);

solution_t = MagTenseLandauLifshitzSolver_mex( prob_struct, solution_t );
toc

figure; sjask = squeeze(solution.M(end,:,:)); quiver(solution.pts(:,1),solution.pts(:,2),sjask(:,1),sjask(:,2)); axis equal; title('Starting state - Fortran')
% figure; sjask = squeeze(solution_t.M(end,:,:)); quiver(solution_t.pts(:,1),solution_t.pts(:,2),sjask(:,1),sjask(:,2)); axis equal;
% figure; hold all; for i=1:3; plot(problem.t,mean(solution_t.M(:,:,i),2),'o'); end;
plot(fig1,problem.t,mean(solution_t.M(:,:,1),2),'rx'); 
plot(fig1,problem.t,mean(solution_t.M(:,:,2),2),'gx'); 
plot(fig1,problem.t,mean(solution_t.M(:,:,3),2),'bx'); 
% figure; hold all; for i=2:4; plot(problem.Hext(:,1),problem.Hext(:,i),'.'); end;

% Run the Matlab version of the micromagnetism code
addpath('..\micromagnetism')
tic
Script_3D_TestDynamicsStdProbl4(fig1);
toc
legend(fig1,'Fortran Mx','Fortran My','Fortran Mz','Matlab Mx','Matlab My','Matlab Mz','NIST Solution Mx','NIST Solution My','NIST Solution Mz');
