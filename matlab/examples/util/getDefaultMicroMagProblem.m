

function problem = getDefaultMicroMagProblem()
%% Returns a struct with the default settings for running a MicroMag
%% problem using the Fortran implementation in MagTense.
%% Note that the naming convention is dictated in the subroutine 
%% getProblemFieldNames in the module MagTenseMicroMagIO defined in
%% the file MagTenseMicroMagIO.f90 in the sub-project MagTenseMicroMag
%% to the main project MagTense.
problem = struct();

%grid resolution
problem.grid_n = [int32(81),int32(21),int32(1)];
%problem.grid_n = [int32(3), int32(3), int32(1)];
ntot = prod(problem.grid_n);
%size of the (rectangular) domain in each direction
problem.grid_L = [500e-9,125e-9,3e-9];%m

%defines the grid type which currently only supports "uniform"
problem.grid_type = getMicroMagGridType('uniform');


%easy axes of each cell
problem.u_ea = zeros( ntot, 3 );
%new or old problem
problem.ProblemMod = getMicroMagProblemMode( 'new' );
%solver type ('Explicit', 'Implicit' or 'Dynamic')
problem.solver = getMicroMagSolver( 'Explicit' );

% Exchange term constant
problem.A0 = 1.3e-11;
% demag magnetization constant
problem.Ms = 8e5; %A/m
%Anisotropy constant
problem.K0 = 0; 

%precession constant
problem.gamma = 0;

%
problem.alpha = -65104e-17;

problem.MaxT0 = 2;

problem.Hext = [1,1,1]./sqrt(3);

%solution times
problem.nt = int32(1000);
problem.t = linspace(0,1,problem.nt);

%initial magnetization (mx = m0(1:n), my = m(n+1:2n), mz = m(2n+1:3n) with
%n = no. of elements )
problem.m0 = zeros(ntot * 3, 1 );
%x-direction
theta = pi * rand(ntot,1);
phi = 2*pi*rand(ntot,1);
problem.m0(1:ntot) = sin(theta) .* cos(phi);
problem.m0(ntot+1:2*ntot) = sin(theta) .* sin(phi);
problem.m0(2*ntot+1:3*ntot) = cos(theta);
end