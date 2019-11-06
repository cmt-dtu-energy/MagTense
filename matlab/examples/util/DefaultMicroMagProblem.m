classdef DefaultMicroMagProblem
    
    %% Defines a class with the default settings for running a MicroMag
    %% problem using the Fortran implementation in MagTense.
    %% Note that the naming convention is dictated in the subroutine 
    %% getProblemFieldNames in the module MagTenseMicroMagIO defined in
    %% the file MagTenseMicroMagIO.f90 in the sub-project MagTenseMicroMag
    %% to the main project MagTense.

properties
    %grid resolution
    grid_n
    %grid_n = [int32(3), int32(3), int32(1)];
    ntot
    %size of the (rectangular) domain in each direction
    grid_L

    %defines the grid type which currently only supports "uniform"
    grid_type


    %easy axes of each cell
    u_ea
    %new or old problem
    ProblemMod
    %solver type ('Explicit', 'Implicit' or 'Dynamic')
    solver

    % Exchange term constant
    A0
    % demag magnetization constant
    Ms
    %Anisotropy constant
    K0

    %precession constant
    gamma

    %
    alpha

    MaxT0

    Hext

    %solution times
    nt
    t

    %initial magnetization (mx = m0(1:n), my = m(n+1:2n), mz = m(2n+1:3n) with
    %n = no. of elements )
    m0
    %x-direction
    theta
    phi    
end
methods
    function obj = DefaultMicroMagProblem(nx,ny,nz)
        nx = int32(nx);
        ny = int32(ny);
        nz = int32(nz);
        %grid resolution
        obj.grid_n = [nx,ny,nz];
        %grid_n = [int32(3), int32(3), int32(1)];
        obj.ntot = prod(obj.grid_n);
        %size of the (rectangular) domain in each direction
        obj.grid_L = [500e-9,125e-9,3e-9];%m

        %defines the grid type which currently only supports "uniform"
        obj.grid_type = getMicroMagGridType('uniform');


        %easy axes of each cell
        obj.u_ea = zeros( obj.ntot, 3 );
        %new or old problem
        obj.ProblemMod = getMicroMagProblemMode( 'new' );
        %solver type ('Explicit', 'Implicit' or 'Dynamic')
        obj.solver = getMicroMagSolver( 'Explicit' );

        % Exchange term constant
        obj.A0 = 1.3e-11;
        % demag magnetization constant
        obj.Ms = 8e5; %A/m
        %Anisotropy constant
        obj.K0 = 0; 

        %precession constant
        obj.gamma = 0;

        %
        obj.alpha = -65104e-17;

        obj.MaxT0 = 2;

        obj.Hext = [1,1,1]./sqrt(3);

        %solution times
        obj.nt = int32(1000);
        obj.t = linspace(0,1,obj.nt);

        %initial magnetization (mx = m0(1:n), my = m(n+1:2n), mz = m(2n+1:3n) with
        %n = no. of elements )
        obj.m0 = zeros(obj.ntot * 3, 1 );
        %x-direction
        obj.theta = pi * rand(obj.ntot,1);
        obj.phi = 2*pi*rand(obj.ntot,1);
        obj.m0(1:obj.ntot) = sin(obj.theta) .* cos(obj.phi);
        obj.m0(obj.ntot+1:2*obj.ntot) = sin(obj.theta) .* sin(obj.phi);
        obj.m0(2*obj.ntot+1:3*obj.ntot) = cos(obj.theta);
    end
end
end