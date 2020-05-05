clear all
%close all
addpath('../../core');
addpath('../../../util');
%setup the problem
%no. of prisms
n = 3;
%size of the prisms
d = 5e-3;%m

tiles(1) = getDefaultMagTile();
tiles(1).tileType = getMagTileType('prism');
tiles(1).abc = [d,d,d];
%boundary conditions (internal) for the left-most tile (it is neighboring
%tile 2, hence n_ind = 2
tiles(1).bdryCdts(1) = struct( 'Type', MagTenseTransientGeometry.FC_INTERNAL, 'l', d/2, 'A', d*d, 'n_ind', 2 );

tiles(2) = tiles(1);
tiles(2).offset = [d,0,0];
tiles(2).bdryCdts(1) = struct( 'Type', MagTenseTransientGeometry.FC_INTERNAL, 'l', d/2, 'A', d*d, 'n_ind', 1 );
tiles(1).bdryCdts(2) = struct( 'Type', MagTenseTransientGeometry.FC_INTERNAL, 'l', d/2, 'A', d*d, 'n_ind', 3 );

tiles(3) = tiles(1);
tiles(3).offset = [2*d,0,0];
tiles(3).bdryCdts(1) = struct( 'Type', MagTenseTransientGeometry.FC_INTERNAL, 'l', d/2, 'A', d*d, 'n_ind', 2 );

%set external condition
tiles(3).bdryCdts(1) = struct( 'Type', MagTenseTransientGeometry.FC_DIRICHLET, 'l', d/2, 'A', d*d, 'n_ind', -1 ); 

debug = false;

%settings
%simple thermal property functions
c = @(T,H,p) ones(size(T)) .* 300;%J/kgK
k = @(T,H,p) ones(size(T)) .* 8;%W/m*K
rho = @(T,H,p) ones(size(T)) .* 7200;%kg/m^3

setts = MagTenseTransientSettings( c, k, rho );
setts.t_tot = 10;%seconds
%fixed timestep
setts.dt = 0.001;%s
%starting time
setts.t = 0;

%the solution object contains the current temperature, field and pressure.
%It also contains the initial condition and may thus be passed directly
%from a file if a simulation is continued
solution = MagTenseTransientSolution();
solution.T = ones(n,1) .* 300;
%solution.T(1)=290;
solution.H = zeros(n,3);
solution.p = zeros(n,1);
%initial thermal properties
solution.c = setts.c(solution.T,solution.H,solution.p);
solution.k = setts.k(solution.T,solution.H,solution.p);
solution.rho = setts.rho(solution.T,solution.H,solution.p);

%setup the geometry. We consider three prisms that have same dimensions and
%are placed adjacent from left to right along the x-axis
geom = MagTenseTransientGeometry( tiles );
%volume of each prism
%geom.dV = ones(n,1) .* (d*d*d);
%thermal resistance (geometrical part) from the center of each prism to the
%center of the face of the currently considered neighboring prism
%geom.R_geom = sparse( n, n );
%prism to the left has only one neighbor. The resistance is the length to
%the center of the mutual face and the total contact area (which, in the
%case of prisms with different dimensions and/or positions may be less than
%the entire surface area of one or both prisms)
%geom.R_geom(1,2) = d/2 / (d*d);
%prism 2 has two neighors
%geom.R_geom(2,1) = d/2 / (d*d);
%geom.R_geom(2,3) = d/2 / (d*d);
%prism 3 has one neighbor
%geom.R_geom(3,2) = d/2 / (d*d);

%run the solution
[solution] = MagTenseTransientSolver( solution, geom, setts, debug );

%do plots

%save data and