clear all
close all

addpath('../../../util/');
addpath('../../core/');
addpath('../../../mex_files/');
clear all
%close all
addpath('../../core');
addpath('../../../util');
%setup the problem
%no. of prisms
n = 3;
%size of the prisms
d = 5e-3;%m

%vacuum permeability
mu0 = 4*pi*1e-7;

%set external condition
%set all x-faces to the left to a fixed temperature (290 K) and all to the
%right to another (310 K)
bdryLeft = @(t) 290;
bdryRight = @(t) 310;

tiles(1) = getDefaultMagTile();
tiles(1).tileType = getMagTileType('prism');
tiles(1).abc = [d,d,d];
%boundary conditions (internal) for the left-most tile (it is neighboring
%tile 2, hence n_ind = 2
tiles(1).bdryCdts(1) = struct( 'Type', MagTenseTransientGeometry.FC_INTERNAL, 'l', d/2, 'A', d*d, 'n_ind', 2, 'FaceID', 2, 'bdryFun', bdryLeft );

tiles(2) = tiles(1);
tiles(2).offset = [d,0,0];
tiles(2).bdryCdts(1) = struct( 'Type', MagTenseTransientGeometry.FC_INTERNAL, 'l', d/2, 'A', d*d, 'n_ind', 1, 'FaceID', 2, 'bdryFun', @dummyFun );
tiles(2).bdryCdts(2) = struct( 'Type', MagTenseTransientGeometry.FC_INTERNAL, 'l', d/2, 'A', d*d, 'n_ind', 3, 'FaceID', 1, 'bdryFun', @dummyFun );

tiles(3) = tiles(1);
tiles(3).offset = [2*d,0,0];
tiles(3).bdryCdts(1) = struct( 'Type', MagTenseTransientGeometry.FC_INTERNAL, 'l', d/2, 'A', d*d, 'n_ind', 2, 'FaceID', 2, 'bdryFun', @dummyFun );

%set external condition
tiles(3).bdryCdts(2) = struct( 'Type', MagTenseTransientGeometry.FC_DIRICHLET, 'l', d/2, 'A', d*d, 'n_ind', -1, 'FaceID', 1, 'bdryFun', bdryRight ); 

debug = false;

%settings
%simple thermal property functions
c = @(T,H,p) ones(size(T)) .* 300;%J/kgK
k = @(T,H,p) ones(size(T)) .* 8;%W/m*K
rho = @(T,H,p) ones(size(T)) .* 7200;%kg/m^3
%applied field
Happ = @(t) [ones(n,1)'; zeros(n,2)']';

%Magnetization state function
%the argument should be an array of size (n,1) containing the indices into
%the state functions for each tile (in this test case only a single state
%function is assumed for all tiles).
stateFct = MagTenseStateFunction( ones( n, 1 ) );

%Magnetization
M = @( H, T, p, Hyst ) H./sqrt(sum(H.^2,2)) .* stateFct.getMnorm( sqrt(sum(H.^2,2))*mu0, T, p, Hyst );


setts = MagTenseTransientSettings( c, k, rho, M, Happ, stateFct );
setts.t_tot = 100.;%seconds
%fixed timestep
setts.dt = 0.1;%s
%starting time
setts.t = 0;


%the solution object contains the current temperature, field and pressure.
%It also contains the initial condition and may thus be passed directly
%from a file if a simulation is continued
solution = MagTenseTransientSolution( n );
solution.T = ones(n,1) .* 300;
solution.H = (zeros(n,3)+1e-3)/mu0;
solution.p = zeros(n,1)+1e5;
solution.hyst = zeros(n,1);
%initial thermal properties
solution.c = setts.c(solution.T,solution.H,solution.p);
solution.k = setts.k(solution.T,solution.H,solution.p);
solution.rho = setts.rho(solution.T,solution.H,solution.p);
solution.hyst = setts.Hyst( solution );
solution.M = setts.M(solution.H, solution.T, solution.p, solution.hyst );

 %find the internal boundary conditions based on the geometric
 %data
 tiles = MagtenseApplyBoundaryConditions( tiles );
%setup the geometry.
geom = MagTenseTransientGeometry( tiles );

%run the solution
[solution] = MagTenseTransientSolver( solution, geom, setts, debug );

%do plots
%get tile centers
%crds=reshape([tiles.offset],[3,n])';
%x=reshape(crds(:,1),[res.ny,res.nx+1]);
%y=reshape(crds(:,2),[res.ny,res.nx+1]);
close all
[fg,ax] = getFigure();
colorbar;
for i=1:solution.nSnap-1
    cla;
    %T = reshape( solution.T_sol{i},[res.ny,res.nx]);
    T = solution.T_sol{i};
    nC = length(unique(T));
    col = parula(nC);
    colInd = round( (T-min(T))./(max(T)-min(T)) * ( nC - 1 ) + 1 );
    for j=1:n
        tiles(j).color = col(colInd(j),:);
    end
    plotTiles(tiles);
    %surf_and_con(x,y,T,ax);
    caxis([290,310]);
    title(['Time = ' num2str(solution.t_sol(i),'%7.6f')]);
    
    drawnow;
end
