clear all
%close all
addpath('../../core');
addpath('../../../util');
%setup the problem

%generate Voronoi map with refinements
nTc = 10;
sigmaTc = 0.5;
ndims = 2;
res = struct('nx',20,'ny',20,'nz',1,'Lx',0.001,'Ly',0.001,'Lz',0.001);
L = 4;
vmap = generateVoronoiRefinedDomain( nTc, sigmaTc, ndims, res, L);

n = length(vmap.pos(:,1));

tl = getDefaultMagTile();
tl.tileType = getMagTileType('prism');
tiles(1) = tl;
tiles(n) = tl;

for i=1:n
    tiles(i) = tl;
    tiles(i).offset = vmap.pos(i,:)*0.001;
    tiles(i).abc = vmap.dims(i,:)*0.001;
end

debug = false;

%settings
%simple thermal property functions
c = @(T,H,p) ones(size(T)) .* 300;%J/kgK
k = @(T,H,p) ones(size(T)) .* 8;%W/m*K
rho = @(T,H,p) ones(size(T)) .* 7200;%kg/m^3

setts = MagTenseTransientSettings( c, k, rho );
setts.t_tot = 0.003;%seconds
%fixed timestep
setts.dt = 0.000001;%s
%starting time
setts.t = 0;

%the solution object contains the current temperature, field and pressure.
%It also contains the initial condition and may thus be passed directly
%from a file if a simulation is continued
solution = MagTenseTransientSolution();
solution.T = ones(n,1) .* 300;
solution.H = zeros(n,3);
solution.p = zeros(n,1);
%initial thermal properties
solution.c = setts.c(solution.T,solution.H,solution.p);
solution.k = setts.k(solution.T,solution.H,solution.p);
solution.rho = setts.rho(solution.T,solution.H,solution.p);


 %find the internal boundary conditions based on the geometric
 %data
 tiles = MagtenseApplyBoundaryConditions( tiles );

%set external condition
%set all x-faces to the left to a fixed temperature (290 K) and all to the
%right to another (310 K)
bdryLeft = @(t) 290;
bdryRight = @(t) 310;
for i=1:length(tiles)
    foundLeft = false;
    foundRight = false;
    for j=1:length(tiles(i).bdryCdts)
        if tiles(i).bdryCdts(j).FaceID == 2 % 2 means x-face to the left
            foundLeft = true;
        elseif tiles(i).bdryCdts(j).FaceID == 1 %means x-face to the right
            foundRight = true;
        end
    end
    if ~foundLeft
        tiles(i).color = [0,1,0];
        tiles(i).bdryCdts(length(tiles(i).bdryCdts)+1) = struct( 'Type', MagTenseTransientGeometry.FC_DIRICHLET, 'l', tiles(i).abc(1)/2, 'A', prod(tiles(i).abc(2:3)), 'n_ind', -1, 'FaceID', 2, 'bdryFun', bdryRight);
    end
    if ~foundRight
        tiles(i).color = [0,0,1];
        tiles(i).bdryCdts(length(tiles(i).bdryCdts)+1) = struct( 'Type', MagTenseTransientGeometry.FC_DIRICHLET, 'l', tiles(i).abc(1)/2, 'A', prod(tiles(i).abc(2:3)), 'n_ind', -1, 'FaceID', 1, 'bdryFun', bdryLeft);
    end
     
end

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
    %caxis([290,310]);
    title(['Time = ' num2str(solution.t_sol(i),'%7.6f')]);
    
    drawnow;
end
%save data and