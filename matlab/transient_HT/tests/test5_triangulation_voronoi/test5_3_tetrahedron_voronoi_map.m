clear all
close all

addpath('../../../util/');
addpath('../../core/');
%generate Voronoi map with refinements
nTc = 100;
sigmaTc = 0.5;
ndims = 2;
res = struct('nx',20,'ny',20,'nz',1,'Lx',0.001,'Ly',0.001,'Lz',0.001);
L = 4;
vmap = generateVoronoiRefinedDomain( nTc, sigmaTc, ndims, res, L);

%add "bounding box"
xl = min(vmap.vor_c(:,1));
xh = max(vmap.vor_c(:,1));

yl = min(vmap.vor_c(:,2));
yh = max(vmap.vor_c(:,2));

pts = [vmap.vor_c; [xl,yl]; [xl,yh]; [xh,yl]; [xh,yh]];
%add z-dimension
pts = [pts';zeros(1,length(pts(:,1)))]';
u = pts;
u(:,3) = 0.1;
pts = [pts; u];


%make the Delaunay triangulation
DT = delaunayTriangulation( pts );
%get the normals of the free boundary triangles
NM = getTetrahedronNormals( DT );


n = length(DT.ConnectivityList(:,1));

%setup the tetrahedra in terms of tiles
tile = getDefaultMagTile();
tile.tileType = getMagTileType('tetrahedron');

tiles(1) = tile;
tiles(n) = tile;

for i=1:n
   tiles(i) = tile;
   tiles(i).vertices = DT.Points(DT.ConnectivityList(i,:),:)';
end

figure;hold on
tetramesh(DT);
figure;hold on;
plotTiles(tiles);


debug = false;

%settings
%simple thermal property functions
c = @(T,H,p) ones(size(T)) .* 300;%J/kgK
k = @(T,H,p) ones(size(T)) .* 8;%W/m*K
rho = @(T,H,p) ones(size(T)) .* 7200;%kg/m^3

setts = MagTenseTransientSettings( c, k, rho );
setts.t_tot = 1000.;%seconds
%fixed timestep
setts.dt = 0.1;%s
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

inds = [[1,2,3];[1,2,4];[1,3,4];[2,3,4]];
for i=1:length(NM.areas)
   found = false;
   if sqrt(sum((NM.normals(i,:) - [-1,0,0]).^2))==0
      bdryF = bdryLeft;
      found = true;
   elseif sqrt(sum((NM.normals(i,:) - [1,0,0]).^2))==0
      bdryF = bdryRight;
      found = true;
   end
   if found
       tiles(NM.normals_list(i)).color = [1,0,1];
       tiles(NM.normals_list(i)).bdryCdts(length(tiles(NM.normals_list(i)).bdryCdts)+1) = struct( 'Type', MagTenseTransientGeometry.FC_DIRICHLET, 'l', sqrt(sum((mean(tiles(NM.normals_list(i)).vertices,2)'-NM.centers(i,:)).^2)), 'A',NM.areas(i) , 'n_ind', -1, 'FaceID', 1, 'bdryFun', bdryF);
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
