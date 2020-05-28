% Rasmus Bj?rk and Kaspar K. Nielsen
%updated 24 May 2020
%Creates and returns a struct containing af Voronoi map based on n_Tc
%points distributed after a normal distribution given sigma_Tc in ndim
%dimensions (1, 2 or 3) with base resolution base_res and refined L times
%according to the criterion where individual refinement cells should belong
%to the closes Voronoi point
function [vor_map_out] = generateVoronoiRefinedDomain( n_Tc, sigma_Tc, ndim, base_res, L, fileOut, distr, doPlots, criteria, vor_map_in )

%--- Initialization of variables
if ~exist( 'distr','var')
    distr = 'norm';
end

if ~exist( 'doPlots','var')
    doPlots = false;
end

if (~exist('criteria','var'))
    criteria.name = 'voronoi_dist';
end

if isempty(distr)
    distr = 'norm';
end

if isempty(doPlots)
    doPlots = false;
end

if (isempty(criteria))
    criteria.name = 'voronoi_dist';
end

%%generate n_Tc voronoi points each with an individual value of Tc drawn
%%from a distribution
if ~exist( 'vor_map_in', 'var' )
    %dTc = random( distr, 0, sigma_Tc, [n_Tc, 1] );
    dTc = randn(n_Tc,1) * sigma_Tc;
    %Generate the Voronoi points
    x = rand([n_Tc,1]);
    y = rand([n_Tc,1]);
    z = rand([n_Tc,1]);
else
    dTc = vor_map_in.dTc;
    x = vor_map_in.x;
    y = vor_map_in.y;
    z = vor_map_in.z;
end

vor_c = zeros( n_Tc, 3 );
vor_c(:,1) = x;
vor_c(:,2) = y;
vor_c(:,3) = z;

if ndim == 2
    vor_c = vor_c(:,1:2);
end

nx = base_res.nx;
ny = base_res.ny;
if ndim == 3
    nz = base_res.nz;
else
    nz = 1;
end

n = nx * ny * nz;

pos = zeros( n, 3 );
dims = zeros( n, 3 );

dx = 1/nx;
dy = 1/ny;
dz = 1/nz;

% xg = dx/2:dx:(nx*dx-dx/2);
% yg = dy/2:dy:(ny*dy-dy/2);
% zg = dz/2:dz:(nz*dz-dz/2);
xg = linspace(dx/2,nx*dx-dx/2,nx);
yg = linspace(dy/2,ny*dy-dy/2,ny);
zg = linspace(dz/2,nz*dz-dz/2,nz);

dims(:,1) = xg(2)-xg(1);
dims(:,2) = yg(2)-yg(1);
if ndim == 3
    dims(:,3) = zg(2)-zg(1);
else
    dims(:,3) = 1-0;
end
[X,Y,Z] = meshgrid(xg,yg,zg);
pos(:,1) = X(:);
pos(:,2) = Y(:);
pos(:,3) = Z(:);

%refine the grid
[pos_out,dims_out,children_out,nn] = refineRectGrid_v2( pos, dims, [x y z], ndim, L, criteria );

Tc_out = zeros( length(pos_out(:,1)), 1 );
for i=1:length(nn)
    Tc_out(i) = dTc( nn(i) );
end
if doPlots
    %--- Visualize the refinement grid
    plot_grid(pos_out, dims_out, ndim, Tc_out, vor_c, dTc)
    colorbar
end
pos = pos_out;
dims = dims_out;

vor_map_out = struct();
vor_map_out.pos = pos;
vor_map_out.dims = dims;
vor_map_out.Tc_distr = Tc_out;
vor_map_out.dTc = dTc;
vor_map_out.vor_c = vor_c;
vor_map_out.x = x;
vor_map_out.y = y;
vor_map_out.z = z;
vor_map_out.L = L;
vor_map_out.n_Tc = n_Tc;
vor_map_out.sigma_Tc = sigma_Tc;
vor_map_out.ndim = ndim;
vor_map_out.base_res = base_res;
%save( [fileOut '_nx_' num2str(nx) '_ny_' num2str(ny) '_nz_' num2str(nz) '_L_' num2str(L) '.mat'],'Tc_out', 'pos', 'vor_c', 'dims', 'L', 'base_res');

end