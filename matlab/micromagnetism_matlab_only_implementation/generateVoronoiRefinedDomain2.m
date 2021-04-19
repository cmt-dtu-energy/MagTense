function [vor_map_out] = generateVoronoiRefinedDomain2( n_Tc, base_res, L, Lsize )

sigma_Tc = 1 ;
ndim = numel(base_res) ;
%%
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

if (~exist( 'Lsize', 'var' )) 
   Lsize = ones(ndim,1) ;
end

%%generate n_Tc voronoi points each with an individual value of Tc drawn
%%from a distribution
if (~exist( 'vor_map_in', 'var' )) 
   GenerateMap = 1 ;
else
    if isempty(vor_map_in)
        GenerateMap = 1 ;
    end
end
if GenerateMap
    dTc = random( distr, 0, sigma_Tc, [n_Tc, 1] );
    %Generate the Voronoi points
    x = Lsize(1).*rand([n_Tc,1]) - Lsize(1)/2;
    y = Lsize(2).*rand([n_Tc,1]) - Lsize(2)/2;
    z = Lsize(3).*rand([n_Tc,1]) - Lsize(3)/2;
    if 1
        fine_res = base_res.*(2^(L-1)) ;
        dimf = Lsize./fine_res ;
        if fine_res(1)>1 ; xf = linspace(-Lsize(1)/2+dimf(1)/2,Lsize(1)/2-dimf(1)/2,fine_res(1)) ; else ; xf = 0 ; end
        if fine_res(2)>1 ; yf = linspace(-Lsize(2)/2+dimf(2)/2,Lsize(2)/2-dimf(2)/2,fine_res(2)) ; else ; yf = 0 ; end
        if fine_res(3)>1 ; zf = linspace(-Lsize(3)/2+dimf(3)/2,Lsize(3)/2-dimf(3)/2,fine_res(3)) ; else ; zf = 0 ; end
        
        [Xf,Yf,Zf]  = ndgrid(xf,yf,zf) ;
        [x,y,z,InVoronoiRegion] = DoLloydIteration(Xf,Yf,Zf,x,y,z) ;
    end
    
    
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

nx = base_res(1);
ny = base_res(2);
if ndim == 3
    nz = base_res(3);
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
if nx>1 ; xg = Lsize(1).*linspace(dx/2,nx*dx-dx/2,nx)-Lsize(1)/2; else ; xg = 0 ; end
if ny>1 ; yg = Lsize(2).*linspace(dy/2,ny*dy-dy/2,ny)-Lsize(2)/2; else ; yg = 0 ; end
if nz>1 ; zg = Lsize(3).*linspace(dz/2,nz*dz-dz/2,nz)-Lsize(3)/2; else ; zg = 0 ; end

if nx>1 ; dims(:,1) = xg(2)-xg(1); else ; dims(:,1) = Lsize(1) ; end
if ny>1 ; dims(:,2) = yg(2)-yg(1); else ; dims(:,2) = Lsize(2) ; end
if ndim == 3
    if nz>1 ; dims(:,3) = zg(2)-zg(1); else ; dims(:,3) = Lsize(3) ; end
else
    dims(:,3) = Lsize(3)-0;
end
if ndim == 3
[X,Y,Z] = ndgrid(xg,yg,zg);     
else
[X,Y,Z] = meshgrid(xg,yg,zg); 
end
% INSERT HERE Cartesian Grid Analysis !!! This is the coarse grid !
% [fNormX,fNormY,fNormZ,AreaFaces,Volumes,TheSigns,X,Y,Z,Xf,Yf,Zf,TheTs] = CartesianGridAnalysis(x,y,z,Lx,Ly,Lz,dx,dy,dz) ;

pos(:,1) = X(:);
pos(:,2) = Y(:);
pos(:,3) = Z(:);

%refine the grid
[pos_out,dims_out,children_out,nn] = refineRectGrid_v4( pos, dims, [x y z], ndim, L, criteria, Lsize );





% children_out is never used !

%%
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
% vor_map_out.GridInfo = GridInfo ; 
try
vor_map_out.Xf = Xf ;
vor_map_out.Yf = Yf ;
vor_map_out.Zf = Zf ;
vor_map_out.InVoronoiRegion = InVoronoiRegion ;

end
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