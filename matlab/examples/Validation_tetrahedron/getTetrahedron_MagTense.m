clear all
close all
addpath('../util/');
addpath('../../mex_files/');
%vertices defining the tetrahedron
v = [[-0.5,0,0];[0.5,0,0];[0,0.5,0];[0,0,-0.5]]';

%points at which the field is to be evaluated
n = 10000;
pts = zeros(n,3);

pts(:,1) = linspace(-0.7,0.7,n);
pts(:,2) = 0.1;
pts(:,3) = 0.01;

tl = getDefaultMagTile();

tl.tileType = getMagTileType('tetrahedron');
tl.M = [0,0,1];

tl.vertices = v;

tic
H_mt = getHFromTiles_mex( tl, pts, int32(length(tl)), int32(length(pts(:,1))) );
toc

tic
H_ml = getHTetrahedron_Matlab( pts', v, [0,0,1]' );
toc

getFigure(true);
plot( pts(:,1), sqrt(sum(H_mt.^2,2)), 'displayname','MagTense','linewidth',2);

plot( pts(:,1), sqrt(sum(H_ml.^2,1)), '--','displayname','Matlab','linewidth',2);
