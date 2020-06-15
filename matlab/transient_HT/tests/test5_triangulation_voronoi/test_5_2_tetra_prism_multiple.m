clear all
close all

%make half the domain a rectangular grid
tileP = getDefaultMagTile();
tileP.tileType = getMagTileType('prism');
tileP.abc = [0.1,0.1,0.1];

res = struct('nx',2,'ny',2,'nz',2);
tiles = refineTiles(tileP,res);

%make points for a tetrahedral mesh
v = [[0.05,0.05,0.05];[0.05,-0.05,0.05];[0.05,-0.05,-0.05];[0.05,0.05,-0.05];
     [0.15,0.05,0.05];[0.15,-0.05,0.05];[0.15,-0.05,-0.05];[0.15,0.05,-0.05]];

DT = delaunayTriangulation(v); 
nT = length(DT.ConnectivityList(:,1));
tl = getDefaultMagTile();
tl.tileType = getMagTileType('tetrahedron');
col = jet(nT);
for i=1:nT
   tl.vertices = DT.Points(DT.ConnectivityList(i,:),:)';
   tl.color = col(i,:);
   tiles = [tiles tl];
end

tiles = MagtenseApplyBoundaryConditions( tiles );

figure;
hold on;
plotTiles(tiles);
alpha 0.2;
%tetramesh(DT);