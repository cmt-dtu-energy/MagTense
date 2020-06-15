function N_tensor = ComputeDemagTensorCartesianMesh(X,Y,Z,dx,dy,dz)
%define a default tile
tl = getDefaultMagTile();
%ensure the tiles to be prisms
tl.tileType = getMagTileType('prism');
%set the dimensions of each tile (to be a cubes)
tl.abc = [dx,dy,dz];

tiles = [];
for i=1:numel(X)
    tl.offset = [X(i),Y(i),Z(i)];
    tl.abc = [dx(i),dy(i),dz(i)];
    tiles = [tiles tl];
end
%pts needs to be (n,3)
pts = reshape([tiles.offset],[3,length(tiles)])';
N_tensor = getNTensor_experiment( tiles, pts );