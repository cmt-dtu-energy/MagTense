clear all
close all

tileA=getDefaultMagTile();
tileA.tileType=getMagTileType('prism');
tileA.offset=[0,0,0];
tileA.abc=[1,1,1];

tileB=getDefaultMagTile();
tileB.tileType=getMagTileType('tetrahedron');
tileB.vertices=[[0.5,0.5,0.5];[0.5,0.5,-0.5];[0.5,-0.5,0];[1.0,0.5,0];]';

face = MagTenseGetSharedFaces( tileA, tileB );

figure;
hold on;
plotTiles([tileA,tileB]);
plot3(0.5,sqrt(face.lA^2-0.5^2),0,'ok');
alpha 0.2;