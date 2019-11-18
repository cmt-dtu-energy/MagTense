
%%This function shows how to use MagTense for a single permanent magnet
%%cube.
function [] = MagTense_Example001_SimpleCube()

%make sure to source the right path for the generic Matlab routines
addpath(genpath('../util/'));
addpath('../../MEX_files/');
%define the vacuum permeability
mu0 = 4*pi*1e-7;

%%Get a default tile from MagTense

tile = getDefaultMagTile();
    
%ensure the tile is a permanent magnet
tile.magnetType = getMagnetType('hard');

%set the geometry to be a rectangular prism
tile.tileType = getMagTileType('prism');

%set the dimensions of the prism (5 cm on either side
tile.abc = [0.05,0.05,0.05];

%set the center position of the prism (centered at Origo)
tile.offset = [0,0,0];

%set the easy axis of the cube. This is expected to be with respect to the global
%coordinate system basis
tile.u_ea = [0,1,0];
%ensure the two hard axes are perpendicular and normalized
tile.u_oa1 = [1,0,0];
tile.u_oa2 = [0,0,1];

%set the relative permeability for the easy axis
tile.mu_r_ea = 1.06;
%and for the two hard axes
tile.mu_r_oa = 1.17;

%set the remanence of the magnet (1.2 T converted to A/m)
tile.Mrem = 1.2 / mu0;

%%Let MagTense find the magnetization vector of the cube by iterating to a
%%self-consistent solution. The two empty arguments ( [] ) ensure that
%%default values are used and they are in any case not relevant for this
%%example. The 1e-6 argument is the maximum relative error (change in
%%magnetization from one iteration to the next) allowed before convergence
%%is reached. Finally, 100 reflects the max. no. of allowed iterations
tl2 = tile;
tl2.offset(2) = 0.1;
tl2.u_ea(2) = -1;
tiles = [tile,tl2];
tile = IterateMagnetization( tiles, [], [], 1e-6, 100 );


%%Now find the field in a set of points
%define a range of points spanning the xy plane at z=0
x = linspace( -0.05,0.05, 20);
y = linspace( -0.06,0.16, 51);
z = 0.0251;

%use meshgrid to fill out the span
[X,Y,Z] = meshgrid(x,y,z);


%get the field
[H,Hnorm] = getHMagTense( tile, X, Y, Z );
H = H .* mu0;
Hnorm = Hnorm .* mu0;
%Plot the solution in a contour map
surf_and_con( X,Y,Hnorm);
xlabel('x');
ylabel('y');

%%plot the magnetic field (H). Notice how the field opposes the
%%magnetization inside the cube.
getFigure();
quiver(X,Y,H(:,:,1),H(:,:,2),2,'linewidth',2);
contour(X,Y,Hnorm,'linewidth',2);
plotTiles(tile,true);
alpha 0.3;
axis equal;
xlabel('x');
ylabel('y');


%define a range of points spanning the xy plane at z=0
x = linspace( -0.05,0.05, 20);
y = 0.025 + 5e-3; %height above the cube
z = linspace( -0.06,0.06, 21);

%use meshgrid to fill out the span
[X,Y,Z] = meshgrid(x,y,z);


%get the field
[H,Hnorm] = getHMagTense( tile, X, Y, Z );
H = H .* mu0;
Hnorm = Hnorm .* mu0;

%Plot the solution in a contour map
surf_and_con( squeeze(X),squeeze(Z),Hnorm);
xlabel('x');
ylabel('z');


getFigure();
quiver(squeeze(X),squeeze(Z),H(:,:,1),H(:,:,3),2,'linewidth',2);
contour(squeeze(X),squeeze(Z),Hnorm,'linewidth',2);

alpha 0.3;
axis equal;
xlabel('x');
ylabel('z');

getFigure();contour(squeeze(X),squeeze(Z),H(:,:,1),'linewidth',2,'showtext','on');
xlabel('x');
ylabel('z');
end