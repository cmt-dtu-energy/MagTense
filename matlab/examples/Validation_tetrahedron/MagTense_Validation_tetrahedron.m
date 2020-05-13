
%%This function compares MagTense to a FEM simulations for a single permanent magnet.
function [] = MagTense_Validation_tetrahedron()

%Make a figure
figure1 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto');
fig1 = axes('Parent',figure1,'Layer','top','FontSize',16);
hold all
grid on
box on

%make sure to source the right path for the generic Matlab routines
addpath('../../util/');
addpath('../../MEX_files/');

%define the vacuum permeability
mu0 = 4*pi*1e-7;

%%Get a default tile from MagTense
tile = getDefaultMagTile();
    
%ensure the tile is a permanent magnet
tile.magnetType = getMagnetType('hard');

%set the geometry to be a rectangular prism
tile.tileType = getMagTileType('tetrahedron');

%set the dimensions of the prism
tile.vertices = [[2.5,3,1];[2,1,4];[1.5,4,3];[4.5,5,2]]';

%set the center position of the prism (centered at Origo)
tile.offset = [0,0,0];

%set the rotation of the prism (centered at Origo)
tile.rotAngles = [0, 0, 0];

%set the easy axis of the cube. This is expected to be with respect to the global
%coordinate system basis
tile.M = 1/(4*pi*1e-7)*[0.324264068, 0.734846928, 0.891545179]';

%set the relative permeability for the easy axis
tile.mu_r_ea = 1.00;
%and for the two hard axes
tile.mu_r_oa = 1.00;

%%Now find the field in a set of points
x = -10:0.01:10;
y = x;
z = x;

offset = [3,3,2.5];

pts = zeros( numel(x)+numel(y)+numel(z), 3 );
pts(1:numel(x),:) = [x; zeros(1,numel(x))+offset(2); zeros(1,numel(x))+offset(3)]';
pts((numel(x)+1):(numel(x)+numel(y)),:) = [zeros(1,numel(y))+offset(1); y; zeros(1,numel(y))+offset(3)]';
pts((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z)),:) = [zeros(1,numel(z))+offset(1); zeros(1,numel(z))+offset(2); z]';

%get the field from the Fortran implementation
tic
H = getHFromTiles_mex( tile, pts, int32(length(tile)), int32(length(pts(:,1))) );
toc

%Find the norm of the field
Hnorm = squeeze( sqrt( sum(H.^2,2) ) );

%Plot the solution
plot(x,4*pi*1e-7*Hnorm(1:numel(x)),'r.');
plot(y,4*pi*1e-7*Hnorm((numel(x)+1):(numel(x)+numel(y))),'g.');
plot(z,4*pi*1e-7*Hnorm((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z))),'b.');


%get the field from the Matlab inplementation
tic
H = getHTetrahedron_Matlab( pts', tile.vertices, tile.M )';
toc

%Find the norm of the field
Hnorm = squeeze( sqrt( sum(H.^2,2) ) );

%Plot the solution
plot(x,4*pi*1e-7*Hnorm(1:numel(x)),'rd','Markersize',5);
plot(y,4*pi*1e-7*Hnorm((numel(x)+1):(numel(x)+numel(y))),'gd','Markersize',5);
plot(z,4*pi*1e-7*Hnorm((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z))),'bd','Markersize',5);


%Load comparison data from FEM simulation
data_FEM_x = load('..\..\..\documentation\examples_FEM_validation\Validation_tetrahedron\Validation_tetrahedron_normH_x.txt');
data_FEM_y = load('..\..\..\documentation\examples_FEM_validation\Validation_tetrahedron\Validation_tetrahedron_normH_y.txt');
data_FEM_z = load('..\..\..\documentation\examples_FEM_validation\Validation_tetrahedron\Validation_tetrahedron_normH_z.txt');

plot(data_FEM_x(:,1),data_FEM_x(:,2),'ro');
plot(data_FEM_y(:,1),data_FEM_y(:,2),'go');
plot(data_FEM_z(:,1),data_FEM_z(:,2),'bo');

h_l = legend('MagTense, x for y,z=offset','MagTense, y for x,z=offset','MagTense, z for x,y=offset','MagTense Matlab, x for y,z=offset','MagTense Matlab, y for x,z=offset','MagTense Matlab, z for x,y=offset','FEM, x for y,z=offset','FEM, y for x,z=offset','FEM, z for x,y=offset','Location','NorthWest');
set(h_l,'fontsize',9);
ylabel('|\mu_0{}H| [T]');
xlabel('x, y or z [m]');
ylim([0 1.2]);
end