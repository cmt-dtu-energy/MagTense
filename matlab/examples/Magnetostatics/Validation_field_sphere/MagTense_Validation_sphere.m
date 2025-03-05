
%%This function compares MagTense to a FEM simulations for a single permanent magnet.
function [] = MagTense_Validation_sphere()

%make sure to source the right path for the generic Matlab routines
addpath(genpath('../../../util/'));
addpath('../../../MEX_files/');
%define the vacuum permeability
mu0 = 4*pi*1e-7;

%%Get a default tile from MagTense
tile = getDefaultMagTile();
    
%ensure the tile is a permanent magnet
tile.magnetType = getMagnetType('hard');

%set the geometry to be a rectangular prism
tile.tileType = getMagTileType('sphere');

%set the dimensions of the prism
tile.abc = [1.6,0,0];

%set the center position of the prism (centered at Origo)
% tile.offset = [0.5,0.4,0.1];
tile.offset = [2,3,4];

%set the easy axis of the cube. This is expected to be with respect to the global
%coordinate system basis
tile.u_ea = [1.1, 0.5, 0.3];
tile.u_ea = tile.u_ea ./ norm(tile.u_ea);

%ensure the two hard axes are perpendicular and normalized
tile.u_oa1 = [1.1, -0.5,  0.3];
tile.u_oa1 = tile.u_oa1 ./ norm(tile.u_oa1);

tile.u_oa2 = [-1.1,  0.5, -0.3 ];
tile.u_oa2 = tile.u_oa2 ./ norm(tile.u_oa2);

%set the relative permeability for the easy axis
tile.mu_r_ea = 1.00;
%and for the two hard axes
tile.mu_r_oa = 1.00;

%set the remanence of the magnet (1.2 T converted to A/m)
tile.Mrem = 1.2 / mu0;

%%Let MagTense find the magnetization vector of the cube by iterating to a
%%self-consistent solution. The two empty arguments ( [] ) ensure that
%%default values are used and they are in any case not relevant for this
%%example. The 1e-6 argument is the maximum relative error (change in
%%magnetization from one iteration to the next) allowed before convergence
%%is reached. Finally, 100 reflects the max. no. of allowed iterations
tile = IterateMagnetization( tile, [], [], 1e-6, 100 );

%%Now find the field in a set of points
x = -10:0.01:10;
y = -10:0.01:10;
z = -10:0.01:10;

pts = zeros( numel(x)+numel(y)+numel(z), 3 );
pts(1:numel(x),:) = [x; zeros(1,numel(y))+tile.offset(2); zeros(1,numel(z))+tile.offset(3)]';
pts((numel(x)+1):(numel(x)+numel(y)),:) = [zeros(1,numel(x))+tile.offset(1); y; zeros(1,numel(z))+tile.offset(3)]';
pts((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z)),:) = [zeros(1,numel(x))+tile.offset(1); zeros(1,numel(y))+tile.offset(2); z]';

%get the field
H = getHFromTiles_mex( tile, pts, int32( length(tile) ), int32( length(pts(:,1)) ) );
         
%Find the norm of the field
Hnorm = squeeze( sqrt( sum(H.^2,2) ) );

%Make a figure
figure1 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold all; grid on; box on
figure2 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig2 = axes('Parent',figure2,'Layer','top','FontSize',16); hold all; grid on; box on
figure3 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig3 = axes('Parent',figure3,'Layer','top','FontSize',16); hold all; grid on; box on

%Plot the solution
plot(fig1,x,H(1:numel(x),1),'r.');
plot(fig1,x,H(1:numel(x),2),'g.');
plot(fig1,x,H(1:numel(x),3),'b.');

plot(fig2,y,H((numel(x)+1):(numel(x)+numel(y)),1),'r.');
plot(fig2,y,H((numel(x)+1):(numel(x)+numel(y)),2),'g.');
plot(fig2,y,H((numel(x)+1):(numel(x)+numel(y)),3),'b.');

plot(fig3,z,H((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z)),1),'r.');
plot(fig3,z,H((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z)),2),'g.');
plot(fig3,z,H((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z)),3),'b.');

%Load comparison data from FEM simulation
data_FEM_Hx_x = load('..\..\..\..\documentation\examples_FEM_validation\Validation_sphere\Validation_sphere_Hx_x.txt');
data_FEM_Hx_y = load('..\..\..\..\documentation\examples_FEM_validation\Validation_sphere\Validation_sphere_Hx_y.txt');
data_FEM_Hx_z = load('..\..\..\..\documentation\examples_FEM_validation\Validation_sphere\Validation_sphere_Hx_z.txt');
data_FEM_Hy_x = load('..\..\..\..\documentation\examples_FEM_validation\Validation_sphere\Validation_sphere_Hy_x.txt');
data_FEM_Hy_y = load('..\..\..\..\documentation\examples_FEM_validation\Validation_sphere\Validation_sphere_Hy_y.txt');
data_FEM_Hy_z = load('..\..\..\..\documentation\examples_FEM_validation\Validation_sphere\Validation_sphere_Hy_z.txt');
data_FEM_Hz_x = load('..\..\..\..\documentation\examples_FEM_validation\Validation_sphere\Validation_sphere_Hz_x.txt');
data_FEM_Hz_y = load('..\..\..\..\documentation\examples_FEM_validation\Validation_sphere\Validation_sphere_Hz_y.txt');
data_FEM_Hz_z = load('..\..\..\..\documentation\examples_FEM_validation\Validation_sphere\Validation_sphere_Hz_z.txt');

plot(fig1,data_FEM_Hx_x(:,1),data_FEM_Hx_x(:,2),'ro')
plot(fig1,data_FEM_Hy_x(:,1),data_FEM_Hy_x(:,2),'go')
plot(fig1,data_FEM_Hz_x(:,1),data_FEM_Hz_x(:,2),'bo')

plot(fig2,data_FEM_Hx_y(:,1),data_FEM_Hx_y(:,2),'ro')
plot(fig2,data_FEM_Hy_y(:,1),data_FEM_Hy_y(:,2),'go')
plot(fig2,data_FEM_Hz_y(:,1),data_FEM_Hz_y(:,2),'bo')

plot(fig3,data_FEM_Hx_z(:,1),data_FEM_Hx_z(:,2),'ro')
plot(fig3,data_FEM_Hy_z(:,1),data_FEM_Hy_z(:,2),'go')
plot(fig3,data_FEM_Hz_z(:,1),data_FEM_Hz_z(:,2),'bo')

h_l = legend(fig1,'MagTense, H_x','MagTense, H_y','MagTense, H_z','FEM, H_x','FEM, H_y','FEM, H_z','Location','NorthWest');
set(h_l,'fontsize',10);
ylabel(fig1,'H [A/m]');
xlabel(fig1,'x [m]');

h_l = legend(fig2,'MagTense, H_x','MagTense, H_y','MagTense, H_z','FEM, H_x','FEM, H_y','FEM, H_z','Location','SouthWest');
set(h_l,'fontsize',10);
ylabel(fig2,'H [A/m]');
xlabel(fig2,'y [m]');

h_l = legend(fig3,'MagTense, H_x','MagTense, H_y','MagTense, H_z','FEM, H_x','FEM, H_y','FEM, H_z','Location','SouthWest');
set(h_l,'fontsize',10);
ylabel(fig3,'H [A/m]');
xlabel(fig3,'z [m]');
end