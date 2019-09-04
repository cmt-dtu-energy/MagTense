
%%This function compares MagTense to a FEM simulations for a single permanent magnet.
function [] = MagTense_Validation_cylindrical_slice()
close all
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
tile.tileType = getMagTileType('cylinder');

%set the dimensions of the prism
tile.r0 = 0.3;
tile.theta0 = pi/2;
tile.z0 = 0.5;

tile.dr = 0.3;
tile.dtheta = pi/4;
tile.dz = 0.1;

%set the center position of the prism (centered at Origo)
tile.offset = [0.8, -0.1, 0.3];
% tile.offset = [0,0,0];

%set the rotation of the prism (centered at Origo)
%the rotation is such that there is first rotation around the z-axis (last
%entry in rotAngles), then around the y-axis (middle entry), and finally
%around the x-axis (first entry in rotAngles). Thus here, the tile is frist
%rotated pi/4 around z, than -pi/3 around y and finally pi/2 around x.
tile.rotAngles = [0,0,0];%[-pi /6, pi /5, pi /2];

%set the easy axis of the cube. This is expected to be with respect to the global
%coordinate system basis
tile.u_ea = [0.35355339, 0.35355339, 0.8660254];
tile.u_ea = tile.u_ea ./ norm(tile.u_ea);

%ensure the two hard axes are perpendicular and normalized
tile.u_oa1 = [0.35355339, -0.35355339,  0.        ];
tile.u_oa1 = tile.u_oa1 ./ norm(tile.u_oa1);

tile.u_oa2 = [0.30618622,  0.30618622, -0.25      ];
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
x = -0.5:0.001:1.5;
y = -0.5:0.001:1.5;
z = -1:0.001:1;

pts = zeros( numel(x)+numel(y)+numel(z), 3 );
pts(1:numel(x),:) = [x; zeros(1,numel(y))+tile.offset(2); zeros(1,numel(z))+tile.offset(3)]';
pts((numel(x)+1):(numel(x)+numel(y)),:) = [zeros(1,numel(x))+tile.offset(1); y; zeros(1,numel(z))+tile.offset(3)]';
pts((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z)),:) = [zeros(1,numel(x))+tile.offset(1); zeros(1,numel(y))+tile.offset(2); z]';

delta = 1e-4;
x = [ linspace(tile.offset(1)*(1-delta),tile.offset(1),10) linspace(tile.offset(1),tile.offset(1)*(1+delta),10)];
y = tile.offset(2);
z = tile.offset(3);

pts = zeros(length(x),3);
pts(:,2) = y;
pts(:,3) = z;
pts(:,1) = x;

%get the field
H = getHFromTiles_mex( tile, pts, int32( length(tile) ), int32( length(pts(:,1)) ) );

N = getNFromTile_mex( tile, pts, int32( length(pts(:,1)) ) );         

%Find the norm of the field
Hnorm = squeeze( sqrt( sum(H.^2,2) ) );
getFigure(true);
plot(x,N(:,1,1),'-o');
plot(x,N(:,2,2),'-o');
plot(x,N(:,3,3),'-o');
plot(x,N(:,1,2),'-o');
plot(x,N(:,1,3),'-o');
plot(x,N(:,2,3),'-o');
getFigure(true);
plot(x,H,'-o');
%Make a figure
figure1 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto');
fig1 = axes('Parent',figure1,'Layer','top','FontSize',16);
hold all
grid on
box on

%Plot the solution
plot(x,4*pi*1e-7*Hnorm(1:numel(x)),'r.');
plot(y,4*pi*1e-7*Hnorm((numel(x)+1):(numel(x)+numel(y))),'g.');
plot(z,4*pi*1e-7*Hnorm((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z))),'b.');

% %Load comparison data from FEM simulation
% data_FEM_x = load('..\..\..\documentation\examples_FEM_validation\Validation_prism\Validation_prism_normH_x.txt');
% data_FEM_y = load('..\..\..\documentation\examples_FEM_validation\Validation_prism\Validation_prism_normH_y.txt');
% data_FEM_z = load('..\..\..\documentation\examples_FEM_validation\Validation_prism\Validation_prism_normH_z.txt');
% 
% plot(data_FEM_x(:,1),data_FEM_x(:,2),'ro');
% plot(data_FEM_y(:,1),data_FEM_y(:,2),'go');
% plot(data_FEM_z(:,1),data_FEM_z(:,2),'bo');

h_l = legend('MagTense, x for y,z=offset','MagTense, y for x,z=offset','MagTense, z for x,y=offset','FEM, x for y,z=offset','FEM, y for x,z=offset','FEM, z for x,y=offset','Location','West');
set(h_l,'fontsize',10);
ylabel('|\mu_0{}H| [T]');
xlabel('x, y or z [m]');
end