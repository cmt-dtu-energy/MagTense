
%%This function compares MagTense to a FEM simulations for a single permanent magnet.
function [] = MagTense_Validation_cylindrical_slice()

%make sure to source the right path for the generic Matlab routines
addpath(genpath('../../util/'));
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
tile.r0 = 5.3984;
tile.theta0 = pi/8;
tile.z0 = 0;

tile.dr = 6.4672-4.3296;
tile.dtheta = pi/8;
tile.dz = 1;

%set the center position of the prism (centered at Origo)
tile.offset = [0, 0, 0];
offset = [5,2,0];

%set the rotation of the prism (centered at Origo)
%the rotation is such that there is first rotation around the z-axis (last
%entry in rotAngles), then around the y-axis (middle entry), and finally
%around the x-axis (first entry in rotAngles). Thus here, the tile is frist
%rotated pi/4 around z, than -pi/3 around y and finally pi/2 around x.
tile.rotAngles = [0,0,0];%[-pi /6, pi /5, pi /2];

%set the easy axis of the cube. This is expected to be with respect to the global
%coordinate system basis
tile.u_ea = [1, 1, 1];
tile.u_ea = tile.u_ea ./ norm(tile.u_ea);

%ensure the two hard axes are perpendicular and normalized
tile.u_oa1 = [1, -1,  1];
tile.u_oa1 = tile.u_oa1 ./ norm(tile.u_oa1);

tile.u_oa2 = [-1,  1, -1 ];
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
% x = -0.2:0.001:1.8;
% y = -0.8:0.001:1.2;
% z = -0.2:0.001:1.8;
% data_FEM_x = load('..\..\..\documentation\examples_FEM_validation\Validation_cylinder\Validation_cylinder_ex2_normH_x.txt');
% data_FEM_y = load('..\..\..\documentation\examples_FEM_validation\Validation_cylinder\Validation_cylinder_ex2_normH_y.txt');
% data_FEM_z = load('..\..\..\documentation\examples_FEM_validation\Validation_cylinder\Validation_cylinder_ex2_normH_z.txt');
data_FEM_x = load('..\..\..\documentation\examples_FEM_validation\Validation_cylinder\Validation_cylinder_example_1_Hx.txt');
data_FEM_y = load('..\..\..\documentation\examples_FEM_validation\Validation_cylinder\Validation_cylinder_example_1_Hy.txt');
data_FEM_z = load('..\..\..\documentation\examples_FEM_validation\Validation_cylinder\Validation_cylinder_example_1_Hz.txt');
x = data_FEM_x(:,1)'+1e-3; %--- Introduce a small offset of 1e-3 to avoid hitting the exact surface
y = data_FEM_y(:,1)'+1e-3;
z = data_FEM_z(:,1)'+1e-3;
% x = data_FEM_x(:,1)'; 
% y = data_FEM_y(:,1)';
% z = data_FEM_z(:,1)';

% pts = zeros( numel(x)+numel(y)+numel(z), 3 );
% pts(1:numel(x),:) = [x; zeros(1,numel(x))+offset(2); zeros(1,numel(x))+offset(3)]';
% pts((numel(x)+1):(numel(x)+numel(y)),:) = [zeros(1,numel(y))+offset(1); y; zeros(1,numel(y))+offset(3)]';
% pts((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z)),:) = [zeros(1,numel(z))+offset(1); zeros(1,numel(z))+offset(2); z]';
pts = [x; y; z]';

tic
%get the field
H = getHFromTiles_mex( tile, pts, int32( length(tile) ), int32( length(pts(:,1)) ) );
% N = getNFromTile_mex( tile, pts, int32( length(pts(:,1)) ) );
toc   

%Find the norm of the field
Hnorm = squeeze( sqrt( sum(H.^2,2) ) );


%Make a figure
figure1 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto');
fig1 = axes('Parent',figure1,'Layer','top','FontSize',16);
hold all
grid on
box on

%The starting point for the line along which the field is calculated
x_start =  2;
y_start = -1;
z_start = -3;

%Plot the solution
plot(sqrt((x-x_start).^2+(y-y_start).^2+(z-z_start).^2),H(:,1),'-','linewidth',2);
plot(sqrt((x-x_start).^2+(y-y_start).^2+(z-z_start).^2),H(:,2),'-','linewidth',2);
plot(sqrt((x-x_start).^2+(y-y_start).^2+(z-z_start).^2),H(:,3),'-','linewidth',2);

%Load comparison data from FEM simulation
FEM_dist = sqrt((data_FEM_x(:,1)-x_start).^2+(data_FEM_y(:,1)-y_start).^2+(data_FEM_z(:,1)-z_start).^2);
plot(FEM_dist,data_FEM_x(:,2),'k--','linewidth',2);
plot(FEM_dist,data_FEM_y(:,2),'k--','linewidth',2);
plot(FEM_dist,data_FEM_z(:,2),'k--','linewidth',2);

h_l = legend('MagTense, H_x','MagTense, H_y','MagTense, H_z','FEM','Location','NorthEast');
set(h_l,'fontsize',10);
ylabel('H_i [A m^{-1}]');
xlabel('Distance from [x,y,z] = (2,-1,-3) to (8,5,3) [m]');
xlim([0 max(FEM_dist)])
end