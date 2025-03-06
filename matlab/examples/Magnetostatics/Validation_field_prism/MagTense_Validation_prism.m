
%%This function compares MagTense to a FEM simulations for a single permanent magnet.
function [int_error] = MagTense_Validation_prism()
%Return the integrated error, to enable check for consistency

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
tile.tileType = getMagTileType('prism');

%set the dimensions of the prism
tile.abc = [0.6,0.1,0.3];

%set the center position of the prism (centered at Origo)
tile.offset = [0.5,0.4,0.1];
% tile.offset = [0,0,0];

%set the rotation of the prism (centered at Origo)
%the rotation is such that there is first rotation around the z-axis (last
%entry in rotAngles), then around the y-axis (middle entry), and finally
%around the x-axis (first entry in rotAngles). Thus here, the tile is frist
%rotated pi/4 around z, than -pi/3 around y and finally pi/2 around x.
tile.rotAngles = [pi/2, -pi/3, pi/4];

%set the easy axis of the cube. This is expected to be with respect to the global
%coordinate system basis
tile.u_ea = [0.35355339, 0.61237244, 0.70710678]/norm([0.35355339, 0.61237244, 0.70710678]);

%ensure the two hard axes are perpendicular and normalized
tile.u_oa1 = [ 0.61237244, -0.35355339,  0.        ]/norm([ 0.61237244, -0.35355339,  0.        ]);
tile.u_oa2 = [ 0.25,       0.4330127, -0.5      ]/norm([ 0.25,       0.4330127, -0.5      ]);

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

%get the field
H = getHFromTiles_mex( tile, pts, int32( length(tile) ), int32( length(pts(:,1)) ) );
         
%Find the norm of the field
Hnorm = squeeze( sqrt( sum(H.^2,2) ) );

%Make a figure
figure1 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto');
fig1 = axes('Parent',figure1,'Layer','top','FontSize',16);
hold all
grid on
box on

%Plot the solution
Hnorm_along_x = 4*pi*1e-7*Hnorm(1:numel(x));
Hnorm_along_y = 4*pi*1e-7*Hnorm((numel(x)+1):(numel(x)+numel(y)));
Hnorm_along_z = 4*pi*1e-7*Hnorm((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z)));
plot(x,Hnorm_along_x,'r.');
plot(y,Hnorm_along_y,'g.');
plot(z,Hnorm_along_z,'b.');

%Load comparison data from FEM simulation
data_FEM_x = load('../../../../documentation/examples_FEM_validation/Validation_prism/Validation_prism_normH_x.txt');
data_FEM_y = load('../../../../documentation/examples_FEM_validation/Validation_prism/Validation_prism_normH_y.txt');
data_FEM_z = load('../../../../documentation/examples_FEM_validation/Validation_prism/Validation_prism_normH_z.txt');

plot(data_FEM_x(:,1),data_FEM_x(:,2),'ro');
plot(data_FEM_y(:,1),data_FEM_y(:,2),'go');
plot(data_FEM_z(:,1),data_FEM_z(:,2),'bo');

h_l = legend('MagTense, x for y,z=offset','MagTense, y for x,z=offset','MagTense, z for x,y=offset','FEM, x for y,z=offset','FEM, y for x,z=offset','FEM, z for x,y=offset','Location','West');
set(h_l,'fontsize',10);
ylabel('|\mu_0{}H| [T]');
xlabel('x, y or z [m]');

% Interpolate the MagTense solution to the NIST-published solutions and
% calculate the difference between the results as an integral.
FEM_interp(:,1) = interp1(data_FEM_x(:,1),data_FEM_x(:,2),x,'linear','extrap');
FEM_interp(:,2) = interp1(data_FEM_y(:,1),data_FEM_y(:,2),y,'linear','extrap');
FEM_interp(:,3) = interp1(data_FEM_z(:,1),data_FEM_z(:,2),z,'linear','extrap');
int_error(1) = trapz(x,abs(FEM_interp(:,1)-Hnorm_along_x));
int_error(2) = trapz(y,abs(FEM_interp(:,2)-Hnorm_along_y));
int_error(3) = trapz(z,abs(FEM_interp(:,3)-Hnorm_along_z));

pointwise_error = [abs(interp1(data_FEM_x(:,1),data_FEM_x(:,2),x,'linear','extrap')'-Hnorm_along_x); abs(interp1(data_FEM_y(:,1),data_FEM_y(:,2),y,'linear','extrap')'-Hnorm_along_y); abs(interp1(data_FEM_z(:,1),data_FEM_z(:,2),z,'linear','extrap')'-Hnorm_along_z)];
disp(['Mean error between MagTense and FEM = ' num2str(mean(pointwise_error)) '+/-' num2str(std(pointwise_error))])
disp(['Integrated error between MagTense and FEM is Mx = ' num2str(int_error(1)) ', My = ' num2str(int_error(2)) ', Mz = ' num2str(int_error(3)) ])

end