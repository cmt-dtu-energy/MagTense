
%%This function compares MagTense to a FEM simulations for a single permanent magnet.
function [] = MagTense_Validation_cylindrical_slice_reproduce_errors()
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

offset_eval = [0.8, 0.2, 0.8];

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
x = 0.675736065:0.000000005:0.6757361;
y = 0.34999:0.000001:0.35001;
z = 0.7499999:0.00000001:0.7500001;

pts = zeros(length(x),3);
pts(:,1) = x;
pts(:,2) = offset_eval(2);
pts(:,3) = offset_eval(3);
% 
% pts = zeros( numel(x)+numel(y)+numel(z), 3 );
% pts(1:numel(x),:) = [x; zeros(1,numel(y))+offset_eval(2); zeros(1,numel(z))+offset_eval(3)]';
% pts((numel(x)+1):(numel(x)+numel(y)),:) = [zeros(1,numel(x))+offset_eval(1); y; zeros(1,numel(z))+offset_eval(3)]';
% pts((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z)),:) = [zeros(1,numel(x))+offset_eval(1); zeros(1,numel(y))+offset_eval(2); z]';
% 
% x = linspace( offset_eval(1)-1e-8,offset_eval(1)+1e-8,100);
% pts = zeros(length(x),3);
% pts(:,1) = x;
% pts(:,2) = offset_eval(2);
% pts(:,3) = offset_eval(3);

%get the field
H = getHFromTiles_mex( tile, pts, int32( length(tile) ), int32( length(pts(:,1)) ) );

N = getNFromTile_mex( tile, pts, int32( length(pts(:,1)) ) );

%%Debugging in Fortran code tells us that the 3rd point has issues in N11
%%and N33 relative to the other tensors. Not exactly equivalent to what is
%%seen in Matlab so something must also be fishy with the rotation
%%transformation. This error seems to be in the transfer of the N tensor
%%(weird!!!) and not the code as such. The code would only produce an odd
%%point in the third value as is also consistent with the returned H field
%%I think the issue is in theta being very close to zero but not quite (ln
%%134 in TileNComponents.f90)

%Find the norm of the field
Hnorm = squeeze( sqrt( sum(H.^2,2) ) );
getFigure(true);
plot(pts(:,1),Hnorm,'o');
getFigure(true);
plot(pts(:,1),N(:,1,1),'-o');
plot(pts(:,1),N(:,2,1),'-o');
plot(pts(:,1),N(:,3,1),'-o');
plot(pts(:,1),N(:,2,2),'-o');
plot(pts(:,1),N(:,2,3),'-o');
plot(pts(:,1),N(:,3,3),'-o');
%Make a figure
figure1 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto');
fig1 = axes('Parent',figure1,'Layer','top','FontSize',16);
hold all
grid on
box on

%Plot the solution
subplot(1,3,1);
plot(x,4*pi*1e-7*Hnorm(1:numel(x)),'r.');
h_l = legend('MagTense, x for y,z=offset','West');
set(h_l,'fontsize',10);
ylabel('|\mu_0{}H| [T]');
xlabel('x [m]');
subplot(1,3,2);
plot(y,4*pi*1e-7*Hnorm((numel(x)+1):(numel(x)+numel(y))),'g.');
h_l = legend('MagTense, y for x,z=offset','West');
set(h_l,'fontsize',10);
ylabel('|\mu_0{}H| [T]');
xlabel('y [m]');
subplot(1,3,3);
plot(z,4*pi*1e-7*Hnorm((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z))),'b.');
h_l = legend('MagTense, z for x,y=offset','West');
set(h_l,'fontsize',10);
ylabel('|\mu_0{}H| [T]');
xlabel('z [m]');
% %Load comparison data from FEM simulation
% data_FEM_x = load('..\..\..\documentation\examples_FEM_validation\Validation_prism\Validation_prism_normH_x.txt');
% data_FEM_y = load('..\..\..\documentation\examples_FEM_validation\Validation_prism\Validation_prism_normH_y.txt');
% data_FEM_z = load('..\..\..\documentation\examples_FEM_validation\Validation_prism\Validation_prism_normH_z.txt');
% 
% plot(data_FEM_x(:,1),data_FEM_x(:,2),'ro');
% plot(data_FEM_y(:,1),data_FEM_y(:,2),'go');
% plot(data_FEM_z(:,1),data_FEM_z(:,2),'bo');

end