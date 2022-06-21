
%%This function shows how to use MagTense for a single permanent magnet
%%cube.
function [] = Test_demag_field()

close all

figure2= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig2 = axes('Parent',figure2,'Layer','top','FontSize',12); hold on; grid on; box on

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

%set the dimensions of the prism (5 cm on either side
tile.abc = [1,1,1];

%set the center position of the prism (centered at Origo)
tile.offset = [0,0,0];

%set the easy axis of the cube. This is expected to be with respect to the global
%coordinate system basis
tile.u_ea = [0,1,0];
%ensure the two hard axes are perpendicular and normalized
tile.u_oa1 = [1,0,0];
tile.u_oa2 = [0,0,1];

%set the relative permeability for the easy axis
tile.mu_r_ea = 1;
%and for the two hard axes
tile.mu_r_oa = 1;

%set the remanence of the magnet (1.2 T converted to A/m)
tile.Mrem = 1.0;

tile.M = [1,1,1];%1.0;

% tile = IterateMagnetization( tile, [], [], 1e-6, 100 );


%%Now find the field in a set of points
%define a range of points spanning the xy plane at z=0
n_pts = 51;

x = 0;
y = linspace( 1, 10000, n_pts);
z = 0;

%use meshgrid to fill out the span
% [X,Y,Z] = meshgrid(x,y,z);
pts = zeros( n_pts, 3 );
pts(:,2) = y;

H = getHFromTiles_mex( tile, pts, int32( length(tile) ), int32( length(pts(:,1)) ) );

%get the field
% [H,Hnorm] = getHMagTense( tile, x, y, z );
H = H .* mu0;

%Find the norm of the field
Hnorm = squeeze( sqrt( sum(H.^2,4) ) );
    
Hnorm = Hnorm .* mu0;
%Plot the solution in a contour map
H = abs(H);
plot(y,H(:,1),'.');
plot(y,H(:,2),'.');
plot(y,H(:,3),'.');
xlabel('y');
ylabel('H_y');

% plot(y,(mu0/(4*pi)*(3*y.*y)./(y.^5)-1./y.^3,'o')
% dipole = mu0*(1/(4*pi)*(2./y.^3));
r = sqrt(x.^2+y.^2+z.^2);
rhat = pts./r';
m = repmat(tile.M,length(pts),1);
% dipole = -mu0*(1/(4*pi)*pts.*dot(m,pts,2)./r'.^5-m./r'.^3);
dipole = mu0*1/(4*pi)*(3*rhat.*dot(m,rhat,2)-m)./r'.^3;
dipole = abs(dipole);
plot(y,dipole(:,1),'o')
plot(y,dipole(:,2),'o')
plot(y,dipole(:,3),'o')

set(gca, 'YScale', 'log');
legend('MagTense','Dipole')

title('Field along y for m=[0,1,0] and box=[1,1,1]')

print('-dpng','MagTense_dipole.png')

figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',12); hold on; grid on; box on
diff_field = (dipole(:,2)-H(:,2))./dipole(:,2);
plot(y,diff_field','.')
set(gca, 'YScale', 'log');
xlabel('y');
ylabel('Error [%]');
print('-dpng','MagTense_dipole_error.png')

end