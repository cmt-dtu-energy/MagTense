
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
tile.abc = [1,1,1]*1e-9;

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
n_pts = 1000;
% tile_dia = sqrt(sum(tile.abc.^2));
tile_dia = max(tile.abc);

x = linspace( 1*tile_dia, 1000*tile_dia, n_pts);
y = linspace( 1*tile_dia, 1000*tile_dia, n_pts);
z = linspace( 1*tile_dia, 1000*tile_dia, n_pts);

%use meshgrid to fill out the span
% [X,Y,Z] = meshgrid(x,y,z);
% [x,y,z] = meshgrid(logspace(0,3,10)*max(tile.abc),logspace(0,3,10)*max(tile.abc),logspace(0,3,10)*max(tile.abc));
r_arr = logspace(0,4,13)*max(tile.abc);
x = []; y = []; z = [];
n_s = 30;
for i = 1:length(r_arr)
    [X2,Y2,Z2] = sphere(n_s);
    x = [x X2(:)*r_arr(i)];
    y = [y Y2(:)*r_arr(i)];
    z = [z Z2(:)*r_arr(i)];
end
% pts = zeros( n_pts, 3 );
pts(:,1) = x(:);
pts(:,2) = y(:);
pts(:,3) = z(:);
% pts(:,1) = X(:);
% pts(:,2) = Y(:);
% pts(:,3) = Z(:);
r = sqrt(pts(:,1).^2+pts(:,2).^2+pts(:,3).^2);
r_s = r./tile_dia; 

H = getHFromTiles_mex( tile, pts, int32( length(tile) ), int32( length(pts(:,1)) ) );


%get the field
% [H,Hnorm] = getHMagTense( tile, x, y, z );
H = H .* mu0;

%Find the norm of the field
Hnorm = squeeze( sqrt( sum(H.^2,2) ) );
    
% Hnorm = Hnorm .* mu0;
%Plot the solution in a contour map
% H = abs(H);
plot(r_s,abs(H(:,1)),'rd');
plot(r_s,abs(H(:,2)),'gd');
plot(r_s,abs(H(:,3)),'bd');

%% Single precision N
N = getNFromTile_mex( tile, pts, int32( length(pts(:,1)) ) );
for i = 1:length(N(:,1,1)); 
    Hs(i,:) = mu0*single(squeeze(N(i,:,:)))*tile.M'; 
end
Hsnorm = squeeze( sqrt( sum(Hs.^2,2) ) );
plot(r_s,abs(Hs(:,1)),'r*');
plot(r_s,abs(Hs(:,2)),'g*');
plot(r_s,abs(Hs(:,3)),'b*');

%% Dipole
% plot(y,(mu0/(4*pi)*(3*y.*y)./(y.^5)-1./y.^3,'o')
% dipole = mu0*(1/(4*pi)*(2./y.^3));

% rhat = pts./r';
rhat = pts./r;
m = repmat(tile.M,length(pts),1)*prod(tile.abc);
% dipole = -mu0*(1/(4*pi)*pts.*dot(m,pts,2)./r'.^5-m./r'.^3);
dipole = mu0*1/(4*pi)*(3*rhat.*dot(m,rhat,2)-m)./r.^3;
% dipole = abs(dipole);
dipolenorm = sqrt(sum(dipole.^2,2));

plot(r_s,abs(dipole(:,1)),'ro')
plot(r_s,abs(dipole(:,2)),'go')
plot(r_s,abs(dipole(:,3)),'bo')

xlabel('r [units of largest tile dimension]');
ylabel('|\mu_0H_i|');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlim(gca,[min(r_s)*1/1.5 max(r_s)*1.5])
legend('MagTense','Dipole')

title('Field along r for m=[1,1,1] and box=[1,1,1]')

print('-dpng','MagTense_dipole.png')

figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',12); hold on; grid on; box on
diff_field = abs(dipolenorm-Hnorm)./dipolenorm;
plot(r_s,diff_field','.')
diff_field_s = abs(dipolenorm-Hsnorm)./dipolenorm;
plot(r_s,diff_field_s','d')

% for i = 1:length(r_arr)
%     indx = ((i-1)*(n_s+1)^2+1):(i*(n_s+1)^2);
%     r_s_red = r_s(indx); 
%     errorbar(r_s_red(1),mean(diff_field(indx)'),std(diff_field(indx)'),'dr')
%     errorbar(r_s_red(1),mean(diff_field_s(indx)'),std(diff_field_s(indx)'),'db')
%     plot(r_s_red(1),max(diff_field(indx)),'ro')
%     plot(r_s_red(1),max(diff_field_s(indx)),'bo')
% end

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlim(gca,[min(r_s)*1/1.5 max(r_s)*1.5])
xlabel('r [units of largest tile dimension]');
ylabel('|H_{analytical}-H_{dipole}|/H_{dipole} [-]');
legend('Double precision','Single precision','Location','SouthEast')
print('-dpng','MagTense_dipole_error.png')

end