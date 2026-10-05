%%A soft magnetic sphere in a uniform applied field.
%%
%%The applied field is entered as a tile of type 'Uniformfield' (102), which is not a geometry
%%but a uniform field source: its M vector holds the applied field H_app in A/m, it magnetizes
%%the other tiles in the iteration and it is included in the field returned at the evaluation
%%points.
%%
%%A sphere with constant relative permeability mu_r magnetizes uniformly with
%%    M = 3 (mu_r - 1) / (mu_r + 2) H_app,
%%and outside it the total field is H_app plus that of a point dipole with moment M V, which on
%%the axis along H_app is H_app + 2 M V / (4 pi r^3). Both are checked here.
function [err_M, err_H] = MagTense_Example005_soft_sphere_in_uniform_field()

addpath(genpath('../../../util/'));
addpath('../../../MEX_files/');
mu0 = 4*pi*1e-7;

mur = 20;
radius = 0.01;                          % m
H_app = [0, 0, 0.1/mu0];                % 0.1 T along z, entered as H in A/m

%%The sphere: soft with a constant permeability
tile = DefaultMagTile();
tile = tile.setMagnetType('Soft_const_mur');
tile = tile.setMagTileType('sphere');
tile.abc = [radius, 0, 0];
tile.offset = [0, 0, 0];
tile.u_ea = [0, 0, 1];
tile.u_oa1 = [1, 0, 0];
tile.u_oa2 = [0, 1, 0];
tile.mu_r_ea = mur;
tile.mu_r_oa = mur;
tile.Mrem = 0;

%%The applied field: a tile that is not a geometry
field = DefaultMagTile();
field = field.setMagTileType('Uniformfield');
field.M = H_app;

tiles = [tile, field];

%%Iterate to find the magnetization of the sphere
tiles = IterateMagnetization( tiles, [], [], 1e-10, 500 );
M = tiles(1).M;
M_exact = 3*(mur-1)/(mur+2) * H_app;
err_M = norm(M - M_exact) / norm(M_exact);

%%The total field on the axis outside the sphere
r = [2, 3, 5, 10]' * radius;
pts = [zeros(size(r)), zeros(size(r)), r];
H = getHFromTiles_mex( tiles, pts, int32(length(tiles)), int32(length(pts(:,1))) );
V = 4/3*pi*radius^3;
H_exact = repmat(H_app, length(r), 1) + 2 * repmat(M_exact*V, length(r), 1) ./ (4*pi*r.^3);
err_H = max( vecnorm(H - H_exact, 2, 2) ./ vecnorm(H_exact, 2, 2) );

disp(['mu0 H_app = ' num2str(mu0*H_app) ' T, mu_r = ' num2str(mur)])
disp(['   mu0 M = ' num2str(mu0*M) ' T, exact ' num2str(mu0*M_exact) ' T, relative error ' num2str(err_M, '%.2e')])
disp(['   total field on the axis at r/a = ' num2str(r'/radius) ': max relative error ' num2str(err_H, '%.2e')])

end
