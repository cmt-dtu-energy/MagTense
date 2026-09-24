%%This function compares the volume-averaged field of a prism with a FEM simulation.
%%
%%The averaged prism (tile type 'Avgprism', 8) returns, at each evaluation point, the field of
%%the prism averaged over a rectangular observation volume centred on that point, which is what
%%a finite-size sensor or a micromagnetic cell sees. The sizes of the observation volumes are
%%the fifth argument of getHFromTiles_mex, one row per point. Python counterpart:
%%python/examples/magnetostatics/Validation_field_avgprism/validation_avgprism.py
function [rel_int_error] = MagTense_Validation_avgprism()
%Return the integrated error of the x and y components, to enable check for consistency

addpath(genpath('../../../util/'));
addpath('../../../MEX_files/');
mu0 = 4*pi*1e-7;

%%A 3 x 5 x 1 prism at the origin, magnetized with 1 A/m along x
tile = DefaultMagTile();
tile = tile.setMagnetType('hard');
tile = tile.setMagTileType('avgprism');
tile.abc = [3, 5, 1];
tile.offset = [0, 0, 0];
tile.u_ea = [1, 0, 0];
tile.Mrem = 1;
%A single hard magnet with mu_r = 1 needs no iteration: its magnetization is its remanence
tile.M = tile.Mrem * tile.u_ea;

%%The FEM reference: the x coordinate of the observation volume and the average of B_x and B_y
%%over it, in tesla
data_FEM = load('../../../../documentation/examples_FEM_validation/Validation_avgprism/Avg_validation_comsol.txt');
x = data_FEM(:,1);

%%Observation volumes of 2 x 4 x 3 centred on a line along x at y = z = 5
pts = [x, 5*ones(size(x)), 5*ones(size(x))];
obs_size = repmat([2, 4, 3], length(x), 1);
H = getHFromTiles_mex( struct(tile), pts, int32(1), int32(length(x)), obs_size );
B = mu0 * H;

rel_int_error(1) = calculate_relative_integral_error(x, data_FEM(:,2), x, B(:,1));
rel_int_error(2) = calculate_relative_integral_error(x, data_FEM(:,3), x, B(:,2));

figure1 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto');
axes('Parent',figure1,'Layer','top','FontSize',16);
hold all; grid on; box on
plot(x, data_FEM(:,2), 'ko');
plot(x, data_FEM(:,3), 'ro');
plot(x, B(:,1), 'k-');
plot(x, B(:,2), 'r-');
legend('<B_x>, FEM','<B_y>, FEM','<B_x>, MagTense','<B_y>, MagTense','Location','best');
xlabel('x [m]');
ylabel('Field averaged over the observation volume [T]');
title('Averaged prism - MagTense vs. FEM');

disp(['Relative integrated error between MagTense and FEM is <B_x> = ' num2str(rel_int_error(1)) ...
      ', <B_y> = ' num2str(rel_int_error(2)) ])

end
