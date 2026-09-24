%%This function compares MagTense to a FEM simulation for a single permanent magnet shaped as a
%%circular piece. Python counterpart:
%%python/examples/magnetostatics/Validation_field_circpiece/validation_circpiece.py
function [rel_int_error] = MagTense_Validation_circpiece()
%Return the integrated error along x, y and z, to enable check for consistency

addpath(genpath('../../../util/'));
addpath('../../../MEX_files/');
mu0 = 4*pi*1e-7;

tile = DefaultMagTile();
tile = tile.setMagnetType('hard');
tile = tile.setMagTileType('circpiece');

%The piece spans r0 +/- dr/2, theta0 +/- dtheta/2 and z0 +/- dz/2 around the offset
tile.r0 = 0.45;
tile.theta0 = pi/3;
tile.z0 = 0.35;
tile.dr = 0.5;
tile.dtheta = pi/7;
tile.dz = 0.15;
tile.offset = [0.1, 0.3, 0.2];
tile.rotAngles = [0, 0, 0];

%Easy axis in the global coordinate system, and two hard axes perpendicular to it
tile.u_ea = [-0.3095974, -0.22493568, 0.92387953];
tile.u_ea = tile.u_ea / norm(tile.u_ea);
tile.u_oa1 = [0.22493568, -0.3095974, 0];
tile.u_oa1 = tile.u_oa1 / norm(tile.u_oa1);
tile.u_oa2 = cross(tile.u_ea, tile.u_oa1);

tile.mu_r_ea = 1.00;
tile.mu_r_oa = 1.00;
tile.Mrem = 1.2 / mu0;

tile = IterateMagnetization( tile, [], [], 1e-6, 100 );

%The FEM field is given along lines through the centre of the piece
eval_offset = [0.406328622087633, 0.8631363102808783, 0.55];
fem_dir = '../../../../documentation/examples_FEM_validation/Validation_circpiece/';
rel_int_error = compareWithFEM(struct(tile), fem_dir, 'Validation_circpiece_normH_%s.txt', eval_offset, 'Circular piece');

end

function rel_int_error = compareWithFEM(tile, fem_dir, file_pattern, eval_offset, name)
    mu0 = 4*pi*1e-7;
    axes_names = 'xyz';
    colors = 'rgb';

    figure1 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto');
    axes('Parent',figure1,'Layer','top','FontSize',16);
    hold all; grid on; box on

    rel_int_error = zeros(1,3);
    for i = 1:3
        data_FEM = load([fem_dir sprintf(file_pattern, axes_names(i))]);
        pts = repmat(eval_offset, size(data_FEM,1), 1);
        pts(:,i) = data_FEM(:,1);
        H = getHFromTiles_mex( tile, pts, int32(length(tile)), int32(size(pts,1)) );
        Hnorm = mu0 * sqrt(sum(H.^2, 2));

        plot(data_FEM(:,1), Hnorm, [colors(i) '.']);
        plot(data_FEM(:,1), data_FEM(:,2), [colors(i) 'o']);
        rel_int_error(i) = calculate_relative_integral_error(data_FEM(:,1), data_FEM(:,2), data_FEM(:,1), Hnorm);
    end

    legend('MagTense, x','FEM, x','MagTense, y','FEM, y','MagTense, z','FEM, z','Location','West');
    ylabel('|\mu_0{}H| [T]');
    xlabel('x, y or z [m]');
    title([name ' - MagTense vs. FEM']);

    disp(['Relative integrated error between MagTense and FEM is x = ' num2str(rel_int_error(1)) ...
          ', y = ' num2str(rel_int_error(2)) ', z = ' num2str(rel_int_error(3)) ])
end
