%%This function compares MagTense to a FEM simulation for a single permanent magnet shaped as an
%%inverted circular piece, i.e. the region between a circular arc and the straight chord through
%%its end points. Python counterpart:
%%python/examples/magnetostatics/Validation_field_circpiece_inverted/validation_circpiece_inverted.py
function [rel_int_error] = MagTense_Validation_circpiece_inverted()
%Return the integrated error along x, y and z, to enable check for consistency

addpath(genpath('../../../util/'));
addpath('../../../MEX_files/');
mu0 = 4*pi*1e-7;

tile = DefaultMagTile();
tile = tile.setMagnetType('hard');
tile = tile.setMagTileType('circpieceinv');

tile.r0 = 0.3;
tile.theta0 = pi/0.55;
tile.z0 = 0.6;
tile.dr = 0.15;
tile.dtheta = pi/6;
tile.dz = 0.4;
tile.offset = [0.3, 0.5, 0.1];
tile.rotAngles = [0, 0, 0];

%Easy axis in the global coordinate system, and two hard axes perpendicular to it
tile.u_ea = [0.41562694, 0.41562694, 0.80901699];
tile.u_ea = tile.u_ea / norm(tile.u_ea);
tile.u_oa1 = [1, -1, 0] / sqrt(2);
tile.u_oa2 = cross(tile.u_ea, tile.u_oa1);

tile.mu_r_ea = 1.00;
tile.mu_r_oa = 1.00;
tile.Mrem = 1.2 / mu0;

tile = IterateMagnetization( tile, [], [], 1e-6, 100 );

%The FEM field is given along lines through the centre of the piece
eval_offset = [0.6271937452259475, 0.27251823835641853, 0.7];
fem_dir = '../../../../documentation/examples_FEM_validation/Validation_circpiece_inverted/';
rel_int_error = compareWithFEM(struct(tile), fem_dir, 'Validation_circpiece_inverted_normH_%s.txt', eval_offset, 'Inverted circular piece');

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
