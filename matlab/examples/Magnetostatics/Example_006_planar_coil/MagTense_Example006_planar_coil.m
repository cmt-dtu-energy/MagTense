%%The field of a flat spiral coil, entered as a tile of type 'Planarcoil' (101).
%%
%%A planar coil tile is a flat coil in the xy-plane, centred on the tile's offset, made of 100
%%concentric circular loops evenly spaced from the inner radius a = abc(1) to the outer radius
%%b = abc(2). Its M vector holds the current in each loop, in ampere, and has to be the same in
%%all three entries: M = [I, I, I]. The coil is not rotated and does not take part in the
%%magnetization iteration.
%%
%%The field is checked against the Biot-Savart law summed over the same 100 loops, on the axis,
%%where each loop contributes I R^2 / (2 (R^2 + z^2)^(3/2)), and along a line off the axis,
%%where the loops are integrated numerically. Returns the maximum relative error along each of
%%the two lines. Python counterpart:
%%python/examples/magnetostatics/Example_006_planar_coil/planar_coil.py
function [err_axis, err_off] = MagTense_Example006_planar_coil()

addpath(genpath('../../../util/'));
addpath('../../../MEX_files/');

a = 0.02;                                   % m, inner radius
b = 0.05;                                   % m, outer radius
I = 1;                                      % A, current in each loop
n_loops = 100;                              % fixed in the Fortran implementation

coil = DefaultMagTile();
coil = coil.setMagTileType('Planarcoil');
coil.abc = [a, b, 0];
coil.offset = [0, 0, 0];
coil.M = [I, I, I];
coil.inclIter = int32(0);
coil = struct(coil);

R = a + (b - a) * (0:(n_loops-1)) / (n_loops - 1);

%%On the axis, from the coil plane to four outer radii away
z = linspace(0.001, 4*b, 100)';
pts = [zeros(size(z)), zeros(size(z)), z];
H_axis = getHFromTiles_mex( coil, pts, int32(1), int32(length(z)) );
Hz_exact = sum( I * R.^2 ./ (2 * (R.^2 + z.^2).^(3/2)), 2 );
err_axis = max( abs(H_axis(:,3) - Hz_exact) ./ Hz_exact );

%%Off the axis, along x at a height of one inner radius above the coil. A point in the xz-plane
%%is used because the planar coil is axisymmetric.
x = linspace(0, 2*b, 41)'; x = x(2:end);
pts_off = [x, zeros(size(x)), a * ones(size(x))];
H_off = getHFromTiles_mex( coil, pts_off, int32(1), int32(length(x)) );
H_off_exact = biotSavartLoops(R, I, pts_off);
err_off = max( vecnorm(H_off - H_off_exact, 2, 2) ./ vecnorm(H_off_exact, 2, 2) );

disp(['Planar coil, a = ' num2str(a) ' m, b = ' num2str(b) ' m, ' num2str(n_loops) ' loops of ' num2str(I) ' A'])
disp(['   max. relative error on the axis             ' num2str(err_axis, '%.2e')])
disp(['   max. relative error along x at z = a        ' num2str(err_off, '%.2e')])

figure1 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto');
subplot(1,2,1); hold all; grid on; box on
plot(z, H_axis(:,3), 'r.');
plot(z, Hz_exact, 'k-');
legend('MagTense','Biot-Savart');
xlabel('z [m]'); ylabel('H_z on the axis [A/m]');
subplot(1,2,2); hold all; grid on; box on
plot(x, H_off(:,1), 'r.');
plot(x, H_off(:,3), 'b.');
plot(x, H_off_exact(:,1), 'r-');
plot(x, H_off_exact(:,3), 'b-');
legend('MagTense, H_x','MagTense, H_z','Biot-Savart, H_x','Biot-Savart, H_z');
xlabel(['x [m], at z = ' num2str(a) ' m']); ylabel('H [A/m]');
sgtitle('Planar coil');

end

function H = biotSavartLoops(R, I, pts)
    %H of coaxial circular loops in the xy-plane through the origin, by direct integration
    phi = linspace(0, 2*pi, 4001); phi = phi(1:end-1);
    dphi = phi(2) - phi(1);
    H = zeros(size(pts));
    for k = 1:length(R)
        src = [R(k)*cos(phi); R(k)*sin(phi); zeros(size(phi))]';
        dl = [-R(k)*sin(phi); R(k)*cos(phi); zeros(size(phi))]' * dphi;
        for j = 1:size(pts,1)
            d = pts(j,:) - src;
            H(j,:) = H(j,:) + sum( cross(dl, d, 2) ./ vecnorm(d, 2, 2).^3, 1 );
        end
    end
    H = I * H / (4*pi);
end
