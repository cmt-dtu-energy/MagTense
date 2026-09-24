function checks = dipole_field_test(options)
%DIPOLE_FIELD_TEST
% Test the far field of a uniformly magnetised cube against the analytical point dipole. This
% is the MATLAB counterpart of
% python/examples/micromagnetism/MagTense_tests/dipole_field_test.py and uses the same geometry,
% the same evaluation points and the same acceptance limits.
%
% A single cubic prism of side L, magnetised along z, is evaluated at points on spheres of
% increasing radius. Far from the cube the field has to approach that of a point dipole of
% moment m = M*L^3, and the way it approaches it is itself a check: a cube has no quadrupole
% moment by symmetry, so the leading correction is the octupole and the relative deviation
% from the dipole has to fall as R^-4.
%
% That makes this two tests in one. The far field checks the absolute accuracy of the
% magnetostatic kernel, and the exponent checks that the near field has the right multipole
% structure, which a plain threshold at one distance would not catch.
%
% Returns a struct array of checks with the fields 'check', 'value', 'limit' and 'passed',
% where a check passes when value < limit. That is the contract used by testMagTenseFunctions.m.

arguments
    options.ShowTheResult {mustBeNumericOrLogical} = true;   % Save the validation figure
    options.use_CUDA {mustBeNumericOrLogical} = false;       % Not used: no micromagnetics here
end

addpath('../../../MEX_files');
addpath('../../../util');

mu0 = 4*pi*1e-7;
L = 1.2e-6;                             % Side length of the cube [m]
M_rem = 1.2/mu0;                        % Remanent magnetisation, along z [A/m]

distances = [5, 10, 20, 30, 50];        % Distances from the centre, in units of L
n_points = 40;                          % Directions per distance
far_field_distance = 30;
far_field_tol = 1e-3;                   % [%]
exponent_tol = 0.15;

%% Directions: evenly spread over the unit sphere (a Fibonacci sphere), so the test is the same
%% in every language and on every run
k = (0:n_points-1)';
z = 1 - 2*(k + 0.5)/n_points;
phi = k * pi*(3 - sqrt(5));
directions = [sqrt(1 - z.^2).*cos(phi), sqrt(1 - z.^2).*sin(phi), z];

tile = DefaultMagTile();
tile = tile.setMagnetType('hard');
tile = tile.setMagTileType('prism');
tile.abc = [L, L, L];
tile.offset = [0, 0, 0];
tile.M = [0, 0, M_rem];
tile = struct(tile);
m_vec = [0, 0, M_rem] * L^3;

fprintf('Comparing the field of a cube with a point dipole at %d directions per distance\n', n_points);
fprintf('%8s %20s %20s\n', 'R / L', 'max deviation [%]', 'mean deviation [%]');
fprintf('%s\n', repmat('-', 1, 50));

max_deviations = zeros(size(distances));
mean_deviations = zeros(size(distances));
for i = 1:numel(distances)
    pts = directions * (L * distances(i));
    H = getHFromTiles_mex(tile, pts, int32(1), int32(n_points));
    %The whole field vector is compared rather than only its magnitude, so a field of the right
    %strength pointing the wrong way still fails
    r = vecnorm(pts, 2, 2);
    r_hat = pts ./ r;
    H_dipole = (3 * r_hat .* (r_hat * m_vec') - m_vec) ./ (4*pi*r.^3);
    deviations = vecnorm(H - H_dipole, 2, 2) ./ vecnorm(H_dipole, 2, 2) * 100;
    max_deviations(i) = max(deviations);
    mean_deviations(i) = mean(deviations);
    fprintf('%8d %20.3e %20.3e\n', distances(i), max_deviations(i), mean_deviations(i));
end
fprintf('%s\n', repmat('-', 1, 50));

[~, far_field_index] = min(abs(distances - far_field_distance));
far_field_deviation = max_deviations(far_field_index);
fprintf('Deviation at R = %dL: %.3e %% (limit %.3e %%)\n', far_field_distance, far_field_deviation, far_field_tol);

p = polyfit(log(distances), log(max_deviations), 1);
exponent = p(1);
exponent_error = abs(exponent + 4);
fprintf('Fitted falloff exponent: %.3f (expected -4, off by %.3f)\n', exponent, exponent_error);

if options.ShowTheResult
    fig = figure('Visible', 'off');
    loglog(distances, max_deviations, 'o', 'Color', [0.86 0.08 0.24], 'MarkerSize', 8); hold on
    loglog(distances, mean_deviations, 's', 'Color', [0.27 0.51 0.71], 'MarkerSize', 6);
    loglog(distances, max_deviations(1) * (distances / distances(1)).^-4, '--', 'Color', [0.13 0.55 0.13]);
    xlabel('R / L'); ylabel('Deviation from point dipole [%]');
    legend('Calculation, largest of the points', 'Calculation, mean of the points', ...
           'Octupole correction, R^{-4}', 'Location', 'best');
    grid on
    exportgraphics(fig, 'dipole_field_test.png', 'Resolution', 300);
    close(fig);
end

checks = struct('check', {sprintf('far field matches a point dipole at R = %dL', far_field_distance), ...
                          'deviation falls as R^-4, as the octupole correction'}, ...
                'value', {far_field_deviation, exponent_error}, ...
                'limit', {far_field_tol, exponent_tol}, ...
                'passed', {far_field_deviation < far_field_tol, exponent_error < exponent_tol});

if nargout == 0
    if all([checks.passed])
        disp('dipole_field_test PASSED')
    else
        disp('dipole_field_test FAILED')
    end
    clear checks
end
end
