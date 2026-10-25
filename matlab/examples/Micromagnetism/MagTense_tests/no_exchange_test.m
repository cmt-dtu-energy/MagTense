function checks = no_exchange_test(options)
%NO_EXCHANGE_TEST
% Test that a micromagnetic problem without exchange, A0 = 0 in every tile, runs on every grid
% type. This is the MATLAB counterpart of
% python/examples/micromagnetism/MagTense_tests/no_exchange_test.py and uses the same problem and
% the same acceptance limits.
%
% With A0 = 0 everywhere the exchange operator is identically zero. MagTense used to normalise A0
% by its largest value before building the operator, which is 0/0 in that case. On the uniform
% grid the result was never read, but on the unstructured meshes the NaN zeroed every
% interpolation weight, the sparse matrices came out empty and the solver segfaulted. MagTense now
% detects the case, uses an explicit zero operator and skips the exchange field altogether.
%
% The same small cube is run on the uniform grid, on an unstructured mesh of prisms and on a
% tetrahedral mesh, with demagnetisation, an applied field and a non-uniform initial state, so
% that the moments do move. For each grid three things are checked:
%
%   1. The exchange field is exactly zero throughout the run.
%   2. The exchange operator handed back by the solver is exactly zero.
%   3. The magnetisation agrees with a run at A0 = 1e-30 J/m in every tile. That run goes through
%      the full construction of the exchange operator, so this checks that the shortcut is the
%      A0 -> 0 limit of the ordinary calculation and not merely something that does not crash.
%
% A fourth check covers the single cell on a uniform grid, the macrospin, whose zero operator is
% built by the same routine: A0 > 0 and A0 = 0 have to give identical results there, as a single
% cell has no neighbour to be exchange coupled to.
%
% A run at a realistic A0 = 1.3e-11 J/m is included in the figure for comparison; it is not
% tested, but it shows that exchange makes a visible difference for this problem, so the
% agreement in check 3 is not a trivial one.
%
% Returns a struct array of checks with the fields 'check', 'value', 'limit' and 'passed',
% where a check passes when value < limit. That is the contract used by testMagTenseFunctions.m.

arguments
    options.ShowTheResult {mustBeNumericOrLogical} = true;   % Save the figure
    options.use_CUDA {mustBeNumericOrLogical} = false;
end

addpath('../../../MEX_files');
addpath('../../../util');

%% Settings

mu0 = 4*pi*1e-7;
p = struct();
p.L = 20e-9;                    % Side length of the cube [m]
p.n_side = 4;                   % Cells along each side on the uniform grid and the prism mesh
p.Ms = 8e5;                     % Saturation magnetisation [A/m]
p.alpha = 4.42e3;               % Damping [m/(A s)]
p.t_end = 2e-9;                 % [s]
p.nt = 21;                      % Number of time steps returned
p.H_applied = 0.1/mu0 * [0.5, 0, sqrt(3)/2];   % 0.1 T at 30 degrees from z in the xz plane
p.use_CUDA = options.use_CUDA;
A_real = 1.3e-11;               % A realistic exchange constant, for the figure only [J/m]
A_tiny = 1e-30;                 % Small enough to make no difference, large enough to build the operator [J/m]

% The exchange field and operator have to vanish exactly, so these are only guards against round-off
field_tol = 1e-12;              % Largest |H_exc| / Ms
operator_tol = 1e-12;           % Largest |entry| of the exchange operator, which is dimensionless
% The A0 = 0 and A0 = A_tiny runs differ by an exchange field of order A_tiny/(mu0 Ms a^2) ~ 1e-14 A/m
limit_tol = 1e-9;               % Largest difference of a moment, in units where |m| = 1
macrospin_tol = 1e-12;          % Same, for the single cell, where the two runs are identical

grids = {'uniform', 'unstructuredPrisms', 'tetrahedron'};
labels = {'uniform grid', 'unstructured mesh', 'tetrahedral mesh'};

%% Test

checks = struct('check', {}, 'value', {}, 'limit', {}, 'passed', {});
curves = cell(numel(grids), 3);
for g = 1:numel(grids)
    fprintf('\n%s\n', labels{g});
    zero = run_case(grids{g}, 0, p);
    tiny = run_case(grids{g}, A_tiny, p);
    realistic = run_case(grids{g}, A_real, p);

    field = max_abs(zero.H_exc) / p.Ms;
    operator = max_abs(zero.exch_v);
    difference = max_abs(zero.M - tiny.M);
    moved = max_abs(zero.M(end,:,:) - zero.M(1,:,:));
    fprintf('  %d cells, largest change of a moment over the run: %.3f\n', zero.ntot, moved);
    fprintf('  max |H_exc| / Ms with A0 = 0:          %.3e (limit %.0e)\n', field, field_tol);
    fprintf('  max |exchange operator| with A0 = 0:   %.3e (limit %.0e)\n', operator, operator_tol);
    fprintf('  max |m(A0 = 0) - m(A0 = %.0e)|:  %.3e (limit %.0e)\n', A_tiny, difference, limit_tol);
    fprintf('  max |m(A0 = 0) - m(A0 = %.1e)|: %.3e (for comparison)\n', A_real, max_abs(zero.M - realistic.M));

    checks(end+1) = make_check([labels{g} ': exchange field is zero with A0 = 0'], field, field_tol); %#ok<AGROW>
    checks(end+1) = make_check([labels{g} ': exchange operator is zero with A0 = 0'], operator, operator_tol); %#ok<AGROW>
    checks(end+1) = make_check([labels{g} ': A0 = 0 matches the limit A0 -> 0'], difference, limit_tol); %#ok<AGROW>
    curves(g,:) = {zero, tiny, realistic};
end

% The macrospin. A single cell has no neighbour, so A0 must not matter at all
fprintf('\nsingle cell\n');
zero = run_case('macrospin', 0, p);
realistic = run_case('macrospin', A_real, p);
difference = max_abs(zero.M - realistic.M);
fprintf('  max |m(A0 = 0) - m(A0 = %.1e)|: %.3e (limit %.0e)\n', A_real, difference, macrospin_tol);
checks(end+1) = make_check('single cell: A0 > 0 and A0 = 0 give the same result', difference, macrospin_tol);

if options.ShowTheResult
    plot_result(curves, labels);
end

if nargout == 0
    if all([checks.passed])
        disp('no_exchange_test PASSED')
    else
        disp('no_exchange_test FAILED')
    end
    clear checks
end
end


function m0 = initial_state(ntot)
% A deterministic, non-uniform initial state, the same in every language.
%
% Cell i (0-based) points along the polar angle 0.2 + 1.2*frac(0.618034*i) and the azimuth 2.4*i,
% so neighbouring cells differ and exchange would have something to act on.

i = (0:ntot-1)';
theta = 0.2 + 1.2 * mod(0.618034 * i, 1);
phi = 2.4 * i;
m0 = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
end


function result = run_case(grid, A0, p)
% Integrate the LL equation for one grid and one value of A0 in every tile.

switch grid
    case 'uniform'
        problem = DefaultMicroMagProblem(p.n_side, p.n_side, p.n_side);
    case 'macrospin'
        problem = DefaultMicroMagProblem(1, 1, 1);
    case 'unstructuredPrisms'
        a = p.L / p.n_side;
        centres = ((0:p.n_side-1) + 0.5) * a - p.L/2;
        [X, Y, Z] = ndgrid(centres, centres, centres);
        pts = [X(:), Y(:), Z(:)];
        problem = DefaultMicroMagProblem(size(pts,1), 1, 1);
        problem = problem.setMicroMagGridType('unstructuredPrisms');
        problem.grid_pts = pts;
        problem.grid_abc = a * ones(size(pts,1), 3);
    case 'tetrahedron'
        % Three cubes along each side, each split into six tetrahedra: 162 tetrahedra
        [nodes, elements] = kuhn_tetra_mesh(3, p.L/3);
        problem = DefaultMicroMagProblem(size(elements,2), 1, 1);
        problem = problem.setMicroMagGridTetrahedron(nodes, elements);
end
ntot = problem.ntot;
problem.grid_L = [p.L, p.L, p.L];
problem = problem.setUseCuda(p.use_CUDA);
problem = problem.setUseCVODE(false);
problem = problem.setMicroMagSolver('Dynamic');

problem.gamma = 0;
problem.alpha = p.alpha;
problem.Ms = p.Ms * ones(ntot,1);
problem.A0 = A0 * ones(ntot,1);
problem.K0 = zeros(ntot,1);
problem.m0 = initial_state(ntot);
problem.ReturnHall = int32(1);

H_applied = p.H_applied;
HextFct = @(t) (t>=0)' * H_applied;
problem = problem.setHext( HextFct, linspace(0, p.t_end, 2) );
problem = problem.setTime( linspace(0, p.t_end, p.nt) );

solution = struct();
prob_struct = struct(problem);
[solution, GridInfo] = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );

nt = size(solution.M, 1);
result = struct();
result.t = solution.t(:);
result.M = reshape(solution.M(:,:,1,:), nt, ntot, 3);      % (time, cell, component)
result.H_exc = solution.H_exc;
result.exch_v = double(GridInfo.ExchMat_v);
result.ntot = ntot;
end


function value = max_abs(x)
% Largest absolute value, NaN if any entry is not finite, so that a NaN fails the check.

x = double(x(:));
if isempty(x)
    value = 0;
elseif ~all(isfinite(x))
    value = NaN;
else
    value = max(abs(x));
end
end


function c = make_check(name, value, limit)
c = struct('check', name, 'value', value, 'limit', limit, 'passed', value < limit);
end


function [nodes, elements] = kuhn_tetra_mesh(n, a)
% A cube of n x n x n cubes of side a, centred on the origin, each split into six tetrahedra by
% the Kuhn subdivision, as create_tetra_mesh in python/src/magtense/utils.py does. Returns the
% nodes as a 3 x M array and the 1-based connectivity as a 4 x N array.

kuhn = [1 2 4 8; 1 2 6 8; 1 3 4 8; 1 3 7 8; 1 5 6 8; 1 5 7 8];
[ii, jj, kk] = ndgrid(0:n, 0:n, 0:n);
nodes = [ii(:) jj(:) kk(:)]' * a - n*a/2;
nid = reshape(1:size(nodes,2), [n+1, n+1, n+1]);

elements = zeros(4, 6*n^3);
e = 0;
for ic = 1:n
    for jc = 1:n
        for kc = 1:n
            corners = zeros(1,8);
            m = 0;
            for k = 0:1
                for j = 0:1
                    for i = 0:1
                        m = m + 1;
                        corners(m) = nid(ic+i, jc+j, kc+k);
                    end
                end
            end
            for t = 1:6
                e = e + 1;
                elements(:,e) = corners(kuhn(t,:))';
            end
        end
    end
end
end


function plot_result(curves, labels)
% <m_i>(t) on each grid for A0 = 0 (lines), A0 = 1e-30 J/m (circles) and A0 = 1.3e-11 J/m (dashed).

results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~isfolder(results_dir)
    mkdir(results_dir);
end
colours = [0.86 0.08 0.24; 0.13 0.55 0.13; 0.27 0.51 0.71];

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 480]);
layout = tiledlayout(fig, 1, size(curves,1), 'TileSpacing', 'compact');
for g = 1:size(curves,1)
    [zero, tiny, realistic] = curves{g,:};
    ax = nexttile(layout);
    hold(ax, 'on')
    t_ns = zero.t * 1e9;
    for c = 1:3
        plot(ax, t_ns, mean(zero.M(:,:,c), 2), '-', 'Color', colours(c,:), 'LineWidth', 1.2);
        plot(ax, t_ns(1:2:end), mean(tiny.M(1:2:end,:,c), 2), 'o', 'Color', colours(c,:));
        plot(ax, t_ns, mean(realistic.M(:,:,c), 2), '--', 'Color', colours(c,:), 'LineWidth', 1.2);
    end
    title(ax, labels{g});
    xlabel(ax, 't [ns]');
    if g == 1
        ylabel(ax, '<m_i>');
    end
    grid(ax, 'on')
    box(ax, 'on')
end

% Colour is the component and the line style the exchange constant, so the legend is split the
% same way and placed below the panels, where it covers no curve
h = gobjects(6,1);
for c = 1:3
    h(c) = plot(ax, NaN, NaN, '-', 'Color', colours(c,:), 'LineWidth', 1.2);
end
h(4) = plot(ax, NaN, NaN, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
h(5) = plot(ax, NaN, NaN, 'o', 'Color', [0.5 0.5 0.5]);
h(6) = plot(ax, NaN, NaN, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
lgd = legend(ax, h, {'<m_x>', '<m_y>', '<m_z>', 'A_0 = 0', 'A_0 = 10^{-30} J/m', 'A_0 = 1.3 10^{-11} J/m'}, ...
             'Orientation', 'horizontal');
lgd.Layout.Tile = 'south';

figure_path = fullfile(results_dir, 'no_exchange_test.png');
exportgraphics(fig, figure_path, 'Resolution', 200);
close(fig);
fprintf('Saved figure to %s\n', figure_path);
end
