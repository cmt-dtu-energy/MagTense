function checks = periodic_exchange_test(options)
%PERIODIC_EXCHANGE_TEST
% Test the implementation of periodic exchange coupling on both grid types of MagTense. This is
% the MATLAB counterpart of python/examples/micromagnetism/periodic_exchange_test.py and uses the
% same geometry, the same material parameters and the same acceptance limits.
%
% MagTense builds the exchange operator in two entirely different ways: from a stencil on a
% uniform grid, and from a mesh analysis on an unstructured mesh of prisms. Both have to honour
% the periodic boundary conditions requested through exchPBC, so both are tested here.
%
% A 6 x 6 x 6 cube with periodic boundary conditions is divided into 8 sections that are
% internally exchange coupled but independent of one another. The test is passed if the
% magnetisation vectors are identical within each connected section but distinct between sections
% due to random initial conditions.
%
% The 8 sections are a 4 x 4 x 4 cube of center tiles, the 8 corners, the xFace (4 tiles on both
% the left and right side), the yFace, the zFace, the xyEdge (8 tiles on the edges perpendicular
% to the xy plane), the xzEdge and the yzEdge. The remaining cells are numerically inert spacers
% and there is no dipole interaction, so the separate sections are effectively decoupled, giving
% 8 tests of the exchange coupling. The center tiles test basic exchange, while the rest check
% the periodic boundary conditions. Only a section that is connected through a periodic boundary
% can become uniform, so this is an end to end test through the actual time integration.
%
% The measure is the largest deviation of a moment from the mean of its section, in units where
% |m| = 1. A section that is not coupled through the boundary keeps its random initial directions
% and lands at a spread of order 1, so the measure separates a working implementation from a
% broken one by orders of magnitude.
%
% Method 2, the exchange operator itself
% --------------------------------------
% A regular lattice of identical cubes is used, so that the assembled exchange matrix can be
% compared with an exact result on both grid types. The matrix is requested as the second output
% of the solver, which returns it as the COO triplets ExchMat_r, ExchMat_c and ExchMat_v. Two
% properties are tested, neither of which depends on the details of the discretisation:
%
%   1. Zero row sums. A uniform magnetisation has to give a vanishing exchange field.
%
%   2. Plane waves are eigenvectors. On a periodic lattice the exchange operator is invariant
%      under a translation by one cell, and a translation invariant operator has the plane waves
%      exp(i*k*x), k = 2*pi*m/L, as exact eigenvectors. The residual |A*phi - lambda*phi| is
%      therefore at round-off level if, and only if, every tile - including the ones at the
%      boundaries - is coupled in the same way. This is the sharp test of the periodic linking.
%
% The eigenvalues are also compared with the Laplacian dispersion relation, lambda(k) -> -k^2 for
% k -> 0, which is reported rather than tested.
%
% Returns a struct array of checks with the fields 'check', 'value', 'limit' and 'passed',
% where a check passes when value < limit. That is the contract used by testMagTenseFunctions.m.

arguments
    options.ShowTheResult {mustBeNumericOrLogical} = true;   % Save the validation figure
    options.use_CUDA {mustBeNumericOrLogical} = false;
end

mu0 = 4*pi*1e-7;

addpath('../../../MEX_files');
addpath('../../../util');

%% Settings

side = 6;                   % The spacer pattern below only works for 6 cells per side
Ntiles = side^3;
a = 1e-8;                   % Side length of an individual micromagnetic cell [m]

% Timestepping
t_end = 1e-8;                            % Total time simulated [s]
t_step = 1e-10;                          % Requested output spacing [s]
nTimesteps = round(t_end/t_step) + 1;    % Include both time endpoints

% Spread of the unit magnetisation vectors within a periodically connected section. A section
% that is not coupled at all keeps its random initial directions and has a spread of order 1, so
% both limits discriminate by more than two orders of magnitude. The unstructured mesh gets the
% looser limit because its stencil reaches every tile sharing a vertex with a face, so the
% section tiles touch the frozen spacers and relax to a spread of about 3e-3 rather than to zero.
section_tol = struct('uniform', 1e-4, 'unstructuredPrisms', 1e-2);

% Exchange operator test: lattice, spacing and thresholds
op_res = [6 6 6];           % Number of cells along x, y and z
op_a = 1e-9;                % Side length of a cell [m]
op_pbc = [1 1 1];           % Periodic exchange boundary conditions along x, y and z
ROW_SUM_TOL = 1e-8;         % Relative to the largest coupling in the matrix
EIGENVECTOR_TOL = 1e-8;     % Relative residual of the plane wave eigenvectors

% Material properties
Ms = 1.0/mu0;               % Saturation magnetisation, 1 T [A/m]
alpha = 1;                  % Gilbert damping [-]
eta = alpha/(1 + alpha^2) * 2.21e5;     % Damping constant [m/(A*s)]
Aex = 1e-11;                % Exchange constant in the coupled sections [J/m]
Aex_spacer = 1e-25;         % Negligible exchange constant in the spacer cells [J/m]

gridTypes = {'uniform', 'unstructuredPrisms'};
gridLabels = {'uniform grid', 'unstructured mesh'};
sectionNames = {'center', 'corner', 'xFace', 'yFace', 'zFace', 'xyEdge', 'xzEdge', 'yzEdge'};

%% Build the sectioned geometry
% Tiles are ordered the way setupGrid in LandauLifshitzEquationSolver.f90 orders them, with x
% running fastest: n = i + 6*(j-1) + 36*(k-1)
[sections, spacers, pts] = build_sections(side, a);

% Ms must remain nonzero because MagTense forms exchange and anisotropy factors by dividing by
% Ms. The tiny spacer exchange decouples the eight sections without introducing 0/0 coefficients.
Aex_n = Aex*ones(Ntiles,1);
Aex_n(spacers) = Aex_spacer;

% Reproducible random initial directions
rng(7);
m0 = 2*rand(Ntiles,3) - 1;
m0 = m0 ./ vecnorm(m0, 2, 2);

%% Run both grid types
deviations = zeros(numel(gridTypes), numel(sectionNames));
operators = cell(numel(gridTypes), 1);
checks = struct('check', {}, 'value', {}, 'limit', {}, 'passed', {});

for g = 1:numel(gridTypes)
    gridType = gridTypes{g};
    fprintf('\n==============================\n%s\n==============================\n', gridLabels{g});

    if strcmp(gridType, 'uniform')
        problem = DefaultMicroMagProblem(side, side, side);
    else
        problem = DefaultMicroMagProblem(Ntiles, 1, 1);
        problem = problem.setMicroMagGridType('unstructuredPrisms');
        problem.grid_pts = pts;
        problem.grid_abc = a*ones(Ntiles,3);
    end
    problem.grid_L = [side*a, side*a, side*a];
    problem = problem.setUseCuda(options.use_CUDA);
    problem = problem.setUseCVODE(false);
    problem = problem.setUseDemag(false);       % No dipole interaction between the sections
    problem = problem.setMicroMagSolver('Dynamic');

    problem.gamma = 2.21e5;
    problem.alpha = eta;
    problem.Ms = Ms*ones(Ntiles,1);
    problem.A0 = Aex_n;
    problem.K0 = zeros(Ntiles,1);
    problem.m0 = m0;

    % The uniform grid has always applied periodic exchange; the unstructured mesh needs to be
    % asked for it explicitly
    problem.exchPBC = int32([1 1 1]);

    HextFct = @(t) (t>=0)' * [0, 0, 0];
    problem = problem.setHext( HextFct, linspace(0, t_end, 2) );
    problem = problem.setTime( linspace(0, t_end, nTimesteps) );

    fprintf('Run simulation on the %s\n', gridLabels{g});
    solution = struct();
    prob_struct = struct(problem);
    solution = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );
    disp('Done running simulation')

    M_npv = squeeze(solution.M(:,:,1,:));
    M_end = squeeze(M_npv(end,:,:));        % (ntot, 3)

    tolerance = section_tol.(gridType);
    for s = 1:numel(sectionNames)
        mask = sections{s};
        Mavg = mean(M_end(mask,:), 1);
        % MagTense returns the magnetisation already normalised to |m| = 1, so the spread is
        % measured in those units and must not be divided by Ms
        maxError = max( sqrt(sum((M_end(mask,:) - Mavg).^2, 2)) );
        deviations(g,s) = maxError;
        if maxError < tolerance
            verdict = 'works';
        else
            verdict = 'FAILED';
        end
        fprintf('  Exchange coupling %s for %s tiles: %.3e\n', verdict, sectionNames{s}, maxError);
        checks(end+1) = struct('check', ...
            sprintf('%s: %s tiles exchange coupled', gridLabels{g}, sectionNames{s}), ...
            'value', maxError, 'limit', tolerance, 'passed', maxError < tolerance); %#ok<AGROW>
    end

    %% Method 2: the assembled exchange operator
    fprintf('\nExchange operator on a %d x %d x %d lattice of %.3e m cubes, exchPBC = [%d %d %d]\n', ...
            op_res(1), op_res(2), op_res(3), op_a, op_pbc(1), op_pbc(2), op_pbc(3));
    A_pbc = exchange_matrix(op_res, op_a, op_pbc, gridType, options.use_CUDA);
    A_free = exchange_matrix(op_res, op_a, [0 0 0], gridType, options.use_CUDA);

    scale = max(abs(A_pbc(:)));
    row_sum = max(abs(sum(A_pbc, 2))) / scale;
    nnzRow = sum(A_pbc ~= 0, 2);
    fprintf('  matrix          : %d x %d, %d nonzeros\n', size(A_pbc,1), size(A_pbc,2), nnz(A_pbc));
    fprintf('  couplings/row   : %d (min) %d (max)\n', min(nnzRow), max(nnzRow));
    fprintf('  max |row sum|   : %.3e [%s]\n', row_sum, passFail(row_sum < ROW_SUM_TOL));
    checks(end+1) = struct('check', ...
        sprintf('%s: exchange operator has vanishing row sums', gridLabels{g}), ...
        'value', row_sum, 'limit', ROW_SUM_TOL, 'passed', row_sum < ROW_SUM_TOL); %#ok<AGROW>

    dirNames = {'x', 'y', 'z'};
    operators{g} = A_pbc; %#ok<AGROW>
    for d = 1:3
        if op_res(d) < 3
            continue
        end
        [lam, residual, kk] = plane_wave_test(A_pbc, op_res, op_a, d, 1);
        % A plane wave is only an eigenvector if the operator is periodic along the direction,
        % so a small residual is a pass for a periodic direction and a failure for a free one
        deviation = abs(real(lam)/(-kk^2) - 1);
        fprintf(['  plane wave %s    : residual %.3e [%s], lambda %+.4e vs -k^2 %+.4e ' ...
                 '(%.1f %% discretisation error)\n'], dirNames{d}, residual, ...
                 passFail(residual < EIGENVECTOR_TOL), real(lam), -kk^2, 100*deviation);
        % The free operator is analysed only as a diagnostic
        [~, residual_free, ~] = plane_wave_test(A_free, op_res, op_a, d, 1);
        fprintf('    (free boundaries residual %.3e, expected to be large)\n', residual_free);
        if op_pbc(d)
            checks(end+1) = struct('check', ...
                sprintf('%s: plane wave along %s is an eigenvector', gridLabels{g}, dirNames{d}), ...
                'value', residual, 'limit', EIGENVECTOR_TOL, ...
                'passed', residual < EIGENVECTOR_TOL); %#ok<AGROW>
        end
    end
end

%% Plot the results
if options.ShowTheResult
    results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
    if ~isfolder(results_dir)
        mkdir(results_dir);
    end

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1300 450]);
    ax = subplot(1, 2, 1, 'Parent', fig);
    b = bar(ax, deviations');
    set(ax, 'YScale', 'log')
    b(1).DisplayName = gridLabels{1};
    b(2).DisplayName = gridLabels{2};
    hold(ax, 'on')
    yline(ax, section_tol.uniform, '--', 'Color', [0 0.39 0], 'LineWidth', 1.2, ...
          'DisplayName', sprintf('Threshold, %s', gridLabels{1}));
    yline(ax, section_tol.unstructuredPrisms, '--', 'Color', [1 0.55 0], 'LineWidth', 1.2, ...
          'DisplayName', sprintf('Threshold, %s', gridLabels{2}));
    % A section that is not coupled through the periodic boundary keeps its random initial
    % directions, which is where the test would land if the periodic linking were missing
    yline(ax, 1.0, ':', 'Color', [0.7 0.13 0.13], 'LineWidth', 1.2, ...
          'DisplayName', 'Uncoupled sections');
    xticks(ax, 1:numel(sectionNames))
    xticklabels(ax, sectionNames)
    xlabel(ax, 'Periodically connected section')
    ylabel(ax, 'Maximum |m_i - <m>|')
    legend(ax, 'Location', 'best')
    title(ax, 'Sectioned simulation')
    box(ax, 'on')

    % ---- Dispersion of the periodic operator, both grids ---------------------------------
    % The two grid types land on top of one another, so the uniform grid is drawn as large open
    % circles with the unstructured mesh marked inside them
    ax2 = subplot(1, 2, 2, 'Parent', fig);
    hold(ax2, 'on')
    markerSpec = {{'o', 11, [0 0 0.5]}, {'x', 7, [0.86 0.08 0.24]}};
    for g = 1:numel(gridTypes)
        for d = 1:3
            if ~op_pbc(d) || op_res(d) < 3
                continue
            end
            modes = 1:floor(op_res(d)/2);
            k_m = zeros(size(modes)); lam_m = zeros(size(modes));
            for mi = 1:numel(modes)
                [lam, ~, kk] = plane_wave_test(operators{g}, op_res, op_a, d, modes(mi));
                k_m(mi) = kk; lam_m(mi) = real(lam);
            end
            spec = markerSpec{g};
            plot(ax2, k_m*op_a, -lam_m*op_a^2, spec{1}, 'MarkerSize', spec{2}, ...
                 'Color', spec{3}, 'LineWidth', 1.4, ...
                 'DisplayName', sprintf('%s, %s', gridLabels{g}, dirNames{d}));
        end
    end
    k_fine = linspace(0, pi, 200);
    plot(ax2, k_fine, k_fine.^2, 'k-', 'DisplayName', 'Continuum, k^2');
    plot(ax2, k_fine, 4*sin(k_fine/2).^2, 'k--', 'DisplayName', '7 point stencil');
    xlabel(ax2, 'k a')
    ylabel(ax2, '-\lambda a^2')
    title(ax2, 'Dispersion of the periodic operator')
    legend(ax2, 'Location', 'best', 'FontSize', 7)
    box(ax2, 'on')

    figure_path = fullfile(results_dir, 'periodic_exchange_test.png');
    exportgraphics(fig, figure_path, 'Resolution', 200);
    close(fig);
    fprintf('\nSaved figure to %s\n', figure_path);
end

if nargout == 0
    if all([checks.passed])
        disp('periodic_exchange_test PASSED')
    else
        disp('periodic_exchange_test FAILED')
    end
    clear checks
end
end


function s = passFail(ok)
if ok
    s = 'pass';
else
    s = 'FAIL';
end
end


function A = exchange_matrix(res, a, pbc, gridType, use_CUDA)
% Run a minimal simulation and return the assembled exchange matrix as a dense array.
%
% Only the setup of the problem is of interest, so the simulation is run for a couple of trivial
% timesteps with the demagnetisation field and the external field switched off. The matrix comes
% back as the second output of the solver, in COO form, exported straight from MKL at exactly
% nnz entries.

ntot = prod(res);
pts = build_prism_mesh(res, a);

if strcmp(gridType, 'uniform')
    problem = DefaultMicroMagProblem(res(1), res(2), res(3));
else
    problem = DefaultMicroMagProblem(ntot, 1, 1);
    problem = problem.setMicroMagGridType('unstructuredPrisms');
    problem.grid_pts = pts;
    problem.grid_abc = a*ones(ntot,3);
end
problem.grid_L = res*a;
problem = problem.setUseCuda(use_CUDA);
problem = problem.setUseCVODE(false);
problem = problem.setUseDemag(false);
problem = problem.setMicroMagSolver('Dynamic');

problem.gamma = 0;
problem.alpha = 4.42e3;
problem.Ms = 8e5*ones(ntot,1);
problem.A0 = 1.3e-11;
problem.K0 = zeros(ntot,1);
problem.m0 = repmat(1/sqrt(3)*[1 1 1], ntot, 1);
problem.exchPBC = int32(pbc);

HextFct = @(t) (t>=0)' * [0, 0, 0];
problem = problem.setHext( HextFct, linspace(0, 1e-12, 2) );
problem = problem.setTime( linspace(0, 1e-12, 2) );

solution = struct();
prob_struct = struct(problem);
[~, GridInfo] = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );

rows = double(GridInfo.ExchMat_r);
cols = double(GridInfo.ExchMat_c);
vals = double(GridInfo.ExchMat_v);
nr = double(GridInfo.ExchMat_nr);
nc = double(GridInfo.ExchMat_nc);
if isempty(rows)
    error('periodic_exchange_test:emptyExchangeMatrix', ...
          'MagTense returned an empty exchange matrix');
end
% The matrix is exported from MKL in COO format, with either zero or one based indexing
if min(rows) >= 1 && min(cols) >= 1
    base = 0;
else
    base = 1;
end
A = full(sparse(rows+base, cols+base, vals, nr, nc));
end


function [lam, residual, k] = plane_wave_test(A, res, a, direction, mode)
% Check that a plane wave along the given direction is an eigenvector of A.
%
% Returns the eigenvalue and the relative residual. Translation invariance along the direction -
% which is what the periodic linking has to produce - is equivalent to a vanishing residual.

nx = res(1); ny = res(2);
ntot = prod(res);
L = res(direction)*a;

% Cell centre coordinate along the direction of the wave, in the same tile ordering as A, which
% setupGrid fills with x running fastest
idx = (0:ntot-1)';
ijk = [mod(idx, nx), mod(floor(idx/nx), ny), floor(idx/(nx*ny))];
x = (ijk(:,direction) + 0.5)*a;

k = 2*pi*mode/L;
phi = exp(1i*k*x);

Aphi = A*phi;
lam = (phi'*Aphi) / (phi'*phi);
residual = norm(Aphi - lam*phi) / max(norm(Aphi), 1e-300);
end


function pts = build_prism_mesh(res, a)
% Centres of a regular lattice of cubes, ordered with x running fastest so that the unstructured
% mesh and the uniform grid put their tiles in the same places

nx = res(1); ny = res(2); nz = res(3);
L = res*a;
pts = zeros(prod(res), 3);
for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            n = i + nx*(j-1) + nx*ny*(k-1);
            pts(n,:) = ([i j k] - 0.5)*a - L/2;
        end
    end
end
end


function [sections, spacers, pts] = build_sections(side, a)
% Boolean masks of the 8 periodically connected sections, the spacer mask, and the cell centres.
% The sectioned geometry happens to be symmetric under exchanging x and z, so getting the tile
% ordering backwards swaps the xFace/zFace and xyEdge/yzEdge labels without changing which eight
% sections are tested.

ntot = side^3;
center = false(ntot,1); corner = false(ntot,1);
xFace = false(ntot,1);  yFace = false(ntot,1);  zFace = false(ntot,1);
xyEdge = false(ntot,1); xzEdge = false(ntot,1); yzEdge = false(ntot,1);
pts = zeros(ntot,3);

n = 0;
for k = 1:side
    for j = 1:side
        for i = 1:side
            n = n + 1;
            pts(n,:) = ([i j k] - 0.5)*a - side*a/2;
            if any(i == [1 side])
                if any(j == [1 side])
                    if any(k == [1 side])
                        corner(n) = true;
                    elseif any(k == [3 4])
                        xyEdge(n) = true;
                    end
                elseif any(j == [3 4])
                    if any(k == [1 side])
                        xzEdge(n) = true;
                    elseif any(k == [3 4])
                        xFace(n) = true;
                    end
                end
            elseif any(i == [3 4])
                if any(j == [1 side])
                    if any(k == [1 side])
                        yzEdge(n) = true;
                    elseif any(k == [3 4])
                        yFace(n) = true;
                    end
                elseif any(j == [3 4])
                    if any(k == [1 side])
                        zFace(n) = true;
                    elseif any(k == [3 4])
                        center(n) = true;
                    end
                end
            end
        end
    end
end

sections = {center, corner, xFace, yFace, zFace, xyEdge, xzEdge, yzEdge};
spacers = ~(center | corner | xFace | yFace | zFace | xyEdge | xzEdge | yzEdge);
end
