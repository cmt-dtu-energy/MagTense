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
% Method 3, the supercell reference
% ---------------------------------
% The plane waves of method 2 are eigenvectors only because the lattice is regular, so that test
% cannot be carried over to an irregular mesh. The same question - is every tile coupled the way
% the periodic linking claims - can be asked without assuming anything about the mesh, by
% replicating the mesh instead of the magnetisation: the exchange operator is assembled on a
% supercell of copies of the mesh with free boundaries, and the rows of the central copy are
% folded back onto the tiles of the original mesh. Every tile of the central copy has exactly the
% neighbourhood that the periodic operator gives the corresponding tile, so the two matrices have
% to agree entry by entry. This is stricter than the plane wave test, which only looks at one
% eigenvector, and it is run on the regular prism mesh as well as on the irregular mesh below.
%
% The irregular mesh
% ------------------
% Both unstructured tests above use a regular lattice of identical cubes, which is the easy case
% for the mesh analysis. They are therefore repeated on an irregular mesh of Voronoi grains of the
% kind used in python/experiments/Grain_perf, where the cells differ in size and a face at one end
% of the domain is linked to several smaller faces at the other end. Method 1 does not carry over,
% since its spacer pattern is built for a 6 x 6 x 6 grid, but its idea does: a slab of cells one
% base cell thick is left out of the mesh, so that what remains falls into two halves which touch
% one another only through the periodic boundary. Leaving the cells out is cleaner than the inert
% spacers of method 1, since there is then nothing in the gap to pin the moments next to it, and
% the halves relax to one common direction if, and only if, the linking works. The same simulation
% is run without exchPBC as a control, which is what shows that the gap really does cut the mesh
% in two.
%
% Returns a struct array of checks with the fields 'check', 'value', 'limit' and 'passed',
% where a check passes when value < limit. That is the contract used by testMagTenseFunctions.m.

arguments
    options.ShowTheResult {mustBeNumericOrLogical} = true;   % Save the validation figure
    options.use_CUDA {mustBeNumericOrLogical} = false;
    options.Grains {mustBeNumericOrLogical} = true;         % Also test the irregular grain mesh
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
SUPERCELL_TOL = 1e-10;      % Deviation from the supercell operator, relative to the largest coupling

% The irregular mesh. A base grid that is refined once wherever a grain boundary passes through a
% cell, which is the mesh type used in python/experiments/Grain_perf. Two sizes are needed: the
% sectioned simulation wants a few cells on either side of the gap, while the supercell reference
% replicates the mesh 27 times and has to stay small enough to assemble quickly.
GRAIN_LABEL = 'grain mesh';
GRAIN_COUNT = 6;            % Number of grains
GRAIN_SEED = 11;            % Seed of the grain positions, so that the mesh is reproducible
GRAIN_REFINEMENTS = 1;      % Number of times a cell at a grain boundary is split into eight
GRAIN_OFFSET = 0.10;        % Half thickness of the refined layer at a boundary, in base cells
GRAIN_RES_SECTION = 5;      % Base cells per side, sectioned simulation (odd, see grain_gap_mask)
GRAIN_RES_OPERATOR = 3;     % Base cells per side, exchange operator test

% The two halves of the split grain mesh, in units where |m| = 1: the distance between their mean
% magnetisations, and the spread of the moments over both of them. Both are far below 1e-2 when
% the halves are linked through the periodic boundary and of order 1 when they are not. Unlike the
% sections above there is no floor from neighbouring spacers, because the gap holds no cells.
SPLIT_TOL = 1e-2;
SPLIT_CONTROL_MIN = 1e-1;

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
supercells = nan(numel(gridTypes), 1);      % Method 3, filled in for the unstructured mesh only
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
    operators{g} = A_pbc;
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

    %% Method 3: the supercell reference
    % The supercell is built from the tiles of a mesh, so it only applies to the unstructured
    % grid type. The uniform grid has no mesh to replicate
    if strcmp(gridType, 'unstructuredPrisms') && any(op_pbc)
        [supercells(g), nCopies] = supercell_deviation(build_prism_mesh(op_res, op_a), ...
            op_a*ones(prod(op_res),3), op_res*op_a, op_pbc, options.use_CUDA);
        fprintf('  supercell       : %.3e away from the operator on %d copies of the mesh [%s]\n', ...
                supercells(g), nCopies, passFail(supercells(g) < SUPERCELL_TOL));
        checks(end+1) = struct('check', ...
            sprintf('%s: periodic operator matches the supercell reference', gridLabels{g}), ...
            'value', supercells(g), 'limit', SUPERCELL_TOL, ...
            'passed', supercells(g) < SUPERCELL_TOL); %#ok<AGROW>
    end
end

%% The irregular mesh of grains
% Method 1 with a gap through the mesh instead of the sections, and methods 2 and 3 on a smaller
% mesh of the same kind
grains = struct();
if options.Grains
    fprintf('\n==============================\n%s\n==============================\n', GRAIN_LABEL);

    [gpts_all, gabc_all, gL, grainIdx] = build_grain_mesh(GRAIN_RES_SECTION, a, GRAIN_COUNT, ...
                                                          GRAIN_REFINEMENTS, GRAIN_OFFSET, GRAIN_SEED);
    gap = grain_gap_mask(gpts_all, a);
    nSizes = numel(unique(round(gabc_all(:,1)/min(gabc_all(:,1)), 9)));
    fprintf('Sectioned mesh: %d cells of %d different sizes, %d of them left out to split the mesh in two\n', ...
            size(gpts_all,1), nSizes, sum(gap));
    % The bounding box, and with it the period, is set by the cells that are kept, and those still
    % reach both ends of the domain along all three directions
    gpts = gpts_all(~gap,:);
    gabc = gabc_all(~gap,:);

    % Reproducible random initial directions, the same ones for both runs
    rng(7);
    m0_grain = 2*rand(size(gpts,1),3) - 1;
    m0_grain = m0_grain ./ vecnorm(m0_grain, 2, 2);
    params = struct('Ms', Ms, 'eta', eta, 'Aex', Aex, 'm0', m0_grain, 't_end', t_end, ...
                    'nTimesteps', nTimesteps);

    % The gap is normal to x, so only the linking along x can connect the two halves. The other
    % two directions are covered by the supercell reference below, which tests all of them at once
    fprintf('Run simulation on the %s, exchPBC = [1 0 0]\n', GRAIN_LABEL);
    [deviation, spread, M_pbc] = split_deviations(gpts, gabc, gL, [1 0 0], options.use_CUDA, params);
    fprintf('Run simulation on the %s, exchPBC = [0 0 0]\n', GRAIN_LABEL);
    [deviation_free, ~, M_free] = split_deviations(gpts, gabc, gL, [0 0 0], options.use_CUDA, params);
    disp('Done running simulation')

    if deviation < SPLIT_TOL
        verdict = 'works';
    else
        verdict = 'FAILED';
    end
    fprintf('  Exchange coupling through the periodic boundary %s: %.3e between the halves, spread %.3e\n', ...
            verdict, deviation, spread);
    if deviation_free > SPLIT_CONTROL_MIN
        verdict = 'as expected';
    else
        verdict = 'FAILED';
    end
    fprintf('  The same halves without periodic boundaries are decoupled, %s: %.3e\n', ...
            verdict, deviation_free);

    checks(end+1) = struct('check', ...
        sprintf('%s: halves coupled through the periodic boundary', GRAIN_LABEL), ...
        'value', deviation, 'limit', SPLIT_TOL, 'passed', deviation < SPLIT_TOL);
    checks(end+1) = struct('check', ...
        sprintf('%s: both halves relax to the same direction', GRAIN_LABEL), ...
        'value', spread, 'limit', SPLIT_TOL, 'passed', spread < SPLIT_TOL);
    % The control has to fail the test above, which is what shows that the gap really does cut the
    % mesh in two. It is written as a ratio so that it fits the value < limit contract
    checks(end+1) = struct('check', ...
        sprintf('%s: halves decoupled without periodic boundaries (control)', GRAIN_LABEL), ...
        'value', SPLIT_CONTROL_MIN/max(deviation_free, 1e-300), 'limit', 1.0, ...
        'passed', deviation_free > SPLIT_CONTROL_MIN);

    % ---- Methods 2 and 3, on a smaller mesh of the same kind --------------------------------
    [opts_pts, opts_abc, opts_L] = build_grain_mesh(GRAIN_RES_OPERATOR, op_a, GRAIN_COUNT, ...
                                                    GRAIN_REFINEMENTS, GRAIN_OFFSET, GRAIN_SEED);
    A_pbc_g = full(exchange_matrix_mesh(opts_pts, opts_abc, opts_L, op_pbc, options.use_CUDA));
    A_free_g = full(exchange_matrix_mesh(opts_pts, opts_abc, opts_L, [0 0 0], options.use_CUDA));
    scale = max(abs(A_pbc_g(:)));
    row_sum = max(abs(sum(A_pbc_g, 2))) / scale;
    nnzRow = sum(A_pbc_g ~= 0, 2);

    fprintf('\nExchange operator on %d irregular cells, exchPBC = [%d %d %d]\n', ...
            size(opts_pts,1), op_pbc(1), op_pbc(2), op_pbc(3));
    fprintf('  matrix          : %d x %d, %d nonzeros\n', size(A_pbc_g,1), size(A_pbc_g,2), nnz(A_pbc_g));
    fprintf('  couplings/row   : %d (min) %d (max)\n', min(nnzRow), max(nnzRow));
    fprintf('  max |row sum|   : %.3e [%s]\n', row_sum, passFail(row_sum < ROW_SUM_TOL));
    % Reported, not tested: how much of the operator the linking changes at all, which is what
    % the supercell reference below has to reproduce. An implementation that quietly ignored
    % exchPBC would leave this at zero and would then miss the supercell by the same amount
    fprintf('  linking changes : %.3e of the largest coupling\n', ...
            max(abs(A_pbc_g(:) - A_free_g(:)))/scale);
    checks(end+1) = struct('check', ...
        sprintf('%s: exchange operator has vanishing row sums', GRAIN_LABEL), ...
        'value', row_sum, 'limit', ROW_SUM_TOL, 'passed', row_sum < ROW_SUM_TOL);

    grain_supercell = NaN;
    if any(op_pbc)
        [grain_supercell, nCopies] = supercell_deviation(opts_pts, opts_abc, opts_L, op_pbc, ...
                                                         options.use_CUDA);
        fprintf('  supercell       : %.3e away from the operator on %d copies of the mesh [%s]\n', ...
                grain_supercell, nCopies, passFail(grain_supercell < SUPERCELL_TOL));
        checks(end+1) = struct('check', ...
            sprintf('%s: periodic operator matches the supercell reference', GRAIN_LABEL), ...
            'value', grain_supercell, 'limit', SUPERCELL_TOL, ...
            'passed', grain_supercell < SUPERCELL_TOL);
    end

    grains = struct('pts', gpts_all, 'abc', gabc_all, 'grain', grainIdx, 'gap', gap, ...
                    'm0', m0_grain, 'M_pbc', M_pbc, 'M_free', M_free, 'deviation', deviation, ...
                    'spread', spread, 'deviation_free', deviation_free, ...
                    'supercell', grain_supercell, 'grid_L', gL, 'a', a, 'nGrains', GRAIN_COUNT);
end

%% Plot the results
if options.ShowTheResult
    results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
    if ~isfolder(results_dir)
        mkdir(results_dir);
    end

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1900 450]);
    ax = subplot(1, 3, 1, 'Parent', fig);
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
    ax2 = subplot(1, 3, 2, 'Parent', fig);
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

    % ---- The irregular mesh and the supercell reference ---------------------------------
    % Every entry here is a dimensionless deviation that has to stay below its own threshold,
    % except for the control, which has to stay above it
    ax3 = subplot(1, 3, 3, 'Parent', fig);
    barLabels = {}; barValues = []; barLimits = []; barAbove = [];
    if options.Grains
        barLabels{end+1} = 'Split mesh, halves';
        barValues(end+1) = grains.deviation;      barLimits(end+1) = SPLIT_TOL;         barAbove(end+1) = false;
        barLabels{end+1} = 'Split mesh, spread';
        barValues(end+1) = grains.spread;         barLimits(end+1) = SPLIT_TOL;         barAbove(end+1) = false;
        barLabels{end+1} = 'Split mesh, free (control)';
        barValues(end+1) = grains.deviation_free; barLimits(end+1) = SPLIT_CONTROL_MIN;  barAbove(end+1) = true;
    end
    for g = 1:numel(gridTypes)
        if ~isnan(supercells(g))
            barLabels{end+1} = sprintf('Supercell, %s', gridLabels{g}); %#ok<AGROW>
            barValues(end+1) = supercells(g);     barLimits(end+1) = SUPERCELL_TOL;      barAbove(end+1) = false; %#ok<AGROW>
        end
    end
    if options.Grains && ~isnan(grains.supercell)
        barLabels{end+1} = sprintf('Supercell, %s', GRAIN_LABEL);
        barValues(end+1) = grains.supercell;      barLimits(end+1) = SUPERCELL_TOL;      barAbove(end+1) = false;
    end

    if isempty(barValues)
        axis(ax3, 'off')
    else
        hold(ax3, 'on')
        xLow = min([barValues barLimits])/10;
        bh = barh(ax3, barValues);
        bh.FaceColor = 'flat';
        bh.BaseValue = xLow;
        for i = 1:numel(barValues)
            passed = xor(logical(barAbove(i)), barValues(i) < barLimits(i));
            if passed
                bh.CData(i,:) = [0.13 0.55 0.13];
            else
                bh.CData(i,:) = [0.7 0.13 0.13];
            end
            plot(ax3, barLimits(i)*[1 1], i + 0.35*[-1 1], 'k-', 'LineWidth', 1.5);
        end
        xHigh = 10*max([barValues barLimits]);
        set(ax3, 'XScale', 'log', 'YDir', 'reverse')
        xlim(ax3, [xLow, xHigh])
        % Without this a range of seventeen decades is given a single tick
        xticks(ax3, 10.^(ceil(log10(xLow)):4:floor(log10(xHigh))))
        yticks(ax3, 1:numel(barValues))
        yticklabels(ax3, barLabels)
        xlabel(ax3, 'Deviation, the black bars are the thresholds')
        title(ax3, 'Irregular mesh and supercell reference')
        box(ax3, 'on')
    end

    figure_path = fullfile(results_dir, 'periodic_exchange_test.png');
    exportgraphics(fig, figure_path, 'Resolution', 200);
    close(fig);
    fprintf('\nSaved figure to %s\n', figure_path);

    if options.Grains
        plot_grain_figure(grains, results_dir);
    end
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
% The exchange matrix of a regular lattice of cubes, as a dense array.

ntot = prod(res);
if strcmp(gridType, 'uniform')
    problem = DefaultMicroMagProblem(res(1), res(2), res(3));
    problem.grid_L = res*a;
    A = full(solve_for_exchange(problem, ntot, pbc, use_CUDA));
else
    A = full(exchange_matrix_mesh(build_prism_mesh(res, a), a*ones(ntot,3), res*a, pbc, use_CUDA));
end
end


function A = exchange_matrix_mesh(pts, abc, grid_L, pbc, use_CUDA)
% The exchange matrix of an unstructured mesh given by its centres and its side lengths.
%
% Kept sparse, because the supercell reference assembles matrices with up to 27 times as many
% rows as the one they are compared with.

ntot = size(pts,1);
problem = DefaultMicroMagProblem(ntot, 1, 1);
problem = problem.setMicroMagGridType('unstructuredPrisms');
problem.grid_pts = pts;
problem.grid_abc = abc;
problem.grid_L = grid_L;
A = solve_for_exchange(problem, ntot, pbc, use_CUDA);
end


function A = solve_for_exchange(problem, ntot, pbc, use_CUDA)
% Run a minimal simulation and return the assembled exchange matrix as a sparse array.
%
% Only the setup of the problem is of interest, so the simulation is run for a couple of trivial
% timesteps with the demagnetisation field and the external field switched off. The matrix comes
% back as the second output of the solver, in COO form, exported straight from MKL at exactly
% nnz entries.

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
A = sparse(rows+base, cols+base, vals, nr, nc);
end


function [deviation, nCopies] = supercell_deviation(pts, abc, grid_L, pbc, use_CUDA)
% Compare the periodic operator of a mesh with the free operator on a supercell of copies.
%
% The mesh is copied once in each direction along every periodic direction, the exchange operator
% is assembled on the resulting supercell without any periodic boundary conditions, and the rows
% of the central copy are folded back onto the tiles of the original mesh. Every tile of the
% central copy has exactly the neighbourhood that the periodic linking gives the corresponding
% tile of the original mesh - the stencil reaches one layer of cells and every copy is at least
% three cells thick - so the two operators have to agree entry by entry. Nothing about the mesh
% enters, so unlike the plane wave test this works on an irregular mesh as well.

ntot = size(pts,1);
A_pbc = exchange_matrix_mesh(pts, abc, grid_L, pbc, use_CUDA);

% All the copies of the mesh that a tile can be coupled to, including the diagonal ones: a tile
% in a corner of the domain is linked across two or three periodic boundaries at once
offsets = cell(1,3);
for d = 1:3
    if pbc(d)
        offsets{d} = [-1 0 1];
    else
        offsets{d} = 0;
    end
end
[X, Y, Z] = ndgrid(offsets{1}, offsets{2}, offsets{3});
shifts = [X(:) Y(:) Z(:)] .* grid_L(:)';
nCopies = size(shifts,1);
central = find(all(shifts == 0, 2), 1);

super_pts = zeros(nCopies*ntot, 3);
for c = 1:nCopies
    super_pts((c-1)*ntot + (1:ntot), :) = pts + shifts(c,:);
end
super_L = grid_L;
super_L(logical(pbc)) = 3*grid_L(logical(pbc));

A_sup = exchange_matrix_mesh(super_pts, repmat(abc, nCopies, 1), super_L, [0 0 0], use_CUDA);

% Fold the rows of the central copy back onto the original mesh. The copies are stored as
% consecutive blocks of ntot tiles, so column c of the supercell is tile mod(c, ntot) of the mesh
block = A_sup((central-1)*ntot + (1:ntot), :);
folded = sparse(ntot, ntot);
for c = 1:nCopies
    folded = folded + block(:, (c-1)*ntot + (1:ntot));
end

deviation = full( max(max(abs(folded - A_pbc))) / max(max(abs(A_pbc))) );   % max of a sparse matrix stays sparse
end


function [pts, abc, grid_L, grainIdx] = build_grain_mesh(res_base, a_base, nGrains, nRefine, offsetFraction, seed)
% Return the centres, the side lengths, the size and the grains of an irregular mesh.
%
% A base grid of res_base^3 cubes is refined recursively: a cell is split into eight whenever it
% is close enough to the boundary between two grains that the boundary can pass through it. The
% grains are the Voronoi cells of a set of random points, so the result is the mesh type used in
% python/experiments/Grain_perf, without the grain bookkeeping. The mesh is centred on the origin
% and is reproducible for a given seed.

rng(seed);
L = res_base*a_base;
% The grains are kept away from the faces of the domain, so that the refined cells at a grain
% boundary do not all end up sitting on the periodic boundaries
seeds = (0.15 + 0.7*rand(nGrains,3))*L;
offset = offsetFraction*a_base;

pts = zeros(0,3);
abc = zeros(0,3);
for k = 1:res_base
    for j = 1:res_base
        for i = 1:res_base
            [p, d] = refine_cell(([i j k] - 0.5)*a_base, a_base*[1 1 1], 0, nRefine, seeds, offset);
            pts = [pts; p]; %#ok<AGROW>
            abc = [abc; d]; %#ok<AGROW>
        end
    end
end

% The grain of a cell is the seed it is closest to, which is what makes the grains the Voronoi
% cells of the seeds. It is not used by the simulation, only to draw the mesh
distance2 = sum((reshape(pts, [], 1, 3) - reshape(seeds, 1, [], 3)).^2, 3);
[~, grainIdx] = min(distance2, [], 2);

pts = pts - L/2;
grid_L = [L L L];
end


function [pts, abc] = refine_cell(centre, sizes, level, maxLevel, seeds, offset)
% Split a cell into eight while it can be cut by a grain boundary

% The corner of a cell reaches half a diagonal beyond its centre, so a cell that is closer to a
% grain boundary than that is cut by the boundary and is refined
if level < maxLevel && grain_boundary_distance(centre, seeds) <= offset + norm(sizes)/2
    corners = [-1 -1 -1; 1 -1 -1; -1 1 -1; 1 1 -1; -1 -1 1; 1 -1 1; -1 1 1; 1 1 1]/4;
    pts = zeros(0,3);
    abc = zeros(0,3);
    for c = 1:8
        [p, d] = refine_cell(centre + corners(c,:).*sizes, sizes/2, level+1, maxLevel, seeds, offset);
        pts = [pts; p]; %#ok<AGROW>
        abc = [abc; d]; %#ok<AGROW>
    end
else
    pts = centre;
    abc = sizes;
end
end


function distance = grain_boundary_distance(centre, seeds)
% Distance from a point to the plane bisecting the two grains that are closest to it

d = sqrt(sum((seeds - centre).^2, 2));
[sorted, order] = sort(d);
distance = (sorted(2)^2 - sorted(1)^2) / (2*norm(seeds(order(1),:) - seeds(order(2),:)));
end


function mask = grain_gap_mask(pts, a_base)
% The middle layer of base cells of a grain mesh, which is left out to split it in two.
%
% The mesh is centred on the origin, so with an odd number of base cells per side the middle layer
% is the one spanning -a_base/2 to +a_base/2. A refined cell lies entirely within its base cell,
% so testing the centres picks out whole cells and leaves a gap that is one base cell wide, i.e.
% wider than any cell in the mesh. The exchange stencil reaches the cells that share a vertex with
% a face and therefore cannot bridge it.

mask = abs(pts(:,1)) < a_base/2;
end


function [deviation, spread, M_end] = split_deviations(pts, abc, grid_L, pbc, use_CUDA, params)
% Simulate a mesh that a gap has split in two and return how far from uniform it ends up.
%
% The two halves touch one another only through the periodic boundary, so they relax to a common
% direction if the linking works and to unrelated directions if it does not. Two measures are
% returned: the distance between the mean magnetisation of the two halves, which is what the
% linking decides, and the spread of the moments over the whole mesh, which is the same measure
% that the sections use. Both are in the units of |m| = 1 that MagTense returns and both are of
% order 1 when the halves are not coupled.

ntot = size(pts,1);
problem = DefaultMicroMagProblem(ntot, 1, 1);
problem = problem.setMicroMagGridType('unstructuredPrisms');
problem.grid_pts = pts;
problem.grid_abc = abc;
problem.grid_L = grid_L;
problem = problem.setUseCuda(use_CUDA);
problem = problem.setUseCVODE(false);
problem = problem.setUseDemag(false);
problem = problem.setMicroMagSolver('Dynamic');

problem.gamma = 2.21e5;
problem.alpha = params.eta;
problem.Ms = params.Ms*ones(ntot,1);
problem.A0 = params.Aex*ones(ntot,1);
problem.K0 = zeros(ntot,1);
problem.m0 = params.m0;
problem.exchPBC = int32(pbc);

HextFct = @(t) (t>=0)' * [0, 0, 0];
problem = problem.setHext( HextFct, linspace(0, params.t_end, 2) );
problem = problem.setTime( linspace(0, params.t_end, params.nTimesteps) );

solution = struct();
prob_struct = struct(problem);
solution = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );

M_npv = squeeze(solution.M(:,:,1,:));
M_end = squeeze(M_npv(end,:,:));

left = pts(:,1) < 0;
right = pts(:,1) > 0;
deviation = norm(mean(M_end(left,:), 1) - mean(M_end(right,:), 1));
spread = max( sqrt(sum((M_end - mean(M_end, 1)).^2, 2)) );
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


function plot_grain_figure(grains, results_dir)
% The mesh itself: the grains it is built from and the states the two simulations end in.
%
% A cross section is enough to see everything that matters, and the plane is placed away from the
% faces of the cells - those sit at multiples of half a base cell - so that it cuts every cell of
% the mesh it passes through exactly once.

pts = grains.pts; abc = grains.abc; gap = grains.gap;
a_base = grains.a; grid_L = grains.grid_L;

% Everything is drawn in nm, so that the panels do not carry an axis exponent that would land on
% top of the titles
nm = 1e-9;
pts = pts/nm; abc = abc/nm; grid_L = grid_L/nm; a_base = a_base/nm;

z0 = a_base/3;
inPlane = abs(pts(:,3) - z0) < abc(:,3)/2;
shown = inPlane & ~gap;
shownOfKept = inPlane(~gap);
cmap = coolwarm_map(256);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1150 1000]);

% ---- The mesh and the grains it was refined around ---------------------------------------
ax = subplot(2, 2, 1, 'Parent', fig);
start_cross_section(ax, grid_L, a_base);
grainColours = lines(max(grains.grain));
draw_cells(ax, pts(shown,:), abc(shown,:), grainColours(grains.grain(shown),:), 'k');
% The cells that are left out are drawn as empty boxes, since they are what splits the mesh
dropped = inPlane & gap;
draw_cells(ax, pts(dropped,:), abc(dropped,:), [], [0.5 0.5 0.5]);
finish_cross_section(ax, grid_L, ...
    sprintf('%d grains, refined at the grain boundaries', grains.nGrains));
text(ax, 0, 0.47*grid_L(2), 'gap', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
for s = [-1 1]
    text(ax, s*0.55*grid_L(1), 0, 'periodic boundary', 'Color', [0 0.39 0], 'Rotation', 90, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 8);
end

% ---- The state the simulations start from, and the two they end in -----------------------
states = {grains.m0,     {'Initial state, random directions', ''}; ...
          grains.M_pbc,  {'Relaxed with exchPBC = [1 0 0]', ...
                          sprintf('Distance between the mean of the two halves: %.1e', ...
                                  grains.deviation)}; ...
          grains.M_free, {'Relaxed with exchPBC = [0 0 0], the control', ...
                          sprintf('Distance between the mean of the two halves: %.1e', ...
                                  grains.deviation_free)}};
for i = 1:3
    ax = subplot(2, 2, i+1, 'Parent', fig);
    M = states{i,1};
    draw_magnetisation(ax, pts(shown,:), abc(shown,:), M(shownOfKept,:), grid_L, a_base, ...
                       states{i,2}, cmap);
    if i == 3
        set(ax, 'CLim', [-1 1]);
        colormap(ax, cmap);
        cb = colorbar(ax);
        cb.Label.String = 'm_z';
    end
end

sgtitle(fig, sprintf(['The irregular mesh of the periodic exchange test, ' ...
                      'cross section at z = %.2f nm'], z0));

figure_path = fullfile(results_dir, 'periodic_exchange_test_grains.png');
exportgraphics(fig, figure_path, 'Resolution', 200);
close(fig);
fprintf('Saved figure to %s\n', figure_path);
end


function draw_cells(ax, pts, abc, colours, edgeColour)
% Draw the cross sections of a set of cells as rectangles in the xy plane

n = size(pts,1);
if n == 0
    return
end
xs = [pts(:,1)-abc(:,1)/2, pts(:,1)+abc(:,1)/2, pts(:,1)+abc(:,1)/2, pts(:,1)-abc(:,1)/2].';
ys = [pts(:,2)-abc(:,2)/2, pts(:,2)-abc(:,2)/2, pts(:,2)+abc(:,2)/2, pts(:,2)+abc(:,2)/2].';
vertices = [xs(:), ys(:)];
faces = reshape(1:4*n, 4, n).';

if isempty(colours)
    patch(ax, 'Faces', faces, 'Vertices', vertices, 'FaceColor', 'none', ...
          'EdgeColor', edgeColour, 'LineWidth', 0.4);
else
    patch(ax, 'Faces', faces, 'Vertices', vertices, 'FaceVertexCData', colours, ...
          'FaceColor', 'flat', 'EdgeColor', edgeColour, 'LineWidth', 0.4);
end
end


function draw_magnetisation(ax, pts, abc, M, grid_L, a_base, titleText, cmap)
% One cross section of the split mesh, with the magnetisation drawn on the cells.
%
% The colour of a cell is the out of plane component and the arrow is the part of the moment that
% lies in the plane, which is the usual way of looking at a micromagnetic state. The two halves of
% the mesh carry the same direction only if they are exchange coupled through the periodic
% boundary, so the panels can be told apart at a glance.

start_cross_section(ax, grid_L, a_base);
mz = min(max(M(:,3), -1), 1);
colours = interp1(linspace(-1, 1, size(cmap,1)), cmap, mz);
draw_cells(ax, pts, abc, colours, [0.4 0.4 0.4]);

% The arrows are scaled by the cell they belong to, so that the small cells at the grain
% boundaries do not end up with arrows reaching across their neighbours, and they are centred on
% the cell rather than starting at it
u = 0.9*M(:,1).*abc(:,1);
v = 0.9*M(:,2).*abc(:,1);
quiver(ax, pts(:,1)-u/2, pts(:,2)-v/2, u, v, 0, 'Color', 'k', 'LineWidth', 0.8, ...
       'MaxHeadSize', 0.4);
finish_cross_section(ax, grid_L, titleText);
end


function start_cross_section(ax, grid_L, a_base)
% The gap that splits the mesh, drawn first so that the cells land on top of it

hold(ax, 'on')
fill(ax, [-1 1 1 -1]*a_base/2, [-1 -1 1 1]*0.52*grid_L(2), [0.85 0.85 0.85], ...
     'EdgeColor', 'none');
end


function finish_cross_section(ax, grid_L, titleText)
% The two boundaries the halves are linked through, and the frame of the panel

for s = [-1 1]
    plot(ax, s*grid_L(1)/2*[1 1], 0.52*grid_L(2)*[-1 1], '--', 'Color', [0 0.39 0], ...
         'LineWidth', 1.5);
end
set(ax, 'DataAspectRatio', [1 1 1])
xlim(ax, 0.62*grid_L(1)*[-1 1])
ylim(ax, 0.52*grid_L(2)*[-1 1])
xlabel(ax, 'x [nm]')
ylabel(ax, 'y [nm]')
title(ax, titleText, 'FontSize', 10)
box(ax, 'on')
end


function cmap = coolwarm_map(n)
% A blue to white to red colour map for the out of plane magnetisation

anchors = [0.230 0.299 0.754; 0.865 0.865 0.865; 0.706 0.016 0.150];
cmap = interp1([0 0.5 1], anchors, linspace(0, 1, n));
end
