function checks = macrogeometry_supercell_test(options)
%MACROGEOMETRY_SUPERCELL_TEST
% Test that the macrogeometry method reproduces an explicitly replicated domain, cell by cell.
% This is the MATLAB counterpart of
% python/examples/micromagnetism/MagTense_tests/macrogeometry_supercell_test.py and uses the same
% geometry, the same cases and the same acceptance limits.
%
% SHORT EXPLANATION : A periodic domain must give exactly the demagnetisation field of a large
% domain built from the same number of copies placed side by side.
%
% LONG EXPLANATION :
% The macrogeometry method adds the field of 2*n_macro copies of the simulated domain, shifted by
% multiples of shiftVec, to the field of the domain itself. When shiftVec equals the size of the
% domain, grid_L, the copies tile space without gaps or overlaps, so the result has to be identical
% to that of a single domain that is (2*n_macro + 1) times larger along every periodic direction
% and carries the same magnetisation pattern in every block. The demagnetisation field in the
% central block of that supercell is computed without any periodic boundary conditions, so the
% comparison tests the placement of the copies directly, without any analytical model in between.
%
% This is the gapless case that macrogeometry_PBC_test does not cover: that test models a sparse
% chain of separated particles, where shiftVec is much larger than the domain. A shift that is off
% by a cell - for instance the distance between the centres of the two end cells, (n - 1)*dx,
% rather than the full length n*dx - makes the copies overlap by one cell layer. MagTense rejects
% such a shiftVec outright, so the last check runs that case in a separate MATLAB process and
% requires it to stop with a message rather than to return a field.
%
% The domain has a different cell size along x, y and z and a random magnetisation, so a mix-up
% between directions or a misplaced copy anywhere shows up. Both the point demagnetisation tensor
% and the tensor averaged over the observation cell are tested. As a control the same comparison
% is made without the macrogeometry: the periodic copies have to change the field by much more
% than the tolerance, otherwise the test could not tell a working implementation from a missing
% one.
%
% Returns a struct array of checks with the fields 'check', 'value', 'limit' and 'passed',
% where a check passes when value < limit. That is the contract used by testMagTenseFunctions.m.

arguments
    options.ShowTheResult {mustBeNumericOrLogical} = true;   % Save the validation figure
    options.use_CUDA {mustBeNumericOrLogical} = false;
end

addpath('../../../MEX_files');
addpath('../../../util');

%% Settings

res = [3 2 2];                      % Cells along x, y and z in the periodic domain
cellSize = [2e-9 3e-9 5e-9];        % Cell size along x, y and z [m], deliberately unequal
grid_L = res .* cellSize;           % Size of the periodic domain [m]

Ms = 8e5;                           % Saturation magnetisation [A/m]
seed = 3;                           % Seed of the random magnetisation

% {label, periodic directions, copies on each side, averaged tensor}. Copies on each side of the
% domain along every periodic direction; the supercell is (2*n + 1) domains long along those.
cases = {
    'x',                            [1 0 0], 4, false
    'y',                            [0 1 0], 4, false
    'z',                            [0 0 1], 4, false
    'x, y and z',                   [1 1 1], 2, false
    'x, y and z, averaged tensor',  [1 1 1], 2, true
    };

% The demagnetisation tensor is stored in single precision, so the two fields agree to about 1e-7
% of the largest field component. The limit leaves room for the different summation order.
match_tol = 1e-5;
% The periodic copies have to change the field by at least this much, relative to the largest
% field component, for the comparison above to mean anything
control_min = 1e-2;

%% Run the cases
fprintf('Periodic domain: %d x %d x %d cells of %g x %g x %g nm, shiftVec = grid_L\n', ...
        res(1), res(2), res(3), cellSize(1)*1e9, cellSize(2)*1e9, cellSize(3)*1e9);

ntot = prod(res);
m0 = random_magnetisation(ntot, seed);

checks = struct('check', {}, 'value', {}, 'limit', {}, 'passed', {});
mismatch = zeros(size(cases,1), 1);
control = zeros(size(cases,1), 1);
for c = 1:size(cases,1)
    [label, pbc, nCopies, useAvg] = cases{c,:};
    pbc = logical(pbc);
    copies = ones(1,3);
    copies(pbc) = 2*nCopies + 1;

    % The periodic domain
    n_macro = zeros(1,3);
    n_macro(pbc) = nCopies;
    shiftVec = zeros(1,3);
    shiftVec(pbc) = grid_L(pbc);
    H_pbc = demag_field(res, grid_L, m0, Ms, useAvg, options.use_CUDA, n_macro, shiftVec);

    % The same domain without any copies, the control
    H_free = demag_field(res, grid_L, m0, Ms, useAvg, options.use_CUDA, [0 0 0], [0 0 0]);

    % The supercell. Cells are numbered with x running fastest, as setupGrid in
    % LandauLifshitzEquationSolver.f90 numbers them, which is MATLAB's own column-major order, so
    % reshaping to (nx, ny, nz, 3) gives the grid itself and repmat lays the copies side by side
    res_super = res .* copies;
    m0_super = repmat(reshape(m0, [res 3]), [copies 1]);
    H_super = demag_field(res_super, grid_L .* copies, reshape(m0_super, [], 3), Ms, useAvg, ...
                          options.use_CUDA, [0 0 0], [0 0 0]);

    % The central block of the supercell coincides with the periodic domain: the supercell is
    % centred on the origin just like the domain, and it is an odd number of domains long
    offset = zeros(1,3);
    offset(pbc) = nCopies * res(pbc);
    H_super_grid = reshape(H_super, [res_super 3]);
    H_central = reshape(H_super_grid(offset(1) + (1:res(1)), offset(2) + (1:res(2)), ...
                                     offset(3) + (1:res(3)), :), ntot, 3);

    scale = max(abs(H_central(:)));
    mismatch(c) = max(abs(H_pbc(:) - H_central(:))) / scale;
    control(c) = max(abs(H_free(:) - H_central(:))) / scale;

    ok = mismatch(c) < match_tol;
    fprintf(['  periodic along %s, %d copies on each side (%d cells in the supercell): ' ...
             'mismatch %.2e [%s], the copies change the field by %.2e\n'], label, nCopies, ...
             prod(res_super), mismatch(c), passFail(ok), control(c));
    checks(end+1) = struct('check', sprintf('periodic along %s: field matches the supercell', label), ...
        'value', mismatch(c), 'limit', match_tol, 'passed', ok); %#ok<AGROW>
    % Written as a ratio so that it fits the value < limit contract
    checks(end+1) = struct('check', sprintf('periodic along %s: copies change the field (control)', label), ...
        'value', control_min / max(control(c), 1e-300), 'limit', 1.0, ...
        'passed', control(c) > control_min); %#ok<AGROW>
end

%% Overlapping copies have to be rejected
[rejected, output] = overlapping_copies_rejected(res, cellSize);
if rejected
    fprintf('  shiftVec = (n - 1)*dx, copies overlapping by one cell: rejected\n');
else
    fprintf('  shiftVec = (n - 1)*dx, copies overlapping by one cell: NOT REJECTED\n');
    disp(output(max(1, end-2000):end));
end
checks(end+1) = struct('check', 'overlapping copies are rejected', ...
    'value', double(~rejected), 'limit', 0.5, 'passed', rejected);

%% Plot the result
if options.ShowTheResult
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 700 450]);
    ax = axes(fig);
    hold(ax, 'on')
    b = bar(ax, [mismatch control]);
    b(1).FaceColor = [0.13 0.55 0.13];
    b(2).FaceColor = [0.5 0.5 0.5];
    b(1).DisplayName = 'Periodic vs supercell';
    b(2).DisplayName = 'No copies vs supercell (control)';
    yline(ax, match_tol, '--', 'Color', [0.13 0.55 0.13], 'DisplayName', 'Tolerance');
    yline(ax, control_min, ':', 'Color', [0.5 0.5 0.5], 'DisplayName', 'Control minimum');
    set(ax, 'YScale', 'log')
    xticks(ax, 1:size(cases,1))
    xticklabels(ax, cases(:,1))
    xlabel(ax, 'Periodic directions')
    ylabel(ax, 'max |\DeltaH| / max |H|')
    legend(ax, 'Location', 'best')
    box(ax, 'on')

    results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
    if ~isfolder(results_dir)
        mkdir(results_dir);
    end
    figure_path = fullfile(results_dir, 'macrogeometry_supercell_test.png');
    exportgraphics(fig, figure_path, 'Resolution', 200);
    close(fig);
    fprintf('Saved figure to %s\n', figure_path);
end

if nargout == 0
    if all([checks.passed])
        disp('macrogeometry_supercell_test PASSED')
    else
        disp('macrogeometry_supercell_test FAILED')
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


function m = random_magnetisation(n, seed)
% Reproducible random unit vectors, one per cell
rng(seed);
m = 2*rand(n, 3) - 1;
m = m ./ vecnorm(m, 2, 2);
end


function H = demag_field(res, grid_L, m0, Ms, useAvg, use_CUDA, n_macro, shiftVec)
% The demagnetisation field of a uniform grid with magnetisation m0, as an (ntot, 3) array.
%
% Only the field of the initial state is needed, so the simulation is run for a vanishing time
% with the precession switched off and the field is read at the first output time.

ntot = prod(res);
problem = DefaultMicroMagProblem(res(1), res(2), res(3));
problem.grid_L = grid_L;
problem = problem.setUseCuda(use_CUDA);
problem = problem.setUseCVODE(false);
problem = problem.setUseDemag(true);
problem = problem.setMicroMagSolver('Dynamic');
problem.useAvgN = int32(useAvg);
% An octree is pointless for grids this small, and the FMM path ignores the macrogeometry
problem.use_fmm = int32(0);
problem.ReturnHall = int32(1);

problem.gamma = 0;
problem.alpha = 1e3;
problem.Ms = Ms*ones(ntot,1);
problem.A0 = 1e-20;
problem.K0 = zeros(ntot,1);
problem.m0 = m0;

problem.n_macro = int32(n_macro);
problem.shiftVec = shiftVec;

HextFct = @(t) (t>=0)' * [0, 0, 0];
problem = problem.setHext( HextFct, linspace(0, 1e-15, 2) );
problem = problem.setTime( linspace(0, 1e-15, 2) );

solution = struct();
prob_struct = struct(problem);
solution = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );

% The returned magnetisation at the first output time, normalised to |m| = 1, has to be the
% initial state, otherwise the fields compared would belong to different magnetisations
M_first = reshape(solution.M(1,:,1,:), ntot, 3);
assert(max(abs(M_first(:) - m0(:))) < 1e-6, ...
       'macrogeometry_supercell_test:initialState', ...
       'The first output time does not hold the initial magnetisation');
H = reshape(solution.H_dem(1,:,1,:), ntot, 3);
end


function [rejected, output] = overlapping_copies_rejected(res, cellSize)
% Run a domain whose copies overlap by one cell layer in a separate MATLAB process.
%
% The shift is the distance between the centres of the two end cells instead of the full length
% of the domain. MagTense has to stop on that rather than return a field, and stopping takes the
% whole process with it, hence the separate process.

here = fileparts(mfilename('fullpath'));
mexDir = fullfile(here, '..', '..', '..', 'MEX_files');
utilDir = fullfile(here, '..', '..', '..', 'util');

script = strjoin({
    sprintf('addpath(''%s''); addpath(''%s'');', mexDir, utilDir)
    sprintf('res = [%d %d %d]; cellSize = [%.17g %.17g %.17g];', res, cellSize)
    'problem = DefaultMicroMagProblem(res(1), res(2), res(3));'
    'problem.grid_L = res .* cellSize;'
    'problem = problem.setUseCuda(false); problem = problem.setUseCVODE(false);'
    'problem = problem.setUseDemag(true); problem = problem.setMicroMagSolver(''Dynamic'');'
    'problem.useAvgN = int32(0); problem.use_fmm = int32(0);'
    'problem.gamma = 0; problem.alpha = 1e3; problem.Ms = 8e5*ones(prod(res),1);'
    'problem.A0 = 1e-20; problem.K0 = zeros(prod(res),1); problem.m0 = repmat([1 0 0], prod(res), 1);'
    'problem.n_macro = int32([2 0 0]); problem.shiftVec = [(res(1)-1)*cellSize(1) 0 0];'
    'problem = problem.setHext(@(t) (t>=0)'' * [0 0 0], linspace(0, 1e-15, 2));'
    'problem = problem.setTime(linspace(0, 1e-15, 2));'
    'solution = problem.MagTenseLandauLifshitzSolver_mex(struct(problem), struct());'
    'disp(''OVERLAP NOT DETECTED'')'
    }, ' ');

scriptFile = [tempname '.m'];
[scriptDir, scriptName] = fileparts(scriptFile);
% A function or script name has to start with a letter
scriptName = ['mtcheck_' regexprep(scriptName, '[^A-Za-z0-9_]', '_')];
scriptFile = fullfile(scriptDir, [scriptName '.m']);
fid = fopen(scriptFile, 'w');
fprintf(fid, '%s\n', script);
fclose(fid);
cleanup = onCleanup(@() delete(scriptFile));

matlabExe = fullfile(matlabroot, 'bin', 'matlab');
command = sprintf('"%s" -batch "cd(''%s''); %s"', matlabExe, scriptDir, scriptName);
[status, output] = system(command);
% The message of the error stop in checkMacrogeometrySpacing, so that the child failing for any
% other reason - a missing licence, a path problem - is not taken for a rejection
rejected = status ~= 0 && contains(output, 'the macrogeometry copies overlap') ...
           && ~contains(output, 'OVERLAP NOT DETECTED');
end
