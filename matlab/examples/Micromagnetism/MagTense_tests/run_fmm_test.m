function [metrics, allPassed] = run_fmm_test(options)
%RUN_FMM_TEST Validate the FMM demag field against the direct calculation.
%
%   The MATLAB counterpart of python/tests/fmm_test. The same unstructured
%   prism problem is solved twice: once with the direct demag tensor
%   (use_fmm = 0) as the reference, and once per FMM configuration with the
%   short-circuit disabled so FMM is genuinely exercised. The FMM fields are
%   then compared against the reference using the same metrics and the same
%   thresholds as the Python test, so the two are directly comparable.
%
%   Errors if any configuration falls outside tolerance, so that
%   "matlab -batch run_fmm_test" returns a non-zero exit code in CI.
%
%   Example:
%       [metrics, ok] = run_fmm_test('USE_CUDA', false);

arguments
    options.USE_CUDA (1,1) logical = true
    %--- Multipole expansion terms to sweep, matching the Python test
    options.nterms (1,:) double = [6 10 15]
    options.nlmax (1,1) double = 2
    %--- Timestep the error metrics are evaluated at, over all cells. The
    %    Python test calls this target_cell_idx = 9, but H_dem is
    %    (nt, ntot, nt_Hext, 3), so index 9 of 10 is the final timestep, not
    %    a cell. 0 means "last".
    options.target_time_idx (1,1) double = 0
    options.nt (1,1) double = 10
    options.t_end (1,1) double = 1e-9
    options.mesh_file char = ''
    options.ShowTheResult (1,1) logical = false
end

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..', '..', 'MEX_files'));
addpath(fullfile(thisDir, '..', '..', '..', 'util'));

if isempty(options.mesh_file)
    options.mesh_file = fullfile(thisDir, '..', ...
        'mumag_micromag_Std_problem_4', ...
        'Std_prob_4_unstructured_mesh_grains_6_res_80_20_ref_2.mat');
end

mu0 = 4 * pi * 1e-7;

%% ---------------------------------------------------------------- Problem
loaded = load(options.mesh_file);
mesh = loaded.mesh;
ntot = size(mesh.pos_out, 1);
[~, meshName, meshExt] = fileparts(options.mesh_file);
fprintf('FMM test: %d unstructured prisms from %s\n', ntot, [meshName meshExt]);

    function problem = buildProblem(useFmm, nterms)
        problem = DefaultMicroMagProblem(ntot, 1, 1);
        problem = problem.setMicroMagGridType('unstructuredPrisms');
        problem.grid_pts = mesh.pos_out;
        problem.grid_abc = mesh.dims_out;

        problem = problem.setMicroMagDemagApproximation('none');
        problem = problem.setUseCuda(options.USE_CUDA);
        problem = problem.setUseCVODE(false);

        %--- Without this the individual field contributions, H_dem among
        %    them, are not returned to MATLAB at all.
        problem.ReturnHall = int32(1);

        %--- Material properties, as in Standard_problem_4_unstructured_cart
        problem.alpha = 4.42e3;
        problem.gamma = 0;
        problem.Ms = 8e5 * ones(ntot, 1);
        problem.K0 = 0 * ones(ntot, 1);
        problem.A0 = 1.3e-11 * ones(ntot, 1);
        problem.exch_weigh = 8;
        problem = problem.setMicroMagExchMethod('DirectLaplacianNeumann');
        problem = problem.setMicroMagExchInterpn('Extended');

        problem.m0(:, 1) = 1 / sqrt(3);
        problem.m0(:, 2) = 1 / sqrt(3);
        problem.m0(:, 3) = 1 / sqrt(3);

        problem = problem.setTime(linspace(0, options.t_end, options.nt));
        problem.setTimeDis = int32(10);

        HystDir = 1 / mu0 * [1, 1, 1];
        HextFct = @(t) (1e-9 - t)' .* HystDir .* (t < 1e-9)';
        problem = problem.setHext(HextFct, linspace(0, options.t_end, 2 * options.nt));

        %--- FMM configuration. fmm_short = 0 disables the short-circuit that
        %    would otherwise fall back to the direct tensor below fmm_min_n.
        if useFmm
            problem.use_fmm = int32(1);
            problem.fmm_short = int32(0);
            problem.nlmax = int32(options.nlmax);
            problem.fmm_nterms = int32(nterms);
        else
            problem.use_fmm = int32(0);
        end
    end

    function sol = solveProblem(problem)
        sol = struct();
        sol = problem.MagTenseLandauLifshitzSolver_mex(struct(problem), sol);
    end

%% -------------------------------------------------------------- Reference
fprintf('Running direct (use_fmm = 0) reference...\n');
refSol = solveProblem(buildProblem(false, 0));

%% ------------------------------------------------------------- FMM sweeps
%--- H_dem and M are (nt, ntot, nt_Hext, 3); slice one timestep over all cells
tIdx = options.target_time_idx;
if tIdx <= 0
    tIdx = size(refSol.H_dem, 1);
end
metrics = struct('Config', {}, 'Hdem_rel_avg', {}, 'Hdem_rel_max', {}, ...
    'Mout_max', {}, 'Hdem_rel_max_alltimes', {});

for nterms = options.nterms
    label = sprintf('FMM_L%d_N%d', options.nlmax, nterms);
    fprintf('Running %s...\n', label);
    sol = solveProblem(buildProblem(true, nterms));

    %--- Metrics at one timestep over all cells, mirroring the Python
    %    calculate_error_metrics
    refH = squeeze(refSol.H_dem(tIdx, :, 1, :));
    fmmH = squeeze(sol.H_dem(tIdx, :, 1, :));
    dH = abs(refH - fmmH);
    normH = mean(abs(refH(:)));

    refM = squeeze(refSol.M(tIdx, :, 1, :));
    fmmM = squeeze(sol.M(tIdx, :, 1, :));

    %--- Over every timestep as well, reported but not gated on: the single
    %    slice above can flatter or unfairly penalise a run.
    dHall = abs(refSol.H_dem - sol.H_dem);
    normHall = mean(abs(refSol.H_dem(:)));

    %--- A bit-identical result means FMM never actually engaged and the run
    %    silently fell back to the direct tensor, in which case every metric
    %    below is 0 and the test would "pass" while testing nothing. This has
    %    happened for real: a MEX built without FMM3D, since buildMagTenseMEX
    %    only rebuilds the solver matching its USE_CUDA setting.
    if max(dHall(:)) == 0
        error('MagTense:fmmNotEngaged', ...
            ['%s produced results bit-identical to the direct reference, so FMM ' ...
             'did not run. Check that the MEX in use was built with USE_FMM3D ' ...
             'and look for "built tree with" in the solver output.'], label);
    end

    m = struct();
    m.Config = label;
    m.Hdem_rel_avg = mean(dH(:)) / normH;
    m.Hdem_rel_max = max(dH(:)) / normH;
    m.Mout_max = max(abs(refM(:) - fmmM(:)));
    m.Hdem_rel_max_alltimes = max(dHall(:)) / normHall;
    metrics(end + 1) = m; %#ok<AGROW>
end

%% ------------------------------------------------------------- Validation
%--- Same thresholds as python/tests/fmm_test/run_fmm_test.py
thresholds = struct( ...
    'L2', struct('Hdem_rel_avg', 1e-3, 'Hdem_rel_max', 0.05, 'Mout_max', 0.03), ...
    'L3', struct('Hdem_rel_avg', 1e-2, 'Hdem_rel_max', 0.40, 'Mout_max', 0.25));

limitKey = sprintf('L%d', options.nlmax);
if ~isfield(thresholds, limitKey)
    limitKey = 'L2';
end
limits = thresholds.(limitKey);
checks = fieldnames(limits);

fprintf('\n%-16s | %-6s | %-10s | %s\n', 'Configuration', 'Status', 'Limit Used', 'Errors Found');
fprintf('%s\n', repmat('-', 1, 85));

allPassed = true;
for k = 1:numel(metrics)
    failures = strings(0);
    for c = 1:numel(checks)
        name = checks{c};
        value = metrics(k).(name);
        if value > limits.(name)
            failures(end + 1) = sprintf('%s(%.2e > %.2e)', name, value, limits.(name)); %#ok<AGROW>
        end
    end

    if isempty(failures)
        status = 'PASS';
        detail = 'Within tolerance';
    else
        status = 'FAIL';
        detail = strjoin(failures, ', ');
        allPassed = false;
    end

    fprintf('%-16s | %-6s | %-10s | %s\n', metrics(k).Config, status, limitKey, detail);
end
fprintf('%s\n', repmat('-', 1, 85));

for k = 1:numel(metrics)
    fprintf('  %-16s Hdem_rel_avg=%.3e  Hdem_rel_max=%.3e  Mout_max=%.3e  (all-timestep Hdem_rel_max=%.3e)\n', ...
        metrics(k).Config, metrics(k).Hdem_rel_avg, metrics(k).Hdem_rel_max, ...
        metrics(k).Mout_max, metrics(k).Hdem_rel_max_alltimes);
end

if allPassed
    fprintf('\nOVERALL FMM TEST: SUCCESS\n');
else
    error('MagTense:fmmTestFailed', ...
        'FMM results fell outside tolerance for at least one configuration.');
end

end
