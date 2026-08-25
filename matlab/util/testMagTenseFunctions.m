%% testMagTenseFunctions
% Combined test suite for the MATLAB version of MagTense.
%
% The magnetostatic validation scripts are compared with FEM simulations and the
% micromagnetic mumag standard problems with their published solutions. Every test
% produces one or more numerical measures, each of which is compared with the acceptance
% limit that belongs to it. All tests are run to the end and the outcome is shown as a
% table in the command window and as a colour coded overview figure, so that a single
% failure does not hide the state of everything else.
%
% This is the counterpart of python/examples/micromagnetism/testMagTenseFunctions.py. The
% micromagnetic tests in examples/Micromagnetism/MagTense_tests mirror the python ones one
% for one, with the same geometries and the same acceptance limits. Beyond those, the
% magnetostatic validations and the standard problem 6 direction and mesh variants only
% exist here, while standard problem 3 only exists in the python suite.
%
% Every test returns a struct array of checks with the fields 'check', 'value', 'limit'
% and 'passed', where a check passes when value < limit.
%
% Author:  Rasmus Bjoerk

clear variables

% mfilename works both when the file is run from the editor and from a batch session,
% unlike matlab.desktop.editor.getActiveFilename which needs an open editor
util_dir = fileparts(mfilename('fullpath'));
addpath(util_dir)

save_figure = true;                                  % Also write the overview to a .png
results_dir = fullfile(util_dir, 'results');
figure_path = fullfile(results_dir, 'testMagTenseFunctions_overview.png');

magnetostatics_dir = fullfile(util_dir, '..', 'examples', 'Magnetostatics');
micromagnetism_dir = fullfile(util_dir, '..', 'examples', 'Micromagnetism');

% The integrated field errors are in percent, so 5 means the MagTense and the FEM curve
% enclose an area of 5 % of the area under the FEM curve
field_error_limit = 5;

% The workflow runs this suite once per CUDA and CVODE variant of the MEX files, so both
% switches are taken from the environment: MAGTENSE_TEST_CUDA and MAGTENSE_TEST_CVODE, each 0
% or 1. A variable that is not set leaves every test with the default it has on its own, which
% is what a plain call from the editor gets. MAGTENSE_TESTS runs only the tests it names and
% MAGTENSE_SKIP runs all but those, for the legs of the workflow that cannot afford the whole
% suite; standard problem 6 alone is four fifths of its running time.
cuda_args = environmentFlag('use_CUDA', 'MAGTENSE_TEST_CUDA');
cvode_args = environmentFlag('use_CVODE', 'MAGTENSE_TEST_CVODE');

%% The tests
% Each row is {name, function returning the checks, one line description}

tests = {
    'magnetostatics_cylindrical_slice_1', ...
        @() validationChecks(fullfile(magnetostatics_dir, 'Validation_field_cylindrical_slice'), ...
                             'MagTense_Validation_cylindrical_slice_example_1', field_error_limit), ...
        'Field of a cylindrical slice vs FEM, example 1'
    'magnetostatics_cylindrical_slice_2', ...
        @() validationChecks(fullfile(magnetostatics_dir, 'Validation_field_cylindrical_slice'), ...
                             'MagTense_Validation_cylindrical_slice_example_2', field_error_limit), ...
        'Field of a cylindrical slice vs FEM, example 2'
    'magnetostatics_prism', ...
        @() validationChecks(fullfile(magnetostatics_dir, 'Validation_field_prism'), ...
                             'MagTense_Validation_prism', field_error_limit), ...
        'Field of a rectangular prism vs FEM'
    'magnetostatics_sphere', ...
        @() validationChecks(fullfile(magnetostatics_dir, 'Validation_field_sphere'), ...
                             'MagTense_Validation_sphere', field_error_limit), ...
        'Field of a sphere vs FEM'
    'magnetostatics_spheroid', ...
        @() validationChecks(fullfile(magnetostatics_dir, 'Validation_field_spheroid'), ...
                             'MagTense_Validation_spheroid', field_error_limit), ...
        'Field of a spheroid vs FEM'
    'magnetostatics_tetrahedron', ...
        @() validationChecks(fullfile(magnetostatics_dir, 'Validation_field_tetrahedron'), ...
                             'MagTense_Validation_tetrahedron', field_error_limit), ...
        'Field of a tetrahedron vs FEM'
    'macrogeometry_PBC_test', ...
        @() micromagTestChecks(fullfile(micromagnetism_dir, 'MagTense_tests'), ...
                               'macrogeometry_PBC_test', cuda_args), ...
        'Periodic boundaries by the macrogeometry method, along x, y and z'
    'periodic_exchange_test', ...
        @() micromagTestChecks(fullfile(micromagnetism_dir, 'MagTense_tests'), ...
                               'periodic_exchange_test', cuda_args), ...
        'Periodic exchange coupling, uniform grid, unstructured mesh and grain mesh'
    'shape_correction_test', ...
        @() micromagTestChecks(fullfile(micromagnetism_dir, 'MagTense_tests'), ...
                               'shape_correction_test', cuda_args), ...
        'Shape correction field of the sample geometry'
    'temperature_test', ...
        @() micromagTestChecks(fullfile(micromagnetism_dir, 'MagTense_tests'), ...
                               'temperature_test', cuda_args), ...
        'Thermal fluctuations against the analytical angular diffusion'
    'std_problem_4', ...
        @() stdProblem4Checks(fullfile(micromagnetism_dir, 'mumag_micromag_Std_problem_4'), ...
                              [cuda_args, cvode_args]), ...
        'mumag standard problem 4 against the published mean solutions'
    'std_problem_6', ...
        @() stdProblem6Checks(fullfile(micromagnetism_dir, 'mumag_micromag_Std_problem_6'), ...
                              [cuda_args, cvode_args]), ...
        'mumag standard problem 6, domain wall depinning fields'
    };

%% Run them

tests = selectTests(tests, getenv('MAGTENSE_TESTS'), getenv('MAGTENSE_SKIP'));
records = runAllTests(tests);

%% Report

printResultTable(records)
plotResultOverview(records, figure_path, save_figure, results_dir)

failed = [records([records.status] ~= "PASS").name];
if ~isempty(failed)
    error('testMagTenseFunctions:testsFailed', ...
          'The following tests did not pass: %s', strjoin(failed, ', '));
end
disp('All tests passed')


%% ----------------------------------------------------------------------------------
%  Running the tests
%  ----------------------------------------------------------------------------------

function args = environmentFlag(name, variable)
    % {name, value} for a name value argument taken from an environment variable, or {} when
    % the variable is not set, so that the test keeps the default it has on its own

    text = strtrim(getenv(variable));
    if isempty(text)
        args = {};
        return
    end

    switch lower(text)
        case {'1', 'true', 'yes', 'on'}
            args = {name, true};
        case {'0', 'false', 'no', 'off'}
            args = {name, false};
        otherwise
            error('testMagTenseFunctions:badFlag', ...
                  '%s has to be 0 or 1, not "%s"', variable, text);
    end
end

function tests = selectTests(tests, selection, skipped)
    % Keep the tests named in a comma separated list and drop the ones named in another, in the
    % order they are listed above. Both lists are normally empty, which runs everything

    names = string(tests(:, 1));
    keep = true(size(names));

    wanted = testNames(selection, names, 'MAGTENSE_TESTS');
    if ~isempty(wanted)
        keep = keep & ismember(names, wanted);
    end

    dropped = testNames(skipped, names, 'MAGTENSE_SKIP');
    if ~isempty(dropped)
        keep = keep & ~ismember(names, dropped);
    end

    if all(keep)
        return
    end

    tests = tests(keep, :);
    fprintf('Running %d of the tests: %s\n', size(tests, 1), strjoin(tests(:, 1)', ', '));
end

function names = testNames(list, known, variable)
    % The names in a comma separated list, checked against the tests that exist

    names = strings(0, 1);
    if isempty(strtrim(list))
        return
    end

    names = strtrim(split(string(list), ','));
    names(names == "") = [];
    unknown = setdiff(names, known);
    if ~isempty(unknown)
        error('testMagTenseFunctions:unknownTest', ...
              '%s names tests that do not exist: %s', variable, ...
              strjoin(cellstr(unknown), ', '));
    end
end

function records = runAllTests(tests)
    % Run every test, catching errors so that one broken test does not hide the rest

    records = struct('name', {}, 'checks', {}, 'status', {}, 'elapsed', {}, 'message', {});

    for i = 1:size(tests, 1)
        name = tests{i, 1};
        fcn = tests{i, 2};

        fprintf('\n%s\n[%d/%d] %s\n%s\n', repmat('=', 1, 78), i, size(tests, 1), ...
                name, repmat('=', 1, 78));

        timer = tic;
        try
            checks = fcn();
            if isempty(checks)
                status = "ERROR";
                message = 'the test returned no checks';
            elseif all([checks.passed])
                status = "PASS";
                message = '';
            else
                status = "FAIL";
                message = '';
            end
        catch err
            checks = emptyChecks();
            status = "ERROR";
            message = err.message;
            fprintf(2, 'ERROR: %s\n', err.message);
        end
        elapsed = toc(timer);

        fprintf('--> %s: %s in %.1f s\n', name, status, elapsed);
        records(end + 1) = struct('name', string(name), 'checks', checks, ...
                                  'status', status, 'elapsed', elapsed, ...
                                  'message', string(message)); %#ok<AGROW>
    end
end

function checks = validationChecks(exampleDir, functionName, limit)
    % Run one of the magnetostatic validation functions from its own directory.
    % They all return the relative integrated error against FEM along x, y and z.

    oldDir = cd(exampleDir);
    cleanupObj = onCleanup(@() cd(oldDir));

    relativeError = feval(functionName);

    componentNames = {'Hx', 'Hy', 'Hz'};
    checks = emptyChecks();
    for i = 1:numel(relativeError)
        if i <= numel(componentNames)
            label = componentNames{i};
        else
            label = sprintf('component %d', i);
        end
        checks(end + 1) = makeCheck(sprintf('rel. integrated error, %s', label), ...
                                    relativeError(i), limit); %#ok<AGROW>
    end
end

function checks = micromagTestChecks(exampleDir, functionName, solverArgs)
    % Run one of the micromagnetic test functions in MagTense_tests from its own directory.
    % Each of them already returns the check struct array this suite expects, and each has a
    % python counterpart in python/examples/micromagnetism with the same limits.

    if nargin < 3
        solverArgs = {};
    end

    oldDir = cd(exampleDir);
    cleanupObj = onCleanup(@() cd(oldDir));

    checks = feval(functionName, solverArgs{:});
end

function checks = stdProblem4Checks(exampleDir, solverArgs)
    % Standard problem 4, first NIST field. The limit of 100 % is loose because <My>
    % and <Mz> average close to zero, so the relative measure has a small denominator.

    if nargin < 2
        solverArgs = {};
    end

    oldDir = cd(exampleDir);
    cleanupObj = onCleanup(@() cd(oldDir));

    if isempty(solverArgs)
        [~, ~, ~, ~, ~, ~, relativeError] = Standard_problem_4(1);
    else
        % The resolution has to be spelled out, since the options can only follow the
        % positional arguments. [36 9 1] is the default of Standard_problem_4
        [~, ~, ~, ~, ~, ~, relativeError] = Standard_problem_4(1, [36 9 1], solverArgs{:});
    end

    componentNames = {'<Mx>', '<My>', '<Mz>'};
    checks = emptyChecks();
    for i = 1:numel(relativeError)
        checks(end + 1) = makeCheck(sprintf('field 1: %s vs mumag', componentNames{i}), ...
                                    relativeError(i), 100); %#ok<AGROW>
    end
end

function checks = stdProblem6Checks(exampleDir, solverArgs)
    % Standard problem 6, comparing the depinning field with the analytical values for
    % the parameter variations, for the three sample orientations, and on an
    % unstructured mesh

    if nargin < 2
        solverArgs = {};
    end

    oldDir = cd(exampleDir);
    cleanupObj = onCleanup(@() cd(oldDir));

    limit = 5;                                       % Allowed deviation from theory [%]
    x_steps = 80;
    field_steps = 201;

    % Theoretical pinning fields [T] for the different parameter variations
    variations = {'akj', 'ak', 'aj', 'a', 'kj', 'k'};
    theory = [1.568, 1.089, 1.206, 0.838, 1.005, 0.565];

    checks = emptyChecks();

    % Parameter variations, all along x
    for i = 1:numel(variations)
        HP = Standard_problem_6(variations{i}, x_steps, field_steps, 'x', ...
                                'ShowTheResult', false, solverArgs{:});
        checks(end + 1) = makeCheck(sprintf('depinning field, variation "%s"', variations{i}), ...
                                    relativeErrorPercent(HP, theory(i)), limit); %#ok<AGROW>
    end

    % The same problem rotated onto each axis. The result must not depend on the
    % orientation, so this tests that the physics is implemented correctly in all three
    % directions
    directions = {'x', 'y', 'z'};
    for i = 1:numel(directions)
        HP = Standard_problem_6('akj', x_steps, field_steps, directions{i}, ...
                                'ShowTheResult', false, solverArgs{:});
        checks(end + 1) = makeCheck(sprintf('depinning field along %s', directions{i}), ...
                                    relativeErrorPercent(HP, theory(1)), limit); %#ok<AGROW>
    end

    % Unstructured mesh, which only works in the x direction
    HP = Standard_problem_6('akj', x_steps, field_steps, 'x', ...
                            'use_uniform_mesh', false, 'ShowTheResult', false, solverArgs{:});
    checks(end + 1) = makeCheck('depinning field, unstructured mesh', ...
                                relativeErrorPercent(HP, theory(1)), limit);
end

function value = relativeErrorPercent(measured, reference)
    % Deviation from the reference in percent. A missing result counts as infinitely bad
    if isempty(measured) || ~isscalar(measured) || ~isfinite(measured)
        value = Inf;
    else
        value = abs(measured - reference) / reference * 100;
    end
end

function check = makeCheck(name, value, limit)
    check = struct('check', name, 'value', double(value), 'limit', double(limit), ...
                   'passed', double(value) < double(limit));
end

function checks = emptyChecks()
    checks = struct('check', {}, 'value', {}, 'limit', {}, 'passed', {});
end


%% ----------------------------------------------------------------------------------
%  Reporting
%  ----------------------------------------------------------------------------------

function rows = buildRows(records)
    % Flatten the records into one row per check, as a table with string columns

    name = strings(0, 1);
    check = strings(0, 1);
    value = strings(0, 1);
    limit = strings(0, 1);
    status = strings(0, 1);
    margin = zeros(0, 1);                            % value/limit, NaN when not applicable

    for i = 1:numel(records)
        record = records(i);
        label = sprintf('%s [%.0f s]', record.name, record.elapsed);

        if isempty(record.checks)
            message = record.message;
            if strlength(message) > 60
                message = extractBefore(message, 61) + "...";
            end
            name(end + 1, 1) = label; %#ok<AGROW>
            check(end + 1, 1) = "raised: " + message; %#ok<AGROW>
            value(end + 1, 1) = "-"; %#ok<AGROW>
            limit(end + 1, 1) = "-"; %#ok<AGROW>
            status(end + 1, 1) = "ERROR"; %#ok<AGROW>
            margin(end + 1, 1) = NaN; %#ok<AGROW>
            continue
        end

        for j = 1:numel(record.checks)
            c = record.checks(j);
            if j == 1
                name(end + 1, 1) = label; %#ok<AGROW>
            else
                name(end + 1, 1) = ""; %#ok<AGROW>
            end
            check(end + 1, 1) = string(c.check); %#ok<AGROW>
            value(end + 1, 1) = formatNumber(c.value); %#ok<AGROW>
            limit(end + 1, 1) = formatNumber(c.limit); %#ok<AGROW>
            if c.passed
                status(end + 1, 1) = "PASS"; %#ok<AGROW>
            else
                status(end + 1, 1) = "FAIL"; %#ok<AGROW>
            end
            if c.limit > 0
                margin(end + 1, 1) = c.value / c.limit; %#ok<AGROW>
            else
                margin(end + 1, 1) = Inf; %#ok<AGROW>
            end
        end
    end

    rows = table(name, check, value, limit, status, margin);
end

function text = formatNumber(value)
    % Compact fixed or exponential notation, whichever reads better
    if isnan(value)
        text = "nan";
    elseif isinf(value)
        text = "inf";
    elseif value == 0
        text = "0";
    elseif abs(value) >= 1e-3 && abs(value) < 1e4
        text = string(sprintf('%.3g', value));
    else
        text = string(sprintf('%.2e', value));
    end
end

function printResultTable(records)
    % Print the overview as plain text in the command window

    rows = buildRows(records);
    header = ["Test", "Check", "Value", "Limit", "Status"];
    columns = {rows.name, rows.check, rows.value, rows.limit, rows.status};

    widths = zeros(1, numel(header));
    for c = 1:numel(header)
        widths(c) = max([strlength(header(c)); strlength(columns{c})]);
    end
    ruleWidth = sum(widths) + 3 * numel(widths);

    fprintf('\n%s\n', repmat('=', 1, ruleWidth));
    fprintf('OVERVIEW\n');
    fprintf('%s\n', repmat('=', 1, ruleWidth));
    for c = 1:numel(header)
        fprintf('%-*s   ', widths(c), header(c));
    end
    fprintf('\n%s\n', repmat('-', 1, ruleWidth));
    for r = 1:height(rows)
        for c = 1:numel(header)
            fprintf('%-*s   ', widths(c), columns{c}(r));
        end
        fprintf('\n');
    end
    fprintf('%s\n', repmat('-', 1, ruleWidth));

    [nPass, nFail, nError, totalTime] = summarise(records);
    fprintf('%d passed, %d failed, %d errored, in %.0f s\n', nPass, nFail, nError, totalTime);
end

function [nPass, nFail, nError, totalTime] = summarise(records)
    status = [records.status];
    nPass = sum(status == "PASS");
    nFail = sum(status == "FAIL");
    nError = sum(status == "ERROR");
    totalTime = sum([records.elapsed]);
end

function plotResultOverview(records, figurePath, saveFigure, resultsDir)
    % Colour coded table of every check, with the margin to its limit beside it

    rows = buildRows(records);
    nRows = height(rows);

    colours = struct('PASS', [0.788 0.906 0.788], ...
                     'FAIL', [0.949 0.722 0.710], ...
                     'ERROR', [0.949 0.722 0.710]);

    figureHeight = min(1400, 120 + 26 * (nRows + 1));
    fig = figure('Color', 'w', 'Position', [80 60 1500 figureHeight], ...
                 'Name', 'MagTense test suite');

    axTable = axes(fig, 'Position', [0.010 0.085 0.640 0.815]);
    axBar = axes(fig, 'Position', [0.685 0.085 0.295 0.815]);

    % ---- Table -------------------------------------------------------------------
    % Row 0 is the header and row r occupies y in [r, r+1], so the bar chart lines up
    % with it by simply using the same y limits
    hold(axTable, 'on')
    axis(axTable, 'off')
    xlim(axTable, [0 1])
    ylim(axTable, [0 nRows + 1])
    set(axTable, 'YDir', 'reverse')

    columnX = [0.015 0.300 0.735 0.825 0.915];
    header = {'Test', 'Check', 'Value', 'Limit', 'Status'};
    columns = {rows.name, rows.check, rows.value, rows.limit, rows.status};

    patch(axTable, [0 1 1 0], [0 0 1 1], [0.25 0.25 0.25], 'EdgeColor', [0.4 0.4 0.4]);
    for c = 1:numel(header)
        text(axTable, columnX(c), 0.5, header{c}, 'Color', 'w', 'FontWeight', 'bold', ...
             'FontSize', 8, 'VerticalAlignment', 'middle');
    end

    for r = 1:nRows
        rowColour = colours.(char(rows.status(r)));
        patch(axTable, [0 1 1 0], [r r r+1 r+1], rowColour, 'EdgeColor', [0.4 0.4 0.4]);
        for c = 1:numel(header)
            text(axTable, columnX(c), r + 0.5, columns{c}(r), 'FontSize', 8, ...
                 'VerticalAlignment', 'middle', 'Interpreter', 'none');
        end
    end

    % ---- Margin bars -------------------------------------------------------------
    % A bar is the measured value divided by its limit, so anything left of 1 passed
    barMin = 1e-10;
    barMax = 1e4;
    hold(axBar, 'on')
    for r = 1:nRows
        if isnan(rows.margin(r))
            continue                                 % Errored, there is nothing to draw
        end
        % Bars outside the plotted range are pulled just inside it to stay visible.
        % The table holds the exact numbers.
        width = min(max(rows.margin(r), 3 * barMin), 0.8 * barMax);
        if rows.status(r) == "PASS"
            barColour = [0.13 0.55 0.13];
        else
            barColour = [0.70 0.13 0.13];
        end
        patch(axBar, [barMin width width barMin], ...
              [r + 0.2, r + 0.2, r + 0.8, r + 0.8], barColour, 'EdgeColor', 'none');
    end
    set(axBar, 'XScale', 'log', 'YDir', 'reverse', 'YTick', [], 'FontSize', 8, ...
        'XGrid', 'on', 'GridLineStyle', ':', 'Box', 'on', ...
        'XMinorGrid', 'off', 'XMinorTick', 'off')
    xlim(axBar, [barMin barMax])
    ylim(axBar, [0 nRows + 1])
    xline(axBar, 1, 'k--', 'LineWidth', 1.2);
    text(axBar, 1, 0.5, ' limit', 'FontSize', 8, 'VerticalAlignment', 'middle');
    xlabel(axBar, 'measured value / acceptance limit', 'FontSize', 9)

    % ---- Headline ----------------------------------------------------------------
    [nPass, nFail, nError, totalTime] = summarise(records);
    if nFail == 0 && nError == 0
        verdict = 'ALL TESTS PASSED';
        titleColour = [0 0.39 0];
    else
        verdict = 'SOME TESTS FAILED';
        titleColour = [0.55 0 0];
    end
    annotation(fig, 'textbox', [0 0.915 1 0.075], 'String', ...
        {sprintf('MagTense MATLAB test suite - %s', verdict), ...
         sprintf('%d passed, %d failed, %d errored, %.0f s', nPass, nFail, nError, totalTime)}, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Color', titleColour, ...
        'EdgeColor', 'none');

    if saveFigure
        if ~isfolder(resultsDir)
            mkdir(resultsDir);
        end
        exportgraphics(fig, figurePath, 'Resolution', 200);
        fprintf('\nSaved overview figure to %s\n', figurePath);
    end
end
