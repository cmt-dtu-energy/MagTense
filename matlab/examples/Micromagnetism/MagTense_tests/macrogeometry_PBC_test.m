function checks = macrogeometry_PBC_test(options)
%MACROGEOMETRY_PBC_TEST
% Test if periodic boundary conditions as modeled by the macrogeometry method are implemented
% correctly in the MagTense micromagnetics solver. This is the MATLAB counterpart of
% python/examples/micromagnetism/macrogeometry_PBC_test.py and uses the same geometry, the same
% material parameters and the same acceptance limits.
%
% SHORT EXPLANATION : If the plotted polar angle is constant at the critical point, the code works.
%
% LONG EXPLANATION :
% The system is a pair of uniformly magnetised cubes of sidelength a, placed side by side along
% the axis of periodicity. There are periodic boundary conditions so that the cubes are repeated
% at regular intervals of A along that entire axis. Both cubes have a magnetocrystalline
% easy-axis perpendicular to it, but the magnetostatic interaction between the tiles and with
% their infinite copies produces an effective easy-axis along the axis of periodicity as well.
% The initial configuration is m1 = (cos pi/4, 0, sin pi/4) and m2 = (cos pi/4, 0, -sin pi/4).
% By symmetry theta_1 = -theta_2 = theta throughout, so the motion boils down to a single angle.
%
% Using the results of "Exact demagnetisation field for periodic one-dimensional array of
% rectangular prisms" by Durhuus et al., one can derive the parameter combination where the
% effective anisotropy is zero. At that critical point the moments stay put; for every other
% spacing they align along the axis of periodicity or anti-align along the easy-axis.
%
% The test runs three simulations at the critical spacing, one per axis of periodicity, and
% requires theta to stay put. Two further simulations at spacings on either side of the critical
% point act as controls: they must move theta towards pi/2 and 0 respectively, in the direction
% predicted by the analytical prefactor. Without those controls a solver that simply froze the
% magnetisation would pass the constant-theta check.
%
% Returns a struct array of checks with the fields 'check', 'value', 'limit' and 'passed',
% where a check passes when value < limit. That is the contract used by testMagTenseFunctions.m.

arguments
    options.ShowTheResult {mustBeNumericOrLogical} = true;   % Save the validation figures
    options.use_CUDA {mustBeNumericOrLogical} = false;
end

mu0 = 4*pi*1e-7;

addpath('../../../MEX_files');
addpath('../../../util');

%% Settings

% Distance between neighbouring domain copies in the macrogeometry
A_crit = 1.18369410e-07;    % Critical point where the torque is zero at all theta
A_wide = 2e-7;              % Should go towards theta = pi/2 (increasing theta)
A_narrow = 0.7e-7;          % Should go towards theta = 0 (decreasing theta)

testAxes = [1, 2, 3];       % Axes to test the PBC along (1 = x, 2 = y, 3 = z)
controlAxis = 1;            % Axis used for the two control simulations

% Timestepping
t_end = 1e-7;                            % Total time simulated [s]
t_step = 1e-10;                          % Requested output spacing [s]
nTimesteps = round(t_end/t_step) + 1;    % Include both time endpoints

% How constant theta has to be at the critical point to count as passed [rad]
angle_tol = 2e-2;
% The controls only approach their fixed point asymptotically, so they get a looser bound [rad]
control_tol = 0.1;

% Material properties
Ms = 1.0/mu0;               % Saturation magnetisation, 1 T [A/m]
alpha = 1;                  % Gilbert damping [-]
eta = alpha/(1 + alpha^2) * 2.21e5;     % Damping constant [m/(A*s)]
K = 2.7e4;                  % Uniaxial anisotropy constant [J/m^3]
Aex = 1e-20;                % Exchange constant [J/m]

a = 1e-8;                   % Side length of an individual micromagnetic cell [m]
n_copies = 1000;            % Copies of the domain on each side along the axis of periodicity

%% Analytical prefactor
% The rate of change of the potential energy with polar angle is dU/dtheta = prefactor *
% sin(theta)*cos(theta), so the sign decides whether the angle goes to pi/2 (negative) or 0
% (positive), and the magnitude decides how quickly. It vanishes at the critical spacing.
% MATLAB's psi(1,x) is the trigamma function, i.e. scipy's polygamma(1,x).
prefactor = @(q) 2*mu0*Ms^2/pi * ( ...
    2*atan(1/sqrt(3)) - atan(3/sqrt(11)) - atan(1/(3*sqrt(11))) ...          % Pair interaction
    + q.^2 .* (3*psi(1,1-q) - 3*psi(1,1+q)) ...                              % Tile and its copies
    + q.^2 .* (1/2*(psi(1,1+q) + psi(1,1-3*q) - psi(1,1+3*q) - psi(1,1-q))) ...  % Cross terms
    - 2*pi*K/(mu0*Ms^2) );                                                   % Anisotropy

disp('Analytical prefactor (>0 drives theta to 0, <0 drives theta to pi/2, 0 means no torque):')
for A = [A_crit, A_wide, A_narrow]
    fprintf('  A = %.4e m : %+.4e\n', A, prefactor(a/(2*A)));
end

%% Run the simulations
runLabels = {'critical', 'critical', 'critical', 'wide', 'narrow'};
runAxes = [testAxes, controlAxis, controlAxis];
runSpacings = [A_crit, A_crit, A_crit, A_wide, A_narrow];

results = cell(numel(runAxes), 1);
disp('Run simulations')
for r = 1:numel(runAxes)
    results{r} = run_case(runAxes(r), runSpacings(r), a, n_copies, Ms, Aex, K, eta, ...
                          t_end, nTimesteps, options.use_CUDA);
end
disp('Done running simulations')

%% Check the results
axis_names = {'x', 'y', 'z'};
checks = struct('check', {}, 'value', {}, 'limit', {}, 'passed', {});

% At the critical spacing theta must stay at its initial value of pi/4 along every axis
for r = 1:numel(testAxes)
    res = results{r};
    drift = max(abs(res.theta1 - pi/4));
    asymmetry = max(abs(res.theta1 - res.theta2));
    ok = (drift < angle_tol) && (asymmetry < angle_tol);
    if ok
        verdict = 'works';
    else
        verdict = 'FAILED';
    end
    fprintf(['Macrogeometry PBC %s along %s: max |theta - pi/4| = %.2e rad, ' ...
             'max |theta_1 + theta_2| = %.2e rad\n'], verdict, axis_names{testAxes(r)}, ...
             drift, asymmetry);
    checks(end+1) = struct('check', ...
        sprintf('PBC along %s: theta constant at critical A', axis_names{testAxes(r)}), ...
        'value', drift, 'limit', angle_tol, 'passed', drift < angle_tol); %#ok<AGROW>
    checks(end+1) = struct('check', ...
        sprintf('PBC along %s: theta_1 = -theta_2', axis_names{testAxes(r)}), ...
        'value', asymmetry, 'limit', angle_tol, 'passed', asymmetry < angle_tol); %#ok<AGROW>
end

% The controls must actually move, and in the direction the analytical prefactor predicts
controlTargets = [pi/2, 0];
controlNames = {'wide', 'narrow'};
controlSpacings = [A_wide, A_narrow];
for c = 1:2
    res = results{numel(testAxes) + c};
    target = controlTargets(c);
    moved = res.theta1(end) - res.theta1(1);
    expected = target - pi/4;
    distance = abs(res.theta1(end) - target);
    right_way = sign(moved) == sign(expected);
    ok = (distance < control_tol) && right_way;
    if ok
        verdict = 'works';
    else
        verdict = 'FAILED';
    end
    fprintf('Control at A = %.4e m %s: theta went from %.4f to %.4f rad, expected %.4f\n', ...
            controlSpacings(c), verdict, res.theta1(1), res.theta1(end), target);
    % A wrong direction is reported as a failure however close the endpoint happens to be
    if right_way
        value = distance;
    else
        value = Inf;
    end
    checks(end+1) = struct('check', ...
        sprintf('control, A %s: theta reaches %.3f rad', controlNames{c}, target), ...
        'value', value, 'limit', control_tol, 'passed', ok); %#ok<AGROW>
end

%% Plot the results
if options.ShowTheResult
    results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
    if ~isfolder(results_dir)
        mkdir(results_dir);
    end

    % ---- The critical spacing against the two controls, periodicity along x ----------------
    fig = figure('Visible', 'off', 'Color', 'w');
    ax = axes(fig); %#ok<LAXES>
    hold(ax, 'on')
    resNarrow = results{numel(testAxes) + 2};
    resWide = results{numel(testAxes) + 1};
    plot(ax, resNarrow.t*1e9, resNarrow.theta1, '--', 'Color', [1 0.55 0], 'LineWidth', 1.5, ...
         'DisplayName', 'A < A_{crit}, should go to 0');
    plot(ax, resWide.t*1e9, resWide.theta1, '--', 'Color', [0.5 0 0.5], 'LineWidth', 1.5, ...
         'DisplayName', 'A > A_{crit}, should go to \pi/2');
    plot(ax, results{controlAxis}.t*1e9, results{controlAxis}.theta1, 'Color', [0.13 0.55 0.13], ...
         'LineWidth', 1.5, 'DisplayName', 'A = A_{crit}, should be constant');
    xlabel(ax, 'Time [ns]')
    ylabel(ax, 'Polar angle')
    yticks(ax, [0 pi/4 pi/2])
    yticklabels(ax, {'0', '\pi/4', '\pi/2'})
    legend(ax, 'Location', 'best')
    box(ax, 'on')
    figure_path = fullfile(results_dir, 'macrogeometry_PBC_test_periodic_along_x.png');
    exportgraphics(fig, figure_path, 'Resolution', 200);
    close(fig);
    fprintf('Saved figure to %s\n', figure_path);

    % ---- The drift at the critical spacing along each axis, against the accepted band ------
    fig = figure('Visible', 'off', 'Color', 'w');
    ax = axes(fig); %#ok<LAXES>
    hold(ax, 'on')
    angle_window = 10^(floor(log10(angle_tol*10)));
    tPlot = results{1}.t*1e9;
    fill(ax, [tPlot(1) tPlot(end) tPlot(end) tPlot(1)], ...
         [-angle_tol -angle_tol angle_tol angle_tol], [0.5 0.5 0.5], ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Accepted interval');
    critical_colours = [0.13 0.55 0.13; 0.86 0.08 0.24; 0 0 0.5];
    for r = 1:numel(testAxes)
        res = results{r};
        plot(ax, res.t*1e9, res.theta1 - res.theta1(1), 'Color', critical_colours(r,:), ...
             'LineWidth', 1.5, 'DisplayName', sprintf('PBCs along %s', axis_names{testAxes(r)}));
    end
    xlabel(ax, 'Time [ns]')
    ylabel(ax, 'Change in polar angle')
    ylim(ax, [-angle_window angle_window])
    yticks(ax, [-angle_window 0 angle_window])
    legend(ax, 'Location', 'best')
    box(ax, 'on')
    figure_path = fullfile(results_dir, 'macrogeometry_PBC_test_periodic_along_x_y_or_z.png');
    exportgraphics(fig, figure_path, 'Resolution', 200);
    close(fig);
    fprintf('Saved figure to %s\n', figure_path);
end

if nargout == 0
    if all([checks.passed])
        disp('macrogeometry_PBC_test PASSED')
    else
        disp('macrogeometry_PBC_test FAILED')
    end
    clear checks
end
end


function res = run_case(testAxis, A, a, n_copies, Ms, Aex, K, eta, t_end, nTimesteps, use_CUDA)
% Simulate the two-cube system with periodic copies spaced A apart along testAxis.
%
% The angle is measured from the axis of periodicity towards the easy-axis, which is z when the
% periodicity is along x or y and x when the periodicity is along z. Both choices keep the
% easy-axis perpendicular to the axis of periodicity, so the three cases are rotations of one
% another.

if testAxis == 3
    perpAxis = 1;
else
    perpAxis = 3;
end

% Two tiles along the axis of periodicity, one along the others
resolution = [1 1 1];
resolution(testAxis) = 2;
grid_L = [a a a];
grid_L(testAxis) = 2*a;

problem = DefaultMicroMagProblem(resolution(1), resolution(2), resolution(3));
problem.grid_L = grid_L;
problem = problem.setUseCuda(use_CUDA);
problem = problem.setUseCVODE(false);
problem = problem.setUseDemag(true);
problem = problem.setMicroMagSolver('Dynamic');
% The analytical prefactor this test is built on does not hold for the averaged demagnetisation
% tensor: with useAvgN = 1 the moments rotate straight past the fixed point and both controls
% overshoot their targets. The python test opts out the same way.
problem.useAvgN = int32(0);

% The precession is switched off so the motion stays in a plane
problem.gamma = 0;
problem.alpha = eta;
problem.Ms = Ms*ones(2,1);
problem.A0 = Aex;
problem.K0 = K*ones(2,1);

% Easy-axis anisotropy perpendicular to the axis of periodicity
problem.u_ea = zeros(2,3);
problem.u_ea(:, perpAxis) = 1;

% Initial magnetisations at +-45 degrees to the axis of periodicity
theta_p = [pi/4; -pi/4];
problem.m0 = zeros(2,3);
problem.m0(:, testAxis) = cos(theta_p);
problem.m0(:, perpAxis) = sin(theta_p);

% Define the macrogeometry. It is a (2nx+1) X (2ny+1) X (2nz+1) grid of domain copies.
n_macro = int32([0 0 0]);
n_macro(testAxis) = int32(n_copies);
shiftVec = [0 0 0];
shiftVec(testAxis) = A;
problem.n_macro = n_macro;
problem.shiftVec = shiftVec;

% Zero applied field. Two samples completely define the constant zero field.
HextFct = @(t) (t>=0)' * [0, 0, 0];
problem = problem.setHext( HextFct, linspace(0, t_end, 2) );
problem = problem.setTime( linspace(0, t_end, nTimesteps) );

solution = struct();
prob_struct = struct(problem);
solution = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );

M_npv = squeeze(solution.M(:,:,1,:));       % (nt, 2, 3)
M1 = squeeze(M_npv(:,1,:));
M2 = squeeze(M_npv(:,2,:));

res.t = solution.t;
res.theta1 = atan2(M1(:,perpAxis), M1(:,testAxis));
% By symmetry theta_2 = -theta_1, so negating the perpendicular component of the second moment
% should reproduce the first angle
res.theta2 = atan2(-M2(:,perpAxis), M2(:,testAxis));
end
