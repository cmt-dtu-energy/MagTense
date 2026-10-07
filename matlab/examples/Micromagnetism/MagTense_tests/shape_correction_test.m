function checks = shape_correction_test(options)
%SHAPE_CORRECTION_TEST
% Test if the shape correction field is implemented correctly in the MagTense micromagnetics
% solver. This is the MATLAB counterpart of
% python/examples/micromagnetism/MagTense_tests/shape_correction_test.py and uses the same geometry, the same
% material parameters and the same acceptance limits.
%
% The system is a uniform grid of micromagnetic tiles forming a rectangular prism. The
% macrogeometry is the simulated domain itself while the sample is a much longer prism, so the
% demagnetisation field of the sample is uniform across the domain. The uniaxial anisotropy is
% set to cancel that field exactly, which means the magnetisation should neither rotate nor
% break up: the polar angle stays where it starts and the magnetisation stays uniform.
%
% A second check uses an unstructured prism mesh running from 0 to L whose first layer of cells
% is refined, so the cells at the two ends differ in size. The macrogeometry prism of the shape
% correction has to be centred on the outer faces of that mesh. For a uniform magnetisation the
% field of the cells and that of the macrogeometry prism then cancel exactly, and what is left is
% the field of the sample. Centring it on the midpoint of the end-cell centres instead shifts it
% by a/8 and leaves a field error of about 20 %.
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

a = 1e-9;                       % Side lengths of individual micromagnetic cell [m]
resolution = [2, 3, 12];        % Number of cells along x, y and z
Ntiles = prod(resolution);
grid_L = resolution*a;          % Size of the simulated domain [m]

% Sample- and macrogeometry shapes (assumes rectangular prisms). The sample is intentionally
% macroscopic so the demagnetisation field becomes uniform across the simulated domain.
macroShape = grid_L;
aSample = 3;                    % Length of sample along x
bSample = 1;                    % Transverse dimension of sample
sampleShape = [aSample, bSample, bSample];

% Timestepping
t_end = 1e-8;                            % Total time simulated [s]
t_step = 5e-11;                          % Requested output spacing [s]
nTimesteps = round(t_end/t_step) + 1;    % Include both time endpoints

% How far the magnetisation may deviate from uniform to count as passed [-]
uniformity_tol = 1e-4;
% How far the mean polar angle may drift over the simulation to count as passed [rad]
drift_tol = 1e-5;
% How far the field on the mesh with unequal end cells may deviate from that of the sample,
% relative to it. Single precision, in which the tensor is stored, sets the floor.
mesh_field_tol = 1e-5;

% Material properties
Ms = 0.1/mu0;                   % Saturation magnetisation [A/m]
alpha = 1;                      % Gilbert damping [-]
eta = alpha/(1 + alpha^2) * 2.21e5;     % Damping constant [m/(A*s)]
% A moderate exchange retains a uniform state without making the ODE solver stiff on the 1 nm
% mesh. It also avoids the zero-exchange 0/0 normalisation case.
Aex = 1e-14;                    % Exchange constant [J/m]

%% Shape dependent anisotropy
% The shape dependent, uniaxial anisotropy constant from the whole sample's demagnetisation
% field. Assumes equal length along y and z, i.e. b = c.
fx = bSample^2 / (aSample * sqrt(aSample^2 + 2*bSample^2));
Nxx = 2/pi * atan(fx);
Ksh = mu0*Ms^2/4 * (1 - 3*Nxx);

%% Setup the problem
problem = DefaultMicroMagProblem(resolution(1), resolution(2), resolution(3));
problem.grid_L = grid_L;
problem = problem.setUseCuda(options.use_CUDA);
problem = problem.setUseCVODE(false);
problem = problem.setUseDemag(true);
problem = problem.setMicroMagSolver('Dynamic');
% The uniaxial anisotropy below is set from the analytical demagnetisation factor of the sample
% prism, which is what the tensor evaluated at the cell centres reproduces. With the averaged
% tensor (the default) the cancellation is only good to 1.1e-4 rad instead of 2.6e-8, so this
% test opts out explicitly to keep the tight drift limit. The python test does the same.
problem.useAvgN = int32(0);

% The precession is switched off so the motion stays in the xz-plane
problem.gamma = 0;
problem.alpha = eta;
problem.Ms = Ms*ones(Ntiles,1);
problem.A0 = Aex;
problem.K0 = -Ksh*ones(Ntiles,1);   % Uniaxial anisotropy equated to the shape anisotropy
problem.u_ea = repmat([1, 0, 0], Ntiles, 1);

% The macrogeometry and the sample it is cut out of
problem.macroShape = macroShape;
problem.sampleShape = sampleShape;

% Initial magnetisation in the xz-plane at 45 degrees to horizontal
m0_v = [1/sqrt(2), 0, 1/sqrt(2)];
problem.m0 = repmat(m0_v/norm(m0_v), Ntiles, 1);

% Zero applied field. Two samples completely define the constant zero field.
HextFct = @(t) (t>=0)' * [0, 0, 0];
problem = problem.setHext( HextFct, linspace(0, t_end, 2) );
problem = problem.setTime( linspace(0, t_end, nTimesteps) );

%% Run the simulation
disp('Run simulation')
solution = struct();
prob_struct = struct(problem);
solution = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );
disp('Done running simulation')

%% Check the results
t_n = solution.t;
M_npv = squeeze(solution.M(:,:,1,:));       % (nt, ntot, 3)

% Is the magnetisation still uniform at the end of the simulation?
M_end = squeeze(M_npv(end,:,:));            % (ntot, 3)
Mavg = mean(M_end, 1);
uniformity = max( sqrt(sum((M_end - Mavg).^2, 2)) ) / Ms;
if uniformity < uniformity_tol
    disp('Magnetisation is uniform')
else
    disp('Magnetisation is NOT uniform')
end

% The shape correction cancels the demagnetisation field of the sample against the uniaxial
% anisotropy, so the mean polar angle should not move over the simulation at all
Mavg_nv = squeeze(mean(M_npv, 2));          % (nt, 3)
theta_n = atan2(Mavg_nv(:,3), Mavg_nv(:,1));
drift = max(abs(theta_n - theta_n(1)));
if drift < drift_tol
    fprintf('Mean polar angle stays put: max change = %.3e rad\n', drift);
else
    fprintf('Mean polar angle DRIFTS: max change = %.3e rad\n', drift);
end

%% The unstructured mesh with unequal end cells
mesh_deviation = mesh_check(a, Ms, Aex, eta, m0_v/norm(m0_v), sampleShape, Nxx, options.use_CUDA);
if mesh_deviation < mesh_field_tol
    verdict = 'is the sample field';
else
    verdict = 'is NOT the sample field';
end
fprintf('Unstructured mesh with unequal end cells: the field %s, max deviation = %.3e\n', ...
        verdict, mesh_deviation);

checks = struct('check', {'magnetisation stays uniform', ...
                          'mean polar angle does not drift', ...
                          'unstructured mesh, unequal end cells: field is the sample field'}, ...
                'value', {uniformity, drift, mesh_deviation}, ...
                'limit', {uniformity_tol, drift_tol, mesh_field_tol}, ...
                'passed', {uniformity < uniformity_tol, drift < drift_tol, ...
                           mesh_deviation < mesh_field_tol});

%% Plot the result
if options.ShowTheResult
    fig = figure('Visible', 'off', 'Color', 'w');
    ax = axes(fig); %#ok<LAXES>
    hold(ax, 'on')
    tPlot_n = t_n*1e9;      % Switch from s to ns
    angle_window = 10^(floor(log10(drift_tol*10)));
    fill(ax, [tPlot_n(1) tPlot_n(end) tPlot_n(end) tPlot_n(1)], ...
         [-drift_tol -drift_tol drift_tol drift_tol], [0.13 0.55 0.13], ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Accepted interval');
    plot(ax, tPlot_n, theta_n - theta_n(1), 'Color', [0.86 0.08 0.24], 'LineWidth', 1.5, ...
         'DisplayName', 'Calculation (should be constant)');
    xlabel(ax, 'Time [ns]')
    ylabel(ax, 'Change in polar angle')
    ylim(ax, [-angle_window angle_window])
    yticks(ax, [-angle_window 0 angle_window])
    legend(ax, 'Location', 'best')
    box(ax, 'on')

    results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
    if ~isfolder(results_dir)
        mkdir(results_dir);
    end
    figure_path = fullfile(results_dir, 'shape_correction_test.png');
    exportgraphics(fig, figure_path, 'Resolution', 200);
    close(fig);
    fprintf('Saved figure to %s\n', figure_path);
end

if nargout == 0
    if all([checks.passed])
        disp('shape_correction_test PASSED')
    else
        disp('shape_correction_test FAILED')
    end
    clear checks
end
end


function deviation = mesh_check(a, Ms, Aex, eta, m0_v, sampleShape, Nxx, use_CUDA)
% Relative deviation of the demagnetisation field from that of the sample, on a mesh whose end
% cells differ in size.
%
% The mesh is a 4 x 2 x 2 block of cubes of side a running from 0 to L, whose first layer along x
% is refined into cells of half the size, so the centres of the end cells are a/4 and a/2 from the
% two faces. For a uniform magnetisation the field of the cells and of the macrogeometry prism
% cancel exactly, by superposition, wherever the prism sits on the domain, provided it sits on the
% domain. What is left is the field of the sample, which is uniform across the nanometre domain.

base = [4 2 2];
grid_L = base*a;
pts = zeros(0,3);
abc = zeros(0,3);
for k = 0:base(3)-1
    for j = 0:base(2)-1
        for i = 0:base(1)-1
            corner = [i j k]*a;
            if i == 0
                [sx, sy, sz] = ndgrid(0:1, 0:1, 0:1);
                sub = [sx(:) sy(:) sz(:)];
                pts = [pts; corner + (sub + 0.5)*a/2]; %#ok<AGROW>
                abc = [abc; repmat(a/2*[1 1 1], 8, 1)]; %#ok<AGROW>
            else
                pts = [pts; corner + 0.5*a]; %#ok<AGROW>
                abc = [abc; a*[1 1 1]]; %#ok<AGROW>
            end
        end
    end
end
ntot = size(pts,1);

problem = DefaultMicroMagProblem(ntot, 1, 1);
problem = problem.setMicroMagGridType('unstructuredPrisms');
problem.grid_pts = pts;
problem.grid_abc = abc;
problem.grid_L = grid_L;
problem = problem.setUseCuda(use_CUDA);
problem = problem.setUseCVODE(false);
problem = problem.setUseDemag(true);
problem = problem.setMicroMagSolver('Dynamic');
% The shape correction is evaluated at the cell centres, so the cells have to be as well for the
% field of the cells and of the macrogeometry prism to cancel
problem.useAvgN = int32(0);
problem.use_fmm = int32(0);
problem.ReturnHall = int32(1);

problem.gamma = 0;
problem.alpha = eta;
problem.Ms = Ms*ones(ntot,1);
problem.A0 = Aex*ones(ntot,1);
problem.K0 = zeros(ntot,1);
problem.m0 = repmat(m0_v, ntot, 1);
problem.macroShape = grid_L;
problem.sampleShape = sampleShape;

% Only the field of the initial state is needed
HextFct = @(t) (t>=0)' * [0, 0, 0];
problem = problem.setHext( HextFct, linspace(0, 1e-15, 2) );
problem = problem.setTime( linspace(0, 1e-15, 2) );

solution = struct();
prob_struct = struct(problem);
solution = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );

H_dem = reshape(solution.H_dem(1,:,1,:), ntot, 3);
% The sample prism has b = c, so its tensor at the centre is diagonal with Nyy = Nzz
Nyy = (1 - Nxx)/2;
H_sample = -Ms * (diag([Nxx Nyy Nyy]) * m0_v(:))';
deviation = max(abs(H_dem - H_sample), [], 'all') / max(abs(H_sample));
end
