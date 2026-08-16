function checks = temperature_test(options)
%TEMPERATURE_TEST
% Test the thermal fluctuation implementation. This is the MATLAB counterpart of
% python/examples/micromagnetism/temperature_test.py and uses the same geometry, the same
% material parameters and the same acceptance limits.
%
% The system is a collection of non-interacting, uniformly magnetised micromagnetic cells whose
% magnetic moments are initially vertical. There is no angular potential, but the temperature is
% finite, so over time the magnetic moment directions randomise until they are uniformly
% distributed on the unit sphere. The angular distribution as function of time is
%
%   P(t, theta) = sin(theta) * sum_{n=0}^inf (2n+1)/2 * Pn(cos theta) * exp(-n(n+1) D_LLG t)
%
% where Pn is the n'th order Legendre polynomial and D_LLG is a diffusion constant specific to
% the LLG equation. The script compares the simulated distribution to the analytical one.
%
% At a finite external field the moments do not randomise but settle into the Langevin
% distribution instead, and the mean magnetisation is compared with coth(xi) - 1/xi. Set the
% Hext option to switch to that mode.
%
% Returns a struct array of checks with the fields 'check', 'value', 'limit' and 'passed',
% where a check passes when value < limit. That is the contract used by testMagTenseFunctions.m.

arguments
    options.ShowTheResult {mustBeNumericOrLogical} = true;   % Save the validation figure
    options.use_CUDA {mustBeNumericOrLogical} = false;
    options.Ncell (1,1) {mustBeNumeric} = 12;   % Cells per direction; 30 gives better statistics
    options.Hext (1,1) {mustBeNumeric} = 0;     % External field along z [A/m]
end

addpath('../../../MEX_files');
addpath('../../../util');

%% Settings

% Constants of nature
mu0 = 4*pi*1e-7;            % Vacuum permeability [N/A^2]
gamma = 2.21e5;             % Gyromagnetic ratio [rad*m/(A*s)]
kB = 1.3806488e-23;         % Boltzmann's constant [J/K]

% Material properties
alpha = 0.1;                % Gilbert damping
T = 300;                    % Temperature [K]
Ms = 1.0/mu0;               % Saturation magnetisation, 1 T [A/m]
% The cells are meant to be independent, so the exchange is negligible rather than exactly zero:
% a zero exchange constant leaves the effective field identically zero and the ODE solver then
% aborts with a hard failure.
Aex = 1e-24;                % Exchange constant [J/m]
eta = alpha/(1 + alpha^2) * gamma;  % Damping constant [m/(A*s)]
Hext = options.Hext;

% Geometry
Ncell = options.Ncell;
NcellTot = Ncell^3;
a = 10e-9;                  % Single-cell sidelength [m]
V = a^3;

% For the analytical formula
nOrder = 20;                % Highest order of Legendre polynomial included
D_LLG = alpha * gamma * kB * T / (mu0 * (1+alpha^2) * Ms * V);   % Diffusion constant [1/s]

% Timestepping
t_end = 1.5/D_LLG;          % Total time interval simulated [s]
nt = 1001;                  % Number of requested output times

% Comparison and plot settings
Nbin = 20;                  % Number of histogram bins
Ndata = 8;                  % Number of times at which to compare simulation and theory

% Acceptance criteria. All are dominated by the Monte-Carlo noise of Ncell^3 samples.
% The first two only apply at Hext = 0, where the analytical solution is free diffusion on the
% sphere; the third only applies at Hext > 0, where the equilibrium is the Langevin distribution.
D_tol = [0.7, 1.4];         % Allowed range of D_simulated / D_LLG
TV_tol = 0.15;              % Allowed total variation distance between the distributions
% MagTense diffuses about 18 % faster than D_LLG, i.e. as if the temperature were 1.18*T, and the
% mean magnetisation in a field is low by the corresponding amount. The tolerance leaves room for
% that known offset plus the sampling noise.
Langevin_tol = 0.10;        % Allowed deviation of <mz> from the Langevin value

%% Setup the problem
problem = DefaultMicroMagProblem(Ncell, Ncell, Ncell);
problem.grid_L = [a*Ncell, a*Ncell, a*Ncell];
problem = problem.setUseCuda(options.use_CUDA);
problem = problem.setUseCVODE(false);
problem = problem.setUseDemag(false);       % The cells are meant to be non-interacting
problem = problem.setMicroMagSolver('Dynamic');

problem.gamma = gamma;
problem.alpha = eta;
problem.Ms = Ms*ones(NcellTot,1);
problem.A0 = Aex;
problem.K0 = zeros(NcellTot,1);
% The temperature has to be given per cell. The DefaultMicroMagProblem default is zeros(ntot),
% i.e. an ntot x ntot matrix, so it must be overwritten with a column vector here.
problem.temperature = T*ones(NcellTot,1);

% All moments start pointing along z
problem.m0 = repmat([0, 0, 1], NcellTot, 1);

% Constant applied field along z, zero by default
HextFct = @(t) (t>=0)' * [0, 0, Hext];
problem = problem.setHext( HextFct, linspace(0, t_end, 2) );
problem = problem.setTime( linspace(0, t_end, nt) );

%% Run the simulation
disp('Run simulation')
solution = struct();
prob_struct = struct(problem);
solution = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );
disp('Done running simulation')

t_n = solution.t;
M_npv = squeeze(solution.M(:,:,1,:));       % (nt, ntot, 3)

%% Analytical result
% Times at which simulation and theory are compared. The first entry is skipped in the
% comparison: at t = 0 the distribution is a delta function and the series does not converge.
compare_i = round(linspace(1, nt, Ndata));
t_cmp = t_n(compare_i);

theta_j = linspace(0, pi, 200);
x_j = cos(theta_j);

% Sum the Legendre series. The polynomials are built from the standard recurrence so that no
% toolbox is needed: (n+1) P_{n+1} = (2n+1) x P_n - n P_{n-1}
P_nj = zeros(Ndata, numel(theta_j));
Pnm1 = ones(size(x_j));     % P_0
Pn = x_j;                   % P_1
for k = 0:(nOrder-1)
    if k == 0
        Pk = Pnm1;
    elseif k == 1
        Pk = Pn;
    else
        Pk = ((2*(k-1)+1)*x_j.*Pn - (k-1)*Pnm1)/k;
        Pnm1 = Pn;
        Pn = Pk;
    end
    decay = exp(-k*(k+1)*D_LLG*t_cmp(:));
    P_nj = P_nj + (2*k+1)/2 * decay * Pk;
end
P_nj = P_nj .* repmat(sin(theta_j), Ndata, 1);

%% Compare analytical and simulated results
% MagTense returns the magnetisation normalised, but renormalise anyway so the polar angle is
% well defined even if the integrator lets |m| drift
M_norm_np = sqrt(sum(M_npv.^2, 3));
thetaSim_np = acos( min(max(M_npv(:,:,3)./M_norm_np, -1), 1) );

binEdges = linspace(0, pi, Nbin+1);
binWidth = mean(diff(binEdges));
thetaSim_nh = (binEdges(1:end-1) + binEdges(2:end))/2;

Psim_nh = zeros(Ndata, Nbin);
for n = 1:Ndata
    counts = histcounts(thetaSim_np(compare_i(n), :), binEdges);
    Psim_nh(n,:) = counts / NcellTot;       % Probability per bin
end

% Bin the analytical distribution the same way so the two can be compared bin by bin
Pana_nh = zeros(Ndata, Nbin);
for n = 1:Ndata
    for h = 1:Nbin
        sel = (theta_j >= binEdges(h)) & (theta_j <= binEdges(h+1));
        Pana_nh(n,h) = trapz(theta_j(sel), P_nj(n,sel));
    end
end

% Total variation distance between the two distributions, skipping the t = 0 comparison
TV_n = 0.5*sum(abs(Psim_nh - Pana_nh), 2);
TV_max = max(TV_n(2:end));

% The first Legendre order alone gives <cos theta> = exp(-2 D t), which is a much less noisy
% estimator of the diffusion constant than the full distribution
cosMean_n = mean(cos(thetaSim_np), 2);
fit_n = (t_n(:) > 0) & (t_n(:) < 0.7/D_LLG) & (cosMean_n > 0);
p = polyfit(t_n(fit_n), log(cosMean_n(fit_n)), 1);
D_sim = -p(1)/2;

fprintf('Number of independent cells        : %d\n', NcellTot);
fprintf('Analytical diffusion constant D_LLG: %.4e 1/s\n', D_LLG);
fprintf('Simulated diffusion constant       : %.4e 1/s  (ratio %.3f, accepted %.1f-%.1f)\n', ...
        D_sim, D_sim/D_LLG, D_tol(1), D_tol(2));
fprintf('Max total variation distance       : %.4f (accepted < %.2f)\n', TV_max, TV_tol);

Langevin_error = [];
if Hext > 0
    xi = mu0 * Ms * V * Hext / (kB * T);
    Langevin_mz = coth(xi) - 1/xi;
    equilibrated = t_n(:) > (2/3)*t_end;
    mzMean_n = mean(M_npv(:,:,3)./M_norm_np, 2);
    mz_sim = mean(mzMean_n(equilibrated));
    Langevin_error = abs(mz_sim - Langevin_mz);
    fprintf('Langevin parameter xi              : %.4f\n', xi);
    fprintf('Mean <mz>, simulated vs Langevin   : %.4f vs %.4f (deviation %.4f, accepted < %.2f)\n', ...
            mz_sim, Langevin_mz, Langevin_error, Langevin_tol);
end

%% Build the checks
if Hext > 0
    % With a field on, the zero-field analytical solution no longer applies
    checks = struct('check', 'mean magnetisation matches the Langevin value', ...
                    'value', Langevin_error, 'limit', Langevin_tol, ...
                    'passed', Langevin_error < Langevin_tol);
else
    D_ratio = D_sim/D_LLG;
    D_centre = mean(D_tol);
    D_halfwidth = (D_tol(2) - D_tol(1))/2;
    checks = struct('check', {'diffusion constant matches D_LLG', ...
                              'angular distribution matches theory'}, ...
                    'value', {abs(D_ratio - D_centre), TV_max}, ...
                    'limit', {D_halfwidth, TV_tol}, ...
                    'passed', {(D_ratio > D_tol(1)) && (D_ratio < D_tol(2)), TV_max < TV_tol});
end

%% Plot the result
if options.ShowTheResult
    results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
    if ~isfolder(results_dir)
        mkdir(results_dir);
    end

    fig = figure('Visible', 'off', 'Color', 'w');
    ax = axes(fig); %#ok<LAXES>
    hold(ax, 'on')
    cmap = colormap(ax, jet(256));
    % Skip the t = 0 comparison, which needs infinite orders of Legendre polynomials
    for n = 2:Ndata
        frac = t_cmp(n)/t_end;
        colour = cmap(max(1, round(frac*size(cmap,1))), :);
        plot(ax, theta_j, P_nj(n,:), '-', 'Color', colour, 'LineWidth', 1.2);
        plot(ax, thetaSim_nh, Psim_nh(n,:)/binWidth, 'o', 'Color', colour, ...
             'MarkerFaceColor', colour, 'MarkerSize', 4);
    end
    caxis(ax, [0 t_end*1e9]);
    cb = colorbar(ax);
    cb.Label.String = 'Time [ns]';
    xlabel(ax, 'Polar angle, \theta')
    ylabel(ax, 'Probability distribution, P(\theta)')
    xticks(ax, [0 pi/4 pi/2 3*pi/4 pi])
    xticklabels(ax, {'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'})
    box(ax, 'on')
    figure_path = fullfile(results_dir, 'temperature_test.png');
    exportgraphics(fig, figure_path, 'Resolution', 200);
    close(fig);
    fprintf('Saved figure to %s\n', figure_path);

    if Hext > 0
        fig = figure('Visible', 'off', 'Color', 'w');
        ax = axes(fig); %#ok<LAXES>
        hold(ax, 'on')
        yline(ax, Langevin_mz, 'Color', [0.5 0 0.5], 'LineWidth', 1.5, ...
              'DisplayName', 'Langevin equilibrium');
        plot(ax, t_n*1e9, mzMean_n, 'x', 'Color', [0 0 0.5], 'DisplayName', 'Calculation');
        xlabel(ax, 'Time, t [ns]')
        ylabel(ax, 'Average magnetisation, <M_z>/M_s')
        legend(ax, 'Location', 'best')
        box(ax, 'on')
        figure_path = fullfile(results_dir, 'temperature_test_langevin.png');
        exportgraphics(fig, figure_path, 'Resolution', 200);
        close(fig);
        fprintf('Saved figure to %s\n', figure_path);
    end
end

if nargout == 0
    if all([checks.passed])
        disp('temperature_test PASSED')
    else
        disp('temperature_test FAILED')
    end
    clear checks
end
end
