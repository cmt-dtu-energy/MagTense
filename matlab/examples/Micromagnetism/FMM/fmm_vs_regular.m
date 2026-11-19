function [results, M_final] = fmm_vs_regular( grid_size, options )
%FMM_VS_REGULAR
%Compare the dip-fmm demagnetization field with the regular dense demagnetization tensor
%
%Syntax:
%------
%   fmm_vs_regular()
%   [results, M_final] = fmm_vs_regular( grid_size, options )
%
%Description:
%------
%   A cube of grid_size^3 cells of 5 nm is relaxed from a slightly tilted uniform state,
%   once with the regular dense demagnetization tensor and once for every requested
%   expansion order and tree depth with dip-fmm. The final states are compared; the
%   relative RMS difference has to fall as the expansion order goes up. This is the
%   MATLAB counterpart of python/examples/micromagnetism/FMM/fmm_vs_regular.py.
%
%   The MEX files have to be built with dip-fmm: the libraries by the root Makefile or a
%   Visual Studio configuration without "_no_FMM", and the MEX files with
%   buildMagTenseMEX('USE_FMM', true, ...), see 'Demag field - FMM' in the TechManual.
%   Otherwise MagTense reports that use_fmm is ignored, both runs use the dense tensor
%   and every difference is zero.
%
%Input arguments:
%------
%   grid_size : Integer (Default value is 15)
%       The number of cells along each side of the cube.
%
%Options:
%-------
%   orders : Array of integers - Default is [2 6 10]
%       The dip-fmm expansion orders to run.
%
%   depths : Array of integers - Default is [2 3]
%       The depths of the dip-fmm uniform tree to run.
%
%   use_CUDA : Interpreted as a logical - Default is true
%       Use CUDA, for the regular tensor as well as for dip-fmm.
%
%   t_end : Double - Default is 40e-9
%       The time to relax for [s].
%
%Output arguments:
%-------
%   results : Array of size (number of runs) x 3
%       One row per dip-fmm run: expansion order, tree depth and the relative RMS
%       difference of the final magnetization from the regular result.
%
%   M_final : Struct
%       The final magnetization, ntot x 3: regular for the dense tensor and fmm, a cell
%       array with one entry per row of results, for dip-fmm.

arguments
    grid_size (1,1) {mustBeInteger,mustBePositive}      = 15
    options.orders (1,:) {mustBeInteger,mustBePositive} = [2 6 10]
    options.depths (1,:) {mustBeInteger,mustBePositive} = [2 3]
    options.use_CUDA {mustBeNumericOrLogical}           = true
    options.t_end (1,1) double                          = 40e-9
end

addpath('../../../MEX_files');
addpath('../../../util');

fprintf('\n=== %d^3 regular reference ===\n', grid_size);
M_regular = run_problem( grid_size, false, 1, 2, options );

results = zeros(numel(options.orders)*numel(options.depths), 3);
M_final = struct('regular', M_regular, 'fmm', {cell(size(results,1), 1)});
k = 0;
for order = options.orders
    for depth = options.depths
        fprintf('\n=== %d^3, order=%d, depth=%d ===\n', grid_size, order, depth);
        M_fmm = run_problem( grid_size, true, order, depth, options );
        k = k + 1;
        results(k,:) = [order, depth, norm(M_fmm(:) - M_regular(:))/norm(M_regular(:))];
        M_final.fmm{k} = M_fmm;
        fprintf('relative RMS: %.3e\n', results(k,3));
    end
end

fprintf('\nGrid     N  Order  Depth  Relative RMS\n');
for k = 1:size(results,1)
    fprintf('%2d^3  %5d  %5d  %5d  %.3e\n', grid_size, grid_size^3, results(k,1), results(k,2), results(k,3));
end

end

function M_end = run_problem( grid_size, use_dip_fmm, order, depth, options )
%Relax one cubic, non-periodic problem and return the final magnetization, ntot x 3

cell_size = 5e-9;
ntot = grid_size^3;

problem = DefaultMicroMagProblem(grid_size, grid_size, grid_size);
problem.grid_L = grid_size*cell_size*[1,1,1];

problem = problem.setMicroMagDemagApproximation('none');
problem = problem.setUseCuda( options.use_CUDA );
problem = problem.setUseCVODE( false );
problem = problem.setMicroMagSolver('Dynamic');

% Material properties; gamma = 0 is pure relaxation, which keeps the run inexpensive
problem.alpha = 4.42e3;
problem.gamma = 0;
problem.Ms = 8e5*ones(ntot,1);
problem.K0 = zeros(ntot,1);
problem.A0 = 1.3e-11*ones(ntot,1);
problem.tol = 1e-3;

% Initial magnetization: along x, tilted a little out of the plane across the cells
m0 = [ones(ntot,1), zeros(ntot,1), linspace(-0.1,0.1,ntot)'];
problem.m0 = m0./vecnorm(m0,2,2);

% No applied field
problem = problem.setTime( linspace(0,options.t_end,201) );
problem = problem.setHext( @(t) zeros(numel(t),3), linspace(0,options.t_end,2) );

% The demagnetization field: dip-fmm or the regular dense tensor
problem = problem.setUseDipFmm( use_dip_fmm );
problem.fmm_order = int32(order);
problem.fmm_depth = int32(depth);

solution = problem.MagTenseLandauLifshitzSolver_mex( struct(problem), struct() );
M_end = squeeze(solution.M(end,:,:));

end
