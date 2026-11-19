function checks = dipfmm_octree_test(options)
%DIPFMM_OCTREE_TEST
% dip-fmm on octree (unstructuredPrisms) meshes against MagTense's dense demagnetisation tensor.
% This is the MATLAB counterpart of
% python/examples/micromagnetism/MagTense_tests/dipfmm_octree_test.py and uses the same problems
% and the same acceptance limits.
%
% A 32 nm cube is split into 8 nm base tiles, and every tile within 3 nm of the plane z = 1.3 nm
% is refined twice (8 -> 4 -> 2 nm), giving a graded mesh whose tiles lie on the dyadic grid of
% the cube. With a random magnetisation the demagnetising field at t = 0 is compared between
%
%   * the dense MagTense tensor (the reference),
%   * dip-fmm on its uniform tree, two levels deep (leaves of 8 nm: every tile inside its leaf),
%   * dip-fmm on its adaptive tree (capacity 32, at most 5 levels, root fixed to the cube).
%
% A uniform 8^3 grid checks that the original uniform-grid path is unchanged and that the
% adaptive tree also runs there, a 10^3 grid and a grid of non-cubic cells check the root that
% is aligned with the cells, and the run without ReturnHall checks that the four field outputs
% are not allocated at full size.
%
% Differences from the python test:
%   * use_CUDA picks the one backend everything runs on, the dense reference included, as in the
%     other MATLAB tests (the CUDA and the NoCUDA solver are separate MEX files). The python test
%     runs the CPU and the GPU in one call against a CPU reference.
%   * The random magnetisation comes from MATLAB's generator, not numpy's, so the errors agree
%     with the python ones in size but not digit for digit.
%   * H_ext_applied is an output of the python interface only and is not checked.
%
% The MEX files have to be built with dip-fmm (buildMagTenseMEX('USE_FMM', true, ...)). One
% that is not uses the dense tensor for every run; the first check catches that.
%
% Returns a struct array of checks with the fields 'check', 'value', 'limit' and 'passed',
% where a check passes when value < limit. That is the contract used by testMagTenseFunctions.m.

arguments
    options.use_CUDA {mustBeNumericOrLogical} = false;
end

addpath('../../../MEX_files');
addpath('../../../util');

%% Settings

p = struct();
p.L = 32e-9;                    % Side length of the cube [m]
p.Ms = 8e5;                     % Saturation magnetisation [A/m]
p.A0 = 1.3e-11;                 % Exchange stiffness [J/m]
p.use_CUDA = options.use_CUDA;
base = 4;                       % Base tiles per side (8 nm)
levels = 2;                     % Refinements around the plane (8 -> 4 -> 2 nm)
plane_z = 1.3e-9;               % Refinement plane [m]
band = 3e-9;                    % Refine tiles reaching within this distance of the plane [m]
root_cube = [0, 0, 0, 0.5*p.L]; % The cube itself as root [centre x, y, z, half-width]
if p.use_CUDA
    tag = 'cuda';
else
    tag = 'cpu';
end

checks = struct('check', {}, 'value', {}, 'limit', {}, 'passed', {});

%% The octree mesh

[pts, abc] = octree_mesh(p.L, base, levels, plane_z, band);
n = size(pts, 1);
[tile_sizes, ~, index] = unique(round(abc(:,1)*1e9, 6));
fprintf('octree mesh: %d prisms, sizes [nm] / number:', n);
fprintf(' %g / %d', [tile_sizes'; accumarray(index, 1)']);
fprintf('\n');

m0 = random_m0(n, 7);
dense = demag_field(octree_problem(pts, abc, m0, p, {}));

orders = [6, 8];
limits = [2e-3, 5e-4];
for i = 1:numel(orders)
    uni = demag_field(octree_problem(pts, abc, m0, p, {'fmm_order', orders(i), 'fmm_depth', 2}));
    if i == 1
        % dip-fmm never reproduces the dense field to the last digit, a MEX without it does
        checks(end+1) = make_check('the MEX files are built with dip-fmm (0 = yes)', ...
                                   double(rel_l2(uni, dense) == 0), 0.5); %#ok<AGROW>
    end
    checks(end+1) = make_check(sprintf('octree %s: dip-fmm uniform tree (depth 2, order %d) vs dense', tag, orders(i)), ...
                               rel_l2(uni, dense), limits(i)); %#ok<AGROW>
    ada = demag_field(octree_problem(pts, abc, m0, p, {'fmm_order', orders(i), 'fmm_tree', 1, ...
                                     'fmm_cap', 32, 'fmm_maxdep', 5, 'fmm_root', root_cube}));
    checks(end+1) = make_check(sprintf('octree %s: dip-fmm adaptive tree (capacity 32, order %d) vs dense', tag, orders(i)), ...
                               rel_l2(ada, dense), limits(i)); %#ok<AGROW>
end
% No fmm_root: MagTense chooses a root that keeps the tiles in their leaves
ada_auto = demag_field(octree_problem(pts, abc, m0, p, {'fmm_order', 6, 'fmm_tree', 1}));
checks(end+1) = make_check(sprintf('octree %s: dip-fmm adaptive tree, root chosen by MagTense, order 6 vs dense', tag), ...
                           rel_l2(ada_auto, dense), 5e-3);

%% Uniform grid: the original path and the adaptive tree on equal cubes

m0u = random_m0(512, 11);
dense_u = demag_field(grid_problem([8 8 8], m0u, p, {}));
uni_u = demag_field(grid_problem([8 8 8], m0u, p, {'fmm_order', 6, 'fmm_depth', 2}));
checks(end+1) = make_check(sprintf('uniform %s: dip-fmm uniform tree vs dense', tag), rel_l2(uni_u, dense_u), 2e-3);
ada_u = demag_field(grid_problem([8 8 8], m0u, p, {'fmm_order', 6, 'fmm_tree', 1, 'fmm_root', root_cube}));
checks(end+1) = make_check(sprintf('uniform %s: dip-fmm adaptive tree vs dense', tag), rel_l2(ada_u, dense_u), 2e-3);
% The original uniform-cuboid plan (bounding-box root) is still reachable
leg = demag_field(grid_problem([8 8 8], m0u, p, {'fmm_order', 6, 'fmm_depth', 2, 'fmm_align', 0}));
checks(end+1) = make_check(sprintf('uniform %s: original uniform plan (fmm_align = 0) vs dense', tag), rel_l2(leg, dense_u), 2e-3);

%% Aligned roots
% A 10^3 grid at depth 2 has 2.5 cells per leaf on its bounding box, so half the cells stick out
% of their leaves; the aligned root (3-cell leaves) keeps them inside

m10 = random_m0(1000, 13);
dense10 = demag_field(grid_problem([10 10 10], m10, p, {}));
al = demag_field(grid_problem([10 10 10], m10, p, {'fmm_order', 6, 'fmm_depth', 2}));
un = demag_field(grid_problem([10 10 10], m10, p, {'fmm_order', 6, 'fmm_depth', 2, 'fmm_align', 0}));
e_al = rel_l2(al, dense10);
e_un = rel_l2(un, dense10);
fprintf('  10^3 grid, depth 2 (%s): aligned %.3e, bounding-box root %.3e\n', tag, e_al, e_un);
checks(end+1) = make_check(sprintf('10^3 %s: uniform tree on the aligned root (3 cells per leaf) vs dense', tag), e_al, 2e-3);
checks(end+1) = make_check(sprintf('10^3 %s: aligned root is more accurate than the bounding box (ratio)', tag), e_al/e_un, 1.0);

% Cells that are not cubes: 10 x 10 x 5 cells of 3.2 x 3.2 x 6.4 nm; the root unit is 6.4 nm
m_aniso = random_m0(500, 17);
dense_an = demag_field(grid_problem([10 10 5], m_aniso, p, {}));
an = demag_field(grid_problem([10 10 5], m_aniso, p, {'fmm_order', 6, 'fmm_depth', 2}));
checks(end+1) = make_check('10x10x5 non-cubic cells: aligned uniform tree vs dense', rel_l2(an, dense_an), 2e-3);
pa = grid_problem([10 10 5], m_aniso, p, {'fmm_depth', 2});
root = pa.alignedDipFmmRoot();
leaf = 2*root(4)/4;
whole = abs(leaf/3.2e-9 - round(leaf/3.2e-9)) < 1e-9 && abs(leaf/6.4e-9 - round(leaf/6.4e-9)) < 1e-9;
checks(end+1) = make_check('10x10x5: aligned leaf holds whole cells along x and z (0 = yes)', double(~whole), 0.5);

% The root sizes of the 80^3 production grid (2.5 nm cells)
p80 = DefaultMicroMagProblem(80, 80, 80);
p80.grid_L = 200e-9*[1, 1, 1];
p80.fmm_tree = int32(0);
p80.fmm_depth = int32(5);
root = p80.alignedDipFmmRoot();
checks(end+1) = make_check('80^3: aligned uniform root at depth 5 is 240 nm (|side - 240 nm| [nm])', abs(2*root(4) - 240e-9)*1e9, 1e-6);
p80.fmm_tree = int32(1);
root = p80.alignedDipFmmRoot();
checks(end+1) = make_check('80^3: aligned adaptive root is 320 nm (|side - 320 nm| [nm])', abs(2*root(4) - 320e-9)*1e9, 1e-6);

%% Outputs without ReturnHall: small field arrays, the magnetisation still in full
% Two applied fields, each relaxed by the energy minimizer

problem = octree_problem(pts, abc, m0, p, {'fmm_order', 6, 'fmm_depth', 2});
problem.ReturnHall = int32(0);
problem = problem.setMicroMagSolver('Explicit');
problem = problem.setTime( linspace(0, 1e-10, 2) );        % LL fallback window of the minimizer
h_applied = [1e4, -2e4, 3e4];
problem = problem.setHext( @(t) (t(:) + 1).*h_applied, [0, 1] );
solution = problem.MagTenseLandauLifshitzSolver_mex( struct(problem), struct() );
sizes_ok = numel(solution.H_exc) == 3 && numel(solution.H_ext) == 3 && ...
           numel(solution.H_dem) == 3 && numel(solution.H_ani) == 3;
checks(end+1) = make_check('field outputs are (1,1,1,3) without ReturnHall (0 = yes)', double(~sizes_ok), 0.5);
m_end = reshape(solution.M(end,:,end,:), [], 3);
checks(end+1) = make_check('magnetisation still returned in full (|m| - 1)', max(abs(vecnorm(m_end, 2, 2) - 1)), 1e-6);

%% Report

for i = 1:numel(checks)
    if checks(i).passed
        verdict = 'PASS';
    else
        verdict = 'FAIL';
    end
    fprintf('  [%s] %s: %.3e (limit %.0e)\n', verdict, checks(i).check, checks(i).value, checks(i).limit);
end
fprintf('\n%d/%d checks passed\n', sum([checks.passed]), numel(checks));

if nargout == 0
    if all([checks.passed])
        disp('dipfmm_octree_test PASSED')
    else
        disp('dipfmm_octree_test FAILED')
    end
    clear checks
end
end


function [centres, abc] = octree_mesh(L, base, levels, plane_z, band)
% Tile centres and full side lengths of the graded cube, centred on the origin.

h0 = L/base;
g = -0.5*L + ((0:base-1) + 0.5)*h0;
[X, Y, Z] = ndgrid(g, g, g);
centres = [X(:), Y(:), Z(:)];
sizes = h0*ones(size(centres,1), 1);
[sx, sy, sz] = ndgrid([-0.25, 0.25]);
children = [sx(:), sy(:), sz(:)];                           % Offsets of the eight children

kept_c = zeros(0, 3);
kept_s = zeros(0, 1);
for level = 0:levels
    refine = abs(centres(:,3) - plane_z) - 0.5*sizes <= band;
    if level == levels
        refine(:) = false;
    end
    kept_c = [kept_c; centres(~refine,:)]; %#ok<AGROW>
    kept_s = [kept_s; sizes(~refine)]; %#ok<AGROW>
    parents = centres(refine,:);
    parent_sizes = sizes(refine);
    centres = repelem(parents, 8, 1) + repmat(children, size(parents,1), 1).*repelem(parent_sizes, 8, 1);
    sizes = repelem(parent_sizes/2, 8, 1);
end
centres = kept_c;
abc = repmat(kept_s, 1, 3);
end


function m0 = random_m0(n, seed)
% Random unit vectors, reproducible from the seed.

stream = RandStream('mt19937ar', 'Seed', seed);
m0 = randn(stream, n, 3);
m0 = m0./vecnorm(m0, 2, 2);
end


function problem = octree_problem(pts, abc, m0, p, fmm)
% The graded mesh as a grid of unstructured prisms.

problem = DefaultMicroMagProblem(size(pts,1), 1, 1);
problem = problem.setMicroMagGridType('unstructuredPrisms');
problem.grid_pts = pts;
problem.grid_abc = abc;
problem = common_settings(problem, m0, p, fmm);
end


function problem = grid_problem(res, m0, p, fmm)
% A uniform grid of the 32 nm cube with res cells (any shape).

problem = DefaultMicroMagProblem(res(1), res(2), res(3));
problem = common_settings(problem, m0, p, fmm);
end


function problem = common_settings(problem, m0, p, fmm)
% Material, solver and dip-fmm settings shared by every problem. fmm is a cell array of
% fmm_* property names and values; an empty one leaves the dense tensor in use.

ntot = problem.ntot;
problem.grid_L = [p.L, p.L, p.L];
problem = problem.setUseCuda(p.use_CUDA);
problem = problem.setUseCVODE(false);
problem = problem.setMicroMagSolver('Dynamic');
problem.gamma = 0;
problem.alpha = 4000;
problem.Ms = p.Ms*ones(ntot,1);
problem.A0 = p.A0*ones(ntot,1);
problem.K0 = zeros(ntot,1);
problem.m0 = m0;
problem.ReturnHall = int32(1);
problem.window_ena = int32(0);

problem = problem.setUseDipFmm(~isempty(fmm));
for i = 1:2:numel(fmm)
    if strcmp(fmm{i}, 'fmm_root')
        problem.fmm_root = fmm{i+1};
    else
        problem.(fmm{i}) = int32(fmm{i+1});
    end
end
end


function H = demag_field(problem)
% The demagnetising field at t = 0, ntot x 3 [A/m].

t = linspace(0, 1e-13, 2);
problem = problem.setTime(t);
problem = problem.setHext( @(t) zeros(numel(t), 3), t );
solution = problem.MagTenseLandauLifshitzSolver_mex( struct(problem), struct() );
H = reshape(solution.H_dem(1,:,1,:), [], 3);
end


function value = rel_l2(a, b)
value = norm(a(:) - b(:))/norm(b(:));
end


function c = make_check(name, value, limit)
c = struct('check', name, 'value', value, 'limit', limit, 'passed', isfinite(value) && value < limit);
end
