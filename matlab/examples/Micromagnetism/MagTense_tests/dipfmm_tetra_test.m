function checks = dipfmm_tetra_test(options)
%DIPFMM_TETRA_TEST
% dip-fmm on tetrahedral meshes against MagTense's dense tetrahedral demagnetisation tensor.
% This is the MATLAB counterpart of
% python/examples/micromagnetism/MagTense_tests/dipfmm_tetra_test.py and uses the same problem
% and the same acceptance limits.
%
% A 32 nm cube is cut into 8 nm cubes and every cube into six tetrahedra (384 elements). With a
% random magnetisation the demagnetising field at t = 0 is compared between
%
%   * the dense MagTense tensor of exact volume-averaged tetrahedron pairs (the reference),
%   * dip-fmm on its uniform tree at depth 1 (16 nm leaves) and 2 (8 nm leaves, one cube per leaf),
%   * dip-fmm on its adaptive tree (capacity 6, at most 3 levels, root fixed to the cube).
%
% A uniformly magnetised mesh checks the exact volume average <H> = -M/3 of a cube however it is
% cut.
%
% dip-fmm assigns elements to leaves by their centroids. When the leaves are much smaller than
% the elements, elements reach far out of their leaves and the far field is no longer accurate
% (dip-fmm warns about this); choose the depth so that a leaf is not much smaller than an element.
%
% Differences from the python test:
%   * use_CUDA picks the one backend everything runs on, the dense reference included, as in the
%     other MATLAB tests (the CUDA and the NoCUDA solver are separate MEX files). The python test
%     runs the CPU and the GPU in one call against a CPU reference.
%   * The random magnetisation comes from MATLAB's generator, not numpy's, so the errors agree
%     with the python ones in size but not digit for digit.
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
h = 8e-9;                       % Cubes of the Kuhn split [m]
root_cube = [0, 0, 0, 0.5*p.L]; % The cube itself as root [centre x, y, z, half-width]
if p.use_CUDA
    tag = 'cuda';
else
    tag = 'cpu';
end

checks = struct('check', {}, 'value', {}, 'limit', {}, 'passed', {});

%% The mesh

[nodes, elements] = kuhn_tetra_mesh(round(p.L/h), h);
vol = tetra_volumes(nodes, elements);
n = size(elements, 2);
fprintf('tetrahedral mesh: %d elements, volume %.6f L^3\n', n, sum(vol)/p.L^3);

%% Random magnetisation

m0 = random_m0(n, 7);
dense = demag_field(tetra_problem(nodes, elements, m0, p, {}));

orders = [6, 8];
limits = [2e-3, 1e-3];
for depth = 1:2
    for i = 1:numel(orders)
        fmm = demag_field(tetra_problem(nodes, elements, m0, p, {'fmm_order', orders(i), 'fmm_depth', depth}));
        if depth == 1 && i == 1
            % dip-fmm never reproduces the dense field to the last digit, a MEX without it does
            checks(end+1) = make_check('the MEX files are built with dip-fmm (0 = yes)', ...
                                       double(rel_l2(fmm, dense) == 0), 0.5); %#ok<AGROW>
        end
        checks(end+1) = make_check(sprintf('%s: dip-fmm uniform tree (depth %d, order %d) vs dense', tag, depth, orders(i)), ...
                                   rel_l2(fmm, dense), limits(i)); %#ok<AGROW>
    end
end
ada = demag_field(tetra_problem(nodes, elements, m0, p, {'fmm_order', 6, 'fmm_tree', 1, 'fmm_cap', 6, ...
                                'fmm_maxdep', 3, 'fmm_root', root_cube}));
checks(end+1) = make_check(sprintf('%s: dip-fmm adaptive tree (capacity 6, root = cube, order 6) vs dense', tag), ...
                           rel_l2(ada, dense), 2e-3);

%% Uniform magnetisation: the volume-averaged field of the cube is exactly -M/3

mu = repmat([0, 0, 1], n, 1);
hu = demag_field(tetra_problem(nodes, elements, mu, p, {'fmm_order', 6, 'fmm_depth', 2}));
avg = sum(vol.*hu, 1)/sum(vol)/p.Ms;
checks(end+1) = make_check(sprintf('%s: uniform M gives <H>/Ms = (0, 0, -1/3) (max deviation)', tag), ...
                           max(abs(avg - [0, 0, -1/3])), 1e-3);

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
        disp('dipfmm_tetra_test PASSED')
    else
        disp('dipfmm_tetra_test FAILED')
    end
    clear checks
end
end


function m0 = random_m0(n, seed)
% Random unit vectors, reproducible from the seed.

stream = RandStream('mt19937ar', 'Seed', seed);
m0 = randn(stream, n, 3);
m0 = m0./vecnorm(m0, 2, 2);
end


function problem = tetra_problem(nodes, elements, m0, p, fmm)
% The tetrahedral mesh with the material, solver and dip-fmm settings of the test. fmm is a cell
% array of fmm_* property names and values; an empty one leaves the dense tensor in use.

problem = DefaultMicroMagProblem(size(elements,2), 1, 1);
problem = problem.setMicroMagGridTetrahedron(nodes, elements);
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


function vol = tetra_volumes(nodes, elements)
% The volumes of the tetrahedra given by the 3 x M nodes and the 1-based 4 x N connectivity,
% as an N x 1 array.

a = nodes(:, elements(2,:)) - nodes(:, elements(1,:));
b = nodes(:, elements(3,:)) - nodes(:, elements(1,:));
c = nodes(:, elements(4,:)) - nodes(:, elements(1,:));
vol = abs(dot(a, cross(b, c, 1), 1))'/6;
end


function [nodes, elements] = kuhn_tetra_mesh(n, a)
% A cube of n x n x n cubes of side a, centred on the origin, each split into six tetrahedra by
% the Kuhn subdivision, as create_tetra_mesh in python/src/magtense/utils.py does. Returns the
% nodes as a 3 x M array and the 1-based connectivity as a 4 x N array.

kuhn = [1 2 4 8; 1 2 6 8; 1 3 4 8; 1 3 7 8; 1 5 6 8; 1 5 7 8];
[ii, jj, kk] = ndgrid(0:n, 0:n, 0:n);
nodes = [ii(:) jj(:) kk(:)]' * a - n*a/2;
nid = reshape(1:size(nodes,2), [n+1, n+1, n+1]);

elements = zeros(4, 6*n^3);
e = 0;
for ic = 1:n
    for jc = 1:n
        for kc = 1:n
            corners = zeros(1,8);
            m = 0;
            for k = 0:1
                for j = 0:1
                    for i = 0:1
                        m = m + 1;
                        corners(m) = nid(ic+i, jc+j, kc+k);
                    end
                end
            end
            for t = 1:6
                e = e + 1;
                elements(:,e) = corners(kuhn(t,:))';
            end
        end
    end
end
end


function value = rel_l2(a, b)
value = norm(a(:) - b(:))/norm(b(:));
end


function c = make_check(name, value, limit)
c = struct('check', name, 'value', value, 'limit', limit, 'passed', isfinite(value) && value < limit);
end
