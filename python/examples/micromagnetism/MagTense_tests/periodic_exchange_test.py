"""Test the implementation of periodic exchange coupling on both grid types of MagTense.

MagTense builds the exchange operator in two entirely different ways: from a stencil on a
uniform grid, and from a mesh analysis on an unstructured mesh of prisms
(UnstructuredMeshAnalysis.f90 and ComputeDifferentialOperatorsFromMesh.f90). Both have to
honour the periodic boundary conditions requested through exchPBC, so both are tested here, with
three independent methods: an end to end simulation, the assembled exchange operator, and a
comparison with the operator on a supercell of copies of the mesh. The two methods that do not
need a regular lattice are then repeated on an irregular mesh of grains.

Method 1, a sectioned simulation
--------------------------------
A 6 X 6 X 6 cube with periodic boundary conditions is divided into 8 sections that are
internally exchange coupled but independent of one another. The test is passed if the
magnetisation vectors are identical within each connected section but distinct between
sections due to random initial conditions.

The 8 sections are a 4 X 4 X 4 cube of center tiles, the 8 corners, the xFace (4 tiles on both
the left and right side), the yFace (4 tiles on both the front and back face), the zFace
(4 tiles on both the top and bottom faces), the xyEdge (8 tiles on the edges perpendicular to
the xy plane), the xzEdge and the yzEdge. The remaining cells are numerically inert spacers and
there is no dipole interaction, so the separate sections are effectively decoupled, giving 8
tests of the exchange coupling. The center tiles test basic exchange, while the rest check the
periodic boundary conditions. Only a section that is connected through a periodic boundary can
become uniform, so this is an end to end test through the actual time integration.

The measure is the largest deviation of a moment from the mean of its section, in units where
|m| = 1, which is how MagTense returns the magnetisation. A section that is not coupled through
the boundary keeps its random initial directions and lands at a spread of order 1, so the
measure separates a working implementation from a broken one by orders of magnitude.

Method 2, the exchange operator itself
--------------------------------------
A regular lattice of identical cubes is used, so that the assembled exchange matrix can be
compared with an exact result on both grid types. Two properties are tested, neither of which
depends on the details of the discretisation:

  1. Zero row sums. A uniform magnetisation has to give a vanishing exchange field.

  2. Plane waves are eigenvectors. On a periodic lattice the exchange operator is invariant
     under a translation by one cell, and a translation invariant operator has the plane waves
     exp(i*k*x), k = 2*pi*m/L, as exact eigenvectors. The residual |A*phi - lambda*phi| is
     therefore at round-off level if, and only if, every tile - including the ones at the
     boundaries - is coupled in the same way. This is the sharp test of the periodic linking.

The eigenvalues are also compared with the Laplacian dispersion relation, lambda(k) -> -k^2 for
k -> 0, so that the operator is not just periodic but is still the correct second derivative,
and the operators are computed without periodic boundaries as well to show what the linking
changes. Those two are reported rather than tested.

Method 3, the supercell reference
---------------------------------
The plane waves of method 2 are eigenvectors only because the lattice is regular, so that test
cannot be carried over to an irregular mesh. The same question - is every tile coupled the way
the periodic linking claims - can be asked without assuming anything about the mesh, by
replicating the mesh instead of the magnetisation: the exchange operator is assembled on a
supercell of copies of the mesh with free boundaries, and the rows of the central copy are
folded back onto the tiles of the original mesh. Every tile of the central copy has exactly the
neighbourhood that the periodic operator gives the corresponding tile, so the two matrices have
to agree entry by entry. This is stricter than the plane wave test, which only looks at one
eigenvector, and it is run on the regular prism mesh as well as on the irregular mesh below.

The irregular mesh
------------------
Both unstructured tests above use a regular lattice of identical cubes, which is the easy case
for the mesh analysis. They are therefore repeated on an irregular mesh of Voronoi grains of the
kind used in python/experiments/Grain_perf, where the cells differ in size and a face at one end
of the domain is linked to several smaller faces at the other end. Method 1 does not carry over,
since its spacer pattern is built for a 6 X 6 X 6 grid, but its idea does: a slab of cells one
base cell thick is left out of the mesh, so that what remains falls into two halves which touch
one another only through the periodic boundary. Leaving the cells out is cleaner than the inert
spacers of method 1, since there is then nothing in the gap to pin the moments next to it, and
the halves relax to one common direction if, and only if, the linking works. The same simulation
is run without exchPBC as a control, which is what shows that the gap really does cut the mesh
in two.

NOTE: the Fortran library has to be rebuilt for the unstructured part of this test to run.

Running the file executes the test and saves two figures, one of the measures and one
of the irregular mesh and the states it relaxes to. ``run_test()`` returns the same result as
a list of checks, which is the contract the combined suite in testMagTenseFunctions.py expects.

Usage:
    python periodic_exchange_test.py                  # 6 x 6 x 6 cells, PBC along x, y and z
    python periodic_exchange_test.py --res 8 4 4 --pbc 1 0 0
    python periodic_exchange_test.py --no-grains     # Skip the irregular mesh
"""

# General modules
import argparse
import itertools
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

# Magtense stuff
from magtense.micromag import MicromagProblem

# Style settings for plots
plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = False

#%% Settings

output_dir = Path(__file__).resolve().parent

# The two grid types that MagTense assembles the exchange operator for
GRID_TYPES = ('uniform', 'unstructuredPrisms')
GRID_LABELS = {'uniform': 'uniform grid', 'unstructuredPrisms': 'unstructured mesh'}

# The irregular mesh. A base grid that is refined once wherever a grain boundary passes through a
# cell, which is the mesh type used in python/experiments/Grain_perf. Two sizes are needed: the
# sectioned simulation wants a few cells on either side of the gap, while the supercell
# reference replicates the mesh 27 times and has to stay small enough to assemble quickly.
GRAIN_LABEL = 'grain mesh'
GRAIN_COUNT = 6           # Number of grains
GRAIN_SEED = 11           # Seed of the grain positions, so that the mesh is reproducible
GRAIN_REFINEMENTS = 1     # Number of times a cell at a grain boundary is split into eight
GRAIN_OFFSET = 0.10       # Half thickness of the refined layer at a boundary, in base cells
GRAIN_RES_SECTION = 5     # Base cells per side, sectioned simulation (odd, see grain_gap_mask)
GRAIN_RES_OPERATOR = 3    # Base cells per side, exchange operator test

# ---- Acceptance thresholds --------------------------------------------------------------
# Spread of the unit magnetisation vectors within a periodically connected section. A section
# that is not coupled at all keeps its random initial directions and has a spread of order 1,
# so both limits discriminate by more than two orders of magnitude.
# The unstructured mesh gets the looser limit because its stencil reaches every tile sharing a
# vertex with a face, so the section tiles touch the frozen spacers and relax to a spread of
# about 3e-3 rather than to zero. That is an equilibrium floor, not slow convergence: it is
# unchanged when the integration time is raised by two orders of magnitude.
SECTION_TOL = {'uniform': 1e-4, 'unstructuredPrisms': 1e-2}
ROW_SUM_TOL = 1e-8      # Relative to the largest coupling in the matrix
EIGENVECTOR_TOL = 1e-8  # Relative residual of the plane wave eigenvectors
SUPERCELL_TOL = 1e-10   # Deviation from the supercell operator, relative to the largest coupling
# The two halves of the split grain mesh, in units where |m| = 1: the distance between their mean
# magnetisations, and the spread of the moments over both of them. Both land near 1e-4 when the
# halves are linked through the periodic boundary and at order 1 when they are not, so the limit
# discriminates by two orders of magnitude in either direction. Unlike the sections above there
# is no floor from neighbouring spacers here, because the gap holds no cells at all.
SPLIT_TOL = 1e-2
SPLIT_CONTROL_MIN = 1e-1

# Number of exchange matrix entries per tile that MagTense reserves for the returned matrix.
# The stencil of the unstructured operator reaches all the tiles sharing a vertex with a face,
# so it is far wider than the seven point stencil of the uniform grid. MagTense silently
# returns zeros if this is too small, which is caught below.
EXCH_PRESIZE = 64

#%% Fixed settings (don't change these)

# Constants of nature
mu0 = 4*np.pi * 1e-7    # Vacuum permeability [N/A^2]
gamma = 2.21*1e5    # Gyromagnetic ratio [m/(A*s)]

# Material properties
Ms = 1.0               # Saturation magnetisation in magnetic region [T]
Ms = Ms/mu0            # [A/m]
alpha = 1              # Gilbert damping [-]
eta = alpha/(1 + alpha**2) * gamma  # Damping constant [m/(A*s)]
K = 0                  # Uniaxial anisotropy constant [J/m^3]
T = 0                  # Temperature [K]
Aex = 1e-11            # Exchange constant in the coupled sections [J/m]
Aex_spacer = 1e-25     # Negligible exchange constant in the spacer cells [J/m]

# Applied field (set to zero)
def fct_h_ext(t) -> np.ndarray:
    return np.zeros([np.size(t), 3])

# Timestepping for the sectioned simulation
t_end = 1e-8                      # Total time simulated [s]
t_step = 1e-10                    # Requested output spacing [s]
nTimesteps = int(round(t_end / t_step)) + 1  # Include both time endpoints

# Side lengths of individual micromagnetic cells
a_section = 1e-8    # Sectioned simulation [m]
a_operator = 1e-9   # Exchange operator test [m]

# The sectioned geometry is fixed: the spacer pattern below only works for 6 cells per side
SECTION_SIDE = 6

sectionNames_s = ['center', 'corner', 'xFace', 'yFace', 'zFace', 'xyEdge', 'xzEdge', 'yzEdge']


#%% Geometry

def random_directions(ntot, seed=7):
    """Reproducible random unit vectors, used as the initial state of every simulation here."""
    rng = np.random.default_rng(seed)
    m0_pv = rng.uniform(-1, 1, [ntot, 3])
    return m0_pv / np.sqrt(np.sum(m0_pv**2, axis=1, keepdims=True))


def build_sections():
    """Return the boolean masks of the 8 periodically connected sections and the spacers.

    Tiles are ordered the way setupGrid in LandauLifshitzEquationSolver.f90 orders them, with x
    running fastest: n = i + 6*j + 36*k. The sectioned geometry happens to be symmetric under
    exchanging x and z, so getting this backwards swaps the xFace/zFace and xyEdge/yzEdge labels
    without changing which eight sections are tested.
    """
    side = SECTION_SIDE
    ntot = side**3

    corner_n, center_n = np.zeros(ntot, dtype=bool), np.zeros(ntot, dtype=bool)
    xFace_n, yFace_n, zFace_n = (np.zeros(ntot, dtype=bool), np.zeros(ntot, dtype=bool),
                                 np.zeros(ntot, dtype=bool))
    xyEdge_n, xzEdge_n, yzEdge_n = (np.zeros(ntot, dtype=bool), np.zeros(ntot, dtype=bool),
                                    np.zeros(ntot, dtype=bool))
    n = 0
    for k in range(side):
        for j in range(side):
            for i in range(side):
                if i in {0, 5}:
                    if j in {0, 5}:
                        if k in {0, 5}:
                            corner_n[n] = True
                        elif k in {2, 3}:
                            xyEdge_n[n] = True
                    elif j in {2, 3}:
                        if k in {0, 5}:
                            xzEdge_n[n] = True
                        elif k in {2, 3}:
                            xFace_n[n] = True
                elif i in {2, 3}:
                    if j in {0, 5}:
                        if k in {0, 5}:
                            yzEdge_n[n] = True
                        elif k in {2, 3}:
                            yFace_n[n] = True
                    elif j in {2, 3}:
                        if k in {0, 5}:
                            zFace_n[n] = True
                        elif k in {2, 3}:
                            center_n[n] = True
                n += 1

    sections_s = [center_n, corner_n, xFace_n, yFace_n, zFace_n, xyEdge_n, xzEdge_n, yzEdge_n]

    # Identify the coupled magnetic sections and the numerically inert spacers.
    nonZeroMagnetisation_n = np.logical_or.reduce(sections_s)
    zeroMagnetisation_n = np.logical_not(nonZeroMagnetisation_n)

    return sections_s, zeroMagnetisation_n


def build_prism_mesh(res, a):
    """Return the centers and the side lengths of a regular lattice of cubes.

    The tiles are ordered with x running fastest, i.e. the same ordering as the uniform grid,
    so that the two exchange matrices can be compared entry by entry.
    """
    nx, ny, nz = res
    L = np.array([nx * a, ny * a, nz * a])

    pts = np.zeros((nx * ny * nz, 3))
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                n = i + nx * j + nx * ny * k
                pts[n, :] = (np.array([i, j, k]) + 0.5) * a - L / 2

    abc = a * np.ones((nx * ny * nz, 3))

    return pts, abc


def build_section_mesh():
    """The sectioned 6 X 6 X 6 geometry as an unstructured prism mesh.

    Uses the tile ordering of build_sections, so the same section masks apply to both grids.
    """
    side = SECTION_SIDE
    ntot = side**3
    pts = np.zeros((ntot, 3))
    n = 0
    for k in range(side):
        for j in range(side):
            for i in range(side):
                pts[n, :] = (np.array([i, j, k]) + 0.5) * a_section - side * a_section / 2
                n += 1
    abc = a_section * np.ones((ntot, 3))
    return pts, abc


def build_grain_mesh(res_base, a_base, n_grains=GRAIN_COUNT, n_refine=GRAIN_REFINEMENTS,
                     seed=GRAIN_SEED):
    """Return the centers, the side lengths, the size and the grains of an irregular mesh.

    A base grid of res_base**3 cubes is refined recursively: a cell is split into eight whenever
    it is close enough to the boundary between two grains that the boundary can pass through it.
    The grains are the Voronoi cells of a set of random points, so the result is the mesh type
    used in python/experiments/Grain_perf, reproduced here without the grain bookkeeping and
    without any dependency beyond numpy. The mesh is centered on the origin and is deterministic
    for a given seed.
    """
    L = res_base * a_base
    rng = np.random.default_rng(seed)
    # The grains are kept away from the faces of the domain, so that the refined cells at a grain
    # boundary do not all end up sitting on the periodic boundaries
    seeds_gd = rng.uniform(0.15, 0.85, (n_grains, 3)) * L
    offset = GRAIN_OFFSET * a_base

    def boundary_distance(center_d):
        """Distance from a point to the plane bisecting the two grains that are closest to it."""
        distance_g = np.linalg.norm(seeds_gd - center_d, axis=1)
        first, second = np.argsort(distance_g)[:2]
        return ((distance_g[second]**2 - distance_g[first]**2)
                / (2 * np.linalg.norm(seeds_gd[first] - seeds_gd[second])))

    pts, abc = [], []

    def refine(center_d, size_d, level):
        # The corner of a cell reaches half a diagonal beyond its center, so a cell that is closer
        # to a grain boundary than that is cut by the boundary and is refined
        if level < n_refine and boundary_distance(center_d) <= offset + np.linalg.norm(size_d) / 2:
            for corner in itertools.product((-0.25, 0.25), repeat=3):
                refine(center_d + np.asarray(corner) * size_d, size_d / 2, level + 1)
        else:
            pts.append(center_d)
            abc.append(size_d)

    for k in range(res_base):
        for j in range(res_base):
            for i in range(res_base):
                refine((np.array([i, j, k]) + 0.5) * a_base, np.full(3, a_base), 0)

    pts, abc = np.array(pts), np.array(abc)
    # The grain of a cell is the seed it is closest to, which is what makes the grains the Voronoi
    # cells of the seeds. It is not used by the simulation, only to draw the mesh
    grain_n = np.argmin(np.linalg.norm(pts[:, np.newaxis, :] - seeds_gd[np.newaxis, :, :],
                                       axis=2), axis=1)

    return pts - L / 2, abc, np.array([L, L, L]), grain_n


def grain_gap_mask(pts, a_base):
    """The middle layer of base cells of a grain mesh, which is left out to split it in two.

    The mesh is centered on the origin, so with an odd number of base cells per side the middle
    layer is the one spanning -a_base/2 to +a_base/2. A refined cell lies entirely within its
    base cell, so testing the centers picks out whole cells and leaves a gap that is one base
    cell wide, i.e. wider than any cell in the mesh. The exchange stencil reaches the cells that
    share a vertex with a face and therefore cannot bridge it.
    """
    return np.abs(pts[:, 0]) < a_base / 2


#%% Method 1: the sectioned simulation

def section_deviations(grid_type, cuda=False, cvode=False):
    """Simulate the sectioned cube and return the spread of the moments in each section.

    A section that is exchange coupled all the way around the periodic boundary relaxes to a
    single common direction, so the spread collapses to the tolerance of the time integration.
    """
    side = SECTION_SIDE
    ntot = side**3
    sections_s, zeroMagnetisation_n = build_sections()

    # Ms must remain nonzero because MagTense forms exchange and anisotropy factors
    # by dividing by Ms.  The tiny spacer exchange decouples the eight sections
    # without introducing 0/0 coefficients in the uniform exchange stencil.
    Ms_n = Ms * np.ones(ntot)
    Aex_n = Aex * np.ones((ntot, 1))
    Aex_n[zeroMagnetisation_n, 0] = Aex_spacer

    m0_pv = random_directions(ntot)

    if grid_type == 'uniform':
        problem_res = [side, side, side]
        grid_pts, grid_abc = None, None
    else:
        problem_res = (ntot, 1, 1)
        grid_pts, grid_abc = build_section_mesh()

    problem = MicromagProblem(
        res=problem_res,
        grid_type=grid_type,
        grid_pts=grid_pts,
        grid_abc=grid_abc,
        grid_L=[side*a_section, side*a_section, side*a_section],
        A0=Aex_n,
        Ms=Ms_n[:, np.newaxis],
        K0=K,
        alpha=eta,
        gamma=gamma,
        m0=m0_pv,
        exchPBC=np.array([1, 1, 1], dtype=np.int32),
        cuda=cuda,
        cvode=cvode,
        usereturnhall=True,
        solver='dynamic',
        T=T,
        usedemag=0,
        exch_presize=EXCH_PRESIZE,
    )

    print(f'Run simulation on the {GRID_LABELS[grid_type]}')
    result = problem.run_simulation(
        t_end=t_end,
        nt=nTimesteps,
        fct_h_ext=fct_h_ext,
        nt_h_ext=2,  # Two samples completely define the constant zero field.
    )
    print('Done running simulation')

    # Test if exchange coupled magnetic moments are indeed identical. MagTense returns the
    # magnetisation already normalised to |m| = 1, so the spread is measured in those units
    # and must not be divided by Ms.
    M_pv = result[1][-1, :, 0, :]
    tolerance = SECTION_TOL[grid_type]
    deviations = []
    for section_n, sectionName in zip(sections_s, sectionNames_s, strict=True):
        Mavg_dv = np.mean(M_pv[section_n, :], axis=0, keepdims=True)
        error_p = np.sqrt(np.sum((M_pv[section_n, :] - Mavg_dv)**2, axis=1))
        max_error = float(np.max(error_p))
        deviations.append(max_error)
        verdict = 'works' if max_error < tolerance else 'FAILED'
        print(f'  Exchange coupling {verdict} for {sectionName} tiles: {max_error:.3e}')

    return deviations


def split_deviations(pts, abc, grid_L, pbc, cuda=False, cvode=False):
    """Simulate a mesh that a gap has split in two and return how far from uniform it ends up.

    The two halves touch one another only through the periodic boundary, so they relax to a
    common direction if the linking works and to unrelated directions if it does not. Two
    measures are returned: the distance between the mean magnetisation of the two halves, which
    is what the linking decides, and the spread of the moments over the whole mesh, which is the
    same measure that the sections above use. Both are in the units of |m| = 1 that MagTense
    returns and both are of order 1 when the halves are not coupled. The magnetisation itself is
    returned as well, so that the figure can show what the two measures are made of.
    """
    ntot = len(pts)
    m0_pv = random_directions(ntot)

    problem = MicromagProblem(
        res=(ntot, 1, 1),
        grid_type='unstructuredPrisms',
        grid_pts=pts,
        grid_abc=abc,
        grid_L=list(grid_L),
        A0=Aex * np.ones((ntot, 1)),
        Ms=Ms * np.ones((ntot, 1)),
        K0=K,
        alpha=eta,
        gamma=gamma,
        m0=m0_pv,
        exchPBC=np.asarray(pbc, dtype=np.int32),
        cuda=cuda,
        cvode=cvode,
        usereturnhall=True,
        solver='dynamic',
        T=T,
        usedemag=0,
        exch_presize=EXCH_PRESIZE,
    )

    print(f'Run simulation on the {GRAIN_LABEL}, exchPBC = {list(pbc)}')
    result = problem.run_simulation(
        t_end=t_end,
        nt=nTimesteps,
        fct_h_ext=fct_h_ext,
        nt_h_ext=2,  # Two samples completely define the constant zero field.
    )
    print('Done running simulation')

    M_pv = result[1][-1, :, 0, :]
    left_n = pts[:, 0] < 0
    right_n = pts[:, 0] > 0

    difference = np.linalg.norm(np.mean(M_pv[left_n, :], axis=0)
                                - np.mean(M_pv[right_n, :], axis=0))
    Mavg_dv = np.mean(M_pv, axis=0, keepdims=True)
    spread = np.max(np.sqrt(np.sum((M_pv - Mavg_dv)**2, axis=1)))

    return float(difference), float(spread), M_pv


#%% Method 2: the exchange operator

def exchange_triplets(problem_res, grid_type, grid_pts, grid_abc, grid_L, pbc, cuda=False):
    """Run a minimal simulation and return the assembled exchange matrix in COO form.

    Only the setup of the problem is of interest here, so the simulation is run for a couple of
    trivial timesteps with the demagnetisation field and the external field switched off. The
    triplets are returned rather than a dense array because the supercell reference assembles
    matrices with up to 27 times as many rows as the one they are compared with.
    """
    ntot = int(np.prod(problem_res))

    problem = MicromagProblem(
        res=problem_res,
        grid_L=list(grid_L),
        grid_type=grid_type,
        grid_pts=grid_pts,
        grid_abc=grid_abc,
        m0=1 / np.sqrt(3),
        A0=1.3e-11,
        Ms=8e5,
        K0=0,
        gamma=0,
        alpha=4.42e3,
        exchPBC=np.asarray(pbc, dtype=np.int32),
        usedemag=0,
        cuda=cuda,
        cvode=False,
        solver='dynamic',
        exch_presize=EXCH_PRESIZE,
    )

    result = problem.run_simulation(
        t_end=1e-12,
        nt=2,
        fct_h_ext=fct_h_ext,
        nt_h_ext=2,
    )

    n_exch, rows, cols, vals, nr, nc = result[7:13]
    if n_exch == 0:
        raise RuntimeError("MagTense returned an empty exchange matrix")
    if n_exch > EXCH_PRESIZE * ntot:
        raise RuntimeError(
            f"The exchange matrix has {n_exch} entries but only {EXCH_PRESIZE * ntot} were "
            f"reserved, so MagTense has returned zeros. Increase EXCH_PRESIZE."
        )

    # The matrix is exported from MKL in COO format with Fortran (one based) indexing
    rows, cols, vals = rows[:n_exch], cols[:n_exch], vals[:n_exch]
    base = 1 if (rows.min() >= 1 and cols.min() >= 1) else 0

    return rows - base, cols - base, vals, nr, nc


def dense_matrix(rows, cols, vals, nr, nc):
    """The triplets returned by exchange_triplets as a dense array."""
    A = np.zeros((nr, nc))
    np.add.at(A, (rows, cols), vals)

    return A


def mesh_exchange_triplets(pts, abc, grid_L, pbc, cuda=False):
    """The exchange matrix of an unstructured mesh given by its centers and its side lengths."""
    return exchange_triplets((len(pts), 1, 1), 'unstructuredPrisms', pts, abc, grid_L, pbc,
                             cuda=cuda)


def exchange_matrix(res, a, pbc, grid_type, cuda=False):
    """The exchange matrix of a regular lattice of cubes, as a dense array."""
    grid_L = [res[0] * a, res[1] * a, res[2] * a]

    if grid_type == 'unstructuredPrisms':
        pts, abc = build_prism_mesh(res, a)
        return dense_matrix(*mesh_exchange_triplets(pts, abc, grid_L, pbc, cuda=cuda))

    return dense_matrix(*exchange_triplets(res, grid_type, None, None, grid_L, pbc, cuda=cuda))


def plane_wave_test(A, res, a, direction, mode=1):
    """Check that a plane wave along the given direction is an eigenvector of A.

    Returns the eigenvalue and the relative residual. Translation invariance along the direction
    - which is what the periodic linking has to produce - is equivalent to a vanishing residual.
    """
    nx, ny, nz = res
    n_d = res[direction]
    L = n_d * a

    # Cell center coordinate along the direction of the wave, in the same tile ordering as A
    idx = np.arange(nx * ny * nz)
    ijk = np.stack([idx % nx, (idx // nx) % ny, idx // (nx * ny)], axis=1)
    x = (ijk[:, direction] + 0.5) * a

    k = 2 * np.pi * mode / L
    phi = np.exp(1j * k * x)

    Aphi = A @ phi
    lam = np.vdot(phi, Aphi) / np.vdot(phi, phi)
    residual = np.linalg.norm(Aphi - lam * phi) / max(np.linalg.norm(Aphi), 1e-300)

    return lam, residual, k


def supercell_deviation(pts, abc, grid_L, pbc, cuda=False):
    """Compare the periodic operator of a mesh with the free operator on a supercell of copies.

    The mesh is copied once in each direction along every periodic direction, the exchange
    operator is assembled on the resulting supercell without any periodic boundary conditions,
    and the rows of the central copy are folded back onto the tiles of the original mesh. Every
    tile of the central copy has exactly the neighbourhood that the periodic linking gives the
    corresponding tile of the original mesh - the stencil reaches one layer of cells and every
    copy is at least three cells thick - so the two operators have to agree entry by entry.
    Nothing about the mesh enters, so unlike the plane wave test this works on an irregular mesh
    as well.

    Returns the largest deviation relative to the largest coupling and the number of copies used.
    """
    ntot = len(pts)
    grid_L = np.asarray(grid_L, dtype=float)
    pbc = np.asarray(pbc)

    A_pbc = dense_matrix(*mesh_exchange_triplets(pts, abc, grid_L, pbc, cuda=cuda))

    # All the copies of the mesh that a tile can be coupled to, including the diagonal ones: a
    # tile in a corner of the domain is linked across two or three periodic boundaries at once
    shifts_cd = np.array(list(itertools.product(*[(-1., 0., 1.) if p else (0.,) for p in pbc])))
    shifts_cd = shifts_cd * grid_L
    central = int(np.flatnonzero((shifts_cd == 0).all(axis=1))[0])

    super_pts = np.concatenate([pts + shift_d for shift_d in shifts_cd], axis=0)
    super_abc = np.concatenate([abc] * len(shifts_cd), axis=0)
    super_L = grid_L * np.where(pbc != 0, 3, 1)

    rows, cols, vals, _, _ = mesh_exchange_triplets(super_pts, super_abc, super_L, [0, 0, 0],
                                                    cuda=cuda)

    # Fold the rows of the central copy back onto the original mesh. The copies are stored as
    # consecutive blocks of ntot tiles, so column c of the supercell is tile c % ntot of the mesh
    in_copy = np.logical_and(rows >= central * ntot, rows < (central + 1) * ntot)
    A_folded = np.zeros((ntot, ntot))
    np.add.at(A_folded, (rows[in_copy] - central * ntot, cols[in_copy] % ntot), vals[in_copy])

    return np.abs(A_folded - A_pbc).max() / np.abs(A_pbc).max(), len(shifts_cd)


def supercell_check(pts, abc, grid_L, pbc, label, cuda=False):
    """Run the supercell reference on a mesh, print the outcome and return it as a check."""
    deviation, n_copies = supercell_deviation(pts, abc, grid_L, pbc, cuda=cuda)
    status = "pass" if deviation < SUPERCELL_TOL else "FAIL"
    print(f"  supercell       : {deviation:.3e} away from the operator on {n_copies} copies of "
          f"the mesh [{status}]")

    return deviation, {
        'check': f'{label}: periodic operator matches the supercell reference',
        'value': deviation,
        'limit': SUPERCELL_TOL,
        'passed': deviation < SUPERCELL_TOL,
    }


def analyse(A, res, a, pbc, label):
    """Run all the checks on a single exchange matrix and print the outcome."""
    scale = np.abs(A).max()
    row_sum = np.abs(A.sum(axis=1)).max() / scale

    print(f"\n--- {label} ---")
    print(f"  matrix          : {A.shape[0]} x {A.shape[1]}, {np.count_nonzero(A)} nonzeros")
    print(f"  couplings/row   : {np.count_nonzero(A, axis=1).min()} (min) "
          f"{np.count_nonzero(A, axis=1).max()} (max)")
    status = "pass" if row_sum < ROW_SUM_TOL else "FAIL"
    print(f"  max |row sum|   : {row_sum:.3e} [{status}]")

    results = {'row_sum': row_sum, 'eigenvalues': {}, 'residuals': {}}
    for direction, name in enumerate('xyz'):
        if res[direction] < 3:
            continue
        lam, residual, k = plane_wave_test(A, res, a, direction)
        results['eigenvalues'][name] = lam
        results['residuals'][name] = residual
        # A plane wave is only an eigenvector if the operator is periodic along the direction, so
        # a small residual is a pass for a periodic direction and a failure for a free one
        is_eigenvector = residual < EIGENVECTOR_TOL
        if pbc[direction]:
            status = "periodic, pass" if is_eigenvector else "periodic, FAIL"
        else:
            status = "free, FAIL" if is_eigenvector else "free, as expected"
        # The relative deviation from the continuum Laplacian is the discretisation error of the
        # stencil at this wavelength. It is of the order (k*a)^2/12 and is reported, not tested
        deviation = abs(lam.real / (-k**2) - 1)
        print(f"  plane wave {name}    : residual {residual:.3e} [{status}], "
              f"lambda {lam.real:+.4e} vs -k^2 {-k**2:+.4e} "
              f"({100*deviation:.1f} % discretisation error)")

    return results


def dispersion(A, res, a, direction):
    """Eigenvalues of the operator for the plane waves that fit in the periodic domain."""
    modes = np.arange(1, res[direction] // 2 + 1)
    k_m = np.zeros(len(modes))
    lam_m = np.zeros(len(modes))
    for i, mode in enumerate(modes):
        lam, _, k = plane_wave_test(A, res, a, direction, mode=mode)
        k_m[i] = k
        lam_m[i] = lam.real

    return k_m, lam_m


def operator_analysis(res, a, pbc, grid_type, cuda=False):
    """Assemble the periodic and the free operator on one grid type and analyse both."""
    label = GRID_LABELS[grid_type]

    A_pbc = exchange_matrix(res, a, pbc, grid_type, cuda=cuda)
    A_free = exchange_matrix(res, a, [0, 0, 0], grid_type, cuda=cuda)

    results = analyse(A_pbc, res, a, pbc, f"{label}, periodic")
    analyse(A_free, res, a, [0, 0, 0], f"{label}, free boundaries")

    # The rows of the tiles that do not touch a periodic boundary must be untouched by the
    # linking. Reported rather than tested, since it is a property of the implementation
    nx, ny, nz = res
    idx = np.arange(nx * ny * nz)
    ijk = np.stack([idx % nx, (idx // nx) % ny, idx // (nx * ny)], axis=1)
    interior = np.ones(len(idx), dtype=bool)
    for direction in range(3):
        if pbc[direction]:
            interior &= (ijk[:, direction] > 0) & (ijk[:, direction] < res[direction] - 1)
    if interior.any():
        difference = np.abs(A_pbc[interior] - A_free[interior]).max() / np.abs(A_pbc).max()
        status = "pass" if difference < 1e-12 else "FAIL"
        print(f"\n  Rows of the {interior.sum()} tiles away from the periodic boundaries are "
              f"unchanged: max difference {difference:.3e} [{status}]")

    return A_pbc, A_free, results


#%% The irregular mesh

def grain_analysis(pbc, cuda=False, cvode=False):
    """Repeat the tests on an irregular mesh of grains, on which nothing is a regular lattice.

    Returns the checks together with the numbers the figure is drawn from. The sectioned
    simulation and the exchange operator use two different meshes: the simulation needs a few
    cells on either side of the gap, while the supercell reference replicates the mesh 27 times
    and has to stay small.
    """
    print(f"\n{'=' * 78}\n{GRAIN_LABEL}\n{'=' * 78}")
    checks = []

    # ---- Method 1, with a gap through the mesh instead of the sections ---------------------
    pts_all, abc_all, grid_L, grain_all = build_grain_mesh(GRAIN_RES_SECTION, a_section)
    gap_n = grain_gap_mask(pts_all, a_section)
    n_sizes = len(np.unique(np.round(abc_all[:, 0] / abc_all[:, 0].min(), 9)))
    print(f"Sectioned mesh: {len(pts_all)} cells of {n_sizes} different sizes, {gap_n.sum()} of "
          f"them left out to split the mesh in two")
    # The bounding box, and with it the period, is set by the cells that are kept, and those still
    # reach both ends of the domain along all three directions
    keep_n = np.logical_not(gap_n)
    pts, abc = pts_all[keep_n, :], abc_all[keep_n, :]

    # The gap is normal to x, so only the linking along x can connect the two halves. The other
    # two directions are covered by the supercell reference below, which tests all of them at once
    deviation, spread, M_pbc = split_deviations(pts, abc, grid_L, [1, 0, 0], cuda=cuda,
                                                cvode=cvode)
    deviation_free, _, M_free = split_deviations(pts, abc, grid_L, [0, 0, 0], cuda=cuda,
                                                 cvode=cvode)

    verdict = 'works' if deviation < SPLIT_TOL else 'FAILED'
    print(f"  Exchange coupling through the periodic boundary {verdict}: {deviation:.3e} between "
          f"the halves, spread {spread:.3e}")
    verdict = 'as expected' if deviation_free > SPLIT_CONTROL_MIN else 'FAILED'
    print(f"  The same halves without periodic boundaries are decoupled, {verdict}: "
          f"{deviation_free:.3e}")

    checks.append({
        'check': f'{GRAIN_LABEL}: halves coupled through the periodic boundary',
        'value': deviation,
        'limit': SPLIT_TOL,
        'passed': deviation < SPLIT_TOL,
    })
    checks.append({
        'check': f'{GRAIN_LABEL}: both halves relax to the same direction',
        'value': spread,
        'limit': SPLIT_TOL,
        'passed': spread < SPLIT_TOL,
    })
    # The control has to fail the test above, which is what shows that the gap really does cut the
    # mesh in two. It is written as a ratio so that it fits the value < limit contract
    checks.append({
        'check': f'{GRAIN_LABEL}: halves decoupled without periodic boundaries (control)',
        'value': SPLIT_CONTROL_MIN / max(deviation_free, 1e-300),
        'limit': 1.0,
        'passed': deviation_free > SPLIT_CONTROL_MIN,
    })

    # ---- Methods 2 and 3, on a smaller mesh of the same kind -------------------------------
    pts, abc, grid_L, _ = build_grain_mesh(GRAIN_RES_OPERATOR, a_operator)
    A_pbc = dense_matrix(*mesh_exchange_triplets(pts, abc, grid_L, pbc, cuda=cuda))
    A_free = dense_matrix(*mesh_exchange_triplets(pts, abc, grid_L, [0, 0, 0], cuda=cuda))
    scale = np.abs(A_pbc).max()
    row_sum = np.abs(A_pbc.sum(axis=1)).max() / scale

    print(f"\nExchange operator on {len(pts)} irregular cells, exchPBC = {pbc}")
    print(f"  matrix          : {A_pbc.shape[0]} x {A_pbc.shape[1]}, "
          f"{np.count_nonzero(A_pbc)} nonzeros")
    print(f"  couplings/row   : {np.count_nonzero(A_pbc, axis=1).min()} (min) "
          f"{np.count_nonzero(A_pbc, axis=1).max()} (max)")
    status = "pass" if row_sum < ROW_SUM_TOL else "FAIL"
    print(f"  max |row sum|   : {row_sum:.3e} [{status}]")
    # Reported, not tested: how much of the operator the linking changes at all, which is what
    # the supercell reference below has to reproduce. An implementation that quietly ignored
    # exchPBC would leave this at zero and would then miss the supercell by the same amount
    print(f"  linking changes : {np.abs(A_pbc - A_free).max() / scale:.3e} of the largest "
          f"coupling")

    checks.append({
        'check': f'{GRAIN_LABEL}: exchange operator has vanishing row sums',
        'value': row_sum,
        'limit': ROW_SUM_TOL,
        'passed': row_sum < ROW_SUM_TOL,
    })

    supercell = None
    if any(pbc):
        supercell, check = supercell_check(pts, abc, grid_L, pbc, GRAIN_LABEL, cuda=cuda)
        checks.append(check)

    return {'checks': checks, 'deviation': deviation, 'spread': spread,
            'deviation_free': deviation_free, 'supercell': supercell, 'A_pbc': A_pbc,
            'pts': pts_all, 'abc': abc_all, 'grain': grain_all, 'gap': gap_n,
            'M_pbc': M_pbc, 'M_free': M_free}


#%% Run the test

def run_test(res=(6, 6, 6), a=a_operator, pbc=(1, 1, 1), cuda=False, cvode=False,
             plotting: bool = True, grains: bool = True) -> list[dict]:
    """Run all three methods on both grid types and return the checks they consist of.

    Each check is a dict with the keys 'check', 'value', 'limit' and 'passed', where the test
    passes when value < limit. This is the contract used by testMagTenseFunctions.py.

    The res, a and pbc arguments configure the exchange operator test only. The sectioned
    simulation uses the fixed 6 X 6 X 6 geometry that its spacer pattern is built for. The
    irregular mesh, which is run unless grains is False, brings its own geometry along but does
    use the same pbc for its exchange operator.
    """
    res, pbc = list(res), list(pbc)
    checks = []
    deviations = {}
    operators = {}
    supercells = {}

    for grid_type in GRID_TYPES:
        label = GRID_LABELS[grid_type]
        print(f"\n{'=' * 78}\n{label}\n{'=' * 78}")

        # ---- Method 1 ----------------------------------------------------------------
        deviations[grid_type] = section_deviations(grid_type, cuda=cuda, cvode=cvode)
        checks.extend(
            {
                'check': f'{label}: {name} tiles exchange coupled',
                'value': error,
                'limit': SECTION_TOL[grid_type],
                'passed': error < SECTION_TOL[grid_type],
            }
            for name, error in zip(sectionNames_s, deviations[grid_type], strict=True)
        )

        # ---- Method 2 ----------------------------------------------------------------
        print(f"\nExchange operator on a {res[0]} x {res[1]} x {res[2]} lattice of "
              f"{a:.3e} m cubes, exchPBC = {pbc}")
        A_pbc, A_free, analysis = operator_analysis(res, a, pbc, grid_type, cuda=cuda)
        operators[grid_type] = (A_pbc, A_free)

        checks.append({
            'check': f'{label}: exchange operator has vanishing row sums',
            'value': analysis['row_sum'],
            'limit': ROW_SUM_TOL,
            'passed': analysis['row_sum'] < ROW_SUM_TOL,
        })
        for direction, name in enumerate('xyz'):
            if not pbc[direction] or name not in analysis['residuals']:
                continue
            residual = analysis['residuals'][name]
            checks.append({
                'check': f'{label}: plane wave along {name} is an eigenvector',
                'value': residual,
                'limit': EIGENVECTOR_TOL,
                'passed': residual < EIGENVECTOR_TOL,
            })

        # ---- Method 3 ----------------------------------------------------------------
        # The supercell is built from the tiles of a mesh, so it only applies to the
        # unstructured grid type. The uniform grid has no mesh to replicate
        if grid_type == 'unstructuredPrisms' and any(pbc):
            pts, abc = build_prism_mesh(res, a)
            grid_L = [res[0] * a, res[1] * a, res[2] * a]
            supercells[grid_type], check = supercell_check(pts, abc, grid_L, pbc, label,
                                                           cuda=cuda)
            checks.append(check)

    grain_results = None
    if grains:
        grain_results = grain_analysis(pbc, cuda=cuda, cvode=cvode)
        checks.extend(grain_results['checks'])

    if plotting:
        plot_results(deviations, operators, res, a, pbc, supercells, grain_results)

    return checks


def plot_results(deviations, operators, res, a, pbc, supercells, grains):
    """Six panels: the sectioned test, the dispersion, the irregular mesh, and the matrices."""
    fig, axes = plt.subplots(2, 3, layout='constrained', figsize=(19, 9))

    # ---- Sectioned simulation, both grids side by side ---------------------------------
    ax = axes[0, 0]
    x = np.arange(len(sectionNames_s))
    width = 0.38
    threshold_colours = {'uniform': 'darkgreen', 'unstructuredPrisms': 'darkorange'}
    for i, grid_type in enumerate(GRID_TYPES):
        errors = deviations[grid_type]
        tolerance = SECTION_TOL[grid_type]
        colours = ['forestgreen' if e < tolerance else 'firebrick' for e in errors]
        ax.bar(x + (i - 0.5) * width, errors, width, color=colours,
               edgecolor='black', linewidth=0.6,
               hatch='' if i == 0 else '//', label=GRID_LABELS[grid_type])
        ax.axhline(tolerance, color=threshold_colours[grid_type], linestyle='--',
                   label=f'Threshold, {GRID_LABELS[grid_type]}')
    # A section that is not coupled through the periodic boundary keeps its random initial
    # directions, which is where the test would land if the periodic linking were missing
    ax.axhline(1.0, color='firebrick', linestyle=':', label='Uncoupled sections')
    ax.set_yscale('log')
    ax.set_xticks(x)
    ax.set_xticklabels(sectionNames_s, rotation=25, ha='right')
    ax.set_xlabel('Periodically connected section')
    ax.set_ylabel(r'Maximum $|\mathbf{m}_i - \langle\mathbf{m}\rangle|$')
    ax.set_title('Sectioned simulation')
    ax.legend(fontsize=8)

    # ---- Dispersion of the periodic operator, both grids -------------------------------
    ax = axes[0, 1]
    # The two grid types land on top of one another, so the uniform grid is drawn as large
    # open circles with the unstructured mesh marked inside them
    styles = {'uniform': dict(marker='o', markersize=11, markerfacecolor='none',
                              color='navy', markeredgewidth=1.6),
              'unstructuredPrisms': dict(marker='x', markersize=7, color='crimson',
                                         markeredgewidth=1.6)}
    for grid_type in GRID_TYPES:
        A_pbc = operators[grid_type][0]
        for direction, name in enumerate('xyz'):
            if not pbc[direction] or res[direction] < 3:
                continue
            k_m, lam_m = dispersion(A_pbc, res, a, direction)
            ax.plot(k_m * a, -lam_m * a**2, linestyle='',
                    label=f'{GRID_LABELS[grid_type]}, {name}', **styles[grid_type])
    k_fine = np.linspace(0, np.pi, 200)
    ax.plot(k_fine, k_fine**2, 'k-', label=r'Continuum, $k^2$')
    ax.plot(k_fine, 4 * np.sin(k_fine / 2)**2, 'k--', label=r'7 point stencil')
    ax.set_xlabel(r'$k a$')
    ax.set_ylabel(r'$-\lambda a^2$')
    ax.set_title('Dispersion of the periodic operator')
    ax.legend(fontsize=8)

    # ---- The irregular mesh and the supercell reference ---------------------------------
    # Every entry here is a dimensionless deviation that has to stay below its own threshold,
    # except for the control, which has to stay above it
    ax = axes[0, 2]
    bars = []
    if grains is not None:
        bars.append(('Split mesh, halves', grains['deviation'], SPLIT_TOL, False))
        bars.append(('Split mesh, spread', grains['spread'], SPLIT_TOL, False))
        bars.append(('Split mesh, free (control)', grains['deviation_free'], SPLIT_CONTROL_MIN,
                     True))
    for grid_type in GRID_TYPES:
        if supercells.get(grid_type) is not None:
            bars.append((f'Supercell, {GRID_LABELS[grid_type]}', supercells[grid_type],
                         SUPERCELL_TOL, False))
    if grains is not None and grains['supercell'] is not None:
        bars.append((f'Supercell, {GRAIN_LABEL}', grains['supercell'], SUPERCELL_TOL, False))

    if bars:
        bar_labels, values, limits, above = zip(*bars)
        y = np.arange(len(bars))
        colours = ['forestgreen' if (v > lim if ab else v < lim) else 'firebrick'
                   for v, lim, ab in zip(values, limits, above)]
        ax.barh(y, values, color=colours, edgecolor='black', linewidth=0.6)
        ax.plot(limits, y, linestyle='', marker='|', markersize=22, markeredgewidth=2,
                color='black', label='Threshold')
        ax.set_xscale('log')
        ax.set_xlim(min(min(values), min(limits)) / 10, 10 * max(max(values), max(limits)))
        ax.set_yticks(y)
        ax.set_yticklabels(bar_labels, fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel('Deviation')
        ax.set_title('Irregular mesh and supercell reference')
        ax.legend(fontsize=8)
    else:
        ax.axis('off')

    # ---- Sparsity patterns of the periodic operators -----------------------------------
    matrices = [(GRID_LABELS[grid_type], operators[grid_type][0]) for grid_type in GRID_TYPES]
    if grains is not None:
        matrices.append((GRAIN_LABEL, grains['A_pbc']))
    for i in range(3):
        ax = axes[1, i]
        if i >= len(matrices):
            ax.axis('off')
            continue
        label, A_pbc = matrices[i]
        ax.spy(np.abs(A_pbc) > 1e-14 * np.abs(A_pbc).max(), markersize=1.2, color='navy')
        ax.set_title(f'Exchange matrix, {label}, periodic {pbc}')
        ax.set_xlabel('Tile index')
        ax.set_ylabel('Tile index')

    figure_path = output_dir / 'periodic_exchange_test.png'
    fig.savefig(figure_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved figure to {figure_path}")

    if grains is not None:
        plot_grain_figure(grains)


def draw_cells(ax, pts, abc, colours, edgecolour='black', linewidth=0.4):
    """Draw the cross sections of a set of cells as rectangles in the xy plane."""
    cells = [Rectangle((x - dx / 2, y - dy / 2), dx, dy) for x, y, dx, dy
             in zip(pts[:, 0], pts[:, 1], abc[:, 0], abc[:, 1])]
    ax.add_collection(PatchCollection(cells, facecolors=colours, edgecolors=edgecolour,
                                      linewidths=linewidth))


def draw_magnetisation(ax, pts, abc, M_pv, grid_L, a_base, title):
    """One cross section of the split mesh, with the magnetisation drawn on the cells.

    The colour of a cell is the out of plane component and the arrow is the part of the moment
    that lies in the plane, which is the usual way of looking at a micromagnetic state. The two
    halves of the mesh carry the same direction only if they are exchange coupled through the
    periodic boundary, so the two panels can be told apart at a glance.
    """
    colours = plt.get_cmap('coolwarm')((M_pv[:, 2] + 1) / 2)
    draw_cells(ax, pts, abc, colours, edgecolour='dimgray', linewidth=0.3)
    # The arrows are scaled by the cell they belong to, so that the small cells at the grain
    # boundaries do not end up with arrows reaching across their neighbours
    ax.quiver(pts[:, 0], pts[:, 1], M_pv[:, 0] * abc[:, 0], M_pv[:, 1] * abc[:, 0],
              angles='xy', scale_units='xy', scale=1 / 0.9, pivot='mid', width=0.004,
              color='black')
    set_up_cross_section(ax, grid_L, a_base, title)


def set_up_cross_section(ax, grid_L, a_base, title):
    """The frame shared by the panels: the gap, the periodic boundaries and the axes."""
    # The gap that splits the mesh, and the two boundaries that the halves are linked through
    ax.axvspan(-a_base / 2, a_base / 2, color='lightgray', zorder=0)
    for sign in (-1, 1):
        ax.axvline(sign * grid_L[0] / 2, color='darkgreen', linestyle='--', linewidth=1.5)
    ax.set_xlim(-0.62 * grid_L[0], 0.62 * grid_L[0])
    ax.set_ylim(-0.52 * grid_L[1], 0.52 * grid_L[1])
    ax.set_aspect('equal')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_title(title, fontsize=11)


def plot_grain_figure(grains):
    """The mesh itself: the grains it is built from and the states the two simulations end in.

    A cross section is enough to see everything that matters, and the plane is placed away from
    the faces of the cells - those sit at multiples of half a base cell - so that it cuts every
    cell of the mesh it passes through exactly once.
    """
    pts, abc, grain_n, gap_n = grains['pts'], grains['abc'], grains['grain'], grains['gap']
    grid_L = [GRAIN_RES_SECTION * a_section] * 3
    z0 = a_section / 3
    in_plane_n = np.abs(pts[:, 2] - z0) < abc[:, 2] / 2
    keep_n = np.logical_not(gap_n)

    fig, axes = plt.subplots(2, 2, layout='constrained', figsize=(13, 11))

    # ---- The mesh and the grains it was refined around ---------------------------------
    ax = axes[0, 0]
    shown_n = np.logical_and(in_plane_n, keep_n)
    colours = plt.get_cmap('tab10')(grain_n[shown_n] % 10)
    draw_cells(ax, pts[shown_n], abc[shown_n], colours)
    # The cells that are left out are drawn as empty boxes, since they are what splits the mesh
    dropped_n = np.logical_and(in_plane_n, gap_n)
    draw_cells(ax, pts[dropped_n], abc[dropped_n], 'none', edgecolour='gray', linewidth=0.4)
    set_up_cross_section(ax, grid_L, a_section,
                         f'{GRAIN_COUNT} grains, refined at the grain boundaries')
    ax.text(0.0, 0.47 * grid_L[1], 'gap', ha='center', va='top', fontsize=10)
    for sign in (-1, 1):
        ax.text(sign * 0.53 * grid_L[0], 0.0, 'periodic boundary', color='darkgreen',
                rotation=90, ha='center', va='center', fontsize=9)

    # ---- The state the simulations start from, and the two they end in -----------------
    m0_pv = random_directions(int(keep_n.sum()))
    shown_of_kept_n = in_plane_n[keep_n]
    states = [
        (axes[0, 1], m0_pv, 'Initial state, random directions'),
        (axes[1, 0], grains['M_pbc'],
         'Relaxed with exchPBC = [1, 0, 0]\nDistance between the mean of the two halves: '
         f"{grains['deviation']:.1e}"),
        (axes[1, 1], grains['M_free'],
         'Relaxed with exchPBC = [0, 0, 0], the control\nDistance between the mean of the two '
         f"halves: {grains['deviation_free']:.1e}"),
    ]
    for ax, M_pv, title in states:
        draw_magnetisation(ax, pts[shown_n], abc[shown_n], M_pv[shown_of_kept_n], grid_L,
                           a_section, title)

    colourbar = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-1, 1), cmap='coolwarm'),
                             ax=axes, shrink=0.4)
    colourbar.set_label(r'$m_z$')
    fig.suptitle(f'The irregular mesh of the periodic exchange test, cross section at '
                 f'z = {z0:.2e} m')

    figure_path = output_dir / 'periodic_exchange_test_grains.png'
    fig.savefig(figure_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved figure to {figure_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--res', type=int, nargs=3, default=[6, 6, 6],
                        help="Number of cells along x, y and z in the exchange operator test")
    parser.add_argument('--a', type=float, default=a_operator,
                        help="Side length of a cell in the exchange operator test [m]")
    parser.add_argument('--pbc', type=int, nargs=3, default=[1, 1, 1],
                        help="Periodic exchange boundary conditions along x, y and z")
    parser.add_argument('--cuda', action='store_true', help="Use the CUDA solver")
    parser.add_argument('--no-grains', action='store_true',
                        help="Skip the irregular mesh of grains")
    args = parser.parse_args()

    checks = run_test(res=args.res, a=args.a, pbc=args.pbc, cuda=args.cuda,
                      grains=not args.no_grains)

    print(f"\n{'=' * 78}")
    for check in checks:
        status = 'pass' if check['passed'] else 'FAIL'
        print(f"  [{status}] {check['check']}: {check['value']:.3e} "
              f"(limit {check['limit']:.1e})")

    passed = all(c['passed'] for c in checks)
    print('\nperiodic_exchange_test ' + ('PASSED' if passed else 'FAILED'))

    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
