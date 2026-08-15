"""Test the implementation of periodic exchange coupling on both grid types of MagTense.

MagTense builds the exchange operator in two entirely different ways: from a stencil on a
uniform grid, and from a mesh analysis on an unstructured mesh of prisms
(UnstructuredMeshAnalysis.f90 and ComputeDifferentialOperatorsFromMesh.f90). Both have to
honour the periodic boundary conditions requested through exchPBC, so both are tested here,
with two independent methods each.

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

NOTE: the Fortran library has to be rebuilt for the unstructured part of this test to run.

Running the file executes the test and saves a figure. ``run_test()`` returns the same result as
a list of checks, which is the contract the combined suite in testMagTenseFunctions.py expects.

Usage:
    python periodic_exchange_test.py                  # 6 x 6 x 6 cells, PBC along x, y and z
    python periodic_exchange_test.py --res 8 4 4 --pbc 1 0 0
"""

# General modules
import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Magtense stuff
from magtense.micromag import MicromagProblem

# Style settings for plots
plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = False

#%% Settings

results_dir = Path(__file__).resolve().parent / "results"

# The two grid types that MagTense assembles the exchange operator for
GRID_TYPES = ('uniform', 'unstructuredPrisms')
GRID_LABELS = {'uniform': 'uniform grid', 'unstructuredPrisms': 'unstructured mesh'}

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

def build_sections():
    """Return the boolean masks of the 8 periodically connected sections and the spacers.

    Tiles are ordered with z running fastest, i.e. n = k + 6*j + 36*i.
    """
    side = SECTION_SIDE
    ntot = side**3

    corner_n, center_n = np.zeros(ntot, dtype=bool), np.zeros(ntot, dtype=bool)
    xFace_n, yFace_n, zFace_n = (np.zeros(ntot, dtype=bool), np.zeros(ntot, dtype=bool),
                                 np.zeros(ntot, dtype=bool))
    xyEdge_n, xzEdge_n, yzEdge_n = (np.zeros(ntot, dtype=bool), np.zeros(ntot, dtype=bool),
                                    np.zeros(ntot, dtype=bool))
    n = 0
    for i in range(side):
        for j in range(side):
            for k in range(side):
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
    for i in range(side):
        for j in range(side):
            for k in range(side):
                pts[n, :] = (np.array([i, j, k]) + 0.5) * a_section - side * a_section / 2
                n += 1
    abc = a_section * np.ones((ntot, 3))
    return pts, abc


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

    # Initial magnetisation directions
    rng = np.random.default_rng(7)
    m0_pv = rng.uniform(-1, 1, [ntot, 3])    # Reproducible random directions
    m0_pv = m0_pv / np.sqrt(np.sum(m0_pv**2, axis=1, keepdims=True))    # Normalise

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
        # The uniform grid has always applied periodic exchange; the unstructured mesh needs
        # to be asked for it explicitly
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


#%% Method 2: the exchange operator

def exchange_matrix(res, a, pbc, grid_type, cuda=False):
    """Run a minimal simulation and return the assembled exchange matrix as a dense array.

    Only the setup of the problem is of interest here, so the simulation is run for a couple of
    trivial timesteps with the demagnetisation field and the external field switched off.
    """
    nx, ny, nz = res
    ntot = nx * ny * nz
    pts, abc = build_prism_mesh(res, a)

    if grid_type == 'unstructuredPrisms':
        problem_res = (ntot, 1, 1)
        grid_pts, grid_abc = pts, abc
    else:
        problem_res = res
        grid_pts, grid_abc = None, None

    problem = MicromagProblem(
        res=problem_res,
        grid_L=[nx * a, ny * a, nz * a],
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
    base = 1 if (rows.min() >= 1 and cols.min() >= 1) else 0
    A = np.zeros((nr, nc))
    np.add.at(A, (rows[:n_exch] - base, cols[:n_exch] - base), vals[:n_exch])

    return A


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


#%% Run the test

def run_test(res=(6, 6, 6), a=a_operator, pbc=(1, 1, 1), cuda=False, cvode=False,
             plotting: bool = True) -> list[dict]:
    """Run both methods on both grid types and return the checks they consist of.

    Each check is a dict with the keys 'check', 'value', 'limit' and 'passed', where the test
    passes when value < limit. This is the contract used by testMagTenseFunctions.py.

    The res, a and pbc arguments configure the exchange operator test only. The sectioned
    simulation uses the fixed 6 X 6 X 6 geometry that its spacer pattern is built for.
    """
    res, pbc = list(res), list(pbc)
    checks = []
    deviations = {}
    operators = {}

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

    if plotting:
        plot_results(deviations, operators, res, a, pbc)

    return checks


def plot_results(deviations, operators, res, a, pbc):
    """Four panels: the sectioned test, the dispersion, and the two exchange matrices."""
    fig, axes = plt.subplots(2, 2, layout='constrained', figsize=(13, 9))

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

    # ---- Sparsity patterns of the two periodic operators -------------------------------
    for i, grid_type in enumerate(GRID_TYPES):
        ax = axes[1, i]
        A_pbc = operators[grid_type][0]
        ax.spy(np.abs(A_pbc) > 1e-14 * np.abs(A_pbc).max(), markersize=1.2, color='navy')
        ax.set_title(f'Exchange matrix, {GRID_LABELS[grid_type]}, periodic {pbc}')
        ax.set_xlabel('Tile index')
        ax.set_ylabel('Tile index')

    results_dir.mkdir(parents=True, exist_ok=True)
    figure_path = results_dir / 'periodic_exchange_test.png'
    fig.savefig(figure_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved figure to {figure_path}")


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
    args = parser.parse_args()

    checks = run_test(res=args.res, a=args.a, pbc=args.pbc, cuda=args.cuda)

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
