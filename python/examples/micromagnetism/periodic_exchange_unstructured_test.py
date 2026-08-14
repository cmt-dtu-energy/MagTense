"""Test the periodic exchange coupling on an unstructured mesh of uniform prisms.

The tiles at the two ends of a periodic direction are linked together in the mesh analysis
(UnstructuredMeshAnalysis.f90), so that they share a face and enter each others interpolation
stencils exactly as an ordinary pair of neighbouring tiles. The exchange operator computed from
the mesh (ComputeDifferentialOperatorsFromMesh.f90) then needs no special treatment of the
periodic boundaries, apart from measuring the distance between a linked pair through the boundary.

A regular lattice of identical cubes is used as the unstructured mesh, so that the resulting
exchange matrix can be compared both with an exact result and with the uniform grid implementation,
which has supported periodic boundaries all along.

Three properties are tested, none of which depend on the details of the discretisation:

  1. Zero row sums. A uniform magnetisation has to give a vanishing exchange field.

  2. Plane waves are eigenvectors. On a periodic lattice the exchange operator is invariant under
     a translation by one cell, and a translation invariant operator has the plane waves
     exp(i*k*x), k = 2*pi*m/L, as exact eigenvectors. The residual |A*phi - lambda*phi| is
     therefore at round-off level if, and only if, every tile - including the ones at the
     boundaries - is coupled in the same way. This is the sharp test of the periodic linking.

  3. The eigenvalues follow the Laplacian dispersion relation, lambda(k) -> -k^2 for k -> 0, so
     that the operator is not just periodic but is still the correct second derivative.

The same quantities are computed for the uniform grid implementation for reference, and the
unstructured operator is also computed without periodic boundaries to show what the linking
changes.

NOTE: the Fortran library has to be rebuilt for this test to run.

Usage:
    python periodic_exchange_unstructured_test.py                  # 6 x 6 x 6 cells, PBC along x, y and z
    python periodic_exchange_unstructured_test.py --res 8 4 4 --pbc 1 0 0
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

# Acceptance thresholds
ROW_SUM_TOL = 1e-8      # Relative to the largest coupling in the matrix
EIGENVECTOR_TOL = 1e-8  # Relative residual of the plane wave eigenvectors

# Number of exchange matrix entries per tile that MagTense reserves for the returned matrix
EXCH_PRESIZE = 64


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


def exchange_matrix(res, a, pbc, unstructured=True, cuda=False):
    """Run a minimal simulation and return the assembled exchange matrix as a dense array.

    Only the setup of the problem is of interest here, so the simulation is run for a couple of
    trivial timesteps with the demagnetisation field and the external field switched off.
    """
    nx, ny, nz = res
    ntot = nx * ny * nz
    pts, abc = build_prism_mesh(res, a)

    if unstructured:
        grid_type = "unstructuredPrisms"
        problem_res = (ntot, 1, 1)
        grid_pts, grid_abc = pts, abc
    else:
        grid_type = "uniform"
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
        # Room for the exchange matrix that is handed back. The stencil of the unstructured
        # operator reaches all the tiles sharing a vertex with a face, so it is far wider than the
        # seven point stencil of the uniform grid. MagTense silently returns zeros if this is
        # too small, which is caught below.
        exch_presize=EXCH_PRESIZE,
    )

    result = problem.run_simulation(
        t_end=1e-12,
        nt=2,
        fct_h_ext=lambda t: np.zeros([np.size(t), 3]),
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


def run_test(res=(6, 6, 6), a=1e-9, pbc=(1, 1, 1), cuda=False, skip_uniform=False,
             plotting: bool = True) -> list[dict]:
    """Run the periodic exchange checks on the unstructured mesh and return them.

    Each check is a dict with the keys 'check', 'value', 'limit' and 'passed', where the test
    passes when value < limit. This is the contract used by testMagTenseFunctions.py. Only the
    periodic operator is checked; the free-boundary and uniform-grid comparisons that are also
    printed are diagnostics rather than tests.
    """
    res, pbc = list(res), list(pbc)
    print(f"Mesh: {res[0]} x {res[1]} x {res[2]} cubes of {a:.3e} m, exchPBC = {pbc}")

    # The exchange operator on the unstructured mesh, with and without the periodic linking
    A_pbc = exchange_matrix(res, a, pbc, unstructured=True, cuda=cuda)
    A_free = exchange_matrix(res, a, [0, 0, 0], unstructured=True, cuda=cuda)

    res_pbc = analyse(A_pbc, res, a, pbc, "Unstructured prisms, periodic")
    analyse(A_free, res, a, [0, 0, 0], "Unstructured prisms, free boundaries")

    # The rows of the tiles that do not touch a periodic boundary must be untouched by the linking
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
        print(f"\nRows of the {interior.sum()} tiles away from the periodic boundaries are "
              f"unchanged: max difference {difference:.3e} [{status}]")

    # Reference: the uniform grid implementation, which has always supported periodic boundaries
    A_unif = None
    if not skip_uniform:
        A_unif = exchange_matrix(res, a, pbc, unstructured=False, cuda=cuda)
        analyse(A_unif, res, a, pbc, "Uniform grid, periodic (reference)")
        difference = np.abs(A_pbc - A_unif).max() / np.abs(A_unif).max()
        print(f"\nLargest difference between the unstructured and the uniform operator: "
              f"{difference:.3e}")
        print("  (the two use different stencils, so this is informative rather than a test)")

    if plotting:
        # Plot the sparsity patterns and the dispersion relation
        fig, axes = plt.subplots(1, 3, layout='constrained', figsize=(14, 4.4))

        for ax, matrix, title in zip(axes[:2], [A_free, A_pbc],
                                     ['Free boundaries', f'Periodic {pbc}'], strict=True):
            ax.spy(np.abs(matrix) > 1e-14 * np.abs(matrix).max(), markersize=1.2, color='navy')
            ax.set_title(f'Exchange matrix, {title}')
            ax.set_xlabel('Tile index')
            ax.set_ylabel('Tile index')

        ax = axes[2]
        for direction, name in enumerate('xyz'):
            if not pbc[direction] or res[direction] < 3:
                continue
            k_m, lam_m = dispersion(A_pbc, res, a, direction)
            ax.plot(k_m * a, -lam_m * a**2, 'o', label=f'MagTense, {name}')
            if A_unif is not None:
                k_u, lam_u = dispersion(A_unif, res, a, direction)
                ax.plot(k_u * a, -lam_u * a**2, 'x', label=f'Uniform grid, {name}')
        k_fine = np.linspace(0, np.pi, 200)
        ax.plot(k_fine, k_fine**2, 'k-', label=r'Continuum, $k^2$')
        ax.plot(k_fine, 4 * np.sin(k_fine / 2)**2, 'k--', label=r'7 point stencil')
        ax.set_xlabel(r'$k a$')
        ax.set_ylabel(r'$-\lambda a^2$')
        ax.set_title('Dispersion of the periodic operator')
        ax.legend(fontsize=9)

        results_dir = Path(__file__).resolve().parent / "results"
        results_dir.mkdir(parents=True, exist_ok=True)
        figure_path = results_dir / 'periodic_exchange_unstructured_test.png'
        fig.savefig(figure_path, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"\nSaved figure to {figure_path}")

    checks = [{
        'check': 'exchange operator has vanishing row sums',
        'value': res_pbc['row_sum'],
        'limit': ROW_SUM_TOL,
        'passed': res_pbc['row_sum'] < ROW_SUM_TOL,
    }]
    for direction, name in enumerate('xyz'):
        if not pbc[direction] or name not in res_pbc['residuals']:
            continue
        residual = res_pbc['residuals'][name]
        checks.append({
            'check': f'plane wave along {name} is an eigenvector',
            'value': residual,
            'limit': EIGENVECTOR_TOL,
            'passed': residual < EIGENVECTOR_TOL,
        })

    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--res', type=int, nargs=3, default=[6, 6, 6],
                        help="Number of cells along x, y and z")
    parser.add_argument('--a', type=float, default=1e-9, help="Side length of a cell [m]")
    parser.add_argument('--pbc', type=int, nargs=3, default=[1, 1, 1],
                        help="Periodic exchange boundary conditions along x, y and z")
    parser.add_argument('--cuda', action='store_true', help="Use the CUDA solver")
    parser.add_argument('--skip-uniform', action='store_true',
                        help="Skip the comparison with the uniform grid implementation")
    args = parser.parse_args()

    checks = run_test(res=args.res, a=args.a, pbc=args.pbc, cuda=args.cuda,
                      skip_uniform=args.skip_uniform)

    passed = all(c['passed'] for c in checks)
    print("\nPeriodic exchange on the unstructured mesh: " + ("PASSED" if passed else "FAILED"))

    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
