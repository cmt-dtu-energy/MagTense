"""
Test that the macrogeometry method reproduces an explicitly replicated domain, cell by cell.

SHORT EXPLANATION : A periodic domain must give exactly the demagnetisation field of a large
domain built from the same number of copies placed side by side.

LONG EXPLANATION :
The macrogeometry method adds the field of 2*n_macro copies of the simulated domain, shifted by
multiples of shiftVec, to the field of the domain itself. When shiftVec equals the size of the
domain, grid_L, the copies tile space without gaps or overlaps, so the result has to be identical
to that of a single domain that is (2*n_macro + 1) times larger along every periodic direction and
carries the same magnetisation pattern in every block. The demagnetisation field in the central
block of that supercell is computed without any periodic boundary conditions, so the comparison
tests the placement of the copies directly, without any analytical model in between.

This is the gapless case that macrogeometry_PBC_test does not cover: that test models a sparse
chain of separated particles, where shiftVec is much larger than the domain. A shift that is off
by a cell - for instance the distance between the centres of the two end cells, (n - 1)*dx, rather
than the full length n*dx - makes the copies overlap by one cell layer. MagTense rejects such a
shiftVec outright, so the last check runs that case in a subprocess and requires it to stop with a
message rather than to return a field.

The domain has a different cell size along x, y and z and a random magnetisation, so a mix-up
between directions or a misplaced copy anywhere shows up. Both the point demagnetisation tensor
and the tensor averaged over the observation cell are tested. As a control the same comparison is
made without the macrogeometry: the periodic copies have to change the field by much more than
the tolerance, otherwise the test could not tell a working implementation from a missing one.

Running the file executes the test and saves a figure. ``run_test()`` returns the same result as a
list of checks, which is the contract the combined suite in testMagTenseFunctions.py expects.
"""

# General modules
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Magtense stuff
from magtense.micromag import MicromagProblem

plt.rcParams['font.size'] = 13
plt.rcParams['text.usetex'] = False

#%% Settings

res = (3, 2, 2)                         # Cells along x, y and z in the periodic domain
cell = np.array([2e-9, 3e-9, 5e-9])     # Cell size along x, y and z [m], deliberately unequal
grid_L = np.array(res) * cell           # Size of the periodic domain [m]

Ms = 8e5                                # Saturation magnetisation [A/m]
seed = 3                                # Seed of the random magnetisation

# (label, periodic directions, copies on each side, averaged tensor). Copies on each side of the
# domain along every periodic direction; the supercell is (2*n + 1) domains long along those.
CASES = (
    ('x', (1, 0, 0), 4, False),
    ('y', (0, 1, 0), 4, False),
    ('z', (0, 0, 1), 4, False),
    ('x, y and z', (1, 1, 1), 2, False),
    ('x, y and z, averaged tensor', (1, 1, 1), 2, True),
)

# The demagnetisation tensor is stored in single precision, so the two fields agree to about
# 1e-7 of the largest field component. The limit leaves room for the different summation order.
match_tol = 1e-5
# The periodic copies have to change the field by at least this much, relative to the largest
# field component, for the comparison above to mean anything
control_min = 1e-2

# Micromagnetic solver settings
cuda = False

output_dir = Path(__file__).resolve().parent

#%% Helpers


def random_magnetisation(n: int) -> np.ndarray:
    """Reproducible random unit vectors, one per cell."""
    rng = np.random.default_rng(seed)
    m = rng.uniform(-1, 1, (n, 3))
    return m / np.linalg.norm(m, axis=1, keepdims=True)


def demag_field(res_: tuple, grid_L_: np.ndarray, m0: np.ndarray, use_avg: bool,
                n_macro: np.ndarray | None = None, shift_vec: np.ndarray | None = None) -> np.ndarray:
    """The demagnetisation field of a uniform grid with magnetisation m0, as an (ntot, 3) array.

    Only the field of the initial state is needed, so the simulation is run for a vanishing time
    with the precession switched off and the field is read at the first output time.
    """
    ntot = int(np.prod(res_))
    problem = MicromagProblem(
        res=list(res_),
        grid_L=list(grid_L_),
        grid_type='uniform',
        A0=1e-20,
        Ms=Ms,
        K0=0,
        alpha=1e3,
        gamma=0,
        m0=m0,
        cuda=cuda,
        cvode=False,
        useavgn=use_avg,
        usereturnhall=True,
        solver='dynamic',
        n_macro=np.zeros(3) if n_macro is None else n_macro,
        shiftVec=np.zeros(3) if shift_vec is None else shift_vec,
    )
    # An octree is pointless for grids this small, and the FMM path ignores the macrogeometry
    problem.use_fmm = 0
    result = problem.run_simulation(
        t_end=1e-15,
        nt=2,
        fct_h_ext=lambda t: np.zeros((len(t), 3)),
        nt_h_ext=2,
    )
    M_out, H_dem = result[1], result[5]
    # The returned magnetisation at the first output time, normalised to |m| = 1, has to be the
    # initial state, otherwise the fields compared below would belong to different magnetisations
    assert np.allclose(M_out[0, :, 0, :], m0, rtol=0, atol=1e-6)
    return np.asarray(H_dem[0, :, 0, :]).reshape(ntot, 3)


def run_case(pbc: tuple, n_copies: int, use_avg: bool) -> dict:
    """Compare the periodic domain with the central block of its supercell."""
    pbc = np.array(pbc, dtype=bool)
    copies = np.where(pbc, 2 * n_copies + 1, 1)
    nx, ny, nz = res
    ntot = nx * ny * nz

    m0 = random_magnetisation(ntot)

    # The periodic domain
    n_macro = np.where(pbc, n_copies, 0)
    shift_vec = np.where(pbc, grid_L, 0.0)
    H_pbc = demag_field(res, grid_L, m0, use_avg, n_macro, shift_vec)

    # The same domain without any copies, the control
    H_free = demag_field(res, grid_L, m0, use_avg)

    # The supercell. Cells are numbered with x running fastest, as setupGrid in
    # LandauLifshitzEquationSolver.f90 numbers them, so an (nz, ny, nx, 3) array in C order is
    # the grid itself and np.tile lays the copies out side by side
    m0_grid = m0.reshape(nz, ny, nx, 3)
    m0_super = np.tile(m0_grid, (copies[2], copies[1], copies[0], 1))
    res_super = tuple(int(r) for r in np.array(res) * copies)
    H_super = demag_field(res_super, grid_L * copies, m0_super.reshape(-1, 3), use_avg)

    # The central block of the supercell coincides with the periodic domain: the supercell is
    # centred on the origin just like the domain, and it is an odd number of domains long
    offset = np.where(pbc, n_copies, 0) * np.array(res)
    H_super_grid = H_super.reshape(res_super[2], res_super[1], res_super[0], 3)
    H_central = H_super_grid[offset[2]:offset[2] + nz,
                             offset[1]:offset[1] + ny,
                             offset[0]:offset[0] + nx].reshape(ntot, 3)

    scale = np.abs(H_central).max()
    return {
        'mismatch': float(np.abs(H_pbc - H_central).max() / scale),
        'control': float(np.abs(H_free - H_central).max() / scale),
        'H_pbc': H_pbc,
        'H_central': H_central,
        'H_free': H_free,
        'n_super': int(np.prod(res_super)),
    }


# The script run in a subprocess for the overlap check. It sets up the domain with the shift that
# results from measuring between the centres of the two end cells instead of between the outer
# faces, so the copies overlap by one cell layer. MagTense has to stop on that rather than return
# a field, and stopping takes the whole interpreter with it, hence the subprocess.
OVERLAP_SCRIPT = textwrap.dedent('''
    import numpy as np
    from magtense.micromag import MicromagProblem
    res = {res}
    cell = np.array({cell})
    grid_L = np.array(res) * cell
    shift = np.array([(res[0] - 1) * cell[0], 0.0, 0.0])
    problem = MicromagProblem(res=list(res), grid_L=list(grid_L), grid_type='uniform', A0=1e-20,
                              Ms=8e5, K0=0, alpha=1e3, gamma=0, cvode=False, cuda=False,
                              useavgn=False, solver='dynamic',
                              n_macro=np.array([2, 0, 0]), shiftVec=shift)
    problem.use_fmm = 0
    problem.run_simulation(t_end=1e-15, nt=2, fct_h_ext=lambda t: np.zeros((len(t), 3)),
                           nt_h_ext=2)
    print('OVERLAP NOT DETECTED')
''').format(res=res, cell=list(cell))


def overlapping_copies_rejected() -> tuple[bool, str]:
    """Run OVERLAP_SCRIPT in a subprocess and report whether MagTense refused it."""
    env = dict(os.environ)
    env['PYTHONPATH'] = os.pathsep.join(p for p in sys.path if p)
    proc = subprocess.run([sys.executable, '-c', OVERLAP_SCRIPT], capture_output=True,
                          text=True, env=env, timeout=600)
    output = proc.stdout + proc.stderr
    # The message of the error stop in checkMacrogeometrySpacing, so that the subprocess failing
    # for any other reason is not taken for a rejection
    rejected = proc.returncode != 0 and 'the macrogeometry copies overlap' in output \
        and 'OVERLAP NOT DETECTED' not in output
    return rejected, output


#%% Run the test


def run_test(plotting: bool = True) -> list[dict]:
    """Run the supercell test and return the checks it consists of.

    Each check is a dict with the keys 'check', 'value', 'limit' and 'passed', where the test
    passes when value < limit. This is the contract used by testMagTenseFunctions.py.
    """
    print(f'Periodic domain: {res[0]} x {res[1]} x {res[2]} cells of '
          f'{cell[0]*1e9:g} x {cell[1]*1e9:g} x {cell[2]*1e9:g} nm, shiftVec = grid_L')

    checks = []
    results = {}
    for label, pbc, n_copies, use_avg in CASES:
        r = run_case(pbc, n_copies, use_avg)
        results[label] = r
        ok = r['mismatch'] < match_tol
        print(f'  periodic along {label}, {n_copies} copies on each side '
              f'({r["n_super"]} cells in the supercell): '
              f'mismatch {r["mismatch"]:.2e} [{"pass" if ok else "FAIL"}], '
              f'the copies change the field by {r["control"]:.2e}')
        checks.append({
            'check': f'periodic along {label}: field matches the supercell',
            'value': r['mismatch'],
            'limit': match_tol,
            'passed': ok,
        })
        # Written as a ratio so that it fits the value < limit contract
        checks.append({
            'check': f'periodic along {label}: copies change the field (control)',
            'value': control_min / max(r['control'], 1e-300),
            'limit': 1.0,
            'passed': r['control'] > control_min,
        })

    rejected, output = overlapping_copies_rejected()
    print(f'  shiftVec = (n - 1)*dx, copies overlapping by one cell: '
          f'{"rejected" if rejected else "NOT REJECTED"}')
    if not rejected:
        print(textwrap.indent(output[-2000:], '    | '))
    checks.append({
        'check': 'overlapping copies are rejected',
        'value': 0.0 if rejected else 1.0,
        'limit': 0.5,
        'passed': rejected,
    })

    if plotting:
        fig, ax = plt.subplots(layout='constrained', figsize=(7, 4.5))
        labels = [c[0] for c in CASES]
        x = np.arange(len(labels))
        mismatch = [results[lbl]['mismatch'] for lbl in labels]
        control = [results[lbl]['control'] for lbl in labels]
        ax.bar(x - 0.2, mismatch, 0.4, color='forestgreen', label='Periodic vs supercell')
        ax.bar(x + 0.2, control, 0.4, color='grey', label='No copies vs supercell (control)')
        ax.axhline(match_tol, color='forestgreen', ls='--', label='Tolerance')
        ax.axhline(control_min, color='grey', ls=':', label='Control minimum')
        ax.set_yscale('log')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=15)
        ax.set_xlabel('Periodic directions')
        ax.set_ylabel(r'max $|\Delta H|$ / max $|H|$')
        ax.legend(loc='best', fontsize=9)
        figure_path = output_dir / 'macrogeometry_supercell_test.png'
        fig.savefig(figure_path, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f'Saved figure to {figure_path}')

    return checks


if __name__ == '__main__':
    checks = run_test()
    print('macrogeometry_supercell_test '
          + ('PASSED' if all(c['passed'] for c in checks) else 'FAILED'))
