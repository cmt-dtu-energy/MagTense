"""
Test that a micromagnetic problem without exchange, A0 = 0 in every tile, runs on every grid type.

With A0 = 0 everywhere the exchange operator is identically zero. MagTense used to normalise A0 by
its largest value before building the operator, which is 0/0 in that case. On the uniform grid the
result was never read, but on the unstructured meshes the NaN zeroed every interpolation weight,
the sparse matrices came out empty and the solver segfaulted. MagTense now detects the case, uses
an explicit zero operator and skips the exchange field altogether.

The same small cube is run on the uniform grid, on an unstructured mesh of prisms and on a
tetrahedral mesh, with demagnetisation, an applied field and a non-uniform initial state, so that
the moments do move. For each grid three things are checked:

1. The exchange field is exactly zero throughout the run.
2. The exchange operator handed back by the solver is exactly zero.
3. The magnetisation agrees with a run at A0 = 1e-30 J/m in every tile. That run goes through the
   full construction of the exchange operator, so this checks that the shortcut is the A0 -> 0
   limit of the ordinary calculation and not merely something that does not crash.

A fourth check covers the single cell on a uniform grid, the macrospin, whose zero operator is
built by the same routine since this change: A0 > 0 and A0 = 0 have to give identical results
there, as a single cell has no neighbour to be exchange coupled to.

A run at a realistic A0 = 1.3e-11 J/m is included in the figure for comparison; it is not tested,
but it shows that exchange makes a visible difference for this problem, so the agreement in check
3 is not a trivial one.

Running the file executes the test and saves a figure. ``run_test()`` returns the same result as a
list of checks, which is the contract the combined suite in testMagTenseFunctions.py expects. The
MATLAB counterpart is matlab/examples/Micromagnetism/MagTense_tests/no_exchange_test.m, with the
same problem and the same limits.
"""

# General modules
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Magtense stuff
from magtense.micromag import MicromagProblem
from magtense.utils import create_tetra_mesh

# Style settings for plots
plt.rcParams['font.size'] = 15
plt.rcParams['text.usetex'] = False

#%% Settings

mu0 = 4 * np.pi * 1e-7

L = 20e-9                   # Side length of the cube [m]
n_side = 4                  # Cells along each side on the uniform grid and the prism mesh
Ms = 8e5                    # Saturation magnetisation [A/m]
A_real = 1.3e-11            # A realistic exchange constant, for the figure only [J/m]
A_tiny = 1e-30              # Small enough to make no difference, large enough to build the operator [J/m]
alpha = 4.42e3              # Damping [m/(A s)]
t_end = 2e-9                # [s]
nt = 21                     # Number of time steps returned

# Applied field, 0.1 T at 30 degrees from z in the xz plane
H_applied = 0.1 / mu0 * np.array([0.5, 0.0, np.sqrt(3) / 2])

# The exchange field and operator have to vanish exactly, so these are only guards against round-off
field_tol = 1e-12           # Largest |H_exc| / Ms
operator_tol = 1e-12        # Largest |entry| of the exchange operator, which is dimensionless
# The A0 = 0 and A0 = A_tiny runs differ by an exchange field of order A_tiny/(mu0 Ms a^2) ~ 1e-14 A/m
limit_tol = 1e-9            # Largest difference of a moment, in units where |m| = 1
macrospin_tol = 1e-12       # Same, for the single cell, where the two runs are identical

output_dir = Path(__file__).resolve().parent


#%% The problem


def initial_state(ntot: int) -> np.ndarray:
    """A deterministic, non-uniform initial state, the same in every language.

    Cell i (0-based) points along the polar angle 0.2 + 1.2*frac(0.618034*i) and the azimuth
    2.4*i, so neighbouring cells differ and exchange would have something to act on.
    """
    i = np.arange(ntot)
    theta = 0.2 + 1.2 * np.mod(0.618034 * i, 1.0)
    phi = 2.4 * i
    return np.column_stack([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])


def grid_settings(grid: str) -> dict:
    """The MicromagProblem arguments that describe each grid."""
    if grid == 'uniform':
        return dict(res=[n_side] * 3, grid_type='uniform', grid_L=[L] * 3)
    if grid == 'macrospin':
        return dict(res=[1, 1, 1], grid_type='uniform', grid_L=[L] * 3)
    if grid == 'unstructuredPrisms':
        a = L / n_side
        centres = (np.arange(n_side) + 0.5) * a - L / 2
        X, Y, Z = np.meshgrid(centres, centres, centres, indexing='ij')
        pts = np.column_stack([X.ravel(order='F'), Y.ravel(order='F'), Z.ravel(order='F')])
        return dict(res=(len(pts), 1, 1), grid_type='unstructuredPrisms', grid_pts=pts,
                    grid_abc=np.full((len(pts), 3), a), grid_L=[L] * 3, exch_presize=64)
    if grid == 'tetrahedron':
        # Three cubes along each side, each split into six tetrahedra: 162 tetrahedra
        nodes, elements, pts = create_tetra_mesh([L] * 3, L / 3)
        return dict(res=(len(pts), 1, 1), grid_type='tetrahedron', grid_nod=nodes,
                    grid_ele=elements.T, grid_pts=pts, grid_L=[L] * 3, exch_presize=64)
    raise ValueError(grid)


def run(grid: str, A0: float, cuda: bool) -> dict:
    """Integrate the LL equation for one grid and one value of A0 in every tile."""
    settings = grid_settings(grid)
    ntot = int(np.prod(settings['res']))
    problem = MicromagProblem(
        A0=A0 * np.ones((ntot, 1)),
        Ms=Ms * np.ones((ntot, 1)),
        K0=0,
        alpha=alpha,
        gamma=0,
        m0=initial_state(ntot),
        cuda=cuda,
        cvode=False,
        usereturnhall=True,
        solver='dynamic',
        **settings,
    )
    problem.use_fmm = 0
    result = problem.run_simulation(
        t_end=t_end, nt=nt,
        fct_h_ext=lambda t: np.tile(H_applied, (len(t), 1)),
        nt_h_ext=2,
    )
    M = np.asarray(result[1])[:, :, 0, :]          # (time, cell, component)
    return {
        't': np.asarray(result[0]),
        'M': M,
        'H_exc': np.asarray(result[3]),
        'exch_v': np.asarray(result[10]),
        'ntot': ntot,
    }


#%% Test

GRIDS = (('uniform', 'uniform grid'),
         ('unstructuredPrisms', 'unstructured mesh'),
         ('tetrahedron', 'tetrahedral mesh'))


def _max_abs(x: np.ndarray) -> float:
    """Largest absolute value, NaN if any entry is not finite, so that a NaN fails the check."""
    x = np.asarray(x, dtype=np.float64)
    if x.size == 0:
        return 0.0
    if not np.all(np.isfinite(x)):
        return float('nan')
    return float(np.max(np.abs(x)))


def run_test(cuda: bool = False, plotting: bool = True) -> list[dict]:
    checks = []
    curves = {}
    for grid, label in GRIDS:
        print(f'\n{label}')
        zero = run(grid, 0.0, cuda)
        tiny = run(grid, A_tiny, cuda)
        real = run(grid, A_real, cuda)

        field = _max_abs(zero['H_exc']) / Ms
        operator = _max_abs(zero['exch_v'])
        difference = _max_abs(zero['M'] - tiny['M'])
        moved = _max_abs(zero['M'][-1] - zero['M'][0])
        print(f'  {zero["ntot"]} cells, largest change of a moment over the run: {moved:.3f}')
        print(f'  max |H_exc| / Ms with A0 = 0:          {field:.3e} (limit {field_tol:.0e})')
        print(f'  max |exchange operator| with A0 = 0:   {operator:.3e} (limit {operator_tol:.0e})')
        print(f'  max |m(A0 = 0) - m(A0 = {A_tiny:.0e})|:  {difference:.3e} (limit {limit_tol:.0e})')
        print(f'  max |m(A0 = 0) - m(A0 = {A_real:.1e})|: {_max_abs(zero["M"] - real["M"]):.3e} (for comparison)')

        checks += [
            {'check': f'{label}: exchange field is zero with A0 = 0',
             'value': field, 'limit': field_tol, 'passed': field < field_tol},
            {'check': f'{label}: exchange operator is zero with A0 = 0',
             'value': operator, 'limit': operator_tol, 'passed': operator < operator_tol},
            {'check': f'{label}: A0 = 0 matches the limit A0 -> 0',
             'value': difference, 'limit': limit_tol, 'passed': difference < limit_tol},
        ]
        curves[label] = (zero, tiny, real)

    # The macrospin. A single cell has no neighbour, so A0 must not matter at all
    print('\nsingle cell')
    zero = run('macrospin', 0.0, cuda)
    real = run('macrospin', A_real, cuda)
    difference = _max_abs(zero['M'] - real['M'])
    print(f'  max |m(A0 = 0) - m(A0 = {A_real:.1e})|: {difference:.3e} (limit {macrospin_tol:.0e})')
    checks.append(
        {'check': 'single cell: A0 > 0 and A0 = 0 give the same result',
         'value': difference, 'limit': macrospin_tol, 'passed': difference < macrospin_tol})

    if plotting:
        colours = ('crimson', 'forestgreen', 'steelblue')
        fig, axes = plt.subplots(1, len(curves), layout='constrained', figsize=(15, 5), sharey=True)
        for ax, (label, (zero, tiny, real)) in zip(axes, curves.items()):
            t_ns = zero['t'] * 1e9
            for c, colour in enumerate(colours):
                ax.plot(t_ns, zero['M'][:, :, c].mean(axis=1), '-', color=colour)
                ax.plot(t_ns[::2], tiny['M'][::2, :, c].mean(axis=1), 'o', color=colour, markerfacecolor='none')
                ax.plot(t_ns, real['M'][:, :, c].mean(axis=1), '--', color=colour)
            ax.set_title(label)
            ax.set_xlabel('t [ns]')
            ax.grid(True, alpha=0.3)
        axes[0].set_ylabel(r'$\langle m_i \rangle$')

        # Colour is the component and the line style the exchange constant, so the legend is
        # split the same way and placed below the panels, where it covers no curve
        handles = [plt.Line2D([], [], color=colour, label=rf'$\langle m_{name} \rangle$')
                   for name, colour in zip('xyz', colours)]
        handles += [plt.Line2D([], [], color='gray', linestyle='-', label=r'$A_0 = 0$'),
                    plt.Line2D([], [], color='gray', linestyle='none', marker='o', markerfacecolor='none',
                               label=r'$A_0 = 10^{-30}$ J/m'),
                    plt.Line2D([], [], color='gray', linestyle='--', label=r'$A_0 = 1.3 \cdot 10^{-11}$ J/m')]
        fig.legend(handles=handles, loc='outside lower center', ncols=len(handles), fontsize='small')

        # Save the figure beside this script and close it so the script does not
        # open an interactive plotting window.
        figure_path = output_dir / 'no_exchange_test.png'
        fig.savefig(figure_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f'Saved figure to {figure_path}')

    return checks


if __name__ == '__main__':
    checks = run_test()
    print('no_exchange_test '
          + ('PASSED' if all(c['passed'] for c in checks) else 'FAILED'))
