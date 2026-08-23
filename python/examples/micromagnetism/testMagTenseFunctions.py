"""Combined micromagnetic test suite for the python version of MagTense.

This is the python counterpart of matlab/util/testMagTenseFunctions.m. It runs the
validation examples in this folder, collects the numerical measure each of them
produces, compares it with the acceptance limit that belongs to it, and prints a table
plus a figure showing which tests passed.

Coverage compared with the MATLAB suite
---------------------------------------
Shared with MATLAB : macrogeometry periodic boundary conditions, periodic exchange on
                     both the uniform grid and the unstructured mesh, shape correction,
                     thermal fluctuations, standard problem 4, standard problem 6. The
                     MATLAB counterparts live in
                     matlab/examples/Micromagnetism/MagTense_tests and use the same
                     geometries and the same acceptance limits.
Only in python     : standard problem 3.
Only in MATLAB     : the six magnetostatic field validations (cylindrical slice, prism,
                     sphere, spheroid, tetrahedron), and the standard problem 6
                     direction and unstructured mesh variants. The python examples for
                     these do not exist or do not expose a comparable error measure, so
                     they are not covered here.

How a test reports its result
-----------------------------
Every example exposes a function that returns a list of checks. A check is a dict

    {'check': str, 'value': float, 'limit': float, 'passed': bool}

and passes when value < limit. Metrics that are naturally two sided, such as the ratio
between a simulated and an analytical diffusion constant, are converted to that form by
the example itself.

Usage
-----
    python testMagTenseFunctions.py                 # everything except the slow tests
    python testMagTenseFunctions.py --include-slow  # add standard problem 3 (very slow)
    python testMagTenseFunctions.py --list          # show the available tests
    python testMagTenseFunctions.py --tests temperature_test,std_problem_4
    python testMagTenseFunctions.py --skip std_problem_6   # all but the slowest one
"""

import argparse
import sys
import time
import traceback
from pathlib import Path

import matplotlib as mpl

# The suite writes figures rather than showing them, so no interactive backend is
# needed. This has to happen before the examples are imported, since they create
# figures at call time.
mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

# The examples live next to this file and are imported as plain modules
EXAMPLE_DIR = Path(__file__).resolve().parent
if str(EXAMPLE_DIR) not in sys.path:
    sys.path.insert(0, str(EXAMPLE_DIR))

RESULTS_DIR = EXAMPLE_DIR / "results"

# Passed on to the examples that accept it. The examples that do not take a cuda
# argument set it themselves; MicromagProblem falls back to the CPU when no GPU is
# present anyway.
USE_CUDA = False

# Single domain limit of the muMag standard problem 3, in units of the exchange length
STD3_SINGLE_DOMAIN_LIMIT = 8.47


#%% Wrappers turning each example into a list of checks


def _macrogeometry_PBC_test() -> list[dict]:
    import macrogeometry_PBC_test as mod
    return mod.run_test()


def _periodic_exchange_test() -> list[dict]:
    import periodic_exchange_test as mod
    return mod.run_test(cuda=USE_CUDA)


def _shape_correction_test() -> list[dict]:
    import shape_correction_test as mod
    return mod.run_test()


def _temperature_test() -> list[dict]:
    import temperature_test as mod
    return mod.run_test()


def _std_problem_4() -> list[dict]:
    """Compare both NIST fields of standard problem 4 with the published mean solutions.

    The limit of 100 % is the one used by the MATLAB suite. It is loose because <My> and
    <Mz> average close to zero, so the relative measure has a small denominator.
    """
    from std_problem_4 import std_prob_4

    checks = []
    for field in (1, 2):
        _, rel_int_error = std_prob_4(
            NIST_field=field,
            cuda=USE_CUDA,
            cvode=False,
            unstructured=False,
            plotting=True,
            figpath=RESULTS_DIR,
        )
        checks.extend(
            {
                'check': f'field {field}: <M{component}> vs mumag',
                'value': error,
                'limit': 100.0,
                'passed': error < 100.0,
            }
            for component, error in zip("xyz", rel_int_error, strict=True)
        )
    return checks


def _std_problem_6() -> list[dict]:
    """Compare the depinning field of standard problem 6 with the analytical values.

    The parameter variations and the 5 % limit are the same as in the MATLAB suite.
    MATLAB additionally runs the field along y and z and on an unstructured mesh; the
    python std_prob_6 does not offer those options, so they are not covered.
    """
    from std_problem_6 import THEORETICAL_PINNING_FIELDS, std_prob_6

    checks = []
    for settings in ("akj", "ak", "aj", "a", "kj", "k"):
        theory = THEORETICAL_PINNING_FIELDS[settings]
        switching_field = std_prob_6(
            settings=settings,
            x_steps=80,
            field_steps=201,
            cuda=USE_CUDA,
            cvode=False,
            plotting=True,
            figpath=RESULTS_DIR,
        )
        if switching_field is None:
            # No switching at all is a failure whatever the limit is
            error = float('inf')
            print(f'settings="{settings}": no switching observed '
                  f'(theory {theory:.3f} T)')
        else:
            error = abs(switching_field - theory) / theory * 100
            print(f'settings="{settings}": switching field = {switching_field:.4f} T '
                  f'(theory {theory:.3f} T, error {error:.1f} %)')
        checks.append({
            'check': f'depinning field, variation "{settings}"',
            'value': error,
            'limit': 5.0,
            'passed': error < 5.0,
        })
    return checks


def _std_problem_3() -> list[dict]:
    """Locate the single domain limit of standard problem 3, compare with the reference.

    The flower and the vortex state swap their role as the ground state at L = 8.47
    exchange lengths. The crossing of the two total energies is found by linear
    interpolation, so the simulated cube sizes have to bracket it.
    """
    from std_problem_3 import std_prob_3

    L_loop = np.linspace(8, 9, 6)
    L_loop, E_arr = std_prob_3(
        L_loop=L_loop,
        cuda=USE_CUDA,
        cvode=False,
        plotting=True,
        figpath=RESULTS_DIR,
    )

    E_flower = np.sum(E_arr[:, :, 0], axis=0)
    E_vortex = np.sum(E_arr[:, :, 1], axis=0)
    # Negative while the flower state is the ground state
    difference = E_flower - E_vortex

    if difference[0] >= 0 or difference[-1] <= 0:
        print("The simulated cube sizes do not bracket the crossing: E_flower - "
              f"E_vortex goes from {difference[0]:.3e} to {difference[-1]:.3e}")
        error = float('inf')
        L_cross = float('nan')
    else:
        L_cross = float(np.interp(0.0, difference, L_loop))
        error = abs(L_cross - STD3_SINGLE_DOMAIN_LIMIT) / STD3_SINGLE_DOMAIN_LIMIT * 100
        print(f"Single domain limit: L = {L_cross:.3f} l_ex "
              f"(accepted value {STD3_SINGLE_DOMAIN_LIMIT}, error {error:.1f} %)")

    return [{
        'check': 'single domain limit L/l_ex',
        'value': error,
        'limit': 5.0,
        'passed': error < 5.0,
    }]


#%% Registry

# name -> (function, slow, one line description)
TESTS = {
    'macrogeometry_PBC_test': (
        _macrogeometry_PBC_test, False,
        'Periodic boundaries by the macrogeometry method, along x, y and z',
    ),
    'periodic_exchange_test': (
        _periodic_exchange_test, False,
        'Periodic exchange coupling, uniform grid, unstructured mesh and grain mesh',
    ),
    'shape_correction_test': (
        _shape_correction_test, False,
        'Shape correction field of the sample geometry',
    ),
    'temperature_test': (
        _temperature_test, False,
        'Thermal fluctuations against the analytical angular diffusion',
    ),
    'std_problem_4': (
        _std_problem_4, False,
        'muMag standard problem 4 against the published mean solutions',
    ),
    'std_problem_6': (
        _std_problem_6, False,
        'muMag standard problem 6, domain wall depinning fields',
    ),
    'std_problem_3': (
        _std_problem_3, True,
        'muMag standard problem 3, single domain limit (slow, tens of minutes)',
    ),
}


#%% Running and reporting


def run_tests(names: list[str]) -> list[dict]:
    """Run the named tests and return one record per test."""
    records = []
    for i, name in enumerate(names, start=1):
        function = TESTS[name][0]
        print(f"\n{'=' * 78}\n[{i}/{len(names)}] {name}\n{'=' * 78}")
        start = time.time()
        try:
            checks = function()
            status = 'PASS' if all(c['passed'] for c in checks) else 'FAIL'
            error = None
        except Exception:
            # A broken example must not stop the rest of the suite
            checks = []
            status = 'ERROR'
            error = traceback.format_exc()
            print(error)
        elapsed = time.time() - start
        print(f"--> {name}: {status} in {elapsed:.1f} s")
        records.append({'name': name, 'checks': checks, 'status': status,
                        'elapsed': elapsed, 'error': error})
    return records


def print_table(records: list[dict], skipped: list[str]) -> None:
    """Print the overview as plain text."""
    rows = _table_rows(records, skipped)
    widths = [max(len(str(row[c])) for row in [_HEADER, *rows])
              for c in range(len(_HEADER))]

    print(f"\n{'=' * (sum(widths) + 3 * len(widths))}")
    print("OVERVIEW")
    print('=' * (sum(widths) + 3 * len(widths)))
    print('   '.join(str(h).ljust(w) for h, w in zip(_HEADER, widths, strict=True)))
    print('-' * (sum(widths) + 3 * len(widths)))
    for row in rows:
        print('   '.join(str(c).ljust(w) for c, w in zip(row, widths, strict=True)))
    print('-' * (sum(widths) + 3 * len(widths)))

    n_pass = sum(r['status'] == 'PASS' for r in records)
    n_fail = sum(r['status'] == 'FAIL' for r in records)
    n_error = sum(r['status'] == 'ERROR' for r in records)
    total_time = sum(r['elapsed'] for r in records)
    print(f"{n_pass} passed, {n_fail} failed, {n_error} errored, "
          f"{len(skipped)} skipped, in {total_time:.0f} s")


_HEADER = ('Test', 'Check', 'Value', 'Limit', 'Status')


def _fmt(value: float) -> str:
    if np.isnan(value):
        return 'nan'
    if np.isinf(value):
        return 'inf'
    if value == 0:
        return '0'
    return f'{value:.3g}' if 1e-3 <= abs(value) < 1e4 else f'{value:.2e}'


def _table_rows(records: list[dict], skipped: list[str]) -> list[tuple]:
    """Flatten the records into table rows, one per check."""
    rows = []
    for record in records:
        label = f"{record['name']} [{record['elapsed']:.0f} s]"
        if record['status'] == 'ERROR':
            first_line = record['error'].strip().splitlines()[-1]
            rows.append((label, f'raised: {first_line[:60]}', '-', '-', 'ERROR'))
            continue
        rows.extend(
            (
                label if j == 0 else '',
                check['check'],
                _fmt(check['value']),
                _fmt(check['limit']),
                'PASS' if check['passed'] else 'FAIL',
            )
            for j, check in enumerate(record['checks'])
        )
        if not record['checks']:
            rows.append((label, 'no checks returned', '-', '-', 'ERROR'))
    rows.extend((name, 'not run', '-', '-', 'SKIPPED') for name in skipped)
    return rows


_STATUS_COLOURS = {
    'PASS': '#c9e7c9',
    'FAIL': '#f2b8b5',
    'ERROR': '#f2b8b5',
    'SKIPPED': '#e0e0e0',
}


def plot_overview(records: list[dict], skipped: list[str], figure_path: Path) -> None:
    """Save a figure with the result table and the margin of every check."""
    rows = _table_rows(records, skipped)
    n_rows = len(rows)

    fig = plt.figure(figsize=(15, 1.6 + 0.30 * (n_rows + 1)))
    grid = fig.add_gridspec(1, 2, width_ratios=[2.6, 1.0], wspace=0.04,
                            left=0.01, right=0.98, top=0.90, bottom=0.06)
    ax_table = fig.add_subplot(grid[0, 0])
    ax_bar = fig.add_subplot(grid[0, 1])

    # ── Table ────────────────────────────────────────────────────────────────
    ax_table.axis('off')
    table = ax_table.table(
        cellText=[[str(c) for c in row] for row in rows],
        colLabels=_HEADER,
        colWidths=[0.30, 0.42, 0.10, 0.10, 0.08],
        cellLoc='left',
        loc='upper left',
        bbox=[0, 0, 1, 1],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    for (row, _col), cell in table.get_celld().items():
        cell.set_linewidth(0.4)
        if row == 0:
            cell.set_facecolor('#404040')
            cell.set_text_props(color='white', fontweight='bold')
        else:
            cell.set_facecolor(_STATUS_COLOURS.get(rows[row - 1][4], 'white'))

    # ── Margin bar chart, aligned row by row with the table ──────────────────
    # A bar is the measured value divided by its limit, so anything left of 1 passed.
    # The rows are laid out so bar i sits at the vertical centre of table row i.
    bar_min, bar_max = 1e-10, 1e4
    for i, row in enumerate(rows):
        status = row[4]
        if status in ('SKIPPED', 'ERROR'):
            continue
        value, limit = float(row[2]), float(row[3])
        margin = value / limit if limit > 0 else np.inf
        # A bar is drawn from the left edge of the axis, so margins outside the plotted
        # range are pulled just inside it to stay visible. The table has the exact
        # numbers.
        margin = min(max(margin, 3 * bar_min), 0.8 * bar_max)
        ax_bar.barh(i + 1.5, margin, height=0.6,
                    color='forestgreen' if status == 'PASS' else 'firebrick')

    ax_bar.axvline(1.0, color='black', linestyle='--', linewidth=1.2)
    ax_bar.set_xscale('log')
    ax_bar.set_xlim(bar_min, bar_max)
    ax_bar.set_ylim(n_rows + 1, 0)
    ax_bar.set_yticks([])
    ax_bar.set_xlabel('measured value / acceptance limit', fontsize=9)
    ax_bar.tick_params(axis='x', labelsize=8)
    ax_bar.grid(axis='x', linestyle=':', alpha=0.5)
    ax_bar.text(1.0, 0.4, ' limit', fontsize=8, va='center')

    n_pass = sum(r['status'] == 'PASS' for r in records)
    n_fail = sum(r['status'] == 'FAIL' for r in records)
    n_error = sum(r['status'] == 'ERROR' for r in records)
    total_time = sum(r['elapsed'] for r in records)
    all_good = n_fail == 0 and n_error == 0
    overall = 'ALL TESTS PASSED' if all_good else 'SOME TESTS FAILED'
    fig.suptitle(
        f"MagTense micromagnetism test suite - {overall}\n"
        f"{n_pass} passed, {n_fail} failed, {n_error} errored, "
        f"{len(skipped)} skipped, {total_time:.0f} s",
        fontsize=13, fontweight='bold',
        color='darkgreen' if (n_fail == 0 and n_error == 0) else 'darkred',
    )

    figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved overview figure to {figure_path}")


#%% Command line


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--list', action='store_true',
                        help="List the available tests and exit")
    parser.add_argument('--tests', type=str, default=None,
                        help="Comma separated subset of tests to run "
                             "(default: all but the slow ones)")
    parser.add_argument('--skip', type=str, default=None,
                        help="Comma separated tests to leave out. Used by the workflow, "
                             "where standard problem 6 is most of the running time")
    parser.add_argument('--include-slow', action='store_true',
                        help="Also run the tests marked slow")
    parser.add_argument('--cuda', action='store_true',
                        help="Ask the examples that support it to use the GPU")
    parser.add_argument('--no-figure', action='store_true',
                        help="Print the overview table but do not save the figure")
    args = parser.parse_args()

    if args.list:
        print(f"{'Test':38s} {'Slow':5s} Description")
        for name, (_, slow, description) in TESTS.items():
            print(f"{name:38s} {'yes' if slow else 'no':5s} {description}")
        return 0

    # A module level switch is the simplest way to reach the wrappers
    global USE_CUDA
    USE_CUDA = args.cuda

    if args.tests is not None:
        requested = [name.strip() for name in args.tests.split(',') if name.strip()]
        unknown = [name for name in requested if name not in TESTS]
        if unknown:
            print(f"Unknown test(s): {', '.join(unknown)}. "
                  "Use --list to see the available ones.")
            return 2
        selected, skipped = requested, []
    else:
        selected = [name for name, (_, slow, _) in TESTS.items()
                    if args.include_slow or not slow]
        skipped = [name for name in TESTS if name not in selected]

    if args.skip is not None:
        dropped = [name.strip() for name in args.skip.split(',') if name.strip()]
        unknown = [name for name in dropped if name not in TESTS]
        if unknown:
            print(f"Unknown test(s) in --skip: {', '.join(unknown)}. "
                  "Use --list to see the available ones.")
            return 2
        selected = [name for name in selected if name not in dropped]
        skipped = [name for name in TESTS if name not in selected]

    print(f"Running {len(selected)} test(s): {', '.join(selected)}")
    if skipped:
        print(f"Skipping {len(skipped)} test(s): {', '.join(skipped)}")

    records = run_tests(selected)
    print_table(records, skipped)

    if not args.no_figure:
        plot_overview(records, skipped,
                      RESULTS_DIR / 'testMagTenseFunctions_overview.png')

    failed = [r['name'] for r in records if r['status'] != 'PASS']
    if failed:
        print(f"\nFAILED: {', '.join(failed)}")
        return 1
    print("\nAll selected tests passed")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
