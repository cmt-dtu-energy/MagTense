"""
@brief Compares the energy minimizer with the Landau-Lifshitz (LL) relaxation.

Two equilibrium problems are solved with both relaxation methods, and the cost is measured in
effective-field evaluations, which is the quantity that scales with the problem size (the
demagnetization product dominates every evaluation):

1. The hysteresis loop of a single uniaxial grain, a Stoner-Wohlfarth particle, swept with the
   adaptive field stepping. The switching field is known analytically, so both methods can be
   checked against the same answer.
2. mumag standard problem 3 at a single cube size: the flower and the vortex state are relaxed
   from their canonical starting states and the energies returned by the two methods are compared.

Run it from a MagTense source checkout with the python extension built. Pass ``--fast`` to reduce
the resolution of the standard problem 3 part.
"""

import argparse
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PYTHON_SOURCE = REPOSITORY_ROOT / "python" / "src"
for search_path in (REPOSITORY_ROOT, PYTHON_SOURCE):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

from magtense.micromag import MicromagProblem  # noqa: E402

MU0 = 4.0 * np.pi * 1.0e-7


def quiet(problem: MicromagProblem) -> None:
    problem.use_fmm = 0
    problem.log_dir = str(Path(tempfile.gettempdir()) / "magtense_minimizer_vs_llg")
    Path(problem.log_dir).mkdir(parents=True, exist_ok=True)
    problem.window_enabled = 0
    problem.trace_enabled = 0
    problem.setTimeDis = 1000


# ----------------------------------------------------------------------------------------------
# 1. Single grain hysteresis with adaptive field stepping
# ----------------------------------------------------------------------------------------------
def single_grain(solver: str, cuda: bool, n: int = 5, tilt_deg: float = 3.0, predictor: bool = False,
                 saddle_check: bool = True) -> dict:
    Bs, K0, A0 = 2.4, 1.0e6, 7.0e-12
    Ms = Bs / MU0
    L = 10.0e-9
    tilt = np.deg2rad(tilt_deg)
    field_direction = np.array([np.sin(tilt), 0.0, np.cos(tilt)])

    problem = MicromagProblem(
        res=[n, n, n],
        grid_L=[L, L, L],
        solver=solver,
        hysteresis_solver="adaptive",
        m0=np.tile(field_direction, (n**3, 1)),
        A0=A0, Ms=Ms, K0=K0,
        alpha=4000.0, gamma=0.0,
        cuda=cuda,
        usereturnhall=True,
        min_predictor=predictor,
        min_saddle_check=saddle_check,
    )
    problem.u_ea[:, :] = [0.0, 0.0, 1.0]
    # Two output times and a convergence check at the second one, so the LL relaxation can stop
    # early when the state no longer changes.
    problem.t = np.linspace(0.0, 1.0e-9, 2)
    problem.nt = 2
    problem.t_conv = problem.t.copy()
    problem.nt_conv = 2
    quiet(problem)

    t0 = time.perf_counter()
    result = problem.run_hysteresis_adaptive(
        H_start=1.0 / MU0 * field_direction,
        H_end=-2.0 / MU0 * field_direction,
        dH_initial=0.5 / MU0, dH_min=0.01 / MU0, dH_max=0.5 / MU0,
        max_steps=512,
        switch_refine_dH=0.01 / MU0,
    )
    wall = time.perf_counter() - t0

    n_acc = int(result[-1])
    H = np.asarray(result[4][0, 0, :n_acc, :]) @ field_direction  # A/m along the sweep direction
    m = np.asarray(result[1][-1, :, :n_acc, :]).mean(axis=0) @ field_direction
    # Switching field: the first field at which the mean projection changes sign
    i_sw = int(np.argmax(m < 0.0))
    H_sw = 0.5 * (H[i_sw - 1] + H[i_sw]) if i_sw > 0 else np.nan

    # Stoner-Wohlfarth switching field of a uniaxial particle: the cube's shape anisotropy is
    # zero by symmetry, so only K0 enters.  H_K = 2 K0 / (mu0 Ms); at angle psi
    # H_sw = H_K / (cos^(2/3) + sin^(2/3))^(3/2)
    HK = 2.0 * K0 / (MU0 * Ms)
    H_sw_SW = HK / (np.cos(tilt) ** (2.0 / 3.0) + np.sin(tilt) ** (2.0 / 3.0)) ** 1.5

    return dict(
        solver=solver, n_fields=n_acc, n_feval=int(problem.n_feval.sum()),
        n_feval_max=int(problem.n_feval.max()), wall=wall,
        H_sw_T=MU0 * H_sw, H_sw_SW_T=MU0 * H_sw_SW,
        n_fallback=int(np.sum(problem.min_status == 1)),
        n_fail=int(np.sum(problem.min_status == 2)),
        E_total=problem.E_out[-1, :n_acc, :].sum(axis=1),
        H_T=MU0 * H, m=m,
    )


# ----------------------------------------------------------------------------------------------
# 2. Standard problem 3 at one cube size
# ----------------------------------------------------------------------------------------------
def std_problem_3(solver: str, cuda: bool, res: int = 10, L_lex: float = 8.5) -> dict:
    A0 = 1.74532925199e-10
    Ms = 1e6
    K0 = 0.1 * 0.5 * MU0 * Ms**2
    lex = np.sqrt(A0 / (0.5 * MU0 * Ms**2))
    Km = 0.5 * MU0 * Ms**2

    out = {}
    for state in ("flower", "vortex"):
        problem = MicromagProblem(
            res=[res, res, res],
            grid_L=[lex * L_lex] * 3,
            solver=solver,
            A0=A0, Ms=Ms, K0=K0,
            alpha=1e3, gamma=0.0,
            cuda=cuda,
            usereturnhall=False,
        )
        problem.u_ea[:, 2] = 1
        quiet(problem)
        if state == "flower":
            problem.m0[:, 0:2] = 0
            problem.m0[:, 2] = 1
            t_end = 10e-9
        else:
            xv = np.linspace(-1, 1, res)
            x, y, z = np.meshgrid(xv, xv, xv, indexing="ij")
            m0 = np.zeros((res**3, 3))
            m0[:, 0] = np.sin(np.arctan2(z, x)).swapaxes(0, 2).reshape(-1)
            m0[:, 2] = -np.cos(np.arctan2(z, x)).swapaxes(0, 2).reshape(-1)
            problem.m0 = m0 / np.linalg.norm(m0, axis=1)[:, None]
            t_end = 200e-9

        # For the LL relaxation the run is a single constant field. A convergence check at every
        # output time lets it stop as soon as the magnetization is stationary.
        nt = 50
        problem.t_conv = np.linspace(0, t_end, nt)
        problem.nt_conv = nt
        problem.conv_tol = np.repeat(1e-6, nt)

        t0 = time.perf_counter()
        problem.run_simulation(
            t_end=t_end, nt=nt, fct_h_ext=lambda t: np.atleast_2d(t).T * np.array([0, 0, 0]),
            nt_h_ext=2,
        )
        wall = time.perf_counter() - t0
        V = (lex * L_lex) ** 3
        E = problem.E_out[-1, 0, :] / (Km * V)  # reduced energies: exc, ext, dem, ani
        out[state] = dict(
            E_red=E, E_tot=E.sum(), n_feval=int(problem.n_feval.sum()), wall=wall,
            iters=int(problem.min_iter.sum()), status=int(problem.min_status[0]),
            torque=float(problem.min_torque[0]),
        )
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cuda", action="store_true")
    parser.add_argument("--fast", action="store_true")
    parser.add_argument("--res3", type=int, default=10, help="cells per side for std problem 3")
    parser.add_argument("--ngrain", type=int, default=5, help="cells per side for the grain")
    args = parser.parse_args()
    res3 = 6 if args.fast else args.res3

    print("=" * 90)
    print("1. Single-grain hysteresis, adaptive field stepping")
    print("=" * 90)
    rows = []
    for solver in ("explicit_ll", "explicit"):
        r = single_grain(solver, args.cuda, n=args.ngrain)
        rows.append(r)
        print(
            f"{solver:>10s}: {r['n_fields']:4d} field steps, {r['n_feval']:7d} field evaluations "
            f"(max {r['n_feval_max']} per step), {r['wall']:7.2f} s, "
            f"mu0*H_sw = {r['H_sw_T']:.4f} T (Stoner-Wohlfarth {r['H_sw_SW_T']:.4f} T), "
            f"fallbacks {r['n_fallback']}, failures {r['n_fail']}"
        )
    # Compare the loops themselves on the common field axis
    Hgrid = np.linspace(min(rows[0]["H_T"].min(), rows[1]["H_T"].min()),
                        max(rows[0]["H_T"].max(), rows[1]["H_T"].max()), 400)
    m_ll = np.interp(Hgrid, rows[0]["H_T"][::-1], rows[0]["m"][::-1])
    m_mn = np.interp(Hgrid, rows[1]["H_T"][::-1], rows[1]["m"][::-1])
    print(f"   max |m_LL - m_min| over the loop: {np.max(np.abs(m_ll - m_mn)):.3e}")
    print(f"   field-evaluation ratio LL / minimizer: {rows[0]['n_feval'] / rows[1]['n_feval']:.1f}")

    print()
    print("=" * 90)
    print(f"2. Standard problem 3, {res3}^3 cells, L = 8.5 l_ex")
    print("=" * 90)
    res = {}
    for solver in ("explicit_ll", "explicit"):
        res[solver] = std_problem_3(solver, args.cuda, res=res3)
        for state in ("flower", "vortex"):
            r = res[solver][state]
            print(
                f"{solver:>10s} {state:>7s}: E/(Km V) = {r['E_tot']:.6f} "
                f"[exc {r['E_red'][0]:.5f} dem {r['E_red'][2]:.5f} ani {r['E_red'][3]:.5f}]  "
                f"{r['n_feval']:6d} field evaluations, {r['wall']:7.2f} s, "
                f"iters {r['iters']}, status {r['status']}, torque {r['torque']:.1e}"
            )
    for state in ("flower", "vortex"):
        dE = res["explicit"][state]["E_tot"] - res["explicit_ll"][state]["E_tot"]
        ratio = res["explicit_ll"][state]["n_feval"] / max(res["explicit"][state]["n_feval"], 1)
        print(f"   {state:>7s}: E_min - E_LL = {dE:+.2e} (reduced units), field-evaluation ratio LL / minimizer: {ratio:.1f}")


if __name__ == "__main__":
    main()
