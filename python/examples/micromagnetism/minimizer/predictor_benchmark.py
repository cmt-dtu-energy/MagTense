"""
@brief Measures what the secant predictor (min_predictor) saves the energy minimizer.

Three hysteresis loops are run with the minimizer, each with the predictor off and on, and the
cost is compared in effective-field evaluations per applied field. The loops must agree: the
predictor only changes the starting state of each relaxation, so the switching fields and the
energies must be the same to the minimizer tolerance.

1. The single uniaxial grain of minimizer_vs_llg.py, swept with the adaptive field stepping
   (non-uniform steps and rejected trials exercise the step-ratio of the predictor).
2. The same grain swept with a fixed table of uniformly spaced fields.
3. A thin permalloy film without anisotropy, fixed table, where the equilibrium is a nonuniform
   state that rotates smoothly with the field between switching events.

Run it from a MagTense source checkout with the python extension built.
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPOSITORY_ROOT = HERE.parents[3]
PYTHON_SOURCE = REPOSITORY_ROOT / "python" / "src"
for search_path in (REPOSITORY_ROOT, PYTHON_SOURCE, HERE):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

from magtense.micromag import MicromagProblem  # noqa: E402
from minimizer_vs_llg import MU0, quiet, single_grain  # noqa: E402


def loop_stats(problem: MicromagProblem, n_fields: int, H_along: np.ndarray, m_along: np.ndarray) -> dict:
    i_sw = int(np.argmax(m_along < 0.0))
    H_sw = 0.5 * (H_along[i_sw - 1] + H_along[i_sw]) if i_sw > 0 else np.nan
    return dict(
        n_fields=n_fields,
        n_feval=int(problem.n_feval[:n_fields].sum()),
        n_feval_per_field=problem.n_feval[:n_fields].astype(int),
        iters=int(problem.min_iter[:n_fields].sum()),
        H_sw_T=MU0 * H_sw,
        n_fallback=int(np.sum(problem.min_status[:n_fields] == 1)),
        n_fail=int(np.sum(problem.min_status[:n_fields] == 2)),
        E_total=problem.E_out[-1, :n_fields, :].sum(axis=1),
        H_T=MU0 * H_along, m=m_along,
    )


# ----------------------------------------------------------------------------------------------
# 1. Single grain, adaptive stepping (reuses minimizer_vs_llg.single_grain)
# ----------------------------------------------------------------------------------------------
def grain_adaptive(cuda: bool, predictor: bool, n: int, saddle: bool = True) -> dict:
    return single_grain("explicit", cuda, n=n, predictor=predictor, saddle_check=saddle)


# ----------------------------------------------------------------------------------------------
# 2. Single grain, fixed table of uniformly spaced fields
# ----------------------------------------------------------------------------------------------
def grain_static(cuda: bool, predictor: bool, n: int, n_steps: int = 121, tilt_deg: float = 3.0, saddle: bool = True) -> dict:
    Bs, K0, A0 = 2.4, 1.0e6, 7.0e-12
    Ms = Bs / MU0
    L = 10.0e-9
    tilt = np.deg2rad(tilt_deg)
    d = np.array([np.sin(tilt), 0.0, np.cos(tilt)])

    problem = MicromagProblem(
        res=[n, n, n], grid_L=[L, L, L], solver="explicit", hysteresis_solver="static",
        m0=np.tile(d, (n**3, 1)), A0=A0, Ms=Ms, K0=K0, alpha=4000.0, gamma=0.0, cuda=cuda,
        usereturnhall=True, min_predictor=predictor, min_saddle_check=saddle,
    )
    problem.u_ea[:, :] = [0.0, 0.0, 1.0]
    problem.t = np.linspace(0.0, 1.0e-9, 2)
    problem.nt = 2
    problem.t_conv = problem.t.copy()
    problem.nt_conv = 2
    quiet(problem)

    B = np.linspace(1.0, -2.0, n_steps)  # T along d
    h_ext = np.zeros((n_steps, 4))
    h_ext[:, 0] = np.arange(n_steps)
    h_ext[:, 1:4] = np.outer(B / MU0, d)

    t0 = time.perf_counter()
    result = problem.run_hysteresis(h_ext)
    wall = time.perf_counter() - t0
    m = np.asarray(result[1][-1, :, :, :]).mean(axis=0) @ d
    out = loop_stats(problem, n_steps, B / MU0, m)
    out["wall"] = wall
    return out


# ----------------------------------------------------------------------------------------------
# 3. Thin film without anisotropy, fixed table
# ----------------------------------------------------------------------------------------------
def film_static(cuda: bool, predictor: bool, res=(20, 5, 1), n_steps: int = 81, tilt_deg: float = 1.0, saddle: bool = True) -> dict:
    Ms, A0 = 8.0e5, 1.3e-11
    L = [200.0e-9, 50.0e-9, 5.0e-9]
    tilt = np.deg2rad(tilt_deg)
    d = np.array([np.cos(tilt), np.sin(tilt), 0.0])
    ntot = int(np.prod(res))

    problem = MicromagProblem(
        res=list(res), grid_L=L, solver="explicit", hysteresis_solver="static",
        m0=np.tile(d, (ntot, 1)), A0=A0, Ms=Ms, K0=0.0, alpha=4000.0, gamma=0.0, cuda=cuda,
        usereturnhall=True, min_predictor=predictor, min_saddle_check=saddle,
    )
    problem.t = np.linspace(0.0, 1.0e-9, 2)
    problem.nt = 2
    problem.t_conv = problem.t.copy()
    problem.nt_conv = 2
    quiet(problem)

    B = np.linspace(0.1, -0.1, n_steps)  # T along d
    h_ext = np.zeros((n_steps, 4))
    h_ext[:, 0] = np.arange(n_steps)
    h_ext[:, 1:4] = np.outer(B / MU0, d)

    t0 = time.perf_counter()
    result = problem.run_hysteresis(h_ext)
    wall = time.perf_counter() - t0
    m = np.asarray(result[1][-1, :, :, :]).mean(axis=0) @ d
    out = loop_stats(problem, n_steps, B / MU0, m)
    out["wall"] = wall
    return out


def compare(name: str, off: dict, on: dict) -> None:
    print("=" * 96)
    print(name)
    print("=" * 96)
    for label, r in (("predictor off", off), ("predictor on", on)):
        n_max = r["n_feval_max"] if "n_feval_max" in r else int(r["n_feval_per_field"].max())
        print(
            f"{label:>14s}: {r['n_fields']:4d} fields, {r['n_feval']:7d} field evaluations "
            f"({r['n_feval'] / r['n_fields']:6.1f} per field, max {n_max:4d}), "
            f"{r['wall']:7.2f} s, mu0*H_sw = {r['H_sw_T']:.4f} T, "
            f"fallbacks {r['n_fallback']}, failures {r['n_fail']}"
        )
    Hgrid = np.linspace(max(off["H_T"].min(), on["H_T"].min()), min(off["H_T"].max(), on["H_T"].max()), 400)
    m_off = np.interp(Hgrid, off["H_T"][::-1], off["m"][::-1])
    m_on = np.interp(Hgrid, on["H_T"][::-1], on["m"][::-1])
    n_common = min(len(off["E_total"]), len(on["E_total"]))
    same_fields = np.allclose(off["H_T"][:n_common], on["H_T"][:n_common])
    dE = (np.max(np.abs(on["E_total"][:n_common] - off["E_total"][:n_common])) / np.max(np.abs(off["E_total"][:n_common]))
          if same_fields else np.nan)
    print(f"   max |m_on - m_off| over the loop: {np.max(np.abs(m_on - m_off)):.3e}")
    print(f"   max relative energy difference at common fields: {dE:.3e}" + ("" if same_fields else " (field tables differ)"))
    print(f"   field-evaluation ratio off / on: {off['n_feval'] / max(on['n_feval'], 1):.2f}")
    if "n_feval_per_field" in off and "n_feval_per_field" in on:
        print("   evaluations per field, off:", " ".join(f"{v:3d}" for v in off["n_feval_per_field"][:40]))
        print("   evaluations per field, on: ", " ".join(f"{v:3d}" for v in on["n_feval_per_field"][:40]))
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cuda", action="store_true")
    parser.add_argument("--ngrain", type=int, default=5, help="cells per side for the grain")
    parser.add_argument("--skip-film", action="store_true")
    parser.add_argument("--no-saddle", action="store_true", help="switch the saddle check off in every run")
    args = parser.parse_args()

    s = not args.no_saddle
    compare("1. Single grain, adaptive field stepping",
            grain_adaptive(args.cuda, False, args.ngrain, s), grain_adaptive(args.cuda, True, args.ngrain, s))
    compare("2. Single grain, fixed field table",
            grain_static(args.cuda, False, args.ngrain, saddle=s), grain_static(args.cuda, True, args.ngrain, saddle=s))
    if not args.skip_film:
        compare("3. Thin film, fixed field table",
                film_static(args.cuda, False, saddle=s), film_static(args.cuda, True, saddle=s))


if __name__ == "__main__":
    main()
