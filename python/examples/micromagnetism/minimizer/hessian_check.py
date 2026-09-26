"""
@brief Validates the Hessian eigenvalue check of the energy minimizer (min_saddle_check = 2).

The minimizer can end on a saddle point, where the torque vanishes as it does at a minimum. With
min_saddle_check = 2 the lowest eigenvalue of the energy Hessian in the tangent space is computed
by Lanczos with matrix-free products (one field evaluation each) and returned as min_eig, in
units of max(Ms). Three checks:

1. A single uniaxial cell, a Stoner-Wohlfarth particle, swept along a fixed field table. The
   demagnetizing field of a single cube is isotropic, so the lowest eigenvalue is the analytic
   in-plane curvature  lambda = H_K cos(2 theta) + H cos(theta - psi)  with theta the angle of m
   and psi the angle of the field from the easy axis, and H the signed field. It goes to zero at
   the switching field.
2. mumag standard problem 3, the vortex started from its canonical (symmetric) state, which is a
   saddle: the eigenvalue check must report a negative eigenvalue there, push off, and reach the
   same energy as the random nudge (mode 1). Mode 0 shows the saddle energy.
3. The single grain of minimizer_vs_llg.py on the adaptive loop: modes 1 and 2 must give the same
   loop, and the cost of the check is compared in field evaluations.

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
from minimizer_vs_llg import MU0, quiet, single_grain, std_problem_3  # noqa: E402


# ----------------------------------------------------------------------------------------------
# 1. Stoner-Wohlfarth cell against the analytic curvature
# ----------------------------------------------------------------------------------------------
def stoner_wohlfarth(cuda: bool, n_steps: int = 61, tilt_deg: float = 3.0, predictor: bool = True) -> dict:
    Bs, K0, A0 = 2.4, 1.0e6, 7.0e-12
    Ms = Bs / MU0
    L = 10.0e-9
    psi = np.deg2rad(tilt_deg)
    d = np.array([np.sin(psi), 0.0, np.cos(psi)])
    HK = 2.0 * K0 / (MU0 * Ms)

    problem = MicromagProblem(
        res=[1, 1, 1], grid_L=[L, L, L], solver="explicit", hysteresis_solver="static",
        m0=np.tile(d, (1, 1)), A0=A0, Ms=Ms, K0=K0, alpha=4000.0, gamma=0.0, cuda=cuda,
        usereturnhall=True, min_saddle_check=2, min_predictor=predictor,
    )
    problem.u_ea[:, :] = [0.0, 0.0, 1.0]
    problem.t = np.linspace(0.0, 1.0e-9, 2)
    problem.nt = 2
    problem.t_conv = problem.t.copy()
    problem.nt_conv = 2
    quiet(problem)

    B = np.linspace(1.0, -2.0, n_steps)
    h_ext = np.zeros((n_steps, 4))
    h_ext[:, 0] = np.arange(n_steps)
    h_ext[:, 1:4] = np.outer(B / MU0, d)
    result = problem.run_hysteresis(h_ext)
    m = np.asarray(result[1][-1, 0, :, :])  # (n_steps, 3)

    # angle of m from the easy axis in the (x, z) plane, signed like the field tilt
    theta = np.arctan2(m[:, 0], m[:, 2])
    H = B / MU0
    lam_analytic = (HK * np.cos(2.0 * theta) + H * np.cos(theta - psi)) / Ms
    lam = problem.min_eig
    err = np.abs(lam - lam_analytic)
    i_sw = int(np.argmax(m @ d < 0.0))
    return dict(
        B=B, lam=lam, lam_analytic=lam_analytic, err=err, i_sw=i_sw,
        n_feval=problem.n_feval.astype(int), theta=theta,
        H_sw_SW_T=MU0 * HK / (np.cos(psi) ** (2.0 / 3.0) + np.sin(psi) ** (2.0 / 3.0)) ** 1.5,
    )


def print_stoner_wohlfarth(r: dict) -> None:
    print("=" * 96)
    print("1. Stoner-Wohlfarth cell: lowest Hessian eigenvalue against the analytic curvature")
    print("=" * 96)
    print(f"   max |lambda - lambda_analytic| / Ms over the loop: {np.nanmax(r['err']):.2e}")
    print(f"   max relative error where |lambda| > 0.05 Ms: "
          f"{np.nanmax(r['err'][np.abs(r['lam_analytic']) > 0.05] / np.abs(r['lam_analytic'][np.abs(r['lam_analytic']) > 0.05])):.2e}")
    print(f"   field evaluations per field: mean {r['n_feval'].mean():.1f}, max {r['n_feval'].max()}")
    print(f"   Stoner-Wohlfarth switching field {r['H_sw_SW_T']:.4f} T; the eigenvalue before and across it:")
    i0 = max(r["i_sw"] - 6, 0)
    print("      mu0 H [T]   lambda/Ms   analytic     theta [deg]")
    for i in range(i0, min(r["i_sw"] + 2, len(r["B"]))):
        print(f"      {r['B'][i]:9.4f}   {r['lam'][i]:9.5f}   {r['lam_analytic'][i]:9.5f}   {np.rad2deg(r['theta'][i]):8.3f}")
    print()


# ----------------------------------------------------------------------------------------------
# 2. Standard problem 3 vortex: saddle detection
# ----------------------------------------------------------------------------------------------
def print_std_problem_3(cuda: bool, res: int) -> None:
    print("=" * 96)
    print(f"2. Standard problem 3, {res}^3 cells: the canonical vortex start is a saddle")
    print("=" * 96)
    rows = {}
    for mode in (0, 1, 2):
        t0 = time.perf_counter()
        r = std_problem_3("explicit", cuda, res=res, saddle_check=mode)
        rows[mode] = r
        for state in ("flower", "vortex"):
            s = r[state]
            eig = f"{s['eig']:+.4e}" if np.isfinite(s["eig"]) else "   n/a    "
            print(f"   mode {mode} {state:>7s}: E/(Km V) = {s['E_tot']:.6f}  {s['n_feval']:5d} field evaluations, "
                  f"status {s['status']}, lambda_min/Ms = {eig}")
    for state in ("flower", "vortex"):
        dE = rows[2][state]["E_tot"] - rows[1][state]["E_tot"]
        print(f"   {state:>7s}: E(mode 2) - E(mode 1) = {dE:+.2e}, field-evaluation ratio mode 1 / mode 2: "
              f"{rows[1][state]['n_feval'] / max(rows[2][state]['n_feval'], 1):.2f}")
    print()


# ----------------------------------------------------------------------------------------------
# 3. Single grain loop: cost of the two checks
# ----------------------------------------------------------------------------------------------
def print_grain(cuda: bool, n: int) -> None:
    print("=" * 96)
    print("3. Single grain, adaptive loop: nudge (mode 1) against eigenvalue (mode 2), predictor on")
    print("=" * 96)
    rows = {}
    for mode in (1, 2):
        rows[mode] = single_grain("explicit", cuda, n=n, predictor=True, saddle_check=mode)
        r = rows[mode]
        print(f"   mode {mode}: {r['n_fields']:4d} fields, {r['n_feval']:6d} field evaluations "
              f"({r['n_feval'] / r['n_fields']:5.1f} per field), mu0*H_sw = {r['H_sw_T']:.4f} T "
              f"(Stoner-Wohlfarth {r['H_sw_SW_T']:.4f} T), fallbacks {r['n_fallback']}, failures {r['n_fail']}")
    Hgrid = np.linspace(max(rows[1]["H_T"].min(), rows[2]["H_T"].min()),
                        min(rows[1]["H_T"].max(), rows[2]["H_T"].max()), 400)
    m1 = np.interp(Hgrid, rows[1]["H_T"][::-1], rows[1]["m"][::-1])
    m2 = np.interp(Hgrid, rows[2]["H_T"][::-1], rows[2]["m"][::-1])
    print(f"   max |m_2 - m_1| over the loop: {np.max(np.abs(m2 - m1)):.3e}")
    print(f"   field-evaluation ratio mode 1 / mode 2: {rows[1]['n_feval'] / max(rows[2]['n_feval'], 1):.2f}")
    lam = rows[2]["eig"]
    H = rows[2]["H_T"]
    print("   lowest eigenvalue / Ms along the descending branch (last fields before switching):")
    i_sw = int(np.argmax(rows[2]["m"] < 0.0))
    for i in range(max(i_sw - 6, 0), min(i_sw + 2, len(H))):
        print(f"      mu0 H = {H[i]:8.4f} T   lambda/Ms = {lam[i]:+.5f}")
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cuda", action="store_true")
    parser.add_argument("--res3", type=int, default=6, help="cells per side for std problem 3")
    parser.add_argument("--ngrain", type=int, default=5, help="cells per side for the grain")
    args = parser.parse_args()

    print_stoner_wohlfarth(stoner_wohlfarth(args.cuda))
    print_std_problem_3(args.cuda, args.res3)
    print_grain(args.cuda, args.ngrain)


if __name__ == "__main__":
    main()
