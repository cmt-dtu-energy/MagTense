"""The adaptive field loop's step metric and the minimizer short-circuit.

Every tile of a graded octree mesh (1632 prisms of 2, 4 and 8 nm, see dipfmm_octree_test.py) is an
independent Stoner-Wohlfarth particle: no exchange, no demagnetisation, uniaxial anisotropy along z
and the field swept from +0.5 T to -1.5 T at 1 degree from the axis. Two phases:

  * soft: the 1536 tiles of 2 nm, Ms / 30 and a low anisotropy field, switching near -0.36 T;
    94 % of the tiles but 2 % of the moment;
  * hard: the 4 and 8 nm tiles, switching near -1.14 T; 6 % of the tiles, 98 % of the moment.

The adaptive loop measures dM on the moment-weighted mean sum_i Ms_i V_i m_i / sum_i Ms_i V_i. The
soft switch moves that by about 0.04 < dM_reject = 0.05 and is taken in stride, while a plain
average over the tiles would move by 1.9 and be bisected down to dH_min. The hard switch is a real
switching event and is refined to dH_min.

Short-circuit: a sweep from -1.0 T (hard phase metastable) with a first step of 0.3 T jumps across
the hard switch at once. With the minimizer short-circuit on (dM_abort = 1.5) and off (0) it must
accept the same fields and states; with it on, the trials across the switch are abandoned early
(status 3, counted in the loop summary) and the sweep needs fewer field evaluations. A small
min_maxrot makes every descent take many iterations, as a large coupled system does.

``run_test()`` returns a list of checks {'check', 'value', 'limit', 'passed'} like the other tests.
"""

from __future__ import annotations

import os
import re
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path

import numpy as np

from magtense.micromag import MicromagProblem

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dipfmm_octree_test import L, octree_mesh  # noqa: E402

MU0 = 4e-7 * np.pi
MS_HARD = 1.0 / MU0            # mu0 Ms = 1 T
K_HARD = 5.0e5                 # mu0 H_K = 2 K / Ms = 1.26 T
MS_SOFT = MS_HARD / 30.0
K_SOFT = 0.25 * MS_SOFT * 0.8  # mu0 H_K = 0.4 T
TILT = np.deg2rad(1.0)          # nearly aligned: little rotation before a switch, so the
                                # controller keeps a large step and a trial jumps across it
DH_MIN_T = 0.01


@contextmanager
def captured_stdout():
    """Collect what the Fortran side writes to file descriptor 1."""
    sys.stdout.flush()
    saved = os.dup(1)
    with tempfile.TemporaryFile(mode="w+b") as tmp:
        os.dup2(tmp.fileno(), 1)
        out = {"text": ""}
        try:
            yield out
        finally:
            sys.stdout.flush()
            os.dup2(saved, 1)
            os.close(saved)
            tmp.seek(0)
            out["text"] = tmp.read().decode(errors="replace")


def summary_counts(text: str) -> dict:
    # list-directed Fortran output wraps at 80 columns with a leading blank on the continuation
    text = text.replace("\n ", "")
    m = re.search(r"Adaptive hysteresis summary: accepted\s+(\d+), rejected \(dM\)\s+(\d+), "
                  r"rejected \(short-circuit\)\s+(\d+), rejected \(sign change\)\s+(\d+), "
                  r"rejected \(not converged\)\s+(\d+)", text)
    if m is None:
        return {}
    keys = ("accepted", "dM", "short_circuit", "sign", "not_converged")
    return dict(zip(keys, map(int, m.groups())))


def sweep(dM_abort: float, h_start_T: float = 0.5, h_end_T: float = -1.5, dh_initial_T: float = 0.1,
          dh_max_T: float = 0.2) -> dict:
    pts, abc = octree_mesh()
    n = len(pts)
    soft = abc[:, 0] < 3e-9
    Ms = np.where(soft, MS_SOFT, MS_HARD)
    K0 = np.where(soft, K_SOFT, K_HARD)
    m0 = np.tile([0.0, 0.0, 1.0], (n, 1))
    # no exchange: A0 tiny rather than 0, which the prism-mesh exchange setup does not survive
    p = MicromagProblem(res=(n, 1, 1), grid_L=[L, L, L], grid_type="unstructuredPrisms", grid_pts=pts, grid_abc=abc,
                        m0=m0, Ms=Ms[:, None], A0=np.full((n, 1), 1e-25), K0=K0[:, None], alpha=4000.0, gamma=0.0, solver="explicit",
                        hysteresis_solver="adaptive", usedemag=False, exch_presize=64, min_saddle_check=0,
                        min_maxrot=0.02, min_predictor=False)
    p.u_ea[:, :] = [0.0, 0.0, 1.0]
    p.window_enabled = 0
    p.t = np.linspace(0.0, 1e-9, 2)
    p.nt = 2
    p.t_conv = np.linspace(0.0, 1e-9, 11)
    p.nt_conv = 11
    p.conv_tol = np.repeat(1e-6, 11)
    d = np.array([np.sin(TILT), 0.0, np.cos(TILT)])
    with captured_stdout() as out:
        res = p.run_hysteresis_adaptive(H_start=h_start_T / MU0 * d, H_end=h_end_T / MU0 * d,
                                        dH_initial=dh_initial_T / MU0, dH_min=DH_MIN_T / MU0,
                                        dH_max=dh_max_T / MU0, max_steps=200, dM_target=0.02,
                                        dM_reject=0.05, dH_grow=2.0, dH_shrink=0.5, dM_abort=dM_abort)
    m = np.asarray(res[1])[-1]                       # (n, fields, 3)
    fields = p.H_ext_applied @ d * MU0
    return {"fields": fields, "mz_soft": m[soft, :, 2].mean(axis=0), "mz_hard": m[~soft, :, 2].mean(axis=0),
            "m": m, "n_feval": np.asarray(p.n_feval, float), "status": np.asarray(p.min_status),
            "counts": summary_counts(out["text"]), "log": out["text"]}


def check(name: str, value: float, limit: float) -> dict:
    ok = bool(np.isfinite(value) and value < limit)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {value:.3e} (limit {limit:.0e})")
    return {"check": name, "value": float(value), "limit": float(limit), "passed": ok}


def step_before_switch(fields: np.ndarray, mz: np.ndarray) -> float:
    """Field step [T] of the first accepted field at which the phase is reversed (mz < 0)."""
    i = int(np.flatnonzero(mz < 0.0)[0])
    return float(fields[i - 1] - fields[i])


def run_test() -> list[dict]:
    checks = []
    full = sweep(1.5)
    print(f"+0.5 -> -1.5 T: {len(full['fields'])} fields, {full['counts']}")

    # moment-weighted metric: the soft switch is accepted at a regular step, the hard one refined
    soft_step = step_before_switch(full["fields"], full["mz_soft"])
    hard_step = step_before_switch(full["fields"], full["mz_hard"])
    checks.append(check("soft phase (2 % of the moment) switches without refinement: dH_min / step [T/T]",
                        DH_MIN_T / soft_step, 0.34))
    checks.append(check("hard phase switch refined to dH_min: step / dH_min - 1",
                        abs(hard_step / DH_MIN_T - 1.0), 0.05))
    checks.append(check("both phases reversed at the end (1 - <mz> reversed)",
                        1.0 - min(-full["mz_soft"][-1], -full["mz_hard"][-1]), 0.05))

    on = sweep(1.5, h_start_T=-1.0, dh_initial_T=0.3, dh_max_T=0.3)
    off = sweep(0.0, h_start_T=-1.0, dh_initial_T=0.3, dh_max_T=0.3)
    for tag, r in (("on ", on), ("off", off)):
        print(f"-1.0 -> -1.5 T, short-circuit {tag}: {len(r['fields'])} fields, {r['counts']}, "
              f"{r['n_feval'].sum():.0f} field evaluations")

    # short-circuit: same accepted fields and states, early rejections, fewer evaluations
    same_fields = len(on["fields"]) == len(off["fields"]) and np.allclose(on["fields"], off["fields"], atol=1e-9)
    checks.append(check("short-circuit on/off accept the same fields (0 = yes)", 0.0 if same_fields else 1.0, 0.5))
    if same_fields:
        checks.append(check("short-circuit on/off accept the same states (max |dm|)",
                            float(np.max(np.abs(on["m"] - off["m"]))), 1e-3))
    checks.append(check("short-circuit rejections with dM_abort = 1.5 (0 = none; value is 1/count)",
                        1.0 / max(on["counts"].get("short_circuit", 0), 1e-9), 1.01))
    checks.append(check("no short-circuit with dM_abort = 0 (count)", float(off["counts"].get("short_circuit", -1)), 0.5))
    checks.append(check("no accepted field carries the short-circuit status 3 (count)",
                        float(np.sum(on["status"] == 3)), 0.5))
    checks.append(check("field evaluations with / without the short-circuit",
                        on["n_feval"].sum() / off["n_feval"].sum(), 1.0))
    return checks


def main() -> int:
    checks = run_test()
    failed = [c for c in checks if not c["passed"]]
    print(f"\n{len(checks) - len(failed)}/{len(checks)} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
