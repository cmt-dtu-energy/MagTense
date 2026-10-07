"""dip-fmm on octree (unstructuredPrisms) meshes against MagTense's dense demagnetisation tensor.

A 32 nm cube is split into 8 nm base tiles, and every tile within 3 nm of the plane z = 1.3 nm is
refined twice (8 -> 4 -> 2 nm), giving a graded mesh whose tiles lie on the dyadic grid of the cube.
With a random magnetisation the demagnetising field at t = 0 is compared between

  * the dense MagTense tensor (the reference),
  * dip-fmm on its uniform tree, two levels deep (leaves of 8 nm: every tile inside its leaf),
  * dip-fmm on its adaptive tree (capacity 32, at most 5 levels, root fixed to the cube),

on the CPU and, when available, on the GPU. A uniform 8^3 grid checks that the original uniform-grid
path is unchanged and that the adaptive tree also runs there. The run without usereturnhall checks
that the four field outputs are not allocated at full size and that the applied fields come back in
``H_ext_applied``.

Usage:
    python dipfmm_octree_test.py            # CPU and, if present, CUDA
    python dipfmm_octree_test.py --cpu-only

``run_test()`` returns a list of checks {'check', 'value', 'limit', 'passed'} like the other tests.
"""

from __future__ import annotations

import argparse
import shutil

import numpy as np

from magtense.micromag import MicromagProblem

L = 32e-9          # cube side [m]
BASE = 4           # base tiles per side (8 nm)
LEVELS = 2         # refinements around the plane (8 -> 4 -> 2 nm)
PLANE_Z = 1.3e-9   # refinement plane
BAND = 3e-9        # refine tiles reaching within this distance of the plane
MS = 8e5
A_EX = 1.3e-11


def octree_mesh() -> tuple[np.ndarray, np.ndarray]:
    """Tile centres and full side lengths of the graded cube, centred on the origin."""
    h0 = L / BASE
    g = -0.5 * L + (np.arange(BASE) + 0.5) * h0
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    centres = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    sizes = np.full(len(centres), h0)
    kept_c, kept_s = [], []
    children = np.array([[sx, sy, sz] for sz in (-0.25, 0.25) for sy in (-0.25, 0.25) for sx in (-0.25, 0.25)])
    for level in range(LEVELS + 1):
        near = np.abs(centres[:, 2] - PLANE_Z) - 0.5 * sizes <= BAND
        refine = near if level < LEVELS else np.zeros(len(centres), dtype=bool)
        kept_c.append(centres[~refine])
        kept_s.append(sizes[~refine])
        r = np.flatnonzero(refine)
        centres = (centres[r][:, None, :] + children[None] * sizes[r][:, None, None]).reshape(-1, 3)
        sizes = np.repeat(sizes[r] / 2, 8)
    centres = np.concatenate(kept_c)
    sizes = np.concatenate(kept_s)
    return centres, np.repeat(sizes[:, None], 3, axis=1)


def random_m0(n: int, seed: int = 7) -> np.ndarray:
    m = np.random.default_rng(seed).normal(size=(n, 3))
    return m / np.linalg.norm(m, axis=1)[:, None]


def demag_field(problem: MicromagProblem) -> np.ndarray:
    res = problem.run_simulation(t_end=1e-13, nt=2, fct_h_ext=lambda t: np.zeros((len(t), 3)), nt_h_ext=2)
    return np.asarray(res[5])[0, :, 0, :]          # H_dem at t = 0, (ntot, 3)


def octree_problem(pts, abc, m0, cuda=False, returnhall=True, **dipfmm) -> MicromagProblem:
    n = len(pts)
    p = MicromagProblem(res=(n, 1, 1), grid_L=[L, L, L], grid_type="unstructuredPrisms", grid_pts=pts, grid_abc=abc,
                        m0=m0, Ms=MS, A0=A_EX, K0=0.0, alpha=4000.0, gamma=0.0, solver="dynamic", cuda=cuda,
                        usereturnhall=returnhall, exch_presize=64, **dipfmm)
    p.window_enabled = 0
    return p


def uniform_problem(m0, cuda=False, **dipfmm) -> MicromagProblem:
    p = MicromagProblem(res=[8, 8, 8], grid_L=[L, L, L], grid_type="uniform", m0=m0, Ms=MS, A0=A_EX, K0=0.0,
                        alpha=4000.0, gamma=0.0, solver="dynamic", cuda=cuda, usereturnhall=True, **dipfmm)
    p.window_enabled = 0
    return p


def rel_l2(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def check(name: str, value: float, limit: float) -> dict:
    ok = bool(np.isfinite(value) and value < limit)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {value:.3e} (limit {limit:.0e})")
    return {"check": name, "value": float(value), "limit": float(limit), "passed": ok}


def run_test(cpu_only: bool = False) -> list[dict]:
    checks = []
    pts, abc = octree_mesh()
    n = len(pts)
    sizes, counts = np.unique(np.round(abc[:, 0] * 1e9, 6), return_counts=True)
    print(f"octree mesh: {n} prisms, sizes [nm] {dict(zip(sizes.tolist(), counts.tolist()))}")
    m0 = random_m0(n)
    dense = demag_field(octree_problem(pts, abc, m0))
    backends = [False] + ([] if cpu_only or not shutil.which("nvidia-smi") else [True])
    for cuda in backends:
        tag = "cuda" if cuda else "cpu"
        for order, limit in ((6, 2e-3), (8, 5e-4)):
            uni = demag_field(octree_problem(pts, abc, m0, cuda=cuda, use_cdfmm=True, cdfmm_order=order, cdfmm_depth=2))
            checks.append(check(f"octree {tag}: dip-fmm uniform tree (depth 2, order {order}) vs dense", rel_l2(uni, dense), limit))
            ada = demag_field(octree_problem(pts, abc, m0, cuda=cuda, use_cdfmm=True, cdfmm_order=order,
                                             cdfmm_tree="adaptive", cdfmm_capacity=32, cdfmm_max_depth=5,
                                             cdfmm_root=[0.0, 0.0, 0.0, 0.5 * L]))
            checks.append(check(f"octree {tag}: dip-fmm adaptive tree (capacity 32, order {order}) vs dense", rel_l2(ada, dense), limit))
        # no cdfmm_root: MagTense uses the bounding cube of the tiles, which keeps them in their leaves
        ada_auto = demag_field(octree_problem(pts, abc, m0, cuda=cuda, use_cdfmm=True, cdfmm_order=6, cdfmm_tree="adaptive"))
        checks.append(check(f"octree {tag}: dip-fmm adaptive tree, bounding-cube root (cdfmm_root=None), order 6 vs dense", rel_l2(ada_auto, dense), 5e-3))

    # uniform grid: the original path and the adaptive tree on equal cubes
    m0u = random_m0(512, seed=11)
    dense_u = demag_field(uniform_problem(m0u))
    for cuda in backends:
        tag = "cuda" if cuda else "cpu"
        uni_u = demag_field(uniform_problem(m0u, cuda=cuda, use_cdfmm=True, cdfmm_order=6, cdfmm_depth=2))
        checks.append(check(f"uniform {tag}: dip-fmm uniform tree vs dense", rel_l2(uni_u, dense_u), 2e-3))
        ada_u = demag_field(uniform_problem(m0u, cuda=cuda, use_cdfmm=True, cdfmm_order=6, cdfmm_tree="adaptive",
                                            cdfmm_root=[0.0, 0.0, 0.0, 0.5 * L]))
        checks.append(check(f"uniform {tag}: dip-fmm adaptive tree vs dense", rel_l2(ada_u, dense_u), 2e-3))

    # outputs without usereturnhall: small field arrays, applied fields still returned
    p = octree_problem(pts, abc, m0, returnhall=False, use_cdfmm=True, cdfmm_order=6, cdfmm_depth=2)
    p.solver = "explicit"                      # field steps relaxed by the energy minimizer
    p.t = np.linspace(0.0, 1e-10, 2)           # LL fallback window of the minimizer
    p.nt = 2
    h_applied = np.array([1e4, -2e4, 3e4])
    res = p.run_hysteresis(np.array([[0.0, *h_applied], [1.0, *(2 * h_applied)]]))
    sizes_ok = all(np.asarray(res[i]).size == 3 for i in (3, 4, 5, 6))
    checks.append(check("field outputs are (1,1,1,3) without usereturnhall (0 = yes)", 0.0 if sizes_ok else 1.0, 0.5))
    applied_err = float(np.max(np.abs(p.H_ext_applied - np.array([h_applied, 2 * h_applied]))))
    checks.append(check("H_ext_applied returns the applied fields [A/m]", applied_err, 1e-6))
    checks.append(check("magnetisation still returned in full (|m| - 1)", float(np.max(np.abs(
        np.linalg.norm(np.asarray(res[1])[-1, :, -1, :], axis=1) - 1.0))), 1e-6))
    return checks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cpu-only", action="store_true")
    args = parser.parse_args()
    checks = run_test(cpu_only=args.cpu_only)
    failed = [c for c in checks if not c["passed"]]
    print(f"\n{len(checks) - len(failed)}/{len(checks)} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
