"""dip-fmm on tetrahedral meshes against MagTense's dense tetrahedral demagnetisation tensor.

A 32 nm cube is cut into 8 nm cubes and every cube into six tetrahedra (magtense.utils
.create_tetra_mesh, 384 elements). With a random magnetisation the demagnetising field at t = 0 is
compared between

  * the dense MagTense tensor of exact volume-averaged tetrahedron pairs (the reference),
  * dip-fmm on its uniform tree at depth 1 (16 nm leaves) and 2 (8 nm leaves, one cube per leaf),
  * dip-fmm on its adaptive tree (capacity 6, at most 3 levels, root fixed to the cube),

on the CPU and, when available, on the GPU. A uniformly magnetised mesh checks the exact volume
average <H> = -M/3 of a cube however it is cut.

dip-fmm assigns elements to leaves by their centroids. When the leaves are much smaller than the
elements, elements reach far out of their leaves and the far field is no longer accurate (dip-fmm
warns about this); choose the depth so that a leaf is not much smaller than an element.

Usage:
    python dipfmm_tetra_test.py            # CPU and, if present, CUDA
    python dipfmm_tetra_test.py --cpu-only

``run_test()`` returns a list of checks {'check', 'value', 'limit', 'passed'} like the other tests.
"""

from __future__ import annotations

import argparse
import shutil

import numpy as np

from magtense.micromag import MicromagProblem
from magtense.utils import create_tetra_mesh, tetra_volumes

L = 32e-9          # cube side [m]
H = 8e-9           # cubes of the Kuhn split [m]
MS = 8e5
A_EX = 1.3e-11


def random_m0(n: int, seed: int = 7) -> np.ndarray:
    m = np.random.default_rng(seed).normal(size=(n, 3))
    return m / np.linalg.norm(m, axis=1)[:, None]


def demag_field(problem: MicromagProblem) -> np.ndarray:
    res = problem.run_simulation(t_end=1e-13, nt=2, fct_h_ext=lambda t: np.zeros((len(t), 3)), nt_h_ext=2)
    return np.asarray(res[5])[0, :, 0, :]          # H_dem at t = 0, (ntot, 3)


def tetra_problem(nodes, elements, pts, m0, cuda=False, **dipfmm) -> MicromagProblem:
    p = MicromagProblem(res=(len(pts), 1, 1), grid_L=[L, L, L], grid_type="tetrahedron", grid_nod=nodes,
                        grid_ele=elements.T, grid_pts=pts, m0=m0, Ms=MS, A0=A_EX, K0=0.0, alpha=4000.0,
                        gamma=0.0, solver="dynamic", cuda=cuda, usereturnhall=True, exch_presize=128, **dipfmm)
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
    nodes, elements, pts = create_tetra_mesh([L, L, L], H)
    vol = tetra_volumes(nodes, elements)
    n = len(pts)
    print(f"tetrahedral mesh: {n} elements, volume {vol.sum() / L**3:.6f} L^3")
    m0 = random_m0(n)
    dense = demag_field(tetra_problem(nodes, elements, pts, m0))
    backends = [False] + ([] if cpu_only or not shutil.which("nvidia-smi") else [True])
    for cuda in backends:
        tag = "cuda" if cuda else "cpu"
        for depth in (1, 2):
            for order, limit in ((6, 2e-3), (8, 1e-3)):
                fmm = demag_field(tetra_problem(nodes, elements, pts, m0, cuda=cuda, use_cdfmm=True,
                                                cdfmm_order=order, cdfmm_depth=depth))
                checks.append(check(f"{tag}: dip-fmm uniform tree (depth {depth}, order {order}) vs dense",
                                    rel_l2(fmm, dense), limit))
        ada = demag_field(tetra_problem(nodes, elements, pts, m0, cuda=cuda, use_cdfmm=True, cdfmm_order=6,
                                        cdfmm_tree="adaptive", cdfmm_capacity=6, cdfmm_max_depth=3,
                                        cdfmm_root=[0.0, 0.0, 0.0, 0.5 * L]))
        checks.append(check(f"{tag}: dip-fmm adaptive tree (capacity 6, root = cube, order 6) vs dense",
                            rel_l2(ada, dense), 2e-3))
        # uniform magnetisation: the volume-averaged field of the cube is exactly -M/3
        mu = np.tile([0.0, 0.0, 1.0], (n, 1))
        hu = demag_field(tetra_problem(nodes, elements, pts, mu, cuda=cuda, use_cdfmm=True, cdfmm_order=6,
                                       cdfmm_depth=2))
        avg = (vol[:, None] * hu).sum(axis=0) / vol.sum() / MS
        checks.append(check(f"{tag}: uniform M gives <H>/Ms = (0, 0, -1/3) (max deviation)",
                            float(np.max(np.abs(avg - [0.0, 0.0, -1.0 / 3.0]))), 1e-3))
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
