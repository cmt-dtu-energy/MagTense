"""Exchange operator on tetrahedral meshes against an analytical exchange field.

A helix m = (cos kx, sin kx, 0) has the exact exchange field H_exc = -(2A/(mu0 Ms)) k^2 m, and it
satisfies the Neumann condition on the faces normal to y and z. MagTense builds the exchange
operator of a tetrahedral mesh itself (mesh analysis plus the Direct Laplacian with a quadratic
face fit, selected explicitly since the default for tetrahedra is now 'cellhessian'); the field it returns at t = 0 is compared with the exact one on the elements at least
10 nm from the two x ends of a 50 nm cube.

Meshes: the Kuhn mesh of magtense.utils.create_tetra_mesh (six tetrahedra per cube) and the same
mesh with its interior nodes moved at random by up to 0.2 of the cube edge (seed 1), which has no
axis aligned interior faces. Each at 10 and 5 nm.

Checks: the interior relative L2 error is below 0.15 at 5 nm on both meshes, and halving the
element size reduces it by at least a factor 1.6 (the scheme is first order on unstructured
tetrahedra). Before the quadratic fit the Kuhn mesh gave errors of order one that did not
decrease under refinement.

Usage:
    python tetra_exchange_test.py

``run_test()`` returns a list of checks {'check', 'value', 'limit', 'passed'} like the other tests.
"""

from __future__ import annotations

import numpy as np

from magtense.micromag import MicromagProblem
from magtense.utils import create_tetra_mesh

MU0 = 4e-7 * np.pi
MS = 8.0e5
A_EX = 1.3e-11
L = 50e-9
K = 2 * np.pi / 50e-9
MARGIN = 10e-9


def jitter(nodes: np.ndarray, h: float, amount: float = 0.2, seed: int = 1) -> np.ndarray:
    """Move every node that is not on the boundary of the cube by up to amount*h per coordinate."""
    rng = np.random.default_rng(seed)
    out = nodes.copy()
    interior = np.all(np.abs(nodes) < L / 2 - 1e-3 * h, axis=1)
    out[interior] += rng.uniform(-amount * h, amount * h, size=(interior.sum(), 3))
    return out


def helix_error(nodes: np.ndarray, elements: np.ndarray, cuda: bool = False) -> float:
    pts = nodes[elements - 1].mean(axis=1)
    m0 = np.stack([np.cos(K * pts[:, 0]), np.sin(K * pts[:, 0]), 0 * pts[:, 0]], axis=1)
    p = MicromagProblem(res=(len(pts), 1, 1), grid_L=[L, L, L], grid_type="tetrahedron",
                        grid_nod=nodes, grid_ele=elements.T, grid_pts=pts, m0=m0, Ms=MS, A0=A_EX,
                        K0=0.0, alpha=4000.0, gamma=0.0, solver="dynamic", cuda=cuda,
                        usedemag=False, usereturnhall=True, exch_presize=256,
                        exch_meth="directlaplacianneumann")
    p.window_enabled = 0
    res = p.run_simulation(t_end=1e-13, nt=2, fct_h_ext=lambda t: np.zeros((len(t), 3)), nt_h_ext=2)
    H = np.asarray(res[3])[0, :, 0, :]
    exact = -(2 * A_EX / (MU0 * MS)) * K * K * m0
    inner = np.abs(pts[:, 0]) <= L / 2 - MARGIN
    return float(np.linalg.norm(H[inner] - exact[inner]) / np.linalg.norm(exact[inner]))


def check(name: str, value: float, limit: float, below: bool = True) -> dict:
    ok = bool(np.isfinite(value) and (value < limit if below else value > limit))
    rel = "<" if below else ">"
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {value:.3e} (limit {rel} {limit:g})")
    return {"check": name, "value": float(value), "limit": float(limit), "passed": ok}


def run_test(cuda: bool = False) -> list[dict]:
    checks = []
    for name, perturb in (("Kuhn", False), ("jittered", True)):
        err = {}
        for h in (10e-9, 5e-9):
            nodes, elements, _ = create_tetra_mesh([L, L, L], h)
            if perturb:
                nodes = jitter(nodes, h)
            err[h] = helix_error(nodes, elements, cuda=cuda)
            print(f"{name} mesh, h = {h * 1e9:g} nm: interior relative L2 error {err[h]:.3e}")
        checks.append(check(f"{name} mesh helix error at 5 nm", err[5e-9], 0.15))
        checks.append(check(f"{name} mesh error ratio 10 nm / 5 nm", err[10e-9] / err[5e-9], 1.6, below=False))
    return checks


if __name__ == "__main__":
    results = run_test()
    print("all passed" if all(c["passed"] for c in results) else "FAILED")
