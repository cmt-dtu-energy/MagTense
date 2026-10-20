"""Direct Laplacian exchange operator on unstructured prism meshes (graded and locally refined).

MagTense builds the exchange operator of an 'unstructuredPrisms' grid with the Direct Laplacian of
Poulsen et al., J. Magn. Magn. Mater. 551, 169093 (2022) when exch_meth='directlaplacianneumann': the
gradient on each face is a weighted linear least-squares fit of the surrounding cells, and the
Laplacian of a cell is the sum of those gradients over its faces. Such an operator is exact for a
linear phi on any mesh, and on a uniform grid it is the 7-point finite-difference Laplacian.

The exchange matrix L (H_exc = (2A/(mu0 Ms)) L m per component) is taken from MagTense and applied to
phi = x, y (exact Laplacian 0) and x^2 + y^2 (exact 4) on a 500 x 125 x 3 nm film, on the cells at
least two base cells from the film edges:
  uniform  40 x 10 cells
  graded   40 x 10 cells, x widths varying smoothly between 0.76 and 1.24 times the mean
  grains   a 40 x 10 base grid in which cells cut by the boundaries of nine Voronoi grains are split
           into four in x and y, twice (hanging faces)

Checks: L x and L y vanish on all three meshes (below 1e-7 of the largest matrix entry: the mesh analysis
rounds the cell bounds to 2^-21 of the smallest cell, which leaves ~1e-9 on the graded grid); L (x^2 + y^2) is 4 to round-off on the
uniform grid and to 1 % on the graded grid, and stays within 4 +- 5 on the grain mesh, where the first
order face fit leaves an O(1) error next to the hanging faces. Before the fix the graded grid lost the
x part of the Laplacian (the mesh analysis rounded the cell centres to half the smallest cell, so the
graded faces did not match), and on the grain mesh L x reached 200 and L (x^2 + y^2) -153 to 149 (the
face fit kept only the gradient component along the face normal).

Usage:
    python prism_directlaplacian_test.py

``run_test()`` returns a list of checks {'check', 'value', 'limit', 'passed'} like the other tests.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from magtense.micromag import MicromagProblem

L = np.array([500e-9, 125e-9, 3e-9])
NX, NY = 40, 10


def film(xb: np.ndarray, ny: int) -> tuple[np.ndarray, np.ndarray]:
    xc, wx = 0.5 * (xb[1:] + xb[:-1]), np.diff(xb)
    yc = (np.arange(ny) + 0.5) * L[1] / ny
    X, Y = np.meshgrid(xc, yc, indexing="ij")
    WX = np.meshgrid(wx, yc, indexing="ij")[0]
    pts = np.stack([X.ravel(), Y.ravel(), np.full(X.size, L[2] / 2)], axis=1)
    abc = np.stack([WX.ravel(), np.full(X.size, L[1] / ny), np.full(X.size, L[2])], axis=1)
    return pts, abc


def grain_mesh(n_grains: int = 9, n_ref: int = 2, seed: int = 2031) -> tuple[np.ndarray, np.ndarray]:
    seeds = np.random.default_rng(seed).uniform(0.1, 0.9, (n_grains, 2)) * L[:2]
    out = []

    def refine(c, a, level):
        d = np.linalg.norm(seeds - c, axis=1)
        i, j = np.argsort(d)[:2]
        dist = (d[j] ** 2 - d[i] ** 2) / (2 * np.linalg.norm(seeds[i] - seeds[j]))
        if level < n_ref and dist <= 0.5 * np.hypot(*a):
            for ox in (-0.25, 0.25):
                for oy in (-0.25, 0.25):
                    refine(c + np.array([ox * a[0], oy * a[1]]), a / 2, level + 1)
        else:
            out.append((c[0], c[1], a[0], a[1]))

    a0 = L[:2] / np.array([NX, NY])
    for ix in range(NX):
        for iy in range(NY):
            refine((np.array([ix, iy]) + 0.5) * a0, a0.copy(), 0)
    o = np.array(out)
    pts = np.column_stack([o[:, 0], o[:, 1], np.full(len(o), L[2] / 2)])
    abc = np.column_stack([o[:, 2], o[:, 3], np.full(len(o), L[2])])
    return pts, abc


def exchange_matrix(pts: np.ndarray, abc: np.ndarray, cuda: bool = False) -> sp.csr_matrix:
    p = MicromagProblem(res=(len(pts), 1, 1), grid_L=list(L), grid_type="unstructuredPrisms", grid_pts=pts,
                        grid_abc=abc, m0=1 / np.sqrt(3), Ms=8e5, A0=1.3e-11, K0=0.0, alpha=4000.0, gamma=0.0,
                        cuda=cuda, usedemag=False, exch_presize=256, exch_meth="directlaplacianneumann")
    p.window_enabled = 0
    res = p.run_simulation(t_end=1e-13, nt=2, fct_h_ext=lambda t: np.zeros((len(t), 3)), nt_h_ext=2)
    nval, rows, cols, vals = res[7:11]
    n = int(nval)
    return sp.csr_matrix((vals[:n], (rows[:n] - 1, cols[:n] - 1)), shape=(len(pts), len(pts)))


def check(name: str, value: float, limit: float) -> dict:
    ok = bool(np.isfinite(value) and value < limit)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {value:.3e} (limit < {limit:g})")
    return {"check": name, "value": float(value), "limit": float(limit), "passed": ok}


def run_test(cuda: bool = False) -> list[dict]:
    s = np.linspace(0, 1, NX + 1)
    xg = s + 0.12 * np.sin(2 * np.pi * s) / np.pi
    meshes = {"uniform": (film(s * L[0], NY), 1e-9), "graded": (film(xg / xg[-1] * L[0], NY), 0.04),
              "grains": (grain_mesh(), 5.0)}
    checks = []
    for name, ((pts, abc), quad_limit) in meshes.items():
        M = exchange_matrix(pts, abc, cuda=cuda) * L[0] ** 2       # Laplacian in units of 1/L0^2
        x, y = pts[:, 0] / L[0], pts[:, 1] / L[0]
        h = 2 * L[0] / NX
        inner = (pts[:, 0] > h) & (pts[:, 0] < L[0] - h) & (pts[:, 1] > h) & (pts[:, 1] < L[1] - h)
        scale = abs(M).max()
        lin = max(np.abs(M @ x)[inner].max(), np.abs(M @ y)[inner].max()) / scale
        quad = np.abs(M @ (x * x + y * y) - 4)[inner].max()
        print(f"{name}: {len(pts)} cells")
        checks.append(check(f"{name}: linear phi, max |L phi| / max |L|", lin, 1e-7))
        checks.append(check(f"{name}: max |L (x^2 + y^2) - 4|", quad, quad_limit))
    return checks


if __name__ == "__main__":
    results = run_test()
    print("all passed" if all(c["passed"] for c in results) else "FAILED")
