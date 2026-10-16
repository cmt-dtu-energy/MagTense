"""Symmetric exchange operator of unstructured grids (exch_meth='sushi'): symmetry, energy, accuracy.

The exchange field of 'sushi' is derived from a discrete exchange energy E(m) = m^T K m with K
symmetric positive semidefinite and K 1 = 0, and the matrix MagTense builds is L = -V^-1 K (V the
cell volumes). Then V L is symmetric, the field is exactly -dE/dm, and the exchange energy is
conserved by undamped and decreases under damped Landau-Lifshitz dynamics. The matrix is returned
as COO (A0 uniform) and checked directly, on the meshes of unstructured_exchange_test (Kuhn,
jittered (0.25 h) and Delaunay tetrahedra; uniform, graded and octree prisms; n = 5 and 10):

  symmetry      ||V L - (V L)^T|| / ||V L|| at round-off
  constants     L 1 = 0
  definiteness  no positive eigenvalue of V^1/2 L V^-1/2 beyond round-off (n = 5)
  accuracy      the solution of (k^2 V - V L) u = 2 k^2 V phi for an exact free-surface eigenfunction
                phi (so u = phi exactly) converges: error at n = 10 below the limit, and it falls
                from n = 5 to n = 10 (finite-volume operators are judged by the solution; the
                pointwise residual of L phi need not converge on irregular meshes)
  regular grid  on uniform prisms the matrix is exactly the 7-point Laplacian with mirror Neumann
                boundaries, i.e. the standard finite-difference operator
  periodic      exchPBC along x: symmetric, and the periodic field is reproduced
  cylinder      curved free surface
  energy        exchange-only LL dynamics (scipy RK45, rtol 1e-10) of a non-uniform state:
                alpha = 0 conserves E to the integrator tolerance, alpha > 0 decreases it at every
                output step; for contrast, the non-symmetric 'cellhessian' operator drifts.

``run_test()`` returns a list of checks {'check', 'value', 'limit', 'passed'} like the other tests.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.integrate import solve_ivp

from magtense.micromag import MicromagProblem
from unstructured_exchange_test import (L, delaunay_box, delaunay_cylinder, f_box, f_cyl, f_per, kuhn,
                                        prisms)

MS, A = 8.0e5, 1.3e-11


def exchange_matrix(grid, method="sushi", pbc=(0, 0, 0), cuda=False):
    if grid[0] == "tet":
        nodes, ele = grid[1], grid[2]
        cen = nodes[ele].mean(axis=1)
        v = nodes[ele]
        vol = np.abs(np.einsum("ij,ij->i", v[:, 1] - v[:, 0], np.cross(v[:, 2] - v[:, 0], v[:, 3] - v[:, 0]))) / 6
        ext = nodes.max(axis=0) - nodes.min(axis=0)
        kw = dict(grid_type="tetrahedron", grid_nod=nodes, grid_ele=(ele + 1).T, grid_pts=cen)
    else:
        cen, abc = grid[1], grid[2]
        vol = np.prod(abc, axis=1)
        ext = (cen + abc / 2).max(axis=0) - (cen - abc / 2).min(axis=0)
        kw = dict(grid_type="unstructuredPrisms", grid_pts=cen, grid_abc=abc)
    N = len(cen)
    p = MicromagProblem(res=(N, 1, 1), grid_L=list(ext), m0=1 / np.sqrt(3), Ms=MS, A0=A, K0=0.0,
                        alpha=4000.0, gamma=0.0, solver="dynamic", cuda=cuda, usedemag=False,
                        exch_presize=400, exch_meth=method, exchPBC=list(pbc), **kw)
    p.window_enabled = 0
    r = p.run_simulation(t_end=1e-13, nt=2, fct_h_ext=lambda t: np.zeros((len(t), 3)), nt_h_ext=2)
    nv = int(r[7])
    M = sp.csr_matrix((r[10][:nv], (r[8][:nv] - 1, r[9][:nv] - 1)), shape=(N, N))
    return M, cen, vol


def asym(M, vol):
    Am = sp.diags(vol) @ M
    return float(sp.linalg.norm(Am - Am.T) / sp.linalg.norm(Am))


def max_eig(M, vol):
    s = np.sqrt(vol)
    ev = np.linalg.eigvalsh((sp.diags(s) @ M @ sp.diags(1 / s)).toarray())
    return float(ev.max() / abs(ev.min()))


def helmholtz(M, cen, vol, f):
    phi, lap = f(cen)
    k2 = -np.sum(vol * phi * lap) / np.sum(vol * phi * phi)
    u = spla.spsolve((k2 * sp.diags(vol) - sp.diags(vol) @ M).tocsc(), 2 * k2 * vol * phi)
    return float(np.sqrt(np.sum(vol * (u - phi) ** 2) / np.sum(vol * phi ** 2)))


def fd_laplacian(cen, abc):
    """7-point Laplacian with mirror Neumann boundaries of a uniform prism grid."""
    h = abc[0, 0]
    N = len(cen)
    idx = {tuple(np.round(c / h * 2).astype(int)): i for i, c in enumerate(cen)}
    rows, cols = [], []
    for i, c in enumerate(cen):
        key = np.round(c / h * 2).astype(int)
        for ax in range(3):
            for s in (-2, 2):
                k = key.copy()
                k[ax] += s
                j = idx.get(tuple(k))
                if j is not None:
                    rows.append(i)
                    cols.append(j)
    Off = sp.csr_matrix((np.ones(len(rows)) / h ** 2, (rows, cols)), shape=(N, N))
    return (Off - sp.diags(np.asarray(Off.sum(axis=1)).ravel())).tocsr()


def ll_energy(M, cen, vol, alpha, t_end=20.0):
    """Exchange-only LL dynamics, h = L m (time in units of 1/(|L| h^2)); returns E(t)/E(0)."""
    rng = np.random.default_rng(7)
    x = (cen - cen.mean(axis=0)) / L * 2 * np.pi
    m = np.stack([np.cos(x[:, 0]) + 0.3 * np.sin(x[:, 2]), np.sin(x[:, 0]) * np.cos(x[:, 1]),
                  0.5 + 0.3 * np.cos(x[:, 1] + x[:, 2])], axis=1) + 0.05 * rng.normal(size=(len(cen), 3))
    m /= np.linalg.norm(m, axis=1, keepdims=True)
    Ls = (M * np.cbrt(vol.mean()) ** 2).tocsr()
    N = len(vol)

    def rhs(t, y):
        mm = y.reshape(N, 3)
        mxh = np.cross(mm, Ls @ mm)
        return (-mxh - alpha * np.cross(mm, mxh)).ravel()

    sol = solve_ivp(rhs, (0, t_end), m.ravel(), t_eval=np.linspace(0, t_end, 101), method="RK45",
                    rtol=1e-10, atol=1e-12)
    E = np.array([-0.5 * np.sum(vol[:, None] * y.reshape(N, 3) * (Ls @ y.reshape(N, 3))) for y in sol.y.T])
    return E / E[0]


def check(name, value, limit, below=True):
    ok = bool(np.isfinite(value) and (value < limit if below else value > limit))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {value:.3e} (limit {'<' if below else '>'} {limit:g})")
    return {"check": name, "value": float(value), "limit": float(limit), "passed": ok}


def run_test(cuda: bool = False) -> list[dict]:
    checks = []
    meshes = {
        "Kuhn tetrahedra": lambda n: ("tet", *kuhn(n)),
        "jittered tetrahedra": lambda n: ("tet", *kuhn(n, jitter=0.25)),
        "Delaunay tetrahedra": lambda n: ("tet", *delaunay_box(n)),
        "uniform prisms": lambda n: ("box", *prisms(n, "uniform")),
        "graded prisms": lambda n: ("box", *prisms(n, "graded")),
        "octree prisms": lambda n: ("box", *prisms(n, "octree")),
    }
    for name, make in meshes.items():
        err = {}
        for n in (5, 10):
            g = make(n)
            M, cen, vol = exchange_matrix(g, cuda=cuda)
            err[n] = helmholtz(M, cen, vol, f_box)
            if n == 5:
                ev = max_eig(M, vol)
            else:
                sym = asym(M, vol)
                const = float(np.abs(M @ np.ones(len(cen))).max() / np.abs(M.diagonal()).max())
                if name == "uniform prisms":
                    fd = fd_laplacian(g[1], g[2])
                    dfd = float(abs(M - fd).max() / abs(fd).max())
        print(f"{name}: solution error n=5 {err[5]:.3e}, n=10 {err[10]:.3e}")
        checks.append(check(f"{name}, symmetric (V L)", sym, 1e-13))
        checks.append(check(f"{name}, constant field annihilated", const, 1e-12))
        checks.append(check(f"{name}, no growing mode (max eig / |min eig|, n=5)", ev, 1e-10))
        checks.append(check(f"{name}, free surfaces, solution error at n=10", err[10], 0.05))
        checks.append(check(f"{name}, error ratio n=5/n=10", err[5] / err[10], 1.8, below=False))
        if name == "uniform prisms":
            checks.append(check("uniform prisms, identical to the 7-point Laplacian (mirror Neumann)", dfd, 1e-12))
    for name, make in (("Kuhn tetrahedra", meshes["Kuhn tetrahedra"]), ("octree prisms", meshes["octree prisms"])):
        M, cen, vol = exchange_matrix(make(10), pbc=(1, 0, 0), cuda=cuda)
        e = helmholtz(M, cen, vol, f_per)
        print(f"{name}: periodic field with exchPBC, solution error {e:.3e}")
        checks.append(check(f"{name}, exchPBC, symmetric", asym(M, vol), 1e-13))
        checks.append(check(f"{name}, exchPBC, periodic field solution error", e, 0.05))
    M, cen, vol = exchange_matrix(("tet", *delaunay_cylinder(5e-9)), cuda=cuda)
    checks.append(check("cylinder, curved free surface, symmetric", asym(M, vol), 1e-13))
    checks.append(check("cylinder, curved free surface, solution error at h=5 nm", helmholtz(M, cen, vol, f_cyl), 0.05))
    for name, make in (("jittered tetrahedra", meshes["jittered tetrahedra"]), ("octree prisms", meshes["octree prisms"])):
        g = make(5)
        M, cen, vol = exchange_matrix(g, cuda=cuda)
        e0 = ll_energy(M, cen, vol, 0.0)
        e1 = ll_energy(M, cen, vol, 0.5)
        Mh, _, _ = exchange_matrix(g, method="cellhessian", cuda=cuda)
        eh = ll_energy(Mh, cen, vol, 0.0)
        print(f"{name}: undamped energy drift sushi {np.abs(e0 - 1).max():.2e}, cellhessian {np.abs(eh - 1).max():.2e}")
        checks.append(check(f"{name}, alpha=0: exchange energy conserved", np.abs(e0 - 1).max(), 1e-6))
        checks.append(check(f"{name}, alpha=0.5: largest energy increase per step", np.diff(e1).max(), 1e-12))
        checks.append(check(f"{name}, alpha=0.5: energy decays", e1[-1], 0.5))
    return checks


if __name__ == "__main__":
    results = run_test()
    print("all passed" if all(c["passed"] for c in results) else "FAILED")
