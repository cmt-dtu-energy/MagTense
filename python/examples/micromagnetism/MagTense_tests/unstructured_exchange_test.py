"""Exchange operator of unstructured grids (exch_meth='cellhessian') against exact Laplacians.

The exchange matrix MagTense builds for a grid (returned as COO; A0 uniform, so it is the
Laplacian with the boundary conditions of the problem) is applied to fields whose Laplacian is
known and which satisfy those boundary conditions exactly, so the error is measured on ALL cells,
the cells at the boundary included:

  box       50 nm cube, free surfaces:   phi = cos(pi X/L) cos(2 pi Y/L), X = x + L/2, Y = y + L/2
  periodic  50 nm cube, exchPBC along x: phi = sin(2 pi x/L) cos(pi Y/L)  (periodic, but NOT a
            free surface at the x faces; without exchPBC the operator must miss it there)
  cylinder  R = 25 nm, H = 30 nm, free surfaces: phi = J0(k r), k R = 3.8317 (first zero of J1),
            a curved boundary; the mesh is the Delaunay tetrahedralisation of points in the cylinder

Meshes (n = 5 and 10 cells per edge): Kuhn tetrahedra, Kuhn tetrahedra with the interior nodes
moved by up to 0.35 h (slivers; a few of them invert, which the cell Hessian fit tolerates),
Delaunay tetrahedra of random points; uniform prisms, prisms with
random widths 0.4-1.6 h (graded), and an octree-like refinement with hanging faces.

Checks: the error over all cells is below the limit at n = 10 and decreases from n = 5 to n = 10
(by at least a factor 1.8), on every mesh; the periodic field is reproduced with exchPBC and not
without it; the operator annihilates a constant field. The previous Direct Laplacian failed all
of these on the irregular meshes (zero exchange on the graded prisms, errors growing under
refinement on the Kuhn, jittered and octree meshes).

``run_test()`` returns a list of checks {'check', 'value', 'limit', 'passed'} like the other tests.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.spatial import Delaunay
from scipy.special import j0

from magtense.micromag import MicromagProblem

L = 50e-9
MS, A = 8.0e5, 1.3e-11


# ------------------------------------------------------------------ meshes
def kuhn(n, jitter=0.0, seed=1):
    h = L / n
    g = np.arange(n + 1) * h - L / 2
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    nodes = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    nid = np.arange(len(nodes)).reshape(n + 1, n + 1, n + 1)
    perms = ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0))
    ele = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for p in perms:
                    c = [i, j, k]
                    t = [nid[tuple(c)]]
                    for ax in p:
                        c[ax] += 1
                        t.append(nid[tuple(c)])
                    ele.append(t)
    if jitter:
        rng = np.random.default_rng(seed)
        inner = np.all(np.abs(nodes) < L / 2 - 1e-3 * h, axis=1)
        nodes = nodes.copy()
        nodes[inner] += rng.uniform(-jitter * h, jitter * h, size=(inner.sum(), 3))
    return nodes, np.array(ele)


def drop_slivers(nodes, ele, h):
    v = nodes[ele]
    vol = np.abs(np.einsum("ij,ij->i", v[:, 1] - v[:, 0], np.cross(v[:, 2] - v[:, 0], v[:, 3] - v[:, 0]))) / 6
    return ele[vol > 1e-9 * h ** 3]


def delaunay_box(n, seed=2):
    rng = np.random.default_rng(seed)
    a = L / 2
    pts = [rng.uniform(-a, a, size=(n ** 3, 3))]
    for ax in range(3):
        for s in (-a, a):
            p = rng.uniform(-a, a, size=(n ** 2, 3))
            p[:, ax] = s
            pts.append(p)
    for ax in range(3):
        o = [i for i in range(3) if i != ax]
        for s1 in (-a, a):
            for s2 in (-a, a):
                p = np.zeros((n, 3))
                p[:, ax] = rng.uniform(-a, a, n)
                p[:, o[0]], p[:, o[1]] = s1, s2
                pts.append(p)
    pts.append(np.array([[x, y, z] for x in (-a, a) for y in (-a, a) for z in (-a, a)]))
    nodes = np.vstack(pts)
    return nodes, drop_slivers(nodes, Delaunay(nodes).simplices, L / n)


def delaunay_cylinder(h, R=25e-9, H=30e-9, seed=5):
    """Convex cylinder: rings of boundary points on the side and the caps, random interior points."""
    rng = np.random.default_rng(seed)
    nz = max(2, int(round(H / h)) + 1)
    nth = max(8, int(round(2 * np.pi * R / h)))
    pts = []
    for z in np.linspace(-H / 2, H / 2, nz):
        th = np.linspace(0, 2 * np.pi, nth, endpoint=False) + rng.uniform(0, 2 * np.pi)
        pts.append(np.stack([R * np.cos(th), R * np.sin(th), np.full(nth, z)], axis=1))
    for z in (-H / 2, H / 2):
        k = int(np.pi * R * R / h ** 2)
        r = R * np.sqrt(rng.uniform(0, 0.95, k))
        th = rng.uniform(0, 2 * np.pi, k)
        pts.append(np.stack([r * np.cos(th), r * np.sin(th), np.full(k, z)], axis=1))
    k = int(np.pi * R * R * H / h ** 3)
    r = R * np.sqrt(rng.uniform(0, 0.95, k))
    th = rng.uniform(0, 2 * np.pi, k)
    pts.append(np.stack([r * np.cos(th), r * np.sin(th), rng.uniform(-0.95, 0.95, k) * H / 2], axis=1))
    nodes = np.vstack(pts)
    return nodes, drop_slivers(nodes, Delaunay(nodes).simplices, h)


def prisms(n, kind, seed=3):
    h = L / n
    rng = np.random.default_rng(seed)
    if kind == "uniform":
        e = [np.linspace(-L / 2, L / 2, n + 1)] * 3
    elif kind == "graded":
        e = []
        for _ in range(3):
            w = rng.uniform(0.4, 1.6, n)
            e.append(np.concatenate([[0], np.cumsum(w / w.sum() * L)]) - L / 2)
    if kind in ("uniform", "graded"):
        cx = [0.5 * (x[1:] + x[:-1]) for x in e]
        wx = [np.diff(x) for x in e]
        C = np.stack(np.meshgrid(*cx, indexing="ij"), axis=-1).reshape(-1, 3)
        W = np.stack(np.meshgrid(*wx, indexing="ij"), axis=-1).reshape(-1, 3)
        return C, W
    pts, abc = prisms(n, "uniform")
    rng = np.random.default_rng(4)
    for _ in range(2):
        split = (rng.random(len(pts)) < 0.4) & (abc[:, 0] <= abc[:, 0].min() * 1.0001)
        P, Q = [pts[~split]], [abc[~split]]
        for p, a in zip(pts[split], abc[split]):
            for d in np.array(np.meshgrid([-.25, .25], [-.25, .25], [-.25, .25])).reshape(3, -1).T:
                P.append((p + d * a)[None])
                Q.append((a / 2)[None])
        pts, abc = np.vstack(P), np.vstack(Q)
    return pts, abc


# ------------------------------------------------------------------ fields
def f_box(c):
    X, Y = c[:, 0] + L / 2, c[:, 1] + L / 2
    a, b = np.pi / L, 2 * np.pi / L
    phi = np.cos(a * X) * np.cos(b * Y)
    return phi, -(a * a + b * b) * phi


def f_per(c):
    a, b = 2 * np.pi / L, np.pi / L
    phi = np.sin(a * c[:, 0]) * np.cos(b * (c[:, 1] + L / 2))
    return phi, -(a * a + b * b) * phi


def f_cyl(c, R=25e-9):
    k = 3.8317059702 / R
    phi = j0(k * np.hypot(c[:, 0], c[:, 1]))
    return phi, -k * k * phi


# ------------------------------------------------------------------ MagTense
def exchange_matrix(grid, pbc=(0, 0, 0), cuda=False):
    if grid[0] == "tet":
        nodes, ele = grid[1], grid[2]
        cen = nodes[ele].mean(axis=1)
        ext = nodes.max(axis=0) - nodes.min(axis=0)
        kw = dict(grid_type="tetrahedron", grid_nod=nodes, grid_ele=(ele + 1).T, grid_pts=cen)
    else:
        cen, abc = grid[1], grid[2]
        ext = (cen + abc / 2).max(axis=0) - (cen - abc / 2).min(axis=0)
        kw = dict(grid_type="unstructuredPrisms", grid_pts=cen, grid_abc=abc)
    N = len(cen)
    p = MicromagProblem(res=(N, 1, 1), grid_L=list(ext), m0=1 / np.sqrt(3), Ms=MS, A0=A, K0=0.0,
                        alpha=4000.0, gamma=0.0, solver="dynamic", cuda=cuda, usedemag=False,
                        exch_presize=400, exch_meth="cellhessian", exchPBC=list(pbc), **kw)
    p.window_enabled = 0
    r = p.run_simulation(t_end=1e-13, nt=2, fct_h_ext=lambda t: np.zeros((len(t), 3)), nt_h_ext=2)
    nv = int(r[7])
    M = sp.csr_matrix((r[10][:nv], (r[8][:nv] - 1, r[9][:nv] - 1)), shape=(N, N))
    return M, cen


def rel_err(M, cen, f):
    phi, lap = f(cen)
    return float(np.linalg.norm(M @ phi - lap) / np.linalg.norm(lap))


def check(name, value, limit, below=True):
    ok = bool(np.isfinite(value) and (value < limit if below else value > limit))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {value:.3e} (limit {'<' if below else '>'} {limit:g})")
    return {"check": name, "value": float(value), "limit": float(limit), "passed": ok}


def run_test(cuda: bool = False) -> list[dict]:
    checks = []
    meshes = {
        "Kuhn tetrahedra": lambda n: ("tet", *kuhn(n)),
        "jittered tetrahedra": lambda n: ("tet", *kuhn(n, jitter=0.35)),
        "Delaunay tetrahedra": lambda n: ("tet", *delaunay_box(n)),
        "uniform prisms": lambda n: ("box", *prisms(n, "uniform")),
        "graded prisms": lambda n: ("box", *prisms(n, "graded")),
        "octree prisms": lambda n: ("box", *prisms(n, "octree")),
    }
    for name, make in meshes.items():
        err = {}
        for n in (5, 10):
            M, cen = exchange_matrix(make(n), cuda=cuda)
            err[n] = rel_err(M, cen, f_box)
            if n == 10:
                const = float(np.abs(M @ np.ones(len(cen))).max() / np.abs(M.diagonal()).max())
        print(f"{name}: free-surface error n=5 {err[5]:.3e}, n=10 {err[10]:.3e}")
        checks.append(check(f"{name}, free surfaces, error at n=10", err[10], 0.15))
        checks.append(check(f"{name}, error ratio n=5/n=10", err[5] / err[10], 1.8, below=False))
        checks.append(check(f"{name}, constant field annihilated", const, 1e-10))
    for name, make in (("Kuhn tetrahedra", meshes["Kuhn tetrahedra"]), ("octree prisms", meshes["octree prisms"])):
        g = make(10)
        M, cen = exchange_matrix(g, pbc=(1, 0, 0), cuda=cuda)
        e_pbc = rel_err(M, cen, f_per)
        M0, _ = exchange_matrix(g, cuda=cuda)
        e_free = rel_err(M0, cen, f_per)
        print(f"{name}: periodic field, with exchPBC {e_pbc:.3e}, without {e_free:.3e}")
        checks.append(check(f"{name}, periodic along x with exchPBC", e_pbc, 0.15))
        checks.append(check(f"{name}, periodic field missed without exchPBC", e_free, 3 * e_pbc, below=False))
    err = {}
    for h in (10e-9, 5e-9):
        M, cen = exchange_matrix(("tet", *delaunay_cylinder(h)), cuda=cuda)
        err[h] = rel_err(M, cen, f_cyl)
    print(f"cylinder (curved free surface): error h=10 nm {err[10e-9]:.3e}, h=5 nm {err[5e-9]:.3e}")
    checks.append(check("cylinder, curved free surface, error at h=5 nm", err[5e-9], 0.15))
    checks.append(check("cylinder, error ratio 10 nm / 5 nm", err[10e-9] / err[5e-9], 1.5, below=False))
    return checks


if __name__ == "__main__":
    results = run_test()
    print("all passed" if all(c["passed"] for c in results) else "FAILED")
