"""Volume-averaged prism-prism demag tensor at all separations (tile type 8 with an observation size).

The averaged tensor of two prisms in closed form (Fukushima et al. 2002) is a corner sum whose terms
grow as the cube of the distance while the tensor falls as its inverse cube, so in double precision it
cancels for well separated tiles (Chernyshenko and Fangohr, JMMM 381, 440, 2015). MagTense evaluates it
in closed form, in units of the tile size, up to 1.5 summed vertex radii and as a Gauss average of the
point tensor of the source over the receiving tile beyond.

Two 5 nm cubes are placed from touching to 200 cells apart along x (where coordinate differences vanish
exactly) and along a generic direction. The reference is a 16-point-per-axis Gauss average of the point
tensor (tile type 2), which agrees with a 12-point rule to 2e-10 or better here. Check: the largest
relative Frobenius difference is below 1e-9 (the previous closed form reached 2.5 at 200 cells along x
and 7.7e-7 for touching cubes).

Usage:
    python avgprism_farfield_test.py

``run_test()`` returns a list of checks {'check', 'value', 'limit', 'passed'} like the other tests.
"""

from __future__ import annotations

import numpy as np

from magtense.magstatics import Tiles, get_demag_tensor

H = 5e-9


def _tile(tile_type: int) -> Tiles:
    return Tiles(n=1, tile_type=tile_type, offset=[0, 0, 0], size=[H, H, H], M_rem=0.0, mu_r_ea=1.0, mu_r_oa=1.0)


def averaged(c: np.ndarray) -> np.ndarray:
    return get_demag_tensor(_tile(8), np.asfortranarray(c[None]), obs_size=np.asfortranarray([[H, H, H]]))[0][0]


def reference(c: np.ndarray, q: int = 16) -> np.ndarray:
    x, w = np.polynomial.legendre.leggauss(q)
    x, w = 0.5 * H * x, 0.5 * w
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    W = np.einsum("i,j,k->ijk", w, w, w).ravel()
    pts = c + np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    return np.einsum("q,qij->ij", W, get_demag_tensor(_tile(2), np.asfortranarray(pts))[0])


def check(name: str, value: float, limit: float) -> dict:
    ok = bool(np.isfinite(value) and value < limit)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {value:.3e} (limit < {limit:g})")
    return {"check": name, "value": float(value), "limit": float(limit), "passed": ok}


def run_test(cuda: bool = False) -> list[dict]:
    checks = []
    for name, u in (("along x", np.array([1.0, 0.0, 0.0])),
                    ("generic direction", np.array([1.0, 0.6, 0.35]) / np.linalg.norm([1.0, 0.6, 0.35]))):
        worst = 0.0
        for d in (1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0):
            if name != "along x" and d < 2:
                continue
            c = d * H * u
            ref = reference(c)
            err = np.linalg.norm(averaged(c) - ref) / np.linalg.norm(ref)
            print(f"{name}, {d:5.0f} cells: {err:.2e}")
            worst = max(worst, err)
        checks.append(check(f"cubes {name}, 1-200 cells: largest relative difference", worst, 1e-9))
    return checks


if __name__ == "__main__":
    results = run_test()
    print("all passed" if all(c["passed"] for c in results) else "FAILED")
