"""Shared helpers for the magnetostatic validations against FEM (COMSOL).

Each Validation_field_* directory holds one script that sets up a single tile and
compares its field with a COMSOL simulation of the same geometry. Every script is a port
of the MATLAB script of the same name in matlab/examples/Magnetostatics, with the same
tile, the same evaluation points and the same error measure, so the two languages run
the same test. The reference data lives in documentation/examples_FEM_validation.

matplotlib is only imported when a figure is actually drawn, so the validations run
without it.
"""

from pathlib import Path

import numpy as np

from magtense.magstatics import Tiles, run_simulation

MU0 = 4 * np.pi * 1e-7

FEM_ROOT = Path(__file__).resolve().parents[3] / "documentation" / "examples_FEM_validation"


def load_fem(relative_path: str) -> np.ndarray:
    """Load a COMSOL text export from documentation/examples_FEM_validation."""
    return np.loadtxt(FEM_ROOT / relative_path, comments="%")


def _interp1_linear_extrap(x_g: np.ndarray, g: np.ndarray, x_f: np.ndarray) -> np.ndarray:
    """MATLAB's interp1(x_g, g, x_f, 'linear', 'extrap')."""
    order = np.argsort(x_g, kind="stable")
    x_g, g = x_g[order], g[order]
    out = np.interp(x_f, x_g, g)
    lo, hi = x_f < x_g[0], x_f > x_g[-1]
    out[lo] = g[0] + (x_f[lo] - x_g[0]) * (g[1] - g[0]) / (x_g[1] - x_g[0])
    out[hi] = g[-1] + (x_f[hi] - x_g[-1]) * (g[-1] - g[-2]) / (x_g[-1] - x_g[-2])
    return out


def calculate_relative_integral_error(x_f, f, x_g, g) -> float:
    """Port of matlab/util/calculate_relative_integral_error.m.

    g(x_g) is interpolated to x_f, points where g is NaN (which happens exactly on a tile
    surface) are left out of the interpolation, and the result is the integral of the
    absolute difference to f divided by the integral of |f|, in percent. It is insensitive
    to the handful of points on a tile surface, where the pointwise difference is large.
    """
    x_f, f = np.asarray(x_f, dtype=np.float64), np.asarray(f, dtype=np.float64)
    x_g, g = np.asarray(x_g, dtype=np.float64), np.asarray(g, dtype=np.float64)
    keep = ~np.isnan(g)
    g_interpolated = _interp1_linear_extrap(x_g[keep], g[keep], x_f)
    return float(np.trapezoid(np.abs(f - g_interpolated), x_f) / np.trapezoid(np.abs(f), x_f) * 100)


def line_points(offset, x, y=None, z=None) -> np.ndarray:
    """Points along x, then y, then z, each line through ``offset``.

    Mirrors the pts arrays the MATLAB validations build. y and z default to x.
    """
    x = np.asarray(x, dtype=np.float64)
    y = x if y is None else np.asarray(y, dtype=np.float64)
    z = x if z is None else np.asarray(z, dtype=np.float64)
    lines = []
    for axis, values in enumerate((x, y, z)):
        pts = np.tile(np.asarray(offset, dtype=np.float64), (len(values), 1))
        pts[:, axis] = values
        lines.append(pts)
    return np.asfortranarray(np.vstack(lines))


def field(tile: Tiles, pts: np.ndarray, obs_size: np.ndarray | None = None) -> np.ndarray:
    """The field H in A/m of ``tile`` at ``pts``."""
    _, H = run_simulation(tile, np.asfortranarray(pts), obs_size=obs_size)
    return np.asarray(H)


def report(title: str, errors: list[float], labels: str = "xyz") -> None:
    print(
        f"{title}: relative integrated error between MagTense and FEM is "
        + ", ".join(f"{a} = {e:.3g} %" for a, e in zip(labels, errors, strict=False))
    )


def plot_comparison(title: str, panels: list[dict], ylabel: str) -> None:
    """One subplot per panel; a panel is {'xlabel': str, 'curves': [(label, x, y, style)]}."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, len(panels), figsize=(5 * len(panels), 4), squeeze=False)
    fig.suptitle(f"{title} - MagTense vs. FEM")
    for axis, panel in zip(ax[0], panels, strict=True):
        for label, x, y, style in panel["curves"]:
            axis.plot(x, y, style, markersize=3, fillstyle="none", label=label)
        axis.set_xlabel(panel["xlabel"])
        axis.set_ylabel(ylabel)
        axis.legend(fontsize=8)
        axis.grid(True)
    fig.tight_layout()
    plt.show()


def norm_lines_validation(
    tile: Tiles,
    title: str,
    fem_files: list[str],
    x_mt: list[np.ndarray],
    offset,
    show_plot: bool,
) -> list[float]:
    """|mu0 H| along x, y and z through ``offset`` against three COMSOL line graphs.

    ``x_mt`` holds the MagTense evaluation coordinates of the three lines, which need not
    be the FEM points: MagTense is interpolated to the FEM points, as in MATLAB.
    """
    pts = line_points(offset, *x_mt)
    Hnorm = MU0 * np.linalg.norm(field(tile, pts), axis=1)
    bounds = np.cumsum([0] + [len(v) for v in x_mt])
    errors, panels = [], []
    for i, fname in enumerate(fem_files):
        fem = load_fem(fname)
        g = Hnorm[bounds[i] : bounds[i + 1]]
        errors.append(calculate_relative_integral_error(fem[:, 0], fem[:, 1], x_mt[i], g))
        panels.append({
            "xlabel": f"{'xyz'[i]} [m]",
            "curves": [("MagTense", x_mt[i], g, "r."), ("FEM", fem[:, 0], fem[:, 1], "bo")],
        })
    report(title, errors)
    if show_plot:
        plot_comparison(title, panels, r"$|\mu_0 H|$ [T]")
    return errors
