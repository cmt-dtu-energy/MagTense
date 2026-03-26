from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib import cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from scipy.sparse import issparse  # <-- NEW


def _get_field(grid_info: Any, name: str) -> Any:
    """Helper to access GridInfo fields as attribute or dict key."""
    if hasattr(grid_info, name):
        return getattr(grid_info, name)
    return grid_info[name]


def cartesian_unstructured_mesh_plot(
    pos: NDArray[np.float64],
    dims: NDArray[np.float64],
    grid_info: Any,
    i_in: Sequence[np.ndarray],
    fname_save: str | Path | None = None,
    fig: Figure | None = None,
    ax: Axes | None = None,
) -> Tuple[Figure, Axes]:
    """
    Plot Cartesian unstructured mesh in 3D.
    """

    # --- Figure / axis setup -------------------------------------------------
    if fig is None or ax is None:
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection="3d")
    else:
        ax.cla()

    ax.view_init(elev=30.0, azim=30.0)

    # --- Extract GridInfo fields ---------------------------------------------
    the_signs = _get_field(grid_info, "TheSigns")  # may be sparse or dense
    Xf: NDArray[np.float64] = np.asarray(_get_field(grid_info, "Xf"), dtype=float)
    Yf: NDArray[np.float64] = np.asarray(_get_field(grid_info, "Yf"), dtype=float)
    Zf: NDArray[np.float64] = np.asarray(_get_field(grid_info, "Zf"), dtype=float)
    fNormX: NDArray[np.float64] = np.asarray(_get_field(grid_info, "fNormX"), dtype=float)
    fNormY: NDArray[np.float64] = np.asarray(_get_field(grid_info, "fNormY"), dtype=float)
    fNormZ: NDArray[np.float64] = np.asarray(_get_field(grid_info, "fNormZ"), dtype=float)

    has_dims_f = hasattr(grid_info, "DimsF") or (
        isinstance(grid_info, Mapping) and "DimsF" in grid_info
    )
    DimsF: NDArray[np.float64] | None = (
        np.asarray(_get_field(grid_info, "DimsF"), dtype=float) if has_dims_f else None
    )

    # --- Determine which faces have a single contributing cell ---------------
    # MATLAB: iOne = sum(abs(TheSigns),1) == 1;
    if issparse(the_signs):
        # the_signs: (Ncells, Nfaces) sparse
        col_sums = np.array(np.abs(the_signs).sum(axis=0)).ravel()
    else:
        ts_dense = np.asarray(the_signs)
        col_sums = np.sum(np.abs(ts_dense), axis=0)

    i_one: NDArray[np.bool_] = col_sums == 1
    i_one_ind: NDArray[np.int_] = np.where(i_one)[0]

    # --- Colour setup --------------------------------------------------------
    n_groups = len(i_in)
    if n_groups > 1:
        base_cols = cm.hsv(np.linspace(0.0, 1.0, n_groups - 1, endpoint=False))[:, :3]
        cols = np.vstack([base_cols, np.array([[0.4, 0.4, 0.4]])])
    else:
        cols = None

    line_width = 1.5
    patches: list[Poly3DCollection] = []

    def _colour_for_n_multi_groups(n_idx: int) -> NDArray[np.float64]:
        assert cols is not None
        for j, arr in enumerate(i_in):
            # ensure 1D comparison
            if np.any(arr == n_idx):
                return cols[j]
        return np.array([0.4, 0.4, 0.4], dtype=float)

    # --- Main face loop ------------------------------------------------------
    for k in i_one_ind:
        # Find n such that TheSigns(n, k) ~= 0
        if issparse(the_signs):
            col = the_signs[:, k]              # sparse column (Ncells, 1)
            n_indices = col.nonzero()[0]       # row indices where != 0
        else:
            ts_dense = np.asarray(the_signs)
            n_indices = np.where(ts_dense[:, k] != 0)[0]

        if n_indices.size == 0:
            continue
        n_idx = int(n_indices[0])

        # Decide colour & style
        if n_groups > 1:
            this_col = _colour_for_n_multi_groups(n_idx)
            line_style = "-"
        else:
            data = i_in[0]
            if data.ndim >= 2 and data.shape[1] == 3:
                this_col = data[n_idx, :, 0] if data.ndim == 3 else data[n_idx, :]
                line_style = "none"
            else:
                max_val = float(np.max(data)) if np.max(data) != 0 else 1.0
                val = float(data[n_idx] / max_val)
                this_col = np.array([val, val, val], dtype=float)
                line_style = "-"

        verts_list: list[list[tuple[float, float, float]]] = []

        def add_quad(xs: NDArray[np.float64],
                     ys: NDArray[np.float64],
                     zs: NDArray[np.float64]) -> None:
            verts_list.append(list(zip(xs.tolist(), ys.tolist(), zs.tolist())))

        # X-normal faces
        if abs(fNormX[k]) > 0:
            if DimsF is not None:
                xs = np.full(4, Xf[k], dtype=float)
                ys = Yf[k] + (np.array([0.0, 0.0, 1.0, 1.0]) - 0.5) * DimsF[k, 1]
                zs = Zf[k] + (np.array([0.0, 1.0, 1.0, 0.0]) - 0.5) * DimsF[k, 2]
            else:
                xs = pos[n_idx, 0] + np.array([1.0, 1.0, 1.0, 1.0]) * fNormX[k] * dims[n_idx, 0] / 2.0
                ys = pos[n_idx, 1] + (np.array([0.0, 0.0, 1.0, 1.0]) - 0.5) * dims[n_idx, 1]
                zs = pos[n_idx, 2] + (np.array([0.0, 1.0, 1.0, 0.0]) - 0.5) * dims[n_idx, 2]
            add_quad(xs, ys, zs)

        # Y-normal faces
        if abs(fNormY[k]) > 0:
            if DimsF is not None:
                xs = Xf[k] + (np.array([0.0, 1.0, 1.0, 0.0]) - 0.5) * DimsF[k, 0]
                ys = np.full(4, Yf[k], dtype=float)
                zs = Zf[k] + (np.array([0.0, 0.0, 1.0, 1.0]) - 0.5) * DimsF[k, 2]
            else:
                xs = pos[n_idx, 0] + (np.array([0.0, 0.0, 1.0, 1.0]) - 0.5) * dims[n_idx, 0]
                ys = pos[n_idx, 1] + np.array([1.0, 1.0, 1.0, 1.0]) * fNormY[k] * dims[n_idx, 1] / 2.0
                zs = pos[n_idx, 2] + (np.array([0.0, 1.0, 1.0, 0.0]) - 0.5) * dims[n_idx, 2]
            add_quad(xs, ys, zs)

        # Z-normal faces
        if abs(fNormZ[k]) > 0:
            if DimsF is not None:
                xs = Xf[k] + (np.array([0.0, 0.0, 1.0, 1.0]) - 0.5) * DimsF[k, 0]
                ys = Yf[k] + (np.array([0.0, 1.0, 1.0, 0.0]) - 0.5) * DimsF[k, 1]
                zs = np.full(4, Zf[k], dtype=float)
            else:
                xs = pos[n_idx, 0] + (np.array([0.0, 1.0, 1.0, 0.0]) - 0.5) * dims[n_idx, 0]
                ys = pos[n_idx, 1] + (np.array([0.0, 0.0, 1.0, 1.0]) - 0.5) * dims[n_idx, 1]
                zs = pos[n_idx, 2] + np.array([1.0, 1.0, 1.0, 1.0]) * fNormZ[k] * dims[n_idx, 2] / 2.0
            add_quad(xs, ys, zs)

        if not verts_list:
            continue

        poly = Poly3DCollection(
            verts_list,
            linewidths=line_width,
            linestyles=line_style,
        )

        poly.set_facecolor(this_col)
        poly.set_edgecolor("black")       # ✔ Black outline around each quad
        poly.set_linewidth(line_width)    # ✔ Same thickness as MATLAB

        ax.add_collection3d(poly)
        patches.append(poly)

    # --- Axis limits & labels -------------------------------------------------
    dd = 0.0
    ax.set_xlim(float(np.min(Xf) - dd / 20.0), float(np.max(Xf) + dd / 20.0))
    ax.set_ylim(float(np.min(Yf) - dd / 20.0), float(np.max(Yf) + dd / 20.0))
    ax.set_zlim(float(np.min(Zf) - dd / 20.0), float(np.max(Zf) + dd / 20.0))

    ax.view_init(elev=30.0, azim=30.0)
    ax.set_box_aspect((1.0, 1.0, 1.0))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    # (GIF stuff omitted here for brevity; can be re-added unchanged.)

    return fig, ax
