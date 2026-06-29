"""Magnetization comparison utilities for the single-grain notebooks.

This module contains the plotting and ipywidgets dashboard code that used to
live directly in compare_magnetization.ipynb.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 - registers 3D projection

RESULTS_DIR = Path("results")
RESULT_PATTERN = "single_grain_*.npz"

try:
    import ipywidgets as widgets
    from IPython.display import clear_output, display
except Exception:  # pragma: no cover - notebooks provide these dependencies.
    widgets = None
    clear_output = None
    display = None

output = widgets.Output() if widgets is not None else None

def load_scalar(data, key, default=np.nan):
    if key not in data.files:
        return default
    value = data[key]
    if np.asarray(value).shape == ():
        return value.item()
    return value


def find_result_file(
    size_factor: float,
    n: int,
    adaptive: bool = False,
    periodic: bool = False,
    adaptive_dh_min_t: float | None = None,
    use_fmm: bool = False,
    fmm_nterms: int | None = None,
    nlmax: int | None = None,
    results_dir: Path = RESULTS_DIR,
) -> Path:
    stem = f'single_grain_sf{size_factor:.2g}_n{n}'

    if periodic:
        stem += '_P'

    if use_fmm:
        if fmm_nterms is None or nlmax is None:
            raise ValueError('fmm_nterms and nlmax must be given when use_fmm=True')
        stem += f'_fmm_on_N{fmm_nterms}_L{nlmax}'

    if adaptive:
        if adaptive_dh_min_t is None:
            raise ValueError('adaptive_dh_min_t must be given for adaptive runs')
        stem += f'_A_FS{abs(adaptive_dh_min_t):.1e}'

    expected_name = f'{stem}.npz'
    path = results_dir / expected_name
    if not path.exists():
        raise FileNotFoundError(f'{path} not found')
    return path


def load_result_data(path: Path):
    with np.load(path, allow_pickle=True) as data:
        H = np.asarray(data['H_array'], dtype=float)
        Mx = np.asarray(data['Mx_array'], dtype=float)
        My = np.asarray(data['My_array'], dtype=float)
        Mz = np.asarray(data['Mz_array'], dtype=float)
        res = data['res'] if 'res' in data.files else None
        return {
            'path': path,
            'H': H,
            'Mx': Mx,
            'My': My,
            'Mz': Mz,
            'res': res,
            'metadata': {
                'size_factor': load_scalar(data, 'size_factor'),
                'n': load_scalar(data, 'n'),
                'adaptive': load_scalar(data, 'adaptive', False),
                'adaptive_dh_min_t': load_scalar(data, 'adaptive_dh_min_t', np.nan),
            },
        }


def find_nearest_field_index(H: np.ndarray, target_H: float) -> int:
    H = np.asarray(H, dtype=float).ravel()
    return int(np.argmin(np.abs(H - target_H)))


def _as_numeric_array(a: object, name: str = "array") -> np.ndarray:
    """
    Convert object-like arrays to a clean float NumPy array.
    """
    try:
        return np.asarray(a, dtype=float)
    except Exception as exc:
        raise ValueError(f"Could not convert {name} to a float array") from exc


def _extract_res_entry(res: object, entry_index: int, name: str) -> np.ndarray:
    """
    Extract an entry from the MagTense-style res object.

    Expected use:
        entry_index = 1 -> magnetisation
        entry_index = 2 -> tile positions
    """
    if res is None:
        raise ValueError("No res object available")

    res_obj = np.asarray(res, dtype=object)

    if res_obj.ndim == 0:
        raise ValueError(f"Unexpected scalar res object while extracting {name}")

    if len(res_obj) <= entry_index:
        raise ValueError(
            f"Could not extract {name}: res has length {len(res_obj)}, "
            f"but requested entry {entry_index}"
        )

    return np.asarray(res_obj[entry_index], dtype=object)


def _extract_vector_quantity_from_candidate(
    candidate: object,
    idx: int | None,
    *,
    name: str,
) -> np.ndarray:
    """
    Extract a vector-valued tile quantity with final dimension 3.

    Supports common layouts:

        (1, n_tiles, n_steps, 3)
        (n_tiles, n_steps, 3)
        (n_steps, n_tiles, 3)
        (n_tiles, 3)

    For magnetisation, idx should normally be given.
    For static positions, idx may be None.
    """
    candidate = np.asarray(candidate, dtype=object)

    if candidate.ndim == 4:
        if candidate.shape[-1] != 3:
            raise ValueError(f"{name}: expected final dimension 3, got shape {candidate.shape}")

        if idx is None:
            values = candidate[0, :, 0, :]
        else:
            values = candidate[0, :, idx, :]

    elif candidate.ndim == 3:
        if candidate.shape[-1] != 3:
            raise ValueError(f"{name}: expected final dimension 3, got shape {candidate.shape}")

        if idx is None:
            values = candidate[:, 0, :]
        else:
            # Prefer layout (n_tiles, n_steps, 3), matching your original extractor.
            if idx < candidate.shape[1]:
                values = candidate[:, idx, :]
            elif idx < candidate.shape[0]:
                values = candidate[idx, :, :]
            else:
                raise IndexError(
                    f"{name}: idx={idx} is out of bounds for candidate with shape {candidate.shape}"
                )

    elif candidate.ndim == 2:
        if candidate.shape[-1] != 3:
            raise ValueError(f"{name}: expected shape (n_tiles, 3), got {candidate.shape}")

        values = candidate

    else:
        raise ValueError(f"{name}: unsupported candidate ndim={candidate.ndim}, shape={candidate.shape}")

    values = _as_numeric_array(values, name=name)

    if values.ndim != 2 or values.shape[1] != 3:
        raise ValueError(f"{name}: expected extracted values with shape (n_tiles, 3), got {values.shape}")

    return values


def extract_tile_magnetization(res: object, idx: int) -> np.ndarray:
    """
    Extract tile-resolved magnetisation from data['res'].

    Parameters
    ----------
    res:
        Full data['res'] object.
    idx:
        Field/time index.

    Returns
    -------
    m_tiles:
        Array with shape (n_tiles, 3), columns are Mx, My, Mz.
    """
    candidate = _extract_res_entry(res, 1, "magnetisation")
    return _extract_vector_quantity_from_candidate(candidate, idx, name="magnetisation")


def extract_tile_positions(res: object, idx: int | None = None) -> np.ndarray:
    """
    Extract tile-centre positions from data['res'].

    Parameters
    ----------
    res:
        Full data['res'] object.
    idx:
        Optional field/time index. Use this if positions are indexed.

    Returns
    -------
    positions:
        Array with shape (n_tiles, 3), columns are x, y, z.
    """
    candidate = _extract_res_entry(res, 2, "tile positions")
    return _extract_vector_quantity_from_candidate(candidate, idx, name="tile positions")


def normalise_vectors(v: np.ndarray, eps: float = 1e-30) -> np.ndarray:
    """
    Row-wise vector normalisation.
    """
    v = np.asarray(v, dtype=float)

    if v.ndim != 2 or v.shape[1] != 3:
        raise ValueError(f"Expected array with shape (n, 3), got {v.shape}")

    norm = np.linalg.norm(v, axis=1)
    norm_safe = np.where(norm > eps, norm, 1.0)

    return v / norm_safe[:, None]


def angle_between_vectors_deg(
    m1: np.ndarray,
    m2: np.ndarray,
    eps: float = 1e-30,
) -> np.ndarray:
    """
    Tile-wise angle between two vector fields in degrees.
    """
    m1 = np.asarray(m1, dtype=float)
    m2 = np.asarray(m2, dtype=float)

    if m1.shape != m2.shape:
        raise ValueError(f"Shape mismatch: m1 has {m1.shape}, m2 has {m2.shape}")

    if m1.ndim != 2 or m1.shape[1] != 3:
        raise ValueError(f"Expected arrays with shape (n_tiles, 3), got {m1.shape}")

    n1 = np.linalg.norm(m1, axis=1)
    n2 = np.linalg.norm(m2, axis=1)

    denom = n1 * n2
    valid = denom > eps

    cosang = np.full(m1.shape[0], np.nan)
    cosang[valid] = np.sum(m1[valid] * m2[valid], axis=1) / denom[valid]
    cosang = np.clip(cosang, -1.0, 1.0)

    return np.degrees(np.arccos(cosang))


def _estimate_spatial_scale(pos: np.ndarray) -> float:
    """
    Estimate a sensible vector length scale from the spatial domain.
    """
    pos = np.asarray(pos, dtype=float)

    extent = np.nanmax(pos, axis=0) - np.nanmin(pos, axis=0)
    max_extent = np.nanmax(extent)

    if not np.isfinite(max_extent) or max_extent <= 0:
        return 1.0

    n_tiles = pos.shape[0]
    spacing = max_extent / max(n_tiles ** (1.0 / 3.0), 1.0)

    if not np.isfinite(spacing) or spacing <= 0:
        spacing = 0.05 * max_extent

    return spacing


def _set_axes_equal_3d(ax, pos: np.ndarray, pad_fraction: float = 0.05):
    """
    Set equal x/y/z limits for a 3D axis based on physical tile positions.
    """
    pos = np.asarray(pos, dtype=float)

    mins = np.nanmin(pos, axis=0)
    maxs = np.nanmax(pos, axis=0)

    centre = 0.5 * (mins + maxs)
    span = np.nanmax(maxs - mins)

    if not np.isfinite(span) or span <= 0:
        span = 1.0

    half = 0.5 * span * (1.0 + pad_fraction)

    ax.set_xlim(centre[0] - half, centre[0] + half)
    ax.set_ylim(centre[1] - half, centre[1] + half)
    ax.set_zlim(centre[2] - half, centre[2] + half)

    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass


def _plot_coloured_quiver_3d(
    ax,
    pos: np.ndarray,
    vec: np.ndarray,
    colour_values: np.ndarray,
    *,
    cmap,
    norm,
    length: float,
    linewidth: float = 0.8,
    alpha: float = 0.95,
    arrow_length_ratio: float = 0.35,
):
    """
    Plot a 3D vector field at physical tile positions.
    """
    pos = np.asarray(pos, dtype=float)
    vec = np.asarray(vec, dtype=float)
    colour_values = np.asarray(colour_values, dtype=float)

    colours = cmap(norm(colour_values))

    ax.quiver(
        pos[:, 0],
        pos[:, 1],
        pos[:, 2],
        vec[:, 0],
        vec[:, 1],
        vec[:, 2],
        length=length,
        normalize=False,
        colors=colours,
        linewidth=linewidth,
        alpha=alpha,
        arrow_length_ratio=arrow_length_ratio,
    )

def draw_magnetization_curve_and_orientation_3d_on_figure(
    fig,
    H: np.ndarray,
    Mx: np.ndarray,
    My: np.ndarray,
    Mz: np.ndarray,
    res: object,
    idx1: int,
    idx2: int,
    *,
    use_unit_vectors: bool = True,
    stride: int = 1,
    vector_length: float = 0.75,
    diff_vector_length: float = 1.0,
    elev: float = 25,
    azim: float = -60,
    cmap_m: str = "coolwarm",
    cmap_angle: str = "viridis",
    positions_idx: int | None = None,
    last_panel: str = "difference",
    plot_title: str | None = None,
):
    """
    Draw the full 2x3 magnetisation plot into an existing figure.

    This function does NOT create a new figure. It clears and redraws the
    existing figure, which prevents widget canvases from stacking in notebooks.
    """
    if stride < 1:
        raise ValueError(f"stride must be >= 1, got {stride}")

    if last_panel not in {"difference", "angle"}:
        raise ValueError(
            f'last_panel must be either "difference" or "angle", got {last_panel!r}'
        )

    H = np.asarray(H, dtype=float)
    Mx = np.asarray(Mx, dtype=float)
    My = np.asarray(My, dtype=float)
    Mz = np.asarray(Mz, dtype=float)

    m1 = extract_tile_magnetization(res, idx1)
    m2 = extract_tile_magnetization(res, idx2)

    if positions_idx is None:
        pos = extract_tile_positions(res, idx=None)
    else:
        pos = extract_tile_positions(res, idx=positions_idx)

    m1 = np.asarray(m1, dtype=float)
    m2 = np.asarray(m2, dtype=float)
    pos = np.asarray(pos, dtype=float)

    if m1.shape != m2.shape:
        raise ValueError(f"Shape mismatch: m1 has {m1.shape}, m2 has {m2.shape}")

    if m1.ndim != 2 or m1.shape[1] != 3:
        raise ValueError(f"Expected m1/m2 with shape (n_tiles, 3), got {m1.shape}")

    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError(f"Expected positions with shape (n_tiles, 3), got {pos.shape}")

    if pos.shape[0] != m1.shape[0]:
        raise ValueError(
            f"Position/magnetisation tile-count mismatch: "
            f"positions have {pos.shape[0]} tiles, magnetisation has {m1.shape[0]}"
        )

    if use_unit_vectors:
        q1 = normalise_vectors(m1)
        q2 = normalise_vectors(m2)

        vector_quantity = r"\hat{m}"
        vector_quantity_label = r"$\hat{m}$"

        colour_1 = q1[:, 2]
        colour_2 = q2[:, 2]
        m_norm = plt.Normalize(vmin=-1.0, vmax=1.0)
    else:
        q1 = m1
        q2 = m2

        vector_quantity = r"m"
        vector_quantity_label = r"$m$"

        colour_1 = q1[:, 2]
        colour_2 = q2[:, 2]

        max_abs_component = np.nanmax(np.abs(np.vstack([q1, q2])))
        if not np.isfinite(max_abs_component) or max_abs_component <= 0:
            max_abs_component = 1.0

        m_norm = plt.Normalize(
            vmin=-max_abs_component,
            vmax=max_abs_component,
        )

    dq = q2 - q1
    angle_deg = angle_between_vectors_deg(m1, m2)

    sl = slice(None, None, stride)

    pos_p = pos[sl]
    q1_p = q1[sl]
    q2_p = q2[sl]
    dq_p = dq[sl]

    colour_1_p = colour_1[sl]
    colour_2_p = colour_2[sl]
    angle_p = angle_deg[sl]

    spatial_scale = _estimate_spatial_scale(pos)

    q_length = vector_length * spatial_scale
    dq_length = diff_vector_length * spatial_scale

    cmap_m_obj = plt.get_cmap(cmap_m)
    cmap_angle_obj = plt.get_cmap(cmap_angle)

    finite_angles = angle_deg[np.isfinite(angle_deg)]

    if finite_angles.size > 0:
        angle_vmax = np.nanpercentile(finite_angles, 99)
        if not np.isfinite(angle_vmax) or angle_vmax <= 0:
            angle_vmax = 1.0
    else:
        angle_vmax = 1.0

    angle_norm = plt.Normalize(vmin=0.0, vmax=angle_vmax)

    mean_angle = np.nanmean(angle_deg)
    max_angle = np.nanmax(angle_deg)

    # This is the important part: clear existing axes/colorbars.
    fig.clf()

    if plot_title:
        # Keep the run information visible without letting a long filename overlap
        # the first row of axes.
        from textwrap import wrap

        wrapped_title = "\n".join(wrap(str(plot_title), width=120))
        fig.suptitle(wrapped_title, fontsize=11, fontweight="bold", y=0.985)

    # Use explicit margins instead of constrained_layout.  This avoids the
    # ipympl canvas becoming very wide and leaving a large blank region between
    # the figure and the sidebar.  The larger top margin also prevents the
    # suptitle from covering the top-row subplot titles.
    gs = fig.add_gridspec(
        nrows=2,
        ncols=3,
        height_ratios=[1.0, 1.38],
        left=0.055,
        right=0.985,
        bottom=0.075,
        top=0.875 if plot_title else 0.945,
        hspace=0.34,
        wspace=0.28,
    )

    ax_mx = fig.add_subplot(gs[0, 0])
    ax_my = fig.add_subplot(gs[0, 1])
    ax_mz = fig.add_subplot(gs[0, 2])

    ax_3d_1 = fig.add_subplot(gs[1, 0], projection="3d")
    ax_3d_2 = fig.add_subplot(gs[1, 1], projection="3d")
    ax_3d_3 = fig.add_subplot(gs[1, 2], projection="3d")

    for ax, data_arr, comp in zip(
        (ax_mx, ax_my, ax_mz),
        (Mx, My, Mz),
        ("Mx [A/m]", "My [A/m]", "Mz [A/m]"),
    ):
        ax.plot(H, data_arr, "-", linewidth=1)

        ax.plot(
            [H[idx1], H[idx2]],
            data_arr[[idx1, idx2]],
            "o",
            color="red",
            markersize=6,
        )

        ax.axvline(H[idx1], color="red", linestyle="--", linewidth=0.8, alpha=0.5)
        ax.axvline(H[idx2], color="red", linestyle="--", linewidth=0.8, alpha=0.5)

        ax.set_title(f"Mean {comp}", fontsize=9)
        ax.set_xlabel("field value [T]", fontsize=8)
        ax.set_ylabel(comp, fontsize=8)
        ax.tick_params(axis="both", labelsize=8)
        ax.grid(True, linestyle="--", linewidth=0.5)

    _plot_coloured_quiver_3d(
        ax_3d_1,
        pos_p,
        q1_p,
        colour_1_p,
        cmap=cmap_m_obj,
        norm=m_norm,
        length=q_length,
        linewidth=0.9,
        alpha=0.95,
    )

    ax_3d_1.set_title(
        f"Index {idx1}, H = {H[idx1]:.4g}\n"
        f"{vector_quantity_label} coloured by z-component"
    )
    ax_3d_1.set_xlabel("x")
    ax_3d_1.set_ylabel("y")
    ax_3d_1.set_zlabel("z")
    ax_3d_1.view_init(elev=elev, azim=azim)
    _set_axes_equal_3d(ax_3d_1, pos)

    _plot_coloured_quiver_3d(
        ax_3d_2,
        pos_p,
        q2_p,
        colour_2_p,
        cmap=cmap_m_obj,
        norm=m_norm,
        length=q_length,
        linewidth=0.9,
        alpha=0.95,
    )

    ax_3d_2.set_title(
        f"Index {idx2}, H = {H[idx2]:.4g}\n"
        f"{vector_quantity_label} coloured by z-component"
    )
    ax_3d_2.set_xlabel("x")
    ax_3d_2.set_ylabel("y")
    ax_3d_2.set_zlabel("z")
    ax_3d_2.view_init(elev=elev, azim=azim)
    _set_axes_equal_3d(ax_3d_2, pos)

    if last_panel == "difference":
        _plot_coloured_quiver_3d(
            ax_3d_3,
            pos_p,
            dq_p,
            angle_p,
            cmap=cmap_angle_obj,
            norm=angle_norm,
            length=dq_length,
            linewidth=1.1,
            alpha=0.95,
        )

        ax_3d_3.set_title(
            rf"Difference vectors: $\Delta {vector_quantity} = "
            rf"{vector_quantity}_2 - {vector_quantity}_1$"
        )

    elif last_panel == "angle":
        ax_3d_3.scatter(
            pos[:, 0],
            pos[:, 1],
            pos[:, 2],
            c=angle_deg,
            s=35,
            cmap=cmap_angle_obj,
            norm=angle_norm,
            depthshade=False,
        )

        ax_3d_3.set_title(
            f"Angular change in 3D space\n"
            f"mean = {mean_angle:.2f}°, max = {max_angle:.2f}°"
        )

    ax_3d_3.set_xlabel("x")
    ax_3d_3.set_ylabel("y")
    ax_3d_3.set_zlabel("z")
    ax_3d_3.view_init(elev=elev, azim=azim)
    _set_axes_equal_3d(ax_3d_3, pos)

    for ax in (ax_3d_1, ax_3d_2, ax_3d_3):
        ax.title.set_fontsize(9)
        ax.xaxis.label.set_size(8)
        ax.yaxis.label.set_size(8)
        ax.zaxis.label.set_size(8)
        ax.tick_params(axis="both", labelsize=7)

    sm_m = plt.cm.ScalarMappable(cmap=cmap_m_obj, norm=m_norm)
    sm_m.set_array([])

    cbar_m = fig.colorbar(
        sm_m,
        ax=[ax_3d_1, ax_3d_2],
        shrink=0.58,
        pad=0.025,
    )

    if use_unit_vectors:
        cbar_m.set_label(r"$\hat{m}_z$")
    else:
        cbar_m.set_label(r"$m_z$")

    sm_angle = plt.cm.ScalarMappable(cmap=cmap_angle_obj, norm=angle_norm)
    sm_angle.set_array([])

    cbar_angle = fig.colorbar(
        sm_angle,
        ax=ax_3d_3,
        shrink=0.58,
        pad=0.045,
    )
    cbar_angle.set_label("angle [deg]")

    axs = {
        "Mx": ax_mx,
        "My": ax_my,
        "Mz": ax_mz,
        "m_idx1": ax_3d_1,
        "m_idx2": ax_3d_2,
        "last": ax_3d_3,
    }

    fig.canvas.draw_idle()

    return axs


def plot_magnetization_curve_and_orientation_3d(
    H: np.ndarray,
    Mx: np.ndarray,
    My: np.ndarray,
    Mz: np.ndarray,
    res: object,
    idx1: int,
    idx2: int,
    *,
    use_unit_vectors: bool = True,
    stride: int = 1,
    vector_length: float = 0.75,
    diff_vector_length: float = 1.0,
    elev: float = 25,
    azim: float = -60,
    cmap_m: str = "coolwarm",
    cmap_angle: str = "viridis",
    positions_idx: int | None = None,
    last_panel: str = "difference",
    plot_title: str | None = None,
):
    """
    Non-interactive convenience wrapper.

    Creates a new figure once and draws the 2x3 layout into it.
    """
    fig = plt.figure(figsize=(15.8, 8.2), constrained_layout=False)

    axs = draw_magnetization_curve_and_orientation_3d_on_figure(
        fig,
        H,
        Mx,
        My,
        Mz,
        res,
        idx1,
        idx2,
        use_unit_vectors=use_unit_vectors,
        stride=stride,
        vector_length=vector_length,
        diff_vector_length=diff_vector_length,
        elev=elev,
        azim=azim,
        cmap_m=cmap_m,
        cmap_angle=cmap_angle,
        positions_idx=positions_idx,
        last_panel=last_panel,
        plot_title=plot_title,
    )

    return fig, axs


def set_3d_view(axs: dict, elev: float = 25, azim: float = -60):
    """
    Set the same 3D view on all 3D axes in an axes dictionary.
    """
    for ax in axs.values():
        if hasattr(ax, "view_init") and hasattr(ax, "get_zlim"):
            ax.view_init(elev=elev, azim=azim)

    fig = next(iter(axs.values())).figure
    fig.canvas.draw_idle()


def add_view_sliders(
    fig,
    axs: dict,
    *,
    elev0: float = 25,
    azim0: float = -60,
    continuous_update: bool = True,
):
    """
    Add simple view sliders for a static figure.
    """
    elev_slider = widgets.FloatSlider(
        value=elev0,
        min=-90,
        max=90,
        step=1,
        description="elev",
        continuous_update=continuous_update,
        readout_format=".0f",
    )

    azim_slider = widgets.FloatSlider(
        value=azim0,
        min=-180,
        max=180,
        step=1,
        description="azim",
        continuous_update=continuous_update,
        readout_format=".0f",
    )

    def update_view(change=None):
        elev = elev_slider.value
        azim = azim_slider.value

        for ax in axs.values():
            if hasattr(ax, "view_init") and hasattr(ax, "get_zlim"):
                ax.view_init(elev=elev, azim=azim)

        fig.canvas.draw_idle()

    elev_slider.observe(update_view, names="value")
    azim_slider.observe(update_view, names="value")

    controls = widgets.HBox([elev_slider, azim_slider])
    display(controls)

    update_view()

    return controls


def interactive_magnetization_curve_and_orientation_3d(
    data: dict,
    *,
    idx1_init: int,
    idx2_init: int,
    use_unit_vectors: bool = True,
    stride_init: int = 1,
    vector_length_init: float = 0.75,
    diff_vector_length_init: float = 1.0,
    elev_init: float = 25,
    azim_init: float = -60,
    last_panel_init: str = "difference",
    positions_idx: int | None = None,
    continuous_update: bool = False,
):
    """
    Fully interactive version.

    Creates one persistent figure and redraws into that same figure when
    sliders change. This avoids stacked notebook widget canvases.
    """
    H = np.asarray(data["H"], dtype=float)
    n_steps = H.size

    if n_steps < 2:
        raise ValueError(f"Need at least two field/time steps, got {n_steps}")

    idx1_init = int(np.clip(idx1_init, 0, n_steps - 1))
    idx2_init = int(np.clip(idx2_init, 0, n_steps - 1))

    if last_panel_init not in {"difference", "angle"}:
        raise ValueError(
            f'last_panel_init must be "difference" or "angle", got {last_panel_init!r}'
        )

    idx1_slider = widgets.IntSlider(
        value=idx1_init,
        min=0,
        max=n_steps - 1,
        step=1,
        description="idx1",
        continuous_update=continuous_update,
        readout=True,
    )

    idx2_slider = widgets.IntSlider(
        value=idx2_init,
        min=0,
        max=n_steps - 1,
        step=1,
        description="idx2",
        continuous_update=continuous_update,
        readout=True,
    )

    elev_slider = widgets.FloatSlider(
        value=elev_init,
        min=-90,
        max=90,
        step=1,
        description="elev",
        continuous_update=continuous_update,
        readout_format=".0f",
    )

    azim_slider = widgets.FloatSlider(
        value=azim_init,
        min=-180,
        max=180,
        step=1,
        description="azim",
        continuous_update=continuous_update,
        readout_format=".0f",
    )

    stride_slider = widgets.IntSlider(
        value=stride_init,
        min=1,
        max=20,
        step=1,
        description="stride",
        continuous_update=continuous_update,
    )

    vector_length_slider = widgets.FloatSlider(
        value=vector_length_init,
        min=0.05,
        max=3.0,
        step=0.05,
        description="vec len",
        continuous_update=continuous_update,
        readout_format=".2f",
    )

    diff_vector_length_slider = widgets.FloatSlider(
        value=diff_vector_length_init,
        min=0.05,
        max=3.0,
        step=0.05,
        description="diff len",
        continuous_update=continuous_update,
        readout_format=".2f",
    )

    last_panel_dropdown = widgets.Dropdown(
        options=[
            ("Difference vectors", "difference"),
            ("Angular change", "angle"),
        ],
        value=last_panel_init,
        description="last",
    )

    index_info = widgets.HTML()

    controls_indices = widgets.HBox(
        [
            idx1_slider,
            idx2_slider,
            last_panel_dropdown,
        ]
    )

    controls_view = widgets.HBox(
        [
            elev_slider,
            azim_slider,
        ]
    )

    controls_style = widgets.HBox(
        [
            stride_slider,
            vector_length_slider,
            diff_vector_length_slider,
        ]
    )

    controls = widgets.VBox(
        [
            controls_indices,
            index_info,
            controls_view,
            controls_style,
        ]
    )

    display(controls)

    # Create one persistent figure.
    fig = plt.figure(figsize=(15.8, 8.2), constrained_layout=False)
    display(fig.canvas)

    state = {
        "fig": fig,
        "axs": None,
        "is_drawing": False,
    }

    def update_index_info():
        i1 = idx1_slider.value
        i2 = idx2_slider.value

        index_info.value = (
            f"<b>idx1:</b> {i1}, H = {H[i1]:.6g} &nbsp;&nbsp; "
            f"<b>idx2:</b> {i2}, H = {H[i2]:.6g} &nbsp;&nbsp; "
            f"<b>ΔH:</b> {H[i2] - H[i1]:.6g}"
        )

    def redraw(change=None):
        if state["is_drawing"]:
            return

        state["is_drawing"] = True

        try:
            update_index_info()

            axs = draw_magnetization_curve_and_orientation_3d_on_figure(
                fig,
                data["H"],
                data["Mx"],
                data["My"],
                data["Mz"],
                data["res"],
                idx1_slider.value,
                idx2_slider.value,
                use_unit_vectors=use_unit_vectors,
                stride=stride_slider.value,
                vector_length=vector_length_slider.value,
                diff_vector_length=diff_vector_length_slider.value,
                elev=elev_slider.value,
                azim=azim_slider.value,
                last_panel=last_panel_dropdown.value,
                positions_idx=positions_idx,
            )

            state["axs"] = axs

        finally:
            state["is_drawing"] = False

    for widget in [
        idx1_slider,
        idx2_slider,
        elev_slider,
        azim_slider,
        stride_slider,
        vector_length_slider,
        diff_vector_length_slider,
        last_panel_dropdown,
    ]:
        widget.observe(redraw, names="value")

    redraw()

    return {
        "idx1": idx1_slider,
        "idx2": idx2_slider,
        "elev": elev_slider,
        "azim": azim_slider,
        "stride": stride_slider,
        "vector_length": vector_length_slider,
        "diff_vector_length": diff_vector_length_slider,
        "last_panel": last_panel_dropdown,
        "index_info": index_info,
        "state": state,
    }

def build_result_file_path(
    size_factor: float,
    n: int,
    adaptive: bool = False,
    periodic: bool = False,
    adaptive_dh_min_t: float | None = None,
    use_fmm: bool = False,
    fmm_nterms: int | None = None,
    nlmax: int | None = None,
    results_dir: Path = RESULTS_DIR,
) -> Path:
    """
    Build the expected result-file path from the run parameters.
    """
    stem = f'single_grain_sf{size_factor:.2g}_n{n}'

    if periodic:
        stem += '_P'

    if use_fmm:
        if fmm_nterms is None or nlmax is None:
            raise ValueError('fmm_nterms and nlmax must be given when use_fmm=True')
        stem += f'_fmm_on_N{fmm_nterms}_L{nlmax}'

    if adaptive:
        if adaptive_dh_min_t is None:
            raise ValueError('adaptive_dh_min_t must be given for adaptive runs')
        stem += f'_A_FS{abs(adaptive_dh_min_t):.1e}'

    return Path(results_dir) / f'{stem}.npz'


def find_result_file(
    size_factor: float,
    n: int,
    adaptive: bool = False,
    periodic: bool = False,
    adaptive_dh_min_t: float | None = None,
    use_fmm: bool = False,
    fmm_nterms: int | None = None,
    nlmax: int | None = None,
    results_dir: Path = RESULTS_DIR,
    raise_on_missing: bool = True,
) -> Path | None:
    """
    Find a result file from run parameters.

    By default this keeps the old behaviour and raises if the file is missing.
    Set raise_on_missing=False when using it from widgets, so the notebook can
    show a warning instead of crashing.
    """
    path = build_result_file_path(
        size_factor=size_factor,
        n=n,
        adaptive=adaptive,
        periodic=periodic,
        adaptive_dh_min_t=adaptive_dh_min_t,
        use_fmm=use_fmm,
        fmm_nterms=fmm_nterms,
        nlmax=nlmax,
        results_dir=results_dir,
    )

    if not path.exists():
        if raise_on_missing:
            raise FileNotFoundError(f'{path} not found')
        return None

    return path


def list_result_files(
    results_dir: Path = RESULTS_DIR,
    result_pattern: str = RESULT_PATTERN,
) -> list[Path]:
    """
    Return all result files matching RESULT_PATTERN.
    """
    return sorted(Path(results_dir).glob(result_pattern))


def _relative_or_full_path_label(path: Path, base_dir: Path = RESULTS_DIR) -> str:
    """
    Compact display label for a result file.
    """
    path = Path(path)

    try:
        return str(path.relative_to(base_dir))
    except ValueError:
        return str(path)


def _result_file_dropdown_options(
    results_dir: Path = RESULTS_DIR,
    result_pattern: str = RESULT_PATTERN,
) -> list[tuple[str, str]]:
    """
    Dropdown options as (label, value). The value is a string path because
    ipywidgets handles that more predictably than Path objects.
    """
    files = list_result_files(results_dir, result_pattern)

    if not files:
        return [('No result files found', '')]

    return [
        (_relative_or_full_path_label(path, results_dir), str(path))
        for path in files
    ]


def _card_title(title: str, subtitle: str | None = None) -> widgets.HTML:
    subtitle_html = f'<div style="color:#666;font-size:12px;margin-top:2px;">{subtitle}</div>' if subtitle else ''
    return widgets.HTML(
        f'<div style="margin-bottom:8px;">'
        f'<div style="font-size:15px;font-weight:700;letter-spacing:0.2px;">{title}</div>'
        f'{subtitle_html}'
        f'</div>'
    )


def _card(children: list[widgets.Widget], width: str = '100%') -> widgets.VBox:
    return widgets.VBox(
        children,
        layout=widgets.Layout(
            width=width,
            padding='12px 14px',
            margin='0 0 12px 0',
            border='1px solid #d8d8d8',
            border_radius='10px',
            overflow='hidden',
            box_sizing='border-box',
        ),
    )


def _html_warning(message: str) -> str:
    return (
        '<div style="background:#fff3f3;border-left:4px solid #b00020;'
        'padding:8px 10px;border-radius:6px;color:#5f000f;font-size:13px;">'
        '<b>Warning:</b> '
        f'{message}'
        '</div>'
    )


def _html_ok(message: str) -> str:
    return (
        '<div style="background:#f1f8f2;border-left:4px solid #1b5e20;'
        'padding:8px 10px;border-radius:6px;color:#174a1c;font-size:13px;">'
        '<b>Loaded:</b> '
        f'{message}'
        '</div>'
    )


def _safe_float_label(value: object, fmt: str = '.4g') -> str:
    try:
        value = float(value)
    except Exception:
        return 'n/a'
    if not np.isfinite(value):
        return 'n/a'
    return format(value, fmt)


def _infer_periodic_from_path(path: Path) -> bool:
    stem = Path(path).stem
    return '_P' in stem.split('_A_')[0]


def _format_run_title(path: Path | None, data: dict | None) -> str:
    if path is None or data is None:
        return 'Magnetisation comparison'

    metadata = data.get('metadata', {})
    size_factor = _safe_float_label(metadata.get('size_factor'))
    n = metadata.get('n', 'n/a')
    adaptive = bool(metadata.get('adaptive', False))
    periodic = _infer_periodic_from_path(path)
    dh_min = _safe_float_label(metadata.get('adaptive_dh_min_t'), '.2e')

    flags = []
    flags.append('periodic' if periodic else 'non-periodic')
    if adaptive:
        flags.append(f'adaptive, dh_min={dh_min} T')
    else:
        flags.append('fixed step')

    return f'{Path(path).name}  |  size_factor={size_factor}  |  n={n}  |  ' + '  |  '.join(flags)


def _format_run_title_html(path: Path | None, data: dict | None) -> str:
    title = _format_run_title(path, data)
    return (
        '<div style="font-size:15px;font-weight:700;line-height:1.35;'
        'padding:6px 8px 8px 8px;white-space:normal;max-width:1260px;">'
        f'{title}'
        '</div>'
    )


def interactive_magnetization_file_selector_and_orientation_3d(
    *,
    results_dir: Path = RESULTS_DIR,
    result_pattern: str = RESULT_PATTERN,
    size_factor_init: float = 30,
    n_init: int = 12,
    adaptive_init: bool = True,
    periodic_init: bool = True,
    adaptive_dh_min_t_init: float = 0.001,
    target_H1_init: float = -0.905,
    target_H2_init: float = -1.0,
    use_unit_vectors: bool = True,
    stride_init: int = 1,
    vector_length_init: float = 0.75,
    diff_vector_length_init: float = 1.0,
    elev_init: float = 60,
    azim_init: float = 25,
    last_panel_init: str = 'difference',
    positions_idx: int | None = None,
    continuous_update: bool = False,
):
    """
    Interactive file selector + magnetisation plot.

    Layout:
      - data selection spans the full dashboard width at the top;
      - persistent plot below on the left;
      - plot controls only in the right sidebar, with the same height as the plot.

    The data can be loaded in two ways:
      1. choose an existing .npz file from the dropdown;
      2. enter run parameters and find the matching file using find_result_file.

    Missing files are reported in the status widget rather than raising.
    """
    results_dir = Path(results_dir)

    narrow = widgets.Layout(width='100%')
    slider_layout = widgets.Layout(width='100%')

    file_dropdown = widgets.Dropdown(
        options=_result_file_dropdown_options(results_dir, result_pattern),
        description='',
        layout=narrow,
    )

    refresh_files_button = widgets.Button(
        description='Refresh',
        tooltip='Refresh the dropdown from RESULTS_DIR',
        button_style='',
        layout=widgets.Layout(width='105px'),
    )

    load_dropdown_button = widgets.Button(
        description='Load selected',
        button_style='primary',
        layout=widgets.Layout(width='150px'),
    )

    size_factor_input = widgets.FloatText(
        value=size_factor_init,
        description='size',
        layout=widgets.Layout(width='155px'),
        style={'description_width': '48px'},
    )

    n_input = widgets.IntText(
        value=n_init,
        description='n',
        layout=widgets.Layout(width='125px'),
        style={'description_width': '24px'},
    )

    adaptive_checkbox = widgets.Checkbox(
        value=adaptive_init,
        description='adaptive',
        indent=False,
        layout=widgets.Layout(width='120px'),
    )

    periodic_checkbox = widgets.Checkbox(
        value=periodic_init,
        description='periodic',
        indent=False,
        layout=widgets.Layout(width='120px'),
    )

    adaptive_dh_min_t_input = widgets.FloatText(
        value=adaptive_dh_min_t_init,
        description='dh_min',
        layout=widgets.Layout(width='170px'),
        style={'description_width': '55px'},
    )

    load_from_parameters_button = widgets.Button(
        description='Find/load',
        button_style='primary',
        layout=widgets.Layout(width='130px'),
    )

    status_html = widgets.HTML(layout=widgets.Layout(width='100%'))
    loaded_file_html = widgets.HTML(layout=widgets.Layout(width='100%'))
    plot_title_html = widgets.HTML(value=_format_run_title_html(None, None))

    idx1_slider = widgets.IntSlider(
        value=0,
        min=0,
        max=1,
        step=1,
        description='idx1',
        continuous_update=continuous_update,
        disabled=True,
        layout=slider_layout,
        style={'description_width': '65px'},
    )

    idx2_slider = widgets.IntSlider(
        value=1,
        min=0,
        max=1,
        step=1,
        description='idx2',
        continuous_update=continuous_update,
        disabled=True,
        layout=slider_layout,
        style={'description_width': '65px'},
    )

    elev_slider = widgets.FloatSlider(
        value=elev_init,
        min=-90,
        max=90,
        step=1,
        description='elev',
        continuous_update=continuous_update,
        readout_format='.0f',
        layout=slider_layout,
        style={'description_width': '65px'},
    )

    azim_slider = widgets.FloatSlider(
        value=azim_init,
        min=-180,
        max=180,
        step=1,
        description='azim',
        continuous_update=continuous_update,
        readout_format='.0f',
        layout=slider_layout,
        style={'description_width': '65px'},
    )

    stride_slider = widgets.IntSlider(
        value=stride_init,
        min=1,
        max=20,
        step=1,
        description='stride',
        continuous_update=continuous_update,
        layout=slider_layout,
        style={'description_width': '65px'},
    )

    vector_length_slider = widgets.FloatSlider(
        value=vector_length_init,
        min=0.05,
        max=3.0,
        step=0.05,
        description='vec len',
        continuous_update=continuous_update,
        readout_format='.2f',
        layout=slider_layout,
        style={'description_width': '65px'},
    )

    diff_vector_length_slider = widgets.FloatSlider(
        value=diff_vector_length_init,
        min=0.05,
        max=3.0,
        step=0.05,
        description='diff len',
        continuous_update=continuous_update,
        readout_format='.2f',
        layout=slider_layout,
        style={'description_width': '65px'},
    )

    last_panel_dropdown = widgets.Dropdown(
        options=[
            ('Difference vectors', 'difference'),
            ('Angular change', 'angle'),
        ],
        value=last_panel_init,
        description='last',
        layout=slider_layout,
        style={'description_width': '65px'},
    )

    index_info = widgets.HTML(layout=widgets.Layout(width='100%'))

    # Compact top-mounted data-selection area.  It is kept above the plot
    # rather than in the right sidebar so the sidebar contains only plot controls.
    file_picker_group = widgets.VBox(
        [
            widgets.HTML('<div style="font-size:12px;font-weight:700;margin-bottom:3px;">Result file</div>'),
            file_dropdown,
            widgets.HBox([load_dropdown_button, refresh_files_button], layout=widgets.Layout(gap='8px')),
        ],
        layout=widgets.Layout(
            flex='1 1 420px',
            min_width='340px',
            max_width='620px',
        ),
    )

    parameter_group = widgets.VBox(
        [
            widgets.HTML('<div style="font-size:12px;font-weight:700;margin-bottom:3px;">Find from run parameters</div>'),
            widgets.HBox([size_factor_input, n_input], layout=widgets.Layout(gap='8px', flex_flow='row wrap')),
            widgets.HBox([adaptive_checkbox, periodic_checkbox], layout=widgets.Layout(gap='12px', flex_flow='row wrap')),
            widgets.HBox([adaptive_dh_min_t_input, load_from_parameters_button], layout=widgets.Layout(gap='8px', flex_flow='row wrap')),
        ],
        layout=widgets.Layout(
            flex='0 1 430px',
            min_width='340px',
        ),
    )

    data_status_group = widgets.VBox(
        [status_html, loaded_file_html],
        layout=widgets.Layout(
            flex='1 1 320px',
            min_width='300px',
        ),
    )

    data_selection = _card(
        [
            _card_title('Data selection', 'Choose an existing result file, or reconstruct the expected filename from run parameters.'),
            widgets.HBox(
                [file_picker_group, parameter_group, data_status_group],
                layout=widgets.Layout(
                    width='100%',
                    gap='14px',
                    align_items='flex-start',
                    flex_flow='row wrap',
                    box_sizing='border-box',
                ),
            ),
        ]
    )

    plot_controls = _card(
        [
            _card_title('Plot controls', 'Controls are kept vertical so the figure remains visible while tuning.'),
            last_panel_dropdown,
            idx1_slider,
            idx2_slider,
            index_info,
            widgets.HTML('<hr style="border:none;border-top:1px solid #e5e5e5;margin:10px 0;">'),
            elev_slider,
            azim_slider,
            stride_slider,
            vector_length_slider,
            diff_vector_length_slider,
        ]
    )

    sidebar_width = '380px'
    plot_height = '820px'

    sidebar = widgets.VBox(
        [plot_controls],
        layout=widgets.Layout(
            width=sidebar_width,
            min_width='360px',
            max_width=sidebar_width,
            height=plot_height,
            min_height=plot_height,
            flex=f'0 0 {sidebar_width}',
            padding='0 0 0 10px',
            overflow='visible',
            box_sizing='border-box',
        ),
    )

    # The figure is intentionally wide.  The data-selection card spans the
    # full dashboard width above both the plot and the plot-control sidebar.
    # Below it, the plot and sidebar sit in one row with matching heights.
    fig = plt.figure(figsize=(15.8, 8.2), dpi=100, constrained_layout=False)
    try:
        fig.canvas.layout = widgets.Layout(
            width='100%',
            height=plot_height,
            min_width='760px',
            flex='1 1 auto',
            overflow='hidden',
        )
        fig.canvas.toolbar_position = 'bottom'
    except Exception:
        pass

    plot_box = widgets.VBox(
        [plot_title_html, fig.canvas],
        layout=widgets.Layout(
            width='100%',
            min_width='760px',
            max_width='none',
            height=plot_height,
            min_height=plot_height,
            flex='1 1 auto',
            overflow='hidden',
            box_sizing='border-box',
        ),
    )

    plot_and_controls = widgets.HBox(
        [plot_box, sidebar],
        layout=widgets.Layout(
            width='100%',
            max_width='100%',
            min_width='0',
            align_items='stretch',
            gap='14px',
            display='flex',
            flex_flow='row nowrap',
            overflow='hidden',
            box_sizing='border-box',
        ),
    )

    dashboard = widgets.VBox(
        [data_selection, plot_and_controls],
        layout=widgets.Layout(
            width='100%',
            max_width='100%',
            min_width='0',
            gap='8px',
            display='flex',
            flex_flow='column nowrap',
            overflow='hidden',
            box_sizing='border-box',
        ),
    )
    display(dashboard)

    state = {
        'fig': fig,
        'axs': None,
        'data': None,
        'path': None,
        'is_drawing': False,
    }

    def clear_plot():
        fig.clf()
        fig.canvas.draw_idle()
        index_info.value = ''

    def set_warning(message: str):
        status_html.value = _html_warning(message)

    def set_loaded_message(path: Path, data: dict):
        title = _format_run_title(path, data)
        loaded_file_html.value = (
            '<div style="font-size:12px;line-height:1.45;margin-top:8px;">'
            f'<b>Current:</b><br><code>{Path(path).name}</code><br>'
            f'<b>Run:</b><br>{title}<br>'
            f'<b>H steps:</b> {np.asarray(data["H"]).size}'
            '</div>'
        )
        status_html.value = _html_ok(Path(path).name)
        plot_title_html.value = _format_run_title_html(path, data)

    def refresh_file_options(preferred_path: Path | None = None):
        old_value = file_dropdown.value
        options = _result_file_dropdown_options(results_dir, result_pattern)
        file_dropdown.options = options

        valid_values = {value for _, value in options}

        if preferred_path is not None and str(preferred_path) in valid_values:
            file_dropdown.value = str(preferred_path)
        elif old_value in valid_values:
            file_dropdown.value = old_value
        else:
            file_dropdown.value = options[0][1]

    def update_index_info():
        data = state['data']

        if data is None:
            index_info.value = ''
            return

        H = np.asarray(data['H'], dtype=float)

        i1 = idx1_slider.value
        i2 = idx2_slider.value

        index_info.value = (
            '<div style="background:#f7f7f7;padding:7px 9px;border-radius:6px;font-size:12px;line-height:1.45;">'
            f'<b>idx1:</b> {i1}, H = {H[i1]:.6g}<br>'
            f'<b>idx2:</b> {i2}, H = {H[i2]:.6g}<br>'
            f'<b>ΔH:</b> {H[i2] - H[i1]:.6g}'
            '</div>'
        )

    def redraw(change=None):
        if state['is_drawing']:
            return

        data = state['data']

        if data is None:
            return

        state['is_drawing'] = True

        try:
            update_index_info()

            axs = draw_magnetization_curve_and_orientation_3d_on_figure(
                fig,
                data['H'],
                data['Mx'],
                data['My'],
                data['Mz'],
                data['res'],
                idx1_slider.value,
                idx2_slider.value,
                use_unit_vectors=use_unit_vectors,
                stride=stride_slider.value,
                vector_length=vector_length_slider.value,
                diff_vector_length=diff_vector_length_slider.value,
                elev=elev_slider.value,
                azim=azim_slider.value,
                last_panel=last_panel_dropdown.value,
                positions_idx=positions_idx,
                plot_title=None,
            )

            state['axs'] = axs

        except Exception as exc:
            set_warning(f'Could not draw plot: {exc}')
        finally:
            state['is_drawing'] = False

    def load_path(path_like):
        if path_like in (None, ''):
            state['data'] = None
            state['path'] = None
            loaded_file_html.value = ''
            plot_title_html.value = _format_run_title_html(None, None)
            set_warning(f'No result files were found in {results_dir.resolve()}')
            clear_plot()
            return

        path = Path(path_like)

        if not path.exists():
            state['data'] = None
            state['path'] = None
            loaded_file_html.value = ''
            plot_title_html.value = _format_run_title_html(None, None)
            set_warning(f'{path} was not found')
            clear_plot()
            return

        try:
            data = load_result_data(path)
        except Exception as exc:
            state['data'] = None
            state['path'] = None
            loaded_file_html.value = ''
            plot_title_html.value = _format_run_title_html(None, None)
            set_warning(f'Could not load {path}: {exc}')
            clear_plot()
            return

        H = np.asarray(data['H'], dtype=float).ravel()

        if H.size < 2:
            state['data'] = None
            state['path'] = None
            loaded_file_html.value = ''
            plot_title_html.value = _format_run_title_html(None, None)
            set_warning(f'{path} has only {H.size} field step(s); at least two are needed')
            clear_plot()
            return

        state['data'] = data
        state['path'] = path

        idx1 = find_nearest_field_index(H, target_H1_init)
        idx2 = find_nearest_field_index(H, target_H2_init)

        # Update slider ranges before setting values.
        idx1_slider.max = H.size - 1
        idx2_slider.max = H.size - 1
        idx1_slider.disabled = False
        idx2_slider.disabled = False
        idx1_slider.value = int(idx1)
        idx2_slider.value = int(idx2)

        refresh_file_options(preferred_path=path)
        set_loaded_message(path, data)
        redraw()

    def on_refresh_files_clicked(button):
        refresh_file_options(preferred_path=state['path'])

    def on_load_dropdown_clicked(button):
        load_path(file_dropdown.value)

    def on_load_from_parameters_clicked(button):
        try:
            expected_path = build_result_file_path(
                size_factor=size_factor_input.value,
                n=n_input.value,
                adaptive=adaptive_checkbox.value,
                periodic=periodic_checkbox.value,
                adaptive_dh_min_t=adaptive_dh_min_t_input.value,
                results_dir=results_dir,
            )
        except Exception as exc:
            set_warning(f'Could not build result-file name from inputs: {exc}')
            return

        path = find_result_file(
            size_factor=size_factor_input.value,
            n=n_input.value,
            adaptive=adaptive_checkbox.value,
            periodic=periodic_checkbox.value,
            adaptive_dh_min_t=adaptive_dh_min_t_input.value,
            results_dir=results_dir,
            raise_on_missing=False,
        )

        if path is None:
            set_warning(f'No matching file found: expected {expected_path}')
            return

        load_path(path)

    refresh_files_button.on_click(on_refresh_files_clicked)
    load_dropdown_button.on_click(on_load_dropdown_clicked)
    load_from_parameters_button.on_click(on_load_from_parameters_clicked)

    for widget in [
        idx1_slider,
        idx2_slider,
        elev_slider,
        azim_slider,
        stride_slider,
        vector_length_slider,
        diff_vector_length_slider,
        last_panel_dropdown,
    ]:
        widget.observe(redraw, names='value')

    # Try to start with the file corresponding to the default parameter values.
    default_path = find_result_file(
        size_factor=size_factor_init,
        n=n_init,
        adaptive=adaptive_init,
        periodic=periodic_init,
        adaptive_dh_min_t=adaptive_dh_min_t_init,
        results_dir=results_dir,
        raise_on_missing=False,
    )

    if default_path is not None:
        load_path(default_path)
    elif file_dropdown.value:
        load_path(file_dropdown.value)
    else:
        set_warning(f'No result files were found in {results_dir.resolve()}')
        clear_plot()

    return {
        'file_dropdown': file_dropdown,
        'refresh_files': refresh_files_button,
        'load_selected_file': load_dropdown_button,
        'size_factor': size_factor_input,
        'n': n_input,
        'adaptive': adaptive_checkbox,
        'periodic': periodic_checkbox,
        'adaptive_dh_min_t': adaptive_dh_min_t_input,
        'load_from_parameters': load_from_parameters_button,
        'idx1': idx1_slider,
        'idx2': idx2_slider,
        'elev': elev_slider,
        'azim': azim_slider,
        'stride': stride_slider,
        'vector_length': vector_length_slider,
        'diff_vector_length': diff_vector_length_slider,
        'last_panel': last_panel_dropdown,
        'index_info': index_info,
        'status': status_html,
        'loaded_file': loaded_file_html,
        'plot_title': plot_title_html,
        'dashboard': dashboard,
        'state': state,
    }


__all__ = [
    "RESULTS_DIR",
    "RESULT_PATTERN",
    "add_view_sliders",
    "angle_between_vectors_deg",
    "build_result_file_path",
    "draw_magnetization_curve_and_orientation_3d_on_figure",
    "extract_tile_magnetization",
    "extract_tile_positions",
    "find_nearest_field_index",
    "find_result_file",
    "interactive_magnetization_curve_and_orientation_3d",
    "interactive_magnetization_file_selector_and_orientation_3d",
    "list_result_files",
    "load_result_data",
    "load_scalar",
    "normalise_vectors",
    "plot_magnetization_curve_and_orientation_3d",
    "set_3d_view",
]
