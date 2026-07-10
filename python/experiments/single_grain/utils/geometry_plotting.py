"""Three-dimensional visualization helpers for adaptive prism meshes."""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from .geometry import PrismMesh


def _cube_faces(center: np.ndarray, dimensions: np.ndarray) -> list[np.ndarray]:
    half = dimensions / 2.0
    vertices = center + np.array(
        [
            [-half[0], -half[1], -half[2]],
            [half[0], -half[1], -half[2]],
            [half[0], half[1], -half[2]],
            [-half[0], half[1], -half[2]],
            [-half[0], -half[1], half[2]],
            [half[0], -half[1], half[2]],
            [half[0], half[1], half[2]],
            [-half[0], half[1], half[2]],
        ]
    )
    return [
        vertices[[0, 1, 2, 3]],
        vertices[[4, 5, 6, 7]],
        vertices[[0, 1, 5, 4]],
        vertices[[2, 3, 7, 6]],
        vertices[[1, 2, 6, 5]],
        vertices[[0, 3, 7, 4]],
    ]


def plot_prism_mesh(
    mesh: PrismMesh,
    *,
    fig: Figure | None = None,
    ax: Axes | None = None,
    color_by_level: bool = True,
    cmap: str = "viridis",
    alpha: float = 0.35,
    wireframe: bool = True,
    hide_axes: bool = False,
    max_plot_tiles: int | None = 5000,
    elevation: float = 25.0,
    azimuth: float = 35.0,
) -> tuple[Figure, Axes]:
    """Plot prism cells in 3D and return the Matplotlib figure and axes."""
    if max_plot_tiles is not None and max_plot_tiles <= 0:
        raise ValueError("max_plot_tiles must be positive or None")
    if fig is None and ax is None:
        fig = plt.figure(figsize=(7, 6))
        ax = fig.add_subplot(111, projection="3d")
    elif fig is None or ax is None:
        raise ValueError("fig and ax must either both be supplied or both be omitted")
    else:
        ax.clear()

    assert fig is not None and ax is not None
    count = mesh.achieved_tiles
    if max_plot_tiles is not None and count > max_plot_tiles:
        indices = np.linspace(0, count - 1, max_plot_tiles, dtype=int)
    else:
        indices = np.arange(count)

    levels = mesh.levels[indices]
    if color_by_level and len(levels):
        minimum = int(np.min(levels))
        span = max(1, int(np.max(levels)) - minimum)
        colors = colormaps[cmap]((levels - minimum) / span)
    else:
        colors = np.tile(np.array([0.2, 0.55, 0.85, 1.0]), (len(indices), 1))

    faces: list[np.ndarray] = []
    face_colors: list[np.ndarray] = []
    for index, color in zip(indices, colors):
        cell_faces = _cube_faces(mesh.centers[index], mesh.dimensions[index])
        faces.extend(cell_faces)
        face_colors.extend([color] * len(cell_faces))

    collection = Poly3DCollection(
        faces,
        facecolors=face_colors,
        edgecolors="black" if wireframe else "none",
        linewidths=0.25 if wireframe else 0.0,
        alpha=alpha,
    )
    ax.add_collection3d(collection)
    lower, upper = mesh.root_bounds
    padding = np.maximum((upper - lower) * 0.03, np.finfo(float).eps)
    ax.set_xlim(lower[0] - padding[0], upper[0] + padding[0])
    ax.set_ylim(lower[1] - padding[1], upper[1] + padding[1])
    ax.set_zlim(lower[2] - padding[2], upper[2] + padding[2])
    ax.set_box_aspect(upper - lower)
    ax.view_init(elev=elevation, azim=azimuth)
    ax.set_title(f"{mesh.shape.replace('_', ' ').title()} ({mesh.achieved_tiles} cells)")
    if hide_axes:
        ax.set_axis_off()
    else:
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
    return fig, ax


def interactive_prism_mesh(
    mesh: PrismMesh,
    *,
    elevation: float = 25.0,
    azimuth: float = 35.0,
    **plot_kwargs: Any,
) -> tuple[Figure, Axes, Any]:
    """Plot a mesh with notebook sliders for elevation and azimuth.

    The returned controls can be displayed directly in a notebook. With an
    interactive Matplotlib backend, the usual mouse-based rotation also remains
    available.
    """
    try:
        import ipywidgets as widgets
    except ImportError as exc:  # pragma: no cover - depends on notebook extras
        raise ImportError("interactive_prism_mesh requires ipywidgets") from exc

    fig, ax = plot_prism_mesh(
        mesh,
        elevation=elevation,
        azimuth=azimuth,
        **plot_kwargs,
    )
    elevation_slider = widgets.FloatSlider(
        value=elevation,
        min=-90.0,
        max=90.0,
        step=1.0,
        description="Elevation",
        continuous_update=True,
    )
    azimuth_slider = widgets.FloatSlider(
        value=azimuth,
        min=-180.0,
        max=180.0,
        step=1.0,
        description="Azimuth",
        continuous_update=True,
    )

    def update_view(change: dict[str, Any]) -> None:
        del change
        ax.view_init(elev=elevation_slider.value, azim=azimuth_slider.value)
        fig.canvas.draw_idle()

    elevation_slider.observe(update_view, names="value")
    azimuth_slider.observe(update_view, names="value")
    controls = widgets.HBox([elevation_slider, azimuth_slider])
    return fig, ax, controls


__all__ = ["interactive_prism_mesh", "plot_prism_mesh"]
