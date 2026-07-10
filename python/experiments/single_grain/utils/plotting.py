"""Static plotting helpers for single-grain material comparisons."""

from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path
from typing import Callable, Iterable

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from .metrics import MU0
from .results import ResultRecord


def plot_hysteresis_curves(
    records: Iterable[ResultRecord],
    *,
    ax=None,
    title: str = "Representative hysteresis curves",
    label_func: Callable[[ResultRecord], str] | None = None,
):
    """Plot mu0*M against mu0*H for the supplied records."""
    if ax is None:
        _, ax = plt.subplots(figsize=(7.2, 4.5))
    for record in records:
        mode = "P" if record.periodic else "non-P"
        label = (
            label_func(record)
            if label_func is not None
            else f"{record.comparison_label}, {record.size_nm:g} nm, {mode}"
        )
        ax.plot(
            record.H_T,
            MU0 * record.Mz_A_per_m,
            ".-",
            markersize=3,
            linewidth=1.2,
            label=label,
        )
    ax.axhline(0.0, color="0.5", linewidth=0.7)
    ax.axvline(0.0, color="0.5", linewidth=0.7)
    ax.set(xlabel=r"$\mu_0 H$ [T]", ylabel=r"$\mu_0 M_z$ [T]", title=title)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(fontsize=8)
    return ax


def plot_second_quadrant(
    records: Iterable[ResultRecord],
    *,
    ax=None,
    title: str = "Second-quadrant demagnetization curves",
):
    """Plot B(H) and mark each calculated maximum energy product."""
    if ax is None:
        _, ax = plt.subplots(figsize=(7.2, 4.5))
    for record in records:
        B = MU0 * (record.H_A_per_m + record.Mz_A_per_m)
        mask = (record.H_A_per_m <= 0.0) & (B >= 0.0)
        mode = "P" if record.periodic else "non-P"
        ax.plot(
            record.H_A_per_m[mask] / 1e3,
            B[mask],
            ".-",
            markersize=3,
            label=(
                f"{record.comparison_label}, {record.size_nm:g} nm, {mode}; "
                f"{record.metrics.BH_max_kJ_per_m3:.1f} kJ/m3"
            ),
        )
    ax.set(xlabel="H [kA/m]", ylabel="B [T]", title=title)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(fontsize=8)
    return ax


_METRICS = {
    "abs_Hc_T": (lambda record: abs(record.metrics.Hc_T), r"$|\mu_0 H_c|$ [T]"),
    "mu0_Mr_T": (lambda record: record.metrics.mu0_Mr_T, r"$\mu_0 M_r$ [T]"),
    "Mr_over_Ms": (lambda record: record.metrics.Mr_over_Ms, r"$M_r/M_s$ [-]"),
    "BH_max_kJ_per_m3": (
        lambda record: record.metrics.BH_max_kJ_per_m3,
        r"$(BH)_{\max}$ [kJ/m$^3$]",
    ),
}


def plot_metric_vs_size(
    records: Iterable[ResultRecord],
    metric: str,
    *,
    ax=None,
    title: str | None = None,
    label_func: Callable[[ResultRecord, bool, int, float], str] | None = None,
):
    """Plot one metric by physical size, preserving each numerical series."""
    if metric not in _METRICS:
        raise ValueError(f"Unknown metric {metric!r}; choose from {sorted(_METRICS)}")
    value, ylabel = _METRICS[metric]
    if ax is None:
        _, ax = plt.subplots(figsize=(6.5, 4.2))

    groups: dict[tuple[str, bool, int, float], list[ResultRecord]] = defaultdict(list)
    for record in records:
        groups[
            (
                record.comparison_label,
                record.periodic,
                record.n,
                record.adaptive_dh_min_t,
            )
        ].append(record)

    for (comparison_label, periodic, n, dh_min), group in sorted(groups.items(), key=str):
        group.sort(key=lambda record: record.size_nm)
        label = (
            label_func(group[0], periodic, n, dh_min)
            if label_func is not None
            else (
                f"{comparison_label}, {'P' if periodic else 'non-P'}, "
                f"n={n}, dh={dh_min:g}"
            )
        )
        ax.plot(
            [record.size_nm for record in group],
            [value(record) for record in group],
            "o-" if periodic else "o--",
            linewidth=1.3,
            markersize=4,
            label=label,
        )
    ax.set(xlabel="grain size [nm]", ylabel=ylabel, title=title or ylabel)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(fontsize=7)
    return ax


def plot_material_metrics(
    records: Iterable[ResultRecord],
    *,
    title: str = "Single-grain permanent-magnet metrics",
    label_func: Callable[[ResultRecord, bool, int, float], str] | None = None,
):
    """Create the standard four-panel material comparison figure."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    for ax, metric in zip(axes.ravel(), _METRICS):
        plot_metric_vs_size(records, metric, ax=ax, label_func=label_func)
    fig.suptitle(title)
    fig.tight_layout()
    return fig, axes


def plot_numerical_convergence(
    records: Iterable[ResultRecord],
    metric: str = "abs_Hc_T",
):
    """Compare raw resolution and adaptive-step results without averaging."""
    records = list(records)
    value, ylabel = _METRICS[metric]
    comparison_labels = sorted({record.comparison_label for record in records})
    n_values = sorted({record.n for record in records})
    all_dh_values = sorted(
        {
            record.adaptive_dh_min_t
            for record in records
            if np.isfinite(record.adaptive_dh_min_t)
        }
    )
    color_cycle = ["black", "red", "blue", "green", "purple", "orange"]
    marker_cycle = ["o", "s", "^", "D", "v", "P"]
    linestyle_cycle = ["-", "--", "-.", ":", (0, (3, 1, 1, 1)), (0, (5, 2))]
    n_colors = {
        n: color_cycle[index % len(color_cycle)]
        for index, n in enumerate(n_values)
    }
    dh_styles = {
        dh_min: (
            marker_cycle[index % len(marker_cycle)],
            linestyle_cycle[index % len(linestyle_cycle)],
        )
        for index, dh_min in enumerate(all_dh_values)
    }
    ncols = 2
    nrows = max(1, math.ceil(len(comparison_labels) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(15.5, 3.8 * nrows),
        squeeze=False,
    )
    for ax, comparison_label in zip(axes.ravel(), comparison_labels):
        subset = [record for record in records if record.comparison_label == comparison_label]
        for periodic in (False, True):
            for n in sorted({record.n for record in subset}):
                subset_dh_values = sorted(
                    {
                        record.adaptive_dh_min_t
                        for record in subset
                        if record.periodic == periodic
                        and record.n == n
                        and np.isfinite(record.adaptive_dh_min_t)
                    }
                )
                for dh_min in subset_dh_values:
                    group = sorted(
                        [
                            record
                            for record in subset
                            if record.periodic == periodic
                            and record.n == n
                            and np.isclose(record.adaptive_dh_min_t, dh_min)
                        ],
                        key=lambda record: record.size_nm,
                    )
                    if not group:
                        continue
                    marker, linestyle = dh_styles[dh_min]
                    ax.plot(
                        [record.size_nm for record in group],
                        [value(record) for record in group],
                        color=n_colors[n],
                        marker=marker,
                        linestyle=linestyle,
                        label=f"{'P' if periodic else 'non-P'}, n={n}, dh={dh_min:g}",
                        markersize=4,
                        linewidth=1.2,
                        alpha=0.85 if periodic else 0.7,
        )
        ax.set(title=comparison_label, xlabel="grain size [nm]", ylabel=ylabel)
        ax.grid(True, linestyle=":", alpha=0.6)
    for ax in axes.ravel()[len(comparison_labels):]:
        ax.set_visible(False)
    spatial_handles = [
        Line2D([0], [0], color=n_colors[n], linewidth=2.0, label=f"n={n}")
        for n in n_values
    ]
    field_handles = [
        Line2D(
            [0],
            [0],
            color="0.35",
            marker=dh_styles[dh_min][0],
            linestyle=dh_styles[dh_min][1],
            linewidth=1.5,
            markersize=5,
            label=f"dh={dh_min:g} T",
        )
        for dh_min in all_dh_values
    ]
    fig.tight_layout(rect=(0.0, 0.0, 0.80, 1.0))
    spatial_legend = fig.legend(
        handles=spatial_handles,
        title="Spatial resolution\n(color)",
        loc="upper left",
        bbox_to_anchor=(0.82, 0.98),
        frameon=False,
        borderaxespad=0.0,
    )
    field_legend = fig.legend(
        handles=field_handles,
        title="Field resolution\n(marker + linestyle)",
        loc="upper left",
        bbox_to_anchor=(0.82, 0.66),
        frameon=False,
        borderaxespad=0.0,
    )
    fig.add_artist(spatial_legend)
    fig.add_artist(field_legend)
    return fig, axes


def plot_property_space(records: Iterable[ResultRecord], *, ax=None):
    """Plot material properties, coloured by maximum energy product."""
    records = list(records)
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5))
    values = [record.metrics.BH_max_kJ_per_m3 for record in records]
    finite_values = [value for value in values if np.isfinite(value)]
    norm = plt.Normalize(
        vmin=min(finite_values) if finite_values else 0.0,
        vmax=max(finite_values) if finite_values else 1.0,
    )
    for periodic, marker in ((False, "o"), (True, "s")):
        points = [record for record in records if record.periodic == periodic]
        if not points:
            continue
        ax.scatter(
            [record.mu0_Ms_T for record in points],
            [record.K0_J_per_m3 / 1e6 for record in points],
            c=[record.metrics.BH_max_kJ_per_m3 for record in points],
            s=[40 + 15 * np.log10(max(record.A0_J_per_m / 1e-12, 1.0)) for record in points],
            cmap="viridis",
            norm=norm,
            marker=marker,
            edgecolor="black",
            linewidth=0.5,
            label="periodic" if periodic else "non-periodic",
        )
    labelled = set()
    for record in records:
        if record.comparison_label in labelled:
            continue
        labelled.add(record.comparison_label)
        ax.annotate(record.comparison_label, (record.mu0_Ms_T, record.K0_J_per_m3 / 1e6), xytext=(4, 4), textcoords="offset points")
    ax.set(xlabel=r"$\mu_0 M_s$ [T]", ylabel=r"$K_0$ [MJ/m$^3$]", title="Material-property space")
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(fontsize=8)
    scalar_map = plt.cm.ScalarMappable(norm=norm, cmap="viridis")
    ax.figure.colorbar(scalar_map, ax=ax, label=r"$(BH)_{\max}$ [kJ/m$^3$]")
    return ax


def save_figure(fig, path: Path | str, *, dpi: int = 200) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path


__all__ = [
    "plot_hysteresis_curves",
    "plot_material_metrics",
    "plot_metric_vs_size",
    "plot_numerical_convergence",
    "plot_property_space",
    "plot_second_quadrant",
    "save_figure",
]
