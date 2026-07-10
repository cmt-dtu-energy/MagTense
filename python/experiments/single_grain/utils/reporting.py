"""Dataset summaries and deterministic selections used by analysis notebooks."""

from __future__ import annotations

from collections import Counter
from typing import Iterable

import numpy as np

from .results import ResultRecord


def material_summary(records: Iterable[ResultRecord]) -> list[dict[str, object]]:
    """Return one property/count row per distinct material and shape variant."""
    grouped: dict[str, list[ResultRecord]] = {}
    for record in records:
        grouped.setdefault(record.comparison_label, []).append(record)
    rows = []
    for comparison_label, group in sorted(grouped.items()):
        first = group[0]
        rows.append(
            {
                "material": first.material,
                "shape": first.shape,
                "shape_variant": first.shape_variant,
                "comparison_label": comparison_label,
                "mu0_Ms_T": first.mu0_Ms_T,
                "A0_J_per_m": first.A0_J_per_m,
                "K0_MJ_per_m3": first.K0_J_per_m3 / 1e6,
                "result_count": len(group),
                "sizes_nm": sorted({round(record.size_nm, 9) for record in group}),
                "resolutions": sorted({record.n for record in group}),
                "dh_min_T": sorted(
                    {
                        record.adaptive_dh_min_t
                        for record in group
                        if np.isfinite(record.adaptive_dh_min_t)
                    }
                ),
                "periodic_modes": sorted({record.periodic for record in group}),
            }
        )
    return rows


def missing_combinations(
    records: Iterable[ResultRecord],
    *,
    materials: Iterable[str] | None = None,
    shapes: Iterable[str] | None = None,
    sizes_nm: Iterable[float] | None = None,
    resolutions: Iterable[int] | None = None,
    dh_min_values: Iterable[float] | None = None,
    periodic_modes: Iterable[bool] = (False, True),
) -> list[tuple[str, str, float, int, float, bool]]:
    """List expected numerical combinations that are not present."""
    records = list(records)
    materials = list(materials or sorted({record.material for record in records}))
    shapes = list(shapes or sorted({record.shape_variant for record in records}))
    sizes_nm = list(sizes_nm or sorted({record.size_nm for record in records}))
    resolutions = list(resolutions or sorted({record.n for record in records}))
    dh_min_values = list(
        dh_min_values
        or sorted(
            {
                record.adaptive_dh_min_t
                for record in records
                if np.isfinite(record.adaptive_dh_min_t)
            }
        )
    )
    present = {
        (
            record.material,
            record.shape_variant,
            round(record.size_nm, 9),
            record.n,
            round(record.adaptive_dh_min_t, 12),
            record.periodic,
        )
        for record in records
        if np.isfinite(record.adaptive_dh_min_t)
    }
    missing = []
    for material in materials:
        for shape in shapes:
            for size_nm in sizes_nm:
                for n in resolutions:
                    for dh_min in dh_min_values:
                        for periodic in periodic_modes:
                            key = (
                                material,
                                shape,
                                round(float(size_nm), 9),
                                n,
                                round(float(dh_min), 12),
                                periodic,
                            )
                            if key not in present:
                                missing.append(
                                    (material, shape, float(size_nm), n, float(dh_min), periodic)
                                )
    return missing


def select_representative_records(
    records: Iterable[ResultRecord],
    *,
    target_size_nm: float = 40.0,
) -> list[ResultRecord]:
    """Choose the nearest size for every material, shape, and periodic mode."""
    grouped: dict[tuple[str, str, bool], list[ResultRecord]] = {}
    for record in records:
        grouped.setdefault(
            (record.material, record.shape_variant, record.periodic), []
        ).append(record)
    return [
        min(group, key=lambda record: (abs(record.size_nm - target_size_nm), str(record.path)))
        for _, group in sorted(grouped.items())
    ]


def diagnostic_counts(records: Iterable[ResultRecord]) -> dict[str, int]:
    """Count metric diagnostics so incomplete curves are visible."""
    counter: Counter[str] = Counter()
    for record in records:
        counter[f"coercivity:{record.metrics.coercivity_status}"] += 1
        counter[f"remanence:{record.metrics.remanence_status}"] += 1
        counter[f"energy_product:{record.metrics.energy_product_status}"] += 1
    return dict(sorted(counter.items()))


__all__ = [
    "diagnostic_counts",
    "material_summary",
    "missing_combinations",
    "select_representative_records",
]
