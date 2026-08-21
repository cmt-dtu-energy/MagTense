"""Derived permanent-magnet metrics for descending hysteresis curves."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

MU0 = 4.0 * np.pi * 1e-7
KJ_PER_M3_PER_MGOE = 7.957747154594767


@dataclass(frozen=True)
class HysteresisMetrics:
    """Scalar metrics and diagnostics derived from one hysteresis curve."""

    Hc_A_per_m: float
    Hc_T: float
    Mr_A_per_m: float
    mu0_Mr_T: float
    Mr_over_Ms: float
    BH_max_J_per_m3: float
    BH_max_kJ_per_m3: float
    BH_max_MGOe: float
    coercivity_status: str
    remanence_status: str
    energy_product_status: str

    def as_dict(self) -> dict[str, float | str]:
        return asdict(self)


def _clean_curve(H: np.ndarray, M: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    H = np.asarray(H, dtype=float).ravel()
    M = np.asarray(M, dtype=float).ravel()
    if H.size != M.size:
        raise ValueError(f"H and M must have the same length, got {H.size} and {M.size}")
    finite = np.isfinite(H) & np.isfinite(M)
    return H[finite], M[finite]


def interpolated_crossing(
    x: np.ndarray,
    y: np.ndarray,
    *,
    target: float = 0.0,
) -> tuple[float, str]:
    """Return the first linearly interpolated x where y crosses ``target``."""
    x, y = _clean_curve(x, y)
    if x.size < 2:
        return float("nan"), "insufficient_points"

    shifted = y - target
    exact = np.flatnonzero(shifted == 0.0)
    if exact.size:
        return float(x[int(exact[0])]), "ok"

    crossings = np.flatnonzero(shifted[:-1] * shifted[1:] < 0.0)
    if crossings.size == 0:
        return float("nan"), "no_crossing"

    index = int(crossings[0])
    x0, x1 = x[index], x[index + 1]
    y0, y1 = shifted[index], shifted[index + 1]
    return float(x0 - y0 * (x1 - x0) / (y1 - y0)), "ok"


def interpolated_coercivity(H: np.ndarray, M: np.ndarray) -> float:
    """Backward-compatible convenience wrapper returning H at M=0."""
    value, _ = interpolated_crossing(H, M)
    return value


def interpolated_remanence(H: np.ndarray, M: np.ndarray) -> tuple[float, str]:
    """Return linearly interpolated M at H=0."""
    value, status = interpolated_crossing(M, H)
    return value, status


def maximum_energy_product(
    H: np.ndarray,
    M: np.ndarray,
) -> tuple[float, str]:
    """Return the standard second-quadrant maximum of ``-B*H`` in J/m^3.

    H and M are in A/m and B is calculated as ``mu0 * (H + M)``. Each
    adjacent pair is treated as a linear B(H) segment. Endpoints, quadrant
    boundaries, and the quadratic segment maximum are considered.
    """
    H, M = _clean_curve(H, M)
    if H.size < 2:
        return float("nan"), "insufficient_points"

    B = MU0 * (H + M)
    best = float("nan")

    for h0, h1, b0, b1 in zip(H[:-1], H[1:], B[:-1], B[1:]):
        dh = h1 - h0
        db = b1 - b0
        lower, upper = 0.0, 1.0

        # Clip the linear segment to H <= 0 and B >= 0.
        for value0, delta, keep_nonnegative in (
            (-h0, -dh, True),
            (b0, db, True),
        ):
            if delta == 0.0:
                if (value0 >= 0.0) != keep_nonnegative:
                    lower, upper = 1.0, 0.0
                    break
                continue
            crossing = -value0 / delta
            if delta > 0.0:
                lower = max(lower, crossing)
            else:
                upper = min(upper, crossing)

        lower = max(0.0, lower)
        upper = min(1.0, upper)
        if lower > upper:
            continue

        candidates = [lower, upper]
        quadratic = -dh * db
        linear = -(h0 * db + b0 * dh)
        if quadratic != 0.0:
            vertex = -linear / (2.0 * quadratic)
            if lower <= vertex <= upper:
                candidates.append(vertex)

        for fraction in candidates:
            h = h0 + fraction * dh
            b = b0 + fraction * db
            product = -h * b
            if h <= 1e-12 and b >= -1e-12 and product >= 0.0:
                best = product if not np.isfinite(best) else max(best, product)

    if not np.isfinite(best):
        return float("nan"), "no_second_quadrant"
    return float(best), "ok"


def calculate_hysteresis_metrics(
    H_A_per_m: np.ndarray,
    M_A_per_m: np.ndarray,
    Ms_A_per_m: float,
) -> HysteresisMetrics:
    """Calculate coercivity, remanence, and maximum energy product."""
    Hc, coercivity_status = interpolated_crossing(H_A_per_m, M_A_per_m)
    Mr, remanence_status = interpolated_remanence(H_A_per_m, M_A_per_m)
    BH_max, energy_product_status = maximum_energy_product(H_A_per_m, M_A_per_m)

    Ms = float(Ms_A_per_m)
    ratio = Mr / Ms if np.isfinite(Mr) and np.isfinite(Ms) and Ms != 0.0 else np.nan
    return HysteresisMetrics(
        Hc_A_per_m=Hc,
        Hc_T=MU0 * Hc if np.isfinite(Hc) else np.nan,
        Mr_A_per_m=Mr,
        mu0_Mr_T=MU0 * Mr if np.isfinite(Mr) else np.nan,
        Mr_over_Ms=float(ratio),
        BH_max_J_per_m3=BH_max,
        BH_max_kJ_per_m3=BH_max / 1000.0 if np.isfinite(BH_max) else np.nan,
        BH_max_MGOe=(BH_max / 1000.0 / KJ_PER_M3_PER_MGOE if np.isfinite(BH_max) else np.nan),
        coercivity_status=coercivity_status,
        remanence_status=remanence_status,
        energy_product_status=energy_product_status,
    )
