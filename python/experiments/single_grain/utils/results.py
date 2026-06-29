"""Discovery and normalized loading of single-grain result files."""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .metrics import MU0, HysteresisMetrics, calculate_hysteresis_metrics

_ADAPTIVE_PATTERN = re.compile(
    r"_A_FS(?P<dh>[+-]?\d+(?:\.\d+)?(?:e[+-]?\d+)?)", re.IGNORECASE
)
_PERIODIC_PATTERN = re.compile(r"_n\d+_P(?:_|$)")


@dataclass
class ResultRecord:
    """Normalized metadata, arrays, and metrics for one result file."""

    path: Path
    material: str
    preset: str | None
    batch: str | None
    shape: str
    shape_variant: str
    shape_metadata: dict[str, Any]
    mu0_Ms_T: float
    Ms_A_per_m: float
    A0_J_per_m: float
    K0_J_per_m3: float
    L_m: float
    size_nm: float
    size_factor: float
    n: int
    ntot: int
    adaptive: bool
    periodic: bool
    adaptive_dh_min_t: float
    use_fmm: bool
    cuda: bool
    cvode: bool
    runtime_s: float
    H_T: np.ndarray
    H_A_per_m: np.ndarray
    M_A_per_m: np.ndarray
    Mx_A_per_m: np.ndarray
    My_A_per_m: np.ndarray
    Mz_A_per_m: np.ndarray
    res: object | None
    metrics: HysteresisMetrics

    @property
    def comparison_label(self) -> str:
        """Return the material/shape label used for comparison plots."""
        return f"{self.material} / {self.shape_variant}"

    def summary_row(self) -> dict[str, object]:
        """Return scalar fields suitable for display or CSV export."""
        row = {
            "file": self.path.name,
            "path": str(self.path),
            "material": self.material,
            "shape": self.shape,
            "shape_variant": self.shape_variant,
            "comparison_label": self.comparison_label,
            "preset": self.preset or "",
            "batch": self.batch or "",
            "mu0_Ms_T": self.mu0_Ms_T,
            "Ms_A_per_m": self.Ms_A_per_m,
            "A0_J_per_m": self.A0_J_per_m,
            "K0_J_per_m3": self.K0_J_per_m3,
            "L_m": self.L_m,
            "size_nm": self.size_nm,
            "size_factor": self.size_factor,
            "n": self.n,
            "ntot": self.ntot,
            "adaptive": self.adaptive,
            "periodic": self.periodic,
            "adaptive_dh_min_t": self.adaptive_dh_min_t,
            "use_fmm": self.use_fmm,
            "cuda": self.cuda,
            "cvode": self.cvode,
            "runtime_s": self.runtime_s,
        }
        row.update(self.metrics.as_dict())
        return row


def scalar(data: np.lib.npyio.NpzFile, key: str, default: object = np.nan) -> object:
    """Read an NPZ scalar while leaving genuine arrays unchanged."""
    if key not in data.files:
        return default
    value = data[key]
    return value.item() if np.asarray(value).shape == () else value


def infer_run_metadata_from_filename(path: Path) -> dict[str, object]:
    """Infer run modes for old files that did not save those fields."""
    stem = path.stem
    adaptive_match = _ADAPTIVE_PATTERN.search(stem)
    return {
        "adaptive": adaptive_match is not None,
        "periodic": bool(_PERIODIC_PATTERN.search(stem)),
        "adaptive_dh_min_t": (
            float(adaptive_match.group("dh")) if adaptive_match else np.nan
        ),
    }


def nearest_manifest(path: Path) -> tuple[Path | None, dict[str, object]]:
    """Return the closest ancestor manifest, if this is a batch result."""
    for parent in path.parents:
        candidate = parent / "manifest.json"
        if candidate.is_file():
            try:
                return candidate, json.loads(candidate.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                return candidate, {}
    return None, {}


def _material_label(
    preset: str | None,
    mu0_ms_t: float,
    a0: float,
    k0: float,
) -> str:
    if preset:
        return preset
    return f"Ms={mu0_ms_t:g}T, A0={a0:.3g}, K0={k0:.3g}"


def _shape_from_path(path: Path) -> str | None:
    """Infer a shape variant from a shape_* result directory."""
    for parent in path.parents:
        if parent.name.startswith("shape_"):
            variant = parent.name.removeprefix("shape_")
            return variant or None
    return None


def _dict_scalar(value: object) -> dict[str, Any]:
    """Normalize optional object-array metadata to a plain dict."""
    if isinstance(value, dict):
        return dict(value)
    return {}


def _shape_metadata(
    data: np.lib.npyio.NpzFile,
    manifest: dict[str, object],
    path: Path,
) -> tuple[str, str, dict[str, Any]]:
    """Resolve shape metadata from NPZ, manifest, path, then legacy default."""
    manifest_shape = manifest.get("shape") if isinstance(manifest, dict) else None
    path_shape = _shape_from_path(path)

    shape = scalar(data, "shape", manifest_shape or path_shape or "cube")
    shape_variant = scalar(data, "shape_variant", manifest_shape or path_shape or shape)
    metadata = scalar(data, "shape_metadata", {})

    shape_label = str(shape) if shape else "cube"
    variant_label = str(shape_variant) if shape_variant else shape_label
    return shape_label, variant_label, _dict_scalar(metadata)


def discover_result_files(
    roots: Path | str | Iterable[Path | str],
    pattern: str = "single_grain*.npz",
) -> list[Path]:
    """Recursively discover legacy and batch result files without duplicates."""
    if isinstance(roots, (str, Path)):
        roots = [roots]
    discovered: dict[Path, Path] = {}
    for root in roots:
        root_path = Path(root).expanduser()
        if root_path.is_file() and root_path.suffix == ".npz":
            discovered[root_path.resolve()] = root_path
        elif root_path.exists():
            for path in root_path.rglob(pattern):
                discovered[path.resolve()] = path
    return sorted(discovered.values(), key=lambda path: str(path))


def load_result(path: Path | str, *, include_res: bool = False) -> ResultRecord:
    """Load one NPZ file and calculate any derived metrics not saved by the run."""
    path = Path(path)
    inferred = infer_run_metadata_from_filename(path)
    manifest_path, manifest = nearest_manifest(path)
    material_manifest = manifest.get("material", {}) if isinstance(manifest, dict) else {}

    with np.load(path, allow_pickle=True) as data:
        H_T = np.asarray(scalar(data, "H_array"), dtype=float).ravel()
        H_A_per_m = np.asarray(
            scalar(data, "H_array_A_per_m", H_T / MU0), dtype=float
        ).ravel()
        Mz = np.asarray(scalar(data, "Mz_array", scalar(data, "M_array")), dtype=float).ravel()
        M_scalar = np.asarray(scalar(data, "M_array", Mz), dtype=float).ravel()
        Mx = np.asarray(scalar(data, "Mx_array", np.full_like(Mz, np.nan)), dtype=float).ravel()
        My = np.asarray(scalar(data, "My_array", np.full_like(Mz, np.nan)), dtype=float).ravel()

        Ms = float(scalar(data, "Ms", np.nan))
        mu0_ms_t = float(
            scalar(
                data,
                "mu0_Ms_T",
                material_manifest.get("mu0_Ms_T", MU0 * Ms),
            )
        )
        if not np.isfinite(Ms) and np.isfinite(mu0_ms_t):
            Ms = mu0_ms_t / MU0
        a0 = float(scalar(data, "A0", material_manifest.get("A0_J_per_m", np.nan)))
        k0 = float(scalar(data, "K0", material_manifest.get("K0_J_per_m3", np.nan)))
        L_m = float(scalar(data, "L"))
        preset_value = manifest.get("preset") if isinstance(manifest, dict) else None
        preset = str(preset_value) if preset_value else None
        batch_data = manifest.get("batch", {}) if isinstance(manifest, dict) else {}
        batch = batch_data.get("label") if isinstance(batch_data, dict) else None
        shape, shape_variant, shape_metadata = _shape_metadata(data, manifest, path)

        metrics = calculate_hysteresis_metrics(H_A_per_m, M_scalar, Ms)
        res = data["res"] if include_res and "res" in data.files else None
        adaptive = bool(scalar(data, "adaptive", inferred["adaptive"]))
        periodic = bool(scalar(data, "periodic", inferred["periodic"]))
        dh_min = float(
            scalar(data, "adaptive_dh_min_t", inferred["adaptive_dh_min_t"])
        )

        return ResultRecord(
            path=path,
            material=_material_label(preset, mu0_ms_t, a0, k0),
            preset=preset,
            batch=str(batch) if batch else (manifest_path.parent.name if manifest_path else None),
            shape=shape,
            shape_variant=shape_variant,
            shape_metadata=shape_metadata,
            mu0_Ms_T=mu0_ms_t,
            Ms_A_per_m=Ms,
            A0_J_per_m=a0,
            K0_J_per_m3=k0,
            L_m=L_m,
            size_nm=L_m * 1e9,
            size_factor=float(scalar(data, "size_factor", np.nan)),
            n=int(scalar(data, "n")),
            ntot=int(scalar(data, "ntot", int(scalar(data, "n")) ** 3)),
            adaptive=adaptive,
            periodic=periodic,
            adaptive_dh_min_t=dh_min,
            use_fmm=bool(scalar(data, "use_fmm", False)),
            cuda=bool(scalar(data, "cuda", False)),
            cvode=bool(scalar(data, "cvode", False)),
            runtime_s=float(scalar(data, "runtime", np.nan)),
            H_T=H_T,
            H_A_per_m=H_A_per_m,
            M_A_per_m=M_scalar,
            Mx_A_per_m=Mx,
            My_A_per_m=My,
            Mz_A_per_m=Mz,
            res=res,
            metrics=metrics,
        )


def load_results(
    roots: Path | str | Iterable[Path | str],
    *,
    include_res: bool = False,
    strict: bool = False,
) -> tuple[list[ResultRecord], list[tuple[Path, str]]]:
    """Load all discovered results and return non-fatal file errors separately."""
    records: list[ResultRecord] = []
    errors: list[tuple[Path, str]] = []
    for path in discover_result_files(roots):
        try:
            records.append(load_result(path, include_res=include_res))
        except Exception as error:
            if strict:
                raise
            errors.append((path, str(error)))
    records.sort(
        key=lambda item: (
            item.material,
            item.shape_variant,
            item.periodic,
            item.n,
            np.inf if not np.isfinite(item.adaptive_dh_min_t) else item.adaptive_dh_min_t,
            item.size_nm,
            str(item.path),
        )
    )
    return records, errors


def filter_results(
    records: Iterable[ResultRecord],
    *,
    n: int | None = None,
    adaptive_dh_min_t: float | None = None,
    periodic: bool | None = None,
    materials: Iterable[str] | None = None,
    shapes: Iterable[str] | None = None,
) -> list[ResultRecord]:
    """Filter without aggregating unlike numerical configurations."""
    selected_materials = set(materials) if materials is not None else None
    selected_shapes = set(shapes) if shapes is not None else None
    selected = []
    for record in records:
        if n is not None and record.n != n:
            continue
        if adaptive_dh_min_t is not None and not np.isclose(
            record.adaptive_dh_min_t, adaptive_dh_min_t
        ):
            continue
        if periodic is not None and record.periodic != periodic:
            continue
        if selected_materials is not None and record.material not in selected_materials:
            continue
        if selected_shapes is not None and record.shape_variant not in selected_shapes:
            continue
        selected.append(record)
    return selected


def write_summary_csv(records: Iterable[ResultRecord], path: Path | str) -> Path:
    """Write scalar result metadata and metrics without requiring pandas."""
    path = Path(path)
    rows = [record.summary_row() for record in records]
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return path
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def format_summary_table(records: Iterable[ResultRecord]) -> str:
    """Return a compact fixed-width table for notebook display."""
    rows = list(records)
    header = "material / shape        size[nm]  n  periodic  dh_min[T]  |mu0 Hc|[T]  mu0 Mr[T]  BHmax[kJ/m3]"
    lines = [header, "-" * len(header)]
    for record in rows:
        lines.append(
            f"{record.comparison_label:<22.22} {record.size_nm:8.3g} {record.n:2d} "
            f"{str(record.periodic):>8} {record.adaptive_dh_min_t:10.3g} "
            f"{abs(record.metrics.Hc_T):12.4g} {record.metrics.mu0_Mr_T:9.4g} "
            f"{record.metrics.BH_max_kJ_per_m3:13.4g}"
        )
    return "\n".join(lines)


__all__ = [
    "ResultRecord",
    "discover_result_files",
    "filter_results",
    "format_summary_table",
    "infer_run_metadata_from_filename",
    "load_result",
    "load_results",
    "nearest_manifest",
    "scalar",
    "write_summary_csv",
]
