from __future__ import annotations

import argparse
from dataclasses import dataclass
import time
from pathlib import Path

import numpy as np

from magtense.micromag import MicromagProblem

from utils.metrics import calculate_hysteresis_metrics
from utils.simulation import (
    DEFAULT_A0,
    DEFAULT_K0,
    DEFAULT_MS,
    DEFAULT_MU0_MS_T,
    DEFAULT_TILT_DEGREES,
    MU0,
    build_hysteresis_field,
    characteristic_length,
    create_single_grain_problem,
    extract_mean_magnetisation,
    interpolated_coercivity,
    output_stem,
    tilted_field_direction,
    z_easy_axis,
)
from utils.geometry import (
    PrismMesh,
    generate_cylinder,
    generate_ellipsoid,
    generate_hexagonal_prism,
    generate_rectangular_prism,
    generate_sphere,
)


SHAPE_ALIASES = {
    "cube": ("cube", "cube"),
    "box": ("cube", "cube"),
    "rect": ("rectangle", "rectangle"),
    "rectangle": ("rectangle", "rectangle"),
    "recantgle": ("rectangle", "rectangle"),
    "rectangular_prism": ("rectangle", "rectangle"),
    "sphere": ("sphere", "sphere"),
    "cylinder": ("cylinder", "cylinder_z"),
    "ellipsoid": ("ellipsoid", "ellipsoid_z"),
    "elipsoid": ("ellipsoid", "ellipsoid_z"),
    "ellipsoid_x": ("ellipsoid", "ellipsoid_x"),
    "elipsoid_x": ("ellipsoid", "ellipsoid_x"),
    "ellipsoid_z": ("ellipsoid", "ellipsoid_z"),
    "elipsoid_z": ("ellipsoid", "ellipsoid_z"),
    "hexagon": ("hexagon", "hexagon"),
    "hexagonal_prism": ("hexagon", "hexagon"),
}


@dataclass(frozen=True)
class ShapeSpec:
    """Resolved geometric shape parameters for one single-grain run."""

    shape: str
    variant: str
    outer_size_m: float
    parameters: dict[str, float | str | tuple[float, ...]]


def resolve_shape_spec(
    shape: str,
    shape_variant: str | None,
    outer_size_m: float,
) -> ShapeSpec:
    """Normalize aliases and convert outer size to concrete shape parameters."""
    shape_key = shape.lower().replace("-", "_")
    if shape_key not in SHAPE_ALIASES:
        valid = ", ".join(sorted(SHAPE_ALIASES))
        raise ValueError(f"Unknown shape {shape!r}. Expected one of: {valid}")
    base_shape, default_variant = SHAPE_ALIASES[shape_key]
    variant = (shape_variant or default_variant).lower().replace("-", "_")
    if variant in SHAPE_ALIASES:
        alias_shape, alias_variant = SHAPE_ALIASES[variant]
        if alias_shape != base_shape and shape_variant is not None:
            raise ValueError(
                f"Shape variant {shape_variant!r} does not belong to shape {shape!r}"
            )
        base_shape, variant = alias_shape, alias_variant

    L = float(outer_size_m)
    if not np.isfinite(L) or L <= 0.0:
        raise ValueError("outer_size_m must be positive and finite")

    if base_shape == "cube":
        variant = "cube"
        parameters = {"side_length": L}
    elif base_shape == "rectangle":
        variant = "rectangle"
        parameters = {"dimensions": (L, 0.7 * L, 0.7 * L)}
    elif base_shape == "sphere":
        variant = "sphere"
        parameters = {"radius": L / 2.0}
    elif base_shape == "cylinder":
        variant = "cylinder_z"
        parameters = {"radius": L / 2.0, "length": L, "axis": "z"}
    elif base_shape == "ellipsoid":
        if variant not in {"ellipsoid_x", "ellipsoid_z", "x", "z"}:
            raise ValueError("--shape-variant for ellipsoid must be x or z")
        variant = "ellipsoid_x" if variant in {"ellipsoid_x", "x"} else "ellipsoid_z"
        if variant == "ellipsoid_x":
            semi_axes = (L / 2.0, 0.35 * L, 0.35 * L)
        else:
            semi_axes = (0.35 * L, 0.35 * L, L / 2.0)
        parameters = {"semi_axes": semi_axes}
    elif base_shape == "hexagon":
        variant = "hexagon"
        parameters = {"side_length": L / 2.0, "height": L, "axis": "z"}
    else:
        raise AssertionError(f"Unhandled shape: {base_shape}")

    return ShapeSpec(
        shape=base_shape,
        variant=variant,
        outer_size_m=L,
        parameters=parameters,
    )


def build_shape_mesh(shape_spec: ShapeSpec, n: int) -> PrismMesh | None:
    """Return a mesh for non-cube shapes, or None for the legacy cube path."""
    if shape_spec.shape == "cube":
        return None
    target_tiles = n**3
    common = {
        "target_tiles": target_tiles,
        "overshoot_policy": "soft",
        "queue_policy": "symmetric_priority",
        "grid_shifts": "half_step",
        "max_depth": 8,
        "eps": 1e-12,
    }
    params = shape_spec.parameters
    if shape_spec.shape == "rectangle":
        dimensions = params["dimensions"]
        return generate_rectangular_prism(dimensions, n, n, n)
    if shape_spec.shape == "sphere":
        return generate_sphere(params["radius"], **common)
    if shape_spec.shape == "cylinder":
        return generate_cylinder(
            params["radius"], params["length"], params["axis"], **common
        )
    if shape_spec.shape == "ellipsoid":
        return generate_ellipsoid(params["semi_axes"], **common)
    if shape_spec.shape == "hexagon":
        return generate_hexagonal_prism(
            params["side_length"], params["height"], params["axis"], **common
        )
    raise AssertionError(f"Unhandled shape: {shape_spec.shape}")


def _json_safe(value):
    """Convert nested numpy-heavy metadata to simple serializable objects."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def run_single_grain_coercivity(
    *,
    L: float,
    size_factor: float,
    n: int,
    use_fmm: bool = False,
    periodic: bool = False,
    cuda: bool = False,
    cvode: bool = False,
    tilt_degrees: float = DEFAULT_TILT_DEGREES,
    field_min_t: float = -7.0,
    field_max_t: float = 1.0,
    field_step_t: float = -0.1,
    output_dir: Path = Path("results"),
    timer_log_dir: Path = Path("timer_logs_single_grain"),
    plotting: bool = True,
    Ms: float = DEFAULT_MS,
    K0: float = DEFAULT_K0,
    A0: float = DEFAULT_A0,
    fmm_cells_per_node: int = 660,
    fmm_eps: float = 1e-4,
    ifunif: int = 1,
    nlmin: int = 1,
    nlmax: int = 5,
    allow_fmm_short_circuit: int = 0,
    fmm_min_n: int = 20000,
    fmm_nterms: int = -1,
    adaptive: bool = False,
    adaptive_max_steps: int | None = None,
    adaptive_dh_initial_t: float | None = None,
    adaptive_dh_min_t: float | None = None,
    adaptive_dh_max_t: float | None = None,
    output_stem_override: str | None = None,
    shape: str = "cube",
    shape_variant: str | None = None,
    material_label: str | None = None,
    sw_coercivity_t: float | None = None,
) -> None:
    """Run one size/resolution hysteresis curve and save compatible arrays."""
    timer_log_dir.mkdir(parents=True, exist_ok=True)
    steps_t = np.arange(field_max_t, field_min_t + 0.5 * field_step_t, field_step_t)
    easy_axis = z_easy_axis()
    field_direction = tilted_field_direction(tilt_degrees)
    H_ext = build_hysteresis_field(steps_t, field_direction)
    shape_spec = resolve_shape_spec(shape, shape_variant, L)
    mesh = build_shape_mesh(shape_spec, n)
    if mesh is None:
        problem = create_single_grain_problem(
            L,
            n,
            use_fmm=use_fmm,
            cuda=cuda,
            cvode=cvode,
            tilt_degrees=tilt_degrees,
            Ms=Ms,
            K0=K0,
            A0=A0,
            fmm_cells_per_node=fmm_cells_per_node,
            fmm_eps=fmm_eps,
            ifunif=ifunif,
            nlmin=nlmin,
            nlmax=nlmax,
            allow_fmm_short_circuit=allow_fmm_short_circuit,
            fmm_min_n=fmm_min_n,
            fmm_nterms=fmm_nterms,
            timer_log_dir=timer_log_dir,
            hysteresis_solver="adaptive" if adaptive else "static",
            easy_axis=easy_axis,
            m0_direction=field_direction,
        )
        grid_L = np.array([L, L, L], dtype=float)
        target_tiles = n**3
        achieved_tiles = n**3
        mesh_metadata = {}
        shape_metadata = _json_safe(shape_spec.parameters)
        root_bounds = np.vstack((-grid_L / 2.0, grid_L / 2.0))
        represented_volume = L**3
    else:
        grid_kwargs = mesh.to_micromag_kwargs()
        grid_L = list(mesh.root_bounds[1] - mesh.root_bounds[0])
        problem = MicromagProblem(
            grid_L=grid_L,
            solver="explicit",
            hysteresis_solver="adaptive" if adaptive else "static",
            m0=np.tile(field_direction, (mesh.achieved_tiles, 1)),
            A0=A0,
            Ms=Ms,
            K0=K0,
            alpha=4000.0,
            gamma=0.0,
            cuda=cuda,
            cvode=cvode,
            useavgn=1,
            **grid_kwargs,
        )
        problem.u_ea[:, :] = easy_axis[np.newaxis, :]
        problem.t = np.linspace(0.0, 1e-9, 2)
        problem.nt = len(problem.t)
        problem.t_conv = np.linspace(0.0, 1e-9, 2)
        problem.nt_conv = len(problem.t_conv)
        problem.exch_presize = problem.ntot * 28
        problem.use_fmm = int(use_fmm)
        problem.fmm_cells_per_node = fmm_cells_per_node if use_fmm else 0
        problem.fmm_eps = fmm_eps
        problem.ifunif = ifunif
        problem.nlmin = nlmin
        problem.nlmax = nlmax
        problem.allow_fmm_short_circuit = allow_fmm_short_circuit
        problem.fmm_min_n = fmm_min_n
        problem.fmm_nterms = fmm_nterms
        problem.log_dir = str(timer_log_dir)
        problem.window_enabled = 1
        problem.window_interval = 30.0
        problem.trace_enabled = 0
        problem.flush_each = 1
        problem.trace_verbose = 1
        target_tiles = n**3
        achieved_tiles = mesh.achieved_tiles
        mesh_metadata = _json_safe(mesh.refinement_metadata)
        shape_metadata = _json_safe(mesh.shape_metadata)
        root_bounds = mesh.root_bounds
        represented_volume = mesh.represented_volume

    # Configure periodic exchange boundary conditions at the script level
    # (do not modify core MicromagProblem implementation here).
    if periodic and shape_spec.shape != "cube":
        raise ValueError("--periodic is only supported for cube runs")
    if periodic:
        try:
            problem.exchPBC = [1, 1, 1]
        except Exception:
            # Best-effort: if attribute not present, set attribute anyway
            setattr(problem, "exchPBC", [1, 1, 1])
    else:
        try:
            problem.exchPBC = [0, 0, 0]
        except Exception:
            setattr(problem, "exchPBC", [0, 0, 0])

    # Also set periodic repetition for the demag field (n_macro):
    # use a single periodic repetition in each dimension when periodic=True.
    if periodic:
        try:
            problem.n_macro = np.array([1, 1, 1], dtype=int)
            problem.shiftVec = np.array(grid_L)
        except Exception:
            setattr(problem, "n_macro", np.array([1, 1, 1], dtype=int))
            setattr(problem, "shiftVec", np.array(grid_L))
    else:
        try:
            problem.n_macro = np.zeros(3, dtype=int)
            problem.shiftVec = np.zeros(3)
        except Exception:
            setattr(problem, "n_macro", np.zeros(3, dtype=int))
            setattr(problem, "shiftVec", np.zeros(3))

    if output_stem_override is None:
        stem = output_stem(size_factor, n, use_fmm, fmm_nterms, nlmax)
        if shape_spec.variant != "cube":
            stem = stem.replace("single_grain_", f"single_grain_{shape_spec.variant}_", 1)
        # Append _P to the stem when periodic boundary conditions are enabled.
        if periodic:
            n_token = f"_n{n}"
            if n_token in stem:
                stem = stem.replace(n_token, n_token + "_P", 1)
            else:
                stem = stem + "_P"
        if adaptive:
            dH_min_t = (
                abs(adaptive_dh_min_t)
                if adaptive_dh_min_t is not None
                else abs(field_step_t / 10.0)
            )
            stem += f"_A_FS{dH_min_t:.1e}"
    else:
        stem = output_stem_override
    problem.timer_log_file = f"{stem}_timer.log"
    problem.trace_log_file = f"{stem}_trace.log"

    print("Single-grain coercivity experiment")
    print(f"  L = {L:.6e} m")
    print(f"  shape = {shape_spec.shape} ({shape_spec.variant})")
    print(f"  size_factor = {size_factor:.2g}")
    print(f"  n = {n} (target {target_tiles} cells, achieved {achieved_tiles})")
    material_length = characteristic_length(A0, Ms)
    print(f"  mu0 Ms = {MU0 * Ms:.6e} T")
    print(f"  Ms = {Ms:.6e} A/m")
    print(f"  A0 = {A0:.6e} J/m")
    print(f"  K0 = {K0:.6e} J/m^3")
    print(f"  characteristic_length = {material_length:.6e} m")
    print(f"  L / characteristic_length = {L / material_length:.6g}")
    print("  easy_axis_tilt = 0.000 degrees")
    print(f"  easy_axis = {easy_axis}")
    print(f"  field_angle = {tilt_degrees:.3f} degrees")
    print(f"  field_direction = {field_direction}")
    print(f"  initial_magnetisation = {field_direction}")
    print(f"  use_fmm = {use_fmm}")
    print(f"  cuda = {cuda}")
    print(f"  cvode = {cvode}")
    print(f"  adaptive = {adaptive}")
    start_time = time.time()
    resolved_adaptive_dh_min_t = np.nan
    if adaptive:
        max_steps = (
            adaptive_max_steps if adaptive_max_steps is not None else len(steps_t)
        )
        dH_initial_t = (
            adaptive_dh_initial_t
            if adaptive_dh_initial_t is not None
            else field_step_t
        )
        dH_min_t = (
            adaptive_dh_min_t
            if adaptive_dh_min_t is not None
            else field_step_t / 10.0
        )
        dH_max_t = (
            adaptive_dh_max_t
            if adaptive_dh_max_t is not None
            else field_step_t * 2.0
        )
        dH_initial = abs(dH_initial_t) / MU0
        dH_min = abs(dH_min_t) / MU0
        dH_max = abs(dH_max_t) / MU0
        resolved_adaptive_dh_min_t = abs(dH_min_t)
        res = problem.run_hysteresis_adaptive(
            H_start=H_ext[0, 1:4],
            H_end=H_ext[-1, 1:4],
            dH_initial=dH_initial,
            dH_min=dH_min,
            dH_max=dH_max,
            max_steps=max_steps,
            switch_refine_dH=dH_min,
        )
        n_fields = res[-1]
        H_vectors_A_per_m = res[4][0, 0, :n_fields, :]
        H_A_per_m = H_vectors_A_per_m @ field_direction
    else:
        res = problem.run_hysteresis(H_ext=H_ext)
        n_fields = len(steps_t)
        H_vectors_A_per_m = H_ext[:, 1:4]
        H_A_per_m = H_vectors_A_per_m @ field_direction
    runtime = time.time() - start_time

    Mx, My, Mz = extract_mean_magnetisation(res, n_fields, Ms)
    M_parallel = (
        np.column_stack((Mx, My, Mz)) @ field_direction
    )
    H_T = MU0 * H_A_per_m
    metrics = calculate_hysteresis_metrics(H_A_per_m, M_parallel, Ms)
    Hc_A_per_m = metrics.Hc_A_per_m
    Hc_T = metrics.Hc_T
    anisotropy_field_A_per_m = 2.0 * K0 / Ms

    output_dir.mkdir(parents=True, exist_ok=True)
    res_path = output_dir / f"{stem}.npz"
    np.savez(
        res_path,
        res=np.array(res, dtype=object),
        H_array=H_T,
        H_array_A_per_m=H_A_per_m,
        Hx_array_A_per_m=H_vectors_A_per_m[:, 0],
        Hy_array_A_per_m=H_vectors_A_per_m[:, 1],
        Hz_array_A_per_m=H_vectors_A_per_m[:, 2],
        M_array=M_parallel,
        Mx_array=Mx,
        My_array=My,
        Mz_array=Mz,
        M_parallel_array=M_parallel,
        Hc=Hc_A_per_m,
        Hc_A_per_m=Hc_A_per_m,
        Hc_T=Hc_T,
        Mr_A_per_m=metrics.Mr_A_per_m,
        mu0_Mr_T=metrics.mu0_Mr_T,
        Mr_over_Ms=metrics.Mr_over_Ms,
        BH_max_J_per_m3=metrics.BH_max_J_per_m3,
        BH_max_kJ_per_m3=metrics.BH_max_kJ_per_m3,
        BH_max_MGOe=metrics.BH_max_MGOe,
        coercivity_status=metrics.coercivity_status,
        remanence_status=metrics.remanence_status,
        energy_product_status=metrics.energy_product_status,
        H_N=anisotropy_field_A_per_m,
        runtime=runtime,
        L=L,
        n=n,
        ntot=achieved_tiles,
        shape=shape_spec.shape,
        shape_variant=shape_spec.variant,
        shape_outer_size_m=shape_spec.outer_size_m,
        shape_parameters=np.array(_json_safe(shape_spec.parameters), dtype=object),
        shape_metadata=np.array(shape_metadata, dtype=object),
        mesh_metadata=np.array(mesh_metadata, dtype=object),
        mesh_root_bounds=np.asarray(root_bounds, dtype=float),
        represented_volume=represented_volume,
        target_tiles=target_tiles,
        achieved_tiles=achieved_tiles,
        use_fmm=use_fmm,
        cuda=cuda,
        cvode=cvode,
        fmm_cells_per_node=fmm_cells_per_node,
        fmm_eps=fmm_eps,
        ifunif=ifunif,
        nlmin=nlmin,
        nlmax=nlmax,
        allow_fmm_short_circuit=allow_fmm_short_circuit,
        fmm_min_n=fmm_min_n,
        fmm_nterms=fmm_nterms,
        characteristic_length=material_length,
        size_factor=L / material_length,
        adaptive=adaptive,
        periodic=periodic,
        adaptive_dh_min_t=resolved_adaptive_dh_min_t,
        easy_axis=easy_axis,
        easy_axis_tilt_degrees=0.0,
        field_direction=field_direction,
        field_angle_degrees=tilt_degrees,
        initial_magnetisation=field_direction,
        magnetisation_projection="field_direction",
        material_label=material_label or "",
        sw_coercivity_t=(
            float(sw_coercivity_t) if sw_coercivity_t is not None else np.nan
        ),
        field_min_t=field_min_t,
        field_max_t=field_max_t,
        mu0_Ms_T=MU0 * Ms,
        Ms=Ms,
        K0=K0,
        A0=A0,
        output_stem=stem,
    )
    print(f"Hysteresis simulation took {runtime:.3f} seconds")
    print(f"Hc = {Hc_A_per_m:.6e} A/m ({Hc_T:.6e} T)")
    print(
        f"Mr = {metrics.Mr_A_per_m:.6e} A/m "
        f"(mu0 Mr = {metrics.mu0_Mr_T:.6e} T)"
    )
    print(
        f"(BH)max = {metrics.BH_max_kJ_per_m3:.6e} kJ/m^3 "
        f"({metrics.BH_max_MGOe:.6e} MGOe)"
    )
    print(f"Saved result object and scalars to: {res_path}")

    if plotting:
        plt = __import__("matplotlib.pyplot", fromlist=["pyplot"])

        fig_path = output_dir / f"{stem}.png"
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(MU0 * H_A_per_m, MU0 * M_parallel, ".-k", linewidth=1.5, markersize=4)
        ax.plot(
            MU0 * Hc_A_per_m,
            0.0,
            "r*",
            markersize=10,
            label=f"Hc = {Hc_T:.3f} T",
        )
        ax.set_xlabel(r"Applied Field, $\mu_0 H$ [ T ]")
        ax.set_ylabel(r"Magnetization parallel to field, $\mu_0 M$ [ T ]")
        ax.set_title(
            f"Single grain hysteresis: size_factor={size_factor:.2g}, n={n}, FMM={use_fmm}"
        )
        ax.grid(True, which="both", linestyle=":", linewidth=0.8, alpha=0.6)
        ax.legend(loc="lower right", frameon=True, framealpha=0.9, edgecolor="inherit")
        fig.tight_layout(pad=0.2)
        fig.savefig(fig_path, dpi=300)
        plt.close(fig)
        print(f"Saved figure to: {fig_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run a simple single-grain grain-size/resolution coercivity test. "
            "With no arguments, this runs size_factor=1, n=10, FMM/CUDA/CVODE off, "
            "and writes the result plus plot to ./results."
        )
    )
    size_group = parser.add_mutually_exclusive_group()
    size_group.add_argument("--L", type=float, help="Cubic domain side length [m].")
    size_group.add_argument(
        "--size-factor",
        type=float,
        help=(
            "Set L = size_factor * sqrt(A0 / (0.5 * mu0 * Ms**2)); "
            "defaults to 1.0 when neither --L nor --size-factor is provided."
        ),
    )
    parser.add_argument(
        "--n", type=int, default=10, help="Cubic grid resolution (default: 10)."
    )
    parser.add_argument(
        "--tilt-degrees",
        type=float,
        default=DEFAULT_TILT_DEGREES,
        help="External-field tilt away from z in the x-z plane [degrees].",
    )
    parser.add_argument("--field-min-t", type=float, default=-3.0)
    parser.add_argument("--field-max-t", type=float, default=1.0)
    parser.add_argument("--field-step-t", type=float, default=-0.1)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument(
        "--output-stem",
        help=(
            "Override the complete result/plot/log stem. Intended for sweep "
            "submitters that use physical-size names."
        ),
    )
    parser.add_argument(
        "--timer-log-dir",
        type=Path,
        default=Path("timer_logs_single_grain"),
        help="Directory for backend timer and trace logs.",
    )
    parser.add_argument(
        "--mu0-ms-t",
        type=float,
        default=DEFAULT_MU0_MS_T,
        help="Saturation magnetization expressed as mu0*Ms [T].",
    )
    parser.add_argument(
        "--a0", type=float, default=DEFAULT_A0, help="Exchange stiffness [J/m]."
    )
    parser.add_argument(
        "--k0",
        type=float,
        default=DEFAULT_K0,
        help="Uniaxial anisotropy constant [J/m^3].",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Disable the default hysteresis plot written to the output directory.",
    )
    parser.add_argument("--use-fmm", action="store_true", help="Enable FMM demag.")
    parser.add_argument("--cuda", action="store_true", help="Enable CUDA backend support.")
    parser.add_argument("--fmm-cpn", type=int, default=660)
    parser.add_argument("--fmm-eps", type=float, default=1e-4)
    parser.add_argument("--ifunif", type=int, default=1)
    parser.add_argument("--nlmin", type=int, default=1)
    parser.add_argument("--nlmax", type=int, default=5)
    parser.add_argument("--allow-fmm-short-circuit", type=int, default=0)
    parser.add_argument("--fmm-min-n", type=int, default=20000)
    parser.add_argument("--fmm-nterms", type=int, default=-1)
    parser.add_argument(
        "--adaptive",
        action="store_true",
        help="Use adaptive backend hysteresis stepping.",
    )
    parser.add_argument("--adaptive-max-steps", type=int, default=500)
    parser.add_argument("--adaptive-dh-initial-t", type=float, default=0.5)
    parser.add_argument("--adaptive-dh-min-t", type=float, default=0.01)
    parser.add_argument("--adaptive-dh-max-t", type=float, default=0.5)
    parser.add_argument(
        "--periodic",
        action="store_true",
        help="Enable periodic exchange boundary conditions (exchPBC=1) and append _P to output names.",
    )
    parser.add_argument(
        "--shape",
        default="cube",
        help=(
            "Shape to simulate: cube, rectangle, sphere, cylinder, ellipsoid, "
            "or hexagon. Common aliases and misspellings are accepted."
        ),
    )
    parser.add_argument(
        "--shape-variant",
        help="Optional variant, currently x/z for ellipsoid orientation.",
    )
    parser.add_argument(
        "--material-label",
        help="Optional metadata label for custom/property-grid sweeps.",
    )
    parser.add_argument(
        "--sw-coercivity-t",
        type=float,
        help="Optional Stoner-Wohlfarth coercivity estimate stored as metadata [T].",
    )
    args = parser.parse_args()

    if not np.isfinite(args.mu0_ms_t) or args.mu0_ms_t <= 0.0:
        parser.error("--mu0-ms-t must be positive")
    if not np.isfinite(args.a0) or args.a0 <= 0.0:
        parser.error("--a0 must be positive")
    if not np.isfinite(args.k0) or args.k0 < 0.0:
        parser.error("--k0 must be non-negative")
    if args.n <= 0:
        parser.error("--n must be positive")
    if args.size_factor is not None and args.size_factor <= 0.0:
        parser.error("--size-factor must be positive")
    if args.L is not None and args.L <= 0.0:
        parser.error("--L must be positive")
    if args.output_stem is not None and (
        not args.output_stem or Path(args.output_stem).name != args.output_stem
    ):
        parser.error("--output-stem must be a non-empty filename stem")
    try:
        shape_spec = resolve_shape_spec(args.shape, args.shape_variant, 1.0)
    except ValueError as error:
        parser.error(str(error))
    if args.periodic and shape_spec.shape != "cube":
        parser.error("--periodic is only supported for cube runs")

    Ms = args.mu0_ms_t / MU0
    material_length = characteristic_length(args.a0, Ms)

    if args.size_factor is None:
        if args.L is None:
            size_factor = 1.0
        else:
            size_factor = args.L / material_length
    else:
        size_factor = args.size_factor
    L = args.L if args.L is not None else size_factor * material_length
    run_single_grain_coercivity(
        L=L,
        size_factor=size_factor,
        n=args.n,
        use_fmm=args.use_fmm,
        cuda=args.cuda,
        cvode=False,
        tilt_degrees=args.tilt_degrees,
        field_min_t=args.field_min_t,
        field_max_t=args.field_max_t,
        field_step_t=args.field_step_t,
        output_dir=args.output_dir,
        timer_log_dir=args.timer_log_dir,
        plotting=not args.no_plot,
        Ms=Ms,
        K0=args.k0,
        A0=args.a0,
        periodic=args.periodic,
        fmm_cells_per_node=args.fmm_cpn,
        fmm_eps=args.fmm_eps,
        ifunif=args.ifunif,
        nlmin=args.nlmin,
        nlmax=args.nlmax,
        allow_fmm_short_circuit=args.allow_fmm_short_circuit,
        fmm_min_n=args.fmm_min_n,
        fmm_nterms=args.fmm_nterms,
        adaptive=args.adaptive,
        adaptive_max_steps=args.adaptive_max_steps,
        adaptive_dh_initial_t=args.adaptive_dh_initial_t,
        adaptive_dh_min_t=args.adaptive_dh_min_t,
        adaptive_dh_max_t=args.adaptive_dh_max_t,
        output_stem_override=args.output_stem,
        shape=args.shape,
        shape_variant=args.shape_variant,
        material_label=args.material_label,
        sw_coercivity_t=args.sw_coercivity_t,
    )


if __name__ == "__main__":
    main()
