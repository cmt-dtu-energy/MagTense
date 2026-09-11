"""
@brief Demonstrates fixed and adaptive MagTense hysteresis field stepping.

The example uses a small 5 x 5 x 5 single-grain problem.  Its 125 magnetic
cells allow spatial variation while keeping all five curves reasonably quick.
The fixed 0.01 T curve is treated as the reference solution.  Adaptive curves
are compared with this reference using coercivity, interpolation error, and
the number of accepted field states.

The code is intentionally explicit and lightly abstracted so the calculation
can be translated to MATLAB one block at a time.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys
import tempfile
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np


# Make direct execution work from either the repository root or this folder.
# An installed MagTense package remains usable because these paths merely point
# at the same source tree and compiled extension in a development checkout.
REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PYTHON_SOURCE = REPOSITORY_ROOT / "python" / "src"
for search_path in (REPOSITORY_ROOT, PYTHON_SOURCE):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

from magtense.micromag import MicromagProblem  # noqa: E402


# -----------------------------------------------------------------------------
# Easily varied physical and numerical parameters
# -----------------------------------------------------------------------------

MU0 = 4.0 * np.pi * 1.0e-7
DEFAULT_OUTPUT_PATH = Path(__file__).with_name("adaptive_hysteresis_comparison.png")


@dataclass(frozen=True)
class DemoConfig:
    """
    @brief Stores all parameters varied by the script and notebook.

    Fields use SI units except angles in degrees and field-step quantities,
    which are specified as mu0*H in tesla for readability.
    """

    # A 5 x 5 x 5 mesh contains 125 cells.  This permits non-uniform magnetic
    # states while remaining small enough for an interactive demonstration.
    grain_size_m: float = 10.0e-9
    cells_per_axis: int = 5

    # These are the standard single-grain material parameters used elsewhere in
    # the repository.  Ms is derived from the specified saturation induction.
    saturation_induction_t: float = 2.4
    anisotropy_j_per_m3: float = 1.0e6
    exchange_j_per_m: float = 7.0e-12
    field_tilt_deg: float = 3.0

    # A descending branch from +1 T to -2 T starts near positive saturation and
    # extends well beyond the reversal of this small uniaxial grain.
    field_start_t: float = 1.0
    field_end_t: float = -2.0
    fixed_steps_t: tuple[float, ...] = (0.5, 0.1, 0.01)

    # Adaptive stepping begins and may grow to 0.5 T.  The two requested minima
    # control the finest step the accept/reject loop may use near switching.
    adaptive_initial_step_t: float = 0.5
    adaptive_max_step_t: float = 0.5
    adaptive_min_steps_t: tuple[float, ...] = (0.1, 0.01)
    adaptive_max_accepted_steps: int = 512


# -----------------------------------------------------------------------------
# Shared problem construction
# -----------------------------------------------------------------------------

def unit_field_direction(tilt_deg: float) -> np.ndarray:
    """
    @brief Returns a unit field direction tilted from z towards x.

    @param tilt_deg Tilt angle in degrees.
    @return Cartesian unit vector [sin(theta), 0, cos(theta)].
    """
    theta = np.deg2rad(float(tilt_deg))
    return np.array([np.sin(theta), 0.0, np.cos(theta)], dtype=float)


def fixed_field_values(start_t: float, end_t: float, step_t: float) -> np.ndarray:
    """
    @brief Constructs an endpoint-inclusive descending field grid.

    @param start_t Initial mu0*H value in tesla.
    @param end_t Final mu0*H value in tesla.
    @param step_t Positive fixed step magnitude in tesla.
    @return Monotonically descending field values including both endpoints.
    @raises ValueError If the interval or step is invalid.
    """
    if not start_t > end_t:
        raise ValueError("field_start_t must be greater than field_end_t")
    if not np.isfinite(step_t) or step_t <= 0.0:
        raise ValueError("step_t must be finite and positive")

    span = float(start_t - end_t)
    interval_count = int(round(span / step_t))
    if not np.isclose(interval_count * step_t, span, rtol=0.0, atol=1.0e-12):
        raise ValueError("the fixed step must divide the selected field interval")
    return np.linspace(start_t, end_t, interval_count + 1)


def static_field_array(field_t: np.ndarray, direction: np.ndarray) -> np.ndarray:
    """
    @brief Builds MagTense's four-column static external-field array.

    Column zero stores the scalar mu0*H value in tesla for bookkeeping.  The
    Fortran solver reads columns one to three as the Cartesian H vector in A/m,
    hence the division by mu0.

    @param field_t Scalar mu0*H samples in tesla.
    @param direction Cartesian field unit vector.
    @return Array with shape (n_fields, 4).
    """
    field_t = np.asarray(field_t, dtype=float)
    array = np.zeros((field_t.size, 4), dtype=float)
    array[:, 0] = field_t
    array[:, 1:4] = (field_t / MU0)[:, np.newaxis] * direction[np.newaxis, :]
    return array


def create_problem(config: DemoConfig, hysteresis_solver: str) -> MicromagProblem:
    """
    @brief Creates a fresh uniformly magnetised problem for one sweep.

    A new problem is constructed for every curve.  This ensures that every
    fixed and adaptive run starts from the same positively saturated state and
    does not inherit the final state of a preceding calculation.

    @param config Physical and numerical demonstration parameters.
    @param hysteresis_solver Either ``"static"`` or ``"adaptive"``.
    @return Initialised MagTense micromagnetic problem.
    """
    n = int(config.cells_per_axis)
    n_tiles = n**3
    direction = unit_field_direction(config.field_tilt_deg)
    easy_axis = np.array([0.0, 0.0, 1.0], dtype=float)
    saturation_magnetisation = config.saturation_induction_t / MU0

    problem = MicromagProblem(
        res=[n, n, n],
        grid_L=[config.grain_size_m] * 3,
        grid_type="uniform",
        solver="explicit",
        hysteresis_solver=hysteresis_solver,
        m0=np.tile(direction, (n_tiles, 1)),
        A0=config.exchange_j_per_m,
        Ms=saturation_magnetisation,
        K0=config.anisotropy_j_per_m3,
        alpha=4000.0,
        gamma=0.0,
        cuda=False,
        cvode=False,
        useavgn=1,
        # This script plots the applied field from result[4], so the individual H
        # components have to be returned.
        usereturnhall=True,
    )

    # The easy axis is stored once per tile even though this example is uniform.
    problem.u_ea[:, :] = easy_axis[np.newaxis, :]

    # Two time samples are sufficient for the equilibrium/static solver used by
    # this demonstration and match the small single-grain reference setup.
    problem.t = np.linspace(0.0, 1.0e-9, 2)
    problem.nt = len(problem.t)
    problem.t_conv = problem.t.copy()
    problem.nt_conv = len(problem.t_conv)
    problem.exch_presize = problem.ntot * 28
    problem.use_fmm = 0
    problem.ifunif = 1
    problem.nlmin = 1
    problem.nlmax = 5

    # Solver log files are placed outside the repository so repeatedly running
    # this educational example does not dirty the working tree.
    log_dir = Path(tempfile.gettempdir()) / "magtense_adaptive_hysteresis"
    log_dir.mkdir(parents=True, exist_ok=True)
    problem.log_dir = str(log_dir)
    problem.window_enabled = 0
    problem.trace_enabled = 0
    return problem


# -----------------------------------------------------------------------------
# Fixed and adaptive solver calls
# -----------------------------------------------------------------------------

def extract_mean_mz_t(
    result: Sequence[object],
    n_fields: int,
    saturation_induction_t: float,
) -> np.ndarray:
    """
    @brief Extracts volume-averaged mu0*Mz from a solver result.

    ``result[1]`` stores magnetisation history.  Index zero selects the final
    returned time state, the next axis enumerates tiles, the third axis is the
    accepted field index, and the final axis contains x/y/z components.  The
    saved magnetisation is normalised, so multiplication by mu0*Ms gives tesla.

    @param result MagTense hysteresis result sequence.
    @param n_fields Number of fixed or accepted field states.
    @param saturation_induction_t Saturation induction mu0*Ms in tesla.
    @return Mean z magnetisation aligned with the field samples.
    """
    normalised_magnetisation = np.asarray(result[1][1, :, :n_fields, :], dtype=float)
    mean_normalised_mz = np.mean(normalised_magnetisation[:, :, 2], axis=0)
    return saturation_induction_t * mean_normalised_mz


def run_fixed_hysteresis(config: DemoConfig, step_t: float) -> dict[str, object]:
    """
    @brief Runs one conventional hysteresis sweep with a uniform field step.

    @param config Physical and numerical demonstration parameters.
    @param step_t Fixed mu0*H step magnitude in tesla.
    @return Simple result dictionary shared with adaptive runs.
    """
    direction = unit_field_direction(config.field_tilt_deg)
    field_t = fixed_field_values(config.field_start_t, config.field_end_t, step_t)
    problem = create_problem(config, "static")
    result = problem.run_hysteresis(H_ext=static_field_array(field_t, direction))
    magnetisation_t = extract_mean_mz_t(
        result,
        field_t.size,
        config.saturation_induction_t,
    )
    return make_run_record("fixed", step_t, field_t, magnetisation_t)


def run_adaptive_hysteresis(config: DemoConfig, minimum_step_t: float) -> dict[str, object]:
    """
    @brief Runs one hysteresis sweep using backend adaptive field stepping.

    The solver proposes a field step, solves the magnetic state, and accepts or
    retries the step according to the magnetisation change.  ``dH_min`` limits
    refinement near reversal; ``switch_refine_dH`` explicitly requests that
    same resolution when the switching region is detected.

    @param config Physical and numerical demonstration parameters.
    @param minimum_step_t Minimum accepted mu0*H step in tesla.
    @return Simple result dictionary containing only accepted field states.
    """
    direction = unit_field_direction(config.field_tilt_deg)
    problem = create_problem(config, "adaptive")

    # The adaptive interface accepts Cartesian H in A/m, not mu0*H in tesla.
    result = problem.run_hysteresis_adaptive(
        H_start=config.field_start_t / MU0 * direction,
        H_end=config.field_end_t / MU0 * direction,
        dH_initial=config.adaptive_initial_step_t / MU0,
        dH_min=minimum_step_t / MU0,
        dH_max=config.adaptive_max_step_t / MU0,
        max_steps=config.adaptive_max_accepted_steps,
        switch_refine_dH=minimum_step_t / MU0,
    )

    # Adaptive output is preallocated in Fortran.  The final result item is the
    # number of accepted states; the Python wrapper has already sliced every
    # field-dependent array to this length.
    n_accepted = int(result[-1])
    field_vectors_a_per_m = np.asarray(result[4][0, 0, :n_accepted, :], dtype=float)
    field_t = MU0 * (field_vectors_a_per_m @ direction)
    magnetisation_t = extract_mean_mz_t(
        result,
        n_accepted,
        config.saturation_induction_t,
    )
    return make_run_record("adaptive", minimum_step_t, field_t, magnetisation_t)


def make_run_record(
    mode: str,
    requested_step_t: float,
    field_t: np.ndarray,
    magnetisation_t: np.ndarray,
) -> dict[str, object]:
    """
    @brief Validates and packages one curve in the common result representation.

    @param mode ``"fixed"`` or ``"adaptive"``.
    @param requested_step_t Fixed step or adaptive minimum step in tesla.
    @param field_t Accepted mu0*H samples in tesla.
    @param magnetisation_t Mean mu0*Mz samples in tesla.
    @return Dictionary used by metrics, tables, and plots.
    @raises ValueError If solver output is empty, misaligned, or non-monotonic.
    """
    field_t = np.asarray(field_t, dtype=float).ravel()
    magnetisation_t = np.asarray(magnetisation_t, dtype=float).ravel()
    if field_t.size == 0 or field_t.size != magnetisation_t.size:
        raise ValueError("field and magnetisation arrays must be non-empty and aligned")
    if not np.all(np.diff(field_t) < 0.0):
        raise ValueError("the accepted hysteresis field must be strictly descending")
    return {
        "mode": mode,
        "requested_step_t": float(requested_step_t),
        "field_t": field_t,
        "magnetisation_t": magnetisation_t,
        "coercivity_t": interpolated_zero_crossing(field_t, magnetisation_t),
        "n_evaluations": int(field_t.size),
    }


def run_all_hysteresis(config: DemoConfig = DemoConfig()) -> list[dict[str, object]]:
    """
    @brief Runs the three fixed and two adaptive comparison curves.

    @param config Physical and numerical demonstration parameters.
    @return Runs ordered as fixed steps followed by adaptive minimum steps.
    """
    runs = [run_fixed_hysteresis(config, step) for step in config.fixed_steps_t]
    runs.extend(
        run_adaptive_hysteresis(config, minimum_step)
        for minimum_step in config.adaptive_min_steps_t
    )
    return runs


# -----------------------------------------------------------------------------
# Quantitative comparison with the fine fixed-step reference
# -----------------------------------------------------------------------------

def interpolated_zero_crossing(field_t: np.ndarray, values: np.ndarray) -> float:
    """
    @brief Linearly interpolates the first zero crossing in sweep order.

    If consecutive samples (H0, M0) and (H1, M1) bracket zero, the estimate is
    Hc = H0 - M0*(H1 - H0)/(M1 - M0).  This is directly portable to MATLAB.

    @param field_t Monotonically descending field samples in tesla.
    @param values Magnetisation samples aligned with ``field_t``.
    @return Signed coercive field in tesla, or NaN when no crossing exists.
    """
    field_t = np.asarray(field_t, dtype=float)
    values = np.asarray(values, dtype=float)
    exact = np.flatnonzero(values == 0.0)
    if exact.size:
        return float(field_t[int(exact[0])])
    brackets = np.flatnonzero(values[:-1] * values[1:] < 0.0)
    if not brackets.size:
        return float("nan")
    index = int(brackets[0])
    h0, h1 = field_t[index], field_t[index + 1]
    m0, m1 = values[index], values[index + 1]
    return float(h0 - m0 * (h1 - h0) / (m1 - m0))


def interpolate_to_reference(
    run: dict[str, object], reference_field_t: np.ndarray
) -> np.ndarray:
    """
    @brief Interpolates one descending curve onto the fine reference grid.

    ``numpy.interp`` requires its sample abscissa to increase, so both stored
    run arrays are reversed.  The reference grid itself may remain descending.

    @param run Result dictionary returned by a fixed or adaptive calculation.
    @param reference_field_t Fine fixed-step field grid in tesla.
    @return Interpolated mean mu0*Mz values in tesla.
    """
    field_t = np.asarray(run["field_t"], dtype=float)
    magnetisation_t = np.asarray(run["magnetisation_t"], dtype=float)
    return np.interp(reference_field_t, field_t[::-1], magnetisation_t[::-1])


def comparison_metrics(
    runs: Sequence[dict[str, object]],
    reference_step_t: float = 0.01,
) -> list[dict[str, float | int | str]]:
    """
    @brief Compares all curves with the finest fixed-step result.

    @param runs Fixed and adaptive result dictionaries.
    @param reference_step_t Fixed step identifying the numerical reference.
    @return One table row per run with cost and agreement metrics.
    @raises ValueError If the requested fixed reference is absent.
    """
    reference = next(
        (
            run
            for run in runs
            if run["mode"] == "fixed"
            and np.isclose(run["requested_step_t"], reference_step_t)
        ),
        None,
    )
    if reference is None:
        raise ValueError(f"fixed {reference_step_t:g} T reference run is missing")

    reference_field_t = np.asarray(reference["field_t"], dtype=float)
    reference_magnetisation_t = np.asarray(reference["magnetisation_t"], dtype=float)
    reference_hc_t = float(reference["coercivity_t"])
    reference_evaluations = int(reference["n_evaluations"])

    rows: list[dict[str, float | int | str]] = []
    for run in runs:
        interpolated = interpolate_to_reference(run, reference_field_t)
        absolute_error = np.abs(interpolated - reference_magnetisation_t)
        evaluations = int(run["n_evaluations"])
        rows.append(
            {
                "mode": str(run["mode"]),
                "step_t": float(run["requested_step_t"]),
                "evaluations": evaluations,
                "reduction_percent": 100.0 * (1.0 - evaluations / reference_evaluations),
                "coercivity_t": float(run["coercivity_t"]),
                "coercivity_difference_t": abs(float(run["coercivity_t"]) - reference_hc_t),
                "mean_absolute_error_t": float(np.mean(absolute_error)),
                "maximum_absolute_error_t": float(np.max(absolute_error)),
            }
        )
    return rows


def print_metrics_table(rows: Sequence[dict[str, float | int | str]]) -> None:
    """
    @brief Prints a dependency-free fixed-width correctness table.

    @param rows Output from :func:`comparison_metrics`.
    @return No value.
    """
    header = (
        "mode      step [T]  fields  reduction  Hc [T]    |dHc| [T]  "
        "mean |dM| [T]  max |dM| [T]"
    )
    print("\nComparison with the fixed 0.01 T reference")
    print(header)
    print("-" * len(header))
    for row in rows:
        print(
            f"{row['mode']:<9} "
            f"{float(row['step_t']):>8.3g}  "
            f"{int(row['evaluations']):>6d}  "
            f"{float(row['reduction_percent']):>8.1f}%  "
            f"{float(row['coercivity_t']):>8.5f}  "
            f"{float(row['coercivity_difference_t']):>10.3e}  "
            f"{float(row['mean_absolute_error_t']):>13.3e}  "
            f"{float(row['maximum_absolute_error_t']):>12.3e}"
        )


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def run_label(run: dict[str, object]) -> str:
    """@brief Returns a concise legend label for one run."""
    if run["mode"] == "fixed":
        return f"fixed dH = {float(run['requested_step_t']):g} T"
    return f"adaptive dH_min = {float(run['requested_step_t']):g} T"


def plot_comparison(
    runs: Sequence[dict[str, object]],
    reference_step_t: float = 0.01,
) -> tuple[plt.Figure, np.ndarray]:
    """
    @brief Plots curves, reference errors, and accepted field-step sizes.

    @param runs Fixed and adaptive result dictionaries.
    @param reference_step_t Fixed step identifying the numerical reference.
    @return Matplotlib figure and three axes.
    """
    reference = next(
        run
        for run in runs
        if run["mode"] == "fixed"
        and np.isclose(run["requested_step_t"], reference_step_t)
    )
    reference_field_t = np.asarray(reference["field_t"], dtype=float)
    reference_magnetisation_t = np.asarray(reference["magnetisation_t"], dtype=float)

    colours = plt.get_cmap("tab10")(np.arange(len(runs)))
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 11.0), sharex=True)

    for colour, run in zip(colours, runs):
        field_t = np.asarray(run["field_t"], dtype=float)
        magnetisation_t = np.asarray(run["magnetisation_t"], dtype=float)
        marker_spacing = max(1, field_t.size // 45)
        axes[0].plot(
            field_t,
            magnetisation_t,
            marker="o",
            markevery=marker_spacing,
            markersize=3.0,
            linewidth=1.3,
            color=colour,
            label=f"{run_label(run)} ({field_t.size} fields)",
        )

        interpolated = interpolate_to_reference(run, reference_field_t)
        axes[1].plot(
            reference_field_t,
            np.abs(interpolated - reference_magnetisation_t),
            linewidth=1.3,
            color=colour,
            label=run_label(run),
        )

        # Associate each accepted increment with its endpoint.  Fixed runs form
        # horizontal lines; adaptive runs reveal refinement around switching.
        axes[2].plot(
            field_t[1:],
            np.abs(np.diff(field_t)),
            marker="o",
            markersize=3.0,
            linewidth=1.2,
            color=colour,
            label=run_label(run),
        )

    axes[0].axhline(0.0, color="0.45", linewidth=0.7)
    axes[0].axvline(0.0, color="0.45", linewidth=0.7)
    axes[0].set_ylabel(r"Mean $\mu_0 M_z$ [T]")
    axes[0].set_title("Fixed and adaptive hysteresis stepping")
    axes[0].legend(fontsize=8, ncol=2)

    axes[1].set_ylabel(r"$|\mu_0 M_z - \mu_0 M_{z,ref}|$ [T]")
    axes[1].set_title("Piecewise-linear error relative to fixed dH = 0.01 T")

    axes[2].set_yscale("log")
    axes[2].set_ylabel(r"Accepted $|\Delta(\mu_0 H)|$ [T]")
    axes[2].set_xlabel(r"Applied field $\mu_0 H$ [T]")
    axes[2].set_title("Adaptive steps become small near magnetisation reversal")

    for axis in axes:
        axis.grid(True, linestyle=":", alpha=0.65)
    fig.tight_layout()
    return fig, axes


def run_demonstration(
    config: DemoConfig = DemoConfig(),
    *,
    output_path: Path | str = DEFAULT_OUTPUT_PATH,
    show: bool = True,
) -> tuple[list[dict[str, object]], list[dict[str, float | int | str]], plt.Figure]:
    """
    @brief Runs, reports, plots, and saves the complete demonstration.

    @param config Physical and numerical demonstration parameters.
    @param output_path Destination PNG path.
    @param show Whether to call ``plt.show()`` after saving.
    @return Runs, metric rows, and the generated figure.
    """
    runs = run_all_hysteresis(config)
    rows = comparison_metrics(runs)
    print_metrics_table(rows)
    figure, _ = plot_comparison(runs)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=180, bbox_inches="tight")
    print(f"\nSaved comparison figure to {output_path}")
    if show:
        plt.show()
    return runs, rows, figure


def main() -> None:
    """@brief Parses display/output options and runs the default example."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help="PNG output path (default: beside this script)",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Save the figure without opening an interactive window",
    )
    arguments = parser.parse_args()
    run_demonstration(output_path=arguments.output, show=not arguments.no_show)


if __name__ == "__main__":
    main()
