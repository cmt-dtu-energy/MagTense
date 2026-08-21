"""
@brief Runs one minimal adaptive hysteresis experiment.

This example deliberately contains no plotting, comparison, interpolation, or
helper functions.  The accepted field values and volume-averaged
magnetisation remain available as ordinary NumPy arrays after the solve.
"""

from pathlib import Path
import sys
import tempfile

import numpy as np


# Allow this file to be run directly from a MagTense source checkout.
REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
PYTHON_SOURCE = REPOSITORY_ROOT / "python" / "src"
for search_path in (REPOSITORY_ROOT, PYTHON_SOURCE):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

from magtense.micromag import MicromagProblem  # noqa: E402


# All material parameters use SI units.  Field values and field-step sizes are
# written below as mu0*H in tesla and converted to H in A/m for MagTense.
MU0 = 4.0 * np.pi * 1.0e-7
BS = 2.4
MS = BS / MU0
K0 = 1.0e6
A0 = 7.0e-12

# The 10 nm cube is divided into 5 x 5 x 5 = 125 magnetic cells.
N = 5
N_TILES = N**3
GRAIN_SIZE = 10.0e-9

# Sweep downwards from positive saturation.  The adaptive solver starts with a
# 0.5 T step and may refine it to 0.01 T near a rapid magnetisation change.
FIELD_START_T = 1.0
FIELD_END_T = -2.0
INITIAL_STEP_T = 0.5
MINIMUM_STEP_T = 0.01
MAXIMUM_STEP_T = 0.5
MAXIMUM_ACCEPTED_STEPS = 512

# A small tilt avoids exact alignment of the field and uniaxial easy axis.
tilt = np.deg2rad(3.0)
field_direction = np.array([np.sin(tilt), 0.0, np.cos(tilt)])

# Initialise every cell along the positive applied-field direction.
problem = MicromagProblem(
    res=[N, N, N],
    grid_L=[GRAIN_SIZE, GRAIN_SIZE, GRAIN_SIZE],
    grid_type="uniform",
    solver="explicit",
    hysteresis_solver="adaptive",
    m0=np.tile(field_direction, (N_TILES, 1)),
    A0=A0,
    Ms=MS,
    K0=K0,
    alpha=4000.0,
    gamma=0.0,
    cuda=False,
    cvode=False,
    useavgn=1,
    # This script reads the applied field back from result[4], so the individual H
    # components have to be returned.
    usereturnhall=True,
)

# Set the uniaxial easy axis to z for every cell.
problem.u_ea[:, :] = np.array([0.0, 0.0, 1.0])

# These solver settings match the small single-grain examples in this project.
problem.t = np.linspace(0.0, 1.0e-9, 2)
problem.nt = len(problem.t)
problem.t_conv = problem.t.copy()
problem.nt_conv = len(problem.t_conv)
problem.exch_presize = problem.ntot * 28
problem.use_fmm = 0
problem.ifunif = 1
problem.nlmin = 1
problem.nlmax = 5
problem.log_dir = str(Path(tempfile.gettempdir()) / "magtense_adaptive_hysteresis")
Path(problem.log_dir).mkdir(parents=True, exist_ok=True)
problem.window_enabled = 0
problem.trace_enabled = 0

# MagTense expects H and dH in A/m.  Dividing mu0*H in tesla by mu0 performs
# that conversion.  dM below is the root-mean-square change of the normalised
# magnetisation components between consecutive field states.
result = problem.run_hysteresis_adaptive(
    # Cartesian field vector at the first, positively saturated state.
    H_start=FIELD_START_T / MU0 * field_direction,
    # Cartesian field vector at which the descending sweep must finish.
    H_end=FIELD_END_T / MU0 * field_direction,
    # Step proposed immediately after solving the initial field state.
    dH_initial=INITIAL_STEP_T / MU0,
    # Smallest permitted step; large changes at this step must be accepted.
    dH_min=MINIMUM_STEP_T / MU0,
    # Largest permitted step, used in slowly varying parts of the curve.
    dH_max=MAXIMUM_STEP_T / MU0,
    # Maximum number of accepted states, including the initial field state.
    max_steps=MAXIMUM_ACCEPTED_STEPS,
    # Grow the next step when dM is below this small-change threshold.
    dM_min=1.0e-3,
    # Shrink the next step when dM exceeds this target threshold.
    dM_target=1.0e-2,
    # Reject and retry a step when dM exceeds this threshold and dH > dH_min.
    dM_reject=5.0e-2,
    # Multiply dH by this factor after a sufficiently small change.
    dH_grow=1.5,
    # Multiply dH by this factor when shrinking or retrying a step.
    dH_shrink=0.75,
    # Reject a sign-changing step larger than this value and retry more finely.
    switch_refine_dH=MINIMUM_STEP_T / MU0,
)

# Fixed-step alternative
# ----------------------
# The corresponding conventional fixed-step solver call is shown below for
# reference.  H_ext has one row per field state.
# Column zero is a scalar field label in tesla; columns one to three contain
# the Cartesian H vector in A/m used by the solver.
#
# FIXED_STEP_T = 0.01
# n_intervals = int(round((FIELD_START_T - FIELD_END_T) / FIXED_STEP_T))
# fixed_field_t = np.linspace(FIELD_START_T, FIELD_END_T, n_intervals + 1)
# H_ext = np.zeros((fixed_field_t.size, 4))
# H_ext[:, 0] = fixed_field_t
# H_ext[:, 1:4] = (
#     fixed_field_t[:, np.newaxis] / MU0 * field_direction[np.newaxis, :]
# )
# result = problem.run_hysteresis(
#     H_ext=H_ext,  # Complete fixed field schedule; no adaptive options apply.
# )
#
# A fixed sweep already knows its number and values of field states:
# n_accepted = fixed_field_t.size
# field_t = fixed_field_t
# The adaptive field extraction immediately below is therefore not required.

# The final result item gives the number of accepted field states.  Convert the
# accepted Cartesian fields back to their scalar mu0*H values in tesla.
n_accepted = int(result[-1])
field_vectors = np.asarray(result[4][0, 0, :n_accepted, :])
field_t = MU0 * (field_vectors @ field_direction)

# result[1] has axes (time, cell, field, component).  Select the final time,
# average the normalised z component over all cells, and multiply by mu0*Ms.
magnetisation = np.asarray(result[1][1, :, :n_accepted, :])
mean_mz_t = BS * np.mean(magnetisation[:, :, 2], axis=0)

print(f"Adaptive sweep completed with {n_accepted} accepted field states.")
