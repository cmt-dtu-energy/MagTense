#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Single-grain coercivity queue submitter
# ============================================================
#
# This script does not run simulations itself. It expands the configured sweep
# into generic scheduler jobs, gives every submission an isolated output folder,
# and adds the jobs to the persistent queue.
#
# IMPORTANT: Keep --dry-run whenever this script is extended. A dry run must do
# all environment/parameter validation and create the exact job JSON files that
# would be submitted, but it must NEVER add jobs to the scheduler queue.
# ============================================================


# ============================================================
# User settings
# ============================================================
#
# These defaults define a normal submission. Command-line flags override them,
# so routine material changes do not require editing this file.

# Runtime used by the queued simulations. The scheduler service does not inherit
# an activated terminal environment, so every job records this Conda runtime.
PYTHON_BIN="${PYTHON_BIN:-python}"
CONDA_ENV="${CONDA_ENV:-magtense-env}"
CONDA_EXE="${CONDA_EXE:-}"

# OpenMP affinity matching the `omp_pcores_1t` preset in ~/.bashrc. These values
# select one hardware thread from each of the eight high-performance P-cores.
# They are written into every job because systemd does not source ~/.bashrc.
OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
OMP_PLACES="${OMP_PLACES:-}"
if [[ -z "$OMP_PLACES" ]]; then
  OMP_PLACES='{0},{2},{4},{6},{8},{10},{12},{14}'
fi
OMP_PROC_BIND="${OMP_PROC_BIND:-close}"

# Material properties. MU0_MS_T is entered as mu0*Ms in tesla; the Python
# experiment converts it to Ms in A/m before constructing MicromagProblem.
MU0_MS_T="2.4"
A0="7e-12"
K0="1e6"
FIELD_MIN_T="-3.0"
FIELD_MAX_T="1.0"

# Run modes applied to every point in the sweep.
USE_ADAPTIVE=true
USE_PERIODIC=true

# Cartesian sweep values. With the defaults this creates
# 3 dh values x 3 resolutions x 8 size factors = 72 jobs.
DH_MIN_VALUES="0.001,0.01,0.1"
RESOLUTIONS="10,12,14"
SIZE_FACTORS="1,5,10,15,20,25,30,40"
DRY_RUN=false


# ============================================================
# Repository and executable paths
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
PY_SCRIPT="${SCRIPT_DIR}/single_grain_coercivity.py"
SCHEDULER="${REPO_ROOT}/tools/job_scheduler/job_scheduler.py"
PYTHON_BIN="$(command -v "$PYTHON_BIN")"
BATCH_ROOT="${BATCH_ROOT:-$SCRIPT_DIR}"


# ============================================================
# Command-line interface
# ============================================================

usage() {
  cat <<'EOF'
Usage: run_script.sh [options]

Material options:
  --mu0-ms-t VALUE       Saturation magnetization mu0*Ms [T]
  --a0 VALUE             Exchange stiffness [J/m]
  --k0 VALUE             Uniaxial anisotropy [J/m^3]
  --field-min-t VALUE    Minimum applied field [T] (default: -3.0)
  --field-max-t VALUE    Maximum applied field [T] (default: 1.0)

Runtime options:
  --conda-env NAME       Conda environment used by every job (default: magtense-env)

Sweep options (comma-separated):
  --resolutions VALUES   Grid resolutions
  --size-factors VALUES  Grain size factors
  --dh-min-values VALUES Adaptive minimum field steps [T]

Modes:
  --adaptive / --no-adaptive
  --periodic / --no-periodic
  --dry-run              Create the batch manifest and print jobs without queueing
  -h, --help

Example:
  ./run_script.sh --mu0-ms-t 1.57 --a0 7e-12 --k0 1e6
EOF
}

while (( $# > 0 )); do
  case "$1" in
    --mu0-ms-t) MU0_MS_T="$2"; shift 2 ;;
    --a0) A0="$2"; shift 2 ;;
    --k0) K0="$2"; shift 2 ;;
    --field-min-t) FIELD_MIN_T="$2"; shift 2 ;;
    --field-max-t) FIELD_MAX_T="$2"; shift 2 ;;
    --conda-env) CONDA_ENV="$2"; shift 2 ;;
    --resolutions) RESOLUTIONS="$2"; shift 2 ;;
    --size-factors) SIZE_FACTORS="$2"; shift 2 ;;
    --dh-min-values) DH_MIN_VALUES="$2"; shift 2 ;;
    --adaptive) USE_ADAPTIVE=true; shift ;;
    --no-adaptive) USE_ADAPTIVE=false; shift ;;
    --periodic) USE_PERIODIC=true; shift ;;
    --no-periodic) USE_PERIODIC=false; shift ;;
    --dry-run) DRY_RUN=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done


# ============================================================
# Validate scheduler and Conda runtime
# ============================================================
#
# All validation happens before batch directories or queue entries are created.
# This prevents a misspelled environment from leaving a half-submitted sweep.

if [[ ! -f "$SCHEDULER" ]]; then
  echo "Scheduler not found: $SCHEDULER" >&2
  exit 1
fi

if [[ -z "$CONDA_EXE" ]]; then
  CONDA_EXE="$(command -v conda || true)"
fi
if [[ -z "$CONDA_EXE" || ! -x "$CONDA_EXE" ]]; then
  echo "Conda executable not found. Set CONDA_EXE or add conda to PATH." >&2
  exit 1
fi
CONDA_EXE="$(realpath "$CONDA_EXE")"

# Resolve the environment name to one absolute prefix. Recording the prefix in
# each job makes execution independent of Conda's currently active environment.
conda_envs_json="$("$CONDA_EXE" env list --json)"
JOB_CONDA_PREFIX="$("$PYTHON_BIN" - "$CONDA_ENV" "$conda_envs_json" <<'PY'
import json
import sys
from pathlib import Path

name = sys.argv[1]
environments = json.loads(sys.argv[2]).get("envs", [])
matches = [path for path in environments if Path(path).name == name]
if len(matches) != 1:
    available = ", ".join(sorted(Path(path).name for path in environments))
    raise SystemExit(
        f"Conda environment {name!r} was not found uniquely. Available: {available}"
    )
print(Path(matches[0]).resolve())
PY
)"

SOURCE_DIR="$(realpath "${SCRIPT_DIR}/../../src")"

# Verify the actual runtime once per batch rather than discovering a missing
# dependency hours later when the first scheduled simulation starts.
if ! PYTHONPATH="${SOURCE_DIR}${PYTHONPATH:+:${PYTHONPATH}}" \
  "$CONDA_EXE" run --no-capture-output --prefix "$JOB_CONDA_PREFIX" \
  python -c 'import numpy; import magtense'; then
  echo "Conda environment '$CONDA_ENV' cannot import numpy and magtense." >&2
  exit 1
fi

# Reject empty affinity settings before creating a batch. Advanced users may
# override these variables when submitting on another computer.
if [[ -z "$OMP_NUM_THREADS" || -z "$OMP_PLACES" || -z "$OMP_PROC_BIND" ]]; then
  echo "OMP_NUM_THREADS, OMP_PLACES, and OMP_PROC_BIND must be non-empty." >&2
  exit 1
fi


# ============================================================
# Parse and validate sweep values
# ============================================================

IFS=',' read -r -a resolution_values <<< "$RESOLUTIONS"
IFS=',' read -r -a size_factor_values <<< "$SIZE_FACTORS"
IFS=',' read -r -a dh_min_values <<< "$DH_MIN_VALUES"

"$PYTHON_BIN" - "$MU0_MS_T" "$A0" "$K0" "$FIELD_MIN_T" "$FIELD_MAX_T" \
  "$RESOLUTIONS" "$SIZE_FACTORS" "$DH_MIN_VALUES" <<'PY'
import math
import sys

mu0_ms_t, a0, k0, field_min_t, field_max_t = map(float, sys.argv[1:6])
if not all(map(math.isfinite, (mu0_ms_t, a0, k0, field_min_t, field_max_t))):
    raise SystemExit("Material values must be finite")
if mu0_ms_t <= 0 or a0 <= 0 or k0 < 0:
    raise SystemExit("mu0*Ms and A0 must be positive; K0 must be non-negative")
if field_min_t >= field_max_t:
    raise SystemExit("field-min-t must be less than field-max-t")
for index, (name, value) in enumerate(zip(
    ("resolutions", "size factors", "dh minimum values"), sys.argv[6:9]
)):
    try:
        parsed = [float(item) for item in value.split(",") if item]
    except ValueError as error:
        raise SystemExit(f"Invalid {name}: {value}") from error
    if not parsed:
        raise SystemExit(f"No {name} were provided")
    if not all(math.isfinite(item) and item > 0 for item in parsed):
        raise SystemExit(f"All {name} must be finite and positive")
    if index == 0 and not all(item.is_integer() for item in parsed):
        raise SystemExit("All resolutions must be integers")
PY


# ============================================================
# Create one isolated batch directory
# ============================================================
#
# Jobs write directly to their final destinations. Nothing needs to be copied
# after execution, and repeated material sweeps cannot overwrite older batches.

timestamp="$(date +%Y%m%d_%H%M%S)"
material_dir="res_Ms_${MU0_MS_T}_A0_${A0}_K0_${K0}"
batch_dir="${BATCH_ROOT}/${material_dir}/batch_${timestamp}_$$"
batch_id="$(realpath -m "$batch_dir")"
batch_label="$(basename "$batch_dir")"
results_dir="${batch_dir}/results"
logs_dir="${batch_dir}/logs"
timer_logs_dir="${batch_dir}/timer_logs_single_grain"
status_dir="${batch_dir}/status"
jobs_dir="${batch_dir}/jobs"
mkdir -p "$results_dir" "$logs_dir" "$timer_logs_dir" "$status_dir" "$jobs_dir"

manifest="${batch_dir}/manifest.json"
summary="${batch_dir}/summary.tsv"

# The manifest describes the complete submission, including the exact Conda
# runtime used by every generated job.
"$PYTHON_BIN" - "$manifest" "$timestamp" "$MU0_MS_T" "$A0" "$K0" \
  "$FIELD_MIN_T" "$FIELD_MAX_T" "$USE_ADAPTIVE" "$USE_PERIODIC" "$RESOLUTIONS" "$SIZE_FACTORS" \
  "$DH_MIN_VALUES" "$DRY_RUN" "$CONDA_ENV" "$CONDA_EXE" \
  "$JOB_CONDA_PREFIX" "$OMP_NUM_THREADS" "$OMP_PLACES" \
  "$OMP_PROC_BIND" "$batch_id" "$batch_label" <<'PY'
import json
import sys
from pathlib import Path

(
    manifest,
    timestamp,
    mu0_ms_t,
    a0,
    k0,
    field_min_t,
    field_max_t,
    adaptive,
    periodic,
    resolutions,
    size_factors,
    dh_min_values,
    dry_run,
    conda_name,
    conda_executable,
    conda_prefix,
    omp_num_threads,
    omp_places,
    omp_proc_bind,
    batch_id,
    batch_label,
) = sys.argv[1:]
data = {
    "schema_version": 1,
    "created_at": timestamp,
    "material": {
        "mu0_Ms_T": float(mu0_ms_t),
        "A0_J_per_m": float(a0),
        "K0_J_per_m3": float(k0),
    },
    "batch": {"id": batch_id, "label": batch_label},
    "field_min_t": float(field_min_t),
    "field_max_t": float(field_max_t),
    "adaptive": adaptive == "true",
    "periodic": periodic == "true",
    "resolutions": [int(value) for value in resolutions.split(",")],
    "size_factors": [float(value) for value in size_factors.split(",")],
    "adaptive_dh_min_t": [float(value) for value in dh_min_values.split(",")],
    "dry_run": dry_run == "true",
    "runtime": {
        "type": "conda",
        "name": conda_name,
        "executable": conda_executable,
        "prefix": conda_prefix,
    },
    "environment": {
        "OMP_NUM_THREADS": omp_num_threads,
        "OMP_PLACES": omp_places,
        "OMP_PROC_BIND": omp_proc_bind,
    },
}
Path(manifest).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
PY

printf 'job_id\tstatus\tadaptive\tperiodic\tdh_min\tn\tsize_factor\tntot\tstdout\tstderr\treceipt\n' > "$summary"


# ============================================================
# Stable output naming
# ============================================================

make_stem() {
  "$PYTHON_BIN" - "$1" "$2" "$USE_PERIODIC" "$USE_ADAPTIVE" "$3" <<'PY'
import sys

size_factor = float(sys.argv[1])
n = int(sys.argv[2])
periodic = sys.argv[3] == "true"
adaptive = sys.argv[4] == "true"
dh_min = sys.argv[5]
stem = f"single_grain_sf{size_factor:.2g}_n{n}"
if periodic:
    stem += "_P"
if adaptive:
    stem += f"_A_FS{abs(float(dh_min)):.1e}"
print(stem)
PY
}


# ============================================================
# Convert one sweep point into one generic scheduler job
# ============================================================

submit_one() {
  local dh_min="$1"
  local n="$2"
  local size_factor="$3"
  local stem job_file stdout_path stderr_path receipt_path job_id status
  stem="$(make_stem "$size_factor" "$n" "$dh_min")"
  job_file="${jobs_dir}/${stem}.json"
  stdout_path="${logs_dir}/${stem}.out"
  stderr_path="${logs_dir}/${stem}.err"
  receipt_path="${status_dir}/${stem}.json"

  # Build JSON with Python instead of shell string concatenation. This preserves
  # argv boundaries and prevents paths or values from being interpreted by a
  # shell when the executor eventually launches the simulation.
  "$PYTHON_BIN" - "$job_file" "$PY_SCRIPT" "$SCRIPT_DIR" \
    "$results_dir" "$timer_logs_dir" "$stdout_path" "$stderr_path" \
    "$receipt_path" "$stem" "$MU0_MS_T" "$A0" "$K0" "$n" \
    "$FIELD_MIN_T" "$FIELD_MAX_T" "$size_factor" "$USE_ADAPTIVE" "$USE_PERIODIC" "$dh_min" \
    "$SOURCE_DIR" "$CONDA_ENV" "$CONDA_EXE" "$JOB_CONDA_PREFIX" \
    "$OMP_NUM_THREADS" "$OMP_PLACES" "$OMP_PROC_BIND" "$batch_id" "$batch_label" <<'PY'
import json
import os
import sys
from pathlib import Path

(
    job_file,
    script,
    cwd,
    results_dir,
    timer_logs_dir,
    stdout_path,
    stderr_path,
    receipt_path,
    stem,
    mu0_ms_t,
    a0,
    k0,
    n,
    field_min_t,
    field_max_t,
    size_factor,
    adaptive,
    periodic,
    dh_min,
    source_dir,
    conda_name,
    conda_executable,
    conda_prefix,
    omp_num_threads,
    omp_places,
    omp_proc_bind,
    batch_id,
    batch_label,
) = sys.argv[1:]
argv = [
    "python",
    script,
    "--mu0-ms-t", mu0_ms_t,
    "--a0", a0,
    "--k0", k0,
    "--n", n,
    "--size-factor", size_factor,
    "--field-min-t", field_min_t,
    "--field-max-t", field_max_t,
    "--output-dir", results_dir,
    "--timer-log-dir", timer_logs_dir,
]
if adaptive == "true":
    argv.extend(["--adaptive", "--adaptive-dh-min-t", dh_min])
if periodic == "true":
    argv.append("--periodic")
existing_pythonpath = os.environ.get("PYTHONPATH", "")
pythonpath = source_dir if not existing_pythonpath else f"{source_dir}:{existing_pythonpath}"
job_environment = {
    "PYTHONPATH": pythonpath,
    "OMP_NUM_THREADS": omp_num_threads,
    "OMP_PLACES": omp_places,
    "OMP_PROC_BIND": omp_proc_bind,
}
job = {
    "batch": {"id": batch_id, "label": batch_label},
    "argv": argv,
    "working_directory": cwd,
    "environment": job_environment,
    "stdout_path": stdout_path,
    "stderr_path": stderr_path,
    "receipt_path": receipt_path,
    "runtime": {
        "type": "conda",
        "name": conda_name,
        "executable": conda_executable,
        "prefix": conda_prefix,
    },
    "metadata": {
        "label": stem,
        "experiment": "single_grain_coercivity",
        "mu0_Ms_T": float(mu0_ms_t),
        "A0": float(a0),
        "K0": float(k0),
        "field_min_t": float(field_min_t),
        "field_max_t": float(field_max_t),
        "n": int(n),
        "size_factor": float(size_factor),
        "adaptive": adaptive == "true",
        "periodic": periodic == "true",
        "adaptive_dh_min_t": float(dh_min) if adaptive == "true" else None,
        "conda_environment": conda_name,
        "conda_prefix": conda_prefix,
        "openmp": {
            "OMP_NUM_THREADS": omp_num_threads,
            "OMP_PLACES": omp_places,
            "OMP_PROC_BIND": omp_proc_bind,
        },
    },
}
Path(job_file).write_text(json.dumps(job, indent=2) + "\n", encoding="utf-8")
PY

  if [[ "$DRY_RUN" == true ]]; then
    # IMPORTANT: print the fully wrapped Conda command, but do not call submit.
    # This branch is the safety check users should run before a large sweep.
    job_id="-"
    status="dry-run"
    echo "DRY RUN: $stem"
    "$PYTHON_BIN" - "$job_file" <<'PY'
import json
import shlex
import sys

job = json.load(open(sys.argv[1], encoding="utf-8"))
runtime = job["runtime"]
command = [
    runtime["executable"],
    "run",
    "--no-capture-output",
    "--prefix",
    runtime["prefix"],
    *job["argv"],
]
print(f"  Conda environment: {runtime['name']} ({runtime['prefix']})")
print(
    "  OpenMP: "
    f"OMP_NUM_THREADS={job['environment']['OMP_NUM_THREADS']}, "
    f"OMP_PLACES={job['environment']['OMP_PLACES']}, "
    f"OMP_PROC_BIND={job['environment']['OMP_PROC_BIND']}"
)
print("  " + shlex.join(command))
PY
  else
    # The scheduler validates and copies the job into its persistent pending
    # directory. The returned ID is recorded in this batch's summary.
    job_id="$("$PYTHON_BIN" "$SCHEDULER" submit --job-file "$job_file")"
    status="pending"
    echo "Queued ${stem}: ${job_id}"
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$job_id" "$status" "$USE_ADAPTIVE" "$USE_PERIODIC" "$dh_min" \
    "$n" "$size_factor" "$((n*n*n))" "$stdout_path" "$stderr_path" \
    "$receipt_path" >> "$summary"
}


# ============================================================
# Expand the Cartesian parameter sweep
# ============================================================

if [[ "$USE_ADAPTIVE" == true ]]; then
  for dh_min in "${dh_min_values[@]}"; do
    for n in "${resolution_values[@]}"; do
      for size_factor in "${size_factor_values[@]}"; do
        submit_one "$dh_min" "$n" "$size_factor"
      done
    done
  done
else
  for n in "${resolution_values[@]}"; do
    for size_factor in "${size_factor_values[@]}"; do
      submit_one "none" "$n" "$size_factor"
    done
  done
fi

# Keep these paths as the final output so they are easy to find in terminal logs.
echo "Batch directory: $batch_dir"
echo "Manifest: $manifest"
echo "Summary: $summary"
