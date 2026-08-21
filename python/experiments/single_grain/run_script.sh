#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Single-grain coercivity queue submitter
# ============================================================
#
# One invocation may contain several named or custom materials. Each material
# receives its own scheduler batch and output directory. Jobs inside that batch
# sweep physical grain size, resolution, adaptive field step, and periodic mode.
#
# IMPORTANT: --dry-run must always build the exact manifests and job JSON files
# without adding anything to the scheduler queue.
# ============================================================


# ============================================================
# User settings
# ============================================================

PYTHON_BIN="${PYTHON_BIN:-python}"
CONDA_ENV="${CONDA_ENV:-magtense-env}"
CONDA_EXE="${CONDA_EXE:-}"

# Explicit service environment matching the local `omp_pcores_1t` preset.
OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
OMP_PLACES="${OMP_PLACES:-}"
if [[ -z "$OMP_PLACES" ]]; then
  OMP_PLACES='{0},{2},{4},{6},{8},{10},{12},{14}'
fi
OMP_PROC_BIND="${OMP_PROC_BIND:-close}"

# Named presets are expanded in this order, which is also their queue order.
MATERIAL_SELECTION="fe16n2,sm2fe17n3"
MATERIALS_EXPLICIT=false
CUSTOM_MATERIAL_SPECS=()

# The singular flags remain available for one custom material submission.
CUSTOM_NAME="custom"
CUSTOM_MU0_MS_T="2.4"
CUSTOM_A0="7e-12"
CUSTOM_K0="1e6"
CUSTOM_FIELD_MAX_T="1.0"
CUSTOM_FIELD_MIN_T="-3.0"
SINGULAR_CUSTOM=false

USE_ADAPTIVE=true
PERIODIC_MODES="both"
DH_MIN_VALUES="0.001,0.01,0.1"
RESOLUTIONS="10,12,14"
SIZES_NM="1,5,10,20,40,60,80,100,200"
SHAPES="cube"
APPEND_EXISTING_BATCHES=false
DRY_RUN=false

PROPERTY_GRID=false
PROPERTY_MS_RANGE_T="1.5:2.4:0.1"
PROPERTY_K0_RANGE="1e6:9e6:5e5"
PROPERTY_A0="7e-12"
PROPERTY_FIELD_MAX_T="2.0"
PLOT_RESULTS="auto"

RESOLUTIONS_EXPLICIT=false
DH_MIN_VALUES_EXPLICIT=false
SIZES_NM_EXPLICIT=false
SHAPES_EXPLICIT=false
PERIODIC_MODES_EXPLICIT=false
FIELD_MIN_EXPLICIT=false
FIELD_MAX_EXPLICIT=false


# ============================================================
# Repository paths
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

Named materials:
  --materials VALUES       Comma-separated presets: fe16n2,sm2fe17n3
                           (default: both, in that order)
  --custom-material SPEC   Repeatable NAME,MU0_MS_T,A0,K0,FIELD_MAX_T,FIELD_MIN_T

Property-grid mode:
  --property-grid          Sweep mu0*Ms and K0 instead of named materials
  --ms-range-t A:B:S       mu0*Ms range in T (default: 1.5:2.4:0.1)
  --k0-range A:B:S         K0 range in J/m^3 (default: 1e6:9e6:5e5)

Custom one-material mode (using any of these disables implicit named defaults):
  --material-name NAME
  --mu0-ms-t VALUE         Saturation magnetization mu0*Ms [T]
  --a0 VALUE               Exchange stiffness [J/m]
  --k0 VALUE               Uniaxial anisotropy [J/m^3]
  --field-min-t VALUE      Minimum applied field [T]
  --field-max-t VALUE      Maximum applied field [T]

Sweep options (comma-separated):
  --sizes-nm VALUES        Physical cubic side lengths [nm]
  --resolutions VALUES     Cubic grid resolutions
  --dh-min-values VALUES   Adaptive minimum field steps [T]
  --shapes VALUES          Comma-separated shapes/variants: cube,sphere,
                           ellipsoid_x,ellipsoid_z,rectangle,cylinder,hexagon
  --periodic-modes MODE    both, periodic, or non-periodic (default: both)

Modes and runtime:
  --adaptive / --no-adaptive
  --periodic               Shortcut for --periodic-modes periodic
  --no-periodic            Shortcut for --periodic-modes non-periodic
  --plot / --no-plot       Enable or disable per-job PNG plots
  --conda-env NAME         Conda environment (default: magtense-env)
  --append-existing-batches
                           Add jobs to the active 20260615 material batches
  --dry-run                Build and print jobs without queueing
  -h, --help

Examples:
  ./run_script.sh --dry-run
  ./run_script.sh --append-existing-batches --dry-run
  ./run_script.sh --materials sm2fe17n3 --periodic --dry-run
  ./run_script.sh --mu0-ms-t 1.57 --a0 7e-12 --k0 1e6 --dry-run
  ./run_script.sh --property-grid --dry-run
EOF
}

while (( $# > 0 )); do
  case "$1" in
    --property-grid) PROPERTY_GRID=true; shift ;;
    --ms-range-t) PROPERTY_MS_RANGE_T="$2"; shift 2 ;;
    --k0-range) PROPERTY_K0_RANGE="$2"; shift 2 ;;
    --materials) MATERIAL_SELECTION="$2"; MATERIALS_EXPLICIT=true; shift 2 ;;
    --custom-material) CUSTOM_MATERIAL_SPECS+=("$2"); shift 2 ;;
    --material-name) CUSTOM_NAME="$2"; SINGULAR_CUSTOM=true; shift 2 ;;
    --mu0-ms-t) CUSTOM_MU0_MS_T="$2"; SINGULAR_CUSTOM=true; shift 2 ;;
    --a0) CUSTOM_A0="$2"; SINGULAR_CUSTOM=true; shift 2 ;;
    --k0) CUSTOM_K0="$2"; SINGULAR_CUSTOM=true; shift 2 ;;
    --field-min-t) CUSTOM_FIELD_MIN_T="$2"; FIELD_MIN_EXPLICIT=true; SINGULAR_CUSTOM=true; shift 2 ;;
    --field-max-t) CUSTOM_FIELD_MAX_T="$2"; PROPERTY_FIELD_MAX_T="$2"; FIELD_MAX_EXPLICIT=true; SINGULAR_CUSTOM=true; shift 2 ;;
    --sizes-nm) SIZES_NM="$2"; SIZES_NM_EXPLICIT=true; shift 2 ;;
    --resolutions) RESOLUTIONS="$2"; RESOLUTIONS_EXPLICIT=true; shift 2 ;;
    --dh-min-values) DH_MIN_VALUES="$2"; DH_MIN_VALUES_EXPLICIT=true; shift 2 ;;
    --shapes) SHAPES="$2"; SHAPES_EXPLICIT=true; shift 2 ;;
    --periodic-modes) PERIODIC_MODES="$2"; PERIODIC_MODES_EXPLICIT=true; shift 2 ;;
    --adaptive) USE_ADAPTIVE=true; shift ;;
    --no-adaptive) USE_ADAPTIVE=false; shift ;;
    --periodic) PERIODIC_MODES="periodic"; PERIODIC_MODES_EXPLICIT=true; shift ;;
    --no-periodic) PERIODIC_MODES="non-periodic"; PERIODIC_MODES_EXPLICIT=true; shift ;;
    --plot) PLOT_RESULTS=true; shift ;;
    --no-plot) PLOT_RESULTS=false; shift ;;
    --conda-env) CONDA_ENV="$2"; shift 2 ;;
    --append-existing-batches) APPEND_EXISTING_BATCHES=true; shift ;;
    --dry-run) DRY_RUN=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ "$PROPERTY_GRID" == true ]]; then
  [[ "$RESOLUTIONS_EXPLICIT" == false ]] && RESOLUTIONS="14"
  [[ "$DH_MIN_VALUES_EXPLICIT" == false ]] && DH_MIN_VALUES="0.001"
  [[ "$SIZES_NM_EXPLICIT" == false ]] && SIZES_NM="20,25,30,35,40,45,50"
  [[ "$SHAPES_EXPLICIT" == false ]] && SHAPES="cube,ellipsoid_z"
  [[ "$PERIODIC_MODES_EXPLICIT" == false ]] && PERIODIC_MODES="non-periodic"
  [[ "$PLOT_RESULTS" == "auto" ]] && PLOT_RESULTS=false
  if [[ "$FIELD_MIN_EXPLICIT" == true ]]; then
    echo "--field-min-t is computed from the SW estimate in --property-grid mode." >&2
    exit 2
  fi
  if [[ "$APPEND_EXISTING_BATCHES" == true ]]; then
    echo "--append-existing-batches cannot be combined with --property-grid." >&2
    exit 2
  fi
else
  [[ "$PLOT_RESULTS" == "auto" ]] && PLOT_RESULTS=true
fi


# ============================================================
# Resolve material records
# ============================================================
#
# Records use a pipe delimiter internally:
# name|mu0_Ms_T|A0_J_per_m|K0_J_per_m3|field_max_T|field_min_T

material_records=()

if [[ "$PROPERTY_GRID" == false ]]; then

add_named_material() {
  case "${1,,}" in
    fe16n2|alpha-fe16n2|alpha_fe16n2)
      material_records+=("fe16n2|2.4|7e-12|1e6|1.0|-3.0")
      ;;
    sm2fe17n3)
      material_records+=("sm2fe17n3|1.54|7e-12|8.9e6|2.0|-15.0")
      ;;
    *)
      echo "Unknown material preset: $1 (expected fe16n2 or sm2fe17n3)" >&2
      exit 2
      ;;
  esac
}

# A singular custom flag means "one custom material" unless --materials was
# explicitly supplied. Repeatable --custom-material entries can be combined
# with explicitly selected presets, or used alone when --materials is omitted.
if [[ "$SINGULAR_CUSTOM" == true ]]; then
  if [[ "$MATERIALS_EXPLICIT" == true || ${#CUSTOM_MATERIAL_SPECS[@]} -gt 0 ]]; then
    echo "Singular material flags cannot be combined with --materials or --custom-material." >&2
    exit 2
  fi
  material_records+=(
    "${CUSTOM_NAME}|${CUSTOM_MU0_MS_T}|${CUSTOM_A0}|${CUSTOM_K0}|${CUSTOM_FIELD_MAX_T}|${CUSTOM_FIELD_MIN_T}"
  )
else
  if [[ ${#CUSTOM_MATERIAL_SPECS[@]} -eq 0 || "$MATERIALS_EXPLICIT" == true ]]; then
    IFS=',' read -r -a selected_materials <<< "$MATERIAL_SELECTION"
    for material in "${selected_materials[@]}"; do
      [[ -n "$material" ]] && add_named_material "$material"
    done
  fi
  for spec in "${CUSTOM_MATERIAL_SPECS[@]}"; do
    IFS=',' read -r name mu0_ms_t a0 k0 field_max_t field_min_t extra <<< "$spec"
    if [[ -n "${extra:-}" || -z "${field_min_t:-}" ]]; then
      echo "Invalid --custom-material: $spec" >&2
      echo "Expected NAME,MU0_MS_T,A0,K0,FIELD_MAX_T,FIELD_MIN_T" >&2
      exit 2
    fi
    material_records+=("${name}|${mu0_ms_t}|${a0}|${k0}|${field_max_t}|${field_min_t}")
  done
fi

if [[ ${#material_records[@]} -eq 0 ]]; then
  echo "At least one material must be selected." >&2
  exit 2
fi

if [[ "$APPEND_EXISTING_BATCHES" == true ]]; then
  for record in "${material_records[@]}"; do
    IFS='|' read -r preset _ <<< "$record"
    case "$preset" in
      fe16n2|sm2fe17n3) ;;
      *)
        echo "--append-existing-batches only supports fe16n2 and sm2fe17n3." >&2
        exit 2
        ;;
    esac
  done
fi

fi

# ============================================================
# Validate runtime and all sweep inputs before creating batches
# ============================================================

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
if ! PYTHONPATH="${SOURCE_DIR}${PYTHONPATH:+:${PYTHONPATH}}" \
  "$CONDA_EXE" run --no-capture-output --prefix "$JOB_CONDA_PREFIX" \
  python -c 'import numpy; import magtense'; then
  echo "Conda environment '$CONDA_ENV' cannot import numpy and magtense." >&2
  exit 1
fi

if [[ -z "$OMP_NUM_THREADS" || -z "$OMP_PLACES" || -z "$OMP_PROC_BIND" ]]; then
  echo "OMP_NUM_THREADS, OMP_PLACES, and OMP_PROC_BIND must be non-empty." >&2
  exit 1
fi

case "$PERIODIC_MODES" in
  both) periodic_values=(false true) ;;
  periodic) periodic_values=(true) ;;
  non-periodic|nonperiodic) periodic_values=(false) ;;
  *) echo "--periodic-modes must be both, periodic, or non-periodic" >&2; exit 2 ;;
esac

IFS=',' read -r -a resolution_values <<< "$RESOLUTIONS"
IFS=',' read -r -a size_nm_values <<< "$SIZES_NM"
IFS=',' read -r -a dh_min_values <<< "$DH_MIN_VALUES"
IFS=',' read -r -a shape_values <<< "$SHAPES"

if [[ "$PERIODIC_MODES" != "non-periodic" && "$PERIODIC_MODES" != "nonperiodic" ]]; then
  for shape in "${shape_values[@]}"; do
    normalized_shape="${shape,,}"
    normalized_shape="${normalized_shape//-/_}"
    if [[ "$normalized_shape" != "cube" ]]; then
      echo "Periodic modes are only supported for cube shape runs." >&2
      exit 2
    fi
  done
fi

"$PYTHON_BIN" - "$RESOLUTIONS" "$SIZES_NM" "$DH_MIN_VALUES" "$SHAPES" \
  "${material_records[@]}" <<'PY'
import math
import re
import sys

resolutions, sizes_nm, dh_values, shapes, *materials = sys.argv[1:]
for index, (name, raw) in enumerate((
    ("resolutions", resolutions),
    ("sizes in nm", sizes_nm),
    ("dh minimum values", dh_values),
)):
    try:
        values = [float(item) for item in raw.split(",") if item]
    except ValueError as error:
        raise SystemExit(f"Invalid {name}: {raw}") from error
    if not values or not all(math.isfinite(value) and value > 0 for value in values):
        raise SystemExit(f"All {name} must be finite and positive")
    if index == 0 and not all(value.is_integer() for value in values):
        raise SystemExit("All resolutions must be integers")

valid_shapes = {
    "cube", "rectangle", "recantgle", "sphere", "cylinder", "ellipsoid",
    "elipsoid", "ellipsoid_x", "ellipsoid_z", "elipsoid_x", "elipsoid_z",
    "hexagon",
}
shape_values = [shape.strip().lower().replace("-", "_") for shape in shapes.split(",") if shape]
if not shape_values:
    raise SystemExit("At least one shape must be selected")
unknown_shapes = sorted(set(shape_values) - valid_shapes)
if unknown_shapes:
    raise SystemExit(f"Unknown shapes: {', '.join(unknown_shapes)}")

for record in materials:
    fields = record.split("|")
    if len(fields) != 6:
        raise SystemExit(f"Invalid material record: {record}")
    name, *numbers = fields
    if not re.fullmatch(r"[A-Za-z0-9_-]+", name):
        raise SystemExit(f"Material name must use letters, numbers, _ or -: {name}")
    mu0_ms_t, a0, k0, field_max_t, field_min_t = map(float, numbers)
    if not all(map(math.isfinite, (mu0_ms_t, a0, k0, field_max_t, field_min_t))):
        raise SystemExit(f"Material values must be finite: {name}")
    if mu0_ms_t <= 0 or a0 <= 0 or k0 < 0:
        raise SystemExit(f"mu0*Ms and A0 must be positive; K0 non-negative: {name}")
    if field_min_t >= field_max_t:
        raise SystemExit(f"field minimum must be below field maximum: {name}")
PY


# ============================================================
# Create and populate one isolated material batch
# ============================================================

submission_timestamp="$(date +%Y%m%d_%H%M%S)"
batch_index=0

create_material_batch() {
  local record="$1"
  local raw_shape="$2"
  local preset mu0_ms_t a0 k0 field_max_t field_min_t
  IFS='|' read -r preset mu0_ms_t a0 k0 field_max_t field_min_t <<< "$record"
  batch_index=$((batch_index + 1))

  local material_dir batch_dir batch_id batch_label results_dir logs_dir
  local timer_logs_dir status_dir jobs_dir manifest summary summary_mode
  local dry_run_dir
  local shape_token shape_arg shape_variant_arg shape_dir
  material_dir="res_${preset}_Ms_${mu0_ms_t}_A0_${a0}_K0_${k0}"
  shape_token="${raw_shape,,}"
  shape_token="${shape_token//-/_}"
  case "$shape_token" in
    recantgle) shape_token="rectangle" ;;
    elipsoid) shape_token="ellipsoid_z" ;;
    ellipsoid|elipsoid_z) shape_token="ellipsoid_z" ;;
    elipsoid_x) shape_token="ellipsoid_x" ;;
  esac
  shape_arg="$shape_token"
  shape_variant_arg=""
  case "$shape_token" in
    ellipsoid_x) shape_arg="ellipsoid"; shape_variant_arg="x" ;;
    ellipsoid_z) shape_arg="ellipsoid"; shape_variant_arg="z" ;;
  esac
  shape_dir="shape_${shape_token}"

  if [[ "$APPEND_EXISTING_BATCHES" == true ]]; then
    case "$preset" in
      fe16n2)
        batch_dir="${SCRIPT_DIR}/res_fe16n2_Ms_2.4_A0_7e-12_K0_1e6/batch_20260615_114706_1_91756"
        ;;
      sm2fe17n3)
        batch_dir="${SCRIPT_DIR}/res_sm2fe17n3_Ms_1.54_A0_7e-12_K0_8.9e6/batch_20260615_114706_2_91756"
        ;;
    esac
    if [[ ! -d "$batch_dir" ]]; then
      echo "Append target batch not found: $batch_dir" >&2
      exit 1
    fi
  else
    batch_dir="${BATCH_ROOT}/${material_dir}/${shape_dir}/batch_${submission_timestamp}_${batch_index}_$$"
  fi

  batch_id="$(realpath -m "$batch_dir")"
  batch_label="$(basename "$batch_dir")"
  results_dir="${batch_dir}/results"
  logs_dir="${batch_dir}/logs"
  timer_logs_dir="${batch_dir}/timer_logs_single_grain"
  status_dir="${batch_dir}/status"
  jobs_dir="${batch_dir}/jobs"
  manifest="${batch_dir}/manifest.json"
  summary="${batch_dir}/summary.tsv"
  summary_mode="create"

  if [[ "$APPEND_EXISTING_BATCHES" == true ]]; then
    for required_dir in "$results_dir" "$logs_dir" "$timer_logs_dir" "$status_dir" "$jobs_dir"; do
      if [[ ! -d "$required_dir" ]]; then
        echo "Append target is missing required directory: $required_dir" >&2
        exit 1
      fi
    done
    if [[ ! -f "$manifest" || ! -f "$summary" ]]; then
      echo "Append target must contain manifest.json and summary.tsv: $batch_dir" >&2
      exit 1
    fi
    summary_mode="append"
    if [[ "$DRY_RUN" == true ]]; then
      dry_run_dir="${TMPDIR:-/tmp}/single_grain_append_dryrun/${material_dir}/${batch_label}"
      mkdir -p "${dry_run_dir}/jobs"
      jobs_dir="${dry_run_dir}/jobs"
      manifest="${dry_run_dir}/manifest.json"
      summary="${dry_run_dir}/summary.tsv"
      summary_mode="create"
    fi
  else
    mkdir -p "$results_dir" "$logs_dir" "$timer_logs_dir" "$status_dir" "$jobs_dir"
  fi

  if [[ "$summary_mode" == "create" ]]; then
    "$PYTHON_BIN" - "$manifest" "$submission_timestamp" "$preset" \
      "$mu0_ms_t" "$a0" "$k0" "$field_min_t" "$field_max_t" \
      "$USE_ADAPTIVE" "$PERIODIC_MODES" "$RESOLUTIONS" "$SIZES_NM" \
      "$DH_MIN_VALUES" "$SHAPES" "$shape_token" "$DRY_RUN" "$CONDA_ENV" "$CONDA_EXE" \
      "$JOB_CONDA_PREFIX" "$OMP_NUM_THREADS" "$OMP_PLACES" \
      "$OMP_PROC_BIND" "$batch_id" "$batch_label" <<'PY'
import json
import sys
from pathlib import Path

(
    manifest, timestamp, preset, mu0_ms_t, a0, k0, field_min_t, field_max_t,
    adaptive, periodic_modes, resolutions, sizes_nm, dh_min_values,
    shapes, shape_token, dry_run, conda_name, conda_executable, conda_prefix,
    omp_num_threads, omp_places, omp_proc_bind, batch_id, batch_label,
) = sys.argv[1:]
data = {
    "schema_version": 2,
    "created_at": timestamp,
    "preset": preset,
    "material": {
        "mu0_Ms_T": float(mu0_ms_t),
        "A0_J_per_m": float(a0),
        "K0_J_per_m3": float(k0),
    },
    "batch": {"id": batch_id, "label": batch_label},
    "field_min_t": float(field_min_t),
    "field_max_t": float(field_max_t),
    "adaptive": adaptive == "true",
    "periodic_modes": (
        [False, True] if periodic_modes == "both"
        else [periodic_modes == "periodic"]
    ),
    "resolutions": [int(value) for value in resolutions.split(",")],
    "sizes_nm": [float(value) for value in sizes_nm.split(",")],
    "shapes": [value for value in shapes.split(",") if value],
    "shape": shape_token,
    "adaptive_dh_min_t": [float(value) for value in dh_min_values.split(",")],
    "dry_run": dry_run == "true",
    "runtime": {
        "type": "conda", "name": conda_name,
        "executable": conda_executable, "prefix": conda_prefix,
    },
    "environment": {
        "OMP_NUM_THREADS": omp_num_threads,
        "OMP_PLACES": omp_places,
        "OMP_PROC_BIND": omp_proc_bind,
    },
}
Path(manifest).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
PY

    printf 'job_id\tstatus\tpreset\tshape\tadaptive\tperiodic\tdh_min\tn\tsize_nm\tL_m\tsize_factor\tntot\tstdout\tstderr\treceipt\n' > "$summary"
  fi

  make_stem() {
    "$PYTHON_BIN" - "$1" "$2" "$3" "$4" "$5" "$6" <<'PY'
import sys

size_nm = float(sys.argv[1])
n = int(sys.argv[2])
periodic = sys.argv[3] == "true"
adaptive = sys.argv[4] == "true"
dh_min = sys.argv[5]
shape = sys.argv[6]
size_token = format(size_nm, ".12g")
if shape == "cube":
    stem = f"single_grain_size{size_token}nm_n{n}"
else:
    stem = f"single_grain_{shape}_size{size_token}nm_n{n}"
if periodic:
    stem += "_P"
if adaptive:
    stem += f"_A_FS{abs(float(dh_min)):.1e}"
print(stem)
PY
  }

  submit_one() {
    local dh_min="$1" n="$2" size_nm="$3" periodic="$4"
    local stem job_file stdout_path stderr_path receipt_path job_id status
    stem="$(make_stem "$size_nm" "$n" "$periodic" "$USE_ADAPTIVE" "$dh_min" "$shape_token")"
    job_file="${jobs_dir}/${stem}.json"
    stdout_path="${logs_dir}/${stem}.out"
    stderr_path="${logs_dir}/${stem}.err"
    receipt_path="${status_dir}/${stem}.json"

    if [[ "$APPEND_EXISTING_BATCHES" == true && "$DRY_RUN" == false && -e "$job_file" ]]; then
      echo "Refusing to overwrite existing append job: $job_file" >&2
      exit 1
    fi

    "$PYTHON_BIN" - "$job_file" "$PY_SCRIPT" "$SCRIPT_DIR" "$results_dir" \
      "$timer_logs_dir" "$stdout_path" "$stderr_path" "$receipt_path" \
      "$stem" "$preset" "$mu0_ms_t" "$a0" "$k0" "$field_min_t" \
      "$field_max_t" "$n" "$size_nm" "$shape_token" "$shape_arg" \
      "$shape_variant_arg" "$USE_ADAPTIVE" "$periodic" "$dh_min" "$SOURCE_DIR" "$CONDA_ENV" "$CONDA_EXE" \
      "$JOB_CONDA_PREFIX" "$OMP_NUM_THREADS" "$OMP_PLACES" \
      "$OMP_PROC_BIND" "$batch_id" "$batch_label" <<'PY'
import json
import math
import os
import sys
from pathlib import Path

(
    job_file, script, cwd, results_dir, timer_logs_dir, stdout_path,
    stderr_path, receipt_path, stem, preset, mu0_ms_t, a0, k0,
    field_min_t, field_max_t, n, size_nm, shape_token, shape_arg,
    shape_variant_arg, adaptive, periodic, dh_min,
    source_dir, conda_name, conda_executable, conda_prefix,
    omp_num_threads, omp_places, omp_proc_bind, batch_id, batch_label,
) = sys.argv[1:]

mu0 = 4.0 * math.pi * 1e-7
ms = float(mu0_ms_t) / mu0
L_m = float(size_nm) * 1e-9
characteristic_length = math.sqrt(float(a0) / (0.5 * mu0 * ms**2))
size_factor = L_m / characteristic_length

argv = [
    "python", script,
    "--mu0-ms-t", mu0_ms_t,
    "--a0", a0,
    "--k0", k0,
    "--n", n,
    "--L", format(L_m, ".17g"),
    "--field-min-t", field_min_t,
    "--field-max-t", field_max_t,
    "--output-stem", stem,
    "--output-dir", results_dir,
    "--timer-log-dir", timer_logs_dir,
    "--shape", shape_arg,
]
if shape_variant_arg:
    argv.extend(["--shape-variant", shape_variant_arg])
if adaptive == "true":
    argv.extend(["--adaptive", "--adaptive-dh-min-t", dh_min])
if periodic == "true":
    argv.append("--periodic")

existing_pythonpath = os.environ.get("PYTHONPATH", "")
pythonpath = source_dir if not existing_pythonpath else f"{source_dir}:{existing_pythonpath}"
openmp = {
    "OMP_NUM_THREADS": omp_num_threads,
    "OMP_PLACES": omp_places,
    "OMP_PROC_BIND": omp_proc_bind,
}
job = {
    "batch": {"id": batch_id, "label": batch_label},
    "argv": argv,
    "working_directory": cwd,
    "environment": {"PYTHONPATH": pythonpath, **openmp},
    "stdout_path": stdout_path,
    "stderr_path": stderr_path,
    "receipt_path": receipt_path,
    "runtime": {
        "type": "conda", "name": conda_name,
        "executable": conda_executable, "prefix": conda_prefix,
    },
    "metadata": {
        "label": stem,
        "experiment": "single_grain_coercivity",
        "preset": preset,
        "shape": shape_token,
        "mu0_Ms_T": float(mu0_ms_t),
        "A0": float(a0),
        "K0": float(k0),
        "field_min_t": float(field_min_t),
        "field_max_t": float(field_max_t),
        "n": int(n),
        "size_nm": float(size_nm),
        "L_m": L_m,
        "size_factor": size_factor,
        "adaptive": adaptive == "true",
        "periodic": periodic == "true",
        "adaptive_dh_min_t": float(dh_min) if adaptive == "true" else None,
        "conda_environment": conda_name,
        "conda_prefix": conda_prefix,
        "openmp": openmp,
    },
}
Path(job_file).write_text(json.dumps(job, indent=2) + "\n", encoding="utf-8")
PY
    dimensions="$("$PYTHON_BIN" - "$job_file" <<'PY'
import json
import sys
job = json.load(open(sys.argv[1], encoding="utf-8"))
metadata = job["metadata"]
print(f"{metadata['L_m']:.17g}\t{metadata['size_factor']:.17g}")
PY
)"
    local L_m size_factor
    IFS=$'\t' read -r L_m size_factor <<< "$dimensions"

    if [[ "$DRY_RUN" == true ]]; then
      job_id="-"
      status="dry-run"
      echo "DRY RUN: ${preset}/${stem}"
      "$PYTHON_BIN" - "$job_file" <<'PY'
import json
import shlex
import sys
job = json.load(open(sys.argv[1], encoding="utf-8"))
runtime = job["runtime"]
command = [
    runtime["executable"], "run", "--no-capture-output", "--prefix",
    runtime["prefix"], *job["argv"],
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
      job_id="$("$PYTHON_BIN" "$SCHEDULER" submit --job-file "$job_file")"
      status="pending"
      echo "Queued ${preset}/${stem}: ${job_id}"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$job_id" "$status" "$preset" "$shape_token" "$USE_ADAPTIVE" "$periodic" \
      "$dh_min" "$n" "$size_nm" "$L_m" "$size_factor" "$((n*n*n))" \
      "$stdout_path" "$stderr_path" "$receipt_path" >> "$summary"
  }

  if [[ "$USE_ADAPTIVE" == true ]]; then
    for dh_min in "${dh_min_values[@]}"; do
      for n in "${resolution_values[@]}"; do
        for size_nm in "${size_nm_values[@]}"; do
          for periodic in "${periodic_values[@]}"; do
            submit_one "$dh_min" "$n" "$size_nm" "$periodic"
          done
        done
      done
    done
  else
    for n in "${resolution_values[@]}"; do
      for size_nm in "${size_nm_values[@]}"; do
        for periodic in "${periodic_values[@]}"; do
          submit_one "none" "$n" "$size_nm" "$periodic"
        done
      done
    done
  fi

  if [[ "$APPEND_EXISTING_BATCHES" == true && "$DRY_RUN" == false ]]; then
    "$PYTHON_BIN" - "$manifest" "$SIZES_NM" <<'PY'
import json
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
new_sizes = [float(value) for value in sys.argv[2].split(",") if value]
data = json.loads(manifest.read_text(encoding="utf-8"))
sizes = {float(value) for value in data.get("sizes_nm", [])}
sizes.update(new_sizes)
data["sizes_nm"] = sorted(sizes)
manifest.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
PY
  fi

  echo "Batch directory: $batch_dir"
  echo "Manifest: $manifest"
  echo "Summary: $summary"
}

create_property_grid_batch() {
  local raw_shape="$1"
  batch_index=$((batch_index + 1))

  local sweep_dir batch_dir batch_id batch_label results_dir logs_dir
  local timer_logs_dir status_dir jobs_dir manifest summary
  local shape_token shape_arg shape_variant_arg shape_dir
  sweep_dir="res_ms_k0_sweep_A0_${PROPERTY_A0}"
  shape_token="${raw_shape,,}"
  shape_token="${shape_token//-/_}"
  case "$shape_token" in
    recantgle) shape_token="rectangle" ;;
    elipsoid) shape_token="ellipsoid_z" ;;
    ellipsoid|elipsoid_z) shape_token="ellipsoid_z" ;;
    elipsoid_x) shape_token="ellipsoid_x" ;;
  esac
  shape_arg="$shape_token"
  shape_variant_arg=""
  case "$shape_token" in
    ellipsoid_x) shape_arg="ellipsoid"; shape_variant_arg="x" ;;
    ellipsoid_z) shape_arg="ellipsoid"; shape_variant_arg="z" ;;
  esac
  shape_dir="shape_${shape_token}"

  batch_dir="${BATCH_ROOT}/${sweep_dir}/${shape_dir}/batch_${submission_timestamp}_${batch_index}_$$"
  batch_id="$(realpath -m "$batch_dir")"
  batch_label="$(basename "$batch_dir")"
  results_dir="${batch_dir}/results"
  logs_dir="${batch_dir}/logs"
  timer_logs_dir="${batch_dir}/timer_logs_single_grain"
  status_dir="${batch_dir}/status"
  jobs_dir="${batch_dir}/jobs"
  manifest="${batch_dir}/manifest.json"
  summary="${batch_dir}/summary.tsv"
  mkdir -p "$results_dir" "$logs_dir" "$timer_logs_dir" "$status_dir" "$jobs_dir"

  "$PYTHON_BIN" - "$manifest" "$summary" "$jobs_dir" "$PY_SCRIPT" "$SCRIPT_DIR" \
    "$results_dir" "$timer_logs_dir" "$logs_dir" "$status_dir" \
    "$SCHEDULER" "$SOURCE_DIR" "$CONDA_ENV" "$CONDA_EXE" "$JOB_CONDA_PREFIX" \
    "$OMP_NUM_THREADS" "$OMP_PLACES" "$OMP_PROC_BIND" "$batch_id" "$batch_label" \
    "$submission_timestamp" "$shape_token" "$shape_arg" "$shape_variant_arg" \
    "$PROPERTY_MS_RANGE_T" "$PROPERTY_K0_RANGE" "$PROPERTY_A0" "$PROPERTY_FIELD_MAX_T" \
    "$USE_ADAPTIVE" "$PERIODIC_MODES" "$RESOLUTIONS" "$SIZES_NM" "$DH_MIN_VALUES" \
    "$DRY_RUN" "$PLOT_RESULTS" <<'PY'
import json
import math
import os
import re
import shlex
import subprocess
import sys
from decimal import Decimal
from pathlib import Path

(
    manifest, summary, jobs_dir, script, cwd, results_dir, timer_logs_dir,
    logs_dir, status_dir, scheduler, source_dir, conda_name, conda_executable,
    conda_prefix, omp_num_threads, omp_places, omp_proc_bind, batch_id,
    batch_label, timestamp, shape_token, shape_arg, shape_variant_arg,
    ms_range_t, k0_range, a0, field_max_t, adaptive, periodic_modes,
    resolutions, sizes_nm, dh_min_values, dry_run, plot_results,
) = sys.argv[1:]

MU0 = 4.0 * math.pi * 1e-7


def parse_range(spec: str, name: str) -> list[Decimal]:
    parts = spec.split(":")
    if len(parts) != 3:
        raise SystemExit(f"{name} must use START:STOP:STEP syntax")
    start, stop, step = (Decimal(part) for part in parts)
    if step <= 0:
        raise SystemExit(f"{name} step must be positive")
    if stop < start:
        raise SystemExit(f"{name} stop must be >= start")
    values: list[Decimal] = []
    value = start
    epsilon = step / Decimal("1000000")
    while value <= stop + epsilon:
        values.append(value)
        value += step
    return values


def parse_positive_csv(raw: str, name: str, *, integer: bool = False) -> list[float]:
    try:
        values = [float(item) for item in raw.split(",") if item]
    except ValueError as error:
        raise SystemExit(f"Invalid {name}: {raw}") from error
    if not values or not all(math.isfinite(value) and value > 0 for value in values):
        raise SystemExit(f"All {name} must be finite and positive")
    if integer and not all(value.is_integer() for value in values):
        raise SystemExit(f"All {name} must be integers")
    return values


def token(value: float) -> str:
    text = format(value, ".12g")
    return text.replace("+", "").replace("-", "m").replace(".", "p")


def size_token(value: float) -> str:
    return format(value, ".12g")


def sw_coercivity_t(mu0_ms_t: float, k0: float) -> float:
    return 2.0 * k0 * MU0 / mu0_ms_t


def field_min_from_sw(sw_t: float) -> float:
    return -(math.ceil(sw_t * 2.0 - 1e-12) / 2.0 + 0.5)


ms_values = [float(value) for value in parse_range(ms_range_t, "ms range")]
k0_values = [float(value) for value in parse_range(k0_range, "K0 range")]
resolution_values = [int(value) for value in parse_positive_csv(resolutions, "resolutions", integer=True)]
size_values = parse_positive_csv(sizes_nm, "sizes in nm")
dh_values = parse_positive_csv(dh_min_values, "dh minimum values")

if periodic_modes not in {"non-periodic", "nonperiodic"}:
    raise SystemExit("--property-grid currently supports only non-periodic runs")
periodic_values = [False]

field_max = float(field_max_t)
if not math.isfinite(field_max):
    raise SystemExit("field max must be finite")
a0_value = float(a0)
if not math.isfinite(a0_value) or a0_value <= 0.0:
    raise SystemExit("A0 must be positive")

plot_enabled = plot_results == "true"
adaptive_enabled = adaptive == "true"
jobs_path = Path(jobs_dir)
logs_path = Path(logs_dir)
status_path = Path(status_dir)

manifest_data = {
    "schema_version": 3,
    "created_at": timestamp,
    "preset": "property_grid",
    "sweep": {
        "type": "ms_k0_grid",
        "mu0_Ms_T": ms_values,
        "K0_J_per_m3": k0_values,
        "sw_field_rule": "field_min_t = -(ceil(SW_T / 0.5) * 0.5 + 0.5)",
    },
    "material": {"A0_J_per_m": a0_value},
    "batch": {"id": batch_id, "label": batch_label},
    "field_max_t": field_max,
    "adaptive": adaptive_enabled,
    "periodic_modes": periodic_values,
    "resolutions": resolution_values,
    "sizes_nm": size_values,
    "shape": shape_token,
    "adaptive_dh_min_t": dh_values,
    "dry_run": dry_run == "true",
    "plot_results": plot_enabled,
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
Path(manifest).write_text(json.dumps(manifest_data, indent=2) + "\n", encoding="utf-8")

summary_columns = [
    "job_id", "status", "preset", "shape", "mu0_ms_t", "k0", "sw_coercivity_t",
    "field_min_t", "field_max_t", "adaptive", "periodic", "dh_min", "n",
    "size_nm", "L_m", "size_factor", "ntot", "stdout", "stderr", "receipt",
]

existing_pythonpath = os.environ.get("PYTHONPATH", "")
pythonpath = source_dir if not existing_pythonpath else f"{source_dir}:{existing_pythonpath}"
openmp = {
    "OMP_NUM_THREADS": omp_num_threads,
    "OMP_PLACES": omp_places,
    "OMP_PROC_BIND": omp_proc_bind,
}

rows: list[list[object]] = []
dry_run_examples: list[list[str]] = []
for mu0_ms_t in ms_values:
    ms_a_per_m = mu0_ms_t / MU0
    characteristic_length = math.sqrt(a0_value / (0.5 * MU0 * ms_a_per_m**2))
    for k0 in k0_values:
        sw_t = sw_coercivity_t(mu0_ms_t, k0)
        field_min = field_min_from_sw(sw_t)
        material_label = f"grid_Ms{token(mu0_ms_t)}T_K0{token(k0)}"
        for dh_min in dh_values:
            for n in resolution_values:
                for size_nm in size_values:
                    for periodic in periodic_values:
                        L_m = size_nm * 1e-9
                        size_factor = L_m / characteristic_length
                        base = (
                            f"single_grain_Ms{token(mu0_ms_t)}T_K0{token(k0)}_"
                            f"size{size_token(size_nm)}nm_n{n}"
                        )
                        stem = (
                            base if shape_token == "cube"
                            else base.replace("single_grain_", f"single_grain_{shape_token}_", 1)
                        )
                        if periodic:
                            stem += "_P"
                        if adaptive_enabled:
                            stem += f"_A_FS{abs(float(dh_min)):.1e}"

                        job_file = jobs_path / f"{stem}.json"
                        stdout_path = logs_path / f"{stem}.out"
                        stderr_path = logs_path / f"{stem}.err"
                        receipt_path = status_path / f"{stem}.json"
                        argv = [
                            "python", script,
                            "--mu0-ms-t", format(mu0_ms_t, ".17g"),
                            "--a0", a0,
                            "--k0", format(k0, ".17g"),
                            "--n", str(n),
                            "--L", format(L_m, ".17g"),
                            "--field-min-t", format(field_min, ".17g"),
                            "--field-max-t", format(field_max, ".17g"),
                            "--output-stem", stem,
                            "--output-dir", results_dir,
                            "--timer-log-dir", timer_logs_dir,
                            "--shape", shape_arg,
                            "--material-label", material_label,
                            "--sw-coercivity-t", format(sw_t, ".17g"),
                        ]
                        if shape_variant_arg:
                            argv.extend(["--shape-variant", shape_variant_arg])
                        if adaptive_enabled:
                            argv.extend(["--adaptive", "--adaptive-dh-min-t", format(dh_min, ".17g")])
                        if periodic:
                            argv.append("--periodic")
                        if not plot_enabled:
                            argv.append("--no-plot")

                        job = {
                            "batch": {"id": batch_id, "label": batch_label},
                            "argv": argv,
                            "working_directory": cwd,
                            "environment": {"PYTHONPATH": pythonpath, **openmp},
                            "stdout_path": str(stdout_path),
                            "stderr_path": str(stderr_path),
                            "receipt_path": str(receipt_path),
                            "runtime": {
                                "type": "conda",
                                "name": conda_name,
                                "executable": conda_executable,
                                "prefix": conda_prefix,
                            },
                            "metadata": {
                                "label": stem,
                                "experiment": "single_grain_coercivity",
                                "preset": "property_grid",
                                "shape": shape_token,
                                "mu0_Ms_T": mu0_ms_t,
                                "A0": a0_value,
                                "K0": k0,
                                "SW_T": sw_t,
                                "field_min_t": field_min,
                                "field_max_t": field_max,
                                "n": n,
                                "size_nm": size_nm,
                                "L_m": L_m,
                                "size_factor": size_factor,
                                "adaptive": adaptive_enabled,
                                "periodic": periodic,
                                "adaptive_dh_min_t": float(dh_min) if adaptive_enabled else None,
                                "conda_environment": conda_name,
                                "conda_prefix": conda_prefix,
                                "openmp": openmp,
                            },
                        }
                        job_file.write_text(json.dumps(job, indent=2) + "\n", encoding="utf-8")

                        if dry_run == "true":
                            job_id = "-"
                            status = "dry-run"
                            if len(dry_run_examples) < 3:
                                command = [
                                    conda_executable, "run", "--no-capture-output",
                                    "--prefix", conda_prefix, *argv,
                                ]
                                dry_run_examples.append([shlex.join(command)])
                        else:
                            job_id = subprocess.check_output(
                                [sys.executable, scheduler, "submit", "--job-file", str(job_file)],
                                text=True,
                            ).strip()
                            status = "pending"

                        rows.append([
                            job_id, status, "property_grid", shape_token, mu0_ms_t, k0, sw_t,
                            field_min, field_max, adaptive_enabled, periodic, dh_min, n,
                            size_nm, L_m, size_factor, n**3, stdout_path, stderr_path,
                            receipt_path,
                        ])

with Path(summary).open("w", encoding="utf-8") as handle:
    handle.write("\t".join(summary_columns) + "\n")
    for row in rows:
        handle.write("\t".join(str(item) for item in row) + "\n")

print(
    f"{'DRY RUN: ' if dry_run == 'true' else 'Queued '}"
    f"property_grid/{shape_token}: {len(rows)} jobs"
)
print(
    f"  Ms values: {len(ms_values)} ({min(ms_values):g}..{max(ms_values):g} T), "
    f"K0 values: {len(k0_values)} ({min(k0_values):g}..{max(k0_values):g} J/m^3)"
)
print(f"  field_min_t range: {min(row[7] for row in rows):g}..{max(row[7] for row in rows):g} T")
print(f"  field_max_t: {field_max:g} T")
if dry_run_examples:
    print("  Example command:")
    print("  " + dry_run_examples[0][0])
PY

  echo "Batch directory: $batch_dir"
  echo "Manifest: $manifest"
  echo "Summary: $summary"
}

if [[ "$PROPERTY_GRID" == true ]]; then
  "$PYTHON_BIN" - "$PROPERTY_MS_RANGE_T" "$PROPERTY_K0_RANGE" "$PROPERTY_A0" \
    "$PROPERTY_FIELD_MAX_T" <<'PY'
import math
import sys
from decimal import Decimal

def parse_range(spec, name):
    parts = spec.split(":")
    if len(parts) != 3:
        raise SystemExit(f"{name} must use START:STOP:STEP syntax")
    start, stop, step = (Decimal(part) for part in parts)
    if step <= 0:
        raise SystemExit(f"{name} step must be positive")
    if stop < start:
        raise SystemExit(f"{name} stop must be >= start")
    values = []
    value = start
    epsilon = step / Decimal("1000000")
    while value <= stop + epsilon:
        values.append(value)
        value += step
    if not values:
        raise SystemExit(f"{name} produced no values")
    return values

ms_range, k0_range, a0, field_max = sys.argv[1:]
ms_values = parse_range(ms_range, "ms range")
k0_values = parse_range(k0_range, "K0 range")
numbers = [Decimal(a0), Decimal(field_max)]
if any(value <= 0 for value in [numbers[0]]) or not math.isfinite(float(field_max)):
    raise SystemExit("A0 must be positive and field max must be finite")
if any(value <= 0 for value in ms_values + k0_values):
    raise SystemExit("Ms and K0 range values must be positive")
print(f"Property grid: {len(ms_values)} Ms values x {len(k0_values)} K0 values")
PY
fi

for material_record in "${material_records[@]}"; do
  for shape in "${shape_values[@]}"; do
    create_material_batch "$material_record" "$shape"
  done
done

if [[ "$PROPERTY_GRID" == true ]]; then
  for shape in "${shape_values[@]}"; do
    create_property_grid_batch "$shape"
  done
fi

if [[ "$DRY_RUN" == true ]]; then
  echo "Dry run complete: no jobs were added to the scheduler queue."
fi
