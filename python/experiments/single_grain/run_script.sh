#!/usr/bin/env bash
set -u

# ============================================================
# Single-grain grain-size/resolution coercivity sweep
# ============================================================
# Physical setup implemented by single_grain_coercivity.py:
#   Lx = Ly = Lz = L, nx = ny = nz = n, one material grain,
#   H_ext along z, and the easy axis tilted slightly away from z.
#
# Grain sizes are defined by SIZE_FACTORS relative to the characteristic
# exchange/demagnetisation scale
#   ell = sqrt(A0 / (0.5 * mu0 * Ms**2))
# using the same SI-unit material defaults as the grain examples
# (A0 = 7.7e-12 J/m, mu0*Ms = 1.61 T).  The Python setup computes
# L = factor * ell to avoid unit conversion mistakes in this shell script.
#
# The default RESOLUTIONS span roughly 1000--10000 cubic cells:
#   10^3, 12^3, 14^3, 16^3, 18^3, 20^3, 22^3.
#
# FMM is disabled by default.  Set USE_FMM=true to add --use-fmm while
# preserving the same physical setup and output naming convention.
# ============================================================

PYTHON_BIN="python"
PY_SCRIPT="single_grain_coercivity.py"

RESOLUTIONS=(10 12 14 16 18 20 22)
SIZE_FACTORS=(0.5 0.75 1.0 1.25 1.5 2.0)
USE_FMM=false
USE_CUDA=false
DO_PLOT=true

TILT_DEGREES=3.0
OUTPUT_ROOT="results"
LOGDIR="logs"

# Hysteresis field sweep in Tesla-equivalent values.  The Python setup converts
# these to A/m for MagTense, matching the existing grain hysteresis convention.
FIELD_MAX_T=1.0
FIELD_MIN_T=-7.0
FIELD_STEP_T=-0.1

# Optional FMM controls for future comparisons when USE_FMM=true.
# Defaults match the grain and grain_perf command-line defaults.
FMM_CPN=660
FMM_EPS="1e-4"
IFUNIF=1
NLMIN=1
NLMAX=5
ALLOW_FMM_SHORT_CIRCUIT=0
FMM_MIN_N=20000
FMM_NTERMS=-1

# Ensure the local package sources are importable when this script is run from
# the experiment directory without installing MagTense as a wheel.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${SCRIPT_DIR}/../../src:${PYTHONPATH:-}"

mkdir -p "$LOGDIR" "$OUTPUT_ROOT"
timestamp="$(date +%Y%m%d_%H%M%S)"
summary="${LOGDIR}/summary_single_grain_${timestamp}.tsv"
echo -e "size_factor\tresolution\tntot\tfmm\texit_code\toutput_dir" > "$summary"

fail_count=0

for size_factor in "${SIZE_FACTORS[@]}"; do
  for n in "${RESOLUTIONS[@]}"; do
    fmm_tag="fmm_off"
    cmd=("$PYTHON_BIN" "$PY_SCRIPT" --size-factor "$size_factor" --n "$n")
    cmd+=(--tilt-degrees "$TILT_DEGREES")
    cmd+=(--field-max-t "$FIELD_MAX_T" --field-min-t "$FIELD_MIN_T" --field-step-t "$FIELD_STEP_T")
    cmd+=(--output-dir "$OUTPUT_ROOT")

    if [[ "$USE_FMM" == true ]]; then
      fmm_tag="fmm_on"
      cmd+=(--use-fmm --fmm-cpn "$FMM_CPN" --fmm-eps "$FMM_EPS")
      cmd+=(--ifunif "$IFUNIF" --nlmin "$NLMIN" --nlmax "$NLMAX")
      cmd+=(--allow-fmm-short-circuit "$ALLOW_FMM_SHORT_CIRCUIT")
      cmd+=(--fmm-min-n "$FMM_MIN_N" --fmm-nterms "$FMM_NTERMS")
    fi

    if [[ "$USE_CUDA" == true ]]; then
      cmd+=(--cuda)
    fi

    if [[ "$DO_PLOT" != true ]]; then
      cmd+=(--no-plot)
    fi

    tag="single_grain_sf${size_factor}_n${n}_${fmm_tag}_${timestamp}"
    out="${LOGDIR}/${tag}.out"
    err="${LOGDIR}/${tag}.err"

    echo "=================================================="
    echo "Size factor : $size_factor"
    echo "Resolution  : $n ($((n*n*n)) cells)"
    echo "FMM         : $USE_FMM"
    echo "Command     : ${cmd[*]}"
    echo "Started at  : $(date)"

    "${cmd[@]}" >"$out" 2>"$err"
    rc=$?

    echo "Finished at : $(date)"
    echo "Exit code   : $rc"
    echo "stdout      : $out"
    echo "stderr      : $err"

    echo -e "${size_factor}\t${n}\t$((n*n*n))\t${USE_FMM}\t${rc}\t${OUTPUT_ROOT}" >> "$summary"
    if (( rc != 0 )); then
      fail_count=$((fail_count + 1))
    fi
  done
done

echo "Summary written to: $summary"
if (( fail_count > 0 )); then
  echo "ERROR: $fail_count run(s) failed." >&2
  exit 1
fi
