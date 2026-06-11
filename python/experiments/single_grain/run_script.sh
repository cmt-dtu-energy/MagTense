#!/usr/bin/env bash
set -u

# ============================================================
# Single-grain coercivity parameter scan
# ============================================================
#
# Scans:
#   --adaptive-dh-min-t   only when USE_ADAPTIVE=true
#   --n
#   --size-factor
#
# Optional flags:
#   USE_ADAPTIVE=true  -> add --adaptive and --adaptive-dh-min-t
#   USE_PERIODIC=true  -> add --periodic
#
# All other parameters are left at the defaults in:
#   single_grain_coercivity.py
#
# Example adaptive + periodic command:
#   python single_grain_coercivity.py \
#       --size-factor 15 \
#       --n 10 \
#       --adaptive \
#       --adaptive-dh-min-t 0.001 \
#       --periodic
#
# Example non-adaptive + non-periodic command:
#   python single_grain_coercivity.py \
#       --size-factor 15 \
#       --n 10
# ============================================================

PYTHON_BIN="python"
PY_SCRIPT="single_grain_coercivity.py"

# Top-level run-mode flags
USE_ADAPTIVE=true
USE_PERIODIC=true

# Smallest dh_min first, so the finest adaptive stepping runs first.
# This array is only used when USE_ADAPTIVE=true.
DH_MIN_VAL=(0.001 0.01 0.1)

#RESOLUTIONS=(10 12 14 16 18 20)
RESOLUTIONS=(10 12 14)
SIZE_FACTORS=(1 5 10 15 20 25 30 40)


LOGDIR="logs"

# Ensure local MagTense sources are importable.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${SCRIPT_DIR}/../../src:${PYTHONPATH:-}"

mkdir -p "$LOGDIR"

timestamp="$(date +%Y%m%d_%H%M%S)"
summary="${LOGDIR}/summary_single_grain_${timestamp}.tsv"

echo -e "adaptive\tperiodic\tdh_min\tn\tsize_factor\tntot\texit_code\tstdout\tstderr" > "$summary"

fail_count=0
run_index=0

if [[ "$USE_ADAPTIVE" == true ]]; then
  n_runs=$((${#DH_MIN_VAL[@]} * ${#RESOLUTIONS[@]} * ${#SIZE_FACTORS[@]}))
else
  n_runs=$((${#RESOLUTIONS[@]} * ${#SIZE_FACTORS[@]}))
fi

make_stem() {
  local size_factor="$1"
  local n="$2"
  local dh_min="${3:-}"

  "$PYTHON_BIN" - "$size_factor" "$n" "$USE_PERIODIC" "$USE_ADAPTIVE" "$dh_min" <<'PY'
import sys
import math

size_factor = float(sys.argv[1])
n = int(sys.argv[2])
periodic = sys.argv[3].lower() == "true"
adaptive = sys.argv[4].lower() == "true"
dh_min_arg = sys.argv[5] if len(sys.argv) > 5 else ""

stem = f"single_grain_sf{size_factor:.2g}_n{n}"

if periodic:
    stem = stem.replace(f"_n{n}", f"_n{n}_P", 1)

if adaptive:
    dh_min = abs(float(dh_min_arg))
    stem += f"_A_FS{dh_min:.1e}"

print(stem)
PY
}

run_one() {
  local dh_min="$1"
  local n="$2"
  local size_factor="$3"

  run_index=$((run_index + 1))

  cmd=(
    "$PYTHON_BIN"
    "$PY_SCRIPT"
    --size-factor "$size_factor"
    --n "$n"
  )

  if [[ "$USE_ADAPTIVE" == true ]]; then
    cmd+=(--adaptive)
    cmd+=(--adaptive-dh-min-t "$dh_min")
  fi

  if [[ "$USE_PERIODIC" == true ]]; then
    cmd+=(--periodic)
  fi

  stem="$(make_stem "$size_factor" "$n" "$dh_min")"
  out="${LOGDIR}/${stem}_${timestamp}.out"
  err="${LOGDIR}/${stem}_${timestamp}.err"

  echo "=================================================="
  echo "Run              : ${run_index}/${n_runs}"
  echo "Adaptive         : $USE_ADAPTIVE"
  echo "Periodic         : $USE_PERIODIC"
  if [[ "$USE_ADAPTIVE" == true ]]; then
    echo "dh_min           : $dh_min"
  fi
  echo "Resolution       : $n ($((n*n*n)) cells)"
  echo "Size factor      : $size_factor"
  echo "Log stem         : $stem"
  echo "Command          : ${cmd[*]}"
  echo "Started at       : $(date)"

  "${cmd[@]}" >"$out" 2>"$err"
  rc=$?

  echo "Finished at      : $(date)"
  echo "Exit code        : $rc"
  echo "stdout           : $out"
  echo "stderr           : $err"

  echo -e "${USE_ADAPTIVE}\t${USE_PERIODIC}\t${dh_min}\t${n}\t${size_factor}\t$((n*n*n))\t${rc}\t${out}\t${err}" >> "$summary"

  if (( rc != 0 )); then
    fail_count=$((fail_count + 1))
  fi
}

if [[ "$USE_ADAPTIVE" == true ]]; then
  for dh_min in "${DH_MIN_VAL[@]}"; do
    for n in "${RESOLUTIONS[@]}"; do
      for size_factor in "${SIZE_FACTORS[@]}"; do
        run_one "$dh_min" "$n" "$size_factor"
      done
    done
  done
else
  for n in "${RESOLUTIONS[@]}"; do
    for size_factor in "${SIZE_FACTORS[@]}"; do
      run_one "none" "$n" "$size_factor"
    done
  done
fi

echo "=================================================="
echo "Summary written to: $summary"

if (( fail_count > 0 )); then
  echo "ERROR: $fail_count run(s) failed." >&2
  exit 1
fi

echo "All runs completed successfully."