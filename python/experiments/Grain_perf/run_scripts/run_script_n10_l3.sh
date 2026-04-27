#!/usr/bin/env bash
set -u

# ============================================================
# Settings (ZIPPED INPUTS)
# Every *_LIST must have the same length.
# Index i defines one run.
# ============================================================

CPU_RANGE="0-15"
PYTHON_BIN="python"
PY_SCRIPT="grain_hysteresis_perf.py"

USE_FMM=1
USE_CUDA=1
# Optional label per run (for human readability in tags)
LABEL_LIST=("n10" "n10" "n10" "n10" "n10" "n10" "n10" "n10" "n10")

# Mesh parameter per run (uniform grid perf: rasBase only)
BASES_LIST=(10 13 16 20 25 31 32 35 40)

# FMM controls per run
CPN_LIST=(10 10 10 10 10 10 10 10 10)
IFUNIF_LIST=(1 1 1 1 1 1 1 1 1)                 # 1 uniform, 0 adaptive
NLMIN_LIST=(1 1 1 1 1 1 1 1 1 )
NLMAX_LIST=(3 3 3 3 3 3 3 3 3)
ALLOW_FMM_SHORT_CIRCUIT_LIST=(0 0 0 0 0 0 0 0 0)  # 0/1 (int-style)
FMM_MIN_N_LIST=(30000 30000 30000 30000 30000 30000 30000 30000 30000)
FMM_NTERMS_LIST=(10 10 10 10 10 10 10 10 10)

# Mesh file naming (uniform perf naming)
MESH_TEMPLATE="Grid_rasBase_%d_nGrains_25_nRef_1_dG_3.75e-09.mat"

LOGDIR="logs_perf"
GPU_POLL_INTERVAL="5.0"
TIMEBIN="/usr/bin/time"

# ============================================================

mkdir -p "$LOGDIR"

timestamp="$(date +%Y%m%d_%H%M%S)"
summary="${LOGDIR}/summary_${timestamp}.tsv"
echo -e "tag\tlabel\tN\tcpn\tifunif\tnlmin\tnlmax\tshort_circuit\tfmm_min_n\tnterms\texit_code\tmax_rss_kb\tmax_gpu_mem_mib\tmesh_file" > "$summary"

if [[ ! -x "$TIMEBIN" ]]; then
  echo "ERROR: $TIMEBIN not found/executable." >&2
  exit 1
fi

# ---- validate zipped list lengths ----
len=${#BASES_LIST[@]}
for name in LABEL_LIST CPN_LIST IFUNIF_LIST NLMIN_LIST NLMAX_LIST ALLOW_FMM_SHORT_CIRCUIT_LIST FMM_MIN_N_LIST FMM_NTERMS_LIST; do
  eval "this_len=\${#${name}[@]}"
  if (( this_len != len )); then
    echo "ERROR: Array length mismatch. BASES_LIST has $len but $name has $this_len." >&2
    exit 2
  fi
done

current_runlog=""

log () {
  if [[ -n "${current_runlog:-}" ]]; then
    echo "$@" | tee -a "$current_runlog"
  else
    echo "$@"
  fi
}

gpu_monitor_pid () {
  local pid="$1"
  local outfile="$2"
  local interval="$3"

  if ! command -v nvidia-smi >/dev/null 2>&1; then
    echo "NA" > "$outfile"
    return 0
  fi

  local memfield=""
  if nvidia-smi --help-query-compute-apps 2>/dev/null | grep -q "used_gpu_memory"; then
    memfield="used_gpu_memory"
  elif nvidia-smi --help-query-compute-apps 2>/dev/null | grep -q "used_memory"; then
    memfield="used_memory"
  else
    memfield="used_gpu_memory"
  fi

  local peak=0
  while kill -0 "$pid" 2>/dev/null; do
    local used_sum
    used_sum=$(
      nvidia-smi --query-compute-apps=pid,${memfield} --format=csv,noheader,nounits 2>/dev/null |
      awk -v P="$pid" -F',' '
        { gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $2);
          if ($1 == P && $2 ~ /^[0-9]+$/) s += $2 }
        END { print s+0 }'
    )
    (( used_sum > peak )) && peak="$used_sum"
    sleep "$interval"
  done

  echo "$peak" > "$outfile"
}

fail_count=0

# ============================================================
# Run zipped combinations
# ============================================================

for ((i=0; i<len; i++)); do
  label="${LABEL_LIST[i]}"
  N="${BASES_LIST[i]}"

  cpn="${CPN_LIST[i]}"
  ifunif="${IFUNIF_LIST[i]}"
  nlmin="${NLMIN_LIST[i]}"
  nlmax="${NLMAX_LIST[i]}"
  allow_sc="${ALLOW_FMM_SHORT_CIRCUIT_LIST[i]}"
  fmm_min_n="${FMM_MIN_N_LIST[i]}"
  fmm_nterms="${FMM_NTERMS_LIST[i]}"

  mesh_file="$(printf "$MESH_TEMPLATE" "$N")"

  tag="perf_${label}_N${N}_cpn${cpn}_ifunif${ifunif}_nl${nlmin}-${nlmax}_sc${allow_sc}_minN${fmm_min_n}_nterms${fmm_nterms}_${timestamp}"

  out="${LOGDIR}/${tag}.out"
  err="${LOGDIR}/${tag}.err"
  timelog="${LOGDIR}/${tag}.time"
  gpulog="${LOGDIR}/${tag}.gpu_peak"
  current_runlog="${LOGDIR}/${tag}.runlog"
  pidfile="${LOGDIR}/${tag}.pid"

  : > "$current_runlog"
  rm -f "$pidfile"

  cmd=(taskset -c "$CPU_RANGE" "$PYTHON_BIN" "$PY_SCRIPT")

  (( USE_FMM == 1 )) && cmd+=(--use-fmm --fmm-cpn "$cpn")
  (( USE_CUDA == 1 )) && cmd+=(--cuda)

  cmd+=(--mesh-file "$mesh_file")
  cmd+=(--ifunif "$ifunif" --nlmin "$nlmin" --nlmax "$nlmax")
  cmd+=(--allow-fmm-short-circuit "$allow_sc" --fmm-min-n "$fmm_min_n" --fmm-nterms "$fmm_nterms")

  log "=================================================="
  log "Index          : $i / $((len-1))"
  log "Tag            : $tag"
  log "Label          : $label"
  log "rasBase (N)    : $N"
  log "CPN            : $cpn"
  log "ifunif         : $ifunif"
  log "nlmin          : $nlmin"
  log "nlmax          : $nlmax"
  log "short_circuit  : $allow_sc"
  log "fmm_min_n      : $fmm_min_n"
  log "fmm_nterms     : $fmm_nterms"
  log "Mesh file      : $mesh_file"
  log "Command        : ${cmd[*]}"
  log "Started at     : $(date)"
  log "--------------------------------------------------"

  "$TIMEBIN" -v -o "$timelog" bash -c '
    pidfile="$1"; shift
    "$@" &
    echo $! > "$pidfile"
    wait $!
  ' _ "$pidfile" "${cmd[@]}" >"$out" 2>"$err" &
  time_pid=$!

  child_pid=""
  for _ in $(seq 1 200); do
    if [[ -s "$pidfile" ]]; then
      child_pid="$(cat "$pidfile" 2>/dev/null || true)"
      break
    fi
    sleep 0.01
  done

  if [[ -n "$child_pid" ]]; then
    log "python_pid     : $child_pid"
    gpu_monitor_pid "$child_pid" "$gpulog" "$GPU_POLL_INTERVAL" &
    mon_pid=$!
  else
    log "WARNING: python PID not found (pidfile missing). GPU peak will be NA."
    echo "NA" > "$gpulog"
    mon_pid=""
  fi

  wait "$time_pid"
  rc=$?

  [[ -n "${mon_pid:-}" ]] && wait "$mon_pid" 2>/dev/null || true

  max_rss_kb="$(awk -F: '/Maximum resident set size/ {gsub(/^[ \t]+/, "", $2); print $2}' "$timelog")"
  [[ -z "$max_rss_kb" ]] && max_rss_kb="NA"

  max_gpu_mem_mib="$(cat "$gpulog" 2>/dev/null || echo "NA")"
  [[ -z "$max_gpu_mem_mib" ]] && max_gpu_mem_mib="NA"

  log "--------------------------------------------------"
  log "Finished at    : $(date)"
  log "Exit code      : $rc"
  log "Max RSS (KB)   : $max_rss_kb"
  log "Max GPU (MiB)  : $max_gpu_mem_mib"

  echo -e "${tag}\t${label}\t${N}\t${cpn}\t${ifunif}\t${nlmin}\t${nlmax}\t${allow_sc}\t${fmm_min_n}\t${rc}\t${max_rss_kb}\t${max_gpu_mem_mib}\t${mesh_file}" >> "$summary"

  if (( rc != 0 )); then
    log "STATUS         : FAILED"
    fail_count=$((fail_count + 1))
  else
    log "STATUS         : OK"
  fi
done

current_runlog=""

log ""
log "All runs finished."
log "Failures : $fail_count"
log "Logs     : $LOGDIR"
log "Summary  : $summary"
