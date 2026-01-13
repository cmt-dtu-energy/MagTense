#!/usr/bin/env bash
# Sweep base resolutions and target FMM tree levels.
# Continues on errors.
# Records peak CPU RAM via /usr/bin/time -v and peak GPU VRAM via nvidia-smi polling.
# Echoes run info to terminal AND saves the same info to a per-run .runlog.

CPU_RANGE="0-15"
PY_SCRIPT="grain_hysteresis_perf.py"
MESH_TEMPLATE="Grid_rasBase_%d_nGrains_25_nRef_1_dG_3.75e-09.mat"

#BASES=(10 13 16 20 25 32 39)
#BASES=(32 40)
#BASES=(31)
#BASES=(10 13 16 20 25 31)
#BASES=(31 25 20 16 13 10)
BASES=(46)
LEVELS=(3)
#LEVELS=(1 2 3)

#BASES=(20 25 32 40)
#LEVELS=(1 2 3)

#BASES=(16)
#LEVELS=(3)

#LOGDIR="logs_CUDA"
LOGDIR="logs_fmm_levels"
#LOGDIR="logs_fmm_large"
mkdir -p "$LOGDIR"

timestamp="$(date +%Y%m%d_%H%M%S)"
summary="${LOGDIR}/summary_${timestamp}.tsv"
echo -e "tag\tN\tL\tcpn\texit_code\tmax_rss_kb\tmax_gpu_mem_mib" > "$summary"

# Prefer GNU time binary path
TIMEBIN="/usr/bin/time"
if [[ ! -x "$TIMEBIN" ]]; then
  echo "ERROR: /usr/bin/time not found/executable. Install GNU time (package often named 'time')." >&2
  exit 1
fi

# GPU polling interval (seconds)
GPU_POLL_INTERVAL="60.0"
#GPU_POLL_INTERVAL="5.0"
# GPU_POLL_INTERVAL="1.0"


pick_cpn () {
  local total="$1"
  local L="$2"

  local denom_low=1
  for ((i=0; i<L; i++)); do denom_low=$((denom_low * 8)); done
  local cpn_min=$(( (total + denom_low - 1) / denom_low ))
  (( cpn_min < 1 )) && cpn_min=1

  local denom_up=1
  for ((i=0; i<L-1; i++)); do denom_up=$((denom_up * 8)); done
  local cpn_max=$(( total / denom_up - 1 ))
  (( cpn_max < 1 )) && cpn_max=1

  if (( cpn_min > cpn_max )); then
    echo "$cpn_min"
  else
    echo $(( (cpn_min + cpn_max) / 2 ))
  fi
}

# Monitor GPU memory for a PID, record peak MiB into a file.
# Auto-detects whether nvidia-smi supports used_gpu_memory or used_memory.
gpu_monitor_pid () {
  local pid="$1"
  local outfile="$2"
  local interval="$3"

  if ! command -v nvidia-smi >/dev/null 2>&1; then
    echo "NA" > "$outfile"
    return 0
  fi

  # Detect supported field name
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
          {
            gsub(/^[ \t]+|[ \t]+$/, "", $1);
            gsub(/^[ \t]+|[ \t]+$/, "", $2);
            if ($1 == P && $2 ~ /^[0-9]+$/) s += $2;
          }
          END { print s+0 }'
    )

    if [[ "$used_sum" =~ ^[0-9]+$ ]] && (( used_sum > peak )); then
      peak="$used_sum"
    fi

    sleep "$interval"
  done

  echo "$peak" > "$outfile"
}

fail_count=0

for N in "${BASES[@]}"; do
  total=$(( N*N*N ))
  mesh_file=$(printf "$MESH_TEMPLATE" "$N")

  for L in "${LEVELS[@]}"; do
    cpn="$(pick_cpn "$total" "$L")"
    tag="N${N}_L${L}_cpn${cpn}_${timestamp}"

    out="${LOGDIR}/${tag}.out"
    err="${LOGDIR}/${tag}.err"
    timelog="${LOGDIR}/${tag}.time"
    gpulog="${LOGDIR}/${tag}.gpu_peak"
    runlog="${LOGDIR}/${tag}.runlog"
    pidfile="${LOGDIR}/${tag}.pid"

    # log() prints to terminal AND appends to runlog
    log () { echo "$@" | tee -a "$runlog"; }

    # command for running with fmm
    cmd=(taskset -c "$CPU_RANGE" python "$PY_SCRIPT" --use-fmm --cuda --fmm-cpn "$cpn" --mesh-file "$mesh_file")

    # command for running without fmm
    #cmd=(taskset -c "$CPU_RANGE" python "$PY_SCRIPT" --cuda --fmm-cpn "$cpn" --mesh-file "$mesh_file")

    # Clear/initialise runlog
    : > "$runlog"

    log "=================================================="
    log "Tag            : $tag"
    log "Base resolution: $N x $N x $N"
    log "Total cells    : $total"
    log "Target FMM lvl : $L"
    log "Chosen cpn     : $cpn"
    log "Mesh file      : $mesh_file"
    log "Stdout         : $out"
    log "Stderr         : $err"
    log "Time log       : $timelog"
    log "GPU peak file  : $gpulog"
    log "Started at     : $(date)"
    log "Command        : ${cmd[*]}"
    log "--------------------------------------------------"

    # Remove stale pidfile
    rm -f "$pidfile"

    # Start timed wrapper in background.
    # IMPORTANT: $! here is the PID of /usr/bin/time, NOT python.
    # We write the python PID to pidfile from inside the wrapper.
    "$TIMEBIN" -v -o "$timelog" bash -c "
      ${cmd[*]} &
      echo \$! > $(printf '%q' "$pidfile")
      wait \$!
    " >"$out" 2>"$err" &
    time_pid=$!

    # Wait until pidfile exists to get the *actual* python PID
    child_pid=""
    for _ in $(seq 1 200); do
      if [[ -s "$pidfile" ]]; then
        child_pid="$(cat "$pidfile" 2>/dev/null)"
        break
      fi
      sleep 0.01
    done

    if [[ -n "$child_pid" ]]; then
      log "time_pid       : $time_pid"
      log "python_pid     : $child_pid"
      gpu_monitor_pid "$child_pid" "$gpulog" "$GPU_POLL_INTERVAL" &
      mon_pid=$!
    else
      log "WARNING: Could not determine python PID (pidfile missing). GPU peak will be NA."
      echo "NA" > "$gpulog"
      mon_pid=""
    fi

    # Wait for timed wrapper to finish and capture exit code
    wait "$time_pid"
    rc=$?

    # Ensure monitor exits too
    if [[ -n "${mon_pid:-}" ]]; then
      wait "$mon_pid" 2>/dev/null
    fi

    # Extract max RSS in KB from GNU time output
    max_rss_kb=$(awk -F: '/Maximum resident set size/ {gsub(/^[ \t]+/, "", $2); print $2}' "$timelog")
    [[ -z "$max_rss_kb" ]] && max_rss_kb="NA"

    # Read GPU peak
    max_gpu_mem_mib=$(cat "$gpulog" 2>/dev/null || echo "NA")
    [[ -z "$max_gpu_mem_mib" ]] && max_gpu_mem_mib="NA"

    log "--------------------------------------------------"
    log "Finished at    : $(date)"
    log "Exit code      : $rc"
    log "Max RSS (KB)   : $max_rss_kb"
    log "Max GPU (MiB)  : $max_gpu_mem_mib"

    echo -e "${tag}\t${N}\t${L}\t${cpn}\t${rc}\t${max_rss_kb}\t${max_gpu_mem_mib}" >> "$summary"

    if (( rc != 0 )); then
      log "STATUS         : FAILED"
      fail_count=$((fail_count + 1))
    else
      log "STATUS         : OK"
    fi
  done
done

echo
echo "All runs finished. Failures: $fail_count"
echo "Logs in: $LOGDIR"
echo "Summary: $summary"
