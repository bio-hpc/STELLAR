#!/usr/bin/env bash
# run_loop.sh [CASE ...]
#
# Detached driver of the STELLAR agentic loop over the improve5Frag_2 cases.
# For each case it runs run_case.py (stage -> pipeline -> diagnose/fix -> success
# check), one at a time, waiting for the case's SLURM MD jobs to drain before the
# next. Safe to launch with nohup; it does not depend on an interactive terminal.
#
#   nohup bash agentic_loop_stellar/run_loop.sh > agentic_loop_stellar/run_loop.out 2>&1 &
#
# Success per case: all_metrics_GN_<CASE>.csv with >=1 rmsd_md != NA.
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
AL="${ROOT}/agentic_loop_stellar"
PY="${AL}/bin/python3"
SUMMARY="${AL}/e2e_summary.txt"
# Source subset dir (override to run a different 50-case subset). Exported so the
# staging step (run_case.py -> stage_case.sh) reads from the same place.
export CASES_DIR="${CASES_DIR:-${ROOT}/improve5Frag_2/cases}"

# Pilot knobs (small = fast). Override via env.
MAX_VALID="${MAX_VALID:-2}"
MAX_ORGANIZE="${MAX_ORGANIZE:-400}"
OVERLAP_TOL="${OVERLAP_TOL:-450}"
# Group A: candidate pool carried to MD so it can retry other geometries when
# some combinations hit GROMACS zero-length constraints (0 = same as MAX_VALID).
MD_POOL="${MD_POOL:-0}"

log(){ echo "[$(date +%H:%M:%S)] $*"; }

# Case list: args, or CASES env (space-separated), or all 50.
if [ "$#" -gt 0 ]; then
  CASES=("$@")
elif [ -n "${CASES:-}" ]; then
  # shellcheck disable=SC2206
  CASES=(${CASES})
else
  mapfile -t CASES < <(ls -1 "${CASES_DIR}" 2>/dev/null)
fi

: > "$SUMMARY"
log "STELLAR agentic loop: ${#CASES[@]} case(s). max_valid=${MAX_VALID} tol=${OVERLAP_TOL}"
echo "# case: result (detail)" >> "$SUMMARY"

drain_queue(){
  for _ in $(seq 1 360); do
    n=$(squeue -u "$USER" -h 2>/dev/null | wc -l)
    [ "$n" -eq 0 ] && return 0
    sleep 10
  done
}

clean_case(){
  local c="$1"
  rm -rf "${ROOT}/${c}_results" "${ROOT}/${c}_GN" \
         "${ROOT}/valid_combinations_GN_${c}" "${ROOT}/valid_GN_${c}_final" \
         "${ROOT}/md_rmsd_peptides_${c}" "${ROOT}"/VS_GR_*"${c}"* \
         "${ROOT}"/*_"${c}".csv "${ROOT}/all_metrics_GN_${c}.csv" 2>/dev/null
  # also match lowercase receptor id in VS_GR_* names
  local lc; lc="$(echo "${c%%_*}" | tr '[:upper:]' '[:lower:]')"
  rm -rf "${ROOT}"/VS_GR_*"${lc}"_* 2>/dev/null
}

drain_queue
for c in "${CASES[@]}"; do
  [ -d "${CASES_DIR}/${c}" ] || { echo "${c}: SKIP (no case dir)" >> "$SUMMARY"; continue; }
  log "=== ${c}: start ==="
  clean_case "$c"
  "$PY" "${AL}/run_case.py" "$c" \
      --max-valid "$MAX_VALID" --md-pool "$MD_POOL" \
      --max-organize "$MAX_ORGANIZE" --overlap-tol "$OVERLAP_TOL" \
      > "${AL}/logs/${c}.out" 2>&1
  rc=$?
  line=$(grep -E "^RESULT ${c}:" "${AL}/logs/${c}.out" | tail -1)
  [ -z "$line" ] && line="RESULT ${c}: FAIL (driver rc=${rc})"
  echo "${line#RESULT }" >> "$SUMMARY"
  log "=== ${c}: ${line#RESULT ${c}: } ==="
  drain_queue
done

log "ALL CASES DONE."
{
  echo "----"
  npass=$(grep -c ": PASS" "$SUMMARY" 2>/dev/null || echo 0)
  nfail=$(grep -c ": FAIL" "$SUMMARY" 2>/dev/null || echo 0)
  echo "PASS=${npass} FAIL=${nfail}"
  echo "finished: $(date)"
} >> "$SUMMARY"
