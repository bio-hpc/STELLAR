#!/usr/bin/env bash
set -u
cd /data4/alejandro/STELLAR_JCIM/STELLAR
for d in improve5Frag_6 improve5Frag_7; do
  echo "############################################################"
  echo "### $(date '+%F %T')  START $d (--no-ai, workers 3)"
  echo "############################################################"
  CASES_DIR="$PWD/$d/cases" agentic_loop_stellar/run_auto_fix.sh --no-ai --workers 3 \
      > "agentic_loop_stellar/run_auto_fix_${d}.out" 2>&1
  echo "### $(date '+%F %T')  END $d (rc=$?)"
  echo "### disco: $(df -h /data4 | tail -1)"
done
echo "### $(date '+%F %T')  ALL DONE"
