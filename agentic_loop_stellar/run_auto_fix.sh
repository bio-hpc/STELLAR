#!/usr/bin/env bash
# Launch the parallel AI-driven STELLAR agentic loop (auto_fix_loop.py).
# Wires the Cursor SDK env el7 needs: 3.11 venv for the orchestrator + a
# CentOS7-compatible node for the SDK bridge. run_case.py stays on the 3.9 shim.
#
#   export CURSOR_API_KEY="cursor_..."
#   agentic_loop_stellar/run_auto_fix.sh --workers 3 --max-fix-attempts 1 3HXI_C 6IFJ_C
#   CASES_DIR=/path/other/cases agentic_loop_stellar/run_auto_fix.sh --workers 4
#   agentic_loop_stellar/run_auto_fix.sh --no-ai 3HXI_C   # smoke test, no key/SDK
set -u
AL="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYSDK="${AL}/.venv-sdk/bin/python"
[ -x "$PYSDK" ] || { echo "ERROR: SDK venv missing ($PYSDK)." >&2; exit 2; }
export CURSOR_SDK_BRIDGE_BIN="${CURSOR_SDK_BRIDGE_BIN:-${AL}/.sdkenv/cursor-sdk-bridge}"
exec "$PYSDK" "${AL}/auto_fix_loop.py" "$@"
