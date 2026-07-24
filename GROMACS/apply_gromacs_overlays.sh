#!/usr/bin/env bash
# Apply STELLAR-local GROMACS patches onto the bundled external_sw tree.
# Run after ./fetch_assets.sh on a new cluster checkout.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="${ROOT}/gromacs_overlays"
DST="${ROOT}/external_sw/gromacs"

if [ ! -d "${DST}" ]; then
  echo "ERROR: ${DST} not found. Run ./fetch_assets.sh first." >&2
  exit 1
fi

echo "Applying GROMACS overlays from ${SRC} -> ${DST}"
cp -f "${SRC}/config_files/nvt.sh" "${DST}/config_files/nvt.sh"
cp -f "${SRC}/config_files/npt.sh" "${DST}/config_files/npt.sh"
cp -f "${SRC}/config_files/simulation.sh" "${DST}/config_files/simulation.sh"
cp -f "${SRC}/topology/fix_gro_zero_coords.py" "${DST}/topology/fix_gro_zero_coords.py"
cp -f "${SRC}/topology/options/components/tools_topology.py" "${DST}/topology/options/components/tools_topology.py"
echo "Done."
