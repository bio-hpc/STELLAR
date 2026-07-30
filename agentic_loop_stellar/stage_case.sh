#!/usr/bin/env bash
# stage_case.sh <CASE>
#
# Stage a case from improve5Frag_2/cases/<CASE>/ into a project-root working dir
# <CASE>_GN/ WITHOUT touching the original case data.
#
# Strategy (write-isolated symlinks):
#   <CASE>_GN/VS_GN_<CASE>_FragN/
#       <every original entry>   -> symlink to original (read-only reads)
#       molecules/               -> REAL dir; each pose file symlinked in.
#                                   New files written by step 1 (Pose*_coordinates_CN.csv)
#                                   therefore land in the staging dir, not the original.
#
# Also symlinks receptor_original.pdb for reference. The crystal (complex/) and
# peptide fragments (peptide_pdb_fragments/<CASE>/) already live at project root.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CASE="${1:?usage: stage_case.sh <CASE>}"
# Source subset dir is configurable so a new 50-case subset can be staged without
# touching the scripts:  CASES_DIR=/path/to/subset/cases bash stage_case.sh <CASE>
CASES_DIR="${CASES_DIR:-${ROOT}/improve5Frag_2/cases}"
SRC="${CASES_DIR}/${CASE}"
DST="${ROOT}/${CASE}_GN"

if [ ! -d "$SRC" ]; then
  echo "ERROR: case source not found: $SRC" >&2
  exit 2
fi

# Fresh staging every time (only removes the staging dir, never the source).
rm -rf "$DST"
mkdir -p "$DST"

shopt -s nullglob
for fragdir in "$SRC"/VS_GN_"${CASE}"_Frag*; do
  [ -d "$fragdir" ] || continue
  name="$(basename "$fragdir")"
  sfrag="$DST/$name"
  mkdir -p "$sfrag"
  for entry in "$fragdir"/* "$fragdir"/.[!.]*; do
    [ -e "$entry" ] || continue
    base="$(basename "$entry")"
    if [ "$base" = "molecules" ] && [ -d "$entry" ]; then
      # Real molecules dir with per-file symlinks (write isolation for step 1).
      # Skip regenerable step-1 outputs (*coordinates_CN.csv) so save_pose can
      # rewrite them into the staging dir instead of following a symlink back
      # into the original data.
      mkdir -p "$sfrag/molecules"
      for mol in "$entry"/*; do
        [ -e "$mol" ] || continue
        mb="$(basename "$mol")"
        case "$mb" in
          *coordinates_CN.csv) continue ;;
        esac
        ln -sf "$mol" "$sfrag/molecules/$mb"
      done
    else
      ln -sf "$entry" "$sfrag/$base"
    fi
  done
done

# Reference receptor (not strictly required; kept for traceability).
if [ -f "$SRC/receptor_original.pdb" ]; then
  ln -sf "$SRC/receptor_original.pdb" "$DST/receptor_original.pdb"
fi

nfrags=$(find "$DST" -maxdepth 1 -type d -name "VS_GN_${CASE}_Frag*" | wc -l)
echo "Staged ${CASE}: ${nfrags} fragment dirs -> ${DST}"
