#!/usr/bin/env python3
"""
Prepare a receptor PDB for GROMACS pdb2gmx.

Removes HETATM records, TER records, and OXT atoms that break topology
generation (e.g. 2B9H: Mg/ADP cofactors and C-terminal OXT).

The receptor is kept as a single chain on purpose: the downstream topology
assembler (ShuttleMol) is only robust for single-chain targets. Spurious
backbone bonds that pdb2gmx creates across missing loops are removed afterwards
by fix_gap_bonds.py (called from generate_topologies.py), which is what keeps
grompp from aborting at NVT with "distance between excluded atoms larger than
the cut-off".
"""

import argparse
import glob
import os
import shutil


def sanitize_target_pdb(src_path, dst_path=None, in_place=False):
    """
    Write a protein-only PDB suitable for pdb2gmx.

    Returns:
        str: path to sanitized PDB
    """
    if in_place:
        dst_path = src_path
        tmp_path = src_path + ".sanitize_tmp"
        out_path = tmp_path
    else:
        out_path = dst_path or (os.path.splitext(src_path)[0] + "_sanitized.pdb")

    kept = 0
    with open(src_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            record = line[:6].strip()
            if record in ("HETATM", "TER"):
                continue
            if record == "ATOM" and line[12:16].strip() == "OXT":
                continue
            fout.write(line)
            if record == "ATOM":
                kept += 1

    if kept == 0:
        if os.path.exists(out_path) and out_path != src_path:
            os.remove(out_path)
        raise ValueError(f"No ATOM records kept from {src_path}")

    if in_place:
        shutil.move(out_path, src_path)
        return src_path

    return out_path


def main():
    parser = argparse.ArgumentParser(description="Sanitize receptor PDB for GROMACS topology")
    parser.add_argument("inputs", nargs="+", help="PDB file(s) or glob patterns")
    parser.add_argument("--in-place", action="store_true", help="Overwrite input files")
    parser.add_argument("--suffix", default="_sanitized", help="Output suffix when not in-place")
    args = parser.parse_args()

    paths = []
    for item in args.inputs:
        if any(ch in item for ch in "*?[]"):
            paths.extend(sorted(glob.glob(item)))
        else:
            paths.append(item)

    if not paths:
        raise SystemExit("No input PDB files found")

    for src in paths:
        if not os.path.isfile(src):
            print(f"Skip (missing): {src}")
            continue
        if args.in_place:
            out = sanitize_target_pdb(src, in_place=True)
        else:
            base, ext = os.path.splitext(src)
            out = sanitize_target_pdb(src, dst_path=f"{base}{args.suffix}{ext}")
        print(f"OK: {out}")


if __name__ == "__main__":
    main()
