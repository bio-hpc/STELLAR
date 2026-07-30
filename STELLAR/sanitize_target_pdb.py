#!/usr/bin/env python3
"""
Prepare a receptor PDB for GROMACS pdb2gmx.

Removes HETATM records, TER records, and OXT atoms that break topology
generation (e.g. 2B9H: Mg/ADP cofactors and C-terminal OXT).

Also truncates protein residues whose side chains are incomplete in the crystal
(missing heavy atoms) down to the largest buildable stub (ALA if CB is present,
otherwise GLY). pdb2gmx cannot rebuild missing heavy atoms and aborts on them
(fatal "residue is missing atoms" without -missing), and its histidine
protonation heuristic aborts even earlier with "Incomplete ring in HISnnn" when
a HIS imidazole ring is truncated (e.g. 3FXV: HIS348 keeps only N/CA/C/O/CB).
Stubbing the incomplete residue keeps the backbone intact for the complex while
letting pdb2gmx build a valid topology.

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

# Backbone heavy atoms common to every amino acid.
_BACKBONE = {"N", "CA", "C", "O"}

# Side-chain heavy atoms expected for each standard amino acid (crystal naming;
# pdb2gmx later renames e.g. ILE CD1 -> CD itself). Used only to decide whether
# a residue is complete enough for pdb2gmx to build.
_SIDECHAIN_HEAVY = {
    "GLY": set(),
    "ALA": {"CB"},
    "SER": {"CB", "OG"},
    "CYS": {"CB", "SG"},
    "THR": {"CB", "OG1", "CG2"},
    "VAL": {"CB", "CG1", "CG2"},
    "LEU": {"CB", "CG", "CD1", "CD2"},
    "ILE": {"CB", "CG1", "CG2", "CD1"},
    "PRO": {"CB", "CG", "CD"},
    "MET": {"CB", "CG", "SD", "CE"},
    "PHE": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
    "TRP": {"CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "ASP": {"CB", "CG", "OD1", "OD2"},
    "ASN": {"CB", "CG", "OD1", "ND2"},
    "GLU": {"CB", "CG", "CD", "OE1", "OE2"},
    "GLN": {"CB", "CG", "CD", "OE1", "NE2"},
    "HIS": {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    "LYS": {"CB", "CG", "CD", "CE", "NZ"},
    "ARG": {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
}

# Protonation-state / naming variants that are all histidine (and other common
# alternates) mapped to the canonical residue name used in _SIDECHAIN_HEAVY.
_RESNAME_ALIASES = {
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
    "HISD": "HIS", "HISE": "HIS", "HISH": "HIS",
    "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "CYX": "CYS", "CYM": "CYS",
    "GLH": "GLU", "ASH": "ASP", "LYN": "LYS",
    "MSE": "MET",
}


def _canonical_resname(resname):
    resname = resname.strip().upper()
    return _RESNAME_ALIASES.get(resname, resname)


def _is_hydrogen(atom_name, element):
    element = element.strip()
    if element:
        return element.upper() == "H"
    # No element column: fall back to the atom name (may be prefixed by a digit).
    name = atom_name.strip().lstrip("0123456789")
    return name[:1].upper() == "H"


def _residue_key(line):
    # (chainID, resSeq, iCode) uniquely identifies a residue within the file.
    return (line[21], line[22:26], line[26])


def _plan_residue_truncations(atom_lines):
    """
    Decide which residues have incomplete side chains and how to stub them.

    Returns a dict: residue_key -> ("ALA"|"GLY", set(atom_names_to_keep)).
    Complete residues, non-standard residues and residues with an incomplete
    backbone are left untouched (absent from the dict).
    """
    heavy_present = {}
    resname_of = {}
    for line in atom_lines:
        key = _residue_key(line)
        atom_name = line[12:16].strip()
        element = line[76:78] if len(line) >= 78 else ""
        if _is_hydrogen(atom_name, element):
            continue
        heavy_present.setdefault(key, set()).add(atom_name)
        resname_of.setdefault(key, _canonical_resname(line[17:20]))

    plan = {}
    for key, present in heavy_present.items():
        resname = resname_of[key]
        sidechain = _SIDECHAIN_HEAVY.get(resname)
        if sidechain is None:
            continue  # unknown/non-standard residue: leave as-is
        if not _BACKBONE.issubset(present):
            continue  # backbone itself incomplete: cannot safely stub
        if sidechain.issubset(present):
            continue  # complete residue
        if "CB" in present:
            plan[key] = ("ALA", _BACKBONE | {"CB"})
        else:
            plan[key] = ("GLY", set(_BACKBONE))
    return plan


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

    # First pass: collect the ATOM records that survive HETATM/TER/OXT filtering
    # so we can assess side-chain completeness per residue.
    with open(src_path) as fin:
        raw_lines = fin.readlines()

    atom_lines = []
    for line in raw_lines:
        record = line[:6].strip()
        if record != "ATOM":
            continue
        if line[12:16].strip() == "OXT":
            continue
        atom_lines.append(line)

    truncation_plan = _plan_residue_truncations(atom_lines)

    kept = 0
    truncated = []
    with open(out_path, "w") as fout:
        for line in raw_lines:
            record = line[:6].strip()
            if record in ("HETATM", "TER"):
                continue
            if record != "ATOM":
                fout.write(line)
                continue
            if line[12:16].strip() == "OXT":
                continue

            key = _residue_key(line)
            plan = truncation_plan.get(key)
            if plan is not None:
                stub_resname, keep_atoms = plan
                if line[12:16].strip() not in keep_atoms:
                    continue  # drop side-chain atoms beyond the stub
                # Rename the residue to the buildable stub (ALA/GLY).
                line = line[:17] + f"{stub_resname:>3}" + line[20:]
                if key not in truncated:
                    truncated.append(key)

            fout.write(line)
            kept += 1

    if truncated:
        summary = ", ".join(
            f"{chain}{resseq.strip()}{icode.strip()}"
            for (chain, resseq, icode) in truncated
        )
        print(
            f"sanitize_target_pdb: stubbed {len(truncated)} incomplete residue(s) "
            f"to ALA/GLY (missing side-chain atoms): {summary}"
        )

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
