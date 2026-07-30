#!/usr/bin/env python3
"""Recompute rmsd_md re-extracting against the CORRECT crystal peptide chain.

Bug fixed (STELLAR/calculate_md_rmsd.py::find_crystal_structure): the case suffix
(e.g. 2V8W_B) was treated as the PROTEIN chain, but in improve5Frag_2 it is the
PEPTIDE chain. When no file had that protein chain, the wrong crystal file was
picked, comparing the docked peptide against a reference peptide bound to a
DIFFERENT receptor copy -- tens of Angstrom away after alignment -- inflating
rmsd_md (2V8W_B: wrong chain F -> 44.6 A vs correct chain B -> 3.7 A).

This script, for every case/combo, reuses the already-saved sim receptor and sim
peptide (no MD re-run), rebuilds the sim complex, re-runs the SAME PyMOL
alignment/extraction against the CORRECT crystal (now selected by the fixed
find_crystal_structure), and recomputes the robust obrms in-place RMSD. It then
rewrites resultados_rmsd_md_<CASE>.csv and the rmsd_md column of
all_metrics_GN_<CASE>.csv.

Run on the HOST with the 3.9 shim (the function spawns Singularity for PyMOL):
  SINGULARITY_BIND=/data4 PYTHONNOUSERSITE=1 \
    agentic_loop_stellar/bin/python3 agentic_loop_stellar/recompute_realign.py <CASE ...>
"""
import csv
import os
import re
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "STELLAR"))
sys.path.insert(0, str(HERE))

import calculate_md_rmsd as C          # noqa: E402
from recompute_obrms import (find_peptide_dir, find_results_csv,  # noqa: E402
                             update_all_metrics)

SIF = "singularity/new_ms.simg"


def sim_combos(pdir):
    """combo -> (protein_sim, peptide_sim) paths that both exist."""
    out = {}
    for f in sorted(pdir.glob("protein_sim_*.pdb")):
        m = re.search(r"protein_sim_(\d+)\.pdb", f.name)
        if not m:
            continue
        combo = m.group(1)
        pep = pdir / f"peptide_sim_{combo}.pdb"
        if pep.is_file():
            out[combo] = (f, pep)
    return out


def build_sim_complex(protein_pdb, peptide_pdb, dest):
    with open(dest, "w") as out:
        for src in (protein_pdb, peptide_pdb):
            with open(src) as fh:
                for l in fh:
                    if l.startswith(("ATOM", "HETATM")):
                        out.write(l)
        out.write("END\n")


def candidate_crystals(case):
    """All crystal files whose PEPTIDE chain matches the case suffix. Crystals
    often contain several lattice copies of the same complex (same peptide chain,
    different protein chain), and only the copy the docking corresponds to aligns
    correctly -- so we try them all and keep the best (see pick-min in
    process_case). Falls back to every file for the PDB, then to
    find_crystal_structure."""
    if "_" not in case:
        return []
    base, suffix = case.split("_", 1)
    base = base.lower()
    suffix = suffix.upper()
    import glob
    all_files = sorted(glob.glob(str(ROOT / "complex" / ("%s_*.pdb" % base))) +
                       glob.glob(str(ROOT / "complex" / ("%s_*.pdb" % base.upper()))))

    def pep_chain(p):
        parts = os.path.basename(p).replace(".pdb", "").split("_")
        return parts[1].upper() if len(parts) >= 2 else None

    by_pep = [f for f in all_files if pep_chain(f) == suffix]
    cands = by_pep if by_pep else all_files
    if not cands:
        c, _ = C.find_crystal_structure("peptide_pdb_fragments", None, str(ROOT), case)
        if c:
            cands = [c if os.path.isabs(c) else str(ROOT / c)]
    # dedup preserving order
    seen, out = set(), []
    for c in cands:
        if c not in seen:
            seen.add(c)
            out.append(c)
    return out


def process_case(case, tmp):
    pdir = find_peptide_dir(case)
    if pdir is None:
        return None
    cands = candidate_crystals(case)
    if not cands:
        print("%-8s  no crystal found" % case)
        return None
    chain = case.split("_", 1)[1].upper() if "_" in case else "A"

    vsmap = {}
    res_csv = find_results_csv(case, "resultados_rmsd_md_%s.csv" % case)
    if res_csv:
        with open(res_csv) as fh:
            for row in csv.DictReader(fh):
                vsmap[row["combination"]] = row.get("vs_folder", "")

    combos = sim_combos(pdir)
    rows = []
    used_crystal = os.path.basename(cands[0])
    for combo, (prot, pep) in combos.items():
        cdir = os.path.join(tmp, "%s_%s" % (case, combo))
        os.makedirs(cdir, exist_ok=True)
        simc = os.path.join(cdir, "sim_complex.pdb")
        build_sim_complex(str(prot), str(pep), simc)
        best_v = None
        best_c = None
        for crystal in cands:
            psim = os.path.join(cdir, "psim.pdb")
            pcry = os.path.join(cdir, "pcry.pdb")
            prot_out = os.path.join(cdir, "prot.pdb")
            ok = C.extract_peptide_with_pymol(simc, crystal, psim, pcry, prot_out,
                                              SIF, chain, force_peptide_chain=True)
            if not ok:
                continue
            cbb = os.path.splitext(pcry)[0] + "_backbone.pdb"
            sbb = os.path.splitext(psim)[0] + "_backbone.pdb"
            if not (os.path.exists(cbb) and os.path.exists(sbb)):
                continue
            v = C.robust_backbone_rmsd(cbb, sbb, cdir, combo, SIF)
            if v is not None and (best_v is None or v < best_v):
                best_v = v
                best_c = os.path.basename(crystal)
        if best_c:
            used_crystal = best_c
        rows.append((combo, vsmap.get(combo, ""), best_v))
    return res_csv, rows, used_crystal, chain


def write_resultados(res_csv, case, rows):
    if res_csv is None:
        res_csv = ROOT / ("%s_results" % case) / "csv_files" / ("resultados_rmsd_md_%s.csv" % case)
        res_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(res_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["combination", "vs_folder", "rmsd"])
        for combo, vs, v in rows:
            w.writerow([combo, vs, "NA" if v is None else "%.4f" % v])


def main():
    cases = sys.argv[1:]
    grand_ok = grand_na = 0
    for case in cases:
        out = process_case(case, tempfile.mkdtemp(prefix="realign_"))
        if out is None:
            print("%-8s  SKIP" % case)
            continue
        res_csv, rows, cname, chain = out
        write_resultados(res_csv, case, rows)
        update_all_metrics(case, {c: v for c, vs, v in rows})
        nok = sum(1 for _, _, v in rows if v is not None)
        nna = sum(1 for _, _, v in rows if v is None)
        grand_ok += nok
        grand_na += nna
        detail = ", ".join("%s=%s" % (c, "NA" if v is None else "%.3f" % v)
                           for c, _, v in rows)
        print("%-8s  crystal=%-14s chain=%s  obrms=%d na=%d  [%s]"
              % (case, cname, chain, nok, nna, detail))
    print("\n=== RESUMEN realineado ===")
    print("poses: obrms=%d NA=%d" % (grand_ok, grand_na))


if __name__ == "__main__":
    main()
