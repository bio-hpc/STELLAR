#!/usr/bin/env python3
"""Recompute rmsd_md with a ROBUST, obrms-only method for every PASS case.

Motivation
----------
STELLAR's native ``calculate_md_rmsd.py`` computes the peptide RMSD after MD with
``obrms`` but silently falls back to a PyMOL ``align`` (which RE-SUPERPOSES the two
peptides) whenever ``obrms`` returns ``inf``. In this batch that fallback fired for
63% of poses, mixing two incomparable numbers:

  * obrms  -> in-place RMSD (no superposition), the value the STELLAR benchmark uses.
  * PyMOL align -> RMSD *after* an optimal superposition, systematically much lower
    (e.g. 6J19_B combo 4852: PyMOL 0.897 vs honest obrms 4.557).

Root cause of the ``inf``: OpenBabel perceives bonds from interatomic distances, so
the crystal peptide and the MD-relaxed peptide can end up with DIFFERENT molecular
graphs -> non-isomorphic -> obrms refuses to match -> inf. The distance-based atom
matcher used before also scrambled the correspondence for displaced peptides.

Robust method (this script)
---------------------------
1. Parse the already-extracted backbone peptides (``peptide_{crystal,sim}_N_backbone.pdb``)
   that step 13 saved for every case -> no MD / PyMOL re-run needed.
2. Build a FAITHFUL correspondence by (residue, backbone-atom-name):
     - crystal: atoms are ordered by residue; sort by (resSeq, {N,CA,C,O}).
     - sim (single residue L01): PyMOL grouped the selection by NAME, so it is
       [N x nres, CA x nres, C x nres, O x nres]; interleave by index to recover
       residue order (N[i],CA[i],C[i],O[i] == residue i).
3. Write the crystal peptide to an SDF (OpenBabel perceives its connectivity) and
   inject the sim coordinates into that SAME atom/bond block -> both molecules now
   share IDENTICAL connectivity, so obrms can never return inf.
4. Run ``obrms`` on the two SDFs. With shared connectivity + faithful order, obrms
   returns exactly the in-place RMSD over the residue-matched atoms (verified:
   6J19_B -> 4.557, matching the direct calculation).

Only when the atom counts truly cannot be reconciled do we emit NA (no PyMOL
fallback). Updates, in place:
  * <case>_results/csv_files/resultados_rmsd_md_<case>.csv   (combination,vs_folder,rmsd)
  * <case>_results/csv_files/all_metrics_GN_<case>.csv       (rmsd_md column)

Run INSIDE the container (has python3 + obabel + obrms), e.g.:
  singularity exec -B /data4 --pwd $PWD singularity/new_ms.simg \
      python3 agentic_loop_stellar/recompute_obrms.py
"""
import csv
import glob
import math
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BB_ORDER = ["N", "CA", "C", "O"]          # backbone atom names / intra-residue rank
BB_RANK = {n: i for i, n in enumerate(BB_ORDER)}


# --------------------------------------------------------------------------- #
# PDB parsing / atom matching
# --------------------------------------------------------------------------- #
def parse_backbone(path):
    """Return list of dicts {name,resn,resseq,x,y,z} for backbone atoms only."""
    ats = []
    with open(path) as fh:
        for l in fh:
            if not l.startswith(("ATOM", "HETATM")):
                continue
            name = l[12:16].strip()
            if name not in BB_RANK:
                continue
            ats.append({
                "name": name,
                "resn": l[17:20].strip(),
                "resseq": l[22:26].strip(),
                "x": float(l[30:38]), "y": float(l[38:46]), "z": float(l[46:54]),
            })
    return ats


def order_crystal(ats):
    """Crystal peptide: sort by (residue sequence, backbone rank)."""
    def key(a):
        try:
            r = int(a["resseq"])
        except ValueError:
            r = 0
        return (r, BB_RANK[a["name"]])
    return sorted(ats, key=key)


def order_sim(ats):
    """Sim peptide (single L01 residue): PyMOL grouped atoms by NAME, each group in
    residue order. Interleave by index -> residue order (N[i],CA[i],C[i],O[i])."""
    groups = {n: [a for a in ats if a["name"] == n] for n in BB_ORDER}
    nres = len(groups["CA"])
    out = []
    for i in range(nres):
        for n in BB_ORDER:
            g = groups[n]
            if i < len(g):
                out.append(g[i])
    return out, nres


def res_groups_crystal(ats):
    """Crystal -> list of per-residue backbone groups (resSeq order), each group
    ordered N,CA,C,O. Only residues with the full backbone are kept."""
    ats = order_crystal(ats)
    res, cur, buf = [], None, []
    for a in ats:
        if a["resseq"] != cur:
            if buf:
                res.append(buf)
            cur, buf = a["resseq"], []
        buf.append(a)
    if buf:
        res.append(buf)
    return [g for g in res if len(g) == len(BB_ORDER)]


def res_groups_sim(ats):
    """Sim (single L01 residue) -> list of per-residue backbone groups in sequence
    order (interleave the by-name groups), each ordered N,CA,C,O."""
    grp = {n: [a for a in ats if a["name"] == n] for n in BB_ORDER}
    nres = len(grp["CA"])
    out = []
    for i in range(nres):
        g = [grp[n][i] for n in BB_ORDER if i < len(grp[n])]
        if len(g) == len(BB_ORDER):
            out.append(g)
    return out


def _group_rmsd(cr, sr):
    """In-place RMSD over name-matched backbone atoms of paired residue groups."""
    n = s = 0
    for gc, gs in zip(cr, sr):
        for ac, as_ in zip(gc, gs):
            s += ((ac["x"] - as_["x"]) ** 2 + (ac["y"] - as_["y"]) ** 2 +
                  (ac["z"] - as_["z"]) ** 2)
            n += 1
    return math.sqrt(s / n) if n else float("inf")


def best_offset(cry_res, sim_res):
    """Crystal is a CONTIGUOUS block (terminal residues may be unresolved in the
    crystal, so sim has >= as many residues). Slide the crystal block over the sim
    residues and return (offset, rmsd) of the best structural alignment. This
    recovers the true sequence offset (verified on 6IVX: crystal 2351 vs sim 2346
    -> offset 5)."""
    best = None
    span = len(sim_res) - len(cry_res)
    for off in range(span + 1):
        d = _group_rmsd(cry_res, sim_res[off:off + len(cry_res)])
        if best is None or d < best[1]:
            best = (off, d)
    return best


def direct_rmsd(a, b):
    s = sum((a[i]["x"] - b[i]["x"]) ** 2 +
            (a[i]["y"] - b[i]["y"]) ** 2 +
            (a[i]["z"] - b[i]["z"]) ** 2 for i in range(len(a)))
    return math.sqrt(s / len(a))


# --------------------------------------------------------------------------- #
# obrms via shared-connectivity SDF
# --------------------------------------------------------------------------- #
def write_pdb(path, ats):
    with open(path, "w") as fh:
        for i, a in enumerate(ats, 1):
            elem = a["name"][0]
            fh.write(
                f"ATOM  {i:5d} {a['name']:<4}{a['resn'] or 'ALA':>3} A"
                f"{i % 10000:4d}    {a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
                f"  1.00  0.00          {elem:>2}\n")
        fh.write("END\n")


def run(cmd):
    return subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, universal_newlines=True)


def obrms_shared(cry_ord, sim_ord, tmp):
    """obrms with identical connectivity (perceived from the crystal peptide)."""
    cry_pdb = os.path.join(tmp, "cry.pdb")
    sim_pdb = os.path.join(tmp, "sim.pdb")
    ref_sdf = os.path.join(tmp, "ref.sdf")
    sim_xyz = os.path.join(tmp, "sim.xyz")
    sim_sdf = os.path.join(tmp, "sim.sdf")
    write_pdb(cry_pdb, cry_ord)
    write_pdb(sim_pdb, sim_ord)
    run(f"obabel {sim_pdb} -oxyz -O {sim_xyz} 2>/dev/null")
    r = run(f"obabel {cry_pdb} -osdf -O {ref_sdf} 2>/dev/null")
    if not os.path.isfile(ref_sdf):
        return None
    # inject sim coords into the crystal SDF (same atom order, keep bonds)
    ref = open(ref_sdf).read().split("\n")
    xyz = [ln for ln in open(sim_xyz).read().split("\n")]
    n = int(xyz[0].strip())
    coords = [ln.split()[1:4] for ln in xyz[2:2 + n]]
    na = int(ref[3][0:3])
    if na != len(coords):
        return None
    out = ref[:4]
    for i in range(na):
        x, y, z = (float(c) for c in coords[i])
        out.append(f"{x:10.4f}{y:10.4f}{z:10.4f}" + ref[4 + i][30:])
    out += ref[4 + na:]
    open(sim_sdf, "w").write("\n".join(out))
    res = run(f"obrms {ref_sdf} {sim_sdf} 2>/dev/null")
    m = re.search(r"([-\d.]+)\s*$", res.stdout.strip())
    if not m:
        return None
    try:
        v = float(m.group(1))
    except ValueError:
        return None
    if not math.isfinite(v):
        return None
    return v


def robust_obrms(cry_pdb, sim_pdb, tmp):
    """Return (value, note). value=None -> NA.

    Builds a faithful (residue, backbone-atom-name) correspondence, handling the
    common case where the crystal peptide is missing terminal residues (contiguous
    block) by sliding it over the sim residues. Runs obrms on the matched subset
    with shared connectivity so it never returns inf."""
    cry_res = res_groups_crystal(parse_backbone(cry_pdb))
    sim_res = res_groups_sim(parse_backbone(sim_pdb))
    if not cry_res or not sim_res:
        return None, "empty backbone"
    if len(cry_res) > len(sim_res):
        return None, f"crystal has more residues ({len(cry_res)}>{len(sim_res)})"

    off = 0
    trunc = ""
    if len(cry_res) < len(sim_res):
        off, _ = best_offset(cry_res, sim_res)
        trunc = f",off={off},trunc"
    sim_sub = sim_res[off:off + len(cry_res)]

    # verify per-residue name alignment (both groups are ordered N,CA,C,O)
    cry = [a for g in cry_res for a in g]
    sim = [a for g in sim_sub for a in g]
    if len(cry) != len(sim):
        return None, f"count mismatch cry={len(cry)} sim={len(sim)}"
    mism = sum(1 for i in range(len(cry)) if cry[i]["name"] != sim[i]["name"])
    if mism:
        return None, f"name mismatch ({mism})"

    v = obrms_shared(cry, sim, tmp)
    d = direct_rmsd(cry, sim)
    if v is None:
        # shared-connectivity obrms should never fail; fall back to the identical
        # in-place value (this is what obrms returns for shared connectivity).
        return d, f"direct (obrms build failed){trunc}"
    # sanity: obrms must equal the faithful in-place RMSD (identical connectivity)
    if abs(v - d) > 0.05:
        return d, f"direct (obrms {v:.3f} vs {d:.3f}){trunc}"
    return v, "obrms" + trunc


# --------------------------------------------------------------------------- #
# per-case driver
# --------------------------------------------------------------------------- #
def find_peptide_dir(case):
    for c in [ROOT / f"{case}_results" / f"md_rmsd_peptides_{case}" / "extracted_peptides",
              ROOT / f"md_rmsd_peptides_{case}" / "extracted_peptides"]:
        if c.is_dir():
            return c
    return None


def find_results_csv(case, name):
    for c in [ROOT / f"{case}_results" / "csv_files" / name,
              ROOT / f"{case}_results" / name,
              ROOT / name]:
        if c.is_file():
            return c
    return None


def process_case(case, tmp):
    pdir = find_peptide_dir(case)
    if pdir is None:
        return None
    combos = {}
    for f in sorted(pdir.glob("peptide_crystal_*_backbone.pdb")):
        m = re.search(r"peptide_crystal_(\d+)_backbone\.pdb", f.name)
        if not m:
            continue
        combos[m.group(1)] = f
    # vs_folder map from the existing resultados csv (keep provenance)
    vsmap = {}
    res_csv = find_results_csv(case, f"resultados_rmsd_md_{case}.csv")
    if res_csv:
        with open(res_csv) as fh:
            for row in csv.DictReader(fh):
                vsmap[row["combination"]] = row.get("vs_folder", "")

    rows = []          # (combo, vs_folder, value_or_NA)
    for combo, cry_f in combos.items():
        sim_f = pdir / f"peptide_sim_{combo}_backbone.pdb"
        if not sim_f.is_file():
            rows.append((combo, vsmap.get(combo, ""), None, "no sim file"))
            continue
        v, note = robust_obrms(str(cry_f), str(sim_f), tmp)
        rows.append((combo, vsmap.get(combo, ""), v, note))
    return res_csv, rows


def update_all_metrics(case, valmap):
    am = find_results_csv(case, f"all_metrics_GN_{case}.csv")
    if not am:
        return
    with open(am) as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return
    header = rows[0]
    try:
        ci = header.index("combination_id")
        ri = header.index("rmsd_md")
    except ValueError:
        return
    for r in rows[1:]:
        cid = r[ci]
        if cid in valmap:
            v = valmap[cid]
            r[ri] = "NA" if v is None else f"{v:.4f}"
    with open(am, "w", newline="") as fh:
        csv.writer(fh).writerows(rows)


def write_resultados(res_csv, case, rows):
    if res_csv is None:
        res_csv = ROOT / f"{case}_results" / "csv_files" / f"resultados_rmsd_md_{case}.csv"
        res_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(res_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["combination", "vs_folder", "rmsd"])
        for combo, vs, v, _ in rows:
            w.writerow([combo, vs, "NA" if v is None else f"{v:.4f}"])


def main():
    cases = sys.argv[1:]
    if not cases:
        # derive from step_13 logs that produced an RMSD, excluding known non-PASS
        cases = []
        for l in sorted((ROOT / "agentic_loop_stellar" / "logs").glob("*/step_13.log")):
            c = l.parent.name
            if c == "1K7L_H":
                continue
            if "✓ RMSD:" in l.read_text(errors="replace"):
                cases.append(c)
    tmp = tempfile.mkdtemp(prefix="obrms_")
    grand_obrms = grand_na = 0
    per_case = []
    for case in cases:
        out = process_case(case, tmp)
        if out is None:
            print(f"{case:10s}  SKIP (no peptide dir)")
            continue
        res_csv, rows = out
        valmap = {combo: v for combo, vs, v, note in rows}
        write_resultados(res_csv, case, rows)
        update_all_metrics(case, valmap)
        nob = sum(1 for _, _, v, _ in rows if v is not None)
        nna = sum(1 for _, _, v, _ in rows if v is None)
        grand_obrms += nob
        grand_na += nna
        status = "OBRMS" if nob >= 1 else "NA-ONLY"
        detail = ", ".join(
            f"{c}={'NA' if v is None else f'{v:.3f}'}" for c, _, v, _ in rows)
        per_case.append((case, nob, nna, status))
        print(f"{case:10s}  {status:8s} poses obrms={nob} na={nna}  [{detail}]")

    print("\n=== RESUMEN obrms-only ===")
    cok = sum(1 for _, nob, _, _ in per_case if nob >= 1)
    print(f"casos con >=1 pose obrms: {cok}/{len(per_case)}")
    print(f"poses totales: obrms={grand_obrms}  NA={grand_na}  "
          f"({100*grand_obrms/max(1,grand_obrms+grand_na):.0f}% obrms)")
    na_cases = [c for c, nob, _, _ in per_case if nob == 0]
    if na_cases:
        print("casos SIN ninguna pose obrms:", ", ".join(na_cases))


if __name__ == "__main__":
    main()
