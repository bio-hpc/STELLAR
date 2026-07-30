#!/usr/bin/env python3
"""
Agentic-loop driver for a single STELLAR case (improve5Frag_2).

Runs the GN pipeline step by step with per-step observability so the loop can
detect exactly where a case fails, feed that back to the diagnosis layer, apply
a targeted fix, and retry.

Pilot-friendly: reduces the number of reconstructed combinations that reach the
(expensive) GROMACS MD to just a couple, instead of the usual ~100.

Success (as defined by the user):
    all_metrics_GN_<CASE>.csv exists AND at least one row has rmsd_md != NA.

Run with Python >= 3.9 (see agentic_loop_stellar/bin/python3 shim).
"""
import argparse
import csv
import os
import re
import shlex
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                      # project root (has STELLAR/, GROMACS/, ...)
SHIM_BIN = HERE / "bin"                 # provides python3 -> python3.9

SIF = ROOT / "singularity" / "STELLAR.simg"   # Ubuntu 22.04: glibc 2.35, full deps
# Top-level mount so the container sees the project path (ROOT is under /data4).
BIND = "/" + str(ROOT).strip("/").split("/")[0]

# Steps that are NOT required for the success criterion (rmsd_md). If they fail
# we log it and keep going so merge_all_metrics can still emit rmsd_md.
NON_FATAL_STEPS = {"0", "15", "15b", "15c", "17"}
# Steps skipped by default (docking already done for these cases; calcium opt).
SKIP_STEPS = {"0", "9"}
# Steps that rely on compiled binaries (filter_fragment_combinations, overlap)
# needing a newer glibc than the CentOS-7 host -> run inside STELLAR.simg.
# These scripts use only the stdlib + the binary, so containerizing is safe
# (no nested singularity).
CONTAINER_STEPS = {"2", "4"}


def containerize(cmd):
    """Wrap a shell command so it runs inside STELLAR.simg at the project root."""
    inner = f"cd {shlex.quote(str(ROOT))} && {cmd}"
    return (
        f"singularity exec -B {BIND} --pwd {shlex.quote(str(ROOT))} "
        f"{shlex.quote(str(SIF))} bash -c {shlex.quote(inner)}"
    )


def log(msg, fh=None):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    if fh:
        fh.write(line + "\n")
        fh.flush()


def prep_environment(master_fh):
    """One-time, idempotent environment adaptations for this (non-BSC) cluster."""
    # The GROMACS topology generator (run inside gr.simg) hardcodes the legacy
    # "ShuttleMol/" directory name (config.cfg, external_sw/, extra_shuttlemol/).
    # This repo ships it as GROMACS/, so expose it under the expected name.
    link = ROOT / "ShuttleMol"
    if not link.exists():
        try:
            link.symlink_to("GROMACS")
            log("created symlink ShuttleMol -> GROMACS", master_fh)
        except OSError as e:
            log(f"WARN: could not create ShuttleMol symlink: {e}", master_fh)


def run_shim_env():
    env = dict(os.environ)
    env["PATH"] = f"{SHIM_BIN}:{env.get('PATH','')}"
    # Make the project path visible inside every Singularity container the
    # pipeline scripts spawn internally (host is CentOS 7; /data4 is not
    # auto-mounted). Without this, containers run from $HOME and relative paths
    # like MetaScreener/.../convert_to.py fail to resolve.
    for var in ("SINGULARITY_BIND", "APPTAINER_BIND"):
        existing = env.get(var, "")
        env[var] = f"{existing},{BIND}" if existing else BIND
    # Stop the host's ~/.local (Python 2.7/3 user site-packages, built for the
    # CentOS-7 host) from leaking into containers via the auto-mounted $HOME and
    # shadowing the container's own numpy/etc. (breaks MGLTools, RDKit, ...).
    env["SINGULARITYENV_PYTHONNOUSERSITE"] = "1"
    env["APPTAINERENV_PYTHONNOUSERSITE"] = "1"
    env["PYTHONNOUSERSITE"] = "1"
    return env


def parse_steps(commands_file):
    """Parse the generated commands file into an ordered list of (step_id, [cmds])."""
    steps = []
    cur_id = None
    cur_cmds = []
    pending = ""  # line-continuation accumulator
    hdr = re.compile(r"^#\s*STEP\s+([0-9]+[a-z]?):", re.IGNORECASE)
    with open(commands_file) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            m = hdr.match(line)
            if m:
                if cur_id is not None:
                    steps.append((cur_id, cur_cmds))
                cur_id = m.group(1)
                cur_cmds = []
                pending = ""
                continue
            if cur_id is None:
                continue
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if stripped.endswith("\\"):
                pending += stripped[:-1] + " "
                continue
            cmd = (pending + stripped).strip()
            pending = ""
            cur_cmds.append(cmd)
    if cur_id is not None:
        steps.append((cur_id, cur_cmds))
    return steps


def apply_pilot_overrides(step_id, cmd, case, max_process, pool, tolerance, md_target):
    """Rewrite commands for the pilot (fewer combinations).

    Overlap filter (step 4): use the production tolerance, process the whole
    (already sampled+capped) organized set, and keep the `pool` lowest-overlap
    combinations as MD candidates.

    MD (step 11) — Group A retry-with-other-combinations:
      The reconstructed peptide of some cases yields GROMACS zero-length
      constraints (assertion ip.constr.dA > 0) for certain fragment geometries.
      Instead of requiring the (few) lowest-overlap combos to all pass MD, we
      carry a larger `pool` of candidates through steps 4-10 and let the MD
      launcher walk them SEQUENTIALLY, stopping once `md_target` of them finish
      successfully (`--max-jobs` counts successes only, so failed combos don't
      consume the budget). `--check-existing` makes it idempotent on resume.
    """
    if step_id == "4" and "check_overlap_combinations.py" in cmd:
        # ... <indir> GN <tol> <max_process> <pool>
        cmd = re.sub(
            r"(check_overlap_combinations\.py\s+\S+\s+GN)\s+\d+\s+\d+\s+\d+",
            rf"\1 {tolerance} {max_process} {pool}",
            cmd,
        )
    if step_id == "11" and "run_md_simulations.py" in cmd:
        if "--max-jobs" not in cmd:
            cmd += f" --max-jobs {md_target}"
        if "--check-existing" not in cmd:
            cmd += " --check-existing"
    return cmd


def count_completed_md(case):
    """Number of MD runs that produced a *_complex_md.pdb for this case.

    MD outputs land in VS_GR_<pdb>_<recchain>_GN_query_combination_*_GPU_<date>/
    at the project root; a completed run writes molecules/*_complex_md.pdb.
    """
    pdb = case.split("_")[0].lower()
    n = 0
    for d in ROOT.glob(f"VS_GR_{pdb}_*_GN_*combination_*"):
        mol = d / "molecules"
        if mol.is_dir() and any(mol.glob("*_complex_md.pdb")):
            n += 1
    return n


def remove_empty_files(base, pattern, master_fh):
    """Delete 0-byte files matching pattern under base (recursively)."""
    if not Path(base).is_dir():
        return
    removed = 0
    for p in Path(base).rglob(pattern):
        try:
            if p.is_file() and p.stat().st_size == 0:
                p.unlink()
                removed += 1
        except OSError:
            pass
    if removed:
        log(f"removed {removed} empty {pattern} under {base}", master_fh)


# Standard bond lengths to a hydrogen (nm), keyed by the heavy-atom element.
# Used to repair ligand bonds that antechamber/GAFF left at 0 (see below).
_H_BOND_R0_NM = {"S": 0.13360, "O": 0.09740, "N": 0.10100,
                 "C": 0.10900, "P": 0.14140}
_H_BOND_K = 3.1400e05  # kJ/mol/nm^2 (ignored while the bond is a LINCS constraint)

# Atomic masses (g/mol) by element, to repair ligand atoms antechamber left as
# massless dummies (atomtype `DU`, mass 0) — see repair_ligand_topologies.
_ELEM_MASS = {"H": 1.00800, "C": 12.01000, "N": 14.01000, "O": 16.00000,
              "S": 32.06000, "P": 30.97000, "F": 19.00000, "CL": 35.45000,
              "BR": 79.90000, "I": 126.90000}


def _atom_element(name):
    """Element guess from a PDB/GAFF atom name (e.g. 'HG'->H, 'SG'->S, '1HB'->H)."""
    if not name:
        return "C"
    c = name[0]
    if c.isdigit() and len(name) > 1:
        c = name[1]
    return c.upper()


def repair_ligand_topologies(case, master_fh):
    """Fix zero-length H bonds in the reconstructed-peptide (ligand) topology.

    Root cause (Group A, systematic): for peptides containing a free cysteine
    thiol, the ligand parameterization run in step 10 fails to assign the S-H
    (SG-HG) bond and writes it as `0.0000e+00 0.0000e+00`. Because the MD uses
    `constraints = h-bonds`, that bond becomes a LINCS constraint of length 0
    and GROMACS aborts at NVT grompp:

        Assertion failed: ip.constr.dA > 0
        We should only have positive constraint lengths here

    This is NOT geometry-dependent, so retrying other combinations never helps.
    We repair it by giving every zero-valued bond-to-hydrogen the standard bond
    length for its heavy atom (only b0 matters once it is a constraint).

    Second defect (also systematic, Group A): for certain reconstructed side
    chains (e.g. a terminal alkyne carbon) antechamber fails to assign a GAFF
    atom type and emits the atom as a massless dummy (atomtype `DU`, mass 0)
    that is left un-bonded to the rest of the molecule.  grompp then aborts at
    the ions/NVT step with:

        ERROR: atom <name> (Res L01-1) has mass 0 (state A) / 0 (state B)

    Since the atom carries no LJ (DU sigma/epsilon = 0) and no bonded terms,
    we make it an inert tracer: give it the real element mass (so grompp is
    happy) and zero its partial charge (so a non-interacting massless-turned-
    massive point can't collapse onto a nearby charge and blow up the MD).  It
    does not participate in any interaction, so it cannot affect the backbone
    RMSD that defines success.

    Neither defect is geometry-dependent, so retrying other combinations never
    helps.  Runs after step 10 and before step 11; edits only the staged working
    copy under valid_GN_<CASE>_final/ (never the original case data).
    """
    base = ROOT / f"valid_GN_{case}_final"
    if not base.is_dir():
        return
    fixed_bonds = 0
    fixed_atoms = 0
    fixed_files = 0
    for itp in base.rglob("*.itp"):
        lines = itp.read_text(errors="replace").splitlines()
        # Map atom index -> element from the [ atoms ] section.
        elem = {}
        sec = None
        for ln in lines:
            s = ln.strip()
            if s.startswith("["):
                sec = s.strip("[] ").lower()
                continue
            if sec == "atoms" and s and not s.startswith(";"):
                p = s.split()
                if len(p) >= 5 and p[0].isdigit():
                    elem[int(p[0])] = _atom_element(p[4])
        # Repair zero-valued bonds that involve a hydrogen.
        out = []
        sec = None
        changed = False
        for ln in lines:
            s = ln.strip()
            if s.startswith("["):
                sec = s.strip("[] ").lower()
                out.append(ln)
                continue
            if sec == "atoms" and s and not s.startswith(";"):
                p = s.split()
                # nr type resi res atom cgnr charge mass [; comment]
                if len(p) >= 8 and p[0].isdigit():
                    try:
                        mass = float(p[7])
                    except ValueError:
                        out.append(ln)
                        continue
                    if mass <= 1e-6:
                        m = _ELEM_MASS.get(_atom_element(p[4]), 12.01000)
                        comment = (" ;" + ln.split(";", 1)[1]) if ";" in ln else ""
                        # Keep the (inert, zero-LJ) atomtype; zero the charge so
                        # this un-bonded point can't collapse onto a charge, and
                        # give it the real element mass so grompp accepts it.
                        out.append(
                            f"{p[0]:>6} {p[1]:>4} {p[2]:>5} {p[3]:>5} "
                            f"{p[4]:>5} {p[5]:>5} {0.0:>12.6f} {m:>12.5f}{comment}")
                        changed = True
                        fixed_atoms += 1
                        continue
            if sec == "bonds" and s and not s.startswith(";"):
                p = s.split()
                if len(p) >= 5 and p[0].isdigit():
                    ai, aj = int(p[0]), int(p[1])
                    try:
                        r, k = float(p[3]), float(p[4])
                    except ValueError:
                        out.append(ln)
                        continue
                    if r <= 1e-6:
                        ei, ej = elem.get(ai, ""), elem.get(aj, "")
                        if "H" in (ei, ej):
                            heavy = ej if ei == "H" else ei
                            r0 = _H_BOND_R0_NM.get(heavy, 0.10900)
                            comment = (" ;" + ln.split(";", 1)[1]) if ";" in ln else ""
                            out.append(
                                f"{ai:>7}{aj:>7}   1    "
                                f"{r0:.4e}    {_H_BOND_K:.4e}{comment}")
                            changed = True
                            fixed_bonds += 1
                            continue
            out.append(ln)
        if changed:
            itp.write_text("\n".join(out) + "\n")
            fixed_files += 1
    if fixed_bonds or fixed_atoms:
        log(f"repaired {fixed_bonds} zero-length H-bond(s) and {fixed_atoms} "
            f"zero-mass dummy atom(s) in {fixed_files} ligand topology file(s)",
            master_fh)


def truncate_combinations_csv(case, n, master_fh):
    """Keep an evenly spread sample of N combinations before organizing.

    Huge speed/disk win for the pilot: organize_valid_combinations would
    otherwise materialize *every* valid combination (often >10k) when we only
    need a couple to reach MD.

    Spread sampling (not top-by-score): the best-scoring combinations tend to
    cluster spatially and clash (all rejected by the overlap filter), so we
    sample across the whole score distribution to include non-overlapping ones.
    """
    path = ROOT / f"{case}_GN" / "valid_fragment_combinations.csv"
    if not path.is_file():
        return
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    if len(rows) <= n + 1:
        return
    header, data = rows[0], rows[1:]
    step = max(1, len(data) // n)
    kept = data[::step][:n]
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(kept)
    log(f"sampled combinations CSV: {len(data)} -> {len(kept)} (every {step}th)", master_fh)


def success(case):
    """True if all_metrics_GN_<CASE>.csv has >=1 row with rmsd_md != NA.

    Step 17 (organize_workflow_results) relocates the CSV into
    <CASE>_results/csv_files/, so check both locations.
    """
    candidates = [
        ROOT / f"all_metrics_GN_{case}.csv",
        ROOT / f"{case}_results" / "csv_files" / f"all_metrics_GN_{case}.csv",
        ROOT / f"{case}_results" / f"all_metrics_GN_{case}.csv",
    ]
    csv_path = next((p for p in candidates if p.is_file()), None)
    if csv_path is None:
        return False, "no all_metrics CSV"
    try:
        with open(csv_path, newline="") as fh:
            reader = csv.DictReader(fh)
            col = None
            for c in reader.fieldnames or []:
                if c.strip().lower() in ("rmsd_md", "md_rmsd", "rmsd_md_backbone"):
                    col = c
                    break
            if col is None:
                return False, f"no rmsd_md column (cols={reader.fieldnames})"
            n_ok = 0
            for row in reader:
                v = (row.get(col) or "").strip()
                if v and v.upper() != "NA" and v.lower() != "nan":
                    n_ok += 1
            if n_ok > 0:
                return True, f"{n_ok} rows with valid {col}"
            return False, f"all {col} are NA"
    except Exception as e:  # noqa: BLE001
        return False, f"error reading CSV: {e}"


def _exec(cmd, sfh, env):
    t0 = time.time()
    proc = subprocess.run(
        cmd, shell=True, cwd=str(ROOT), env=env,
        stdout=sfh, stderr=subprocess.STDOUT,
    )
    return proc.returncode, time.time() - t0


def run_step(step_id, cmds, case, logdir, master_fh, env):
    """Run all commands of a step. Return (ok, failing_cmd, log_path).

    Applies the diagnosis/feedback rules:
      - Steps with compiled binaries run inside STELLAR.simg.
      - Any command that fails with a GLIBC version error is retried once
        inside STELLAR.simg.
    """
    step_log = logdir / f"step_{step_id}.log"
    with open(step_log, "w") as sfh:
        for raw_cmd in cmds:
            cmd = containerize(raw_cmd) if step_id in CONTAINER_STEPS else raw_cmd
            log(f"STEP {step_id}: $ {cmd}", master_fh)
            sfh.write(f"$ {cmd}\n")
            sfh.flush()
            rc, dt = _exec(cmd, sfh, env)
            sfh.write(f"[exit {rc} in {dt:.1f}s]\n")
            sfh.flush()
            if rc != 0:
                # Feedback rule: GLIBC mismatch -> retry inside container.
                logtxt = tail(step_log, 60)
                if "GLIBC" in logtxt and "singularity exec" not in cmd:
                    fixed = containerize(raw_cmd)
                    log(f"STEP {step_id}: GLIBC error -> retry in container", master_fh)
                    sfh.write(f"# DIAGNOSE: GLIBC -> retry containerized\n$ {fixed}\n")
                    sfh.flush()
                    rc, dt = _exec(fixed, sfh, env)
                    sfh.write(f"[exit {rc} in {dt:.1f}s]\n")
                    sfh.flush()
                if rc != 0:
                    return False, raw_cmd, step_log
    return True, None, step_log


def tail(path, n=25):
    try:
        lines = Path(path).read_text(errors="replace").splitlines()
        return "\n".join(lines[-n:])
    except Exception:  # noqa: BLE001
        return "(no log)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--max-valid", type=int, default=2,
                    help="target number of combinations that must COMPLETE MD "
                         "(step 11 stops after this many successes)")
    ap.add_argument("--md-pool", type=int, default=0,
                    help="Group A: candidate combinations carried through steps "
                         "4-10 so MD can retry other geometries when some hit "
                         "GROMACS zero-length constraints (0 = same as max-valid)")
    ap.add_argument("--max-process", type=int, default=0,
                    help="combinations processed by the overlap filter (0=all organized)")
    ap.add_argument("--max-organize", type=int, default=400,
                    help="spread-sample cap before step 3 (pilot speed/disk)")
    ap.add_argument("--overlap-tol", type=int, default=450,
                    help="overlap tolerance for step 4 (production value)")
    ap.add_argument("--from-step", default=None,
                    help="resume from this step id (skip earlier)")
    ap.add_argument("--stop-on-fail", action="store_true", default=True)
    args = ap.parse_args()
    case = args.case

    logdir = HERE / "logs" / case
    logdir.mkdir(parents=True, exist_ok=True)
    master = logdir / "run.log"
    env = run_shim_env()

    with open(master, "a") as mfh:
        log(f"=== run_case {case} (max_valid={args.max_valid}) ===", mfh)
        prep_environment(mfh)

        # 1) stage (write-isolated). Only (re)stage when we will run the early
        # steps that consume the staged inputs; re-staging wipes <CASE>_GN and
        # the intermediate CSVs stored there, so skip it when resuming later.
        need_stage = args.from_step is None or args.from_step in ("1", "2")
        if need_stage:
            r = subprocess.run(["bash", str(HERE / "stage_case.sh"), case],
                               cwd=str(ROOT), env=env)
            if r.returncode != 0:
                log(f"STAGE FAILED for {case}", mfh)
                return 3
        else:
            log(f"skip staging (resume from step {args.from_step})", mfh)

        # 2) generate commands
        cmds_file = logdir / f"commands_{case}.txt"
        r = subprocess.run(
            ["python3", "generate_workflow_commands.py", case, str(cmds_file)],
            cwd=str(ROOT), env=env,
            stdout=open(logdir / "generate.log", "w"), stderr=subprocess.STDOUT,
        )
        if r.returncode != 0 or not cmds_file.is_file():
            log(f"GENERATE FAILED for {case}", mfh)
            return 4

        steps = parse_steps(cmds_file)
        log(f"parsed {len(steps)} steps: {[s for s,_ in steps]}", mfh)

        started = args.from_step is None
        for step_id, cmds in steps:
            if not started:
                if step_id == args.from_step:
                    started = True
                else:
                    continue
            if step_id in SKIP_STEPS:
                log(f"STEP {step_id}: skipped (policy)", mfh)
                continue
            # Pilot: cap combinations right before organizing them.
            if step_id == "3" and args.max_organize:
                truncate_combinations_csv(case, args.max_organize, mfh)
            # Robustness: steps skip-if-exists on filename only. Remove 0-byte
            # artifacts left by a previous failed conversion so they regenerate.
            if step_id in ("5", "6"):
                remove_empty_files(
                    ROOT / f"valid_combinations_GN_{case}" / "valid_no_overlap",
                    "*.mol2", mfh)
            # Group A root-cause fix: repair zero-length S-H (and other H) bonds
            # the ligand parameterizer left in the step-10 topology, which
            # otherwise crash NVT grompp with a zero-length LINCS constraint.
            if step_id == "11":
                repair_ligand_topologies(case, mfh)
            pool = args.md_pool if args.md_pool else args.max_valid
            cmds = [apply_pilot_overrides(step_id, c, case,
                                          args.max_process, pool,
                                          args.overlap_tol, args.max_valid)
                    for c in cmds]
            ok, bad, slog = run_step(step_id, cmds, case, logdir, mfh, env)
            if not ok and step_id == "11":
                # Group A: run_md_simulations exits 1 if ANY combination fails,
                # but the pipeline only needs >=1 completed MD to produce a
                # non-NA rmsd_md (step 13 skips combos without MD output). Treat
                # the step as OK when at least one combination finished.
                n_md = count_completed_md(case)
                if n_md >= 1:
                    log(f"STEP 11: {n_md} combination(s) completed MD despite "
                        f"launcher errors -> continue", mfh)
                    ok = True
                else:
                    log("STEP 11: 0 combinations completed MD", mfh)
            if not ok:
                if step_id in NON_FATAL_STEPS:
                    log(f"STEP {step_id}: FAILED (non-fatal), continuing. cmd: {bad}", mfh)
                    continue
                log(f"STEP {step_id}: FAILED (fatal). cmd: {bad}", mfh)
                log(f"--- tail {slog} ---\n{tail(slog)}", mfh)
                print(f"\nRESULT {case}: FAIL at step {step_id}")
                return 1

        ok, why = success(case)
        log(f"success check: {ok} ({why})", mfh)
        print(f"\nRESULT {case}: {'PASS' if ok else 'FAIL'} ({why})")
        return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
