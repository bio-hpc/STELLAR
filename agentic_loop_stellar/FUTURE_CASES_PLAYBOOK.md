# Context and playbook for future STELLAR cases

> Goal of this document: whoever picks this up (person or agent) should **know
> what to do without having to rediscover it**. It captures the environment, the
> recipe, the already-solved failures with their fix, and the **constraints that
> must be respected**.

---

## 0. TL;DR (quick recipe)

```bash
cd /data4/alejandro/STELLAR_JCIM/STELLAR

# One new case (reduced pilot: 2 combinations that complete MD):
agentic_loop_stellar/bin/python3 agentic_loop_stellar/run_case.py <CASE> \
    --max-valid 2 --md-pool 8 --max-organize 400 --overlap-tol 450

# Success = all_metrics_GN_<CASE>.csv exists with >=1 row where rmsd_md != NA (via obrms).
```

`<CASE>` = `<PDBID>_<CHAIN>` (e.g. `6J19_B`). Input data is read from
`improve5Frag_2/cases/<CASE>/` (or the equivalent dataset) and is **never modified**.

---

## 0.b. Running ANOTHER subset of 50 cases

The harness generalizes (all 12+ fixes are integrated). For a new subset:

**Data requirements** (same convention as `improve5Frag_2`):
1. Cases under `<DIR>/cases/<CASE>/VS_GN_<CASE>_Frag*/…` (fragments, `molecules/`,
   `receptor_original.pdb`). `<CASE>` = `<PDBID>_<PEPTIDE_CHAIN>`.
2. Reference crystals in `complex/<pdbid>_<pep>_<prot>.pdb` (mandatory for
   `rmsd_md`; the chain/copy fix depends on all of them being present).
3. Peptide fragments in `peptide_pdb_fragments/<CASE>/` (backup lookup).
4. Same containers in `singularity/` (already present).

**Launch** (point at the subset with `CASES_DIR`):
```bash
CASES_DIR=/path/to/new_subset/cases MAX_VALID=2 MD_POOL=8 \
  nohup bash agentic_loop_stellar/run_loop.sh \
  > agentic_loop_stellar/run_loop_new.out 2>&1 &
```
`run_loop.sh` with no arguments picks up every case in `CASES_DIR`; or pass a
list: `... bash run_loop.sh 1ABC_B 2XYZ_C`. Per-case result in `e2e_summary.txt`
(PASS/FAIL) and logs in `agentic_loop_stellar/logs/<CASE>.out`.

**If a new case fails**, look at the failing step (§3) — usually the existing
fixes already cover it; if a new failure appears, diagnose and integrate it just
like in the first subset (inference + feedback loop).

**Recompute `rmsd_md`** (obrms + correct chain) after the batch:
```bash
SINGULARITY_BIND=/data4 PYTHONNOUSERSITE=1 \
  agentic_loop_stellar/bin/python3 agentic_loop_stellar/recompute_realign.py <CASE ...>
```

## 1. Environment (what the cluster needs)

- CentOS 7 host: **glibc 2.17** and system `python3` = **3.6** → do NOT use the
  system `python3` for the harness. Use the shim `agentic_loop_stellar/bin/python3`
  (→ 3.9).
- SLURM (partition `all`), Singularity/Apptainer.
- Containers (in `singularity/`):
  - `STELLAR.simg` — Ubuntu 22.04 (glibc 2.35), Python deps; for compiled
    binaries (steps 2 and 4).
  - `gr.simg` — GROMACS (topology and MD).
  - `new_ms.simg` / `metascreener*.simg` — MGLTools, RDKit, OpenBabel/**obrms**, PyMOL.
- The project lives under `/data4`, which is **NOT** auto-mounted inside the containers.

## 2. Fixes ALREADY integrated (do not rediscover)

The harness (`run_case.py`) and the pipeline already apply all of this automatically:

1. **Python:** shim `bin/python3 → 3.9` prepended to `PATH`.
2. **glibc:** steps 2 and 4 (compiled binaries) run inside `STELLAR.simg`; any
   `GLIBC` error retries inside the container.
3. **Valid combinations:** **spread** sampling (`--max-organize`) + production
   overlap tolerance **450** (not top-by-score, which clashes).
4. **`/data4` bind:** `SINGULARITY_BIND=/data4` (+APPTAINER) across the whole
   environment, including nested containers.
5. **Host Python isolation:** `PYTHONNOUSERSITE=1` propagated into the containers
   (prevents `~/.local` from breaking MGLTools/RDKit/numpy).
6. **0-byte files:** deleted before steps 5/6 (otherwise they are skipped for
   "already existing" and break RDKit).
7. **`ShuttleMol/`:** symlink `ShuttleMol → GROMACS` (the topology generator uses
   the legacy name).
8. **Resume:** staging is only redone when running from ≤ step 2; the success
   check also inspects `<CASE>_results/csv_files/`.
9. **Group A (ligand topology):** `repair_ligand_topologies` repairs hydrogen
   bonds with `r=0` (cysteine S–H, etc.) after step 10.
10. **Group B (calcium ion):** removal of `resn CA`/`resn CA2` when extracting the
    peptide (`STELLAR/calculate_md_rmsd.py`).
11. **Robust `obrms` (obrms-only):** `backbone` mode pairs by (residue+name),
    handles an offset-truncated crystal and injects identical connectivity → never
    `inf`, no PyMOL fallback.
12. **Correct reference chain/copy selection:** the case suffix (`2V8W_B`) is the
    **peptide** chain; `find_crystal_structure` uses it to pick the crystal, and
    `find_crystal_candidates` + `process_simulation` try all crystal copies and
    keep the **minimum RMSD** (avoids comparing against another copy's peptide,
    tens of Å away).
    - **[FIX 2026-07-28] The guard was wrong:** multi-copy iteration only fired
      with `chain_override is None AND combo_dir is None`. But `combo_dir` is just
      the combination folder (for topology), **not** the crystal source (which
      always comes from `complex/`). Since the batch runner **always** sets
      `combo_dir`, iteration never fired and it kept **one** copy (often the wrong
      receptor's). Fixed to `if chain_override is None:` in `process_simulation`.
      Iteration only acts when `complex/` has >1 copy, so single-copy cases are
      unchanged.
    - **Measured effect (50 improve5Frag cases):** removed the 2 outliers:
      `3KUT_D` (combos 23.61/23.94 Å → 5.46/5.38 Å, copy `3kut_D_A`) and
      `1QD6_B` (combos 21.52/23.16 Å → 5.33/5.88 Å, copies `1qd6_B_C`/`B_D`;
      different combos match different copies, which is why the per-combo `min` is
      key). After the fix, **global max = 6.44 Å** (1ZAF_D) and **0 cases > 7 Å**.
    - **Recompute after the fix** (without redoing MD): re-run step 13 → 16 → 17
      with the shim + bind (see §5, "rmsd_md far too high" entry).
    - **[FIX 2026-07-30] Sibling copies in OTHER chains (lattice symmetry):**
      sometimes the crystal stores identical peptide copies under **different chain
      ids** (e.g. 6g0y as chains H/I/J). The docked pose can land on any copy, so
      comparing only against the case's own chain leaves combos at ~40 Å in `obrms`
      despite `RMSD_ALIGN_PYMOL ≈ 0.4 Å` (correct conformation, different lattice
      position). `find_crystal_candidates` now **extends** the candidates with the
      sibling copies whose **peptide sequence is identical** (guarded by sequence:
      backbone `obrms` is sequence-agnostic and would be fooled by a genuinely
      different peptide that happens to sit nearby; e.g. `6g0y_G_E` has 16 residues
      vs 13 → **excluded**). Effect on `improve5Frag_5`: `6G0Y_I`/`6G0Y_J` combos
      4089/4162 (39–40 Å → 5.29/5.63 Å, copy `6g0y_H_C`) and one combo that was NA
      became 4.19 Å. Subset global max after the fix: **6.52 Å** (`1HQQ_E`),
      0 cases > 7 Å.
13. **[FIX 2026-07-29] `.pdbqt`-only receptor → correct `VS_GR_*` naming**
    (`generate_workflow_commands.py`, `find_receptor_pdb`): some cases ship the
    receptor only as `.pdbqt` in `best_scores` (no `.pdb`), and the pipeline falls
    back to `receptor_original.pdb`. Do **NOT** pass that basename as-is: it
    propagates as the receptor id into `VS_GR_receptor_original_*`, which (a) makes
    `run_md_simulations --check-existing` match **other** cases' folders (same
    `combination_N`) and **skip MD**, and (b) makes `calculate_md_rmsd` (filters
    `VS_GR_*<case>*`) find nothing → **`rmsd_md` NA**. The fix copies
    `receptor_original.pdb` under the real `.pdbqt` basename (e.g. `4p3w_A.pdb`),
    yielding a per-case `VS_GR_4p3w_A_*`.
    - **Typical symptom:** step 11 says `N combinations already have PDB (will skip)`
      / `0 combinations will run` on a clean case, and step 13 aborts with
      `No VS_GR_* folders found for case <CASE>`. `VS_GR_receptor_original_*`
      folders appear.
    - **Fixed case:** `4P3W_G` (receptor `4p3w_A`, `.pdbqt` only) → combos
      4.29/4.35 Å. `improve5Frag_4` ends **50/50** with valid `rmsd_md`.
14. **[FIX 2026-07-30] `VS_GR_*` (MD) archived under `<CASE>_results/md_simulations/`**
    (`organize_workflow_results.py`, step 17): they used to remain in the project
    **root** (the `VS_GR_*<case>*` patterns didn't match because the folder is
    named by the RECEPTOR, `VS_GR_<pdbid>_<recchain>_...`, not by the case chain)
    and saturated the root / the disk. They are now identified by the `vs_folder`
    column of `resultados_rmsd_md_<CASE>.csv` (fallback: pdbid glob) and moved into
    the `md_simulations/` subfolder. It is a `move` within `/data4` (rename, no disk
    cost). `calculate_md_rmsd.find_vs_folders()` **also searches
    `*_results/md_simulations/`**, so step 13 recompute (§5) still works after
    organizing. Note: cases sharing the same pdbid+peptide (e.g. `6G0Y_I`/`6G0Y_J`)
    **share** the same `VS_GR_*`; the first one to organize moves them and the rest
    skip (find_vs_folders still locates them).

## 3. Quick per-step diagnosis (if a new case fails)

| Step | Expected output | If it fails, check… |
|---|---|---|
| 2, 4 | filtered combinations | GLIBC error? → must run in `STELLAR.simg`. 0 valid? → raise `--max-organize` or check `--overlap-tol`. |
| 5 | `.mol2` from `.pdbqt` | `convert_to.py` not found? → bind `/data4`. 0-byte `.mol2`? → cleanup hook. |
| 6 | RDKit merge | empty/corrupt `.mol2`? → delete and regenerate. |
| 7 | receptor located | **Group C**: if the receptor is not found, check names/chain (see §5). |
| 10 | GROMACS topology | HETATM in receptor → already sanitized by `sanitize_target_pdb`. `ShuttleMol/` → symlink. |
| 11 | MD `VS_GR_*/molecules/*_complex_md.pdb` | **Group A**: `Assertion ip.constr.dA > 0` → `r=0` bond in `.itp` (repaired). |
| 13 | `resultados_rmsd_md_<CASE>.csv` with `rmsd` | **Group B**: all NA due to `obrms=inf` → ion/solvent leaked in or different graphs (see §4). |
| 16 | `all_metrics_GN_<CASE>.csv` with `rmsd_md` | consolidates; if 13 produced non-NA, so does this. |

## 4. How to compute `rmsd_md` correctly (mandatory)

- **Method:** `obrms` (OpenBabel) over the **backbone** (N, CA, C, O), **in-place**
  (no re-superposition). Implemented in `robust_backbone_rmsd`.
- **Pairing:** by **(residue, atom name)**.
  - Crystal: ordered by `resSeq` then N,CA,C,O.
  - Sim (single residue `L01`): PyMOL returns it **grouped by name** (all N, then
    all CA, …) → **interleave by index** to recover residue order.
- **Truncated crystal** (unresolved terminal residues missing): the crystal is a
  shorter contiguous block → **slide** it over the sim residues and pick the best
  overlap (recovers the real sequence offset).
- **Avoid `inf`:** perceive connectivity **once from the crystal** (SDF) and
  **inject** the sim coordinates into that same atom/bond block → both share a
  graph → `obrms` always returns a value.
- **Recompute without redoing MD/PyMOL:** use `agentic_loop_stellar/recompute_obrms.py`
  on the already-extracted peptides.

## 5. Known failures and their fix (reference)

- **`GLIBC_2.34 not found`** → run the step in `STELLAR.simg`.
- **0 valid combinations** → spread sampling + `--overlap-tol 450`.
- **`can't open MetaScreener/.../convert_to.py`** → bind `/data4`.
- **broken `numpy`/MGLTools** → `PYTHONNOUSERSITE=1`.
- **0-byte `.mol2`** → delete before 5/6.
- **"Error en el fichero de configuracion" (step 10)** → symlink `ShuttleMol → GROMACS`.
- **`Assertion failed: ip.constr.dA > 0` (MD)** → `r=0` bond in `.itp` (Group A);
  repair H-bond length.
- **`obrms` = `inf` / all `rmsd_md` NA** → leaked ion (Group B: add to the removal
  list) and/or use the robust pairing with shared connectivity.
- **Receptor not found (step 7, Group C)** → covered by fix #13 (`.pdbqt`-only
  receptor fallback with correct naming).
- **`rmsd_md` far too high (20–45 Å) with the receptor well superposed** → wrong
  reference chain/copy: comparing against the peptide of another chain or another
  crystal copy. Tell-tale sign: `RMSD_ALIGN_PYMOL ≈ 0.3–0.6 Å` (that is the
  **receptor** fit, not the peptide) while `obrms` gives ~20+ Å. **Already fixed**
  (§2 point 12): all `complex/` copies are tried and the minimum per-combo `obrms`
  is taken. To **recompute a case** already simulated without redoing MD:
  ```bash
  cd /data4/alejandro/STELLAR_JCIM/STELLAR
  export SINGULARITY_BIND=/data4 APPTAINER_BIND=/data4 PYTHONNOUSERSITE=1
  export PATH="$PWD/agentic_loop_stellar/bin:$PATH"   # python3.9 shim
  # 1) recompute RMSD (tries all copies, picks the minimum):
  python3 STELLAR/calculate_md_rmsd.py --case <CASE> \
      --output-csv resultados_rmsd_md_<CASE>.csv --output-dir md_rmsd_peptides_<CASE>
  # 2) dump into results and re-merge (inputs are already in <CASE>_results/):
  cp -f resultados_rmsd_md_<CASE>.csv <CASE>_results/csv_files/
  cd <CASE>_results/csv_files && python3 ../../STELLAR/merge_all_metrics.py --prefix-type GN \
      --combinations-file valid_fragment_combinations_GN_no_overlap_<CASE>.csv \
      --rmsd-fragments-file ../valid_combinations_GN_<CASE>/valid_no_overlap/rmsd_results.csv \
      --md-rmsd-file resultados_rmsd_md_<CASE>.csv --score-only-file score_only_results_<CASE>.csv \
      --mmpbsa-file mmpbsa_results_<CASE>.csv --fragment-energies-file fragment_energies_<CASE>.csv \
      --output all_metrics_GN_<CASE>.csv
  ```
  **Caution:** without `SINGULARITY_BIND=/data4` the container cannot see the
  project and the PyMOL extraction comes out **empty** (obrms fails silently).
  And use the 3.9 shim, not the host `python3` (3.6 lacks
  `subprocess.capture_output`).

## 6. CONSTRAINTS (must be respected)

1. **Do not modify the original data** in `improve5Frag_2/cases/` (or other
   datasets). All work goes in `agentic_loop_stellar/` / the `<CASE>_GN/` staging
   directories with symlinks and write isolation.
2. **`rmsd_md` ONLY via `obrms`.** PyMOL `align` fallback is forbidden (it
   re-superposes → incomparable values). If `obrms` cannot reconcile the peptides
   → **NA**, not a PyMOL value.
3. **Do not use the host `python3` (3.6).** Use the `bin/python3` shim (3.9).
4. **Compiled-binary steps (2, 4) run in `STELLAR.simg`** (glibc).
5. **Always propagate `SINGULARITY_BIND=/data4` and `PYTHONNOUSERSITE=1`** to the
   containers.
6. **Topology repairs (Group A) edit ONLY the working copy**
   (`valid_GN_<CASE>_final/`), never the original `.itp` files.
7. **Pilot = 1–2 combinations** (`--max-valid 2`). For real production, raise to
   ~100 combinations and expect much longer MD times.
8. **Success = `all_metrics_GN_<CASE>.csv` with ≥1 `rmsd_md` ≠ NA** (via `obrms`).
   A high `rmsd_md` (e.g. >20 Å) is a poor pose but **still meets the criterion**;
   structural quality is a separate analysis.

## 7. Code pointers

- RMSD pipeline: `STELLAR/calculate_md_rmsd.py`
  (`robust_backbone_rmsd`, `_bb_res_groups_*`, `_bb_best_offset`,
  `_obrms_shared_connectivity`, `calculate_rmsd_obrms`).
- Per-case driver: `agentic_loop_stellar/run_case.py`
  (`prep_environment`, `repair_ligand_topologies`, `remove_empty_files`,
  `apply_pilot_overrides`, `count_completed_md`, `success`).
- Staging: `agentic_loop_stellar/stage_case.sh`.
- obrms recompute: `agentic_loop_stellar/recompute_obrms.py`.
- Detailed historical log: `agentic_loop_stellar/LOG_AGENTIC_LOOP_STELLAR.md`.
- Full report: `agentic_loop_stellar/INFORME_STELLAR.md`.
- Parallel AI loop (Cursor SDK): `agentic_loop_stellar/auto_fix_loop.py`.

## 8. Parallel AI agentic loop (`auto_fix_loop.py`)

Orchestrator that runs cases **in parallel** and, when one fails, launches a
Cursor agent (Python SDK, **local** runtime) that receives THIS playbook as rules
+ the tail of the failing step's log, **diagnoses and applies a targeted fix**,
then re-runs the case (resuming from the failed step). At the end it saves an
**MD + PDF** report with what each agent considered and the diff it applied.

Requirements and setup on this host (CentOS 7 / el7) — ALREADY INSTALLED:
- The SDK requires **Python ≥3.10**, but the host only has 3.9. A **standalone
  Python 3.11 was installed with `uv`** in `agentic_loop_stellar/.venv-sdk/`
  (with `cursor-sdk` inside). The 3.9 shim (`bin/python3`) is still used for
  `run_case.py`; only the orchestrator runs under 3.11.
- `cursor-sdk` ships its own `node` that needs GLIBC_2.27+/CXXABI_1.3.11 (el7 has
  neither). A launcher `.sdkenv/cursor-sdk-bridge` was created that redirects the
  bridge to the **cursor-server node** (v20+, compatible), via
  `CURSOR_SDK_BRIDGE_BIN`. It auto-detects the node even if the hash changes after
  updates (override with `CURSOR_SDK_NODE`).
- The SDK's default local store uses `node:sqlite` (requires Node ≥22.13, which
  does NOT run on el7 due to glibc). `auto_fix_loop.py` forces a **JSONL store**
  (`LocalAgentStoreConfig(type="jsonl", root_dir=.sdkenv/agent_store)`) that runs
  inside the bridge without sqlite. Without this, the agent fails to start.
- `export CURSOR_API_KEY="cursor_..."` (user or service-account key).
- Model Opus 4.8 = id `claude-opus-4-8` (`models.list` exposes it; `pick_model`
  prefers it automatically).

**Parallelism at step 10 (ACPYPE) — RESOLVED:** `generate_topologies.py` runs
ACPYPE, which creates temp dirs with a FIXED name
(`fragmento_final_charge_drift` / `.acpype_tmp_...`) in the project root (shared
CWD). With `--workers >1`, two cases reaching step 10 at the same time
**collided** and one failed with
`ACPYPE FAILED: ... fragmento_final_charge_drift_AC.prmtop` (seen with
1FIW_L+1NRL_C and with 2B9I_C+1QD6_A). It cannot be isolated by changing the CWD
because ShuttleMol's inner pipeline assumes CWD == project root (it resolves
`config.cfg`, `external_sw`, `check_target.py`… relative to CWD). **Applied fix:**
an inter-process `flock` (`.topo_acpype.lock`) in `generate_topologies.py`
serializes only the ACPYPE invocation, keeping the CWD at the root. Validated:
2B9I_C + 1QD6_A with `--workers 2` → both PASS (2 valid `rmsd_md` each). It is now
safe to use `--workers >1`.

To reinstall the SDK environment from scratch (if `.venv-sdk/` is deleted):
```bash
cd agentic_loop_stellar
UV=.sdkenv/uv/uv
UV_PYTHON_INSTALL_DIR=$PWD/.sdkenv/pythons UV_CACHE_DIR=$PWD/.sdkenv/cache \
  $UV venv .venv-sdk --python 3.11
$UV pip install --python .venv-sdk/bin/python cursor-sdk
```

Usage (ALWAYS via the wrapper, which sets venv 3.11 + the bridge node):
```bash
export CURSOR_API_KEY="cursor_..."
agentic_loop_stellar/run_auto_fix.sh --workers 3 --max-fix-attempts 1 3HXI_C 6IFJ_C
# all cases in a subset:
CASES_DIR=/path/other_subset/cases agentic_loop_stellar/run_auto_fix.sh --workers 4
# smoke test without AI (parallel pipeline only, no key needed):
agentic_loop_stellar/run_auto_fix.sh --no-ai CASE
```

Model (**Opus 4.8**):
- Without `--model`, the script queries `Cursor.models.list()` and picks an Opus
  (prefers 4.8); if the key sees none, it uses `auto`.
- Force it: `--model <id>` or `export STELLAR_FIX_MODEL=<id>`.

Design notes / limits:
- Cases run in parallel, but the **AI code-editing phase is serialized** (an
  `EDIT_LOCK`) so two agents don't touch STELLAR's shared code at once. Full
  isolation = one git worktree per case.
- The real bottleneck is SLURM/MD; `--workers` should not exceed what the queue
  tolerates.
- The playbook constraints are respected: do not touch original data, `rmsd_md`
  only via `obrms` (never PyMOL fallback).
- Outputs in `agentic_loop_stellar/auto_fix_reports/RUN_<timestamp>/`:
  `informe_auto_fix.md`, `.tex`, `.pdf` and `<CASE>_attemptN.patch` per fix.
