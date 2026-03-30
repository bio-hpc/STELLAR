# External dependencies (STELLAR_github)

Only paths and binaries **actually referenced** by the scripts in this folder. You do **not** need the full upstream GROMACS-side tree or MetaScreener repository—only the subtrees below (plus the Singularity images you choose to use).

**Layout:** most pipeline scripts live under **`STELLAR/`** (e.g. `STELLAR/calculate_mmpbsa.py`). Entry points at the **repository root** are `run_propedia_commands_single.py`, `run_md_simulations.py`, and `sm.sh`. Paths below are “used by” those scripts (with the `STELLAR/` prefix where applicable).

## Scope: this folder vs the full pipeline

**Yes — for actually running docking and MD you usually need more** than the paths listed here. This document answers: “what do the Python/bash helpers in `STELLAR_github` call?” It does **not** enumerate everything a full STELLAR/GROMACS checkout uses to **submit** Propedia docking, generate MD inputs, or drive `execute_jobs`-style workflows.

- **Docking (e.g. Propedia):** preparing complexes, running the docking engine, and producing `Propedia_*` / PDBQT outputs typically relies on **MetaScreener and cluster scripts outside this repo** (and data trees like `Propedia_pdbqt_5Frag/`). `STELLAR_github` assumes those outputs exist and processes them onward.
- **Simulation:** `run_md_simulations.py` only tracks SLURM jobs and log patterns. **`sm.sh`** is a thin entry point into the **`GROMACS/login_node/` shell stack** (`help.sh`, `parameters.sh`, `folders_experiment.sh`, …) and the rest of the GROMACS job/submission layout — that is **much larger** than `generate_topology.py` + `g_mmpbsa`, and is required if you launch MD through that workflow.

So: keep using a **full or near-full STELLAR/GROMACS/MetaScreener deployment** for launch; use the table below to know which **minimal subtrees** the *scripts in this folder* must be able to resolve when you merge them into that tree.

---

## GROMACS (relative to `GROMACS/`)

| Path | Used by | Role |
|------|---------|------|
| `external_sw/gromacs/topology/generate_topology.py` | `generate_topologies.py` | Invoked inside `gr.simg` to build GROMACS topologies. |
| `analyze_results/Simulation_gromacs/analyze_trajectory/extra/g_mmpbsa` | `calculate_mmpbsa.py` | Default binary for MM-PBSA (`--g-mmpbsa` overrides). |

**Notes**

- `run_md_simulations.py` only mentions GROMACS-style **log patterns** (comments); it does not import files from the `GROMACS/` tree.
- `sm.sh` assumes a GROMACS-like directory layout (`login_node`, `execute_jobs.sh`, etc.) on the cluster; it is not a Python import.

---

## MetaScreener (relative to `MetaScreener/`)

| Path | Used by | Role |
|------|---------|------|
| `extra_metascreener/used_by_metascreener/convert_to.py` | `convert_combinations_pdbqt_to_mol2.py` | Converts PDBQT → MOL2 (called via Singularity + Python). |
| `extra_metascreener/used_by_metascreener/` (same folder as `convert_to.py`) | `aggregate_gn_pose_coords.py` | Adds this directory to `PYTHONPATH` and imports `save_pose_CN_coordinates` (that module must live next to `convert_to.py`, as in upstream MetaScreener). |
| `external_sw/mgltools/` (MGLTools: `pythonsh`, `prepare_ligand4.py`, …) | `calculate_score_only.py` | Receptor/ligand prep for AutoDock Vina when not using precomputed PDBQT. |
| `external_sw/gnina/gnina` | `calculate_score_only.py` | GNINA rescoring (default path; `--gnina` overrides). |
| `external_sw/overlap/overlap` | `check_overlap_combinations.py` | Steric clash check between peptide and ligand. |

**Not** required for these scripts: Propedia trees, the rest of `extra_metascreener/`, or other `external_sw` tools unless you extend the pipeline.

---

## Singularity images (paths as in code defaults / search order)

| Image(s) | Used by | Role |
|----------|---------|------|
| `gr.simg` (project root or `PATH`) | `generate_topologies.py`, `calculate_mmpbsa.py` | GROMACS / topology / g_mmpbsa runtime. |
| `singularity/metascreener.simg` | `convert_combinations_pdbqt_to_mol2.py` | RDKit + `convert_to.py` (default `--singularity-image`). |
| `singularity/metascreener_22.04.simg` | `calculate_score_only.py`, `merge_all_combinations.py` (candidate) | MGLTools + GNINA + RDKit stack. |
| `singularity/STELLAR.simg`, `singularity/BIOFRAGMEN.simg` | `merge_all_combinations.py` (candidates) | Fallback images if the default image is missing. |
| `new_ms.simg` (project root) | `calculate_rmsd_combinations.py` | Open Babel `obrms` + embedded PyMOL script. |

Override paths with each script’s CLI flags where available.

---

## Project-local (not GROMACS/MetaScreener, but required by score workflow)

| Path | Used by |
|------|---------|
| `conversion_targets/pdbqtconvert.sh` (+ layout expected by that script) | `STELLAR/calculate_score_only.py` (protein → PDBQT when needed) |

---

## Cluster / data (environment)

- **SLURM**: `sbatch` / job IDs — `run_md_simulations.py`, `sm.sh`.
- **Per-case inputs**: `complex/`, `peptide_pdb_fragments/`, `conversion_targets/`, etc., as described in `README.md`.
