# STELLAR — minimal GN scripts

This directory contains **only the scripts and utilities** needed to run the STELLAR pipeline in **GN** mode, from Propedia docking through **`all_metrics_*.csv`**.

Large data and Singularity images are **not** included. You only need the **GROMACS / MetaScreener subtrees and binaries** listed in **DEPENDENCIES.md** (not whole upstream repos).

## Layout

- **Repository root (this folder):** scripts that **launch** external commands or the MD flow — `run_propedia_commands_single.py`, `run_md_simulations.py`, `sm.sh`.
- **`STELLAR/`:** all other pipeline steps (metrics, conversion, merge, topology, etc.). Invoke them from the project root, e.g. `python3 STELLAR/merge_all_metrics.py ...`.

**Deployment:** merge the contents of `STELLAR_github/` into the root of a full STELLAR checkout; do not replace the entire repository.

Scripts assume **paths relative to the project root** where `GROMACS/`, `MetaScreener/`, `singularity/`, `conversion_targets/`, `complex/`, `peptide_pdb_fragments/`, etc. coexist.

## GN pipeline (logical order)

| # | Script |
|---|--------|
| 0 | `run_propedia_commands_single.py` |
| 1 | `STELLAR/save_pose_CN_coordinates.py`, `STELLAR/aggregate_gn_pose_coords.py` |
| 2 | `STELLAR/filter_fragment_combinations` |
| 3 | `STELLAR/organize_valid_combinations.py` |
| 4 | `STELLAR/check_overlap_combinations.py` |
| 5 | `STELLAR/convert_combinations_pdbqt_to_mol2.py` |
| 6 | `STELLAR/merge_all_combinations.py` (+ `STELLAR/relax_merge_mol2.py`, `STELLAR/fix_charge_drift.py`) |
| 7 | `STELLAR/prepare_final_combinations.py` |
| 8 | `STELLAR/fix_zero_charge_atoms.py` |
| 9 | `STELLAR/generate_topologies.py` |
| 10 | `run_md_simulations.py` + `sm.sh` |
| 11 | `STELLAR/calculate_rmsd_combinations.py` |
| 12 | `STELLAR/calculate_md_rmsd.py` |
| 13 | `STELLAR/filter_valid_combinations_csv.py` |
| 14 | `STELLAR/calculate_score_only.py` |
| 15 | `STELLAR/calculate_mmpbsa.py` |
| 16 | `STELLAR/calculate_fragment_energies.py` |
| 17 | **`STELLAR/merge_all_metrics.py`** → `all_metrics_GN_<case>.csv` |

Command reference: `docs/WORKFLOW_GN_timed_reference.txt`.
