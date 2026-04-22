#!/usr/bin/env python3
"""
Generate a SLURM shell script that runs the full GN workflow sequentially.

Reuses detection helpers from generate_workflow_commands.py.
Also writes run_md_step_<CASE>.sh (step 11 only) for separate MD submission.

Layout: pipeline scripts under STELLAR/; launchers (run_propedia_commands_single.py,
run_md_simulations.py, sm.sh) at project root.

Environment (optional overrides for #SBATCH):
  SLURM_ACCOUNT, SLURM_QOS, SLURM_TIME, SLURM_PARTITION
"""

import os
import sys
from pathlib import Path
from typing import Optional

_HERE = Path(__file__).resolve().parent
_PROJECT = _HERE
if str(_PROJECT) not in sys.path:
    sys.path.insert(0, str(_PROJECT))

from generate_workflow_commands import (
    detect_case_and_base_dir,
    detect_fragment_count,
    detect_propedia_folder,
    extract_fragment_count_from_propedia_name,
    find_receptor_pdb,
    get_case_crystal_pdb,
    get_existing_fragment_numbers,
    project_root,
)

ST = "STELLAR"


def generate_workflow_job(
    case_name: str,
    output_file: Optional[str] = None,
    base_dir: str = ".",
    start_from_step: int = 0,
) -> str:
    """Write main SLURM job script and run_md_step_<CASE>.sh under project root."""
    if output_file is None:
        output_file = f"job_workflow_{case_name}.sh"

    case_upper = case_name.upper()
    case_lower = case_name.lower()
    proj = str(project_root())

    propedia_folder = detect_propedia_folder(case_name, base_dir)
    print(f"  Propedia folder: {propedia_folder}")

    crystal_pdb = get_case_crystal_pdb(case_upper, base_dir)
    if crystal_pdb:
        print(f"  Crystal (complex/): {crystal_pdb}")
    else:
        print(f"  Warning: no crystal in complex/ for {case_upper}")
    crystal_comment = (
        f"# Crystal for this case: {crystal_pdb}"
        if crystal_pdb
        else "# Crystal: not found under complex/"
    )

    num_fragments, frag_base_dir = detect_fragment_count(case_name, base_dir)
    if num_fragments is None:
        print(f"  Warning: no VS_GN_{case_upper}_Frag* (step 0 may not have run)")
        n_prop = extract_fragment_count_from_propedia_name(propedia_folder)
        if n_prop is not None:
            num_fragments = n_prop
            print(f"  Using {num_fragments} fragments from Propedia name: {propedia_folder}")
        else:
            num_fragments = 5
            print("  Using default 5 fragments")
        frag_base_dir = f"{case_upper}_GN"
    else:
        print(f"  Fragments: {num_fragments}")
        if frag_base_dir and frag_base_dir != base_dir:
            print(f"  Location: {frag_base_dir}")

    frag_path_prefix = "" if frag_base_dir == "." else f"{frag_base_dir}/"

    pdb_file = find_receptor_pdb(case_name, base_dir, frag_base_dir)
    if pdb_file:
        pdb_arg = pdb_file
        pdb_search_suffix = ""
    else:
        pdb_code = case_lower.split("_")[0] if "_" in case_lower else case_lower
        pdb_arg = f"{pdb_code}_*.pdb"
        best_scores = os.path.join(frag_base_dir or ".", f"VS_GN_{case_upper}_Frag1", "results", "best_scores")
        pdb_search_suffix = f" --pdb-search-dir {best_scores}"

    frag_nums = get_existing_fragment_numbers(case_name, frag_base_dir)
    if not frag_nums and num_fragments and num_fragments > 0:
        frag_nums = list(range(1, num_fragments + 1))

    frag_lines = []
    for i in frag_nums:
        frag_dir = f"{frag_path_prefix}VS_GN_{case_upper}_Frag{i}"
        frag_lines.append(f"python3 {ST}/save_pose_CN_coordinates.py {frag_dir}/molecules/")
    frag_block = "\n  ".join(frag_lines)

    if frag_nums:
        num_fragments = max(frag_nums)

    aggregate_base = frag_base_dir if frag_base_dir != "." else base_dir
    if aggregate_base == ".":
        aggregate_cmd = f"python3 {ST}/aggregate_gn_pose_coords.py --case-filter {case_upper}"
    else:
        aggregate_cmd = f"python3 {ST}/aggregate_gn_pose_coords.py {aggregate_base} --case-filter {case_upper}"

    filter_base = aggregate_base
    if filter_base.startswith("./"):
        filter_base = filter_base[2:]
    if filter_base == ".":
        filter_csv = "valid_fragment_combinations.csv"
        organize_base = "."
    else:
        filter_csv = f"{filter_base}/valid_fragment_combinations.csv"
        organize_base = filter_base

    base_dir_rmsd_arg = f" --base-dir {case_upper}" if base_dir != "." else ""
    case_base = case_upper.split("_")[0] if "_" in case_upper else case_upper

    account = os.environ.get("SLURM_ACCOUNT", "ucam37")
    qos = os.environ.get("SLURM_QOS", "gp_resa")
    walltime = os.environ.get("SLURM_TIME", "3-00:00:00")

    start_val = int(start_from_step) if start_from_step is not None else 0
    g_mmpbsa = "./GROMACS/analyze_results/Simulation_gromacs/analyze_trajectory/extra/g_mmpbsa"

    job_content = f"""#!/bin/sh
#SBATCH --output={proj}/output_{case_upper}.out
#SBATCH --error={proj}/output_{case_upper}.err
#SBATCH -J STELLAR_{case_upper}
#SBATCH --time={walltime}
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --qos={qos}
#SBATCH -A {account}
#
# Each step skips work if outputs already exist.
# Resume from step N: START_FROM_STEP=10 sbatch this_script.sh  (skips steps 0–9)
# {crystal_comment}

set -e

echo "Starting workflow for {case_upper}"
date

START_FROM_STEP=${{START_FROM_STEP:-{start_val}}}

cd {proj}

# STEP 0: Propedia docking
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 0 ]; then
  echo "STEP 0: skipped (START_FROM_STEP=$START_FROM_STEP)."
else
echo "STEP 0: Propedia docking..."
python3 run_propedia_commands_single.py --complex {case_upper} --propedia-folder {propedia_folder} --skip-if-exists
fi

# STEP 1: C/N coordinates
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 1 ]; then
  echo "STEP 1: skipped."
else
if [ -f "{filter_csv}" ]; then
  echo "STEP 1: {filter_csv} exists; skipping pose CN + aggregate."
else
  echo "STEP 1: C/N coordinates..."
  {frag_block}
  echo "Aggregating..."
  {aggregate_cmd}
fi
fi

# STEP 2: Distance filter
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 2 ]; then
  echo "STEP 2: skipped."
else
if [ -f "{filter_csv}" ]; then
  echo "STEP 2: {filter_csv} exists; skipping distance filter."
else
  echo "STEP 2: distance filter..."
  ./{ST}/filter_fragment_combinations {filter_base} VS_GN {case_upper} --threshold-gn 1.8
fi
fi

# STEP 3: Organize
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 3 ]; then
  echo "STEP 3: skipped."
else
if [ -d "valid_combinations_GN_{case_upper}" ] && [ -n "$(ls -A valid_combinations_GN_{case_upper} 2>/dev/null)" ]; then
  echo "STEP 3: valid_combinations_GN_{case_upper}/ exists; skipping."
else
  echo "STEP 3: organize combinations..."
  python3 {ST}/organize_valid_combinations.py {filter_csv} {organize_base} valid_combinations_GN_{case_upper} GN {case_upper}
fi
fi

# STEP 4: Overlap
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 4 ]; then
  echo "STEP 4: skipped."
else
if [ -d "valid_combinations_GN_{case_upper}/valid_no_overlap" ] && [ -n "$(ls -A valid_combinations_GN_{case_upper}/valid_no_overlap 2>/dev/null)" ]; then
  echo "STEP 4: valid_no_overlap exists; skipping."
else
  echo "STEP 4: overlap filter..."
  python3 {ST}/check_overlap_combinations.py valid_combinations_GN_{case_upper} GN 450 0 100
fi
fi

# STEP 5: PDBQT -> MOL2
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 5 ]; then
  echo "STEP 5: skipped."
else
if [ -n "$(ls valid_combinations_GN_{case_upper}/valid_no_overlap/*.mol2 2>/dev/null)" ] || [ -n "$(ls valid_combinations_GN_{case_upper}/valid_no_overlap/combination_*/*.mol2 2>/dev/null)" ]; then
  echo "STEP 5: MOL2 present; skipping conversion."
else
  echo "STEP 5: PDBQT to MOL2..."
  python3 {ST}/convert_combinations_pdbqt_to_mol2.py GN --base-dir valid_combinations_GN_{case_upper}/valid_no_overlap
fi
fi

# STEP 6: Merge fragments
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 6 ]; then
  echo "STEP 6: skipped."
else
if [ -n "$(ls -d valid_combinations_GN_{case_upper}/valid_no_overlap/combination_*/query_combination_* 2>/dev/null | head -1)" ]; then
  echo "STEP 6: query_combination_* exists; skipping merge."
else
  echo "STEP 6: merge fragments..."
  python3 {ST}/merge_all_combinations.py GN --base-dir-gn valid_combinations_GN_{case_upper}/valid_no_overlap
fi
fi

# STEP 7: Final structures
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 7 ]; then
  echo "STEP 7: skipped."
else
if [ -d "valid_GN_{case_upper}_final" ] && [ -n "$(ls -A valid_GN_{case_upper}_final/combination_* 2>/dev/null)" ]; then
  echo "STEP 7: valid_GN_{case_upper}_final exists; skipping."
else
  echo "STEP 7: prepare final..."
  python3 {ST}/prepare_final_combinations.py GN --pdb-file {pdb_arg}{pdb_search_suffix} --base-dir-gn valid_combinations_GN_{case_upper}/valid_no_overlap --output-dir-gn valid_GN_{case_upper}_final
fi
fi

# STEP 8: Fix charges
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 8 ]; then
  echo "STEP 8: skipped."
else
if [ -n "$(ls valid_GN_{case_upper}_final/combination_*/*.top 2>/dev/null)" ]; then
  echo "STEP 8: topologies exist; skipping charge fix."
else
  echo "STEP 8: fix charges..."
  python3 {ST}/fix_zero_charge_atoms.py GN --base-dir-gn valid_GN_{case_upper}_final --fix-aromatic
fi
fi

# STEP 10: Topologies (step 9 calcium script omitted)
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 10 ]; then
  echo "STEP 10: skipped."
else
if python3 {ST}/generate_topologies.py GN --base-dir-gn valid_GN_{case_upper}_final --check-only 2>/dev/null; then
  echo "STEP 10: topologies complete; skipping."
else
  echo "STEP 10: generate topologies..."
  python3 {ST}/generate_topologies.py GN --base-dir-gn valid_GN_{case_upper}_final
fi
fi

# STEP 11: MD (separate sbatch script to limit memory on main job)
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 11 ]; then
  echo "STEP 11: skipped."
else
if [ -n "$(ls valid_GN_{case_upper}_final/combination_*/*.tpr 2>/dev/null)" ] || [ -n "$(ls valid_GN_{case_upper}_final/combination_*/*.edr 2>/dev/null)" ]; then
  echo "STEP 11: MD outputs exist; skipping launch."
elif [ -n "$(find . -maxdepth 4 -type f -path '*/molecules/*.pdb' 2>/dev/null | grep -i "VS_GR.*{case_base}" | head -1)" ]; then
  echo "STEP 11: VS_GR PDBs found; skipping MD launch."
else
  echo "STEP 11: submitting MD sub-job..."
  MD_JOBID=$(sbatch --parsable run_md_step_{case_upper}.sh)
  echo "  MD job $MD_JOBID submitted; waiting..."
  while squeue -j $MD_JOBID -h 2>/dev/null | grep -q .; do sleep 30; done
  echo "  MD job finished."
fi
fi

# STEP 12: Fragment RMSD
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 12 ]; then
  echo "STEP 12: skipped."
else
if [ -f "valid_combinations_GN_{case_upper}/valid_no_overlap/rmsd_results.csv" ]; then
  echo "STEP 12: rmsd_results.csv exists; skipping."
else
  echo "STEP 12: fragment RMSD..."
  python3 {ST}/calculate_rmsd_combinations.py valid_combinations_GN_{case_upper}/valid_no_overlap peptide_pdb_fragments/{case_upper} GN{base_dir_rmsd_arg}
fi
fi

# STEP 13: Peptide RMSD after MD
if [ -d "md_rmsd_peptides_{case_upper}" ] && [ -n "$(ls -A md_rmsd_peptides_{case_upper} 2>/dev/null)" ]; then
  echo "STEP 13: md_rmsd_peptides exists; skipping."
else
  echo "STEP 13: MD peptide RMSD..."
  python3 {ST}/calculate_md_rmsd.py --output-csv resultados_rmsd_md_{case_upper}.csv --case {case_upper} --output-dir md_rmsd_peptides_{case_upper}
fi

# STEP 14
if [ -f "valid_fragment_combinations_GN_no_overlap_{case_upper}.csv" ]; then
  echo "STEP 14: filtered CSV exists; skipping."
else
  echo "STEP 14: filter valid combinations CSV..."
  python3 {ST}/filter_valid_combinations_csv.py GN --combinations-dir valid_combinations_GN_{case_upper}/valid_no_overlap --input-csv {filter_csv} --output-csv valid_fragment_combinations_GN_no_overlap_{case_upper}.csv
fi

# STEP 15
if [ -f "score_only_results_{case_upper}.csv" ]; then
  echo "STEP 15: skipping score_only."
else
  echo "STEP 15: score_only..."
  python3 {ST}/calculate_score_only.py --peptides-dir md_rmsd_peptides_{case_upper}/extracted_peptides --output-csv score_only_results_{case_upper}.csv --singularity-image singularity/metascreener_22.04.simg --singularity-bind /gpfs/
fi

# STEP 15b
if [ -f "mmpbsa_results_{case_upper}.csv" ]; then
  echo "STEP 15b: skipping MM/PBSA."
else
  echo "STEP 15b: MM/PBSA..."
  python3 {ST}/calculate_mmpbsa.py --case {case_upper} --type GN --gr-simg singularity/gr.simg --g-mmpbsa-bin {g_mmpbsa} --output-csv mmpbsa_results_{case_upper}.csv
fi

# STEP 15c
if [ -f "fragment_energies_{case_upper}.csv" ]; then
  echo "STEP 15c: skipping fragment energies."
else
  echo "STEP 15c: fragment energies..."
  python3 {ST}/calculate_fragment_energies.py --case {case_upper} --type GN --combinations-file valid_fragment_combinations_GN_no_overlap_{case_upper}.csv --output-csv fragment_energies_{case_upper}.csv
fi

# STEP 16
if [ -f "all_metrics_GN_{case_upper}.csv" ]; then
  echo "STEP 16: all_metrics exists; skipping."
else
  echo "STEP 16: merge_all_metrics..."
  python3 {ST}/merge_all_metrics.py --prefix-type GN \\
    --combinations-file valid_fragment_combinations_GN_no_overlap_{case_upper}.csv \\
    --rmsd-fragments-file valid_combinations_GN_{case_upper}/valid_no_overlap/rmsd_results.csv \\
    --md-rmsd-file resultados_rmsd_md_{case_upper}.csv \\
    --score-only-file score_only_results_{case_upper}.csv \\
    --mmpbsa-file mmpbsa_results_{case_upper}.csv \\
    --fragment-energies-file fragment_energies_{case_upper}.csv \\
    --output all_metrics_GN_{case_upper}.csv
fi

# STEP 17
if [ -d "{case_upper}_results" ] && [ -n "$(ls -A {case_upper}_results 2>/dev/null)" ]; then
  echo "STEP 17: results folder exists; skipping."
else
  echo "STEP 17: organize results..."
  python3 organize_workflow_results.py {case_upper} --target-dir {case_upper}_results
fi

echo "Workflow finished for {case_upper}"
date
"""

    out_path = Path(proj) / output_file
    out_path.write_text(job_content, encoding="utf-8")
    out_path.chmod(0o755)

    md_script = f"run_md_step_{case_upper}.sh"
    md_content = f"""#!/bin/sh
# MD only (step 11); submitted by main workflow to isolate memory use.
#SBATCH --output={proj}/output_md_{case_upper}.out
#SBATCH --error={proj}/output_md_{case_upper}.err
#SBATCH -J MD_{case_upper}
#SBATCH --time={walltime}
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --qos={qos}
#SBATCH -A {account}

set -e
cd {proj}
python3 run_md_simulations.py GN --base-dir-gn valid_GN_{case_upper}_final --sequential --check-existing
"""
    md_path = Path(proj) / md_script
    md_path.write_text(md_content, encoding="utf-8")
    md_path.chmod(0o755)

    print(f"Written: {out_path}")
    print(f"  MD helper: {md_path}")
    print(f"  Case: {case_upper}  Fragments: {num_fragments}")
    print(f"  Submit: sbatch {output_file}")
    return str(out_path)


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: python3 generate_workflow_job.py <CASE_OR_FOLDER> [output_job.sh]")
        print("Example: python3 generate_workflow_job.py 1CJR_A")
        sys.exit(1)

    input_path = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else None
    case_name, base_dir = detect_case_and_base_dir(input_path)
    print(f"Case: {case_name}")
    if base_dir != ".":
        print(f"Base directory: {base_dir}")
    os.chdir(project_root())
    generate_workflow_job(case_name, out, base_dir)


if __name__ == "__main__":
    main()
