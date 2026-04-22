#!/bin/sh
#SBATCH --output=/gpfs/scratch/ucam37/alejandro/stellar_github_Def/output_1AQD_C.out
#SBATCH --error=/gpfs/scratch/ucam37/alejandro/stellar_github_Def/output_1AQD_C.err
#SBATCH -J STELLAR_1AQD_C
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --qos=gp_resa
#SBATCH -A ucam37
#
# Each step skips work if outputs already exist.
# Resume from step N: START_FROM_STEP=10 sbatch this_script.sh  (skips steps 0–9)
# # Crystal for this case: /gpfs/scratch/ucam37/alejandro/stellar_github_Def/complex/1aqd_C_B.pdb

set -e

echo "Starting workflow for 1AQD_C"
date

START_FROM_STEP=${START_FROM_STEP:-0}
CASE_ID="1AQD_C"
DOCKING_PREFIX="VS_GN"
DOCKING_CASE_DIR="${CASE_ID}_GN"

cd /gpfs/scratch/ucam37/alejandro/stellar_github_Def

# Docking outputs may live in root (legacy) or inside ${CASE_ID}_GN (new organizer behavior).
if [ -d "${DOCKING_CASE_DIR}/${DOCKING_PREFIX}_${CASE_ID}_Frag1" ]; then
  DOCKING_BASE="${DOCKING_CASE_DIR}"
else
  DOCKING_BASE="."
fi
VALID_COMBINATIONS_CSV="${DOCKING_BASE}/valid_fragment_combinations.csv"
echo "Resolved docking base: ${DOCKING_BASE}"

# STEP 0: Propedia docking
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 0 ]; then
  echo "STEP 0: skipped (START_FROM_STEP=$START_FROM_STEP)."
else
echo "STEP 0: Propedia docking..."
python3 run_propedia_commands_single.py --complex 1AQD_C --propedia-folder Propedia_pdbqt_5Frag --skip-if-exists
fi

# STEP 1: C/N coordinates
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 1 ]; then
  echo "STEP 1: skipped."
else
if [ -f "${VALID_COMBINATIONS_CSV}" ]; then
  echo "STEP 1: valid_fragment_combinations.csv exists; skipping pose CN + aggregate."
else
  echo "STEP 1: C/N coordinates..."
  python3 STELLAR/save_pose_CN_coordinates.py "${DOCKING_BASE}/VS_GN_1AQD_C_Frag1/molecules/"
  python3 STELLAR/save_pose_CN_coordinates.py "${DOCKING_BASE}/VS_GN_1AQD_C_Frag2/molecules/"
  python3 STELLAR/save_pose_CN_coordinates.py "${DOCKING_BASE}/VS_GN_1AQD_C_Frag3/molecules/"
  python3 STELLAR/save_pose_CN_coordinates.py "${DOCKING_BASE}/VS_GN_1AQD_C_Frag4/molecules/"
  python3 STELLAR/save_pose_CN_coordinates.py "${DOCKING_BASE}/VS_GN_1AQD_C_Frag5/molecules/"
  echo "Aggregating..."
  python3 STELLAR/aggregate_gn_pose_coords.py "${DOCKING_BASE}" --case-filter "${CASE_ID}"
fi
fi

# STEP 2: Distance filter
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 2 ]; then
  echo "STEP 2: skipped."
else
if [ -f "${VALID_COMBINATIONS_CSV}" ]; then
  echo "STEP 2: valid_fragment_combinations.csv exists; skipping distance filter."
else
  echo "STEP 2: distance filter..."
  ./STELLAR/filter_fragment_combinations "${DOCKING_BASE}" "${DOCKING_PREFIX}" "${CASE_ID}" --threshold-gn 1.8
fi
fi

# STEP 3: Organize
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 3 ]; then
  echo "STEP 3: skipped."
else
if [ -d "valid_combinations_GN_1AQD_C" ] && [ -n "$(ls -A valid_combinations_GN_1AQD_C 2>/dev/null)" ]; then
  echo "STEP 3: valid_combinations_GN_1AQD_C/ exists; skipping."
else
  echo "STEP 3: organize combinations..."
  python3 STELLAR/organize_valid_combinations.py "${VALID_COMBINATIONS_CSV}" "${DOCKING_BASE}" valid_combinations_GN_1AQD_C GN "${CASE_ID}"
fi
fi

# STEP 4: Overlap
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 4 ]; then
  echo "STEP 4: skipped."
else
if [ -d "valid_combinations_GN_1AQD_C/valid_no_overlap" ] && [ -n "$(ls -A valid_combinations_GN_1AQD_C/valid_no_overlap 2>/dev/null)" ]; then
  echo "STEP 4: valid_no_overlap exists; skipping."
else
  echo "STEP 4: overlap filter..."
  python3 STELLAR/check_overlap_combinations.py valid_combinations_GN_1AQD_C GN 450 0 100
fi
fi

# STEP 5: PDBQT -> MOL2
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 5 ]; then
  echo "STEP 5: skipped."
else
if [ -n "$(ls valid_combinations_GN_1AQD_C/valid_no_overlap/*.mol2 2>/dev/null)" ] || [ -n "$(ls valid_combinations_GN_1AQD_C/valid_no_overlap/combination_*/*.mol2 2>/dev/null)" ]; then
  echo "STEP 5: MOL2 present; skipping conversion."
else
  echo "STEP 5: PDBQT to MOL2..."
  python3 STELLAR/convert_combinations_pdbqt_to_mol2.py GN --base-dir valid_combinations_GN_1AQD_C/valid_no_overlap
fi
fi

# STEP 6: Merge fragments
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 6 ]; then
  echo "STEP 6: skipped."
else
if [ -n "$(ls -d valid_combinations_GN_1AQD_C/valid_no_overlap/combination_*/query_combination_* 2>/dev/null | head -1)" ]; then
  echo "STEP 6: query_combination_* exists; skipping merge."
else
  echo "STEP 6: merge fragments..."
  python3 STELLAR/merge_all_combinations.py GN --base-dir-gn valid_combinations_GN_1AQD_C/valid_no_overlap
fi
fi

# STEP 7: Final structures
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 7 ]; then
  echo "STEP 7: skipped."
else
if [ -d "valid_GN_1AQD_C_final" ] && [ -n "$(ls -A valid_GN_1AQD_C_final/combination_* 2>/dev/null)" ]; then
  echo "STEP 7: valid_GN_1AQD_C_final exists; skipping."
else
  echo "STEP 7: prepare final..."
  python3 STELLAR/prepare_final_combinations.py GN --pdb-file 1aqd_*.pdb --pdb-search-dir "${DOCKING_BASE}/VS_GN_1AQD_C_Frag1/results/best_scores" --base-dir-gn valid_combinations_GN_1AQD_C/valid_no_overlap --output-dir-gn valid_GN_1AQD_C_final
fi
fi

# STEP 8: Fix charges
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 8 ]; then
  echo "STEP 8: skipped."
else
if [ -n "$(ls valid_GN_1AQD_C_final/combination_*/*.top 2>/dev/null)" ]; then
  echo "STEP 8: topologies exist; skipping charge fix."
else
  echo "STEP 8: fix charges..."
  python3 STELLAR/fix_zero_charge_atoms.py GN --base-dir-gn valid_GN_1AQD_C_final --fix-aromatic
fi
fi

# STEP 10: Topologies (step 9 calcium script omitted)
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 10 ]; then
  echo "STEP 10: skipped."
else
if python3 STELLAR/generate_topologies.py GN --base-dir-gn valid_GN_1AQD_C_final --check-only 2>/dev/null; then
  echo "STEP 10: topologies complete; skipping."
else
  echo "STEP 10: generate topologies..."
  python3 STELLAR/generate_topologies.py GN --base-dir-gn valid_GN_1AQD_C_final
fi
fi

# STEP 11: MD (separate sbatch script to limit memory on main job)
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 11 ]; then
  echo "STEP 11: skipped."
else
if [ -n "$(ls valid_GN_1AQD_C_final/combination_*/*.tpr 2>/dev/null)" ] || [ -n "$(ls valid_GN_1AQD_C_final/combination_*/*.edr 2>/dev/null)" ]; then
  echo "STEP 11: MD outputs exist; skipping launch."
elif [ -n "$(find . -maxdepth 4 -type f -path '*/molecules/*.pdb' 2>/dev/null | grep -i "VS_GR.*1AQD" | head -1)" ]; then
  echo "STEP 11: VS_GR PDBs found; skipping MD launch."
else
  echo "STEP 11: submitting MD sub-job..."
  MD_JOBID=$(sbatch --parsable run_md_step_1AQD_C.sh)
  echo "  MD job $MD_JOBID submitted; waiting..."
  while squeue -j $MD_JOBID -h 2>/dev/null | grep -q .; do sleep 30; done
  echo "  MD job finished."
fi
fi

# STEP 12: Fragment RMSD
if [ -n "$START_FROM_STEP" ] && [ "$START_FROM_STEP" -gt 12 ]; then
  echo "STEP 12: skipped."
else
if [ -f "valid_combinations_GN_1AQD_C/valid_no_overlap/rmsd_results.csv" ]; then
  echo "STEP 12: rmsd_results.csv exists; skipping."
else
  echo "STEP 12: fragment RMSD..."
  python3 STELLAR/calculate_rmsd_combinations.py valid_combinations_GN_1AQD_C/valid_no_overlap peptide_pdb_fragments/1AQD_C GN
fi
fi

# STEP 13: Peptide RMSD after MD
if [ -d "md_rmsd_peptides_1AQD_C" ] && [ -n "$(ls -A md_rmsd_peptides_1AQD_C 2>/dev/null)" ]; then
  echo "STEP 13: md_rmsd_peptides exists; skipping."
else
  echo "STEP 13: MD peptide RMSD..."
  python3 STELLAR/calculate_md_rmsd.py --output-csv resultados_rmsd_md_1AQD_C.csv --case 1AQD_C --output-dir md_rmsd_peptides_1AQD_C
fi

# STEP 14
if [ -f "valid_fragment_combinations_GN_no_overlap_1AQD_C.csv" ]; then
  echo "STEP 14: filtered CSV exists; skipping."
else
  echo "STEP 14: filter valid combinations CSV..."
  python3 STELLAR/filter_valid_combinations_csv.py GN --combinations-dir valid_combinations_GN_1AQD_C/valid_no_overlap --input-csv "${VALID_COMBINATIONS_CSV}" --output-csv valid_fragment_combinations_GN_no_overlap_1AQD_C.csv
fi

# STEP 15
if [ -f "score_only_results_1AQD_C.csv" ]; then
  echo "STEP 15: skipping score_only."
else
  echo "STEP 15: score_only..."
  python3 STELLAR/calculate_score_only.py --peptides-dir md_rmsd_peptides_1AQD_C/extracted_peptides --output-csv score_only_results_1AQD_C.csv --singularity-image singularity/metascreener_22.04.simg --singularity-bind /gpfs/
fi

# STEP 15b
if [ -f "mmpbsa_results_1AQD_C.csv" ]; then
  echo "STEP 15b: skipping MM/PBSA."
else
  echo "STEP 15b: MM/PBSA..."
  python3 STELLAR/calculate_mmpbsa.py --case 1AQD_C --type GN --gr-simg singularity/gr.simg --g-mmpbsa-bin ./GROMACS/analyze_results/Simulation_gromacs/analyze_trajectory/extra/g_mmpbsa --output-csv mmpbsa_results_1AQD_C.csv
fi

# STEP 15c
if [ -f "fragment_energies_1AQD_C.csv" ]; then
  echo "STEP 15c: skipping fragment energies."
else
  echo "STEP 15c: fragment energies..."
  python3 STELLAR/calculate_fragment_energies.py --case 1AQD_C --type GN --combinations-file valid_fragment_combinations_GN_no_overlap_1AQD_C.csv --output-csv fragment_energies_1AQD_C.csv
fi

# STEP 16
if [ -f "all_metrics_GN_1AQD_C.csv" ]; then
  echo "STEP 16: all_metrics exists; skipping."
else
  echo "STEP 16: merge_all_metrics..."
  python3 STELLAR/merge_all_metrics.py --prefix-type GN \
    --combinations-file valid_fragment_combinations_GN_no_overlap_1AQD_C.csv \
    --rmsd-fragments-file valid_combinations_GN_1AQD_C/valid_no_overlap/rmsd_results.csv \
    --md-rmsd-file resultados_rmsd_md_1AQD_C.csv \
    --score-only-file score_only_results_1AQD_C.csv \
    --mmpbsa-file mmpbsa_results_1AQD_C.csv \
    --fragment-energies-file fragment_energies_1AQD_C.csv \
    --output all_metrics_GN_1AQD_C.csv
fi

# STEP 17
if [ -d "1AQD_C_results" ] && [ -n "$(ls -A 1AQD_C_results 2>/dev/null)" ]; then
  echo "STEP 17: results folder exists; skipping."
else
  echo "STEP 17: organize results..."
  if [ -f "organize_workflow_results.py" ]; then
    python3 organize_workflow_results.py 1AQD_C --target-dir 1AQD_C_results
  else
    echo "STEP 17: organize_workflow_results.py not found; skipping final organization."
  fi
fi

echo "Workflow finished for 1AQD_C"
date
