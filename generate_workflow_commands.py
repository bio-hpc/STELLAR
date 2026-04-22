#!/usr/bin/env python3
"""
Generate a commented workflow command list for any GN case.

Auto-detects fragment count from VS_GN_{case}_Frag* folders or, if missing, from the
Propedia folder name (e.g. Propedia_pdbqt_6Frag -> 6).

Pipeline scripts live under STELLAR/; launchers stay at project root (see README.md).
"""

import glob
import os
import re
import sys
from pathlib import Path
from typing import Optional, Tuple

# Directory containing this script and the STELLAR/ package (see README merge instructions).
_HERE = Path(__file__).resolve().parent
# Unmerged: .../STELLAR_github/{this file, STELLAR/}. Merged into project root: .../{this file, STELLAR/}.
_PROJECT = _HERE
if str(_PROJECT) not in sys.path:
    sys.path.insert(0, str(_PROJECT))

# Pipeline script prefix (STELLAR_github layout)
ST = "STELLAR"

try:
    from STELLAR.calculate_md_rmsd import find_crystal_structure
except ImportError:
    find_crystal_structure = None


def project_root() -> Path:
    """Directory that contains STELLAR/, MetaScreener/, GROMACS/, etc. (same as this file after merge)."""
    return _PROJECT


def py(script: str) -> str:
    """Return `python3 STELLAR/<script>` (script includes .py)."""
    return f"python3 {ST}/{script}"


def get_case_crystal_pdb(case_name: str, base_dir: str = ".") -> Optional[str]:
    """Return path to crystal (complex) PDB for the case, same logic as calculate_md_rmsd."""
    if find_crystal_structure is None:
        return None
    project_root = os.path.abspath(base_dir) if base_dir != "." else os.getcwd()
    case_upper = case_name.upper() if case_name else None
    crystal_file, _ = find_crystal_structure("complex", None, project_root, case_upper)
    return crystal_file


def detect_case_and_base_dir(input_path: str) -> tuple[str, str]:
    """
    Parse folder path or bare case name.

    Returns:
        (case_name, base_dir) where base_dir is "." if only a case name was given.
    """
    input_path = input_path.rstrip("/")
    if os.path.isdir(input_path) or "/" in input_path or os.path.sep in input_path:
        case_name = os.path.basename(input_path)
        base_dir = os.path.dirname(input_path) or "."
        return case_name, base_dir
    return input_path, "."


def detect_fragment_count(case_name: str, base_dir: str = ".") -> Tuple[Optional[int], Optional[str]]:
    """
    Detect fragment count from VS_GN_{case}_Frag* (or VS_LF_*) under likely locations.

    Returns:
        (max consecutive fragment index or max found, directory containing VS_* folders),
        or (None, None) if none found.
    """
    case_upper = case_name.upper()
    frag_dirs: list[str] = []
    frag_base_dir = None

    search_paths = [
        os.path.join(base_dir, f"{case_upper}_GN"),
        os.path.join(base_dir, f"{case_name.lower()}_GN"),
        base_dir,
    ]
    if "_" in case_name:
        case_short = case_name.split("_")[0]
        search_paths.append(os.path.join(base_dir, case_short))
    search_paths.append(os.path.join(base_dir, case_upper))

    for search_path in search_paths:
        pattern = os.path.join(search_path, f"VS_GN_{case_upper}_Frag*")
        found_dirs = glob.glob(pattern)
        if found_dirs:
            frag_dirs.extend(found_dirs)
            frag_base_dir = search_path
            break

    if not frag_dirs:
        for search_path in search_paths:
            pattern = os.path.join(search_path, f"VS_LF_{case_upper}_Frag*")
            found_dirs = glob.glob(pattern)
            if found_dirs:
                frag_dirs.extend(found_dirs)
                frag_base_dir = search_path
                break

    if not frag_dirs:
        return None, None

    frag_numbers: list[int] = []
    for frag_dir in frag_dirs:
        basename = os.path.basename(frag_dir)
        if os.path.isdir(frag_dir):
            match = re.search(r"_Frag(\d+)", basename)
            if match:
                frag_numbers.append(int(match.group(1)))

    if not frag_numbers:
        return None, None

    max_frag = max(frag_numbers)
    frag_set = set(frag_numbers)
    last_consecutive = 0
    for i in range(1, max_frag + 1):
        if i in frag_set:
            last_consecutive = i
        else:
            break

    if last_consecutive > 0:
        return last_consecutive, frag_base_dir
    return max_frag, frag_base_dir


def get_existing_fragment_numbers(case_name: str, frag_base_dir: Optional[str]) -> list[int]:
    """Sorted list of fragment indices that exist as VS_GN_*_FragN directories."""
    case_upper = case_name.upper()
    frag_numbers: list[int] = []
    if not frag_base_dir or frag_base_dir == ".":
        return frag_numbers

    pattern = os.path.join(frag_base_dir, f"VS_GN_{case_upper}_Frag*")
    for frag_dir in glob.glob(pattern):
        if os.path.isdir(frag_dir):
            match = re.search(r"_Frag(\d+)", os.path.basename(frag_dir))
            if match:
                frag_numbers.append(int(match.group(1)))
    return sorted(frag_numbers)


def detect_propedia_folder(case_name: str, base_dir: str = ".") -> str:
    """Return Propedia_pdbqt_* folder name containing this case with Commands/, or default."""
    case_upper = case_name.upper()
    for propedia_dir in glob.glob(os.path.join(base_dir, "Propedia_pdbqt_*")):
        if not os.path.isdir(propedia_dir):
            continue
        case_dir = os.path.join(propedia_dir, case_upper)
        if os.path.isdir(case_dir) and os.path.isdir(os.path.join(case_dir, "Commands")):
            return os.path.basename(propedia_dir)
    return "Propedia_pdbqt_5Frag"


def extract_fragment_count_from_propedia_name(propedia_folder_name: str) -> Optional[int]:
    """Parse N from Propedia_pdbqt_NFrag."""
    if not propedia_folder_name:
        return None
    match = re.search(r"Propedia_pdbqt_(\d+)Frag", propedia_folder_name, re.IGNORECASE)
    return int(match.group(1)) if match else None


def find_receptor_pdb(case_name: str, base_dir: str = ".", frag_base_dir: Optional[str] = None) -> Optional[str]:
    """Search receptor PDB under VS_* Frag*/results/best_scores and Propedia_Final_mol2."""
    case_upper = case_name.upper()
    case_lower = case_name.lower()
    pdb_code = case_lower.split("_")[0] if "_" in case_lower else case_lower

    if frag_base_dir and frag_base_dir != "." and os.path.isdir(frag_base_dir):
        for prefix in ("VS_GN", "VS_LF"):
            for pat in (f"{pdb_code}_*.pdb", f"{pdb_code}.pdb"):
                p = os.path.join(
                    frag_base_dir, f"{prefix}_{case_upper}_Frag1", "results", "best_scores", pat
                )
                g = [f for f in glob.glob(p) if os.path.isfile(f)]
                if g:
                    return sorted(g)[0]
        for prefix in ("VS_GN", "VS_LF"):
            for pat in (f"{pdb_code}_*.pdb", f"{pdb_code}.pdb"):
                g = glob.glob(
                    os.path.join(
                        frag_base_dir, f"{prefix}_{case_upper}_Frag*", "results", "best_scores", pat
                    )
                )
                g = [f for f in g if os.path.isfile(f)]
                if g:
                    return sorted(g)[0]
        for prefix in ("VS_GN", "VS_LF"):
            for pat in (f"{pdb_code}_*.pdb", f"{pdb_code}.pdb"):
                g = glob.glob(os.path.join(frag_base_dir, f"{prefix}_{case_upper}_Frag*", pat))
                g = [f for f in g if os.path.isfile(f)]
                if g:
                    return sorted(g)[0]
        for pat in (f"{pdb_code}_*.pdb", f"{pdb_code}.pdb"):
            found = [f for f in glob.glob(os.path.join(frag_base_dir, pat)) if os.path.isfile(f)]
            if found:
                return sorted(found)[0]

    receptor_dir = os.path.join(base_dir, "Propedia_Final_mol2", case_upper, "Receptor")
    if os.path.isdir(receptor_dir):
        pdb_files = glob.glob(os.path.join(receptor_dir, "*.pdb"))
        no_pr = [f for f in pdb_files if not f.endswith("_Pr.pdb")]
        if no_pr:
            return sorted(no_pr)[0]
        if pdb_files:
            return sorted(pdb_files)[0]
    return None


def generate_commands(case_name: str, output_file: Optional[str] = None, base_dir: str = ".") -> str:
    """Write workflow commands for case_name to output_file."""
    if output_file is None:
        output_file = f"commands_{case_name}_workflow.txt"

    case_upper = case_name.upper()
    case_lower = case_name.lower()

    propedia_folder = detect_propedia_folder(case_name, base_dir)
    print(f"  Propedia folder: {propedia_folder}")

    crystal_pdb = get_case_crystal_pdb(case_upper, base_dir)
    if crystal_pdb:
        print(f"  Crystal (complex/): {crystal_pdb}")
    else:
        print(f"  Warning: no crystal in complex/ for {case_upper} (step 13 uses --case)")
    crystal_comment = (
        f"# Crystal (complex) for this case: {crystal_pdb}"
        if crystal_pdb
        else f"# Crystal: not found under complex/ (step 13 will use --case {case_upper})"
    )

    num_fragments, frag_base_dir = detect_fragment_count(case_name, base_dir)
    if num_fragments is None:
        print(f"  Warning: no VS_GN_{case_upper}_Frag* folders (step 0 may not have run yet)")
        n_from_propedia = extract_fragment_count_from_propedia_name(propedia_folder)
        if n_from_propedia is not None:
            num_fragments = n_from_propedia
            print(f"  Using {num_fragments} fragments from Propedia folder name: {propedia_folder}")
        else:
            num_fragments = 5
            print("  Using default 5 fragments")
        frag_base_dir = f"{case_upper}_GN"
    else:
        print(f"  Fragments detected: {num_fragments}")
        if frag_base_dir and frag_base_dir != base_dir:
            print(f"  Location: {frag_base_dir}")

    path_prefix = "" if base_dir == "." else f"{base_dir}/"
    frag_path_prefix = "" if frag_base_dir == "." else f"{frag_base_dir}/"

    pdb_file = find_receptor_pdb(case_name, base_dir, frag_base_dir)
    if pdb_file:
        pdb_comment = f"# Receptor PDB found: {pdb_file}"
        pdb_arg = pdb_file
        pdb_search_suffix = ""
    else:
        receptor_path = (
            os.path.join(base_dir, "Propedia_Final_mol2", case_upper, "Receptor")
            if base_dir != "."
            else os.path.join("Propedia_Final_mol2", case_upper, "Receptor")
        )
        pdb_comment = f"# WARNING: no PDB under {receptor_path}"
        pdb_comment += "\n# prepare_final_combinations will search VS_GN_{case}_Frag1/results/best_scores"
        pdb_code = case_lower.split("_")[0] if "_" in case_lower else case_lower
        pdb_arg = f"{pdb_code}_*.pdb"
        best_scores = os.path.join(frag_base_dir or ".", f"VS_GN_{case_upper}_Frag1", "results", "best_scores")
        pdb_search_suffix = f" --pdb-search-dir {best_scores}"

    frag_numbers = get_existing_fragment_numbers(case_name, frag_base_dir)
    if not frag_numbers and num_fragments and num_fragments > 0:
        frag_numbers = list(range(1, num_fragments + 1))

    frag_lines = []
    for i in frag_numbers:
        frag_dir = f"{frag_path_prefix}VS_GN_{case_upper}_Frag{i}"
        frag_lines.append(f"{py('save_pose_CN_coordinates.py')} {frag_dir}/molecules/")
    frag_block = "\n".join(frag_lines)

    if frag_numbers:
        num_fragments = max(frag_numbers)
        print(f"  Fragment dirs: {frag_numbers}")

    aggregate_base = frag_base_dir if frag_base_dir != "." else base_dir
    if aggregate_base == ".":
        aggregate_cmd = f"{py('aggregate_gn_pose_coords.py')} --case-filter {case_upper}"
    else:
        aggregate_cmd = f"{py('aggregate_gn_pose_coords.py')} {aggregate_base} --case-filter {case_upper}"

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

    g_mmpbsa = "./GROMACS/analyze_results/Simulation_gromacs/analyze_trajectory/extra/g_mmpbsa"

    content = f"""# ====================================================================
# GN WORKFLOW FOR {case_upper} ({num_fragments} FRAGMENTS) — COMMANDS ONLY
# ====================================================================
# Base directory: {base_dir if base_dir != "." else "project root"}
{crystal_comment}

# ------------------------------------------------------------------------------
# STEP 0: Propedia docking (sequential SLURM)
# ------------------------------------------------------------------------------
# Runs commands from Propedia for this complex; waits for each job.
# Propedia folder: {propedia_folder}
python3 run_propedia_commands_single.py --complex {case_upper} --propedia-folder {propedia_folder} --skip-if-exists

# ------------------------------------------------------------------------------
# STEP 1: C/N coordinates CSV
# ------------------------------------------------------------------------------
{frag_block}

# Aggregate:
{aggregate_cmd}

# ------------------------------------------------------------------------------
# STEP 2: Distance filter
# ------------------------------------------------------------------------------
# First argument: folder containing VS_GN_* (e.g. {case_upper}_GN). Adjust --threshold-gn if needed.
./{ST}/filter_fragment_combinations {filter_base} VS_GN {case_upper}
# ./{ST}/filter_fragment_combinations {filter_base} VS_GN {case_upper} --threshold-gn 50.0

# ------------------------------------------------------------------------------
# STEP 3: Organize combinations into folders
# ------------------------------------------------------------------------------
# CSV written by filter: {filter_csv}
{py('organize_valid_combinations.py')} {filter_csv} {organize_base} valid_combinations_GN_{case_upper} GN {case_upper}

# ------------------------------------------------------------------------------
# STEP 4: Overlap filter
# ------------------------------------------------------------------------------
{py('check_overlap_combinations.py')} valid_combinations_GN_{case_upper} GN 350 0 100

# ------------------------------------------------------------------------------
# STEP 5: PDBQT -> MOL2 (GN)
# ------------------------------------------------------------------------------
{py('convert_combinations_pdbqt_to_mol2.py')} GN --base-dir valid_combinations_GN_{case_upper}/valid_no_overlap

# ------------------------------------------------------------------------------
# STEP 6: Merge fragments (RDKit via Singularity inside script)
# ------------------------------------------------------------------------------
{py('merge_all_combinations.py')} GN --base-dir-gn valid_combinations_GN_{case_upper}/valid_no_overlap

# ------------------------------------------------------------------------------
# STEP 7: Final structures
# ------------------------------------------------------------------------------
{pdb_comment}
{py('prepare_final_combinations.py')} GN --pdb-file {pdb_arg}{pdb_search_suffix} --base-dir-gn valid_combinations_GN_{case_upper}/valid_no_overlap --output-dir-gn valid_GN_{case_upper}_final

# ------------------------------------------------------------------------------
# STEP 8: Fix zero / aromatic charges in MOL2
# ------------------------------------------------------------------------------
{py('fix_zero_charge_atoms.py')} GN --base-dir-gn valid_GN_{case_upper}_final --fix-aromatic

# ------------------------------------------------------------------------------
# STEP 9: Remove calcium from PDBs (optional; usually commented)
# ------------------------------------------------------------------------------
# {py('remove_calcium_from_pdb.py')} valid_GN_{case_upper}_final

# ------------------------------------------------------------------------------
# STEP 10: GROMACS topologies
# ------------------------------------------------------------------------------
{py('generate_topologies.py')} GN --base-dir-gn valid_GN_{case_upper}_final

# ------------------------------------------------------------------------------
# STEP 11: MD simulations (launcher at project root; uses ./sm.sh and GROMACS/login_node)
# ------------------------------------------------------------------------------
python3 run_md_simulations.py GN --base-dir-gn valid_GN_{case_upper}_final --sequential
# python3 run_md_simulations.py GN --base-dir-gn valid_GN_{case_upper}_final --sequential -o VS -s GR -n 1000

# ------------------------------------------------------------------------------
# STEP 12: Ligand RMSD vs crystal
# ------------------------------------------------------------------------------
{py('calculate_rmsd_combinations.py')} valid_combinations_GN_{case_upper}/valid_no_overlap peptide_pdb_fragments/{case_upper} GN{base_dir_rmsd_arg}

# ------------------------------------------------------------------------------
# STEP 13: Peptide RMSD after MD
# ------------------------------------------------------------------------------
# Crystal: {crystal_pdb if crystal_pdb else "(resolve via complex/ and --case)"}
{py('calculate_md_rmsd.py')} --output-csv resultados_rmsd_md_{case_upper}.csv --case {case_upper} --output-dir md_rmsd_peptides_{case_upper}

# ------------------------------------------------------------------------------
# STEP 14: Filter valid combinations CSV (optional but needed for 15c/16)
# ------------------------------------------------------------------------------
{py('filter_valid_combinations_csv.py')} GN --combinations-dir valid_combinations_GN_{case_upper}/valid_no_overlap --input-csv {filter_csv} --output-csv valid_fragment_combinations_GN_no_overlap_{case_upper}.csv

# ------------------------------------------------------------------------------
# STEP 15: GNINA score_only
# ------------------------------------------------------------------------------
{py('calculate_score_only.py')} --peptides-dir md_rmsd_peptides_{case_upper}/extracted_peptides --output-csv score_only_results_{case_upper}.csv --singularity-image singularity/metascreener_22.04.simg --singularity-bind /gpfs/

# ------------------------------------------------------------------------------
# STEP 15b: MM/PBSA (needs completed VS_GR_* MD)
# ------------------------------------------------------------------------------
{py('calculate_mmpbsa.py')} --case {case_upper} --type GN --gr-simg singularity/gr.simg --g-mmpbsa-bin {g_mmpbsa} --output-csv mmpbsa_results_{case_upper}.csv

# ------------------------------------------------------------------------------
# STEP 15c: Per-fragment energies
# ------------------------------------------------------------------------------
{py('calculate_fragment_energies.py')} --case {case_upper} --type GN --combinations-file valid_fragment_combinations_GN_no_overlap_{case_upper}.csv --output-csv fragment_energies_{case_upper}.csv

# ------------------------------------------------------------------------------
# STEP 16: Merge all metrics
# ------------------------------------------------------------------------------
{py('merge_all_metrics.py')} --prefix-type GN \\
    --combinations-file valid_fragment_combinations_GN_no_overlap_{case_upper}.csv \\
    --rmsd-fragments-file valid_combinations_GN_{case_upper}/valid_no_overlap/rmsd_results.csv \\
    --md-rmsd-file resultados_rmsd_md_{case_upper}.csv \\
    --score-only-file score_only_results_{case_upper}.csv \\
    --mmpbsa-file mmpbsa_results_{case_upper}.csv \\
    --fragment-energies-file fragment_energies_{case_upper}.csv \\
    --output all_metrics_GN_{case_upper}.csv

# ------------------------------------------------------------------------------
# STEP 17: Organize outputs (if script exists at project root)
# ------------------------------------------------------------------------------
python3 organize_workflow_results.py {case_upper} --target-dir {case_upper}_results

# ==============================================================================
# Notes
# ==============================================================================
# Case: {case_upper}
# Data base dir: {base_dir if base_dir != "." else "project root"}
# GN-only pipeline; tune overlap (350) and distance threshold as needed
# Peptide/crystal inputs: {path_prefix}peptide_pdb_fragments/{case_upper}/
# Receptor PDB: {path_prefix}Propedia_Final_mol2/{case_upper}/Receptor/
"""
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(content)

    print(f"Written: {output_file}")
    print(f"  Case: {case_upper}  Fragments: {num_fragments}")
    return output_file


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: python3 generate_workflow_commands.py <CASE_OR_FOLDER> [output.txt]")
        print("Example: python3 generate_workflow_commands.py 1CJR_A")
        print("Example: python3 generate_workflow_commands.py STELLAR_4Frag/1AB9 commands_1AB9.txt")
        sys.exit(1)

    input_path = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else None
    case_name, base_dir = detect_case_and_base_dir(input_path)
    print(f"Case: {case_name}")
    if base_dir != ".":
        print(f"Base directory: {base_dir}")
    os.chdir(project_root())
    generate_commands(case_name, out, base_dir)


if __name__ == "__main__":
    main()
