#!/usr/bin/env python3
"""
Compute per-fragment RMSD for valid fragment combinations against the crystal structure.
"""

import os
import csv
import glob
import subprocess
import sys
import tempfile
from collections import Counter

def calculate_rmsd(reference_file, query_file, singularity_image="singularity/new_ms.simg"):
    """
    Compute RMSD with obrms inside Singularity.

    Args:
        reference_file: Reference structure (crystal)
        query_file: Structure to compare
        singularity_image: Singularity image path

    Returns:
        RMSD value (float) or None if error
    """
    cmd = f'singularity exec {singularity_image} obrms "{reference_file}" "{query_file}"'

    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            check=True,
            cwd=os.getcwd()
        )
        # obrms usually prints RMSD on stdout
        rmsd_str = result.stdout.strip()
        try:
            return float(rmsd_str)
        except ValueError:
            # If not a plain number, try to parse it
            # obrms may print text like "RMSD: 1.23"
            for word in rmsd_str.split():
                try:
                    return float(word)
                except ValueError:
                    continue
            return None
    except subprocess.CalledProcessError as e:
        print(f"  Error computing RMSD: {e.stderr}")
        return None
    except Exception as e:
        print(f"  Unexpected error: {e}")
        return None

def process_combinations(combinations_dir, crystal_dir, output_csv, prefix_type="GN", singularity_image="singularity/new_ms.simg"):
    """
    Process all combinations and compute RMSD per fragment.
    Runs all obrms calls in a single Singularity session.

    Args:
        combinations_dir: Directory containing combination_* folders
        crystal_dir: Directory with crystal fragment PDBs
        output_csv: Output CSV path
        prefix_type: 'GN' or 'LF'
        singularity_image: Singularity image path
    """
    # Infer case name from crystal_dir (e.g. peptide_pdb_fragments/4WVI_D -> 4WVI_D)
    case_name = os.path.basename(crystal_dir.rstrip('/'))
    if not case_name:
        # If crystal_dir ends in /, use parent basename
        case_name = os.path.basename(os.path.dirname(crystal_dir.rstrip('/')))

    # Locate crystal PDBs (names may vary)
    crystal_files = {}
    max_fragments = 20  # search up to fragment 20

    for frag_num in range(1, max_fragments + 1):
        patterns = [
            os.path.join(crystal_dir, "Fragments", f"Fragmento{frag_num}", f"{case_name}_Frag{frag_num}.pdb"),
            os.path.join(crystal_dir, "Fragments", f"Fragmento{frag_num}", f"{case_name.upper()}_Frag{frag_num}.pdb"),
            os.path.join(crystal_dir, "Fragments", f"Fragmento{frag_num}", f"{case_name.lower()}_Frag{frag_num}.pdb"),
            os.path.join(crystal_dir, "Fragments", f"Fragmento{frag_num}", "*.pdb"),  # any PDB
        ]

        for pattern in patterns:
            matches = glob.glob(pattern)
            if matches:
                crystal_files[frag_num] = os.path.abspath(matches[0])
                break
        # Do not stop if a fragment is missing: keep scanning up to max_fragments
        # (allows Frag1..Frag5 without requiring every index to exist)

    if not crystal_files:
        print(f"Error: no crystal PDBs found under {crystal_dir}")
        return False

    for frag_num, crystal_file in crystal_files.items():
        if not os.path.exists(crystal_file):
            print(f"Error: crystal PDB missing for fragment {frag_num}: {crystal_file}")
            return False
        print(f"Fragment {frag_num} crystal: {crystal_file}")

    combo_dirs = sorted(glob.glob(os.path.join(combinations_dir, "combination_*")))

    if not combo_dirs:
        print(f"Error: no combination_* folders in {combinations_dir}")
        print("  Possible causes:")
        print("  - Step 2 (filter_fragment_combinations): no distance-valid combinations.")
        print("  - Step 3 (organize_valid_combinations): no combination_* folders created.")
        print("  - Step 4 (check_overlap_combinations): all combinations overlapped; none moved to valid_no_overlap.")
        print("  Check the step 2 CSV and the logs from steps 3 and 4.")
        return False

    print(f"\nProcessing {len(combo_dirs)} combinations...")
    print(f"Type: {prefix_type}")
    print(f"Batch mode: all obrms runs in one Singularity session\n")

    # Build work list (.mol2 first, then .pdbqt, etc.)
    calculations = []  # (combo_id, frag_num, crystal_file, query_file)

    for combo_dir in combo_dirs:
        combo_name = os.path.basename(combo_dir)
        combo_id = combo_name.replace("combination_", "")

        for frag_num in sorted(crystal_files.keys()):
            if frag_num not in crystal_files or not crystal_files[frag_num]:
                continue
            query_file = None
            # Include .pdb if poses are PDB; try frag and Frag prefixes
            for name_prefix in (f"frag{frag_num}_pose", f"Frag{frag_num}_pose"):
                for ext in (".pdbqt", ".mol2", ".pdb", ".PDB"):
                    pattern = os.path.join(combo_dir, f"{name_prefix}*{ext}")
                    matches = glob.glob(pattern)
                    if matches:
                        query_file = os.path.abspath(matches[0])
                        break
                if query_file:
                    break
            if query_file:
                calculations.append((combo_id, frag_num, crystal_files[frag_num], query_file))

    frag_counts = Counter(frag_num for _, frag_num, _, _ in calculations)
    print(f"Total RMSD evaluations: {len(calculations)}")
    for fn in sorted(crystal_files.keys()):
        print(f"  Fragment {fn}: {frag_counts.get(fn, 0)} evaluations")

    temp_output = os.path.join(combinations_dir, "rmsd_temp_output.txt")
    temp_output_abs = os.path.abspath(temp_output)

    script_content = "#!/bin/bash\n"
    script_content += "# Generated script for batch RMSD\n"
    script_content += f'OUTPUT_FILE="{temp_output_abs}"\n'
    script_content += '> "$OUTPUT_FILE"\n\n'  # truncate output file

    for combo_id, frag_num, crystal_file, query_file in calculations:
        script_content += f'echo "{combo_id}|{frag_num}|$(obrms "{crystal_file}" "{query_file}" 2>/dev/null || echo "ERROR")" >> "$OUTPUT_FILE"\n'

    with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as f:
        f.write(script_content)
        temp_script = f.name

    os.chmod(temp_script, 0o755)

    print("Running RMSD batch in Singularity...")
    cmd = f'singularity exec {singularity_image} bash "{temp_script}"'

    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True, cwd=os.getcwd())
        print("  Batch completed")
    except subprocess.CalledProcessError as e:
        print(f"Error running batch script: {e.stderr}")
        os.unlink(temp_script)
        return False

    rmsd_results = {}  # {(combo_id, frag_num): rmsd_value}

    if os.path.exists(temp_output):
        with open(temp_output, 'r') as f:
            for line in f:
                line = line.strip()
                if '|' in line:
                    parts = line.split('|')
                    if len(parts) == 3:
                        combo_id = parts[0]
                        frag_num = int(parts[1])
                        rmsd_str = parts[2].strip()

                        if rmsd_str != "ERROR" and rmsd_str:
                            try:
                                rmsd_value = float(rmsd_str)
                                rmsd_results[(combo_id, frag_num)] = rmsd_value
                            except ValueError:
                                for word in rmsd_str.split():
                                    try:
                                        rmsd_value = float(word)
                                        rmsd_results[(combo_id, frag_num)] = rmsd_value
                                        break
                                    except ValueError:
                                        continue

    os.unlink(temp_script)
    if os.path.exists(temp_output):
        os.unlink(temp_output)

    all_combo_ids = []
    for combo_dir in combo_dirs:
        combo_name = os.path.basename(combo_dir)
        combo_id = combo_name.replace("combination_", "")
        all_combo_ids.append(combo_id)

    results = []
    max_fragments = max(crystal_files.keys()) if crystal_files else 0

    for combo_id in sorted(all_combo_ids, key=lambda x: int(x) if x.isdigit() else 0):
        result_row = {'combination_id': combo_id}
        for frag_num in range(1, max_fragments + 1):
            result_row[f'rmsd_frag{frag_num}'] = rmsd_results.get((combo_id, frag_num), '')

        for key in result_row:
            if result_row[key] is None:
                result_row[key] = ''
            elif str(result_row[key]).lower() == 'inf' or str(result_row[key]).lower() == 'error':
                result_row[key] = ''  # treat inf/error as empty
        results.append(result_row)

    fieldnames = ['combination_id'] + [f'rmsd_frag{i}' for i in range(1, max_fragments + 1)]
    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"\n✓ Done:")
    print(f"  Combinations written: {len(results)}")
    print(f"  Successful RMSD values: {len(rmsd_results)}")
    print(f"  CSV: {output_csv}")

    valid_rmsd = [r for r in results if any(r.get(f'rmsd_frag{i}', '') != '' for i in range(1, max_fragments + 1))]
    print(f"  Combinations with at least one RMSD: {len(valid_rmsd)}")

    return True

def main():
    # Optional --base-dir: project root containing docking and crystal data
    base_dir = None
    singularity_image = "singularity/new_ms.simg"
    pos_args = []
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == "--base-dir" and i + 1 < len(sys.argv):
            base_dir = sys.argv[i + 1]
            i += 2
            continue
        if sys.argv[i] == "--singularity-image" and i + 1 < len(sys.argv):
            singularity_image = sys.argv[i + 1]
            i += 2
            continue
        pos_args.append(sys.argv[i])
        i += 1

    if len(pos_args) < 2:
        print("Usage: python STELLAR/calculate_rmsd_combinations.py <combinations_dir> <crystal_dir> [prefix_type] [--base-dir DIR] [--singularity-image IMG]")
        print("  combinations_dir: e.g. valid_combinations_GN/valid_no_overlap")
        print("  crystal_dir: e.g. peptide_pdb_fragments/1A1M_C")
        print("  prefix_type: 'GN' or 'LF' (default: GN)")
        print("  --base-dir: project root with docking and crystal data (e.g. STELLAR_4Frag)")
        print("  --singularity-image: image path (default: singularity/new_ms.simg)")
        sys.exit(1)

    combinations_dir = pos_args[0]
    crystal_dir = pos_args[1]
    prefix_type = pos_args[2] if len(pos_args) > 2 else "GN"

    if prefix_type not in ['GN', 'LF']:
        print("Error: prefix_type must be 'GN' or 'LF'")
        sys.exit(1)

    if base_dir and os.path.isdir(base_dir):
        for candidate in [
            os.path.join(base_dir, crystal_dir),
            os.path.join(base_dir, "peptide_pdb_fragments", os.path.basename(crystal_dir.rstrip("/"))),
        ]:
            c = os.path.normpath(candidate)
            if os.path.isdir(c):
                crystal_dir = os.path.abspath(c)
                print(f"Using crystal directory from base_dir: {crystal_dir}")
                break

    if not os.path.exists(singularity_image):
        legacy_image = "new_ms.simg"
        if os.path.exists(legacy_image):
            print(f"Using legacy image path: {legacy_image}")
            singularity_image = legacy_image
        else:
            print(f"Error: Singularity image not found: {singularity_image}")
            sys.exit(1)

    output_csv = os.path.join(combinations_dir, "rmsd_results.csv")

    print(f"Computing RMSD for {prefix_type} combinations")
    print(f"Combinations directory: {combinations_dir}")
    print(f"Crystal directory: {crystal_dir}")
    if base_dir:
        print(f"Base dir (docking/crystal): {base_dir}")
    print(f"Output CSV: {output_csv}\n")

    success = process_combinations(combinations_dir, crystal_dir, output_csv, prefix_type, singularity_image)

    if not success:
        sys.exit(1)

if __name__ == '__main__':
    main()
