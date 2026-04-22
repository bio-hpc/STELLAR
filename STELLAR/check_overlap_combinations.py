#!/usr/bin/env python3
"""
Check overlap between poses in valid fragment combinations.
Moves combinations without problematic overlap into a subfolder.
"""

import os
import sys
import subprocess
import shutil
import glob

def get_volumen(ligand1, ligand2=''):
    """
    Get overlap-related volume for two ligands using the overlap binary.
    """
    overlap_cmd = "MetaScreener/external_sw/overlap/overlap {} {} | grep -i volume | tail -1 | awk '{{print $2}}'"
    
    cmd = overlap_cmd.format(ligand1, ligand2)
    try:
        ret = subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL)
        return float(ret.strip())
    except (subprocess.CalledProcessError, ValueError) as e:
        print(f"  Error computing volume for {ligand1} and {ligand2}: {e}")
        return None

def check_overlap(ligand1, ligand2, tolerancia):
    """
    Return True if there is problematic overlap between two ligands, False otherwise.
    """
    vol_combined = get_volumen(ligand1, ligand2)
    vol1 = get_volumen(ligand1)
    vol2 = get_volumen(ligand2)
    
    if vol_combined is None or vol1 is None or vol2 is None:
        return True  # treat error as overlap
    
    # Combined volume below sum minus tolerance => overlap
    return vol_combined < vol1 + vol2 - tolerancia

def get_fragment_files(combo_dir, file_extension):
    """
    List fragment files in a combination folder, sorted by fragment index.
    """
    pattern = os.path.join(combo_dir, f"frag*{file_extension}")
    files = glob.glob(pattern)
    
    def get_frag_number(filename):
        basename = os.path.basename(filename)
        # Index after "frag"
        try:
            parts = basename.split('_')
            for part in parts:
                if part.startswith('frag'):
                    return int(part.replace('frag', ''))
        except:
            return 0
        return 0
    
    files.sort(key=get_frag_number)
    return files

def check_combination_overlap(combo_dir, file_extension, tolerancia, return_score=False):
    """
    Check overlap between consecutive fragments in one combination.
    Returns True if there is no problematic overlap, False if there is overlap.
    If return_score=True, returns (is_valid, overlap_score) where overlap_score is max overlap_vol
    (lower is better); overlap_score is None if invalid.
    """
    fragment_files = get_fragment_files(combo_dir, file_extension)

    if len(fragment_files) < 2:
        return (True, 0.0) if return_score else True

    max_overlap = -float("inf")
    for i in range(len(fragment_files) - 1):
        frag1 = fragment_files[i]
        frag2 = fragment_files[i + 1]

        if not os.path.exists(frag1) or not os.path.exists(frag2):
            print(f"  Warning: file not found in {combo_dir}")
            return (False, None) if return_score else False

        vol_combined = get_volumen(frag1, frag2)
        vol1 = get_volumen(frag1)
        vol2 = get_volumen(frag2)

        if vol_combined is None or vol1 is None or vol2 is None:
            return (False, None) if return_score else False

        overlap_vol = vol1 + vol2 - vol_combined
        if overlap_vol > max_overlap:
            max_overlap = overlap_vol

        has_overlap = overlap_vol > tolerancia  # same as vol_combined < vol1 + vol2 - tolerancia
        if has_overlap:
            print(f"  Overlap between frag{i+1} and frag{i+2}: {overlap_vol:.2f} (tolerance: {tolerancia})")
            return (False, None) if return_score else False

    return (True, max_overlap) if return_score else True

def process_combinations(
    input_dir,
    output_subdir,
    file_extension,
    tolerancia,
    prefix_type,
    max_combinations=None,
    max_valid=None,
):
    """
    Process all combinations and move valid ones (no overlap) into a subfolder.
    If max_valid is set and there are more than max_valid valid combos, only the
    max_valid with lowest overlap scores are kept.

    Args:
        max_combinations: If set, only process the first N combinations (for testing).
        max_valid: If more valid combinations exist, keep only max_valid with lowest overlap.
    """
    if not os.path.exists(input_dir):
        print(f"Error: directory not found: {input_dir}")
        return

    # Output subfolder for valid combinations
    output_dir = os.path.join(input_dir, output_subdir)
    os.makedirs(output_dir, exist_ok=True)

    # combination_* folders
    combo_dirs = [
        d
        for d in os.listdir(input_dir)
        if os.path.isdir(os.path.join(input_dir, d)) and d.startswith("combination_")
    ]
    combo_dirs.sort(key=lambda x: int(x.split("_")[1]) if x.split("_")[1].isdigit() else 0)

    if not combo_dirs:
        print(f"Error: no combination_* folders in {input_dir}")
        print("  Check that step 3 (organize_valid_combinations) ran successfully")
        print("  and step 2 CSV (filter_fragment_combinations) has valid combinations.")
        return

    if max_combinations:
        combo_dirs = combo_dirs[:max_combinations]
        print(f"[TEST MODE] Processing only the first {max_combinations} combinations...")

    if max_valid is not None:
        print(f"If more than {max_valid} are valid, only the {max_valid} with lowest overlap are kept.")

    print(f"Processing {len(combo_dirs)} {prefix_type} combinations...")
    print(f"Overlap tolerance: {tolerancia}")
    print(f"File extension: {file_extension}\n")

    valid_with_scores = []  # (combo_path, combo_name, overlap_score)
    invalid_count = 0

    for i, combo_name in enumerate(combo_dirs, 1):
        combo_path = os.path.join(input_dir, combo_name)

        if i % 100 == 0:
            print(
                f"  Processed {i}/{len(combo_dirs)} combinations... "
                f"(valid: {len(valid_with_scores)}, invalid: {invalid_count})"
            )

        is_valid, overlap_score = check_combination_overlap(
            combo_path, file_extension, tolerancia, return_score=True
        )

        if is_valid and overlap_score is not None:
            valid_with_scores.append((combo_path, combo_name, overlap_score))
        else:
            invalid_count += 1

    # Keep only max_valid with lowest overlap if needed
    num_valid_total = len(valid_with_scores)
    if max_valid is not None and num_valid_total > max_valid:
        valid_with_scores.sort(key=lambda x: x[2])  # lower overlap_score is better
        valid_with_scores = valid_with_scores[:max_valid]
        print(
            f"\n  Keeping {max_valid} combinations with lowest overlap (of {num_valid_total} valid)."
        )

    for combo_path, combo_name, _ in valid_with_scores:
        dest_path = os.path.join(output_dir, combo_name)
        shutil.move(combo_path, dest_path)

    valid_count = len(valid_with_scores)

    print(f"\n✓ Done:")
    print(f"  Valid combinations (no overlap) moved: {valid_count}")
    print(f"  Invalid combinations (overlap): {invalid_count}")
    print(f"  Valid combinations moved to: {output_dir}")

def main():
    if len(sys.argv) < 3:
        print("Usage: python STELLAR/check_overlap_combinations.py <input_dir> <prefix_type> [tolerance] [max_combinations] [max_valid]")
        print("  input_dir: Directory with combination folders (e.g. valid_combinations_GN)")
        print("  prefix_type: 'GN' or 'LF'")
        print("  tolerance: Overlap tolerance (optional; default 225 for GN, 350 for LF)")
        print("  max_combinations: Max combinations to process (optional; for testing)")
        print("  max_valid: If more valid combinations exist, keep only max_valid with lowest overlap (e.g. 100)")
        print("\nExamples:")
        print("  python STELLAR/check_overlap_combinations.py valid_combinations_GN GN")
        print("  python STELLAR/check_overlap_combinations.py valid_combinations_GN GN 225  # tolerance 225")
        print("  python STELLAR/check_overlap_combinations.py valid_combinations_GN GN 225 0 100  # at most 100 valid, lowest overlap")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    prefix_type = sys.argv[2].upper()
    
    if prefix_type not in ['GN', 'LF']:
        print("Error: prefix_type must be 'GN' or 'LF'")
        sys.exit(1)
    
    # Default tolerances
    default_tolerances = {'GN': 225, 'LF': 350}
    tolerancia = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3].isdigit() else default_tolerances[prefix_type]
    
    # Max combinations to process (testing); 0 = no limit
    max_combinations = None
    if len(sys.argv) > 4:
        try:
            n = int(sys.argv[4])
            max_combinations = n if n > 0 else None
        except ValueError:
            pass

    # Optional cap: keep only max_valid with lowest overlap
    max_valid = None
    if len(sys.argv) > 5:
        try:
            max_valid = int(sys.argv[5])
        except ValueError:
            pass

    # Extensión de archivo según el tipo
    file_extension = '.pdbqt' if prefix_type == 'GN' else '.mol2'
    
    # Output subfolder name
    output_subdir = 'valid_no_overlap'
    
    # Require overlap binary
    overlap_binary = "MetaScreener/external_sw/overlap/overlap"
    if not os.path.exists(overlap_binary):
        print(f"Error: overlap binary not found at {overlap_binary}")
        print("Ensure MetaScreener/external_sw/overlap/overlap exists")
        sys.exit(1)
    
    process_combinations(
        input_dir, output_subdir, file_extension, tolerancia, prefix_type,
        max_combinations=max_combinations,
        max_valid=max_valid,
    )

if __name__ == '__main__':
    main()

