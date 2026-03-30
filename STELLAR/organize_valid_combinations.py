#!/usr/bin/env python3
"""
Organize valid fragment combinations into separate folders.
Each folder holds the pose files for one valid combination.
"""

import os
import sys
import csv
import shutil
import glob


def find_molecule_file(fragment_dir, pose_number, file_type='pdbqt'):
    """
    Find the molecule file for a given pose number.

    Args:
        fragment_dir: Fragment directory (e.g. VS_GN_1A1M_C_Frag1).
        pose_number: Pose number from the CSV.
        file_type: 'pdbqt' for VS_GN, 'mol2' for VS_LF.

    Returns:
        Path to the file, or None.
    """
    molecules_dir = os.path.join(fragment_dir, 'molecules')
    if not os.path.exists(molecules_dir):
        return None
    
    # Match files by pattern
    pattern = os.path.join(molecules_dir, f'*.{file_type}')
    all_files = glob.glob(pattern)
    
    if not all_files:
        return None
    
    # Try numbering formats in preference order
    # Format 1: _pose_number (e.g. _1, _2, _10) — most common
    pattern_file = os.path.join(molecules_dir, f'*_{pose_number}.{file_type}')
    matches = glob.glob(pattern_file)
    if matches:
        return matches[0]
    
    # Format 2: _pose_number with 2-digit padding (e.g. _01, _02, _10)
    pattern_file = os.path.join(molecules_dir, f'*_{pose_number:02d}.{file_type}')
    matches = glob.glob(pattern_file)
    if matches:
        return matches[0]
    
    # Format 3: no trailing number (for pose 1, may be the base file)
    if pose_number == 1:
        # Base prefix from first file
        base_name = os.path.basename(all_files[0])
        # If last segment after _ is not a digit, treat as base name
        if '_' in base_name:
            parts = base_name.rsplit('_', 1)
            last_part = parts[1].replace(f'.{file_type}', '')
            # Last segment not a digit: may be base file
            if not last_part.isdigit():
                # File ending exactly with prefix + extension
                base_prefix = parts[0]
                pattern_file = os.path.join(molecules_dir, f'{base_prefix}.{file_type}')
                if os.path.exists(pattern_file):
                    return pattern_file
        
        # Also try file with no trailing number (prefix only)
        # Common prefix across files
        all_basenames = [os.path.basename(f) for f in all_files]
        # Common prefix up to last _
        if all_basenames:
            first_name = all_basenames[0]
            if '_' in first_name:
                base_prefix = first_name.rsplit('_', 1)[0]
                pattern_file = os.path.join(molecules_dir, f'{base_prefix}.{file_type}')
                if os.path.exists(pattern_file):
                    return pattern_file
    
    # Fallback: match by sorted index (less reliable)
    all_files_sorted = sorted(all_files)
    # Files may be ordered; padding vs no-padding must be handled carefully
    if 1 <= pose_number <= len(all_files_sorted):
        # Unnumbered file may be pose 0 or 1; assume sorted list: pose N -> index N-1
        return all_files_sorted[pose_number - 1]
    
    return None


def get_fragment_info(base_dir, fragment_name, file_type):
    """
    Return paths and naming prefix for a fragment.

    Returns:
        (fragment_dir, molecules_dir, prefix_pattern)
    """
    fragment_dir = os.path.join(base_dir, fragment_name)
    molecules_dir = os.path.join(fragment_dir, 'molecules')
    
    if not os.path.exists(molecules_dir):
        return None, None, None
    
    # Sample file to infer prefix
    pattern = os.path.join(molecules_dir, f'*.{file_type}')
    example_files = glob.glob(pattern)
    
    if not example_files:
        return fragment_dir, molecules_dir, None
    
    # Base prefix
    example_name = os.path.basename(example_files[0])
    # Prefix: up to last _ before digits, or whole stem
    if '_' in example_name:
        parts = example_name.rsplit('_', 1)
        if parts[1].replace(f'.{file_type}', '').isdigit():
            prefix = parts[0]
        else:
            prefix = example_name.replace(f'.{file_type}', '')
    else:
        prefix = example_name.replace(f'.{file_type}', '')
    
    return fragment_dir, molecules_dir, prefix

def extract_case_from_csv(csv_file):
    """
    Extract case code from the CSV filename.
    Examples:
    - valid_fragment_combinations_GN_2W10_C_final.csv -> 2W10_C
    - valid_fragment_combinations_LF_1A1M_C.csv -> 1A1M_C
    """
    import re
    basename = os.path.basename(csv_file)
    patterns = [
        r'([0-9A-Z]+_[A-Z])',  # e.g. 2W10_C
        r'([0-9][A-Z0-9]+_[A-Z])',  # alternate
    ]
    for pattern in patterns:
        match = re.search(pattern, basename)
        if match:
            return match.group(1)
    return None

def organize_combinations(csv_file, base_dir, output_dir, prefix_type='GN', case_filter=None):
    """
    Copy valid combinations into per-combination folders.

    Args:
        csv_file: CSV with valid combinations.
        base_dir: Root containing VS_GN or VS_LF fragment folders.
        output_dir: Where combination_* folders are created.
        prefix_type: 'GN' or 'LF'.
        case_filter: Case code (e.g. '2W10_C') to restrict fragments; if None,
            tries to infer from the CSV filename.
    """
    file_type = 'pdbqt' if prefix_type == 'GN' else 'mol2'
    prefix_str = f'VS_{prefix_type}'
    
    # Auto-detect case if not given
    if case_filter is None:
        detected_case = extract_case_from_csv(csv_file)
        if detected_case:
            case_filter = detected_case
            print(f"ℹ Case inferred from CSV filename: {case_filter}")
        else:
            print("⚠ Warning: could not infer case from CSV. Using ALL available fragment folders.")
            print("  This may mix different cases. Pass case_filter to avoid that.")
    
    # Output root
    os.makedirs(output_dir, exist_ok=True)
    
    # Read combinations CSV
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        combinations = list(reader)
    
    print(f"Processing {len(combinations)} combinations...")
    
    # Fragment paths and prefixes
    fragment_dirs = {}
    fragment_prefixes = {}
    
    # Discover fragment folders
    all_fragment_folders = sorted([d for d in os.listdir(base_dir) 
                                   if d.startswith(prefix_str) and 'Frag' in d])
    
    # Strict case filter
    if case_filter:
        # Normalize filter (2W10_C, 2w10_c, etc.)
        case_filter_upper = case_filter.upper()
        case_filter_lower = case_filter.lower()
        # Filtrar de manera estricta: el caso debe aparecer en el nombre de la carpeta
        fragment_folders = [d for d in all_fragment_folders 
                           if case_filter_upper in d.upper() or case_filter_lower in d.lower()]
        
        if not fragment_folders:
            print(f"✗ ERROR: no fragment folders match filter '{case_filter}'")
            print(f"  Sample folders: {all_fragment_folders[:10]}...")
            print(f"  Total fragment folders: {len(all_fragment_folders)}")
            sys.exit(1)
        
        # Report excluded folders
        excluded_folders = [d for d in all_fragment_folders if d not in fragment_folders]
        if excluded_folders:
            print(f"✓ Filter: {len(fragment_folders)} folders for {case_filter}, {len(excluded_folders)} excluded")
            other_cases = set()
            for folder in excluded_folders[:5]:
                # Guess case labels from excluded folder names
                parts = folder.split('_')
                if len(parts) >= 3:
                    possible_case = '_'.join(parts[2:4]) if len(parts) >= 4 else parts[2]
                    if possible_case.upper() != case_filter_upper:
                        other_cases.add(possible_case)
            if other_cases:
                print(f"  Excluded (other cases): {', '.join(sorted(other_cases)[:5])}")
    else:
        fragment_folders = all_fragment_folders
        print(f"⚠ WARNING: no case filter; using {len(fragment_folders)} fragment folders.")
        print("  This may mix different cases. Pass case_filter to avoid that.")
    
    for frag_folder in fragment_folders:
        frag_num = int(''.join(filter(str.isdigit, frag_folder.split('Frag')[-1])))
        frag_dir, mol_dir, prefix = get_fragment_info(base_dir, frag_folder, file_type)
        if frag_dir:
            fragment_dirs[frag_num] = frag_dir
            fragment_prefixes[frag_num] = prefix
            print(f"Fragment {frag_num}: {frag_folder} (prefix: {prefix})")
    
    # Process each combination
    successful = 0
    failed = 0
    
    for i, combo in enumerate(combinations, 1):
        combo_id = combo.get('combination_id', str(i))
        combo_dir = os.path.join(output_dir, f'combination_{combo_id}')
        os.makedirs(combo_dir, exist_ok=True)
        
        all_found = True
        
        # Copy each fragment's pose file
        for frag_num in sorted(fragment_dirs.keys()):
            pose_key = f'frag{frag_num}_pose'
            if pose_key not in combo:
                continue
            
            pose_number = int(combo[pose_key])
            frag_dir = fragment_dirs[frag_num]
            
            # Locate molecule file
            mol_file = find_molecule_file(frag_dir, pose_number, file_type)
            
            if mol_file and os.path.exists(mol_file):
                # Extra check: filename should contain case when filtering
                if case_filter:
                    case_filter_upper = case_filter.upper()
                    filename = os.path.basename(mol_file)
                    if case_filter_upper not in filename.upper():
                        print(f"  ✗ ERROR: file {filename} lacks '{case_filter}' in name but folder was case-filtered")
                        print(f"    Folder: {frag_dir}")
                        all_found = False
                        continue
                
                # Copy into combination folder
                dest_file = os.path.join(combo_dir, 
                                       f'frag{frag_num}_pose{pose_number}.{file_type}')
                shutil.copy2(mol_file, dest_file)
            else:
                print(f"  Warning: no file for {frag_dir}, pose {pose_number}")
                all_found = False
        
        if all_found:
            successful += 1
        else:
            failed += 1
        
        if i % 100 == 0:
            print(f"  Processed {i}/{len(combinations)} combinations...")
    
    print(f"\nDone:")
    print(f"  Successful combinations: {successful}")
    print(f"  Failed combinations: {failed}")
    print(f"  Output directory: {output_dir}")

def main():
    if len(sys.argv) < 3:
        print("Usage: python STELLAR/organize_valid_combinations.py <csv_file> <base_dir> <output_dir> <prefix_type> [case_filter]")
        print("Example: python STELLAR/organize_valid_combinations.py valid_combinations.csv . valid_combinations_GN GN")
        print("With case filter: python STELLAR/organize_valid_combinations.py valid_combinations.csv . valid_combinations_GN GN 2W10_C")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    base_dir = sys.argv[2]
    output_dir = sys.argv[3]
    prefix_type = sys.argv[4] if len(sys.argv) > 4 else 'GN'
    case_filter = sys.argv[5] if len(sys.argv) > 5 else None
    
    if prefix_type not in ['GN', 'LF']:
        print("Error: prefix_type must be 'GN' or 'LF'")
        sys.exit(1)
    
    if not os.path.exists(csv_file):
        print(f"Error: file not found: {csv_file}")
        sys.exit(1)
    
    if not os.path.exists(base_dir):
        print(f"Error: directory not found: {base_dir}")
        sys.exit(1)
    
    organize_combinations(csv_file, base_dir, output_dir, prefix_type, case_filter)

if __name__ == '__main__':
    main()

