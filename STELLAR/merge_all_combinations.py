#!/usr/bin/env python3
"""
Merge fragments for every valid GN and LF combination.

For each combination:
1. Merge fragments with relax_merge_mol2.py
2. Fix charge drift with fix_charge_drift.py

Processes:
- valid_combinations_GN/valid_no_overlap/combination_*
- valid_combinations_LF/valid_no_overlap/combination_*
"""

import os
import glob
import subprocess
import sys
import argparse
from pathlib import Path

_STELLAR_DIR = Path(__file__).resolve().parent


def find_singularity_image():
    """
    Find an available Singularity image with RDKit.

    Returns:
        Path to the image, or None.
    """
    possible_images = [
        "singularity/STELLAR.simg",
        "singularity/BIOFRAGMEN.simg",
        "singularity/metascreener_22.04.simg",
        "singularity/metascreener.simg",
    ]
    
    for img_path in possible_images:
        if os.path.exists(img_path):
            # Check RDKit import inside image
            try:
                result = subprocess.run(
                    ["singularity", "exec", img_path, "python3", "-c", "import rdkit"],
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                if result.returncode == 0:
                    return img_path
            except:
                continue
    
    return None


def find_fragment_files(combo_dir):
    """
    Find fragment .mol2 files in a combination folder.
    Fragment count is auto-detected (not limited to 3).

    Args:
        combo_dir: Combination directory.

    Returns:
        {1: path1, 2: path2, ...} or None if invalid.
    """
    import re
    
    def extract_fragment_number(filename):
        """Parse fragment index from filename."""
        patterns = [
            r'[Ff]rag(?:ment)?(\d+)',
            r'[Ff](\d+)',
        ]
        for pattern in patterns:
            match = re.search(pattern, filename)
            if match:
                return int(match.group(1))
        return None
    
    # Collect .mol2 files
    mol2_files = glob.glob(os.path.join(combo_dir, "*.mol2"))
    
    frag_files = {}
    for mol2_file in mol2_files:
        frag_num = extract_fragment_number(os.path.basename(mol2_file))
        if frag_num and frag_num >= 1:  # no upper limit
            frag_files[frag_num] = mol2_file
    
    # Require at least one fragment and consecutive numbering from 1
    if len(frag_files) > 0:
        max_frag = max(frag_files.keys())
        # Fragments 1..max_frag must all be present
        expected_frags = set(range(1, max_frag + 1))
        if expected_frags.issubset(set(frag_files.keys())):
            return frag_files
    
    return None


def merge_fragments(combo_dir, output_file, dry_run=False, singularity_img=None):
    """
    Merge fragments with relax_merge_mol2.py.

    Args:
        combo_dir: Directory with fragment .mol2 files.
        output_file: Output path.
        dry_run: If True, only print the command.
        singularity_img: Singularity image path (optional).

    Returns:
        True on success.
    """
    script_path = os.path.abspath("relax_merge_mol2.py")
    combo_dir_abs = os.path.abspath(combo_dir)
    output_abs = os.path.abspath(output_file)
    
    # Auto-detect Singularity image if not given
    if singularity_img is None:
        singularity_img = find_singularity_image()
    
    if singularity_img and os.path.exists(singularity_img):
        # Run inside Singularity
        cmd = [
            "singularity", "exec",
            singularity_img,
            "python3",
            script_path,
            "--folder", combo_dir_abs,
            "-o", output_abs
        ]
    else:
        # Host Python
        cmd = [
            sys.executable,
            script_path,
            "--folder", combo_dir_abs,
            "-o", output_abs
        ]
    
    if dry_run:
        print(f"  [DRY RUN] {' '.join(cmd)}")
        return True
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False
        )
        
        if result.returncode == 0:
            if os.path.exists(output_abs):
                return True
            else:
                print(f"    ✗ Error: output file was not created: {output_abs}")
                if result.stderr:
                    print(f"      stderr: {result.stderr[-500:]}")
                return False
        else:
            print(f"    ✗ Error merging fragments")
            if result.stderr:
                print(f"      stderr: {result.stderr[-500:]}")
            return False
    
    except Exception as e:
        print(f"    ✗ Exception: {e}")
        return False


def fix_charge_drift(input_file, output_file, target=2.0, tolerance=0.008, dry_run=False):
    """
    Fix charge drift with fix_charge_drift.py.

    Args:
        input_file: Input .mol2.
        output_file: Output .mol2.
        target: Target total charge (default: 2.0).
        tolerance: Max tolerance (default: 0.008).
        dry_run: If True, only print the command.

    Returns:
        True on success.
    """
    script_path = str(_STELLAR_DIR / "fix_charge_drift.py")
    input_abs = os.path.abspath(input_file)
    output_abs = os.path.abspath(output_file)
    
    cmd = [
        sys.executable,
        script_path,
        input_abs,
        "-o", output_abs,
        "--target", str(target),
        "--tolerance", str(tolerance)
    ]
    
    if dry_run:
        print(f"  [DRY RUN] {' '.join(cmd)}")
        return True
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False
        )
        
        if result.returncode == 0:
            if os.path.exists(output_abs):
                return True
            else:
                print(f"    ✗ Error: output file was not created: {output_abs}")
                if result.stderr:
                    print(f"      stderr: {result.stderr[-500:]}")
                return False
        else:
            print(f"    ✗ Error fixing charge drift")
            if result.stderr:
                print(f"      stderr: {result.stderr[-500:]}")
            return False
    
    except Exception as e:
        print(f"    ✗ Exception: {e}")
        return False


def process_combination(combo_dir, target_charge=2.0, tolerance=0.008, dry_run=False, singularity_img=None):
    """
    Merge fragments and fix charge drift for one combination.

    Args:
        combo_dir: Combination directory.
        target_charge: Target charge for fix_charge_drift.
        tolerance: Tolerance for fix_charge_drift.
        dry_run: If True, only print actions.
        singularity_img: Singularity image path (optional).

    Returns:
        (success, error_message or None)
    """
    combo_name = os.path.basename(combo_dir)
    
    # Verificar que existen los archivos de fragmentos
    frag_files = find_fragment_files(combo_dir)
    if not frag_files:
        num_found = len(frag_files) if frag_files else 0
        return False, f"No valid fragments in {combo_name} (found: {num_found})"
    
    num_fragments = len(frag_files)
    print(f"  → Found {num_fragments} fragments: {sorted(frag_files.keys())}")
    
    # Output paths
    merged_file = os.path.join(combo_dir, "fragmento_final.mol2")
    final_file = os.path.join(combo_dir, "fragmento_final_charge_drift.mol2")
    
    if os.path.exists(final_file) and not dry_run:
        return True, None  # already done
    
    if not merge_fragments(combo_dir, merged_file, dry_run, singularity_img):
        return False, f"Error merging fragments in {combo_name}"
    
    if not fix_charge_drift(merged_file, final_file, target_charge, tolerance, dry_run):
        return False, f"Error fixing charge drift in {combo_name}"
    
    return True, None


def process_all_combinations(base_dir, prefix_type="GN", target_charge=2.0, tolerance=0.008, dry_run=False, singularity_img=None):
    """
    Process all combinations under base_dir.

    Args:
        base_dir: Root with combination_* folders.
        prefix_type: 'GN' or 'LF' (informational).
        target_charge: Target charge for fix_charge_drift.
        tolerance: Tolerance for fix_charge_drift.
        dry_run: If True, only print actions.
        singularity_img: Singularity image path (optional).

    Returns:
        (total_processed, total_errors, errors_list)
    """
    combo_dirs = sorted(glob.glob(os.path.join(base_dir, "combination_*")))
    
    if not combo_dirs:
        print(f"⚠ No combination_* folders in {base_dir}")
        return 0, 0, []
    
    print(f"Found {len(combo_dirs)} combinations")
    if dry_run:
        print("⚠ DRY RUN — no files will be modified")
    if singularity_img:
        print(f"📦 Singularity image: {singularity_img}")
    print()
    
    total_processed = 0
    total_errors = 0
    errors_list = []
    
    for i, combo_dir in enumerate(combo_dirs, 1):
        combo_name = os.path.basename(combo_dir)
        print(f"[{i}/{len(combo_dirs)}] {combo_name}")
        
        success, error = process_combination(combo_dir, target_charge, tolerance, dry_run, singularity_img)
        
        if success:
            total_processed += 1
            if not dry_run:
                print(f"  ✓ Processed successfully")
        else:
            total_errors += 1
            errors_list.append((combo_name, error))
            print(f"  ✗ {error}")
        
        if (i % 10 == 0) and not dry_run:
            print(f"  Progress: {total_processed} ok, {total_errors} errors...")
    
    return total_processed, total_errors, errors_list


def main():
    parser = argparse.ArgumentParser(
        description="Merge fragments for all valid combinations"
    )
    parser.add_argument(
        'type',
        choices=['GN', 'LF', 'all'],
        help='Combination type: GN, LF, or all (both)'
    )
    parser.add_argument(
        '--base-dir-gn',
        help='GN base directory (default: valid_combinations_GN/valid_no_overlap)'
    )
    parser.add_argument(
        '--base-dir-lf',
        help='LF base directory (default: valid_combinations_LF/valid_no_overlap)'
    )
    parser.add_argument(
        '--target-charge',
        type=float,
        default=2.0,
        help='Carga objetivo para fix_charge_drift (default: 2.0)'
    )
    parser.add_argument(
        '--tolerance',
        type=float,
        default=0.008,
        help='Max tolerance for fix_charge_drift (default: 0.008)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print actions only; do not modify files'
    )
    
    args = parser.parse_args()
    
    singularity_img = find_singularity_image()
    if singularity_img:
        print(f"✓ Singularity image: {singularity_img}")
    else:
        print("⚠ No Singularity image with RDKit found")
        print("  Trying host Python (may fail if RDKit is not installed)")
    
    # Establecer directorios base por defecto (buscar en ambos lugares posibles)
    if args.base_dir_gn is None:
        if os.path.exists("valid_combinations_GN/valid_no_overlap"):
            args.base_dir_gn = "valid_combinations_GN/valid_no_overlap"
        elif os.path.exists("valid_combinations_GN"):
            args.base_dir_gn = "valid_combinations_GN"
        else:
            args.base_dir_gn = "valid_combinations_GN/valid_no_overlap"
    if args.base_dir_lf is None:
        if os.path.exists("valid_combinations_LF/valid_no_overlap"):
            args.base_dir_lf = "valid_combinations_LF/valid_no_overlap"
        elif os.path.exists("valid_combinations_LF"):
            args.base_dir_lf = "valid_combinations_LF"
        else:
            args.base_dir_lf = "valid_combinations_LF/valid_no_overlap"
    
    print("=" * 70)
    print("Merging fragments for valid combinations")
    print("=" * 70)
    print(f"Carga objetivo: {args.target_charge}")
    print(f"Tolerancia: {args.tolerance}")
    print()
    
    total_processed = 0
    total_errors = 0
    all_errors = []
    
    if args.type in ['GN', 'all']:
        if os.path.exists(args.base_dir_gn):
            print(f"\n📁 Processing GN: {args.base_dir_gn}")
            processed, errors, errors_list = process_all_combinations(
                args.base_dir_gn, 'GN', args.target_charge, args.tolerance, args.dry_run, singularity_img
            )
            total_processed += processed
            total_errors += errors
            all_errors.extend(errors_list)
        else:
            print(f"\n⚠ No existe el directorio GN: {args.base_dir_gn}")
    
    if args.type in ['LF', 'all']:
        if os.path.exists(args.base_dir_lf):
            print(f"\n📁 Processing LF: {args.base_dir_lf}")
            processed, errors, errors_list = process_all_combinations(
                args.base_dir_lf, 'LF', args.target_charge, args.tolerance, args.dry_run, singularity_img
            )
            total_processed += processed
            total_errors += errors
            all_errors.extend(errors_list)
        else:
            print(f"\n⚠ LF directory does not exist: {args.base_dir_lf}")
    
    print("\n" + "=" * 70)
    print(f"Summary: {total_processed} combinations processed, {total_errors} errors")
    print("=" * 70)
    
    if all_errors:
        print(f"\nErrors ({len(all_errors)}):")
        for combo_name, error in all_errors[:10]:
            print(f"  - {combo_name}: {error}")
        if len(all_errors) > 10:
            print(f"  ... and {len(all_errors) - 10} more")
    
    if total_errors > 0:
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())


