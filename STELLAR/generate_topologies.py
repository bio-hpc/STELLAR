#!/usr/bin/env python3
"""
Generate GROMACS topologies for all valid fragment combinations.

Runs:
singularity exec gr.simg python GROMACS/external_sw/gromacs/topology/generate_topology.py \
    -t {combination_dir}/1a1m_A.pdb \
    -q {combination_dir}/query_combination_{num}/
"""

import os
import glob
import subprocess
import sys
import argparse
from pathlib import Path


def query_combination_has_valid_topology(combo_dir):
    """
    True if query_combination_{N}/ contains *complex.top and *query.itp (existing GROMACS topology).

    Args:
        combo_dir: Combination folder (e.g. valid_GN_final/combination_102)

    Returns:
        bool
    """
    combo_name = os.path.basename(combo_dir)
    combo_number = combo_name.replace("combination_", "")
    query_dir = os.path.join(combo_dir, f"query_combination_{combo_number}")
    if not os.path.isdir(query_dir):
        return False
    try:
        files = os.listdir(query_dir)
    except OSError:
        return False
    has_complex_top = any(f.endswith("complex.top") for f in files)
    has_query_itp = any(f.endswith("query.itp") for f in files)
    return has_complex_top and has_query_itp


def all_combinations_have_valid_topology(base_dir):
    """
    Whether every combination_* under base_dir already has valid topology
    (complex.top + query.itp in each query_combination_*).

    Returns:
        tuple: (all_ok: bool, total: int, with_topology: int)
    """
    combo_dirs = sorted(glob.glob(os.path.join(base_dir, "combination_*")))
    if not combo_dirs:
        return False, 0, 0
    with_topology = sum(1 for d in combo_dirs if query_combination_has_valid_topology(d))
    return with_topology == len(combo_dirs), len(combo_dirs), with_topology


def generate_topology(combo_dir, prefix_type, singularity_image="gr.simg", dry_run=False):
    """
    Generate topology for one combination.

    Args:
        combo_dir: Combination directory (e.g. valid_GN_final/combination_102)
        prefix_type: 'GN' or 'LF'
        singularity_image: Singularity image (default: gr.simg)
        dry_run: If True, only print the command

    Returns:
        tuple: (success, error_message or None)
    """
    combo_name = os.path.basename(combo_dir)

    combo_number = combo_name.replace("combination_", "")

    pdb_files = glob.glob(os.path.join(combo_dir, "*.pdb"))

    pdb_file = None

    for pdb_path in pdb_files:
        if pdb_path.endswith(f"_{prefix_type}.pdb") and "2w10" in os.path.basename(pdb_path).lower():
            pdb_file = pdb_path
            break

    if pdb_file is None:
        for pdb_path in pdb_files:
            if pdb_path.endswith(f"_{prefix_type}.pdb"):
                pdb_file = pdb_path
                break

    if pdb_file is None:
        if pdb_files:
            pdb_file = pdb_files[0]
        else:
            pdb_file = os.path.join(combo_dir, f"*_{prefix_type}.pdb")

    query_dir = os.path.join(combo_dir, f"query_combination_{combo_number}")

    if not os.path.exists(pdb_file):
        return False, f"PDB file not found: {pdb_file}"

    if not os.path.exists(query_dir):
        return False, f"Query directory not found: {query_dir}"

    if not os.path.exists(singularity_image) and not dry_run:
        return False, f"Singularity image not found: {singularity_image}"

    script_path = "GROMACS/external_sw/gromacs/topology/generate_topology.py"
    current_dir = os.getcwd()

    pdb_rel = os.path.relpath(pdb_file, current_dir)
    query_rel = os.path.relpath(query_dir, current_dir)
    # Trailing slash on query path helps ligand detection
    if not query_rel.endswith('/'):
        query_rel = query_rel + '/'

    cmd = [
        "singularity", "exec",
        singularity_image,
        "python", script_path,
        "-t", pdb_rel,
        "-q", query_rel
    ]

    if dry_run:
        print(f"  [DRY RUN] {' '.join(cmd)}")
        return True, None

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            cwd=current_dir
        )

        if result.returncode == 0:
            return True, None
        else:
            error_msg = result.stderr.strip() if result.stderr else result.stdout.strip()
            if error_msg:
                error_lines = error_msg.split('\n')
                error_preview = '\n'.join(error_lines[-5:]) if len(error_lines) > 5 else error_msg
                return False, f"Execution error:\n{error_preview}"
            else:
                return False, f"Unknown error (exit code: {result.returncode})"

    except Exception as e:
        return False, f"Exception: {e}"


def count_expected_topology_files(base_dir):
    """
    Count combination_* dirs that contain at least one .top file.

    Returns:
        tuple: (count_with_top, total_combination_dirs)
    """
    combo_dirs = sorted(glob.glob(os.path.join(base_dir, "combination_*")))
    count_with_top = 0
    for combo_dir in combo_dirs:
        top_files = glob.glob(os.path.join(combo_dir, "*.top"))
        if top_files:
            count_with_top += 1
    return count_with_top, len(combo_dirs)


def process_all_combinations(base_dir, prefix_type, singularity_image="gr.simg", dry_run=False, max_combinations=None):
    """
    Process all combinations under base_dir.

    Args:
        base_dir: Root with combination_* folders
        prefix_type: 'GN' or 'LF'
        singularity_image: Singularity image path
        dry_run: If True, only show what would run
        max_combinations: If set, only first N combinations (e.g. 1 for smoke test)

    Returns:
        tuple: (total_processed, total_errors, errors_list)
    """
    combo_dirs = sorted(glob.glob(os.path.join(base_dir, "combination_*")))

    if not combo_dirs:
        print(f"⚠ No combination_* folders in {base_dir}")
        return 0, 0, []

    if max_combinations is not None and max_combinations > 0:
        combo_dirs = combo_dirs[:max_combinations]
        print(f"Found {len(combo_dirs)} combinations (limit: first {max_combinations})")
    else:
        print(f"Found {len(combo_dirs)} combinations")
    if dry_run:
        print("⚠ DRY RUN — no commands will be executed")
    print()

    total_processed = 0
    total_errors = 0
    errors_list = []

    for i, combo_dir in enumerate(combo_dirs, 1):
        combo_name = os.path.basename(combo_dir)
        print(f"[{i}/{len(combo_dirs)}] {combo_name}")

        if query_combination_has_valid_topology(combo_dir):
            if not dry_run:
                print(f"  ✓ Valid topology already present (complex.top + query.itp); skipping.")
            total_processed += 1
            continue

        success, error = generate_topology(combo_dir, prefix_type, singularity_image, dry_run)

        if success:
            total_processed += 1
            if not dry_run:
                print(f"  ✓ Topology generated successfully")
        else:
            total_errors += 1
            errors_list.append((combo_name, error))
            print(f"  ✗ {error}")

        if (i % 10 == 0) and not dry_run:
            print(f"  Progress: {total_processed} ok, {total_errors} failed...")

    return total_processed, total_errors, errors_list


def main():
    parser = argparse.ArgumentParser(
        description="Generate GROMACS topologies for all valid combinations"
    )
    parser.add_argument(
        'type',
        choices=['GN', 'LF', 'all'],
        help='Combination type: GN, LF, or all (both)'
    )
    parser.add_argument(
        '--base-dir-gn',
        help='Base directory for GN (default: valid_GN_final)'
    )
    parser.add_argument(
        '--base-dir-lf',
        help='Base directory for LF (default: valid_LF_final)'
    )
    parser.add_argument(
        '--singularity-image',
        default='gr.simg',
        help='Singularity image (default: gr.simg)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print commands only; do not run'
    )
    parser.add_argument(
        '--check-only',
        action='store_true',
        help='Exit 0 if every combination has complex.top and query.itp under query_combination_*; else 1 (for workflow jobs)'
    )
    parser.add_argument(
        '--max-combinations',
        type=int,
        default=None,
        metavar='N',
        help='Process only the first N combinations (e.g. --max-combinations 1 for a smoke test)'
    )

    args = parser.parse_args()

    if args.base_dir_gn is None:
        args.base_dir_gn = "valid_GN_final"
    if args.base_dir_lf is None:
        args.base_dir_lf = "valid_LF_final"

    if args.check_only:
        base_to_check = args.base_dir_gn if args.type in ('GN', 'all') else args.base_dir_lf
        if args.type == 'all':
            gn_ok, gn_tot, gn_has = all_combinations_have_valid_topology(args.base_dir_gn) if os.path.exists(args.base_dir_gn) else (True, 0, 0)
            lf_ok, lf_tot, lf_has = all_combinations_have_valid_topology(args.base_dir_lf) if os.path.exists(args.base_dir_lf) else (True, 0, 0)
            all_ok = gn_ok and lf_ok
        else:
            all_ok, total, with_top = all_combinations_have_valid_topology(base_to_check) if os.path.exists(base_to_check) else (False, 0, 0)
        sys.exit(0 if all_ok else 1)

    print("=" * 70)
    print("GROMACS topology generation")
    print("=" * 70)
    print(f"Singularity image: {args.singularity_image}")
    print()

    total_processed = 0
    total_errors = 0
    all_errors = []

    max_comb = getattr(args, 'max_combinations', None)
    if args.type in ['GN', 'all']:
        if os.path.exists(args.base_dir_gn):
            print(f"\n📁 GN: {args.base_dir_gn}")
            processed, errors, errors_list = process_all_combinations(
                args.base_dir_gn, 'GN', args.singularity_image, args.dry_run, max_combinations=max_comb
            )
            total_processed += processed
            total_errors += errors
            all_errors.extend(errors_list)
        else:
            print(f"\n⚠ GN directory missing: {args.base_dir_gn}")

    if args.type in ['LF', 'all']:
        if os.path.exists(args.base_dir_lf):
            print(f"\n📁 LF: {args.base_dir_lf}")
            processed, errors, errors_list = process_all_combinations(
                args.base_dir_lf, 'LF', args.singularity_image, args.dry_run, max_combinations=max_comb
            )
            total_processed += processed
            total_errors += errors
            all_errors.extend(errors_list)
        else:
            print(f"\n⚠ LF directory missing: {args.base_dir_lf}")

    # Abort only if zero valid topologies (e.g. all acpype failures). If some succeeded, continue.
    # Skip this check when --max-combinations was used.
    if not args.dry_run and max_comb is None:
        if args.type in ['GN', 'all'] and os.path.exists(args.base_dir_gn):
            all_ok, total_combo, with_top = all_combinations_have_valid_topology(args.base_dir_gn)
            if total_combo > 0 and not all_ok:
                if with_top == 0:
                    print("\n" + "=" * 70)
                    print("ERROR: No valid GN topology was generated.")
                    print(f"With complex.top + query.itp: {with_top}/{total_combo}.")
                    print("Aborting job.")
                    print("=" * 70)
                    return 1
                print("\n" + "=" * 70)
                print("WARNING: Not all GN combinations have valid topology.")
                print(f"With complex.top + query.itp: {with_top}/{total_combo}.")
                print("Continuing with valid ones (others will not run MD).")
                print("=" * 70)
        if args.type in ['LF', 'all'] and os.path.exists(args.base_dir_lf):
            all_ok, total_combo, with_top = all_combinations_have_valid_topology(args.base_dir_lf)
            if total_combo > 0 and not all_ok:
                if with_top == 0:
                    print("\n" + "=" * 70)
                    print("ERROR: No valid LF topology was generated.")
                    print(f"With complex.top + query.itp: {with_top}/{total_combo}.")
                    print("Aborting job.")
                    print("=" * 70)
                    return 1
                print("\n" + "=" * 70)
                print("WARNING: Not all LF combinations have valid topology.")
                print(f"With complex.top + query.itp: {with_top}/{total_combo}.")
                print("Continuing with valid ones (others will not run MD).")
                print("=" * 70)

    print("\n" + "=" * 70)
    if args.dry_run:
        print(f"Summary (DRY RUN): {total_processed} combinations ready to process")
    else:
        print(f"Summary: {total_processed} topologies ok, {total_errors} failures")
    print("=" * 70)

    if all_errors:
        print(f"\nErrors ({len(all_errors)}):")
        for combo_name, error in all_errors[:10]:
            print(f"  - {combo_name}: {error}")
            if len(error.split('\n')) > 1:
                print()
        if len(all_errors) > 10:
            print(f"  ... and {len(all_errors) - 10} more")

    if total_errors > 0 and total_processed == 0:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
