#!/usr/bin/env python3
"""
Prepare the final combination layout for GN and LF.

Creates:
- valid_GN_final/combination_X/query_combination_X/fragmento_final_charge_drift.mol2
- valid_GN_final/combination_X/1a1m_A_GN.pdb

- valid_LF_final/combination_X/query_combination_X/fragmento_final_charge_drift.mol2
- valid_LF_final/combination_X/1a1m_A_LF.pdb
"""

import os
import re
import glob
import fnmatch
import shutil
import sys
import argparse


def find_pdb_file(pdb_name, search_dirs=None):
    """
    Find a PDB under the given dirs or VS_GN/VS_LF.
    pdb_name may be a glob pattern (e.g. 1ab9_a_*.pdb).

    Args:
        pdb_name: PDB filename or glob (e.g. 1a1m_A.pdb, 1ab9_a_*.pdb).
        search_dirs: Base directories to search (default: VS_GN*, VS_LF*).

    Returns:
        Path to the PDB, or None if not found.
    """
    if search_dirs is None:
        search_dirs = []
        for pattern in ["VS_GN*", "VS_LF*"]:
            search_dirs.extend(glob.glob(pattern))
    
    use_fnmatch = '*' in pdb_name or '?' in pdb_name
    
    for base_dir in search_dirs:
        pdb_path = os.path.join(base_dir, "results", "best_scores", pdb_name)
        if not use_fnmatch and os.path.exists(pdb_path):
            return pdb_path
        if use_fnmatch and os.path.isdir(os.path.dirname(pdb_path)):
            for f in (os.listdir(os.path.dirname(pdb_path)) or []):
                if fnmatch.fnmatch(f, pdb_name) and f.endswith('.pdb'):
                    return os.path.join(os.path.dirname(pdb_path), f)
        
        for root, _dirs, files in os.walk(base_dir):
            for f in files:
                if not f.endswith('.pdb'):
                    continue
                if use_fnmatch and fnmatch.fnmatch(f, pdb_name):
                    return os.path.join(root, f)
                if not use_fnmatch and f == pdb_name:
                    return os.path.join(root, f)
    
    return None


def resolve_pdb_file(pdb_file, pdb_search_dirs=None):
    """
    Resolve PDB path: literal path, glob, or name search.
    If --pdb-file is a glob (e.g. 1ab9_a_*.pdb), searches cwd, --pdb-search-dir,
    and common locations (Propedia_Final_mol2, Propedia_mol2_4Frag, STELLAR_4Frag, {case}).

    Returns:
        (absolute_path, basename_for_output) or (None, None) if not found.
    """
    pdb_search_dirs = list(pdb_search_dirs or [])
    
    if os.path.exists(pdb_file):
        return (os.path.abspath(pdb_file), os.path.basename(pdb_file))
    
    base = os.path.basename(pdb_file)
    is_glob = '*' in pdb_file or '?' in base
    
    if is_glob:
        # 1) Current directory
        matches = glob.glob(pdb_file)
        if matches:
            m = [x for x in matches if os.path.isfile(x)]
            if m:
                return (os.path.abspath(m[0]), os.path.basename(m[0]))
        
        # 2) Infer case (e.g. 1ab9_a_*.pdb -> 1AB9_A) and common dirs
        m = re.match(r'^([0-9][a-zA-Z0-9]+_[a-zA-Z])', base, re.IGNORECASE)
        extra = ['.']
        if m:
            case = m.group(1).upper()
            extra.extend([
                case,
                os.path.join("Propedia_Final_mol2", case, "Receptor"),
                os.path.join("Propedia_mol2_4Frag", case),
                os.path.join("STELLAR_4Frag", case),
            ])
        # User --pdb-search-dir first, then inferred dirs
        search_dirs = list(dict.fromkeys(pdb_search_dirs + extra))
        
        # 3) Looser patterns: 1ab9_a_*.pdb -> also 1ab9_a*.pdb and 1ab9_a.pdb
        patterns = [base]
        if re.search(r'_[*?]', base):
            patterns.append(re.sub(r'_[*?][^*?]*\.pdb$', '*.pdb', base))
            patterns.append(re.sub(r'_[*?][^*?]*\.pdb$', '.pdb', base))
        
        for d in search_dirs:
            if not os.path.isdir(d):
                continue
            for pat in patterns:
                for p in glob.glob(os.path.join(d, pat)):
                    if os.path.isfile(p):
                        return (os.path.abspath(p), os.path.basename(p))
                # Recursive search only under user --pdb-search-dir (avoids huge trees)
                if d in (pdb_search_dirs or []):
                    for p in glob.glob(os.path.join(d, "**", pat)):
                        if os.path.isfile(p):
                            return (os.path.abspath(p), os.path.basename(p))
    
    # 4) Name search in VS_GN*, VS_LF*, and pdb_search_dirs
    all_dirs = pdb_search_dirs + [x for x in glob.glob("VS_GN*") + glob.glob("VS_LF*") if os.path.isdir(x)]
    found = find_pdb_file(base, all_dirs if all_dirs else None)
    if found:
        return (os.path.abspath(found), os.path.basename(found))
    
    return (None, None)


def prepare_combination(combo_dir, output_base_dir, pdb_file, pdb_name, prefix_type, dry_run=False):
    """
    Prepare one combination: copy files into the final layout.

    Args:
        combo_dir: Source combination directory.
        output_base_dir: Output root (valid_GN_final or valid_LF_final).
        pdb_file: Path to PDB to copy.
        pdb_name: PDB basename (e.g. 1a1m_A.pdb).
        prefix_type: Prefix type ('GN' or 'LF').
        dry_run: If True, only print what would be done.

    Returns:
        (success, error_message)
    """
    combo_name = os.path.basename(combo_dir)
    
    # Combination index (e.g. combination_102 -> 102)
    combo_number = combo_name.replace("combination_", "")
    
    # Directory layout
    output_combo_dir = os.path.join(output_base_dir, combo_name)
    query_dir = os.path.join(output_combo_dir, f"query_combination_{combo_number}")
    
    # PDB name with _GN or _LF suffix
    pdb_base = os.path.splitext(pdb_name)[0]  # stem
    pdb_ext = os.path.splitext(pdb_name)[1]  # .pdb
    pdb_final_name = f"{pdb_base}_{prefix_type}{pdb_ext}"
    
    # Files to copy
    mol2_source = os.path.join(combo_dir, "fragmento_final_charge_drift.mol2")
    mol2_dest = os.path.join(query_dir, "fragmento_final_charge_drift.mol2")
    pdb_dest = os.path.join(output_combo_dir, pdb_final_name)
    
    if not os.path.exists(mol2_source):
        return False, f"Not found: {mol2_source}"
    
    if not os.path.exists(pdb_file):
        return False, f"PDB file not found: {pdb_file}"
    
    if dry_run:
        print(f"  [DRY RUN] Create: {query_dir}/")
        print(f"  [DRY RUN] Copy: {mol2_source} -> {mol2_dest}")
        print(f"  [DRY RUN] Copy: {pdb_file} -> {pdb_dest} (as {pdb_final_name})")
        return True, None
    
    try:
        os.makedirs(query_dir, exist_ok=True)
        shutil.copy2(mol2_source, mol2_dest)
        shutil.copy2(pdb_file, pdb_dest)
        
        return True, None
    
    except Exception as e:
        return False, f"Error copying files: {e}"


def process_all_combinations(base_dir, output_base_dir, prefix_type, pdb_file, pdb_name, dry_run=False):
    """
    Process all combinations under base_dir.

    Args:
        base_dir: Root containing combination_* folders.
        output_base_dir: Output root.
        prefix_type: Prefix type ('GN' or 'LF').
        pdb_file: Path to PDB to copy.
        pdb_name: PDB basename.
        dry_run: If True, only print actions.

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
    print()
    
    total_processed = 0
    total_errors = 0
    errors_list = []
    
    for i, combo_dir in enumerate(combo_dirs, 1):
        combo_name = os.path.basename(combo_dir)
        print(f"[{i}/{len(combo_dirs)}] {combo_name}")
        
        success, error = prepare_combination(combo_dir, output_base_dir, pdb_file, pdb_name, prefix_type, dry_run)
        
        if success:
            total_processed += 1
            if not dry_run:
                print(f"  ✓ Prepared successfully")
        else:
            total_errors += 1
            errors_list.append((combo_name, error))
            print(f"  ✗ {error}")
        
        if (i % 10 == 0) and not dry_run:
            print(f"  Progress: {total_processed} ok, {total_errors} errors...")
    
    return total_processed, total_errors, errors_list


def main():
    parser = argparse.ArgumentParser(
        description="Prepare final combination layout for GN and LF"
    )
    parser.add_argument(
        'type',
        choices=['GN', 'LF', 'all'],
        help='Combination type: GN, LF, or all (both)'
    )
    parser.add_argument(
        '--pdb-file',
        required=True,
        help='PDB path, glob (e.g. 1ab9_a_*.pdb), or name; if missing, search --pdb-search-dir and Propedia/STELLAR/VS_GN/VS_LF'
    )
    parser.add_argument(
        '--pdb-name',
        help='PDB basename in output layout (default: from resolved file)'
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
        '--output-dir-gn',
        help='GN output directory (default: valid_GN_final)'
    )
    parser.add_argument(
        '--output-dir-lf',
        help='LF output directory (default: valid_LF_final)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print actions only; do not modify files'
    )
    parser.add_argument(
        '--pdb-search-dir',
        action='append',
        default=[],
        dest='pdb_search_dirs',
        metavar='DIR',
        help='Extra directories to search for PDB when --pdb-file is a glob or missing (repeatable)'
    )
    
    args = parser.parse_args()
    
    # Default base directories
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
    if args.output_dir_gn is None:
        args.output_dir_gn = "valid_GN_final"
    if args.output_dir_lf is None:
        args.output_dir_lf = "valid_LF_final"
    
    # Resolve PDB: literal path, glob, or search
    pdb_file, pdb_name_resolved = resolve_pdb_file(args.pdb_file, getattr(args, 'pdb_search_dirs', []))
    if pdb_file is None:
        print(f"✗ Error: PDB not found: '{args.pdb_file}'")
        print("  Searched: ., Propedia_Final_mol2/<case>/Receptor, Propedia_mol2_4Frag, STELLAR_4Frag, VS_GN*, VS_LF*")
        print("  Use --pdb-search-dir DIR to add search paths.")
        return 1
    args.pdb_name = pdb_name_resolved if args.pdb_name is None else args.pdb_name
    print(f"✓ PDB file: {pdb_file}")
    
    print("=" * 70)
    print("Final combination layout preparation")
    print("=" * 70)
    print(f"PDB file: {pdb_file}")
    print(f"Output PDB name: {args.pdb_name}")
    print()
    
    total_processed = 0
    total_errors = 0
    all_errors = []
    
    # Process by type
    if args.type in ['GN', 'all']:
        if os.path.exists(args.base_dir_gn):
            # Remove existing output to regenerate clean
            if os.path.exists(args.output_dir_gn) and not args.dry_run:
                import shutil
                print(f"\n🗑️  Removing existing directory: {args.output_dir_gn}")
                shutil.rmtree(args.output_dir_gn)
            
            print(f"\n📁 Processing GN:")
            print(f"   Source: {args.base_dir_gn}")
            print(f"   Output: {args.output_dir_gn}")
            processed, errors, errors_list = process_all_combinations(
                args.base_dir_gn, args.output_dir_gn, 'GN', pdb_file, args.pdb_name, args.dry_run
            )
            total_processed += processed
            total_errors += errors
            all_errors.extend(errors_list)
        else:
            print(f"\n⚠ GN directory does not exist: {args.base_dir_gn}")
    
    if args.type in ['LF', 'all']:
        if os.path.exists(args.base_dir_lf):
            # Eliminar directorio de salida si existe (para regenerar limpio)
            if os.path.exists(args.output_dir_lf) and not args.dry_run:
                import shutil
                print(f"\n🗑️  Removing existing directory: {args.output_dir_lf}")
                shutil.rmtree(args.output_dir_lf)
            
            print(f"\n📁 Processing LF:")
            print(f"   Source: {args.base_dir_lf}")
            print(f"   Output: {args.output_dir_lf}")
            processed, errors, errors_list = process_all_combinations(
                args.base_dir_lf, args.output_dir_lf, 'LF', pdb_file, args.pdb_name, args.dry_run
            )
            total_processed += processed
            total_errors += errors
            all_errors.extend(errors_list)
        else:
            print(f"\n⚠ No existe el directorio LF: {args.base_dir_lf}")
    
    print("\n" + "=" * 70)
    print(f"Summary: {total_processed} combinations prepared, {total_errors} errors")
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


