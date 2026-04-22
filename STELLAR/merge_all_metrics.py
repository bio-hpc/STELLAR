#!/usr/bin/env python3
"""
Merge all metrics into one table:
- Per-fragment RMSD (rmsd_results.csv)
- Docking poses and scores (valid_fragment_combinations_*_no_overlap.csv)
- MD peptide RMSD (resultados_rmsd_md.csv)
"""

import csv
import argparse
import sys
from pathlib import Path


def normalize_combo_id(key):
    """
    Normalize combination_id to canonical string integer form
    so keys like "1", "1.0", " 1 ", "001" match across files.
    """
    if key is None or (isinstance(key, str) and not key.strip()):
        return ''
    s = str(key).strip()
    try:
        return str(int(float(s)))
    except (ValueError, TypeError):
        return s


def load_csv_dict(csv_file, key_column, encoding='utf-8', normalize_key=True):
    """
    Load a CSV into a dict keyed by key_column.

    Args:
        csv_file: Path to CSV
        key_column: Column name to use as key
        encoding: File encoding
        normalize_key: If True, normalize keys (e.g. "1.0" -> "1") for joins

    Returns:
        dict: Rows keyed by key_column
    """
    data = {}
    if not Path(csv_file).exists():
        print(f"⚠ Warning: file not found: {csv_file}")
        return data
    
    with open(csv_file, 'r', encoding=encoding) as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = row[key_column]
            k = normalize_combo_id(key) if normalize_key else key
            if k:
                data[k] = row
    
    return data


def merge_metrics_single_type(
    combinations_file,
    rmsd_fragments_file,
    md_rmsd_file,
    output_file,
    prefix_type="GN",
    score_only_file=None,
    mmpbsa_file=None,
    fragment_energies_file=None,
):
    """
    Merge all metrics into one CSV for a single type (GN or LF).

    Args:
        combinations_file: CSV with poses, scores, distances
        rmsd_fragments_file: CSV with per-fragment RMSD
        md_rmsd_file: CSV with MD RMSD
        output_file: Output CSV path
        prefix_type: 'GN' or 'LF'
        score_only_file: Optional score_only CSV
        mmpbsa_file: Optional MM/PBSA CSV
        fragment_energies_file: Optional per-fragment energy CSV
    """
    print("=" * 70)
    print(f"Merging fragment-combination metrics - {prefix_type}")
    print("=" * 70)
    print(f"Combinations file: {combinations_file}")
    print(f"Fragment RMSD file: {rmsd_fragments_file}")
    print(f"MD RMSD file: {md_rmsd_file}")
    if score_only_file:
        print(f"Score-only file: {score_only_file}")
    if mmpbsa_file:
        print(f"MM/PBSA file: {mmpbsa_file}")
    if fragment_energies_file:
        print(f"Fragment energies file: {fragment_energies_file}")
    print(f"Output file: {output_file}")
    print()
    
    # Load inputs
    print("Loading data...")
    combinations_data = load_csv_dict(combinations_file, 'combination_id')
    rmsd_fragments_data = load_csv_dict(rmsd_fragments_file, 'combination_id')
    
    # Load score_only when available
    score_only_data = {}
    if score_only_file and Path(score_only_file).exists():
        score_only_data = load_csv_dict(score_only_file, 'combination_id')
        print(f"  Score_only rows loaded: {len(score_only_data)}")

    # Load MM/PBSA when available
    mmpbsa_data_all = {}
    mmpbsa_filtered = {}
    if mmpbsa_file and Path(mmpbsa_file).exists():
        mmpbsa_data_all = load_csv_dict(mmpbsa_file, 'combination')
        print(f"  MM/PBSA rows loaded: {len(mmpbsa_data_all)}")

    # Load per-fragment energy contributions when available
    fragment_energies_data = {}
    fragment_energies_columns = []
    if fragment_energies_file and Path(fragment_energies_file).exists():
        fragment_energies_data = load_csv_dict(fragment_energies_file, 'combination_id')
        if fragment_energies_data:
            first_row = next(iter(fragment_energies_data.values()))
            fragment_energies_columns = [k for k in first_row.keys() if k != 'combination_id']
        print(f"  Fragment energies loaded: {len(fragment_energies_data)} ({len(fragment_energies_columns)} columns)")
    # Load MD with multiple rows per combination (GN and LF)
    md_rmsd_data_all = {}
    if Path(md_rmsd_file).exists():
        with open(md_rmsd_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                combo_id = row.get('combination', '')
                if combo_id:
                    # Unique key: combination + type (normalized ID)
                    combo_norm = normalize_combo_id(combo_id)
                    vs_folder = row.get('vs_folder', '')
                    is_lf = 'LF' in vs_folder or '_LF_' in vs_folder
                    is_gn = 'GN' in vs_folder or '_GN_' in vs_folder
                    if is_lf:
                        key = f"{combo_norm}_LF"
                    elif is_gn:
                        key = f"{combo_norm}_GN"
                    else:
                        key = combo_norm
                    md_rmsd_data_all[key] = row
    
    print(f"  Combinations loaded: {len(combinations_data)}")
    print(f"  Fragment RMSD rows loaded: {len(rmsd_fragments_data)}")
    print(f"  MD RMSD rows loaded (total): {len(md_rmsd_data_all)}")
    if mmpbsa_data_all:
        print(f"  MM/PBSA (total): {len(mmpbsa_data_all)}")
    if fragment_energies_data:
        print(f"  Fragment energies (total): {len(fragment_energies_data)}")
    print()
    
    # Filter MD by type (GN or LF) using vs_folder
    filtered_md_data = {}
    for key, row in md_rmsd_data_all.items():
        vs_folder = row.get('vs_folder', '')
        is_lf = 'LF' in vs_folder or '_LF_' in vs_folder
        is_gn = 'GN' in vs_folder or '_GN_' in vs_folder
        # Filter explicitly by type
        if (prefix_type == 'LF' and is_lf) or (prefix_type == 'GN' and is_gn):
            # Use only the combination number as key (normalized)
            combo_id = row.get('combination', '')
            if combo_id:
                filtered_md_data[normalize_combo_id(combo_id)] = row
    
    print(f"  MD RMSD filtered ({prefix_type}): {len(filtered_md_data)}")
    print()
    
    # Filtrar MM/PBSA por tipo
    if mmpbsa_data_all:
        for combo, row in mmpbsa_data_all.items():
            vs_folder = row.get('vs_folder', '')
            row_type = row.get('type', '')
            is_lf = 'LF' in vs_folder or '_LF_' in vs_folder or row_type == 'LF'
            is_gn = 'GN' in vs_folder or '_GN_' in vs_folder or row_type == 'GN'
            if (prefix_type == 'LF' and is_lf) or (prefix_type == 'GN' and is_gn):
                mmpbsa_filtered[combo] = row

    # All unique combinations for this type
    all_combinations = set(combinations_data.keys())
    all_combinations.update(rmsd_fragments_data.keys())
    all_combinations.update(filtered_md_data.keys())
    if score_only_data:
        all_combinations.update(score_only_data.keys())
    if mmpbsa_filtered:
        all_combinations.update(mmpbsa_filtered.keys())
    if fragment_energies_data:
        all_combinations.update(fragment_energies_data.keys())
    all_combinations = sorted([int(c) for c in all_combinations if c.isdigit()])
    
    print(f"Total unique combinations: {len(all_combinations)}")
    print()
    
    # Infer max fragment count from data
    max_fragments = 0
    # Scan combinations_data
    for combo_row in combinations_data.values():
        for key in combo_row.keys():
            if key.startswith('frag') and '_pose' in key:
                try:
                    frag_num = int(key.replace('frag', '').replace('_pose', ''))
                    max_fragments = max(max_fragments, frag_num)
                except ValueError:
                    pass
    # Scan rmsd_fragments_data
    for rmsd_row in rmsd_fragments_data.values():
        for key in rmsd_row.keys():
            if key.startswith('rmsd_frag'):
                try:
                    frag_num = int(key.replace('rmsd_frag', ''))
                    max_fragments = max(max_fragments, frag_num)
                except ValueError:
                    pass
    
    if max_fragments == 0:
        # No fragments found in data: use a reasonable default
        max_fragments = 4
        print(f"⚠ No fragments detected in data, using default {max_fragments}")
    else:
        print(f"✓ Detected {max_fragments} fragments")
    print()
    
    # Build output columns dynamically
    fieldnames = ['combination_id']
    for i in range(1, max_fragments + 1):
        fieldnames.extend([f'frag{i}_pose', f'frag{i}_score', f'rmsd_frag{i}'])
    fieldnames.append('total_score')
    for i in range(1, max_fragments):
        fieldnames.append(f'distance_frag{i}_to_frag{i+1}')
    fieldnames.extend(['rmsd_md', 'vs_folder', 'successfully_generated'])
    if score_only_data:
        fieldnames.append('score_only')
    if mmpbsa_filtered:
        fieldnames.extend(['mmpbsa_vdw', 'mmpbsa_elec', 'mmpbsa_total', 'mmpbsa_frames'])
    if fragment_energies_columns:
        fieldnames.extend(fragment_energies_columns)
    
    # Write merged CSV
    print("Writing merged CSV...")
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        merged_count = 0
        for combo_id in all_combinations:
            combo_str = str(combo_id)
            row = {'combination_id': combo_str}
            
            # Combination data (poses, scores, distances)
            if combo_str in combinations_data:
                combo_row = combinations_data[combo_str]
                for i in range(1, max_fragments + 1):
                    row[f'frag{i}_pose'] = combo_row.get(f'frag{i}_pose', '')
                    row[f'frag{i}_score'] = combo_row.get(f'frag{i}_score', '')
                row['total_score'] = combo_row.get('total_score', '')
                for i in range(1, max_fragments):
                    row[f'distance_frag{i}_to_frag{i+1}'] = combo_row.get(f'distance_frag{i}_to_frag{i+1}', '')
            else:
                # Fill with empty values
                for i in range(1, max_fragments + 1):
                    row[f'frag{i}_pose'] = ''
                    row[f'frag{i}_score'] = ''
                row['total_score'] = ''
                for i in range(1, max_fragments):
                    row[f'distance_frag{i}_to_frag{i+1}'] = ''
            
            # Fragment RMSD data
            if combo_str in rmsd_fragments_data:
                rmsd_row = rmsd_fragments_data[combo_str]
                for i in range(1, max_fragments + 1):
                    row[f'rmsd_frag{i}'] = rmsd_row.get(f'rmsd_frag{i}', '')
            else:
                for i in range(1, max_fragments + 1):
                    row[f'rmsd_frag{i}'] = ''
            
            # Datos de RMSD de MD (ya filtrados por tipo)
            if combo_str in filtered_md_data:
                md_row = filtered_md_data[combo_str]
                row['rmsd_md'] = md_row.get('rmsd', '')
                row['vs_folder'] = md_row.get('vs_folder', '')
            else:
                row['rmsd_md'] = 'NA'
                row['vs_folder'] = 'No generated'
            # successfully_generated: 1 if rmsd_md is a valid number, 0 if NA/empty/non-numeric
            rmsd_val = (row.get('rmsd_md') or '').strip().upper()
            if rmsd_val in ('', 'NA', 'N/A', 'NAN'):
                row['successfully_generated'] = '0'
            else:
                try:
                    float(rmsd_val)
                    row['successfully_generated'] = '1'
                except ValueError:
                    row['successfully_generated'] = '0'
            
            # score_only data
            if score_only_data:
                if combo_str in score_only_data:
                    score_row = score_only_data[combo_str]
                    row['score_only'] = score_row.get('score_only', '')
                else:
                    row['score_only'] = ''

            # MM/PBSA data
            if mmpbsa_filtered:
                if combo_str in mmpbsa_filtered:
                    mmpbsa_row = mmpbsa_filtered[combo_str]
                    row['mmpbsa_vdw'] = mmpbsa_row.get('mmpbsa_vdw', '')
                    row['mmpbsa_elec'] = mmpbsa_row.get('mmpbsa_elec', '')
                    row['mmpbsa_total'] = mmpbsa_row.get('mmpbsa_total', '')
                    row['mmpbsa_frames'] = mmpbsa_row.get('frames', '')
                else:
                    row['mmpbsa_vdw'] = ''
                    row['mmpbsa_elec'] = ''
                    row['mmpbsa_total'] = ''
                    row['mmpbsa_frames'] = ''

            # Per-fragment energy contribution data
            if fragment_energies_columns:
                if combo_str in fragment_energies_data:
                    fe_row = fragment_energies_data[combo_str]
                    for col in fragment_energies_columns:
                        row[col] = fe_row.get(col, '')
                else:
                    for col in fragment_energies_columns:
                        row[col] = ''
            
            writer.writerow(row)
            merged_count += 1
    
    print(f"✓ Merged CSV written: {output_file}")
    print(f"  Total rows: {merged_count}")
    print()
    
    # Statistics (keys normalized; combo_str = str(combo_id) is "1", "2", ...)
    with_combinations = sum(1 for c in all_combinations if str(c) in combinations_data)
    with_rmsd_fragments = sum(1 for c in all_combinations if str(c) in rmsd_fragments_data)
    with_md_rmsd = sum(1 for c in all_combinations if str(c) in filtered_md_data)
    with_score_only = sum(1 for c in all_combinations if str(c) in score_only_data) if score_only_data else 0
    with_mmpbsa = sum(1 for c in all_combinations if str(c) in mmpbsa_filtered) if mmpbsa_filtered else 0
    with_fragment_energies = sum(1 for c in all_combinations if str(c) in fragment_energies_data) if fragment_energies_data else 0

    def has_all_metrics(combo_str):
        if combo_str not in combinations_data or combo_str not in rmsd_fragments_data or combo_str not in filtered_md_data:
            return False
        if score_only_data and combo_str not in score_only_data:
            return False
        if mmpbsa_filtered and combo_str not in mmpbsa_filtered:
            return False
        if fragment_energies_data and combo_str not in fragment_energies_data:
            return False
        return True

    n_all_metrics = sum(1 for c in all_combinations if has_all_metrics(str(c)))
    n_comb_rmsd_md = sum(1 for c in all_combinations if (
        str(c) in combinations_data and str(c) in rmsd_fragments_data and str(c) in filtered_md_data
    ))

    print("Statistics:")
    print(f"  With combination data: {with_combinations}")
    print(f"  With fragment RMSD: {with_rmsd_fragments}")
    print(f"  With MD RMSD: {with_md_rmsd}")
    print(f"  With combinations + frag RMSD + MD RMSD: {n_comb_rmsd_md}")
    if score_only_data:
        print(f"  With score_only: {with_score_only}")
    if mmpbsa_filtered:
        print(f"  With MM/PBSA: {with_mmpbsa}")
    if fragment_energies_data:
        print(f"  With fragment energies: {with_fragment_energies}")
    print(f"  With all metrics (all loaded sources): {n_all_metrics}")

    # Warn if the three base sources should align but mostly do not
    min_core = min(with_combinations, with_rmsd_fragments, with_md_rmsd)
    if min_core > 0 and n_comb_rmsd_md < min_core * 0.5:
        print()
        print("⚠ WARNING: Very few rows have all three base sources (combinations + frag RMSD + MD RMSD).")
        print("  This usually means the combinations CSV and the RMSD CSVs do not come from the same")
        print("  valid_no_overlap. Fix: regenerate the combinations CSV (step 14) using the same")
        print("  directory as for fragment RMSD (step 12) and MD RMSD (step 13):")
        print("  python3 STELLAR/filter_valid_combinations_csv.py GN --combinations-dir <valid_no_overlap> \\")
        print("    --input-csv valid_fragment_combinations.csv --output-csv valid_fragment_combinations_GN_no_overlap_8YP8_B.csv")


def main():
    parser = argparse.ArgumentParser(
        description="Merge fragment-combination metrics into one table"
    )
    parser.add_argument(
        '--prefix-type',
        choices=['GN', 'LF', 'all'],
        default='GN',
        help='Combination type: GN, LF, or all to write both (default: GN)'
    )
    parser.add_argument(
        '--combinations-file',
        help='CSV with combinations (poses, scores, distances). Default path follows prefix-type if omitted'
    )
    parser.add_argument(
        '--rmsd-fragments-file',
        help='CSV with per-fragment RMSD. Default path follows prefix-type if omitted'
    )
    parser.add_argument(
        '--md-rmsd-file',
        default='resultados_rmsd_md.csv',
        help='CSV with MD RMSD (default: resultados_rmsd_md.csv)'
    )
    parser.add_argument(
        '--score-only-file',
        help='CSV with score_only (optional; default: score_only_results.csv)'
    )
    parser.add_argument(
        '--mmpbsa-file',
        help='CSV with MM/PBSA results (optional; default: mmpbsa_results.csv)'
    )
    parser.add_argument(
        '--fragment-energies-file',
        help='CSV with per-fragment energy contributions (optional)'
    )
    parser.add_argument(
        '--output',
        help='Output CSV path. Default follows prefix-type if omitted'
    )
    
    args = parser.parse_args()
    
    # Default paths by type
    if not args.combinations_file:
        args.combinations_file = f'valid_fragment_combinations_{args.prefix_type}_no_overlap.csv'
    
    if not args.rmsd_fragments_file:
        args.rmsd_fragments_file = f'valid_combinations_{args.prefix_type}/valid_no_overlap/rmsd_results.csv'
    
    if not args.output:
        args.output = f'all_metrics_{args.prefix_type}.csv'
    
    # Default score_only file
    if not args.score_only_file:
        args.score_only_file = 'score_only_results.csv'
    # Default MM/PBSA file
    if not args.mmpbsa_file:
        args.mmpbsa_file = 'mmpbsa_results.csv'
    
    # Resolve alternate paths when files live under a case folder
    def resolve_path(path: str, kind: str) -> str:
        if Path(path).exists():
            return path
        # Infer case id from --output (e.g. all_metrics_GN_1AQD_C.csv -> 1AQD_C)
        case_from_output = None
        if args.output:
            stem = Path(args.output).stem
            for pfx in ("GN_", "LF_"):
                if pfx in stem:
                    case_from_output = stem.split(pfx)[-1].strip()
                    break
        if case_from_output:
            if kind == "combinations":
                alts = [
                    Path(case_from_output) / Path(path).name,
                    Path(case_from_output) / "valid_fragment_combinations.csv",
                ]
            elif kind == "rmsd_fragments":
                alts = [
                    Path(f"valid_combinations_{args.prefix_type}_{case_from_output}") / "valid_no_overlap" / "rmsd_results.csv",
                ]
            elif kind == "md_rmsd":
                alts = [
                    Path(f"resultados_rmsd_md_{case_from_output}.csv"),
                ]
            else:
                alts = []
            for p in alts:
                if Path(p).exists():
                    print(f"Using {kind}: {p}")
                    return str(p)
        return path

    args.combinations_file = resolve_path(args.combinations_file, "combinations")
    args.rmsd_fragments_file = resolve_path(args.rmsd_fragments_file, "rmsd_fragments")
    args.md_rmsd_file = resolve_path(args.md_rmsd_file, "md_rmsd")

    # Require core files (score_only is optional)
    missing_files = []
    if not Path(args.combinations_file).exists():
        missing_files.append(args.combinations_file)
    if not Path(args.rmsd_fragments_file).exists():
        missing_files.append(args.rmsd_fragments_file)
    if not Path(args.md_rmsd_file).exists():
        missing_files.append(args.md_rmsd_file)
    
    if missing_files:
        print("Error: The following files were not found:")
        for f in missing_files:
            print(f"  - {f}")
        print("  Hint: run step 14 to generate valid_fragment_combinations_GN_no_overlap.csv")
        sys.exit(1)
    
    # score_only: optional; warn if missing
    if not Path(args.score_only_file).exists():
        print(f"⚠ Warning: score_only file not found: {args.score_only_file}")
        print("  Continuing without score_only...")
        args.score_only_file = None

    # MM/PBSA: optional
    if not Path(args.mmpbsa_file).exists():
        print(f"⚠ Warning: MM/PBSA file not found: {args.mmpbsa_file}")
        print("  Continuing without MM/PBSA...")
        args.mmpbsa_file = None

    # Fragment energies (optional)
    fragment_energies_file = getattr(args, 'fragment_energies_file', None) or None
    
    # Run merge for the requested type
    merge_metrics_single_type(
        args.combinations_file,
        args.rmsd_fragments_file,
        args.md_rmsd_file,
        args.output,
        args.prefix_type,
        args.score_only_file,
        args.mmpbsa_file,
        getattr(args, 'fragment_energies_file', None),
    )
    
    # If 'all', also write LF output
    if args.prefix_type == 'all':
        print("\n" + "=" * 70)
        print("Also generating LF output file...")
        print("=" * 70)
        
        lf_combinations_file = 'valid_fragment_combinations_LF_no_overlap.csv'
        lf_rmsd_fragments_file = 'valid_combinations_LF/valid_no_overlap/rmsd_results.csv'
        lf_output = 'all_metrics_LF.csv'
        
        merge_metrics_single_type(
            lf_combinations_file,
            lf_rmsd_fragments_file,
            args.md_rmsd_file,
            lf_output,
            'LF',
            args.score_only_file,
            args.mmpbsa_file,
            getattr(args, 'fragment_energies_file', None),
        )


if __name__ == "__main__":
    main()

