#!/usr/bin/env python3
"""
Filter the valid-combinations CSV using folders that passed the overlap filter.

Extracts from the original CSV only rows whose combinations exist under valid_no_overlap/.
"""

import os
import csv
import glob
import argparse


def get_valid_combination_ids(valid_combinations_dir):
    """
    Collect valid combination IDs from folder names.

    Args:
        valid_combinations_dir: Directory containing combination_* folders.

    Returns:
        set of valid combination IDs.
    """
    valid_ids = set()
    
    if not os.path.exists(valid_combinations_dir):
        print(f"⚠ Warning: directory does not exist: {valid_combinations_dir}")
        return valid_ids
    
    for combo_dir in glob.glob(os.path.join(valid_combinations_dir, "combination_*")):
        combo_id = os.path.basename(combo_dir).replace("combination_", "")
        if combo_id.isdigit():
            valid_ids.add(int(combo_id))
    
    return valid_ids


def filter_csv(input_csv, output_csv, valid_ids, id_column='combination_id'):
    """
    Filter the input CSV to rows whose combination ID is in valid_ids.

    Args:
        input_csv: Path to the source CSV.
        output_csv: Path to the output CSV.
        valid_ids: Set of valid combination IDs.
        id_column: Column name for the combination ID.

    Returns:
        Number of rows written.
    """
    if not os.path.exists(input_csv):
        print(f"✗ Error: file not found: {input_csv}")
        return 0
    
    # Read source CSV and filter
    valid_rows = []
    with open(input_csv, 'r') as f_in:
        reader = csv.DictReader(f_in)
        header = reader.fieldnames
        
        if id_column not in header:
            print(f"✗ Error: column '{id_column}' not found in CSV")
            print(f"  Available columns: {', '.join(header)}")
            return 0
        
        for row in reader:
            try:
                combo_id = int(row[id_column])
                if combo_id in valid_ids:
                    valid_rows.append(row)
            except (ValueError, KeyError) as e:
                print(f"⚠ Warning: error processing row: {e}")
                continue
    
    # Write filtered CSV
    if valid_rows:
        with open(output_csv, 'w', newline='') as f_out:
            writer = csv.DictWriter(f_out, fieldnames=header)
            writer.writeheader()
            writer.writerows(valid_rows)
        
        print(f"✓ Wrote {len(valid_rows)} valid combinations to {output_csv}")
        return len(valid_rows)
    else:
        print(f"⚠ Warning: no valid combinations found")
        return 0


def main():
    parser = argparse.ArgumentParser(
        description="Filter valid-combinations CSV using valid_no_overlap folders"
    )
    parser.add_argument(
        'type',
        choices=['GN', 'LF'],
        help='Combination type: GN or LF'
    )
    parser.add_argument(
        '--input-csv',
        help='Input CSV (default: valid_fragment_combinations_{type}.csv)'
    )
    parser.add_argument(
        '--output-csv',
        help='Output CSV (default: valid_fragment_combinations_{type}_no_overlap.csv)'
    )
    parser.add_argument(
        '--combinations-dir',
        help='Combinations directory (default: valid_combinations_{type}/valid_no_overlap)'
    )
    parser.add_argument(
        '--id-column',
        default='combination_id',
        help='Column name for combination ID (default: combination_id)'
    )
    
    args = parser.parse_args()
    
    # Defaults
    if args.input_csv is None:
        args.input_csv = "valid_fragment_combinations.csv"
    
    if args.output_csv is None:
        args.output_csv = f"valid_fragment_combinations_{args.type}_no_overlap.csv"
    
    if args.combinations_dir is None:
        args.combinations_dir = f"valid_combinations_{args.type}/valid_no_overlap"
    
    print("=" * 70)
    print(f"Filtering valid combinations for {args.type}")
    print("=" * 70)
    print(f"Combinations directory: {args.combinations_dir}")
    print(f"Input CSV: {args.input_csv}")
    print(f"Output CSV: {args.output_csv}")
    print()
    
    # Collect valid IDs
    print("Looking for valid combinations...")
    valid_ids = get_valid_combination_ids(args.combinations_dir)
    
    if not valid_ids:
        print(f"✗ No valid combinations found in {args.combinations_dir}")
        return 1
    
    print(f"✓ Found {len(valid_ids)} valid combinations")
    
    # Filter CSV
    print(f"\nFiltering CSV...")
    num_written = filter_csv(args.input_csv, args.output_csv, valid_ids, args.id_column)
    
    if num_written == 0:
        return 1
    
    print("\n" + "=" * 70)
    print(f"✓ Done: {num_written} valid combinations written")
    print("=" * 70)
    
    return 0


if __name__ == "__main__":
    exit(main())


