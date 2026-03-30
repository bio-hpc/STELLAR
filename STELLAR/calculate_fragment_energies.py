#!/usr/bin/env python3
"""
Extract energy terms (graph_global_field / graph_global_score) for each pose
listed in the combinations CSV, from JSON files under VS_GN_*/energies/.

For each combination and fragment, reads the pose JSON and collects:
Gauss1, Gauss2, Repulsion, Hydrophobic, Hydrogen_Bonds, Rotational, Total_Affinity.

Output: CSV with combination_id and frag{N}_{field} columns (e.g. frag1_Gauss1, frag1_Total_Affinity).
"""

import argparse
import csv
import glob
import json
import os
import re
import sys
from typing import Dict, Optional, Tuple


# Expected keys in graph_global_field / graph_global_score
ENERGY_FIELDS = [
    "Gauss1", "Gauss2", "Repulsion", "Hydrophobic",
    "Hydrogen_Bonds", "Rotational", "Total_Affinity"
]


def find_json_for_pose(energies_dir: str, pose_number: int) -> Optional[str]:
    """
    Find a JSON in energies_dir matching the pose (suffix _<pose>.json).
    """
    if not os.path.isdir(energies_dir):
        return None
    pattern = os.path.join(energies_dir, f"*_{pose_number}.json")
    matches = glob.glob(pattern)
    return matches[0] if matches else None


def extract_global_scores(json_path: str) -> Optional[Dict[str, str]]:
    """
    Read graph_global_field and graph_global_score from JSON;
    returns {field_name: value_str}.
    """
    try:
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        return None
    fields = data.get("graph_global_field") or []
    scores = data.get("graph_global_score") or []
    if len(fields) != len(scores):
        return None
    return dict(zip(fields, scores))


def resolve_energies_root(
    base_dir: str,
    case_upper: str,
    prefix_type: str,
) -> Optional[Tuple[str, str]]:
    """
    Detect which complex folder exists on disk and return (parent_folder, effective_case)
    for consistent use across the run (avoids mixing A/B chains of the same PDB).

    Search order:
    1) {case}_GN (e.g. 3C3O_A_GN) with VS_{prefix}_{case}_Frag1/energies inside
    2) Same PDB, other chains (3C3O_B_GN, etc.) when case looks like PDB_X
    3) {case}/ (no _GN)
    4) Flat layout: VS_{prefix}_{case}_Frag1/energies under base_dir

    Returns:
        (parent_folder, effective_case) relative to base_dir, or None.
        parent_folder may be "" for a flat layout.
    """
    prefix = "VS_GN" if prefix_type.upper() == "GN" else "VS_LF"
    case_upper = case_upper.strip().upper()

    def has_energies(parent: str, case_val: str) -> bool:
        frag1 = f"{prefix}_{case_val}_Frag1"
        path = os.path.join(base_dir, parent, frag1, "energies") if parent else os.path.join(base_dir, frag1, "energies")
        return os.path.isdir(path)

    # 1) Case as given with _GN (priority: --case value)
    if has_energies(f"{case_upper}_GN", case_upper):
        return (f"{case_upper}_GN", case_upper)
    if has_energies(case_upper, case_upper):
        return (case_upper, case_upper)
    if has_energies("", case_upper):
        return ("", case_upper)

    # 2) Same PDB, other single-letter chains: try A, B, C
    if "_" in case_upper and len(case_upper.split("_")[-1]) == 1:
        base_pdb = case_upper.rsplit("_", 1)[0]
        for ch in "ABC":
            alt_case = f"{base_pdb}_{ch}"
            if has_energies(f"{alt_case}_GN", alt_case):
                return (f"{alt_case}_GN", alt_case)
            if has_energies(alt_case, alt_case):
                return (alt_case, alt_case)
            if has_energies("", alt_case):
                return ("", alt_case)

    return None


def get_fragment_energies_row(
    case: str,
    prefix_type: str,
    combo_row: dict,
    base_dir: str,
    max_frag: int,
    resolved_root: Optional[Tuple[str, str]] = None,
) -> Dict[str, str]:
    """
    For one combination row, fill a dict with combination_id and
    frag{N}_{field} for each fragment and energy field.

    If resolved_root is set, only that (parent_folder, effective_case) is used,
    avoiding mixed A/B chains for the same complex.
    """
    row = {"combination_id": combo_row.get("combination_id", "")}
    case_upper = case.upper().strip()
    prefix = "VS_GN" if prefix_type.upper() == "GN" else "VS_LF"

    for frag_num in range(1, max_frag + 1):
        pose_key = f"frag{frag_num}_pose"
        pose_str = combo_row.get(pose_key, "")
        if not pose_str:
            for field in ENERGY_FIELDS:
                row[f"frag{frag_num}_{field}"] = ""
            continue
        try:
            pose_number = int(pose_str)
        except ValueError:
            for field in ENERGY_FIELDS:
                row[f"frag{frag_num}_{field}"] = ""
            continue

        if resolved_root:
            parent_folder, effective_case = resolved_root
            frag_folder = f"{prefix}_{effective_case}_Frag{frag_num}"
            energies_dir = os.path.join(base_dir, parent_folder, frag_folder, "energies") if parent_folder else os.path.join(base_dir, frag_folder, "energies")
            json_path = find_json_for_pose(energies_dir, pose_number)
        else:
            frag_folder_tpl = f"{prefix}_{case_upper}_Frag{{}}"
            frag_folder = frag_folder_tpl.format(frag_num)
            candidates = [
                os.path.join(base_dir, case_upper, frag_folder, "energies"),
                os.path.join(base_dir, f"{case_upper}_GN", frag_folder, "energies"),
                os.path.join(base_dir, frag_folder, "energies"),
            ]
            if "_" in case_upper and len(case_upper.split("_")[-1]) == 1:
                base_pdb = case_upper.rsplit("_", 1)[0]
                for alt_chain in "ABC":
                    if alt_chain == case_upper.rsplit("_", 1)[1]:
                        continue
                    alt_case = f"{base_pdb}_{alt_chain}"
                    alt_folder = f"{prefix}_{alt_case}_Frag{frag_num}"
                    candidates.append(os.path.join(base_dir, f"{alt_case}_GN", alt_folder, "energies"))
                    candidates.append(os.path.join(base_dir, alt_case, alt_folder, "energies"))
            json_path = None
            for energies_dir in candidates:
                json_path = find_json_for_pose(energies_dir, pose_number)
                if json_path:
                    break
        if not json_path:
            for field in ENERGY_FIELDS:
                row[f"frag{frag_num}_{field}"] = ""
            continue
        scores = extract_global_scores(json_path)
        if not scores:
            for field in ENERGY_FIELDS:
                row[f"frag{frag_num}_{field}"] = ""
            continue
        for field in ENERGY_FIELDS:
            row[f"frag{frag_num}_{field}"] = scores.get(field, "")
    return row


def main():
    parser = argparse.ArgumentParser(
        description="Extract per-fragment energy terms from JSON under VS_GN_*/energies/"
    )
    parser.add_argument("--case", required=True, help="Case (e.g. 1AQD_C, 1CJR_B)")
    parser.add_argument("--type", choices=["GN", "LF"], default="GN", help="GN or LF (default: GN)")
    parser.add_argument(
        "--combinations-file",
        default="valid_fragment_combinations_GN_no_overlap.csv",
        help="CSV with combination_id and frag{N}_pose (default: valid_fragment_combinations_GN_no_overlap.csv)",
    )
    parser.add_argument(
        "--output-csv",
        default="fragment_energies_results.csv",
        help="Output CSV (default: fragment_energies_results.csv)",
    )
    parser.add_argument(
        "--base-dir",
        default=".",
        help="Project root containing {case}/VS_GN_* (default: .)",
    )
    args = parser.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    combos_path = os.path.join(base_dir, args.combinations_file) if not os.path.isabs(args.combinations_file) else args.combinations_file
    if not os.path.isfile(combos_path):
        # Alternate paths (case subfolder or cwd)
        case_dir = args.case.strip().upper()
        no_overlap_name = "valid_fragment_combinations_GN_no_overlap.csv" if "GN" in args.type.upper() else "valid_fragment_combinations_LF_no_overlap.csv"
        alternatives = [
            os.path.join(base_dir, args.combinations_file),
            os.path.join(base_dir, case_dir, args.combinations_file),
            args.combinations_file,
            os.path.join(case_dir, args.combinations_file),
            os.path.join(base_dir, case_dir, "valid_fragment_combinations.csv"),
            os.path.join(case_dir, "valid_fragment_combinations.csv"),
        ]
        combos_path = None
        for p in alternatives:
            if os.path.isfile(p):
                combos_path = p
                if no_overlap_name not in p and "valid_fragment_combinations.csv" in p:
                    print(f"⚠ Using unfiltered CSV (not no_overlap): {p}")
                    print("  For overlap-filtered combinations only, run step 14 then re-run.")
                else:
                    print(f"Using combinations: {p}")
                break
        if combos_path is None:
            print(f"Error: combinations CSV not found.")
            print(f"  Looked for: {os.path.join(base_dir, args.combinations_file)}")
            print(f"  Run step 14 first: python3 STELLAR/filter_valid_combinations_csv.py GN --combinations-dir valid_combinations_GN_{case_dir}/valid_no_overlap --input-csv {case_dir}/valid_fragment_combinations.csv")
            sys.exit(1)

    with open(combos_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        combinations = list(reader)

    if not combinations:
        print("No rows in combinations CSV.")
        sys.exit(0)

    max_frag = 0
    for row in combinations:
        for key in row:
            m = re.match(r"frag(\d+)_pose", key)
            if m:
                max_frag = max(max_frag, int(m.group(1)))
    if max_frag == 0:
        print("No frag{N}_pose columns found in CSV.")
        sys.exit(1)

    case_dir = args.case.strip().upper()
    resolved_root = resolve_energies_root(base_dir, case_dir, args.type)
    if resolved_root is None:
        print("Error: no energies folder found for this case.")
        print(f"  Case: {case_dir}")
        print(f"  Base directory: {base_dir}")
        print("  Expected layouts: {case}_GN/VS_GN_{case}_Frag1/energies or VS_GN_{case}_Frag1/energies")
        print("  (Chains A, B, C are tried for the same PDB when applicable.)")
        sys.exit(1)
    parent_folder, effective_case = resolved_root
    if parent_folder:
        print(f"Using energies under: {parent_folder}/ (effective case: {effective_case})")
    else:
        print(f"Using energies at tree root (effective case: {effective_case})")

    fieldnames = ["combination_id"]
    for i in range(1, max_frag + 1):
        for field in ENERGY_FIELDS:
            fieldnames.append(f"frag{i}_{field}")

    results = []
    for combo_row in combinations:
        row = get_fragment_energies_row(
            args.case,
            args.type,
            combo_row,
            base_dir,
            max_frag,
            resolved_root=resolved_root,
        )
        results.append(row)

    out_path = os.path.join(base_dir, args.output_csv) if not os.path.isabs(args.output_csv) else args.output_csv
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    print(f"✓ Wrote {len(results)} rows to {out_path}")


if __name__ == "__main__":
    main()
