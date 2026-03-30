#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Aggregate C/N coordinates for all GN poses and attach docking scores.

Writes one CSV per VS_GN_* folder:
    <VS_GN_dir>/<VS_GN_name>_coordinates_CN_with_score.csv

Columns:
    pose, x_C_last_residue, y_C_last_residue, z_C_last_residue,
    x_N_first_residue, y_N_first_residue, z_N_first_residue,
    DockingScore

Score sources:
    energies/<pose>.en  -> line "Affinity: <value>"
    (fallback) energies/<pose>.json -> "global_score"
"""
from pathlib import Path
import csv
import json
import re

# Local copy in STELLAR/ (same functions as MetaScreener used_by_metascreener)
from save_pose_CN_coordinates import extract_CN_coordinates, extract_pose_number  # type: ignore


def safe_pose_number(name: str):
    """
    Strict pose index: last _<digits>.<ext> only (avoids decimals like 15.774).
    """
    m = re.search(r"_(\d+)\.(pdbqt|mol2)$", name, re.IGNORECASE)
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            pass
    # Fallback al extractor previo (por compatibilidad)
    return extract_pose_number(name)


def parse_affinity_from_en(en_path: Path):
    """Parse Affinity from a .en file."""
    try:
        with en_path.open() as f:
            for line in f:
                if line.startswith("Affinity:"):
                    parts = line.split()
                    if len(parts) >= 2:
                        return float(parts[1])
    except Exception:
        pass
    return None


def parse_score_from_json(json_path: Path):
    """Fallback: read global_score from JSON."""
    try:
        data = json.loads(json_path.read_text())
        if "global_score" in data:
            return float(data["global_score"])
    except Exception:
        pass
    return None


def process_vs_gn_dir(vs_dir: Path):
    molecules_dir = vs_dir / "molecules"
    energies_dir = vs_dir / "energies"

    if not molecules_dir.is_dir():
        print(f"[WARN] missing molecules/ under {vs_dir}")
        return

    poses_rows = []
    poses = {}

    # Valid pose indices from _<n>.pdbqt filenames
    valid_poses = set()
    for pdbqt_file in molecules_dir.glob("*.pdbqt"):
        if "ligand" in pdbqt_file.name.lower():
            continue
        m = re.search(r"_(\d+)\.pdbqt$", pdbqt_file.name, re.IGNORECASE)
        if m:
            try:
                valid_poses.add(int(m.group(1)))
            except ValueError:
                pass

    # -- Caso 1: usar directamente los PoseX_coordinates_CN.csv existentes (evita parsear nombres raros)
    pose_csv_files = sorted(molecules_dir.glob("Pose*_coordinates_CN.csv"))
    use_pose_csv = bool(pose_csv_files)

    if use_pose_csv:
        for csv_file in pose_csv_files:
            m = re.search(r"Pose(\d+)_coordinates_CN\.csv$", csv_file.name, re.IGNORECASE)
            if not m:
                continue
            pose_num = int(m.group(1))
            # Only valid pose indices
            if pose_num not in valid_poses:
                continue
            # First data row
            try:
                with csv_file.open() as f:
                    reader = csv.reader(f)
                    next(reader, None)  # header
                    row = next(reader, None)
                    if not row or len(row) < 6:
                        continue
                    coords = tuple(float(x) for x in row[:6])
            except Exception:
                continue

            # Score from energies/*_<pose>*.en or .json
            score = None
            for en_path in energies_dir.glob(f"*_{pose_num}*.en"):
                score = parse_affinity_from_en(en_path)
                if score is not None:
                    break
            if score is None:
                for js in energies_dir.glob(f"*_{pose_num}*.json"):
                    score = parse_score_from_json(js)
                    if score is not None:
                        break

            poses_rows.append((pose_num,) + coords + (score,))

    # Case 2: build from pdbqt if no Pose CSVs
    if not use_pose_csv:
        pose_files = sorted(p for p in molecules_dir.glob("*.pdbqt") if "ligand" not in p.name.lower())

        for pose_file in pose_files:
            pose_num = safe_pose_number(pose_file.name)
            if pose_num is None:
                continue
            # Only valid pose indices
            if pose_num not in valid_poses:
                continue

            coords = extract_CN_coordinates(str(pose_file))
            if coords is None:
                continue  # no coords, descartar

            base_name = pose_file.stem  # sin .pdbqt
            en_path = energies_dir / f"{base_name}.en"
            json_path = energies_dir / f"{base_name}.json"

            score = None
            if en_path.exists():
                score = parse_affinity_from_en(en_path)
            if score is None and json_path.exists():
                score = parse_score_from_json(json_path)

            current = poses.get(pose_num)
            # Prefer row with score; keep first if both have scores
            if current is None:
                poses[pose_num] = {"coords": coords, "score": score}
            else:
                if current["score"] is None and score is not None:
                    poses[pose_num] = {"coords": coords, "score": score}

    # Final rows
    if poses_rows:
        rows = sorted(poses_rows, key=lambda r: r[0])
    else:
        rows = []
        for pose_num in sorted(poses.keys()):
            coords = poses[pose_num]["coords"]
            score = poses[pose_num]["score"]
            rows.append((pose_num,) + coords + (score,))

    # Output at VS_GN root (parent of molecules/)
    prefix = vs_dir.name
    out_csv = vs_dir / f"{prefix}_coordinates_CN_with_score.csv"
    header = [
        "pose",
        "x_C_last_residue",
        "y_C_last_residue",
        "z_C_last_residue",
        "x_N_first_residue",
        "y_N_first_residue",
        "z_N_first_residue",
        "DockingScore",
    ]
    with out_csv.open("w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)

    print(f"[OK] {vs_dir.name}: {len(rows)} poses -> {out_csv}")


def find_vs_gn_dirs(root: Path, case_filter=None):
    """Find VS_GN_* dirs, optionally filtered by case string."""
    all_dirs = [p for p in root.iterdir() if p.is_dir() and p.name.startswith("VS_GN_")]
    if case_filter:
        case_filter_upper = case_filter.upper()
        return [d for d in all_dirs if case_filter_upper in d.name.upper()]
    return all_dirs


def main(case_filter=None):
    import argparse

    parser = argparse.ArgumentParser(
        description="Aggregate C/N coordinates and docking score for VS_GN_* folders"
    )
    parser.add_argument(
        "dirs",
        nargs="*",
        help="Specific VS_GN folders to process; if omitted, all VS_GN_* in cwd.",
    )
    parser.add_argument(
        "--case-filter",
        type=str,
        default=None,
        help="Case substring (e.g. 2W10_C) to restrict folders",
    )
    args = parser.parse_args()

    root = Path.cwd()
    if args.dirs:
        targets = []
        for d in args.dirs:
            p = (root / d).resolve()
            if p.is_dir():
                if p.name.startswith("VS_GN_"):
                    targets.append(p)
                else:
                    # Container: search VS_GN_* inside
                    vs_gn_dirs = find_vs_gn_dirs(p, case_filter=args.case_filter)
                    targets.extend(vs_gn_dirs)
            else:
                print(f"[WARN] Folder not found: {d}")
    else:
        targets = find_vs_gn_dirs(root, case_filter=args.case_filter)
        if args.case_filter:
            excluded = [d for d in find_vs_gn_dirs(root) if d not in targets]
            if excluded:
                print(f"✓ Filtrado: {len(targets)} carpetas de {args.case_filter}, {len(excluded)} excluidas")

    if not targets:
        print("No VS_GN folders to process.")
        return

    for vs_dir in targets:
        process_vs_gn_dir(vs_dir)


if __name__ == "__main__":
    main()

