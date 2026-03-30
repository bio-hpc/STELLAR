#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Save coordinates of backbone C in the last residue and N in the first residue
for each generated pose.

Requires: Python 3.4+

Usage:
    python3 STELLAR/save_pose_CN_coordinates.py <pose_file> [--output-dir <dir>]
    python3 STELLAR/save_pose_CN_coordinates.py <directory> [--output-dir <dir>] [--format pdbqt|mol2]

May be called from scriptGN.sh and scriptLF.sh after pose generation.
"""

import os
import sys
import argparse
import csv
import re
from pathlib import Path


def read_pdbqt(file_path):
    """Read a PDBQT file and return a list of atom records.
    Supports:
    1. With chain: ATOM id name resname chain resnum x y z ...
    2. Without chain: ATOM id name resname resnum x y z ...
    """
    atoms = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.split()
                if len(parts) >= 8:
                    atom_id = int(parts[1])
                    atom_name = parts[2]
                    residue_name = parts[3]
                    
                    # Chain: if parts[4] is not numeric and parts[5] is numeric, parts[4] is chain
                    has_chain = False
                    if len(parts) >= 9:
                        try:
                            float(parts[4])
                            has_chain = False
                        except ValueError:
                            try:
                                float(parts[5])
                                has_chain = True
                            except ValueError:
                                has_chain = False
                    
                    try:
                        if has_chain and len(parts) >= 9:
                            # With chain
                            residue_number = parts[5]
                            x = float(parts[6])
                            y = float(parts[7])
                            z = float(parts[8])
                        else:
                            # Without chain
                            residue_number = parts[4]
                            x = float(parts[5])
                            y = float(parts[6])
                            z = float(parts[7])
                    except (ValueError, IndexError) as e:
                        # Alternate layout
                        try:
                            residue_number = parts[4] if len(parts) > 4 else "1"
                            x = float(parts[5])
                            y = float(parts[6])
                            z = float(parts[7])
                        except (ValueError, IndexError):
                            continue  # skip unparseable line
                    
                    atoms.append((atom_id, atom_name, x, y, z, residue_number, residue_name))
    return atoms


def read_mol2(file_path):
    """Read a MOL2 file and return a list of atom records."""
    atoms = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        in_atom_section = False
        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                continue
            elif line.startswith("@<TRIPOS>"):
                in_atom_section = False
                continue
            
            if in_atom_section and line.strip():
                parts = line.split()
                if len(parts) >= 6:
                    atom_id = int(parts[0])
                    atom_name = parts[1]
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                    atom_type = parts[5]
                    residue_number = parts[6] if len(parts) > 6 else "1"
                    residue_name = parts[7] if len(parts) > 7 else "UNK"
                    atoms.append((atom_id, atom_name, x, y, z, residue_number, residue_name))
    return atoms


def extract_first_residue(atoms):
    """Atoms of the first residue (minimum residue number).
    Same logic as check_overlap_VS_GN.py and check_overlap_VS_LF_Def.py."""
    if not atoms:
        return []
    # Convert the 6th element (index 5) to integer and find the minimum value
    try:
        min_value = min(int(t[5]) for t in atoms)
        # Extract and return the tuples with the minimum 6th element value
        return [t for t in atoms if int(t[5]) == min_value]
    except (ValueError, IndexError):
        return atoms


def extract_last_residue(atoms):
    """Atoms of the last residue (maximum residue number).
    Same logic as check_overlap_VS_GN.py and check_overlap_VS_LF_Def.py."""
    if not atoms:
        return []
    # Convert the 6th element (index 5) to integer and find the maximum value
    try:
        max_value = max(int(t[5]) for t in atoms)
        # Extract and return the tuples with the maximum 6th element value
        return [t for t in atoms if int(t[5]) == max_value]
    except (ValueError, IndexError):
        return atoms


def find_atom_coordinates(residue_atoms, atom_name):
    """Coordinates of a named atom in a residue.
    Same as check_overlap_VS_GN.py / check_overlap_VS_LF_Def.py: exact atom name match."""
    matching_atoms = [(t[2], t[3], t[4]) for t in residue_atoms if t[1] == atom_name]
    if matching_atoms:
        return matching_atoms[0]
    return None


def extract_CN_coordinates(file_path, file_format='auto'):
    """
    C coordinates of the last residue and N of the first residue.

    Args:
        file_path: Pose file (pdbqt or mol2).
        file_format: 'pdbqt', 'mol2', or 'auto'.

    Returns:
        (x_C, y_C, z_C, x_N, y_N, z_N) or None.
    """
    if file_format == 'auto':
        if file_path.endswith('.pdbqt'):
            file_format = 'pdbqt'
        elif file_path.endswith('.mol2'):
            file_format = 'mol2'
        else:
            print(f"Error: cannot infer file format for {file_path}")
            return None
    
    # Read atoms
    try:
        if file_format == 'pdbqt':
            atoms = read_pdbqt(file_path)
        elif file_format == 'mol2':
            atoms = read_mol2(file_path)
        else:
            print(f"Error: unsupported format: {file_format}")
            return None
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None
    
    if not atoms:
        print(f"Warning: no atoms in {file_path}")
        return None
    
    # First and last residue
    first_residue = extract_first_residue(atoms)
    last_residue = extract_last_residue(atoms)
    
    if not first_residue:
        print(f"Warning: first residue not found in {file_path}")
        return None
    
    if not last_residue:
        print(f"Warning: last residue not found in {file_path}")
        return None
    
    # C and N atoms
    coord_C = find_atom_coordinates(last_residue, "C")
    coord_N = find_atom_coordinates(first_residue, "N")
    
    if coord_C is None:
        print(f"Warning: atom C not found in last residue of {file_path}")
        return None
    
    if coord_N is None:
        print(f"Warning: atom N not found in first residue of {file_path}")
        return None
    
    return (coord_C[0], coord_C[1], coord_C[2], coord_N[0], coord_N[1], coord_N[2])


def save_coordinates_to_csv(coords, output_path):
    """
    Write coordinates to CSV.

    Args:
        coords: (x_C, y_C, z_C, x_N, y_N, z_N)
        output_path: Output CSV path.
    """
    header = ['x_C_last_residue', 'y_C_last_residue', 'z_C_last_residue', 
              'x_N_first_residue', 'y_N_first_residue', 'z_N_first_residue']
    
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerow(coords)


def extract_pose_number(filename):
    """Parse pose index from filename."""
    patterns = [
        r'_(\d+)\.(pdbqt|mol2)$',
        r'Pose(\d+)',              # Pose1
        r'pose(\d+)',               # pose1
        r'(\d+)\.(pdbqt|mol2)$',    # 1.pdbqt
    ]
    
    for pattern in patterns:
        match = re.search(pattern, filename, re.IGNORECASE)
        if match:
            return int(match.group(1))
    
    return None


def process_single_file(file_path, output_dir=None, pose_number=None):
    """
    Process one pose file and write coordinates.

    Args:
        file_path: Pose file path.
        output_dir: CSV directory (None = same as file).
        pose_number: Pose index (None = parse from filename).
    """
    file_path = Path(file_path)
    
    if not file_path.exists():
        print(f"Error: file not found: {file_path}")
        return False
    
    coords = extract_CN_coordinates(str(file_path))
    if coords is None:
        return False
    
    if pose_number is None:
        pose_number = extract_pose_number(file_path.name)
        if pose_number is None:
            pose_number = 1
    
    # Output directory
    if output_dir is None:
        output_dir = file_path.parent
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    # CSV filename
    csv_filename = f"Pose{pose_number}_coordinates_CN.csv"
    csv_path = output_dir / csv_filename
    
    save_coordinates_to_csv(coords, csv_path)
    print(f"✓ Saved: {csv_path}")
    print(f"  Coordinates: C({coords[0]:.3f}, {coords[1]:.3f}, {coords[2]:.3f}), "
          f"N({coords[3]:.3f}, {coords[4]:.3f}, {coords[5]:.3f})")
    
    return True


def process_directory(directory, file_format='auto', output_dir=None):
    """
    Process all pose files in a directory.

    Args:
        directory: Directory with pose files.
        file_format: 'pdbqt', 'mol2', or 'auto'.
        output_dir: CSV output directory (None = same as poses).
    """
    directory = Path(directory)
    
    if not directory.exists() or not directory.is_dir():
        print(f"Error: not a directory or missing: {directory}")
        return
    
    # Extensions by format
    if file_format == 'auto':
        extensions = ['.pdbqt', '.mol2']
    elif file_format == 'pdbqt':
        extensions = ['.pdbqt']
    elif file_format == 'mol2':
        extensions = ['.mol2']
    else:
        print(f"Error: invalid format: {file_format}")
        return
    
    # Collect pose files
    pose_files = []
    for ext in extensions:
        pose_files.extend(directory.glob(f"*{ext}"))
        # auto: also match by extension pattern
        if file_format == 'auto':
            pose_files.extend(directory.glob(f"*.{ext[1:]}"))
    
    # Skip temp files named with "ligand"
    pose_files = [f for f in pose_files if "ligand" not in f.name.lower()]
    
    if not pose_files:
        print(f"No pose files found in {directory}")
        return
    
    pose_files.sort()
    
    print(f"Processing {len(pose_files)} pose file(s)...")
    
    success_count = 0
    for pose_file in pose_files:
        if process_single_file(pose_file, output_dir):
            success_count += 1
    
    print(f"\n✓ Processed {success_count}/{len(pose_files)} file(s) successfully")


def main():
    parser = argparse.ArgumentParser(
        description='Extract and save C (last residue) and N (first residue) coordinates from poses'
    )
    parser.add_argument('input', help='Pose file or directory of pose files')
    parser.add_argument('--output-dir', '-o', help='Directory for CSV output (default: same as input)')
    parser.add_argument('--format', '-f', choices=['pdbqt', 'mol2', 'auto'], default='auto',
                       help='File format (pdbqt, mol2, or auto)')
    
    args = parser.parse_args()
    
    input_path = Path(args.input)
    
    if input_path.is_file():
        process_single_file(input_path, args.output_dir)
    elif input_path.is_dir():
        process_directory(input_path, args.format, args.output_dir)
    else:
        print(f"Error: {args.input} is not a valid file or directory")
        sys.exit(1)


if __name__ == "__main__":
    main()

