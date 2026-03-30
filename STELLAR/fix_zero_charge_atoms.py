#!/usr/bin/env python3
"""
Detect and fix atoms with exactly 0.0000 charge before building GROMACS topologies.

Issue: acpype assigns type "DU" (dummy) with zero mass to atoms it cannot map to the
AMBER force field, especially aromatic carbons at 0.0000 charge.

Fix: Adjust 0.0000 charges by spreading to bonded neighbors or assigning a small
minimum so acpype does not label them as dummy.
"""

import argparse
import sys
import os
from pathlib import Path
import re
from collections import defaultdict


def parse_mol2_file(mol2_path):
    """Parse MOL2 file and return structured data."""
    atoms = []
    bonds = []
    header_lines = []
    footer_lines = []  # Lines after @<TRIPOS>BOND
    molecule_info = {}
    
    in_atom_section = False
    in_bond_section = False
    in_molecule_section = False
    past_bond_section = False  # True after BOND section
    
    with open(mol2_path, 'r') as f:
        for line in f:
            if line.startswith("@<TRIPOS>MOLECULE"):
                in_molecule_section = True
                in_atom_section = False
                in_bond_section = False
                past_bond_section = False
                header_lines.append(line)
                continue
            elif line.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                in_molecule_section = False
                in_bond_section = False
                past_bond_section = False
                header_lines.append(line)
                continue
            elif line.startswith("@<TRIPOS>BOND"):
                in_bond_section = True
                in_atom_section = False
                in_molecule_section = False
                past_bond_section = False
                header_lines.append(line)
                continue
            elif line.startswith("@<TRIPOS>"):
                # Section after BOND (e.g. SUBSTRUCTURE)
                in_atom_section = False
                in_bond_section = False
                in_molecule_section = False
                past_bond_section = True
                footer_lines.append(line)
                continue
            
            if in_molecule_section:
                header_lines.append(line)
                # Parse molecule info
                parts = line.strip().split()
                if len(parts) >= 3:
                    try:
                        molecule_info['natoms'] = int(parts[0])
                        molecule_info['nbonds'] = int(parts[1])
                    except:
                        pass
            elif in_atom_section:
                # Parse atom line: atom_id atom_name x y z atom_type subst_id subst_name charge
                parts = line.strip().split()
                if len(parts) >= 8:
                    try:
                        atom_id = int(parts[0])
                        atom_name = parts[1]
                        x = float(parts[2])
                        y = float(parts[3])
                        z = float(parts[4])
                        atom_type = parts[5]
                        subst_id = int(parts[6])
                        subst_name = parts[7]
                        charge = float(parts[8])
                        
                        atoms.append({
                            'id': atom_id,
                            'name': atom_name,
                            'x': x, 'y': y, 'z': z,
                            'type': atom_type,
                            'subst_id': subst_id,
                            'subst_name': subst_name,
                            'charge': charge,
                            'line': line.rstrip(),  # keep original line for formatting
                            'original_charge': charge
                        })
                    except (ValueError, IndexError):
                        pass
            elif in_bond_section:
                # Parse bond line: bond_id atom1 atom2 bond_type
                parts = line.strip().split()
                if len(parts) >= 3:
                    try:
                        bond_id = int(parts[0])
                        atom1 = int(parts[1])
                        atom2 = int(parts[2])
                        bond_type = parts[3] if len(parts) > 3 else "1"
                        bonds.append({
                            'id': bond_id,
                            'atom1': atom1,
                            'atom2': atom2,
                            'type': bond_type,
                            'line': line.rstrip()
                        })
                    except (ValueError, IndexError):
                        pass
            elif past_bond_section:
                # Lines after @<TRIPOS>BOND (SUBSTRUCTURE, comments, etc.)
                footer_lines.append(line)
            else:
                header_lines.append(line)
    
    return {
        'atoms': atoms,
        'bonds': bonds,
        'header': header_lines,
        'footer': footer_lines,
        'molecule_info': molecule_info
    }


def build_bond_graph(bonds):
    """Build adjacency graph from bonds."""
    graph = defaultdict(list)
    for bond in bonds:
        a1 = bond['atom1']
        a2 = bond['atom2']
        graph[a1].append(a2)
        graph[a2].append(a1)
    return graph


def fix_zero_charges(atoms, bonds, min_charge=0.0001, max_charge=0.001, fix_aromatic=True):
    """
    Fix atoms with exactly 0.0000 charge by distributing charge to neighbors.
    Also fixes atoms with very small charges (< 0.01) if they are aromatic (C.ar)
    to prevent acpype from assigning DU (dummy) type.
    
    Args:
        atoms: List of atom dictionaries
        bonds: List of bond dictionaries
        min_charge: Minimum charge to assign (default: 0.0001)
        max_charge: Maximum charge adjustment (default: 0.001)
        fix_aromatic: If True, also fix C.ar atoms with charge < 0.01 (default: True)
    
    Returns:
        Modified atoms list and list of fixes applied
    """
    # Build bond graph
    bond_graph = build_bond_graph(bonds)
    
    # Find atoms with zero charge or problematic small charges
    zero_charge_atoms = [a for a in atoms if abs(a['charge']) < 1e-6]
    
    # Also find aromatic carbons with very small charges that acpype might mark as DU
    if fix_aromatic:
        small_charge_aromatic = [a for a in atoms 
                                if a['type'] == 'C.ar' and 0 < abs(a['charge']) < 0.01]
        zero_charge_atoms.extend(small_charge_aromatic)
    
    if not zero_charge_atoms:
        return atoms, []
    
    fixes = []
    modified_atoms = atoms.copy()
    
    # Create atom index map
    atom_dict = {a['id']: a for a in modified_atoms}
    
    # Calculate total charge change needed
    total_charge_change = 0.0
    charge_adjustments = {}
    
    for zero_atom in zero_charge_atoms:
        atom_id = zero_atom['id']
        atom_name = zero_atom['name']
        atom_type = zero_atom['type']
        old_charge = zero_atom['charge']
        
        # Get neighbors
        neighbors = bond_graph.get(atom_id, [])
        
        # For aromatic carbons with small charges, use higher minimum charge
        if atom_type == 'C.ar' and abs(old_charge) < 0.01:
            # Use higher charge for aromatic carbons to avoid DU assignment
            new_charge = 0.01 if old_charge >= 0 else -0.01
            charge_change = new_charge - old_charge
            charge_adjustments[atom_id] = {
                'atom_name': atom_name,
                'atom_type': atom_type,
                'old_charge': old_charge,
                'new_charge': new_charge,
                'method': 'aromatic_small_charge_fix',
                'change': charge_change
            }
            total_charge_change += charge_change
        elif abs(old_charge) < 1e-6:
            # Zero charge atoms
            if not neighbors:
                # No neighbors, assign minimum charge
                new_charge = min_charge
                charge_change = new_charge - old_charge
                charge_adjustments[atom_id] = {
                    'atom_name': atom_name,
                    'atom_type': atom_type,
                    'old_charge': 0.0,
                    'new_charge': new_charge,
                    'method': 'no_neighbors_min_charge',
                    'change': charge_change
                }
                total_charge_change += charge_change
            else:
                # Distribute charge to neighbors
                charge_per_neighbor = min_charge / len(neighbors)
                
                # Adjust neighbor charges
                for neighbor_id in neighbors:
                    if neighbor_id in atom_dict:
                        neighbor = atom_dict[neighbor_id]
                        # Small adjustment to neighbor charge
                        adjustment = charge_per_neighbor * 0.1  # Small adjustment
                        neighbor['charge'] += adjustment
                        total_charge_change += adjustment
                
                # Assign small positive charge to zero-charge atom
                new_charge = min_charge
                charge_change = new_charge - old_charge
                charge_adjustments[atom_id] = {
                    'atom_name': atom_name,
                    'atom_type': atom_type,
                    'old_charge': 0.0,
                    'new_charge': new_charge,
                    'neighbors': len(neighbors),
                    'method': 'distributed_to_neighbors',
                    'change': charge_change
                }
                total_charge_change += charge_change
    
    # Apply charge adjustments and compensate for total charge change
    if total_charge_change != 0.0 and abs(total_charge_change) > 1e-6:
        # Distribute compensation across all non-fixed atoms (prefer atoms with larger absolute charges)
        all_atom_ids = [a['id'] for a in atoms if a['id'] not in charge_adjustments]
        if all_atom_ids:
            # Sort by absolute charge (largest first) to distribute compensation
            sorted_atoms = sorted(all_atom_ids, 
                                key=lambda aid: abs(atom_dict[aid]['charge']), 
                                reverse=True)
            
            # Distribute compensation across top atoms
            compensation_per_atom = -total_charge_change / min(len(sorted_atoms), 20)
            for atom_id in sorted_atoms[:20]:
                if abs(compensation_per_atom) < abs(atom_dict[atom_id]['charge']) * 0.1:  # Don't change more than 10% of original charge
                    atom_dict[atom_id]['charge'] += compensation_per_atom
    
    # Apply all charge adjustments (preserve original atom order)
    for atom_id, adj_info in charge_adjustments.items():
        atom_dict[atom_id]['charge'] = adj_info['new_charge']
        fixes.append({
            'atom_id': atom_id,
            'atom_name': adj_info['atom_name'],
            'atom_type': adj_info['atom_type'],
            'old_charge': adj_info['old_charge'],
            'new_charge': adj_info['new_charge'],
            'method': adj_info['method'],
            'neighbors': adj_info.get('neighbors', None)
        })
    
    # Return atoms in original order (consecutive IDs matter for PyMOL)
    result_atoms = [atom_dict[atom['id']] for atom in atoms]
    return result_atoms, fixes


def write_mol2_file(mol2_path, data, atoms):
    """Write MOL2 file with updated charges, preserving all sections."""
    with open(mol2_path, 'w') as f:
        # Write header up to @<TRIPOS>ATOM (do not duplicate @<TRIPOS>ATOM)
        atom_section_found = False
        for line in data['header']:
            if line.startswith("@<TRIPOS>ATOM"):
                f.write(line)
                atom_section_found = True
                break
            f.write(line)
        
        # If @<TRIPOS>ATOM was not in header, write it now
        if not atom_section_found:
            f.write("@<TRIPOS>ATOM\n")
        
        # Write atoms; preserve line format, replace charge only if changed
        import re
        for atom in atoms:
            original_line = atom.get('line', '')
            original_charge = atom.get('original_charge', atom['charge'])
            new_charge = atom['charge']
            
            if abs(original_charge - new_charge) > 1e-6 and original_line:
                # Replace trailing float (charge) in original line
                pattern = r'([\s]*)(-?\d+\.\d+)(\s*)$'
                match = re.search(pattern, original_line)
                if match:
                    # Replace charge number, keep spacing
                    prefix = original_line[:match.start()]
                    new_line = prefix + match.group(1) + f'{new_charge:7.4f}' + match.group(3)
                    f.write(new_line + '\n')
                else:
                    # Fallback: replace last float
                    pattern2 = r'(-?\d+\.\d+)(\s*)$'
                    new_line = re.sub(pattern2, f'{new_charge:7.4f}\\2', original_line)
                    f.write(new_line + '\n')
            else:
                if original_line:
                    f.write(original_line + '\n')
                else:
                    # Standard MOL2 line if original missing
                    subst_name = atom['subst_name']
                    f.write(f"{atom['id']:7d} {atom['name']:<8s} "
                           f"{atom['x']:9.4f} {atom['y']:9.4f} {atom['z']:9.4f} "
                           f"{atom['type']:<6s} {atom['subst_id']:5d} "
                           f"{subst_name}  {atom['charge']:7.4f}\n")
        
        # Write bonds
        f.write("@<TRIPOS>BOND\n")
        for bond in data['bonds']:
            f.write(f"{bond['id']:6d} {bond['atom1']:6d} {bond['atom2']:6d} {bond['type']}\n")
        
        # Write footer (secciones después de @<TRIPOS>BOND: SUBSTRUCTURE, etc.)
        for line in data.get('footer', []):
            f.write(line)


def process_mol2_file(mol2_path, min_charge=0.0001, dry_run=False, fix_aromatic=True):
    """Process a single MOL2 file to fix zero charges."""
    print(f"  Processing: {mol2_path}")
    
    # Parse MOL2 file
    data = parse_mol2_file(mol2_path)
    atoms = data['atoms']
    bonds = data['bonds']
    
    if not atoms:
        print(f"    ⚠️  No atoms found in file")
        return False, []
    
    # Check for zero charges and problematic aromatic atoms
    zero_charge_count = sum(1 for a in atoms if abs(a['charge']) < 1e-6)
    small_aromatic_count = sum(1 for a in atoms if a['type'] == 'C.ar' and 0 < abs(a['charge']) < 0.01) if fix_aromatic else 0
    
    if zero_charge_count == 0 and small_aromatic_count == 0:
        print(f"    ✓ No problematic atoms")
        return False, []
    
    if zero_charge_count > 0:
        print(f"    ⚠️  Found {zero_charge_count} atoms with charge 0.0000")
    if small_aromatic_count > 0:
        print(f"    ⚠️  Found {small_aromatic_count} C.ar atoms with charge < 0.01")
    
    # Fix zero charges and problematic aromatic atoms
    fixed_atoms, fixes = fix_zero_charges(atoms, bonds, min_charge=min_charge, fix_aromatic=fix_aromatic)
    
    if not fixes:
        return False, []
    
    # Report fixes
    total_charge_before = sum(a['charge'] for a in atoms)
    total_charge_after = sum(a['charge'] for a in fixed_atoms)
    charge_change = total_charge_after - total_charge_before
    
    print(f"    ✓ Fixed {len(fixes)} atoms")
    print(f"    ✓ Total charge: {total_charge_before:.6f} -> {total_charge_after:.6f} (Δ={charge_change:.6f})")
    
    for fix in fixes:
        print(f"      - Atom {fix['atom_id']} ({fix['atom_name']}, {fix['atom_type']}): "
              f"{fix['old_charge']:.4f} -> {fix['new_charge']:.4f} ({fix['method']})")
    
    if not dry_run:
        # Update data with fixed atoms
        data['atoms'] = fixed_atoms
        write_mol2_file(mol2_path, data, fixed_atoms)
        print(f"    ✓ File updated")
    
    return True, fixes


def find_mol2_files(base_dir, pattern="fragmento_final_charge_drift.mol2"):
    """Find all MOL2 files matching the pattern."""
    mol2_files = []
    base_path = Path(base_dir)
    
    if not base_path.exists():
        return mol2_files
    
    # Search in query_combination_* subdirectories
    for query_dir in base_path.rglob("query_combination_*"):
        mol2_file = query_dir / pattern
        if mol2_file.exists():
            mol2_files.append(mol2_file)
    
    return sorted(mol2_files)


def main():
    parser = argparse.ArgumentParser(
        description="Detect and fix atoms with 0.0000 charge in MOL2 files"
    )
    parser.add_argument(
        "prefix_type",
        choices=["GN", "LF", "all"],
        help="Combination type to process (GN, LF, or all)"
    )
    parser.add_argument(
        "--base-dir-gn",
        default="valid_GN_final",
        help="GN base directory (default: valid_GN_final)"
    )
    parser.add_argument(
        "--base-dir-lf",
        default="valid_LF_final",
        help="LF base directory (default: valid_LF_final)"
    )
    parser.add_argument(
        "--min-charge",
        type=float,
        default=0.0001,
        help="Minimum charge to assign (default: 0.0001)"
    )
    parser.add_argument(
        "--fix-aromatic",
        action="store_true",
        default=True,
        help="Also fix C.ar atoms with charge < 0.01 (default: True)"
    )
    parser.add_argument(
        "--pattern",
        default="fragmento_final_charge_drift.mol2",
        help="MOL2 filename pattern (default: fragmento_final_charge_drift.mol2)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print actions only; do not modify files"
    )
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("Fix atoms with 0.0000 charge")
    print("=" * 70)
    print(f"Type: {args.prefix_type}")
    print(f"Minimum charge: {args.min_charge}")
    print(f"Mode: {'DRY RUN' if args.dry_run else 'UPDATE FILES'}")
    print()
    
    all_fixes = []
    total_files = 0
    total_fixed = 0
    
    # Process GN
    if args.prefix_type in ["GN", "all"]:
        print(f"▶ Processing GN: {args.base_dir_gn}")
        print("-" * 70)
        gn_files = find_mol2_files(args.base_dir_gn, args.pattern)
        total_files += len(gn_files)
        
        if not gn_files:
            print(f"  ⚠️  No MOL2 files found in {args.base_dir_gn}")
        else:
            print(f"  Found {len(gn_files)} MOL2 files")
            print()
            
            for mol2_file in gn_files:
                fixed, fixes = process_mol2_file(mol2_file, args.min_charge, args.dry_run, args.fix_aromatic)
                if fixed:
                    total_fixed += 1
                    all_fixes.extend(fixes)
                print()
        
        print()
    
    # Process LF
    if args.prefix_type in ["LF", "all"]:
        print(f"▶ Processing LF: {args.base_dir_lf}")
        print("-" * 70)
        lf_files = find_mol2_files(args.base_dir_lf, args.pattern)
        total_files += len(lf_files)
        
        if not lf_files:
            print(f"  ⚠️  No MOL2 files found in {args.base_dir_lf}")
        else:
            print(f"  Found {len(lf_files)} MOL2 files")
            print()
            
            for mol2_file in lf_files:
                fixed, fixes = process_mol2_file(mol2_file, args.min_charge, args.dry_run, args.fix_aromatic)
                if fixed:
                    total_fixed += 1
                    all_fixes.extend(fixes)
                print()
        
        print()
    
    # Summary
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"Files scanned: {total_files}")
    print(f"Files modified: {total_fixed}")
    print(f"Total atom fixes: {len(all_fixes)}")
    
    if all_fixes:
        print()
        print("Fixed atoms by type:")
        atom_types = defaultdict(int)
        for fix in all_fixes:
            atom_types[fix['atom_type']] += 1
        for atom_type, count in sorted(atom_types.items()):
            print(f"  {atom_type}: {count}")
    
    print()
    if args.dry_run:
        print("⚠️  DRY RUN: no files were modified")
        print("   Run without --dry-run to apply changes")
    else:
        print("✓ Fixes applied")
        print("  You can build topologies without dummy-atom issues")


if __name__ == "__main__":
    main()

