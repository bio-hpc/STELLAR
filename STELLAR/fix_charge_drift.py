#!/usr/bin/env python3
"""
Fix Charge Drift in MOL2 File for ACPYPE Compatibility
======================================================

ACPYPE requires the net charge drift to be < 0.01. This script:
1. Reads a MOL2 file
2. Calculates the net charge
3. Adjusts charges to match target charge (default: 2.0)
4. Writes corrected MOL2 file

Usage:
    python fix_charge_drift.py input.mol2 -o output.mol2 --target 2.0
    python fix_charge_drift.py input.mol2 -o output.mol2 --target 2.0 --tolerance 0.009
"""

import argparse
import sys
import re

def parse_mol2_file(filename):
    """Parse MOL2 file and extract atom information."""
    atoms = []
    bonds = []
    header_lines = []
    in_atom_section = False
    in_bond_section = False
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                in_bond_section = False
                header_lines.append(line)
                continue
            elif line.startswith("@<TRIPOS>BOND"):
                in_atom_section = False
                in_bond_section = True
                header_lines.append(line)
                continue
            elif line.startswith("@<TRIPOS>"):
                in_atom_section = False
                in_bond_section = False
                header_lines.append(line)
                continue
            
            if in_atom_section and line.strip():
                parts = line.split()
                if len(parts) >= 8:
                    atom_id = int(parts[0])
                    atom_name = parts[1]
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                    atom_type = parts[5]
                    subst_id = int(parts[6])
                    subst_name = parts[7]
                    charge = float(parts[8]) if len(parts) > 8 else 0.0
                    
                    atoms.append({
                        'id': atom_id,
                        'name': atom_name,
                        'x': x,
                        'y': y,
                        'z': z,
                        'type': atom_type,
                        'subst_id': subst_id,
                        'subst_name': subst_name,
                        'charge': charge,
                        'original_line': line
                    })
            elif in_bond_section and line.strip():
                parts = line.split()
                if len(parts) >= 4:
                    bonds.append(line)
            elif not in_atom_section and not in_bond_section:
                header_lines.append(line)
    
    return atoms, bonds, header_lines

def calculate_net_charge(atoms):
    """Calculate net charge from atom charges."""
    return sum(atom['charge'] for atom in atoms)

def adjust_charges(atoms, target_charge, tolerance=0.009):
    """
    Adjust charges to match target charge within tolerance.
    Uses weighted distribution to preserve charge distribution.
    """
    current_net = calculate_net_charge(atoms)
    drift = current_net - target_charge
    
    print(f"[Info] Current net charge: {current_net:.6f}")
    print(f"[Info] Target charge: {target_charge:.2f}")
    print(f"[Info] Charge drift: {abs(drift):.6f}")
    print(f"[Info] Tolerance: {tolerance:.6f}")
    
    if abs(drift) <= tolerance:
        print(f"[Info] ✓ Charge drift already within tolerance!")
        return atoms
    
    print(f"[Info] Adjusting charges to fix drift...")
    
    # Calculate adjustment needed
    adjustment = target_charge - current_net
    
    # Calculate weights based on absolute charge (preserve charge distribution)
    abs_charges = [abs(atom['charge']) for atom in atoms]
    total_abs = sum(abs_charges) if sum(abs_charges) > 0.001 else len(atoms)
    
    # First pass: distribute adjustment weighted by absolute charge
    adjusted_atoms = []
    for atom in atoms:
        if total_abs > 0.001:
            weight = abs_charges[atoms.index(atom)] / total_abs
        else:
            weight = 1.0 / len(atoms)
        
        new_charge = atom['charge'] + adjustment * weight
        atom_copy = atom.copy()
        atom_copy['charge'] = new_charge
        adjusted_atoms.append(atom_copy)
    
    # Verify first pass
    new_net = calculate_net_charge(adjusted_atoms)
    new_drift = abs(new_net - target_charge)
    
    print(f"[Info] After first pass: net={new_net:.6f}, drift={new_drift:.6f}")
    
    if new_drift <= tolerance:
        print(f"[Info] ✓ Charge drift fixed!")
        return adjusted_atoms
    
    # Second pass: fine-tune with very small adjustments
    print(f"[Info] Fine-tuning charges...")
    remaining = target_charge - new_net
    
    # Sort by absolute charge (highest first)
    sorted_indices = sorted(range(len(adjusted_atoms)), 
                           key=lambda i: abs(adjusted_atoms[i]['charge']), 
                           reverse=True)
    
    # Apply small adjustments to top atoms
    adjustment_per_atom = remaining / min(20, len(adjusted_atoms))
    
    for idx in sorted_indices[:min(20, len(sorted_indices))]:
        if abs(remaining) < 0.0001:
            break
        
        current = adjusted_atoms[idx]['charge']
        # Make tiny adjustment proportional to remaining
        adjustment_small = remaining * 0.05  # Very conservative
        new_val = current + adjustment_small
        adjusted_atoms[idx]['charge'] = new_val
        remaining -= adjustment_small
    
    # Final verification
    final_net = calculate_net_charge(adjusted_atoms)
    final_drift = abs(final_net - target_charge)
    
    print(f"[Info] Final net charge: {final_net:.6f}")
    print(f"[Info] Final charge drift: {final_drift:.6f}")
    
    if final_drift <= tolerance:
        print(f"[Info] ✓ Charge drift fixed after fine-tuning!")
        return adjusted_atoms
    else:
        print(f"[WARN] Charge drift still high: {final_drift:.6f} > {tolerance:.6f}")
        print(f"[WARN] Trying more aggressive adjustment...")
        
        # More aggressive: adjust all atoms slightly
        final_remaining = target_charge - final_net
        per_atom = final_remaining / len(adjusted_atoms)
        
        for atom in adjusted_atoms:
            atom['charge'] += per_atom
        
        very_final_net = calculate_net_charge(adjusted_atoms)
        very_final_drift = abs(very_final_net - target_charge)
        
        print(f"[Info] Very final net charge: {very_final_net:.6f}")
        print(f"[Info] Very final charge drift: {very_final_drift:.6f}")
        
        if very_final_drift <= tolerance:
            print(f"[Info] ✓ Charge drift fixed with aggressive adjustment!")
        else:
            print(f"[ERROR] Could not reduce drift below tolerance!")
            print(f"[ERROR] You may need to check the input charges manually.")
        
        return adjusted_atoms

def write_mol2_file(filename, atoms, bonds, header_lines):
    """Write corrected MOL2 file."""
    with open(filename, 'w') as f:
        # Write header (everything before @<TRIPOS>ATOM)
        for line in header_lines:
            if not line.startswith("@<TRIPOS>ATOM"):
                f.write(line)
            else:
                break
        
        # Write atoms
        f.write("@<TRIPOS>ATOM\n")
        for atom in atoms:
            # Format: atom_id atom_name x y z atom_type subst_id subst_name charge
            # %7d %-8s %9.4f %9.4f %9.4f %-6s %5d %-4s %7.4f
            f.write(f"{atom['id']:7d} {atom['name']:8s} "
                   f"{atom['x']:9.4f} {atom['y']:9.4f} {atom['z']:9.4f} "
                   f"{atom['type']:6s} {atom['subst_id']:5d} {atom['subst_name']:4s} "
                   f"{atom['charge']:7.4f}\n")
        
        # Write bonds
        f.write("@<TRIPOS>BOND\n")
        for bond in bonds:
            f.write(bond)
        
        # Write any remaining header lines after bonds
        in_bond_section = False
        for line in header_lines:
            if line.startswith("@<TRIPOS>BOND"):
                in_bond_section = True
            elif line.startswith("@<TRIPOS>") and in_bond_section:
                f.write(line)
                in_bond_section = False

def main():
    ap = argparse.ArgumentParser(
        description="Fix charge drift in MOL2 file for ACPYPE compatibility",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Fix charge drift to target 2.0
  python fix_charge_drift.py merged.mol2 -o merged_fixed.mol2 --target 2.0
  
  # Use stricter tolerance
  python fix_charge_drift.py merged.mol2 -o merged_fixed.mol2 --target 2.0 --tolerance 0.009
        """
    )
    ap.add_argument("input", help="Input MOL2 file")
    ap.add_argument("-o", "--output", required=True, help="Output MOL2 file")
    ap.add_argument("--target", type=float, default=2.0, help="Target net charge (default: 2.0)")
    ap.add_argument("--tolerance", type=float, default=0.009, help="Maximum allowed drift (default: 0.009, stricter than ACPYPE's 0.01)")
    
    args = ap.parse_args()
    
    print(f"[Info] Reading MOL2 file: {args.input}")
    try:
        atoms, bonds, header_lines = parse_mol2_file(args.input)
        print(f"[Info] Found {len(atoms)} atoms and {len(bonds)} bonds")
    except Exception as e:
        print(f"[ERROR] Failed to read MOL2 file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Adjust charges
    adjusted_atoms = adjust_charges(atoms, args.target, args.tolerance)
    
    # Write corrected file
    print(f"\n[Info] Writing corrected MOL2 file: {args.output}")
    try:
        write_mol2_file(args.output, adjusted_atoms, bonds, header_lines)
        print(f"[OK] Successfully wrote {args.output}")
        
        # Verify final charge
        final_net = calculate_net_charge(adjusted_atoms)
        final_drift = abs(final_net - args.target)
        print(f"\n[Info] Verification:")
        print(f"  Final net charge: {final_net:.6f}")
        print(f"  Charge drift: {final_drift:.6f}")
        print(f"  Tolerance: {args.tolerance:.6f}")
        
        if final_drift <= args.tolerance:
            print(f"[OK] ✓ File is ACPYPE-compatible!")
        else:
            print(f"[WARN] ⚠ Charge drift still exceeds tolerance")
            print(f"[WARN] ACPYPE may still fail. Consider using --tolerance 0.008")
    except Exception as e:
        print(f"[ERROR] Failed to write MOL2 file: {e}", file=sys.stderr)
        sys.exit(2)

if __name__ == "__main__":
    main()

