#!/usr/bin/env python3
# Merge & Light-Minimize MOL2 Fragments (RDKit)
# ============================================
#
# Quickly merge MOL2 fragments (3, 4, 5, or more), (optionally) add inter-fragment bonds, and do a gentle geometry
# minimization to relieve small steric clashes—without drastically changing the geometry.
#
# Requirements
# ------------
# - Python 3.8+
# - RDKit (2022.09+ recommended)
#
# Install RDKit (conda):
#     conda install -c conda-forge rdkit
#
# Usage
# -----
# Basic (no new bonds, just merge & relax):
#     python relax_merge_mol2.py --folder combo_dir -o merged_relaxed.mol2
#     # O especificar fragmentos manualmente:
#     python relax_merge_mol2.py frag1.mol2 frag2.mol2 frag3.mol2 -o merged_relaxed.mol2
#
# Specify bonds to connect fragments (atom indices are 0-based in each original fragment,
# but the script will map them to the combined indexing automatically when given via flags):
#
#     # Ejemplo: conectar atom 12 en frag1 a atom 7 en frag2, y atom 3 en frag3 a atom 25 en frag1:
#     python relax_merge_mol2.py frag1.mol2 frag2.mol2 frag3.mol2 \
#         -o merged_relaxed.mol2 \
#         --bond f1:12-f2:7 --bond f3:3-f1:25
#
# Control minimization strength (defaults are gentle):
#     python relax_merge_mol2.py --folder combo_dir -o out.mol2 --iters 100 --ff UFF
#
# Notes
# -----
# - The script prints a brief clash report (pairs closer than scaled VDW sum) before and after minimization.
# - If MMFF fails (e.g., missing parameters), it falls back to UFF automatically.
# - If your fragments overlap heavily, consider pre-aligning in PyMOL/ChimeraX or giving a few bonds so the force field "knows" connectivity.
#
import argparse
import sys
import re
import os
import glob
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
import math

VDW = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
    'P': 1.80, 'S': 1.80, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98
}
DEFAULT_SCALE = 0.85  # how strict clash threshold is vs VDW sum

def parse_mol2_atoms(path):
    """Parse MOL2 file to extract atom info (name, type, residue) for peptide terminus detection."""
    atoms = []
    in_atom_section = False
    with open(path, 'r') as f:
        for line in f:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                continue
            elif line.startswith("@<TRIPOS>"):
                in_atom_section = False
                continue
            
            if in_atom_section and line.strip():
                parts = line.split()
                if len(parts) >= 8:
                    atom_id = int(parts[0]) - 1  # Convert to 0-based
                    atom_name = parts[1]
                    atom_type = parts[5]
                    residue_num = int(parts[6])
                    residue_name = parts[7]
                    atoms.append({
                        'id': atom_id,
                        'name': atom_name,
                        'type': atom_type,
                        'residue_num': residue_num,
                        'residue_name': residue_name
                    })
    return atoms

def find_cterminus(mol2_atoms):
    """Find C-terminus atom (C with type C.2) in the last residue."""
    if not mol2_atoms:
        return None
    
    # Find the last residue number
    last_residue = max(atom['residue_num'] for atom in mol2_atoms)
    
    # Find C atom with type C.2 in the last residue (C-terminus carbon)
    for atom in reversed(mol2_atoms):
        if (atom['residue_num'] == last_residue and 
            atom['name'] == 'C' and 
            atom['type'].startswith('C.2')):
            return atom['id']
    
    # Fallback: any C in last residue
    for atom in reversed(mol2_atoms):
        if atom['residue_num'] == last_residue and atom['name'] == 'C':
            return atom['id']
    
    return None

def find_nterminus(mol2_atoms):
    """Find N-terminus atom (N) in the first residue."""
    if not mol2_atoms:
        return None
    
    # Find the first residue number
    first_residue = min(atom['residue_num'] for atom in mol2_atoms)
    
    # Find N atom in the first residue
    for atom in mol2_atoms:
        if (atom['residue_num'] == first_residue and 
            atom['name'] == 'N' and 
            atom['type'].startswith('N')):
            return atom['id']
    
    return None

def extract_fragment_number(filename):
    """Extract fragment number from filename (e.g., Frag1.mol2 -> 1)."""
    basename = os.path.basename(filename)
    # Try different patterns
    patterns = [
        r'[Ff]rag\s*(\d+)',  # Frag1, Frag2, etc.
        r'fragment\s*(\d+)',  # fragment1, fragment2, etc.
        r'_(\d+)\.mol2',     # _1.mol2, _2.mol2, etc.
        r'(\d+)\.mol2'       # 1.mol2, 2.mol2, etc.
    ]
    for pattern in patterns:
        match = re.search(pattern, basename, re.IGNORECASE)
        if match:
            return int(match.group(1))
    return None

def load_mol2_with_names(path):
    """Load MOL2 and extract atom names, residue info, and charges to store as properties."""
    mol = Chem.MolFromMol2File(path, sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError(f"Failed to read MOL2: {path}")
    if mol.GetNumConformers() == 0:
        raise ValueError(f"No 3D coordinates found in: {path}")
    
    # Parse atom names, residue info, and charges from MOL2 file
    in_atom_section = False
    atom_idx = 0
    with open(path, 'r') as f:
        for line in f:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                continue
            elif line.startswith("@<TRIPOS>"):
                in_atom_section = False
                continue
            
            if in_atom_section and line.strip():
                parts = line.split()
                if len(parts) >= 2 and atom_idx < mol.GetNumAtoms():
                    atom_name = parts[1].strip()
                    atom = mol.GetAtomWithIdx(atom_idx)
                    atom.SetProp('_Name', atom_name)
                    
                    # Store residue information if available
                    if len(parts) >= 8:
                        residue_num = parts[6]
                        residue_name = parts[7]
                        atom.SetProp('_ResidueNum', residue_num)
                        atom.SetProp('_ResidueName', residue_name)
                        
                        # Store charge if available
                        if len(parts) >= 9:
                            try:
                                charge = float(parts[8])
                                atom.SetProp('_Charge', f"{charge:.4f}")
                            except:
                                atom.SetProp('_Charge', "0.0000")
                        else:
                            atom.SetProp('_Charge', "0.0000")
                    
                    atom_idx += 1
    
    return mol

def load_mol2(path):
    return load_mol2_with_names(path)

def combine_with_coords(mols):
    """Combine molecules preserving 3D coordinates and residue information; returns (RWMol, atom_offsets, residue_offsets).
    Supports any number of fragments (not only three)."""
    # Combine all fragments
    if len(mols) == 0:
        raise ValueError("No molecules to combine")
    elif len(mols) == 1:
        combo = mols[0]
    else:
        # Combinar de dos en dos: primero combinar los últimos dos, luego con el anterior, etc.
        combo = mols[0]
        for i in range(1, len(mols)):
            combo = Chem.CombineMols(combo, mols[i])
    rw = Chem.RWMol(combo)
    conf = Chem.Conformer(rw.GetNumAtoms())
    offsets = []
    residue_offsets = []  # Track residue number offsets for each fragment
    idx = 0
    max_residue = 0
    
    for m in mols:
        offsets.append(idx)
        
        # Find max residue number in this fragment
        frag_max_res = 0
        for a in range(m.GetNumAtoms()):
            atom = m.GetAtomWithIdx(a)
            if atom.HasProp('_ResidueNum'):
                try:
                    res_num = int(atom.GetProp('_ResidueNum'))
                    frag_max_res = max(frag_max_res, res_num)
                except:
                    pass
        
        residue_offsets.append(max_residue)
        max_residue += frag_max_res
        
        # Copy coordinates
        c = m.GetConformer()
        for a in range(m.GetNumAtoms()):
            p = c.GetAtomPosition(a)
            conf.SetAtomPosition(idx, Point3D(p.x, p.y, p.z))
            
            # Copy residue properties to combined molecule
            atom_orig = m.GetAtomWithIdx(a)
            atom_new = rw.GetAtomWithIdx(idx)
            
            if atom_orig.HasProp('_Name'):
                atom_new.SetProp('_Name', atom_orig.GetProp('_Name'))
            if atom_orig.HasProp('_ResidueName'):
                atom_new.SetProp('_ResidueName', atom_orig.GetProp('_ResidueName'))
            if atom_orig.HasProp('_ResidueNum'):
                # Update residue number to be sequential across fragments
                orig_res_num = int(atom_orig.GetProp('_ResidueNum'))
                new_res_num = max_residue - frag_max_res + orig_res_num
                atom_new.SetProp('_ResidueNum', str(new_res_num))
            if atom_orig.HasProp('_Charge'):
                atom_new.SetProp('_Charge', atom_orig.GetProp('_Charge'))
            
            idx += 1
    
    rw.AddConformer(conf, assignId=True)
    return rw, offsets, residue_offsets

def guess_vdw(symbol):
    return VDW.get(symbol, 1.75)

def distance(conf, i, j):
    pi = conf.GetAtomPosition(i)
    pj = conf.GetAtomPosition(j)
    dx, dy, dz = (pi.x - pj.x), (pi.y - pj.y), (pi.z - pj.z)
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def clash_pairs(mol, scale=DEFAULT_SCALE):
    conf = mol.GetConformer()
    pairs = []
    nat = mol.GetNumAtoms()
    for i in range(nat):
        ai = mol.GetAtomWithIdx(i)
        si = ai.GetSymbol()
        ri = guess_vdw(si)
        for j in range(i+1, nat):
            aj = mol.GetAtomWithIdx(j)
            if mol.GetBondBetweenAtoms(i, j):
                continue
            sj = aj.GetSymbol()
            rj = guess_vdw(sj)
            cutoff = scale*(ri + rj)
            d = distance(conf, i, j)
            if d < cutoff:
                pairs.append((i, j, d, cutoff))
    return pairs

def check_atoms_too_close(mol, min_distance=0.5):
    """
    Check for atoms that are too close together (< min_distance, default 0.5 Å).
    ACPYPE will fail on these. Returns list of (atom_idx1, atom_idx2, distance, is_bonded).
    """
    conf = mol.GetConformer()
    too_close = []
    nat = mol.GetNumAtoms()
    
    for i in range(nat):
        for j in range(i+1, nat):
            d = distance(conf, i, j)
            if d < min_distance:
                is_bonded = mol.GetBondBetweenAtoms(i, j) is not None
                too_close.append((i, j, d, is_bonded))
    
    return too_close

def fix_atoms_too_close(mol, min_distance=0.5, step_size=0.4, max_iters=150):
    """
    Fix atoms that are too close together by moving them apart.
    For bonded atoms: moves them to reasonable bond lengths (~1.0-1.5 Å).
    For non-bonded: moves them apart more aggressively.
    """
    conf = mol.GetConformer()
    fixed_count = 0
    
    for iteration in range(max_iters):
        too_close = check_atoms_too_close(mol, min_distance)
        
        if len(too_close) == 0:
            if iteration > 0:
                print(f"[Info] Fixed {fixed_count} atom pair(s) that were too close in {iteration} iterations")
            break
        
        if iteration == 0:
            print(f"[Info] Found {len(too_close)} atom pair(s) too close (< {min_distance} Å), fixing...")
            # Show first few
            for i, j, d, is_bonded in too_close[:10]:
                atom_i = mol.GetAtomWithIdx(i)
                atom_j = mol.GetAtomWithIdx(j)
                bond_type = "BONDED" if is_bonded else "non-bonded"
                print(f"  Atoms {i+1} ({atom_i.GetSymbol()}) - {j+1} ({atom_j.GetSymbol()}): {d:.4f} Å ({bond_type})")
        
        # Fix each pair
        for atom_idx1, atom_idx2, dist, is_bonded in too_close:
            pi = conf.GetAtomPosition(atom_idx1)
            pj = conf.GetAtomPosition(atom_idx2)
            
            # Vector from atom1 to atom2
            dx = pj.x - pi.x
            dy = pj.y - pi.y
            dz = pj.z - pi.z
            
            # Normalize
            norm = math.sqrt(dx*dx + dy*dy + dz*dz)
            if norm < 0.001:
                # Atoms are at same position - separate in a direction
                dx, dy, dz = 1.0, 0.0, 0.0
                norm = 1.0
            
            dx /= norm
            dy /= norm
            dz /= norm
            
            # Target distance based on atom types
            atom1 = mol.GetAtomWithIdx(atom_idx1)
            atom2 = mol.GetAtomWithIdx(atom_idx2)
            symbols = sorted([atom1.GetSymbol(), atom2.GetSymbol()])
            
            if is_bonded:
                # For bonded atoms, use appropriate bond lengths
                if 'H' in symbols:
                    target_dist = 1.1  # H-X bonds
                elif 'C' in symbols and 'N' in symbols:
                    target_dist = 1.4  # C-N bond
                elif 'C' in symbols and 'O' in symbols:
                    target_dist = 1.4  # C-O bond
                elif symbols == ['C', 'C']:
                    target_dist = 1.5  # C-C bond
                else:
                    target_dist = 1.5  # Default bond length
            else:
                # For non-bonded, move apart more aggressively
                target_dist = min_distance * 3.0  # At least 1.5 Å for non-bonded
            
            # How much to separate
            move_amount = (target_dist - dist) * step_size
            
            # Weight by atomic number (heavier atoms move less)
            wi = 1.0 / (1.0 + atom1.GetAtomicNum())
            wj = 1.0 / (1.0 + atom2.GetAtomicNum())
            total_w = wi + wj
            
            move_i = move_amount * (wj / total_w)
            move_j = move_amount * (wi / total_w)
            
            # Limit movement per iteration
            max_move = 0.5
            if move_i > max_move:
                move_i = max_move
            if move_j > max_move:
                move_j = max_move
            
            # Move atoms apart (in opposite directions along the vector)
            new_pos_i = Point3D(pi.x - dx * move_i, pi.y - dy * move_i, pi.z - dz * move_i)
            new_pos_j = Point3D(pj.x + dx * move_j, pj.y + dy * move_j, pj.z + dz * move_j)
            
            conf.SetAtomPosition(atom_idx1, new_pos_i)
            conf.SetAtomPosition(atom_idx2, new_pos_j)
            fixed_count += 2
        
        # Check progress
        too_close_after = check_atoms_too_close(mol, min_distance)
        if len(too_close_after) == 0:
            break
        elif len(too_close_after) >= len(too_close) and iteration > 20:
            if iteration > 0:
                print(f"[Info] Fixed {fixed_count} atom pair(s), {len(too_close_after)} remain")
            break
    
    return fixed_count

def check_atoms_near_newly_formed_bonds(mol, new_bonds, threshold=1.0):
    """
    Check for atoms that are too close (< threshold Å) to newly formed bonds.
    This ensures proper geometry around inter-fragment bonds.
    INCLUDES HYDROGENS: All atoms including hydrogens are checked, even if they're
    bonded to the bond atoms (hydrogens can still be too close to the bond vector).
    Returns list of (atom_idx, bond_atom1, bond_atom2, distance_to_bond).
    """
    conf = mol.GetConformer()
    problematic = []
    
    for bond_atom1, bond_atom2 in new_bonds:
        p1 = conf.GetAtomPosition(bond_atom1)
        p2 = conf.GetAtomPosition(bond_atom2)
        
        # Bond vector
        bond_vec = Point3D(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z)
        bond_length = math.sqrt(bond_vec.x**2 + bond_vec.y**2 + bond_vec.z**2)
        
        if bond_length < 0.001:
            continue
        
        # Normalize bond vector
        bond_vec.x /= bond_length
        bond_vec.y /= bond_length
        bond_vec.z /= bond_length
        
        # Check all other atoms (including hydrogens)
        for atom_idx in range(mol.GetNumAtoms()):
            if atom_idx == bond_atom1 or atom_idx == bond_atom2:
                continue
            
            atom = mol.GetAtomWithIdx(atom_idx)
            is_hydrogen = atom.GetSymbol() == 'H'
            
            # For hydrogens: always check them (they can be too close even if bonded)
            # For other atoms: skip if directly bonded to either bond atom (expected to be close)
            if not is_hydrogen:
                if mol.GetBondBetweenAtoms(atom_idx, bond_atom1) or mol.GetBondBetweenAtoms(atom_idx, bond_atom2):
                    continue
            
            p_atom = conf.GetAtomPosition(atom_idx)
            
            # Vector from bond_atom1 to the atom
            vec_to_atom = Point3D(p_atom.x - p1.x, p_atom.y - p1.y, p_atom.z - p1.z)
            
            # Project onto bond vector to find closest point on bond
            proj = vec_to_atom.x * bond_vec.x + vec_to_atom.y * bond_vec.y + vec_to_atom.z * bond_vec.z
            
            # Clamp to bond segment
            proj = max(0.0, min(bond_length, proj))
            
            # Closest point on bond
            closest_on_bond = Point3D(
                p1.x + proj * bond_vec.x,
                p1.y + proj * bond_vec.y,
                p1.z + proj * bond_vec.z
            )
            
            # Distance from atom to closest point on bond
            dx = p_atom.x - closest_on_bond.x
            dy = p_atom.y - closest_on_bond.y
            dz = p_atom.z - closest_on_bond.z
            dist_to_bond = math.sqrt(dx*dx + dy*dy + dz*dz)
            
            if dist_to_bond < threshold:
                problematic.append((atom_idx, bond_atom1, bond_atom2, dist_to_bond))
    
    return problematic

def fix_atoms_near_newly_formed_bonds(mol, new_bonds, threshold=1.0, step_size=0.3, max_iters=100):
    """
    Fix atoms that are too close to newly formed bonds by moving them away.
    """
    conf = mol.GetConformer()
    fixed_count = 0
    
    for iteration in range(max_iters):
        problematic = check_atoms_near_newly_formed_bonds(mol, new_bonds, threshold)
        
        if len(problematic) == 0:
            if iteration > 0:
                print(f"[Info] Fixed {fixed_count} atom(s) too close to newly formed bonds in {iteration} iterations")
            break
        
        if iteration == 0:
            print(f"[Info] Found {len(problematic)} atom(s) within {threshold} Å of newly formed bonds, fixing...")
            for atom_idx, bond_atom1, bond_atom2, dist in problematic[:10]:
                atom = mol.GetAtomWithIdx(atom_idx)
                atom_name = atom.GetProp('_Name') if atom.HasProp('_Name') else f"{atom.GetSymbol()}{atom_idx+1}"
                print(f"  Atom {atom_idx+1} ({atom_name}): {dist:.4f} Å from bond {bond_atom1+1}-{bond_atom2+1}")
        
        # Move each problematic atom away from the bond
        for atom_idx, bond_atom1, bond_atom2, dist_to_bond in problematic:
            p_atom = conf.GetAtomPosition(atom_idx)
            p1 = conf.GetAtomPosition(bond_atom1)
            p2 = conf.GetAtomPosition(bond_atom2)
            
            atom = mol.GetAtomWithIdx(atom_idx)
            is_hydrogen = atom.GetSymbol() == 'H'
            
            # Find closest point on bond
            bond_vec = Point3D(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z)
            bond_length = math.sqrt(bond_vec.x**2 + bond_vec.y**2 + bond_vec.z**2)
            
            if bond_length < 0.001:
                continue
            
            bond_vec.x /= bond_length
            bond_vec.y /= bond_length
            bond_vec.z /= bond_length
            
            vec_to_atom = Point3D(p_atom.x - p1.x, p_atom.y - p1.y, p_atom.z - p1.z)
            proj = vec_to_atom.x * bond_vec.x + vec_to_atom.y * bond_vec.y + vec_to_atom.z * bond_vec.z
            proj = max(0.0, min(bond_length, proj))
            
            closest_on_bond = Point3D(
                p1.x + proj * bond_vec.x,
                p1.y + proj * bond_vec.y,
                p1.z + proj * bond_vec.z
            )
            
            # Vector from closest point on bond to atom (direction to move away)
            dx = p_atom.x - closest_on_bond.x
            dy = p_atom.y - closest_on_bond.y
            dz = p_atom.z - closest_on_bond.z
            
            norm = math.sqrt(dx*dx + dy*dy + dz*dz)
            if norm < 0.001:
                # Atom is on the bond line - move perpendicular
                # Use a perpendicular vector
                if abs(bond_vec.x) < 0.9:
                    dx, dy, dz = 1.0, 0.0, 0.0
                else:
                    dx, dy, dz = 0.0, 1.0, 0.0
                norm = 1.0
            
            dx /= norm
            dy /= norm
            dz /= norm
            
            # Move amount - adjust for atom type
            target_dist = threshold * 1.2  # Move to 1.2x threshold
            move_amount = (target_dist - dist_to_bond) * step_size
            
            # Hydrogens can move more easily (they're lighter)
            if is_hydrogen:
                move_amount *= 2.0  # Allow hydrogens to move 50% more
                max_move = 0.7  # Higher max movement for hydrogens
            else:
                max_move = 0.5  # Standard max movement for heavy atoms
            
            if move_amount > max_move:
                move_amount = max_move
            if move_amount < 0.05:
                move_amount = 0.05
            
            # Move atom away from bond
            new_pos = Point3D(
                p_atom.x + dx * move_amount,
                p_atom.y + dy * move_amount,
                p_atom.z + dz * move_amount
            )
            conf.SetAtomPosition(atom_idx, new_pos)
            fixed_count += 1
        
        # Check progress
        problematic_after = check_atoms_near_newly_formed_bonds(mol, new_bonds, threshold)
        if len(problematic_after) == 0:
            break
        elif len(problematic_after) >= len(problematic) and iteration > 20:
            break
    
    return fixed_count

def check_newly_formed_bond_lengths(mol, new_bonds, max_length=3.5):
    """
    Check bond lengths of newly formed bonds.
    Returns list of (bond_atom1, bond_atom2, bond_length) for bonds > max_length.
    """
    conf = mol.GetConformer()
    too_long = []
    
    for bond_atom1, bond_atom2 in new_bonds:
        p1 = conf.GetAtomPosition(bond_atom1)
        p2 = conf.GetAtomPosition(bond_atom2)
        
        dx = p2.x - p1.x
        dy = p2.y - p1.y
        dz = p2.z - p1.z
        bond_length = math.sqrt(dx*dx + dy*dy + dz*dz)
        
        if bond_length > max_length:
            too_long.append((bond_atom1, bond_atom2, bond_length))
    
    return too_long

def fix_newly_formed_bond_lengths(mol, new_bonds, max_length=3.5, target_length=1.4, step_size=0.3, max_iters=100):
    """
    Fix newly formed bonds that are too long by moving atoms closer together.
    For peptide bonds, target_length is typically 1.4 Å (C-N bond).
    """
    conf = mol.GetConformer()
    fixed_count = 0
    
    for iteration in range(max_iters):
        too_long = check_newly_formed_bond_lengths(mol, new_bonds, max_length)
        
        if len(too_long) == 0:
            if iteration > 0:
                print(f"[Info] Fixed {fixed_count} newly formed bond(s) that were too long in {iteration} iterations")
            break
        
        if iteration == 0:
            print(f"[Info] Found {len(too_long)} newly formed bond(s) longer than {max_length} Å, fixing...")
            for bond_atom1, bond_atom2, bond_length in too_long:
                atom1 = mol.GetAtomWithIdx(bond_atom1)
                atom2 = mol.GetAtomWithIdx(bond_atom2)
                atom1_name = atom1.GetProp('_Name') if atom1.HasProp('_Name') else f"{atom1.GetSymbol()}{bond_atom1+1}"
                atom2_name = atom2.GetProp('_Name') if atom2.HasProp('_Name') else f"{atom2.GetSymbol()}{bond_atom2+1}"
                print(f"  Bond {bond_atom1+1} ({atom1_name}) - {bond_atom2+1} ({atom2_name}): {bond_length:.4f} Å")
        
        # Fix each bond that's too long
        for bond_atom1, bond_atom2, bond_length in too_long:
            p1 = conf.GetAtomPosition(bond_atom1)
            p2 = conf.GetAtomPosition(bond_atom2)
            
            # Vector from atom1 to atom2
            dx = p2.x - p1.x
            dy = p2.y - p1.y
            dz = p2.z - p1.z
            
            # Normalize
            norm = math.sqrt(dx*dx + dy*dy + dz*dz)
            if norm < 0.001:
                continue
            
            dx /= norm
            dy /= norm
            dz /= norm
            
            # How much to move atoms closer together
            # We want to reduce the bond length to target_length
            move_amount = (bond_length - target_length) * step_size
            
            # Weight by atomic number (heavier atoms move less)
            atom1 = mol.GetAtomWithIdx(bond_atom1)
            atom2 = mol.GetAtomWithIdx(bond_atom2)
            wi = 1.0 / (1.0 + atom1.GetAtomicNum())
            wj = 1.0 / (1.0 + atom2.GetAtomicNum())
            total_w = wi + wj
            
            move_i = move_amount * (wj / total_w)
            move_j = move_amount * (wi / total_w)
            
            # Limit movement per iteration
            max_move = 0.8
            if move_i > max_move:
                move_i = max_move
            if move_j > max_move:
                move_j = max_move
            
            # Move atoms closer together
            # Atom1 moves towards atom2, atom2 moves towards atom1
            new_pos1 = Point3D(
                p1.x + dx * move_i,
                p1.y + dy * move_i,
                p1.z + dz * move_i
            )
            new_pos2 = Point3D(
                p2.x - dx * move_j,
                p2.y - dy * move_j,
                p2.z - dz * move_j
            )
            
            conf.SetAtomPosition(bond_atom1, new_pos1)
            conf.SetAtomPosition(bond_atom2, new_pos2)
            fixed_count += 2
        
        # Check progress
        too_long_after = check_newly_formed_bond_lengths(mol, new_bonds, max_length)
        if len(too_long_after) == 0:
            break
        elif len(too_long_after) >= len(too_long) and iteration > 20:
            break
    
    return fixed_count

def write_mol2_with_residues(mol, filename):
    """
    Write MOL2 file manually preserving residue information.
    This ensures proper residue names and numbers are maintained.
    """
    conf = mol.GetConformer()
    nat = mol.GetNumAtoms()
    nbonds = mol.GetNumBonds()
    
    with open(filename, 'w') as f:
        # Header
        f.write("@<TRIPOS>MOLECULE\n")
        f.write("merged_molecule\n")
        f.write(f"  {nat} {nbonds} 0 0 0\n")
        f.write("SMALL\n")
        f.write("USER_CHARGES\n\n")
        
        # Atoms
        f.write("@<TRIPOS>ATOM\n")
        for i in range(nat):
            atom = mol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            
            # Get atom name
            atom_name = atom.GetProp('_Name') if atom.HasProp('_Name') else f"{atom.GetSymbol()}{i+1}"
            
            # Get residue information
            if atom.HasProp('_ResidueNum') and atom.HasProp('_ResidueName'):
                residue_num = int(atom.GetProp('_ResidueNum'))
                residue_name = atom.GetProp('_ResidueName')
            else:
                # Default if residue info not available
                residue_num = 1
                residue_name = "UNK"
            
            # Get charge
            if atom.HasProp('_Charge'):
                try:
                    charge = float(atom.GetProp('_Charge'))
                except:
                    charge = 0.0
            else:
                charge = 0.0
            
            # Map RDKit atom types to MOL2 atom types
            symbol = atom.GetSymbol()
            degree = atom.GetDegree()
            
            if symbol == 'C':
                if degree == 4:
                    mol2_type = "C.3"
                elif degree == 3:
                    mol2_type = "C.2"
                elif degree == 2:
                    mol2_type = "C.1"
                else:
                    mol2_type = "C.ar"
            elif symbol == 'N':
                if degree == 3:
                    mol2_type = "N.3"
                elif degree == 2:
                    mol2_type = "N.2"
                elif degree == 1:
                    mol2_type = "N.1"
                else:
                    mol2_type = "N.am"
            elif symbol == 'O':
                if degree == 2:
                    mol2_type = "O.2"
                elif degree == 1:
                    mol2_type = "O.3"
                else:
                    mol2_type = "O.co2"
            elif symbol == 'S':
                mol2_type = "S.3"
            elif symbol == 'P':
                mol2_type = "P.3"
            elif symbol == 'H':
                mol2_type = "H"
            else:
                mol2_type = symbol
            
            # Write atom line: atom_id atom_name x y z atom_type subst_id subst_name charge
            f.write(f"{i+1:7d} {atom_name:8s} {pos.x:9.4f} {pos.y:9.4f} {pos.z:9.4f} {mol2_type:6s} {residue_num:5d} {residue_name:4s} {charge:7.4f}\n")
        
        # Bonds
        f.write("@<TRIPOS>BOND\n")
        bond_idx = 1
        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx() + 1  # MOL2 is 1-based
            j = bond.GetEndAtomIdx() + 1
            order = bond.GetBondTypeAsDouble()
            if order == 1.0:
                bond_type = "1"
            elif order == 2.0:
                bond_type = "2"
            elif order == 3.0:
                bond_type = "3"
            else:
                bond_type = "am"  # amide
            
            f.write(f"{bond_idx:6d} {i:6d} {j:6d} {bond_type}\n")
            bond_idx += 1

def add_bond_from_flag(rw, flag, offsets):
    # flag like "f1:12-f2:7" or "f2:5-f3:10" (any fragment count)
    left, right = flag.split('-')
    fA, aA = left.split(':')
    fB, aB = right.split(':')
    # Extraer número de fragmento: f1 -> 1, f2 -> 2, etc. (convertir a índice 0-based)
    i_frag = int(fA.lower().replace('f', '')) - 1
    j_frag = int(fB.lower().replace('f', '')) - 1
    if i_frag < 0 or i_frag >= len(offsets) or j_frag < 0 or j_frag >= len(offsets):
        raise ValueError(f"Fragment indices out of range: f{i_frag+1} or f{j_frag+1} (total fragments: {len(offsets)})")
    i_local = int(aA)
    j_local = int(aB)
    i_global = offsets[i_frag] + i_local
    j_global = offsets[j_frag] + j_local
    rw.AddBond(i_global, j_global, Chem.BondType.SINGLE)
    return (i_global, j_global)

def optimize_gently(mol, iters=75, ff_name="MMFF"):
    confId = mol.GetConformer().GetId()
    
    # Ensure ring info is initialized
    try:
        mol.GetRingInfo().NumRings()
    except:
        try:
            Chem.GetSSSR(mol)
        except Exception as e:
            print(f"[WARN] Could not initialize ring info: {e}, skipping minimization")
            return "SKIPPED"
    
    if ff_name.upper() == "MMFF":
        try:
            props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
            if props is not None:
                try:
                    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=confId)
                    ff.Initialize()
                    ff.Minimize(maxIts=iters)
                    return "MMFF"
                except Exception as e:
                    print(f"[Info] MMFF optimization failed: {e}, trying UFF...")
                    pass
        except Exception as e:
            print(f"[Info] MMFF setup failed: {e}, trying UFF...")
            pass
    
    # Fall back to UFF
    try:
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=confId)
        ff.Initialize()
        ff.Minimize(maxIts=iters)
        return "UFF"
    except RuntimeError as e:
        error_str = str(e)
        if "bad params pointer" in error_str or "UFFTYPER" in error_str or "RingInfo" in error_str:
            print(f"[WARN] Force field cannot handle this structure (likely due to invalid valences from fragment bonding)")
            print(f"[WARN] Skipping minimization - output structure may have geometric issues")
            return "SKIPPED"
        raise
    except (Chem.AtomValenceException, ValueError) as e:
        error_str = str(e)
        if "valence" in error_str.lower() or "Explicit valence" in error_str:
            print(f"[WARN] Invalid valences detected (N-terminus may need H removal after peptide bond formation)")
            print(f"[WARN] Skipping minimization - output structure may have geometric issues")
            print(f"[WARN] Note: This is expected when connecting peptide fragments - valences need post-processing")
            return "SKIPPED"
        raise
    except Exception as e:
        error_str = str(e)
        if "RingInfo" in error_str:
            print(f"[WARN] Ring info issue prevents optimization: {e}")
            return "SKIPPED"
        raise

def main():
    ap = argparse.ArgumentParser(description="Merge MOL2 fragments (3, 4, 5, or more) and do a gentle minimization to relieve clashes.")
    ap.add_argument("frag1", nargs='?', help="First fragment MOL2 file (or use --folder)")
    ap.add_argument("frag2", nargs='?', help="Second fragment MOL2 file")
    ap.add_argument("frag3", nargs='?', help="Third fragment MOL2 file (opcional si se usa --folder)")
    ap.add_argument("-o", "--out", required=True, help="Output MOL2 path")
    ap.add_argument("--folder", help="Folder containing Frag1.mol2, Frag2.mol2, Frag3.mol2, ... (auto-detects all fragments and order from filenames)")
    ap.add_argument("--bond", action="append", default=[], help="Inter-fragment bond like f1:12-f2:7 (0-based indices within each input)")
    ap.add_argument("--no-auto", dest="auto", action="store_false", help="Disable auto-linking C-terminus→N-terminus when no --bond flags given")
    ap.set_defaults(auto=True)
    ap.add_argument("--iters", type=int, default=75, help="Max iterations for minimization")
    ap.add_argument("--ff", choices=["UFF","MMFF"], default="MMFF", help="Force field preference (MMFF tries first, falls back to UFF)")
    ap.add_argument("--clash_scale", type=float, default=DEFAULT_SCALE, help="Clash threshold scale * (vdW_i + vdw_j)")
    args = ap.parse_args()

    # Auto-detect fragment order from filenames if folder is provided
    frag_files = []
    if args.folder:
        folder_path = args.folder
        frag_files_found = {}
        for f in glob.glob(os.path.join(folder_path, "*.mol2")):
            # Excluir archivos de salida (fragmento_final*.mol2)
            basename = os.path.basename(f)
            if basename.startswith("fragmento_final"):
                continue
            frag_num = extract_fragment_number(f)
            if frag_num:
                frag_files_found[frag_num] = f
        
        if len(frag_files_found) < 1:
            print(f"[ERROR] Could not find any fragments in {folder_path}", file=sys.stderr)
            sys.exit(1)
        
        # Obtener todos los fragmentos en orden (1, 2, 3, 4, 5, ...)
        max_frag = max(frag_files_found.keys())
        for i in range(1, max_frag + 1):
            if i not in frag_files_found:
                print(f"[ERROR] Fragment {i} not found in {folder_path} (found: {sorted(frag_files_found.keys())})", file=sys.stderr)
                sys.exit(1)
            frag_files.append(frag_files_found[i])
        print(f"[Info] Auto-detected {len(frag_files)} fragments: {[os.path.basename(f) for f in frag_files]}")
    elif args.frag1 and args.frag2 and args.frag3:
        # Manual order: try to detect and sort by fragment number
        frag_list = [(args.frag1, extract_fragment_number(args.frag1)),
                     (args.frag2, extract_fragment_number(args.frag2)),
                     (args.frag3, extract_fragment_number(args.frag3))]
        
        if all(num is not None for _, num in frag_list):
            frag_list.sort(key=lambda x: x[1])
            frag_files = [f for f, _ in frag_list]
            print(f"[Info] Auto-sorted fragments by filename order")
        else:
            frag_files = [args.frag1, args.frag2, args.frag3]
    else:
        print("[ERROR] Either provide --folder or specify fragment files", file=sys.stderr)
        sys.exit(1)

    # Load MOL2 files and parse atom information for peptide terminus detection
    mols = [load_mol2(f) for f in frag_files]
    mol2_atoms_list = [parse_mol2_atoms(f) for f in frag_files]
    
    rw, offsets, residue_offsets = combine_with_coords(mols)

    added = []
    for flag in args.bond:
        try:
            pair = add_bond_from_flag(rw, flag, offsets)
            added.append(pair)
        except Exception as e:
            print(f"[ERROR] {e}", file=sys.stderr)
            sys.exit(2)

    # Automatic bonding rule when no explicit --bond supplied:
    # Connect C-terminus of previous fragment → N-terminus of next fragment (preserves peptide sequence)
    # Any fragment count: Frag1→Frag2, Frag2→Frag3, Frag3→Frag4, ...
    if not args.bond and args.auto:
        for prev_idx in range(len(mols) - 1):
            next_idx = prev_idx + 1
            # Find C-terminus in previous fragment (C atom with type C.2 in last residue)
            cterm_local = find_cterminus(mol2_atoms_list[prev_idx])
            # Find N-terminus in next fragment (N atom in first residue)
            nterm_local = find_nterminus(mol2_atoms_list[next_idx])
            
            if cterm_local is not None and nterm_local is not None:
                i_global = offsets[prev_idx] + cterm_local
                j_global = offsets[next_idx] + nterm_local
                
                # Get residue names for better reporting
                prev_atom_info = next((a for a in mol2_atoms_list[prev_idx] if a['id'] == cterm_local), None)
                next_atom_info = next((a for a in mol2_atoms_list[next_idx] if a['id'] == nterm_local), None)
                prev_res = prev_atom_info['residue_name'] if prev_atom_info else "UNK"
                next_res = next_atom_info['residue_name'] if next_atom_info else "UNK"
                
                # Check if bond already exists
                if rw.GetBondBetweenAtoms(i_global, j_global) is None:
                    rw.AddBond(i_global, j_global, Chem.BondType.SINGLE)
                    added.append((i_global, j_global))
                    print(f"[Info] Auto-linked peptide bond: {i_global} (C-term {prev_res} of Frag{prev_idx+1}) → {j_global} (N-term {next_res} of Frag{next_idx+1})")
                else:
                    print(f"[Info] Peptide bond already exists between {i_global} and {j_global}")
            else:
                if cterm_local is None:
                    print(f"[WARN] Could not find C-terminus in Frag{prev_idx+1}")
                if nterm_local is None:
                    print(f"[WARN] Could not find N-terminus in Frag{next_idx+1}")
                print(f"[WARN] Auto-link skipped for Frag{prev_idx+1}→Frag{next_idx+1}")
                
                # Fallback: try simple last C → first N approach
                full_mol = rw.GetMol()
                prev_start = offsets[prev_idx]
                next_start = offsets[next_idx]
                prev_end = offsets[prev_idx] + mols[prev_idx].GetNumAtoms()
                next_end = offsets[next_idx] + mols[next_idx].GetNumAtoms()
                
                prev_c_candidates = [i for i in range(prev_start, prev_end) 
                                    if full_mol.GetAtomWithIdx(i).GetSymbol() == 'C']
                next_candidates = [j for j in range(next_start, next_end) 
                                 if full_mol.GetAtomWithIdx(j).GetSymbol() == 'N']
                
                if prev_c_candidates and next_candidates:
                    i_global = prev_c_candidates[-1]
                    j_global = next_candidates[0]
                    if rw.GetBondBetweenAtoms(i_global, j_global) is None:
                        rw.AddBond(i_global, j_global, Chem.BondType.SINGLE)
                        added.append((i_global, j_global))
                        print(f"[Info] Fallback: linked {i_global} (last C of Frag{prev_idx+1}) → {j_global} (first N of Frag{next_idx+1})")

    # Sanitize after adding bonds - skip valence checking to allow bonds that create "invalid" valences
    # (this is common when connecting fragments, especially terminal atoms)
    mol = rw.GetMol()
    try:
        # Try full sanitization first
        Chem.SanitizeMol(mol)
    except Exception:
        # If that fails, try without properties checking only (valence errors are expected when connecting fragments)
        try:
            sanitize_ops = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
            Chem.SanitizeMol(mol, sanitizeOps=sanitize_ops)
            print("[Info] Sanitized with properties checking disabled (fragment bonding may create temporary invalid valences)")
        except Exception as e:
            # Even partial sanitization failed - try to continue anyway
            print(f"[WARN] Sanitization issues (may be OK for fragment merging): {e}")
            # Try to at least update property caches manually
            try:
                mol.UpdatePropertyCache()
            except:
                pass
    
    # Initialize ring info (needed for force field calculations)
    try:
        mol.GetRingInfo().NumRings()
    except:
        # If ring info initialization fails, try to compute it manually
        try:
            Chem.GetSSSR(mol)
        except:
            print("[WARN] Could not initialize ring info - optimization may fail")

    # Check and fix bond lengths of newly formed bonds (< 3.5 Å requirement)
    if len(added) > 0:
        print(f"\n[Info] Checking bond lengths of newly formed bonds (< 3.5 Å requirement)...")
        too_long_before = check_newly_formed_bond_lengths(mol, added, max_length=3.5)
        if len(too_long_before) > 0:
            print(f"[Info] Found {len(too_long_before)} newly formed bond(s) longer than 3.5 Å:")
            for bond_atom1, bond_atom2, bond_length in too_long_before:
                atom1 = mol.GetAtomWithIdx(bond_atom1)
                atom2 = mol.GetAtomWithIdx(bond_atom2)
                atom1_name = atom1.GetProp('_Name') if atom1.HasProp('_Name') else f"{atom1.GetSymbol()}{bond_atom1+1}"
                atom2_name = atom2.GetProp('_Name') if atom2.HasProp('_Name') else f"{atom2.GetSymbol()}{bond_atom2+1}"
                print(f"  Bond {bond_atom1+1} ({atom1_name}) - {bond_atom2+1} ({atom2_name}): {bond_length:.4f} Å")
            
            fixed_count = fix_newly_formed_bond_lengths(mol, added, max_length=3.5, target_length=1.4)
            too_long_after = check_newly_formed_bond_lengths(mol, added, max_length=3.5)
            
            if len(too_long_after) == 0:
                print(f"[Info] ✓ All newly formed bonds are now < 3.5 Å!")
            else:
                print(f"[WARN] {len(too_long_after)} newly formed bond(s) still longer than 3.5 Å after fixing:")
                for bond_atom1, bond_atom2, bond_length in too_long_after:
                    atom1 = mol.GetAtomWithIdx(bond_atom1)
                    atom2 = mol.GetAtomWithIdx(bond_atom2)
                    atom1_name = atom1.GetProp('_Name') if atom1.HasProp('_Name') else f"{atom1.GetSymbol()}{bond_atom1+1}"
                    atom2_name = atom2.GetProp('_Name') if atom2.HasProp('_Name') else f"{atom2.GetSymbol()}{bond_atom2+1}"
                    print(f"  Bond {bond_atom1+1} ({atom1_name}) - {bond_atom2+1} ({atom2_name}): {bond_length:.4f} Å")
        else:
            # Report actual bond lengths
            conf = mol.GetConformer()
            print(f"[Info] Newly formed bond lengths:")
            for bond_atom1, bond_atom2 in added:
                p1 = conf.GetAtomPosition(bond_atom1)
                p2 = conf.GetAtomPosition(bond_atom2)
                dx = p2.x - p1.x
                dy = p2.y - p1.y
                dz = p2.z - p1.z
                bond_length = math.sqrt(dx*dx + dy*dy + dz*dz)
                atom1 = mol.GetAtomWithIdx(bond_atom1)
                atom2 = mol.GetAtomWithIdx(bond_atom2)
                atom1_name = atom1.GetProp('_Name') if atom1.HasProp('_Name') else f"{atom1.GetSymbol()}{bond_atom1+1}"
                atom2_name = atom2.GetProp('_Name') if atom2.HasProp('_Name') else f"{atom2.GetSymbol()}{bond_atom2+1}"
                print(f"  Bond {bond_atom1+1} ({atom1_name}) - {bond_atom2+1} ({atom2_name}): {bond_length:.4f} Å")
            print(f"[Info] ✓ All newly formed bonds are < 3.5 Å!")

    pre = clash_pairs(mol, scale=args.clash_scale)
    print(f"[Info] Potential clashes before minimization: {len(pre)}")
    if len(pre) > 0:
        print("  (showing up to 10) i  j   dist  < cutoff")
        for (i,j,d,cut) in pre[:10]:
            print(f"    {i:4d} {j:4d}  {d:6.3f} < {cut:6.3f}")

    used_ff = optimize_gently(mol, iters=args.iters, ff_name=args.ff)
    print(f"[Info] Minimization done with {used_ff}, iters={args.iters}")

    post = clash_pairs(mol, scale=args.clash_scale)
    print(f"[Info] Potential clashes after minimization: {len(post)}")
    if len(post) > 0:
        print(f"  (showing up to 10) i  j   dist  < {args.clash_scale}*(vdW_i+vdW_j)")
        for (i,j,d,cut) in post[:10]:
            print(f"    {i:4d} {j:4d}  {d:6.3f} < {cut:6.3f}")
    else:
        print(f"[Info] ✓ No clashes remaining!")
    
    # Check and fix atoms that are too close (< 0.5 Å) - ACPYPE will fail on these
    print(f"\n[Info] Checking for atoms too close together (< 0.5 Å - ACPYPE requirement)...")
    too_close_before = check_atoms_too_close(mol, min_distance=0.5)
    if len(too_close_before) > 0:
        print(f"[Info] Found {len(too_close_before)} atom pair(s) too close (< 0.5 Å):")
        for atom_idx1, atom_idx2, dist, is_bonded in too_close_before[:15]:  # Show first 15
            atom1 = mol.GetAtomWithIdx(atom_idx1)
            atom2 = mol.GetAtomWithIdx(atom_idx2)
            bond_type = "BONDED" if is_bonded else "non-bonded"
            # Try to get atom names if available
            atom1_name = atom1.GetProp('_Name') if atom1.HasProp('_Name') else f"{atom1.GetSymbol()}{atom_idx1+1}"
            atom2_name = atom2.GetProp('_Name') if atom2.HasProp('_Name') else f"{atom2.GetSymbol()}{atom_idx2+1}"
            print(f"  Atom {atom_idx1+1} ({atom1_name}) - Atom {atom_idx2+1} ({atom2_name}): {dist:.5f} Å ({bond_type})")
        
        fixed_count = fix_atoms_too_close(mol, min_distance=0.5)
        too_close_after = check_atoms_too_close(mol, min_distance=0.5)
        
        if len(too_close_after) == 0:
            print(f"[Info] ✓ All atoms too close fixed!")
        else:
            print(f"[WARN] {len(too_close_after)} atom pair(s) still too close after fixing:")
            for atom_idx1, atom_idx2, dist, is_bonded in too_close_after:
                atom1 = mol.GetAtomWithIdx(atom_idx1)
                atom2 = mol.GetAtomWithIdx(atom_idx2)
                bond_type = "BONDED" if is_bonded else "non-bonded"
                atom1_name = atom1.GetProp('_Name') if atom1.HasProp('_Name') else f"{atom1.GetSymbol()}{atom_idx1+1}"
                atom2_name = atom2.GetProp('_Name') if atom2.HasProp('_Name') else f"{atom2.GetSymbol()}{atom_idx2+1}"
                print(f"  Atom {atom_idx1+1} ({atom1_name}) - Atom {atom_idx2+1} ({atom2_name}): {dist:.5f} Å ({bond_type})")
    else:
        print(f"[Info] ✓ No atoms too close (< 0.5 Å) - ACPYPE-safe!")
    
    # Check and fix atoms too close to newly formed bonds (< 1.0 Å threshold)
    if len(added) > 0:
        print(f"\n[Info] Checking for atoms too close to newly formed bonds (< 1.0 Å threshold)...")
        problematic_before = check_atoms_near_newly_formed_bonds(mol, added, threshold=1.0)
        if len(problematic_before) > 0:
            print(f"[Info] Found {len(problematic_before)} atom(s) within 1.0 Å of newly formed bonds:")
            for atom_idx, bond_atom1, bond_atom2, dist in problematic_before[:15]:
                atom = mol.GetAtomWithIdx(atom_idx)
                atom_name = atom.GetProp('_Name') if atom.HasProp('_Name') else f"{atom.GetSymbol()}{atom_idx+1}"
                print(f"  Atom {atom_idx+1} ({atom_name}): {dist:.4f} Å from bond {bond_atom1+1}-{bond_atom2+1}")
            
            fixed_count = fix_atoms_near_newly_formed_bonds(mol, added, threshold=1.0)
            problematic_after = check_atoms_near_newly_formed_bonds(mol, added, threshold=1.0)
            
            if len(problematic_after) == 0:
                print(f"[Info] ✓ All atoms too close to newly formed bonds fixed!")
            else:
                print(f"[WARN] {len(problematic_after)} atom(s) still too close to newly formed bonds after fixing:")
                for atom_idx, bond_atom1, bond_atom2, dist in problematic_after:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    atom_name = atom.GetProp('_Name') if atom.HasProp('_Name') else f"{atom.GetSymbol()}{atom_idx+1}"
                    print(f"  Atom {atom_idx+1} ({atom_name}): {dist:.4f} Å from bond {bond_atom1+1}-{bond_atom2+1}")
        else:
            print(f"[Info] ✓ No atoms too close to newly formed bonds (< 1.0 Å) - GROMACS-safe!")

    # Write output: always use custom MOL2 writer to preserve residue information
    try:
        if args.out.lower().endswith(".mol2"):
            write_mol2_with_residues(mol, args.out)
            print(f"[OK] Wrote {args.out} with preserved residue information")
        else:
            # If not .mol2, try RDKit writer
            mol2_writer = getattr(Chem, "MolToMol2File", None)
            if mol2_writer is not None:
                result = mol2_writer(mol, args.out)
                if result is None:
                    print(f"[OK] Wrote {args.out}")
                else:
                    print(f"[WARN] Writer returned {result}; check output.")
            else:
                alt_out = args.out
                if alt_out.lower().endswith(".mol2"):
                    alt_out = alt_out[:-5] + ".mol"
                Chem.MolToMolFile(mol, alt_out)
                print(f"[OK] MOL2 writer unavailable in this RDKit build; wrote MOL to {alt_out}")
    except Exception as e:
        print(f"[ERROR] Failed to write output: {e}", file=sys.stderr)
        sys.exit(3)

if __name__ == "__main__":
    main()
