#!/usr/bin/env python3
"""
Compute RMSD of simulated (MD) peptides against the crystal structure.

For each MD simulation:
1. Extract the peptide from the last frame (L01 or last available frame)
2. Align crystal to simulation using PyMOL
3. Extract peptides from both structures
4. Compute RMSD with obrms (or PyMOL if obrms returns inf):
   - --peptide-rmsd-atoms backbone (default): N, CA, C, O (intermediate choice)
   - --peptide-rmsd-atoms heavy: all heavy atoms
   - --peptide-rmsd-atoms ca: CA atoms only
5. Write a CSV with the results
"""

import os
import re
import glob
import csv
import subprocess
import sys
import argparse
import tempfile
from pathlib import Path


def find_vs_folders(project_root="."):
    """
    Find all VS_GR_* folders under the project root.

    Args:
        project_root: Project root directory

    Returns:
        list: Paths to VS_GR_* folders
    """
    pattern = os.path.join(project_root, "VS_GR_*")
    vs_folders = [f for f in glob.glob(pattern) if os.path.isdir(f)]
    # Also search the organized location: step 17 (organize_workflow_results.py) moves
    # each case's MD runs into <CASE>_results/md_simulations/ to keep the project root
    # uncluttered. Include them so re-running step 13 (recompute) after organize still
    # finds the simulations. Cases that share a pdbid land in one case's folder, so we
    # scan every *_results/md_simulations and let the per-case pdbid filter select.
    organized = os.path.join(project_root, "*_results", "md_simulations", "VS_GR_*")
    vs_folders += [f for f in glob.glob(organized) if os.path.isdir(f)]
    return sorted(set(vs_folders))


def find_simulation_peptide(vs_folder):
    """
    Find the simulated peptide PDB (last frame).
    Searches in:
    1. molecules/ under VS_GR
    2. query_combination_X/ under VS_GR
    3. Any subdirectory of VS_GR (last resort)

    Args:
        vs_folder: VS_GR_* folder

    Returns:
        str: Path to simulated peptide PDB or None
    """
    # 1. Search in molecules/
    molecules_dir = os.path.join(vs_folder, "molecules")
    if os.path.exists(molecules_dir):
        pdb_files = glob.glob(os.path.join(molecules_dir, "*.pdb"))
        if pdb_files:
            # Prefer files with "no_solvent" or "md" in the name
            preferred = [f for f in pdb_files if "no_solvent" in f or "_md" in f]
            if preferred:
                return preferred[0]
            return pdb_files[0]
    
    # 2. Search under query_combination_X/ inside VS_GR
    query_dirs = glob.glob(os.path.join(vs_folder, "query_combination_*"))
    for query_dir in query_dirs:
        # Look for *_complex_md*.pdb or *_md_no_solvent.pdb
        patterns = [
            os.path.join(query_dir, "*_complex_md_no_solvent.pdb"),
            os.path.join(query_dir, "*_complex_md.pdb"),
            os.path.join(query_dir, "*_md_no_solvent.pdb"),
            os.path.join(query_dir, "*_md.pdb"),
            os.path.join(query_dir, "*.pdb"),
        ]
        for pattern in patterns:
            matches = glob.glob(pattern)
            if matches:
                # Prefer no_solvent
                preferred = [f for f in matches if "no_solvent" in f]
                if preferred:
                    return preferred[0]
                return matches[0]
    
    # 3. Recursive search under VS_GR (last resort)
    pdb_files = glob.glob(os.path.join(vs_folder, "**", "*.pdb"), recursive=True)
    if pdb_files:
        # Prefer files with "no_solvent", "md", or "complex" in the name
        preferred = [f for f in pdb_files if any(x in f for x in ["no_solvent", "_md", "complex"])]
        if preferred:
            return preferred[0]
        return pdb_files[0]
    
    return None


def find_crystal_structure(crystal_base_dir, combo_dir=None, project_root=".", case_name=None, preferred_peptide_chain=None):
    """
    Locate the full crystal structure file.
    Search order:
    1. complex/ under project root (complex/{case}_*.pdb)
    2. query_combination_X/complex/ inside the combination folder
    3. peptide_pdb_fragments/{CASE}/

    Args:
        crystal_base_dir: Base directory for crystal inputs
        combo_dir: Combination directory (optional, for local crystal)
        project_root: Project root
        case_name: Case id (e.g. "1CJR_A") to pick the right PDB
                   Filename pattern: {case_base}_{peptide_chain}_{protein_chain}.pdb
        preferred_peptide_chain: If --chain was set, prefer PDB where peptide is on this chain

    Returns:
        tuple: (path to full crystal PDB, peptide chain id) or (None, None)
    """
    # 1. complex/ at project root (highest priority)
    complex_root = os.path.join(project_root, "complex")
    if os.path.exists(complex_root):
        # With case name, search for that PDB specifically
        if case_name:
            case_lower = case_name.lower()
            case_upper = case_name.upper()
            # Case may be "1CJR_A" or "1cjr_B"
            if "_" in case_lower:
                case_parts = case_lower.split("_")
                case_base = case_parts[0]  # 1cjr
                case_chain = case_parts[1] if len(case_parts) > 1 else None  # A or B
            else:
                case_base = case_lower
                case_chain = None

            # Glob all PDBs starting with case base; pattern {case_base}_{*}_*.pdb
            patterns = [
                os.path.join(complex_root, f"{case_base.upper()}_*.pdb"),  # 1CJR_*.pdb
                os.path.join(complex_root, f"{case_base}_*.pdb"),  # 1cjr_*.pdb
            ]
            
            all_matches = []
            for pattern in patterns:
                all_matches.extend(glob.glob(pattern))
            
            if all_matches:
                # {case_base}_{peptide_chain}_{protein_chain}.pdb → parts[2] = protein chain
                def peptide_chain_from_path(path):
                    parts = os.path.basename(path).replace(".pdb", "").split("_")
                    return parts[1] if len(parts) >= 2 else None
                def protein_chain_from_path(path):
                    parts = os.path.basename(path).replace(".pdb", "").split("_")
                    return parts[2].upper() if len(parts) >= 3 else None
                # Priority 1: match the case suffix to the correct crystal file.
                # The suffix is AMBIGUOUS across datasets: it can be the PEPTIDE
                # chain (e.g. improve5Frag_2: 2V8W_B -> 2v8w_B_A.pdb) or the
                # PROTEIN chain (legacy: 1AQD_A -> 1aqd_C_A.pdb). Disambiguate by
                # what actually exists: prefer a file whose PEPTIDE chain matches
                # the suffix (that is the reference peptide we must compare), and
                # only fall back to protein-chain matching. Picking the wrong file
                # here compares against a peptide bound to a DIFFERENT receptor
                # copy, which sits far away after alignment and inflates rmsd_md to
                # tens of Angstrom (e.g. 2V8W_B: wrong chain F -> 44.6 A vs correct
                # chain B -> 3.7 A).
                if case_chain:
                    cc = case_chain.upper()
                    by_peptide_case = [f for f in all_matches
                                       if (peptide_chain_from_path(f) or "").upper() == cc]
                    by_protein = [f for f in all_matches
                                  if protein_chain_from_path(f) == cc]
                    if by_peptide_case:
                        pool = by_peptide_case
                    elif by_protein:
                        pool = by_protein
                    else:
                        pool = all_matches
                else:
                    pool = all_matches
                # Priority 2: if user passed --chain X, prefer PDB with peptide on chain X
                if preferred_peptide_chain:
                    by_peptide = [f for f in pool if peptide_chain_from_path(f) == preferred_peptide_chain.upper()]
                    pool = by_peptide if by_peptide else pool
                # Priority 3: else prefer peptide on chain C (e.g. 1aqd_C_A.pdb)
                preferred_c = [f for f in pool if peptide_chain_from_path(f) == "C"]
                match = preferred_c[0] if preferred_c else pool[0]
                filename = os.path.basename(match)
                parts = filename.replace(".pdb", "").split("_")
                if len(parts) >= 2:
                    peptide_chain = parts[1]
                    return match, peptide_chain
                return match, "A"
        
        # Only without case_name: any PDB (avoid wrong crystal when case is set but file missing)
        if not case_name:
            patterns = [
                os.path.join(complex_root, "*.pdb"),
            ]
            for pattern in patterns:
                matches = glob.glob(pattern)
                if matches:
                    filename = os.path.basename(matches[0])
                    parts = filename.replace(".pdb", "").split("_")
                    if len(parts) >= 2:
                        peptide_chain = parts[1]
                        return matches[0], peptide_chain
                    return matches[0], "A"

    # 2. Combination folder (if present)
    if combo_dir:
        combo_number = os.path.basename(combo_dir).replace("combination_", "")
        query_dir = os.path.join(combo_dir, f"query_combination_{combo_number}")
        # Use case_name to prefer case PDB (e.g. 1AQD_C -> 1aqd*.pdb)
        case_base = case_name.split("_")[0].lower() if case_name and "_" in case_name else (case_name.lower() if case_name else "")
        complex_patterns = [
            os.path.join(query_dir, "complex", "*.pdb"),
            os.path.join(query_dir, "complex", "*.mol2"),
            os.path.join(query_dir, "*complex*.pdb"),
            os.path.join(query_dir, "*complex*.mol2"),
        ]
        if case_base:
            complex_patterns.insert(0, os.path.join(query_dir, "complex", f"{case_base}*.pdb"))
        for pattern in complex_patterns:
            matches = glob.glob(pattern)
            if matches:
                by_case = [f for f in matches if case_base and case_base in os.path.basename(f).lower()]
                pool = by_case if by_case else matches
                def chain_from_path(p):
                    parts = os.path.basename(p).replace(".pdb", "").replace(".mol2", "").split("_")
                    return parts[1] if len(parts) >= 2 else None
                def protein_chain_from_path_combo(p):
                    parts = os.path.basename(p).replace(".pdb", "").replace(".mol2", "").split("_")
                    return parts[2].upper() if len(parts) >= 3 else None
                case_chain_combo = case_name.split("_")[1].upper() if case_name and "_" in case_name else None
                # Priority 1: match suffix to PEPTIDE chain first (improve5Frag_2),
                # else PROTEIN chain (legacy). See note in the complex/ branch.
                if case_chain_combo:
                    by_peptide_case = [f for f in pool if (chain_from_path(f) or "").upper() == case_chain_combo]
                    by_protein = [f for f in pool if protein_chain_from_path_combo(f) == case_chain_combo]
                    if by_peptide_case:
                        pool = by_peptide_case
                    elif by_protein:
                        pool = by_protein
                # Priority 2: if --chain X, prefer PDB with peptide on chain X
                if preferred_peptide_chain:
                    by_peptide = [f for f in pool if chain_from_path(f) == preferred_peptide_chain.upper()]
                    pool = by_peptide if by_peptide else pool
                preferred_c = [f for f in pool if chain_from_path(f) == "C"]
                selected_file = preferred_c[0] if preferred_c else pool[0]
                filename = os.path.basename(selected_file)
                parts = filename.replace(".pdb", "").replace(".mol2", "").split("_")
                if len(parts) >= 2:
                    peptide_chain = parts[1]
                    return selected_file, peptide_chain
                return selected_file, "A"
    
    # 3. Full complex in peptide_pdb_fragments using case_name
    if case_name:
        case_upper = case_name.upper()
        crystal_patterns = [
            os.path.join(crystal_base_dir, case_upper, "*complete*.pdb"),
            os.path.join(crystal_base_dir, case_upper, "*full*.pdb"),
            os.path.join(crystal_base_dir, case_upper, "*.pdb"),
        ]
        for pattern in crystal_patterns:
            matches = glob.glob(pattern)
            if matches:
                # Prefer non-fragment files (exclude Frag/FragmentoN)
                non_frag = [f for f in matches if "Frag" not in f and "Fragment" not in f]
                selected_file = non_frag[0] if non_frag else matches[0]
                filename = os.path.basename(selected_file)
                parts = filename.replace(".pdb", "").split("_")
                if len(parts) >= 2:
                    peptide_chain = parts[1]
                    return selected_file, peptide_chain
                return selected_file, "A"
    
    # 4. Fallback: peptide_pdb_fragments/1A1M_C (legacy)
    crystal_patterns = [
        os.path.join(crystal_base_dir, "1A1M_C", "*complete*.pdb"),
        os.path.join(crystal_base_dir, "1A1M_C", "*full*.pdb"),
        os.path.join(crystal_base_dir, "1A1M_C", "*.pdb"),
    ]
    for pattern in crystal_patterns:
        matches = glob.glob(pattern)
        if matches:
            non_frag = [f for f in matches if "Frag" not in f and "Fragment" not in f]
            selected_file = non_frag[0] if non_frag else matches[0]
            filename = os.path.basename(selected_file)
            parts = filename.replace(".pdb", "").split("_")
            if len(parts) >= 2:
                peptide_chain = parts[1]
                return selected_file, peptide_chain
            return selected_file, "A"
    
    return None, None


def _crystal_peptide_sequence(path, chain):
    """Ordered residue names of `chain` in a crystal PDB (its peptide sequence).

    Used to confirm two crystal files hold the SAME peptide before treating them as
    interchangeable lattice copies (see find_crystal_candidates). Returns None if the
    chain has no residues or the file can't be read."""
    chain = (chain or "").upper()
    seq, seen = [], set()
    try:
        with open(path) as fh:
            for line in fh:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                if line[21:22].upper() != chain:
                    continue
                res_key = line[22:27]  # resSeq + insertion code
                if res_key in seen:
                    continue
                seen.add(res_key)
                seq.append(line[17:20].strip())
    except OSError:
        return None
    return seq or None


def find_crystal_candidates(project_root, case_name):
    """All complex/ crystal files that could be the reference for this case.

    A crystal frequently contains several lattice copies of the same complex: the
    same peptide chain bound to different protein chains (e.g. 1k7l_D_A.pdb and
    1k7l_D_C.pdb). Only the copy the docking corresponds to superposes correctly;
    the others sit tens of Angstrom away and inflate the RMSD. We therefore return
    ALL files whose PEPTIDE chain matches the case suffix so the caller can try
    each and keep the lowest RMSD. Falls back to protein-chain match, then to every
    file for the PDB. Each returned path follows {base}_{peptide}_{protein}.pdb."""
    complex_root = os.path.join(project_root, "complex")
    if not case_name or "_" not in case_name or not os.path.isdir(complex_root):
        return []
    base, suffix = case_name.split("_", 1)
    base_l, suffix_u = base.lower(), suffix.upper()
    files = sorted(set(glob.glob(os.path.join(complex_root, base_l + "_*.pdb")) +
                       glob.glob(os.path.join(complex_root, base_l.upper() + "_*.pdb"))))

    def _pep(p):
        parts = os.path.basename(p).replace(".pdb", "").split("_")
        return parts[1].upper() if len(parts) >= 2 else None

    def _prot(p):
        parts = os.path.basename(p).replace(".pdb", "").split("_")
        return parts[2].upper() if len(parts) >= 3 else None

    by_pep = [f for f in files if _pep(f) == suffix_u]
    if by_pep:
        # Symmetry-sibling copies: a crystal lattice often stores identical peptide
        # copies under DIFFERENT chain labels (e.g. 6g0y peptide as chains H/I/J).
        # The docking may land on any copy, so a pose that actually matched a sibling
        # copy scores tens of A on obrms against the case's own chain despite a ~0.4 A
        # PyMOL align (correct conformation, different lattice position). Include the
        # siblings so the caller can pick the closest copy. Guard by EXACT peptide
        # sequence identity: backbone obrms is sequence-agnostic and would otherwise
        # be fooled into matching a genuinely different peptide chain that happens to
        # sit nearby. Single-copy cases (no sibling with identical sequence) keep the
        # original single-file behaviour untouched.
        ref_seq = _crystal_peptide_sequence(by_pep[0], _pep(by_pep[0]))
        if ref_seq:
            siblings = [
                f for f in files
                if f not in by_pep and _crystal_peptide_sequence(f, _pep(f)) == ref_seq
            ]
            if siblings:
                return sorted(set(by_pep + siblings))
        return by_pep
    by_prot = [f for f in files if _prot(f) == suffix_u]
    if by_prot:
        return by_prot
    return files


def crystal_peptide_chain(path):
    """Peptide chain id from a {base}_{peptide}_{protein}.pdb crystal filename."""
    parts = os.path.basename(path).replace(".pdb", "").split("_")
    return parts[1] if len(parts) >= 2 else "A"


def extract_peptide_with_pymol(simulation_file, crystal_file, output_peptide_sim, output_peptide_crystal, output_protein_sim=None, singularity_image="singularity/new_ms.simg", peptide_chain="D", force_peptide_chain=False):
    """
    Use PyMOL to align crystal to simulation and extract peptide and receptor.

    Args:
        simulation_file: Simulation PDB
        crystal_file: Crystal PDB
        output_peptide_sim: Output path for simulation peptide
        output_peptide_crystal: Output path for crystal peptide (aligned)
        output_protein_sim: Optional output path for simulation receptor
        singularity_image: Singularity image path
        peptide_chain: Peptide chain id in crystal (from filename or --chain)
        force_peptide_chain: If True, do not auto-detect chain (--chain)

    Returns:
        bool: True on success
    """
    # Temporary PyMOL driver script; use absolute paths
    sim_abs = os.path.abspath(simulation_file)
    crystal_abs = os.path.abspath(crystal_file)
    out_sim_abs = os.path.abspath(output_peptide_sim)
    out_crystal_abs = os.path.abspath(output_peptide_crystal)
    out_protein_sim_abs = os.path.abspath(output_protein_sim) if output_protein_sim else None
    
    # PyMOL driver script (.format inserts paths; brace-doubling for pymol iterate)
    python_script = """# -*- coding: utf-8 -*-
import sys
import pymol
pymol.finish_launching(['pymol', '-c', '-Q'])

pymol.cmd.load('{sim_abs}', 'sim')
pymol.cmd.load('{crystal_abs}', 'crystal')

num_models = pymol.cmd.count_states('sim')
if num_models > 1:
    pymol.cmd.frame(num_models)

solvent_ion = 'resn HOH or resn WAT or resn SOL or resn TIP or resn TIP3 or resn TIP4 or resn NA or resn CL or resn K or resn MG or resn ZN or resn MN or resn CAL or resn CA or resn CA2'
pymol.cmd.remove('sim and (' + solvent_ion + ')')
pymol.cmd.remove('crystal and (' + solvent_ion + ')')

align_result = pymol.cmd.align('sim', 'crystal', cycles=1)
if align_result is not None and len(align_result) > 0:
    print('RMSD_ALIGN_PYMOL:', round(align_result[0], 4))
sys.stdout.flush()

peptide_sim_selector = 'sim and resn L01'
if pymol.cmd.count_atoms(peptide_sim_selector) == 0:
    peptide_sim_selector = 'sim and chain L'

force_chain = {force_peptide_chain}
peptide_crystal_chain = '{peptide_chain}'
peptide_crystal_selector = 'crystal and chain ' + peptide_crystal_chain
if pymol.cmd.count_atoms(peptide_crystal_selector) == 0:
    if force_chain:
        sys.stderr.write('ERROR: Requested chain (' + peptide_crystal_chain + ') not found or empty in crystal. Auto-detection disabled (--chain).\\n')
        sys.exit(1)
    chain_counts_dict = {{}}
    pymol.cmd.iterate('crystal and name CA', 'chain_counts_dict[chain] = chain_counts_dict.get(chain, 0) + 1', space={{'chain_counts_dict': chain_counts_dict}})
    if len(chain_counts_dict) > 1:
        peptide_crystal_chain = min(chain_counts_dict, key=chain_counts_dict.get)
    elif len(chain_counts_dict) == 1:
        peptide_crystal_chain = list(chain_counts_dict.keys())[0]
    else:
        if pymol.cmd.count_atoms('crystal and chain L') > 0:
            peptide_crystal_chain = 'L'
        elif pymol.cmd.count_atoms('crystal and chain C') > 0:
            peptide_crystal_chain = 'C'
        else:
            sys.stderr.write('ERROR: Could not identify peptide chain in crystal\\n')
            sys.exit(1)
    peptide_crystal_selector = 'crystal and chain ' + peptide_crystal_chain

pymol.cmd.select('peptide_sim_sel', peptide_sim_selector)
pymol.cmd.select('peptide_crystal_sel', peptide_crystal_selector)

if pymol.cmd.count_atoms('peptide_sim_sel') == 0 or pymol.cmd.count_atoms('peptide_crystal_sel') == 0:
    sys.stderr.write('ERROR: Could not select peptides after alignment\\n')
    sys.exit(1)

pymol.cmd.create('peptide_sim', 'peptide_sim_sel')
pymol.cmd.create('peptide_crystal', 'peptide_crystal_sel')

pymol.cmd.select('peptide_sim_ca_sel', 'peptide_sim and name CA')
pymol.cmd.select('peptide_crystal_ca_sel', 'peptide_crystal and name CA')

if pymol.cmd.count_atoms('peptide_sim_ca_sel') == 0 or pymol.cmd.count_atoms('peptide_crystal_ca_sel') == 0:
    sys.stderr.write('ERROR: No CA atoms found in peptides\\n')
    sys.exit(1)

pymol.cmd.create('peptide_sim_ca', 'peptide_sim_ca_sel')
pymol.cmd.create('peptide_crystal_ca', 'peptide_crystal_ca_sel')

pymol.cmd.save('{out_sim_abs}', 'peptide_sim')
pymol.cmd.save('{out_crystal_abs}', 'peptide_crystal')

import os as os_module
out_sim_ca = os_module.path.splitext('{out_sim_abs}')[0] + '_ca.pdb'
out_crystal_ca = os_module.path.splitext('{out_crystal_abs}')[0] + '_ca.pdb'
pymol.cmd.save(out_sim_ca, 'peptide_sim_ca')
pymol.cmd.save(out_crystal_ca, 'peptide_crystal_ca')

# Save backbone (N, CA, C, O) versions for obrms (backbone vs backbone)
bb_sel = 'name N or name CA or name C or name O'
pymol.cmd.select('peptide_sim_bb_sel', 'peptide_sim and (' + bb_sel + ')')
pymol.cmd.select('peptide_crystal_bb_sel', 'peptide_crystal and (' + bb_sel + ')')
if pymol.cmd.count_atoms('peptide_sim_bb_sel') > 0 and pymol.cmd.count_atoms('peptide_crystal_bb_sel') > 0:
    pymol.cmd.create('peptide_sim_bb', 'peptide_sim_bb_sel')
    pymol.cmd.create('peptide_crystal_bb', 'peptide_crystal_bb_sel')
    out_sim_bb = os_module.path.splitext('{out_sim_abs}')[0] + '_backbone.pdb'
    out_crystal_bb = os_module.path.splitext('{out_crystal_abs}')[0] + '_backbone.pdb'
    pymol.cmd.save(out_sim_bb, 'peptide_sim_bb')
    pymol.cmd.save(out_crystal_bb, 'peptide_crystal_bb')

# Simulation receptor only (exclude peptide, water, ions)
if '{out_protein_sim_abs}' and '{out_protein_sim_abs}' != 'None':
    solvent_ion_resn = 'resn HOH or resn WAT or resn SOL or resn TIP or resn TIP3 or resn TIP4 or resn NA or resn CL or resn K or resn MG or resn ZN or resn MN or resn CAL or resn CA or resn CA2'
    pymol.cmd.select('protein_sim_full', 'sim and not (' + peptide_sim_selector + ') and not (' + solvent_ion_resn + ')')
    if pymol.cmd.count_atoms('protein_sim_full') > 0:
        pymol.cmd.save('{out_protein_sim_abs}', 'protein_sim_full')
""".format(
        sim_abs=sim_abs,
        crystal_abs=crystal_abs,
        out_sim_abs=out_sim_abs,
        out_crystal_abs=out_crystal_abs,
        out_protein_sim_abs=out_protein_sim_abs if out_protein_sim_abs else 'None',
        peptide_chain=peptide_chain,
        force_peptide_chain=force_peptide_chain
    )
    
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False, encoding='utf-8') as f:
            f.write(python_script)
            python_script_path = f.name

        python_script_abs = os.path.abspath(python_script_path)
        cmd = f'singularity exec {singularity_image} python3 {python_script_abs}'

        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            check=False,
            timeout=300
        )

        if os.path.exists(python_script_path):
            os.remove(python_script_path)

        if os.path.exists(output_peptide_sim) and os.path.exists(output_peptide_crystal):
            if result.stdout:
                for line in result.stdout.strip().splitlines():
                    if 'RMSD_ALIGN_PYMOL:' in line:
                        print(f"  {line.strip()}")
            return True
        else:
            error_msg = result.stderr.strip() if result.stderr else ""
            stdout_msg = result.stdout.strip() if result.stdout else ""
            
            all_output = (error_msg + "\n" + stdout_msg).strip()
            error_lines = [line for line in all_output.split('\n') 
                          if line and "INFO:" not in line and "Converting SIF" not in line 
                          and "Cleaning up" not in line and "PyMOL" not in line]
            
            if error_lines:
                error_preview = '\n'.join(error_lines[-5:])
                print(f"  PyMOL error: {error_preview}")
            elif result.returncode != 0:
                print(f"  PyMOL error: exit code {result.returncode}")

            return False

    except subprocess.TimeoutExpired:
        print(f"  PyMOL error: timeout (> 5 minutes)")
        if 'python_script_path' in locals() and os.path.exists(python_script_path):
            os.remove(python_script_path)
        return False
    except Exception as e:
        print(f"  PyMOL run error: {e}")
        if 'python_script_path' in locals() and os.path.exists(python_script_path):
            os.remove(python_script_path)
        return False


def _parse_pdb_atoms(path):
    """Parse PDB into a list of atom dicts (record, name, resn, resi, chain, element, x, y, z)."""
    atoms = []
    elem_map = {'C': 'C', 'N': 'N', 'O': 'O', 'S': 'S', 'P': 'P', 'H': 'H', 'CA': 'C'}
    with open(path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            record = line[:6].strip()
            try:
                serial = int(line[6:11])
                name = line[12:16].strip()
                resn = line[17:21].strip()
                chain = line[21:22].strip() if len(line) > 21 else ''
                resi = line[22:26].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except (ValueError, IndexError):
                continue
            elem = line[76:78].strip() if len(line) >= 78 else ''
            if not elem and name:
                elem = elem_map.get(name[:1], name[:1]) if name[0].isalpha() else 'C'
            atoms.append({
                'record': record, 'name': name, 'resn': resn, 'resi': resi, 'chain': chain,
                'element': elem or 'C', 'x': x, 'y': y, 'z': z
            })
    return atoms


def _dist2(a, b):
    return (a['x'] - b['x']) ** 2 + (a['y'] - b['y']) ** 2 + (a['z'] - b['z']) ** 2


def _match_atoms_by_distance(ref_atoms, full_atoms, same_element=True):
    """Map each ref atom to nearest in full_atoms (same element); 1:1 without reuse."""
    used = set()
    mapping = []
    for ra in ref_atoms:
        best_j = None
        best_d2 = float('inf')
        for j, fa in enumerate(full_atoms):
            if j in used:
                continue
            if same_element and ra['element'] != fa['element']:
                continue
            d2 = _dist2(ra, fa)
            if d2 < best_d2:
                best_d2 = d2
                best_j = j
        if best_j is not None:
            used.add(best_j)
            mapping.append(best_j)
        else:
            mapping.append(None)
    return mapping


def _write_pdb(path, atoms_template, coords_from, indices):
    """Escribe PDB con nombres/residuos de atoms_template y coordenadas de coords_from[indices[i]].
    Standard PDB columns: chain 22, resi 23-26, element 77-78 (Open Babel).
    """
    with open(path, 'w') as f:
        for i, t in enumerate(atoms_template):
            j = indices[i] if i < len(indices) else i
            if j is None or j >= len(coords_from):
                continue
            c = coords_from[j]
            rec = t.get('record', 'ATOM')
            if len(rec) < 6:
                rec = rec.ljust(6)
            name = (t['name'][:4] if len(t['name']) <= 4 else t['name']).ljust(4)
            resn = (t['resn'][:3] if len(t['resn']) >= 3 else t['resn']).ljust(3)
            chain = t.get('chain', '')[:1] or ' '
            resi = str(t.get('resi', ''))[:4].rjust(4)
            elem = (c.get('element') or t.get('element') or 'C')[:2].strip().rjust(2)
            # PDB: cols 18-20 resn, 22 chain, 23-26 resi, 77-78 element
            occ, tf = 1.00, 0.00
            line = (f"{rec:6s}{i+1:5d} {name}{resn} {chain}{resi}   "
                    f"{c['x']:8.3f}{c['y']:8.3f}{c['z']:8.3f}{occ:6.2f}{tf:6.2f}            {elem}\n")
            f.write(line)
        f.write("END\n")


def _calculate_rmsd_pymol_selector(ref_pdb, query_pdb, singularity_image, selector):
    """
    Calcula RMSD usando PyMOL align con un selector dado.
    selector: string PyMOL (ej. 'backbone', 'not elem H', 'name CA')
    """
    ref_abs = os.path.abspath(ref_pdb)
    query_abs = os.path.abspath(query_pdb)
    script = f"""# -*- coding: utf-8 -*-
import sys
import pymol
pymol.finish_launching(['pymol', '-c', '-Q'])
pymol.cmd.load('{ref_abs}', 'ref')
pymol.cmd.load('{query_abs}', 'query')
ref_sel = pymol.cmd.count_atoms('ref and ({selector})')
query_sel = pymol.cmd.count_atoms('query and ({selector})')
if ref_sel == 0 or query_sel == 0:
    sys.stderr.write('No atoms matching selector')
    sys.exit(1)
r = pymol.cmd.align('query and ({selector})', 'ref and ({selector})', cycles=0)
if r and len(r) > 0:
    print('RMSD:', round(float(r[0]), 4))
""".format(selector=selector)
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False, encoding='utf-8') as f:
            f.write(script)
            script_path = f.name
        result = subprocess.run(
            f'singularity exec {singularity_image} python3 {script_path}',
            shell=True, capture_output=True, text=True, timeout=60
        )
        if os.path.exists(script_path):
            os.remove(script_path)
        if result.returncode != 0:
            return None
        for line in (result.stdout or '').splitlines():
            if 'RMSD:' in line:
                try:
                    return float(line.split(':', 1)[1].strip())
                except (ValueError, IndexError):
                    pass
        return None
    except Exception:
        return None


BACKBONE_NAMES = frozenset(['N', 'CA', 'C', 'O'])


def _filter_atoms_by_mode(atoms, mode):
    """Filter atoms: 'heavy' (non-H), 'backbone' (N, CA, C, O), 'ca' (CA only)."""
    if mode == 'heavy':
        return [a for a in atoms if a['element'] != 'H']
    if mode == 'backbone':
        return [a for a in atoms if a['name'] in BACKBONE_NAMES]
    if mode == 'ca':
        return [a for a in atoms if a['name'] == 'CA']
    return atoms


# --------------------------------------------------------------------------- #
# Robust backbone obrms (residue+name matching + shared connectivity)
# --------------------------------------------------------------------------- #
# Why: OpenBabel perceives bonds from interatomic distances, so a MD-relaxed
# peptide and its crystal reference can end up with DIFFERENT molecular graphs
# and obrms returns 'inf'. The previous nearest-neighbour matcher also scrambled
# the atom correspondence for displaced peptides. We instead build a FAITHFUL
# (residue, backbone-atom-name) correspondence and force IDENTICAL connectivity
# (perceived once from the crystal) so obrms always returns the honest in-place
# RMSD -- no PyMOL fallback needed. See agentic_loop_stellar/recompute_obrms.py.
_BB_ORDER = ['N', 'CA', 'C', 'O']
_BB_RANK = {n: i for i, n in enumerate(_BB_ORDER)}


def _bb_res_groups_crystal(atoms):
    """Crystal -> per-residue backbone groups (resi order), each ordered N,CA,C,O.
    Only residues with the complete backbone are kept."""
    sel = [a for a in atoms if a['name'] in _BB_RANK]

    def key(a):
        try:
            r = int(a['resi'])
        except (ValueError, TypeError):
            r = 0
        return (r, _BB_RANK[a['name']])

    sel.sort(key=key)
    res, cur, buf = [], None, []
    for a in sel:
        if a['resi'] != cur:
            if buf:
                res.append(buf)
            cur, buf = a['resi'], []
        buf.append(a)
    if buf:
        res.append(buf)
    return [g for g in res if len(g) == len(_BB_ORDER)]


def _bb_res_groups_sim(atoms):
    """Sim (single L01 residue) -> per-residue backbone groups in sequence order.
    PyMOL emits the backbone selection grouped by NAME (all N, then all CA, ...),
    so interleave the by-name groups to recover residue order."""
    grp = {n: [a for a in atoms if a['name'] == n] for n in _BB_ORDER}
    nres = len(grp['CA'])
    out = []
    for i in range(nres):
        g = [grp[n][i] for n in _BB_ORDER if i < len(grp[n])]
        if len(g) == len(_BB_ORDER):
            out.append(g)
    return out


def _bb_group_rmsd(cry_res, sim_res):
    """In-place RMSD over name-matched backbone atoms of paired residue groups."""
    n = s = 0
    for gc, gs in zip(cry_res, sim_res):
        for ac, as_ in zip(gc, gs):
            s += _dist2(ac, as_)
            n += 1
    return (s / n) ** 0.5 if n else float('inf')


def _bb_best_offset(cry_res, sim_res):
    """Crystal peptides frequently miss terminal residues (unresolved density), so
    the crystal is a CONTIGUOUS block shorter than the sim. Slide it over the sim
    residues and pick the best structural alignment; this recovers the true
    sequence offset (verified on 6IVX: crystal 2351 vs sim 2346 -> offset 5)."""
    best = None
    for off in range(len(sim_res) - len(cry_res) + 1):
        d = _bb_group_rmsd(cry_res, sim_res[off:off + len(cry_res)])
        if best is None or d < best[1]:
            best = (off, d)
    return best if best else (0, float('inf'))


def _bb_align_blocks(cry_res, sim_res):
    """Align crystal vs sim backbone residue blocks of possibly DIFFERENT length.

    Either side may be the shorter contiguous block: the crystal can miss
    unresolved terminal residues (crystal shorter), and the MD-built peptide can
    be one/few residues short of the crystal reference (sim shorter, e.g. 4P3W_G:
    crystal 15 res vs sim 14 res). Slide the shorter block over the longer one and
    keep the equal-length window with the best structural overlap. Returns
    (cry_sub, sim_sub) of identical length, ordered residue-by-residue N,CA,C,O."""
    L = min(len(cry_res), len(sim_res))
    if L == 0:
        return [], []
    best = None  # (cry_offset, sim_offset, rmsd)
    if len(cry_res) <= len(sim_res):
        for off in range(len(sim_res) - L + 1):
            d = _bb_group_rmsd(cry_res, sim_res[off:off + L])
            if best is None or d < best[2]:
                best = (0, off, d)
    else:
        for off in range(len(cry_res) - L + 1):
            d = _bb_group_rmsd(cry_res[off:off + L], sim_res)
            if best is None or d < best[2]:
                best = (off, 0, d)
    coff, soff, _ = best
    return cry_res[coff:coff + L], sim_res[soff:soff + L]


def _obrms_shared_connectivity(cry_atoms, sim_atoms, output_dir, combo_number, singularity_image):
    """obrms with connectivity perceived once from the crystal and shared with the
    sim (same atom order) -> can never return inf. Returns float or None."""
    cry_pdb = os.path.join(output_dir, f"peptide_{combo_number}_bbmatch_cry.pdb")
    sim_pdb = os.path.join(output_dir, f"peptide_{combo_number}_bbmatch_sim.pdb")
    ref_sdf = os.path.join(output_dir, f"peptide_{combo_number}_bbmatch_ref.sdf")
    sim_xyz = os.path.join(output_dir, f"peptide_{combo_number}_bbmatch_sim.xyz")
    sim_sdf = os.path.join(output_dir, f"peptide_{combo_number}_bbmatch_sim.sdf")
    idx = list(range(len(cry_atoms)))
    _write_pdb(cry_pdb, cry_atoms, cry_atoms, idx)
    _write_pdb(sim_pdb, sim_atoms, sim_atoms, idx)
    try:
        subprocess.run(f'singularity exec {singularity_image} obabel "{cry_pdb}" -osdf -O "{ref_sdf}"',
                       shell=True, capture_output=True, text=True, timeout=60)
        subprocess.run(f'singularity exec {singularity_image} obabel "{sim_pdb}" -oxyz -O "{sim_xyz}"',
                       shell=True, capture_output=True, text=True, timeout=60)
    except Exception:
        return None
    if not (os.path.isfile(ref_sdf) and os.path.isfile(sim_xyz)):
        return None
    ref = open(ref_sdf).read().split('\n')
    xyz = open(sim_xyz).read().split('\n')
    try:
        n = int(xyz[0].strip())
        coords = [ln.split()[1:4] for ln in xyz[2:2 + n]]
        na = int(ref[3][0:3])
    except (ValueError, IndexError):
        return None
    if na != len(coords):
        return None
    out = ref[:4]
    for i in range(na):
        try:
            x, y, z = (float(c) for c in coords[i])
        except (ValueError, IndexError):
            return None
        out.append("%10.4f%10.4f%10.4f" % (x, y, z) + ref[4 + i][30:])
    out += ref[4 + na:]
    open(sim_sdf, 'w').write('\n'.join(out))
    return calculate_rmsd_obrms(ref_sdf, sim_sdf, singularity_image)


def robust_backbone_rmsd(crystal_bb_path, sim_bb_path, output_dir, combo_number, singularity_image):
    """Faithful backbone RMSD via obrms (residue+name matching, terminal-truncation
    aware, shared connectivity). Returns float, or None if the peptides cannot be
    reconciled (obrms-only policy: no PyMOL fallback)."""
    cry_res = _bb_res_groups_crystal(_parse_pdb_atoms(crystal_bb_path))
    sim_res = _bb_res_groups_sim(_parse_pdb_atoms(sim_bb_path))
    if not cry_res or not sim_res:
        return None
    # Either side may be the shorter contiguous block (crystal missing terminal
    # residues, or MD peptide short of the crystal reference); align by sliding
    # the shorter over the longer and taking the best equal-length window.
    cry_sub, sim_sub = _bb_align_blocks(cry_res, sim_res)
    cry = [a for g in cry_sub for a in g]
    sim = [a for g in sim_sub for a in g]
    if len(cry) != len(sim):
        return None
    if any(cry[i]['name'] != sim[i]['name'] for i in range(len(cry))):
        return None
    # obrms with shared connectivity == the faithful in-place RMSD.
    direct = (sum(_dist2(cry[i], sim[i]) for i in range(len(cry))) / len(cry)) ** 0.5
    v = _obrms_shared_connectivity(cry, sim, output_dir, combo_number, singularity_image)
    if v is None or abs(v - direct) > 0.05:
        # obrms build hiccup: fall back to the identical in-place value (this is
        # exactly what obrms returns for shared connectivity), still no PyMOL.
        return round(direct, 4)
    return v


def calculate_rmsd_peptide_obrms_atoms(peptide_crystal_pdb, peptide_sim_pdb, output_peptides_dir, combo_number, singularity_image, atom_mode='heavy'):
    """
    Peptide RMSD on selected atoms: filter crystal/sim, match by distance, rewrite sim with crystal
    topology, run obrms; if inf, fall back to PyMOL align.
    atom_mode: 'heavy', 'backbone' (N, CA, C, O), or 'ca'
    obrms needs equal atom counts; crystal_ref uses selected atoms only.
    For backbone, use _backbone.pdb files from PyMOL when present.
    """
    # Backbone (default): robust residue+name matching + shared-connectivity obrms.
    # This never returns 'inf' and needs no PyMOL fallback (obrms-only policy).
    if atom_mode == 'backbone':
        crystal_bb_path = os.path.splitext(peptide_crystal_pdb)[0] + '_backbone.pdb'
        sim_bb_path = os.path.splitext(peptide_sim_pdb)[0] + '_backbone.pdb'
        if os.path.exists(crystal_bb_path) and os.path.exists(sim_bb_path):
            return robust_backbone_rmsd(
                crystal_bb_path, sim_bb_path, output_peptides_dir,
                combo_number, singularity_image)
        return None

    crystal_atoms = _parse_pdb_atoms(peptide_crystal_pdb)
    sim_atoms = _parse_pdb_atoms(peptide_sim_pdb)
    if not crystal_atoms or not sim_atoms:
        return None
    crystal_sel = _filter_atoms_by_mode(crystal_atoms, atom_mode)
    sim_sel = _filter_atoms_by_mode(sim_atoms, atom_mode)
    if not crystal_sel or not sim_sel:
        return None
    crystal_to_sim = _match_atoms_by_distance(crystal_sel, sim_sel)
    if None in crystal_to_sim:
        return None
    suffix = {'ca': 'ca', 'backbone': 'backbone', 'heavy': 'heavy'}.get(atom_mode, 'heavy')
    # Backbone: distinct filename from ca/heavy
    sim_as_crystal_path = os.path.join(
        output_peptides_dir,
        f"peptide_sim_{combo_number}_backbone_as_crystal.pdb" if atom_mode == 'backbone' else f"peptide_sim_{combo_number}_as_crystal.pdb"
    )

    _write_pdb(sim_as_crystal_path, crystal_sel, sim_sel, crystal_to_sim)

    # Backbone: PyMOL crystal_backbone + sim_as_crystal (obrms ~4 Å)
    if atom_mode == 'backbone':
        crystal_bb_path = os.path.splitext(peptide_crystal_pdb)[0] + '_backbone.pdb'
        if os.path.exists(crystal_bb_path):
            rmsd = calculate_rmsd_obrms(crystal_bb_path, sim_as_crystal_path, singularity_image)
            if rmsd is not None:
                return rmsd

    # Para ca/heavy o fallback backbone: crystal_ref nuestro
    crystal_ref_path = os.path.join(output_peptides_dir, f"peptide_crystal_{combo_number}_{suffix}_ref.pdb")
    _write_pdb(crystal_ref_path, crystal_sel, crystal_sel, list(range(len(crystal_sel))))
    rmsd = calculate_rmsd_obrms(crystal_ref_path, sim_as_crystal_path, singularity_image)
    if rmsd is not None:
        return rmsd
    # Fallback: obrms inf → PyMOL align
    if atom_mode == 'heavy':
        selector = "not elem H"
        desc = "heavy atoms"
    elif atom_mode == 'backbone':
        # "backbone" keyword fails for L01 ligand; use explicit atom names
        selector = "name N or name CA or name C or name O"
        desc = "backbone (N, CA, C, O)"
    else:
        selector = "name CA"
        desc = "CA"
    print(f"    (obrms inf → PyMOL align on {desc})")
    return _calculate_rmsd_pymol_selector(peptide_crystal_pdb, peptide_sim_pdb, singularity_image, selector)


def calculate_rmsd_obrms(reference_file, query_file, singularity_image="singularity/new_ms.simg"):
    """
    RMSD with obrms.

    Args:
        reference_file: Reference (aligned crystal)
        query_file: Query (simulation)
        singularity_image: Singularity image path

    Returns:
        float RMSD or None on error
    """
    cmd = f'singularity exec {singularity_image} obrms "{reference_file}" "{query_file}"'
    
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            check=True,
            timeout=30
        )
        
        rmsd_str = result.stdout.strip()

        rmsd_lower = rmsd_str.lower()
        if 'inf' in rmsd_lower or 'infinity' in rmsd_lower:
            print(f"    ⚠ obrms returned 'inf' — possible structure mismatch")
            return None

        try:
            rmsd_value = float(rmsd_str)
            if rmsd_value == float('inf') or rmsd_value == float('-inf'):
                print(f"    ⚠ RMSD is infinite")
                return None
            if rmsd_value < 0 or rmsd_value > 1000:
                print(f"    ⚠ RMSD out of reasonable range: {rmsd_value}")
                return None
            return rmsd_value
        except ValueError:
            for word in rmsd_str.split():
                try:
                    rmsd_value = float(word)
                    if rmsd_value == float('inf') or rmsd_value == float('-inf'):
                        continue
                    if rmsd_value >= 0 and rmsd_value <= 1000:
                        return rmsd_value
                except ValueError:
                    continue
            print(f"    ⚠ Could not parse RMSD from: {rmsd_str}")
            return None

    except subprocess.TimeoutExpired:
        print(f"    ⚠ obrms timeout (30s)")
        return None
    except subprocess.CalledProcessError as e:
        if e.stderr:
            print(f"  obrms error: {e.stderr.strip()}")
        return None
    except Exception as e:
        print(f"  Unexpected error: {e}")
        return None


def process_simulation(vs_folder, crystal_file, output_dir, combo_dir=None, singularity_image="singularity/new_ms.simg", dry_run=False, case_override=None, chain_override=None, peptide_rmsd_atoms="heavy"):
    """
    Process one MD run: extract peptides, align, compute RMSD.

    Args:
        vs_folder: VS_GR_* folder
        crystal_file: Full crystal PDB or None to search
        output_dir: Working directory for intermediates
        combo_dir: Combination folder (local crystal lookup)
        singularity_image: Singularity image path
        dry_run: Print actions only
        case_override: --case for crystal lookup
        chain_override: --chain peptide chain in crystal
        peptide_rmsd_atoms: 'backbone', 'heavy', or 'ca'

    Returns:
        tuple: (combo_number, rmsd) or (combo_number, None) on failure
    """
    vs_name = os.path.basename(vs_folder)
    
    # Combination id from folder name: VS_GR_*query_combination_{num}*
    import re
    match = re.search(r'query_combination_(\d+)', vs_name)
    if not match:
        return None, None
    
    combo_number = match.group(1)
    
    if dry_run:
        print(f"  [DRY RUN] {vs_name} (combination {combo_number})")
        return combo_number, None
    
    sim_file = find_simulation_peptide(vs_folder)
    if not sim_file:
        print(f"  ✗ {vs_name}: simulation PDB not found")
        return combo_number, None
    
    # Case id from folder e.g. VS_GR_1CJR_A_GN_query_combination_13 -> 1CJR_A
    case_match = re.search(r'VS_GR_([a-zA-Z0-9]+)_([A-Z])_', vs_name)
    case_name = None
    peptide_chain = None
    if case_match:
        case_base = case_match.group(1)  # 1CJR
        protein_chain = case_match.group(2)  # A
        case_name = f"{case_base}_{protein_chain}"
    # --case overrides folder-derived case name
    if case_override:
        case_name = case_override.strip().upper()
    
    if not crystal_file and combo_dir:
        crystal_base = "peptide_pdb_fragments"
        project_root = os.getcwd()
        crystal_file, peptide_chain = find_crystal_structure(crystal_base, combo_dir, project_root, case_name, preferred_peptide_chain=chain_override)
    elif not crystal_file:
        # Try complex/ then peptide_pdb_fragments/{case}
        project_root = os.getcwd()
        crystal_file, peptide_chain = find_crystal_structure("complex", None, project_root, case_name, preferred_peptide_chain=chain_override)
        if not crystal_file and case_name:
            crystal_file, peptide_chain = find_crystal_structure("peptide_pdb_fragments", None, project_root, case_name, preferred_peptide_chain=chain_override)
    
    if not crystal_file:
        print(f"  ✗ {vs_name}: crystal PDB not found")
        return combo_number, None
    
    # --chain overrides inferred peptide chain
    if chain_override is not None:
        peptide_chain = str(chain_override).strip().upper() or "A"
    elif not peptide_chain:
        peptide_chain = "A"
    
    peptides_dir = os.path.join(output_dir, "extracted_peptides")
    os.makedirs(peptides_dir, exist_ok=True)

    peptide_sim = os.path.join(peptides_dir, f"peptide_sim_{combo_number}.pdb")
    peptide_crystal = os.path.join(peptides_dir, f"peptide_crystal_{combo_number}.pdb")
    
    # Simulation receptor PDB
    protein_sim = os.path.join(peptides_dir, f"protein_sim_{combo_number}.pdb")

    desc = {'ca': 'CA atoms', 'backbone': 'backbone (N, CA, C, O)', 'heavy': 'heavy atoms'}.get(peptide_rmsd_atoms, peptide_rmsd_atoms)

    def _extract_and_rmsd(cryst, pchain):
        if not extract_peptide_with_pymol(sim_file, cryst, peptide_sim, peptide_crystal,
                                          protein_sim, singularity_image, pchain,
                                          force_peptide_chain=True):
            return None
        return calculate_rmsd_peptide_obrms_atoms(
            peptide_crystal, peptide_sim, peptides_dir, combo_number,
            singularity_image, atom_mode=peptide_rmsd_atoms)

    print(f"  sim:   {os.path.abspath(sim_file)}")

    # Multi-copy crystals: try every candidate copy and keep the lowest RMSD (the
    # docking only matches one lattice copy; others sit tens of A away and inflate
    # rmsd_md). Enabled whenever the peptide chain isn't pinned with --chain. Note:
    # combo_dir is only the per-combination folder used for topology lookup, NOT a
    # crystal source (the reference is always resolved from complex/), so it must
    # NOT gate this. find_crystal_candidates only returns >1 when complex/ actually
    # holds several lattice copies for this case (e.g. 3kut_D_A vs 3kut_D_B), so
    # single-copy cases keep the original single path untouched.
    candidates = None
    if chain_override is None:
        cand = find_crystal_candidates(os.getcwd(), case_name)
        if len(cand) > 1:
            candidates = cand

    if candidates:
        best = None
        for cryst in candidates:
            pchain = crystal_peptide_chain(cryst)
            r = _extract_and_rmsd(cryst, pchain)
            if r is not None and r >= 0 and r < 1e10 and (best is None or r < best[0]):
                best = (r, cryst, pchain)
        if best is None:
            print(f"  ✗ {vs_name}: RMSD computation failed (all {len(candidates)} crystal copies)")
            return combo_number, None
        rmsd, crystal_file, peptide_chain = best
        # Re-extract with the winning copy so the saved peptide files match it.
        _extract_and_rmsd(crystal_file, peptide_chain)
        print(f"  crystal: {os.path.abspath(crystal_file)} (best of {len(candidates)} copies, chain {peptide_chain})")
        print(f"    (obrms on {desc}; sim rewritten with crystal topology)")
    else:
        print(f"  crystal: {os.path.abspath(crystal_file)}")
        print(f"  [1/3] Extracting peptide and receptor with PyMOL... (peptide chain: {peptide_chain})")
        if not extract_peptide_with_pymol(sim_file, crystal_file, peptide_sim, peptide_crystal, protein_sim, singularity_image, peptide_chain, force_peptide_chain=(chain_override is not None)):
            print(f"  ✗ {vs_name}: peptide/receptor extraction failed")
            return combo_number, None
        print(f"  [2/3] Computing RMSD with obrms...")
        rmsd = calculate_rmsd_peptide_obrms_atoms(
            peptide_crystal, peptide_sim, peptides_dir, combo_number, singularity_image,
            atom_mode=peptide_rmsd_atoms
        )
        if rmsd is not None:
            print(f"    (obrms on {desc}; sim rewritten with crystal topology)")
    
    if rmsd is None:
        print(f"  ✗ {vs_name}: RMSD computation failed")
        print(f"    Peptide/receptor PDBs saved under: {peptides_dir}")
        return combo_number, None
    
    if rmsd == float('inf') or rmsd == float('-inf') or (isinstance(rmsd, float) and not (rmsd >= 0 and rmsd < 1e10)):
        print(f"  ⚠ {vs_name}: invalid RMSD (inf or huge): {rmsd}")
        print(f"    Peptide/receptor PDBs saved under: {peptides_dir}")
        return combo_number, None
    
    print(f"  [3/3] ✓ RMSD: {rmsd:.3f} Å")
    print(f"    Peptide/receptor PDBs saved under: {peptides_dir}")
    
    return combo_number, rmsd


def main():
    parser = argparse.ArgumentParser(
        description="Compute MD peptide RMSD vs crystal structure"
    )
    parser.add_argument(
        '--crystal-dir',
        default='peptide_pdb_fragments',
        help='Base directory for crystal inputs (default: peptide_pdb_fragments)'
    )
    parser.add_argument(
        '--output-csv',
        default='md_rmsd_results.csv',
        help='Output CSV path (default: md_rmsd_results.csv)'
    )
    parser.add_argument(
        '--singularity-image',
        default='singularity/new_ms.simg',
        help='Singularity image (default: singularity/new_ms.simg)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print planned actions only; do not run calculations'
    )
    parser.add_argument(
        '--project-root',
        default='.',
        help='Project root directory (default: .)'
    )
    parser.add_argument(
        '--output-dir',
        default='md_rmsd_peptides',
        help='Directory for extracted peptide/receptor PDBs (default: md_rmsd_peptides)'
    )
    parser.add_argument(
        '--case',
        default=None,
        metavar='CASE',
        help='Only process this case (e.g. 1AQD_A, 1AQD_C). Omit to process all VS_GR_* folders'
    )
    parser.add_argument(
        '--chain',
        default=None,
        metavar='CHAIN',
        help='Peptide chain id in crystal (e.g. A, B, C, H). If omitted, inferred from filename.'
    )
    parser.add_argument(
        '--peptide-rmsd-atoms',
        choices=['ca', 'backbone', 'heavy'],
        default='backbone',
        help="'backbone' = N, CA, C, O (default); 'heavy' = all heavy atoms; 'ca' = CA only"
    )

    args = parser.parse_args()
    if not os.path.exists(args.singularity_image):
        legacy_image = "new_ms.simg"
        if os.path.exists(legacy_image):
            print(f"ℹ Using legacy image path: {legacy_image}")
            args.singularity_image = legacy_image

    print("=" * 70)
    print("MD peptide RMSD vs crystal")
    print("=" * 70)
    print(f"Crystal directory: {args.crystal_dir}")
    print(f"Singularity image: {args.singularity_image}")
    print(f"Output CSV: {args.output_csv}")
    rmsd_desc = {'backbone': 'backbone (N, CA, C, O)', 'heavy': 'heavy atoms', 'ca': 'CA atoms'}.get(args.peptide_rmsd_atoms, 'backbone')
    print(f"RMSD mode: {rmsd_desc}")
    print()

    crystal_file = None
    print("Crystal PDB is resolved per combination from the case name.")
    print()

    vs_folders = find_vs_folders(args.project_root)
    if args.case:
        case_upper = args.case.strip().upper()
        case_base = case_upper.split("_")[0] if "_" in case_upper else case_upper
        def _case_base_from_folder(path):
            m = re.search(r'VS_GR_([a-zA-Z0-9]+)_', os.path.basename(path), re.IGNORECASE)
            return m.group(1).upper() if m else None
        vs_folders = [p for p in vs_folders if _case_base_from_folder(p) == case_base]
        print(f"Case filter: {case_upper} (base {case_base}) → {len(vs_folders)} runs")
    def _combo_sort_key(path):
        m = re.search(r'query_combination_(\d+)', os.path.basename(path))
        return (int(m.group(1)) if m else 0, path)
    vs_folders = sorted(vs_folders, key=_combo_sort_key)
    if not vs_folders:
        msg = "No VS_GR_* folders found" + (f" for case {args.case}" if args.case else "")
        print(f"Error: {msg} under {args.project_root}")
        sys.exit(1)

    print(f"Found {len(vs_folders)} MD runs")
    if args.dry_run:
        print("⚠ DRY RUN — no calculations will be executed")
    print()

    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "extracted_peptides"), exist_ok=True)
    peptides_out = os.path.join(output_dir, "extracted_peptides")
    print(f"Output directory: {output_dir}")
    print(f"  PDBs: {peptides_out}")
    print()
    
    results = []
    
    for i, vs_folder in enumerate(vs_folders, 1):
        vs_name = os.path.basename(vs_folder)
        print(f"[{i}/{len(vs_folders)}] {vs_name}")
        
        match = re.search(r'query_combination_(\d+)', vs_name)
        combo_dir = None
        if match:
            combo_num = match.group(1)
            case_upper = args.case.strip().upper() if args.case else None
            if case_upper:
                base_dirs = [
                    f"valid_GN_{case_upper}_final",
                    f"valid_LF_{case_upper}_final",
                    "valid_GN_final",
                    "valid_LF_final",
                ]
            else:
                base_dirs = ["valid_GN_final", "valid_LF_final"]
            for base_dir in base_dirs:
                potential_dir = os.path.join(base_dir, f"combination_{combo_num}")
                if os.path.exists(potential_dir):
                    combo_dir = potential_dir
                    break
        
        case_override = args.case.strip().upper() if args.case else None
        chain_override = args.chain.strip().upper() if args.chain else None
        combo_number, rmsd = process_simulation(
            vs_folder, crystal_file, output_dir, combo_dir, args.singularity_image, args.dry_run,
            case_override=case_override, chain_override=chain_override,
            peptide_rmsd_atoms=args.peptide_rmsd_atoms
        )
        
        if combo_number:
            results.append({
                'combination': combo_number,
                'vs_folder': vs_name,
                'rmsd': rmsd if rmsd is not None else 'NA'
            })
    
    print()
    print("=" * 70)
    print("Saving results...")

    with open(args.output_csv, 'w', newline='') as csvfile:
        fieldnames = ['combination', 'vs_folder', 'rmsd']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"✓ Results written to: {args.output_csv}")
    print(f"  Total rows: {len(results)}")
    print(f"  With valid RMSD: {sum(1 for r in results if r['rmsd'] != 'NA')}")
    print(f"  With NA: {sum(1 for r in results if r['rmsd'] == 'NA')}")
    print(f"  Extracted structures: {output_dir}")


if __name__ == "__main__":
    main()


