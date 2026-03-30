#!/usr/bin/env python3
"""
Compute score_only with gnina for peptide and protein structures extracted from MD.

For each combination:
1. Convert peptide and protein to PDBQT
2. Run gnina with --score_only
3. Parse the Affinity value
4. Write a CSV with results
"""

import os
import glob
import csv
import subprocess
import sys
import argparse
import re
from pathlib import Path


def _wrap_singularity(cmd, singularity_image, singularity_bind):
    """Prefix cmd with singularity exec --bind when an image is set."""
    if not singularity_image or not singularity_bind:
        return cmd
    return ["singularity", "exec", "--bind", singularity_bind, singularity_image] + cmd


def convert_ligand_to_pdbqt(peptide_pdb, output_pdbqt, mgltools_path="MetaScreener/external_sw/mgltools", singularity_image=None, singularity_bind=None, timeout=180):
    """
    Convert a peptide PDB to PDBQT using prepare_ligand4.py

    Args:
        peptide_pdb: Input peptide PDB path
        output_pdbqt: Output PDBQT path
        mgltools_path: Path to mgltools
        timeout: Max seconds for conversion (default 180; Singularity often needs more)

    Returns:
        bool: True on success, False otherwise
    """
    if singularity_image:
        mgltools_path = os.path.abspath(mgltools_path)
    prepare_ligand = os.path.join(mgltools_path, "bin/pythonsh")
    script_path = os.path.join(mgltools_path, "MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py")

    if not os.path.exists(prepare_ligand):
        print(f"  ✗ Error: not found: {prepare_ligand}")
        return False

    if not os.path.exists(script_path):
        print(f"  ✗ Error: not found: {script_path}")
        return False

    # With Singularity, use absolute paths so they resolve inside the container
    if singularity_image:
        peptide_pdb = os.path.abspath(peptide_pdb)
        output_pdbqt = os.path.abspath(output_pdbqt)
    cmd = [
        prepare_ligand,
        script_path,
        "-l", peptide_pdb,
        "-o", output_pdbqt,
        "-A", "bond_hydrogens"
    ]
    cmd = _wrap_singularity(cmd, singularity_image, singularity_bind)

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=os.getcwd()
        )

        if result.returncode == 0 and os.path.exists(output_pdbqt):
            return True
        else:
            out = (result.stdout or "").strip()
            err = (result.stderr or "").strip()
            if err:
                print(f"  ⚠ prepare_ligand4 error: {err}")
            if out and result.returncode != 0:
                for line in out.splitlines()[-3:]:
                    print(f"    {line}")
            return False
    except subprocess.TimeoutExpired:
        print(f"  ✗ Timeout converting ligand (limit {timeout}s)")
        return False
    except Exception as e:
        print(f"  ✗ Error converting ligand: {e}")
        return False


def convert_protein_to_pdbqt(protein_pdb, output_dir, conversion_script="conversion_targets/pdbqtconvert.sh", singularity_image=None, singularity_bind=None, timeout=300):
    """
    Convert a protein PDB to PDBQT using pdbqtconvert.sh

    Args:
        protein_pdb: Input protein PDB path
        output_dir: Directory for output (PDBQT uses same basename as PDB)
        conversion_script: Path to conversion script
        singularity_image: Optional Singularity image
        singularity_bind: Optional Singularity bind mount
        timeout: Max seconds (default 300)

    Returns:
        str: Path to generated PDBQT or None on failure
    """
    if not os.path.exists(conversion_script):
        print(f"  ✗ Error: not found: {conversion_script}")
        return None

    # With Singularity, use absolute paths
    if singularity_image:
        conversion_script = os.path.abspath(conversion_script)
        protein_pdb = os.path.abspath(protein_pdb)
    cmd = ["bash", conversion_script, protein_pdb]
    cmd = _wrap_singularity(cmd, singularity_image, singularity_bind)
    # pdbqtconvert.sh writes PDBQT next to the PDB
    pdbqt_file = protein_pdb.replace(".pdb", ".pdbqt")
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=os.getcwd()
        )

        if result.returncode == 0 and os.path.exists(pdbqt_file):
            return pdbqt_file
        else:
            err = (result.stderr or "").strip()
            out = (result.stdout or "").strip()
            if err:
                print(f"  ⚠ pdbqtconvert error: {err}")
            if out:
                for line in out.splitlines()[-5:]:
                    print(f"    {line}")
            return None
    except subprocess.TimeoutExpired:
        print(f"  ✗ Timeout converting protein (limit {timeout}s)")
        return None
    except Exception as e:
        print(f"  ✗ Error converting protein: {e}")
        return None


def calculate_score_only(protein_pdbqt, ligand_pdbqt, gnina_path="MetaScreener/external_sw/gnina/gnina", singularity_image=None, singularity_bind=None, timeout=180):
    """
    Run gnina --score_only.

    Args:
        protein_pdbqt: Protein PDBQT path
        ligand_pdbqt: Ligand PDBQT path
        gnina_path: Path to gnina executable
        singularity_image: Optional Singularity image
        singularity_bind: Optional Singularity bind mount
        timeout: Max seconds (default 180)

    Returns:
        float: Affinity in kcal/mol, or None on error
    """
    if not singularity_image and not os.path.exists(gnina_path):
        print(f"  ✗ Error: not found: {gnina_path}")
        return None

    # With Singularity, use absolute paths
    if singularity_image:
        protein_pdbqt = os.path.abspath(protein_pdbqt)
        ligand_pdbqt = os.path.abspath(ligand_pdbqt)
        gnina_path = os.path.abspath(gnina_path) if os.path.exists(gnina_path) else gnina_path
    cmd = [
        gnina_path,
        "-r", protein_pdbqt,
        "-l", ligand_pdbqt,
        "--score_only",
        "--cnn_scoring=none"  # Avoid OOM without GPU; gnina recommends this when no GPU
    ]
    cmd = _wrap_singularity(cmd, singularity_image, singularity_bind)

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=os.getcwd()
        )

        if result.returncode != 0:
            err = (result.stderr or "").strip()
            out = (result.stdout or "").strip()
            print(f"  ⚠ gnina exited with code {result.returncode}")
            print(f"  Command: {' '.join(cmd)}")
            if err:
                for line in err.splitlines()[-15:]:
                    print(f"    stderr: {line}")
            if out:
                for line in out.splitlines()[-15:]:
                    print(f"    stdout: {line}")
            return None

        # Look for "Affinity: X.XXXXX (kcal/mol)"
        output = result.stdout + result.stderr
        affinity_match = re.search(r'Affinity:\s+([-\d.]+)\s+\(kcal/mol\)', output)

        if affinity_match:
            try:
                affinity_value = float(affinity_match.group(1))
                return affinity_value
            except ValueError:
                print(f"  ⚠ Could not parse Affinity as float: {affinity_match.group(1)}")
                return None
        else:
            print(f"  ⚠ No 'Affinity' line found in gnina output")
            print(f"  Last lines of output:")
            for line in output.strip().splitlines()[-15:]:
                print(f"    {line}")
            return None

    except subprocess.TimeoutExpired:
        print(f"  ✗ Timeout running gnina (limit {timeout}s)")
        return None
    except Exception as e:
        print(f"  ✗ Error running gnina: {e}")
        return None


def process_combination(combo_number, peptides_dir, output_dir, mgltools_path, conversion_script, gnina_path, singularity_image=None, singularity_bind=None, timeout_ligand=180, timeout_protein=300, timeout_gnina=180):
    """
    Process one combination: convert inputs and compute score_only

    Args:
        combo_number: Combination index
        peptides_dir: Directory with extracted peptide/protein PDBs
        output_dir: Directory for temporary PDBQT files
        mgltools_path: Path to mgltools
        conversion_script: Protein conversion script path
        gnina_path: Path to gnina
        singularity_image: Optional Singularity image
        singularity_bind: Optional Singularity bind mount
        timeout_ligand: Max seconds for peptide conversion
        timeout_protein: Max seconds for protein conversion
        timeout_gnina: Max seconds for gnina score_only

    Returns:
        float: score_only value or None on error
    """
    peptide_pdb = os.path.join(peptides_dir, f"peptide_sim_{combo_number}.pdb")
    protein_pdb = os.path.join(peptides_dir, f"protein_sim_{combo_number}.pdb")

    if not os.path.exists(peptide_pdb):
        print(f"  ✗ Not found: {peptide_pdb}")
        return None

    if not os.path.exists(protein_pdb):
        print(f"  ✗ Not found: {protein_pdb}")
        return None

    os.makedirs(output_dir, exist_ok=True)

    peptide_pdbqt = os.path.join(output_dir, f"peptide_sim_{combo_number}.pdbqt")
    print(f"  [1/3] Converting peptide to PDBQT...")
    if not convert_ligand_to_pdbqt(peptide_pdb, peptide_pdbqt, mgltools_path, singularity_image, singularity_bind, timeout=timeout_ligand):
        return None

    print(f"  [2/3] Converting protein to PDBQT...")
    protein_pdbqt = convert_protein_to_pdbqt(protein_pdb, output_dir, conversion_script, singularity_image, singularity_bind, timeout=timeout_protein)
    if not protein_pdbqt:
        return None

    print(f"  [3/3] Running gnina score_only...")
    score = calculate_score_only(protein_pdbqt, peptide_pdbqt, gnina_path, singularity_image, singularity_bind, timeout=timeout_gnina)

    return score


def main():
    parser = argparse.ArgumentParser(
        description="Compute score_only with gnina for peptide and protein structures from MD"
    )
    parser.add_argument(
        '--peptides-dir',
        default='md_rmsd_peptides/extracted_peptides',
        help='Directory with extracted peptide/protein PDBs (default: md_rmsd_peptides/extracted_peptides)'
    )
    parser.add_argument(
        '--output-csv',
        default='score_only_results.csv',
        help='Output CSV path (default: score_only_results.csv)'
    )
    parser.add_argument(
        '--output-dir',
        default='md_rmsd_peptides/pdbqt_files',
        help='Directory for temporary PDBQT files (default: md_rmsd_peptides/pdbqt_files)'
    )
    parser.add_argument(
        '--mgltools-path',
        default='MetaScreener/external_sw/mgltools',
        help='Path to mgltools (default: MetaScreener/external_sw/mgltools)'
    )
    parser.add_argument(
        '--conversion-script',
        default='conversion_targets/pdbqtconvert.sh',
        help='Protein conversion script (default: conversion_targets/pdbqtconvert.sh)'
    )
    parser.add_argument(
        '--gnina-path',
        default='MetaScreener/external_sw/gnina/gnina',
        help='Path to gnina executable (default: MetaScreener/external_sw/gnina/gnina)'
    )
    parser.add_argument(
        '--singularity-image',
        default='singularity/metascreener_22.04.simg',
        help='Singularity image for mgltools/gnina (default: singularity/metascreener_22.04.simg). Empty or --no-singularity to disable.'
    )
    parser.add_argument(
        '--singularity-bind',
        default='/gpfs/',
        help='Bind path for singularity exec (default: /gpfs/)'
    )
    parser.add_argument(
        '--no-singularity',
        action='store_true',
        help='Do not use Singularity; run commands on the host'
    )
    parser.add_argument(
        '--timeout-ligand',
        type=int,
        default=180,
        help='Timeout (seconds) for peptide→PDBQT (default: 180; Singularity may need more)'
    )
    parser.add_argument(
        '--timeout-protein',
        type=int,
        default=300,
        help='Timeout (seconds) for protein→PDBQT (default: 300)'
    )
    parser.add_argument(
        '--timeout-gnina',
        type=int,
        default=180,
        help='Timeout (seconds) for gnina --score_only (default: 180)'
    )

    args = parser.parse_args()

    singularity_image = None if args.no_singularity else args.singularity_image
    singularity_bind = args.singularity_bind if singularity_image else None
    if singularity_image:
        singularity_image = os.path.abspath(singularity_image)
        if not os.path.exists(singularity_image):
            print(f"Warning: not found {singularity_image}; running without Singularity")
            singularity_image = None
            singularity_bind = None

    print("=" * 70)
    print("score_only calculation with gnina")
    print("=" * 70)
    if singularity_image:
        print(f"Singularity: {singularity_image} (--bind {singularity_bind})")
    print(f"Peptides directory: {args.peptides_dir}")
    print(f"PDBQT output directory: {args.output_dir}")
    print(f"Output CSV: {args.output_csv}")
    print(f"Timeouts: ligand={args.timeout_ligand}s, protein={args.timeout_protein}s, gnina={args.timeout_gnina}s")
    print()

    if not os.path.exists(args.peptides_dir):
        print(f"Error: directory does not exist: {args.peptides_dir}")
        sys.exit(1)

    peptide_files = glob.glob(os.path.join(args.peptides_dir, "peptide_sim_*.pdb"))

    if not peptide_files:
        print(f"No peptide_sim_*.pdb files in {args.peptides_dir}; writing empty CSV.")
        with open(args.output_csv, 'w', newline='') as csvfile:
            csvfile.write("combination_id,score_only\n")
        print(f"✓ Created {args.output_csv} (no rows).")
        sys.exit(0)

    combo_numbers = []
    for peptide_file in peptide_files:
        filename = os.path.basename(peptide_file)
        match = re.search(r'peptide_sim_(\d+)\.pdb', filename)
        if match:
            combo_numbers.append(match.group(1))

    combo_numbers = sorted(combo_numbers, key=lambda x: int(x))

    print(f"Found {len(combo_numbers)} combinations")
    print()

    results = []

    for i, combo_number in enumerate(combo_numbers, 1):
        print(f"[{i}/{len(combo_numbers)}] Combination {combo_number}")

        score = process_combination(
            combo_number,
            args.peptides_dir,
            args.output_dir,
            args.mgltools_path,
            args.conversion_script,
            args.gnina_path,
            singularity_image,
            singularity_bind,
            timeout_ligand=args.timeout_ligand,
            timeout_protein=args.timeout_protein,
            timeout_gnina=args.timeout_gnina
        )

        if score is not None:
            print(f"  ✓ Score_only: {score:.5f} kcal/mol")
            results.append({
                'combination_id': combo_number,
                'score_only': score
            })
        else:
            print(f"  ✗ Error computing score_only")
            results.append({
                'combination_id': combo_number,
                'score_only': 'NA'
            })
        print()

    print("=" * 70)
    print("Saving results...")

    with open(args.output_csv, 'w', newline='') as csvfile:
        fieldnames = ['combination_id', 'score_only']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"✓ Results saved to: {args.output_csv}")
    print(f"  Total processed: {len(results)}")
    print(f"  With valid score_only: {sum(1 for r in results if r['score_only'] != 'NA')}")
    print(f"  With error (NA): {sum(1 for r in results if r['score_only'] == 'NA')}")


if __name__ == "__main__":
    main()
