#!/usr/bin/env python3
"""
Convert .pdbqt files to .mol2 under valid combination folders.

Uses singularity and MetaScreener's convert_to.py for all .pdbqt files in combination_*.
"""

import os
import glob
import subprocess
import sys
import argparse


def convert_pdbqt_to_mol2_batch(pdbqt_files, mol2_files, singularity_image="singularity/metascreener.simg", dry_run=False):
    """
    Convert multiple pdbqt to mol2 in one singularity session.

    Args:
        pdbqt_files: List of .pdbqt paths.
        mol2_files: List of output .mol2 paths (same order as pdbqt_files).
        singularity_image: Singularity image path.
        dry_run: If True, only print commands.

    Returns:
        (total_processed, total_errors)
    """
    if not pdbqt_files:
        return 0, 0
    
    # Verificar que la imagen de singularity existe
    if not os.path.exists(singularity_image):
        print(f"  ✗ Error: No se encuentra la imagen de singularity: {singularity_image}")
        return 0, len(pdbqt_files)
    
    # Skip already converted
    files_to_convert = []
    for pdbqt_file, mol2_file in zip(pdbqt_files, mol2_files):
        pdbqt_abs = os.path.abspath(pdbqt_file)
        mol2_abs = os.path.abspath(mol2_file)
        
        if not os.path.exists(pdbqt_abs):
            print(f"  ⚠ Archivo no encontrado: {pdbqt_abs}")
            continue
        
        if os.path.exists(mol2_abs) and not dry_run:
            continue  # skip existing
        
        files_to_convert.append((pdbqt_abs, mol2_abs))
    
    if not files_to_convert:
        return 0, 0
    
    convert_script = "MetaScreener/extra_metascreener/used_by_metascreener/convert_to.py"
    
    if dry_run:
        print(f"  [DRY RUN] Would convert {len(files_to_convert)} files in one singularity session")
        for pdbqt_abs, mol2_abs in files_to_convert[:5]:
            print(f"    {os.path.basename(pdbqt_abs)} -> {os.path.basename(mol2_abs)}")
        if len(files_to_convert) > 5:
            print(f"    ... and {len(files_to_convert) - 5} more")
        return len(files_to_convert), 0
    
    # Bash script chaining all conversions
    bash_commands = []
    for pdbqt_abs, mol2_abs in files_to_convert:
        bash_commands.append(f'python {convert_script} "{pdbqt_abs}" "{mol2_abs}"')
    
    # Join with &&
    bash_script = ' && '.join(bash_commands)
    
    # Ejecutar todo en una sola sesión de singularity
    cmd = f'singularity exec {singularity_image} bash -c \'{bash_script}\''
    
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            check=False
        )
        
        # Count outputs
        successful = 0
        failed = 0
        
        for pdbqt_abs, mol2_abs in files_to_convert:
            if os.path.exists(mol2_abs):
                successful += 1
            else:
                failed += 1
                print(f"  ✗ Error: not created: {os.path.basename(mol2_abs)}")
        
        if result.returncode != 0 and result.stderr:
            print(f"  ⚠ Algunos errores durante la conversión:")
            print(f"    {result.stderr[:500]}")  # Mostrar primeros 500 caracteres
        
        return successful, failed
    
    except Exception as e:
        print(f"  ✗ Exception during batch conversion: {e}")
        return 0, len(files_to_convert)


def process_combinations_directory(base_dir, prefix_type="GN", singularity_image="singularity/metascreener.simg", dry_run=False, batch_size=100):
    """
    Process all combination folders under base_dir.

    Args:
        base_dir: Root with combination_* folders.
        prefix_type: 'GN' or 'LF' (unused except for logging context).
        singularity_image: Singularity image path.
        dry_run: If True, only print actions.
        batch_size: Files per batch (default: 100).

    Returns:
        (total_processed, total_errors)
    """
    combo_dirs = sorted(glob.glob(os.path.join(base_dir, "combination_*")))
    
    if not combo_dirs:
        print(f"⚠ No combination_* folders in {base_dir}")
        return 0, 0
    
    print(f"Found {len(combo_dirs)} combination folders")
    if dry_run:
        print("⚠ DRY RUN — no files will be modified")
    print()
    
    # Collect all conversions
    all_pdbqt_files = []
    all_mol2_files = []
    
    for combo_dir in combo_dirs:
        pdbqt_files = glob.glob(os.path.join(combo_dir, "*.pdbqt"))
        for pdbqt_file in pdbqt_files:
            mol2_file = pdbqt_file.replace(".pdbqt", ".mol2")
            
            if os.path.exists(mol2_file) and not dry_run:
                continue  # skip existing .mol2
            
            all_pdbqt_files.append(pdbqt_file)
            all_mol2_files.append(mol2_file)
    
    if not all_pdbqt_files:
        print("✓ All .mol2 files already exist")
        return 0, 0
    
    print(f"Total de archivos a convertir: {len(all_pdbqt_files)}")
    print()
    
    total_processed = 0
    total_errors = 0
    
    # Batches to avoid overly long shell commands
    for batch_start in range(0, len(all_pdbqt_files), batch_size):
        batch_end = min(batch_start + batch_size, len(all_pdbqt_files))
        batch_pdbqt = all_pdbqt_files[batch_start:batch_end]
        batch_mol2 = all_mol2_files[batch_start:batch_end]
        
        batch_num = (batch_start // batch_size) + 1
        total_batches = (len(all_pdbqt_files) + batch_size - 1) // batch_size
        
        print(f"Batch {batch_num}/{total_batches} ({len(batch_pdbqt)} files)...")
        
        processed, errors = convert_pdbqt_to_mol2_batch(
            batch_pdbqt, batch_mol2, singularity_image, dry_run
        )
        
        total_processed += processed
        total_errors += errors
        
        if not dry_run:
            print(f"  ✓ Lote {batch_num} completado: {processed} exitosos, {errors} errores")
    
    return total_processed, total_errors


def main():
    parser = argparse.ArgumentParser(
        description="Convert .pdbqt to .mol2 under valid combination folders"
    )
    parser.add_argument(
        'type',
        choices=['GN', 'LF'],
        help='Combination type: GN or LF'
    )
    parser.add_argument(
        '--base-dir',
        help='Directorio base con las combinaciones (default: valid_combinations_{type}/valid_no_overlap)'
    )
    parser.add_argument(
        '--singularity-image',
        default='singularity/metascreener.simg',
        help='Singularity image path (default: singularity/metascreener.simg)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Solo mostrar lo que haría sin ejecutar'
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        default=100,
        help='Files per singularity batch (default: 100)'
    )
    
    args = parser.parse_args()
    
    # Establecer directorio base por defecto
    if args.base_dir is None:
        args.base_dir = f"valid_combinations_{args.type}/valid_no_overlap"
    
    if not os.path.exists(args.base_dir):
        print(f"✗ Error: directory does not exist: {args.base_dir}")
        return 1
    
    print("=" * 70)
    print(f".pdbqt to .mol2 conversion for {args.type}")
    print("=" * 70)
    print(f"Directorio base: {args.base_dir}")
    print(f"Imagen singularity: {args.singularity_image}")
    print()
    
    total_processed, total_errors = process_combinations_directory(
        args.base_dir,
        args.type,
        args.singularity_image,
        args.dry_run,
        args.batch_size
    )
    
    print("\n" + "=" * 70)
    print(f"Summary: {total_processed} files processed, {total_errors} errors")
    print("=" * 70)
    
    if total_errors > 0:
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

