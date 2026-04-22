#!/usr/bin/env python3
"""
Script para organizar todos los archivos generados durante el workflow en una carpeta.

Mueve:
- Carpetas VS_GN_* y VS_LF_*
- Carpetas valid_combinations_*
- Carpetas valid_GN_*_final y valid_LF_*_final
- Carpetas VS_GR_* (simulaciones MD)
- Archivos CSV generados
- Carpetas md_rmsd_peptides
- Otros archivos relacionados
"""

import os
import glob
import shutil
import sys
import argparse
from pathlib import Path


def find_and_move_files(case_name, target_dir, dry_run=False):
    """
    Encuentra y mueve todos los archivos relacionados con el caso.
    
    Args:
        case_name: Nombre del caso (ej: "1CJR_A")
        target_dir: Directorio destino
        dry_run: Si es True, solo muestra lo que haría sin mover
    
    Returns:
        dict: Estadísticas de archivos movidos
    """
    case_upper = case_name.upper()
    case_lower = case_name.lower()
    
    stats = {
        'folders_moved': 0,
        'files_moved': 0,
        'total_size': 0
    }
    
    # Crear directorio destino si no existe
    if not dry_run:
        os.makedirs(target_dir, exist_ok=True)
    
    # Patrones de archivos y carpetas a mover
    patterns_to_move = []
    
    # 1. Carpetas VS_GN_* y VS_LF_* del caso
    patterns_to_move.extend([
        f"VS_GN_{case_upper}_*",
        f"VS_GN_{case_lower}_*",
        f"VS_LF_{case_upper}_*",
        f"VS_LF_{case_lower}_*",
    ])
    
    # 2. Carpetas valid_combinations_*
    patterns_to_move.extend([
        f"valid_combinations_GN_{case_upper}",
        f"valid_combinations_LF_{case_upper}",
        f"valid_combinations_GN_{case_lower}",
        f"valid_combinations_LF_{case_lower}",
    ])
    
    # 3. Carpetas valid_GN_*_final y valid_LF_*_final
    patterns_to_move.extend([
        f"valid_GN_{case_upper}_final",
        f"valid_LF_{case_upper}_final",
        f"valid_GN_{case_lower}_final",
        f"valid_LF_{case_lower}_final",
    ])
    
    # 4. Carpeta {case}_GN (ej: 1BXP_B_GN) que puede contener VS_GN_*_Frag*
    patterns_to_move.extend([
        f"{case_upper}_GN",
        f"{case_lower}_GN",
    ])
    
    # 5. Carpetas VS_GR_* del caso (simulaciones MD)
    # Buscar con diferentes variaciones del nombre del caso
    patterns_to_move.extend([
        f"VS_GR_{case_upper}_*",
        f"VS_GR_{case_lower}_*",
        f"VS_GR_*{case_upper}*",
        f"VS_GR_*{case_lower}*",
    ])
    
    # 6. Archivos CSV relacionados (incluye nombres por caso para ejecución paralela)
    csv_patterns = [
        f"valid_fragment_combinations_GN_no_overlap.csv",
        f"valid_fragment_combinations_GN_no_overlap_{case_upper}.csv",
        f"valid_fragment_combinations_GN_no_overlap_{case_lower}.csv",
        f"valid_fragment_combinations_LF_no_overlap.csv",
        f"valid_fragment_combinations_LF_no_overlap_{case_upper}.csv",
        f"valid_fragment_combinations_LF_no_overlap_{case_lower}.csv",
        f"valid_fragment_combinations.csv",
        f"resultados_rmsd_md_{case_upper}.csv",
        f"resultados_rmsd_md_{case_lower}.csv",
        f"score_only_results.csv",
        f"score_only_results_{case_upper}.csv",
        f"score_only_results_{case_lower}.csv",
        f"all_metrics_GN_{case_upper}.csv",
        f"all_metrics_LF_{case_upper}.csv",
        f"all_metrics_GN_{case_lower}.csv",
        f"all_metrics_LF_{case_lower}.csv",
        f"mmpbsa_results_{case_upper}.csv",
        f"mmpbsa_results_{case_lower}.csv",
        f"fragment_energies_{case_upper}.csv",
        f"fragment_energies_{case_lower}.csv",
    ]
    
    # 7. Carpetas de resultados MD (incluye nombres por caso para ejecución paralela)
    patterns_to_move.extend([
        "md_rmsd_peptides",
        f"md_rmsd_peptides_{case_upper}",
        f"md_rmsd_peptides_{case_lower}",
    ])
    
    print("=" * 70)
    print(f"Organizando archivos del caso: {case_upper}")
    print("=" * 70)
    print(f"Directorio destino: {target_dir}")
    if dry_run:
        print("⚠ MODO DRY RUN - No se moverán archivos")
    print()
    
    # Mover carpetas
    print("Buscando carpetas...")
    folders_found = []
    for pattern in patterns_to_move:
        matches = glob.glob(pattern)
        for match in matches:
            if os.path.isdir(match) and match not in folders_found:
                folders_found.append(match)
    
    # Mover archivos CSV
    print("Buscando archivos CSV...")
    csv_files_found = []
    for pattern in csv_patterns:
        if os.path.exists(pattern) and os.path.isfile(pattern):
            csv_files_found.append(pattern)
    
    # Mostrar lo que se va a mover
    print(f"\nCarpetas encontradas: {len(folders_found)}")
    for folder in sorted(folders_found):
        size = get_dir_size(folder)
        stats['total_size'] += size
        print(f"  - {folder} ({format_size(size)})")
    
    print(f"\nArchivos CSV encontrados: {len(csv_files_found)}")
    for csv_file in sorted(csv_files_found):
        size = os.path.getsize(csv_file)
        stats['total_size'] += size
        print(f"  - {csv_file} ({format_size(size)})")
    
    if dry_run:
        print("\n⚠ DRY RUN: No se moverán archivos")
        return stats
    
    # Mover carpetas
    print(f"\nMoviendo {len(folders_found)} carpetas...")
    for folder in sorted(folders_found):
        try:
            dest = os.path.join(target_dir, os.path.basename(folder))
            if os.path.exists(dest):
                print(f"  ⚠ {folder} -> Ya existe en destino, omitiendo")
            else:
                shutil.move(folder, dest)
                print(f"  ✓ {folder} -> {dest}")
                stats['folders_moved'] += 1
        except Exception as e:
            print(f"  ✗ Error moviendo {folder}: {e}")
    
    # Mover archivos CSV
    print(f"\nMoviendo {len(csv_files_found)} archivos CSV...")
    csv_dest_dir = os.path.join(target_dir, "csv_files")
    os.makedirs(csv_dest_dir, exist_ok=True)
    
    for csv_file in sorted(csv_files_found):
        try:
            dest = os.path.join(csv_dest_dir, os.path.basename(csv_file))
            if os.path.exists(dest):
                print(f"  ⚠ {csv_file} -> Ya existe en destino, omitiendo")
            else:
                shutil.move(csv_file, dest)
                print(f"  ✓ {csv_file} -> {dest}")
                stats['files_moved'] += 1
        except Exception as e:
            print(f"  ✗ Error moviendo {csv_file}: {e}")
    
    return stats


def get_dir_size(path):
    """Calcula el tamaño total de un directorio."""
    total = 0
    try:
        for dirpath, dirnames, filenames in os.walk(path):
            for filename in filenames:
                filepath = os.path.join(dirpath, filename)
                if os.path.exists(filepath):
                    total += os.path.getsize(filepath)
    except Exception:
        pass
    return total


def format_size(size_bytes):
    """Formatea el tamaño en bytes a formato legible."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.2f} TB"


def main():
    parser = argparse.ArgumentParser(
        description="Organizar todos los archivos generados durante el workflow en una carpeta"
    )
    parser.add_argument(
        'case_name',
        help='Nombre del caso (ej: 1CJR_A)'
    )
    parser.add_argument(
        '--target-dir',
        help='Directorio destino (default: {case_name}_results)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Solo mostrar lo que haría sin mover archivos'
    )
    
    args = parser.parse_args()
    
    case_name = args.case_name
    if args.target_dir:
        target_dir = args.target_dir
    else:
        target_dir = f"{case_name}_results"
    
    stats = find_and_move_files(case_name, target_dir, args.dry_run)
    
    if not args.dry_run:
        print("\n" + "=" * 70)
        print("Resumen:")
        print("=" * 70)
        print(f"  Carpetas movidas: {stats['folders_moved']}")
        print(f"  Archivos movidos: {stats['files_moved']}")
        print(f"  Tamaño total: {format_size(stats['total_size'])}")
        print(f"  Archivos organizados en: {target_dir}")
        print("=" * 70)


if __name__ == "__main__":
    main()

