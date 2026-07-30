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


def find_case_vs_folders(case_name):
    """VS_GR_* (MD) folders that belong to this case.

    VS_GR folders are named by the RECEPTOR id (VS_GR_<pdbid>_<recchain>_...), not by
    the case's peptide chain, so a case-name glob (e.g. VS_GR_*4P3W_G*) never matches
    and the folders pile up at the project root. Identify them precisely from the
    vs_folder column of resultados_rmsd_md_<CASE>.csv (the exact runs this case used);
    fall back to a pdbid glob, which mirrors how calculate_md_rmsd selects runs for a
    case. Note: cases that share a pdbid+peptide (e.g. 6G0Y_I/6G0Y_J) share the same
    VS_GR folders on disk; the first case to organize moves them and the rest skip.
    """
    case_upper = case_name.upper()
    case_lower = case_name.lower()
    pdbid = case_upper.split("_")[0]
    found = []

    for csvname in (f"resultados_rmsd_md_{case_upper}.csv",
                    f"resultados_rmsd_md_{case_lower}.csv"):
        if not os.path.isfile(csvname):
            continue
        try:
            with open(csvname) as fh:
                header = fh.readline().rstrip("\n").split(",")
                idx = header.index("vs_folder") if "vs_folder" in header else None
                if idx is not None:
                    for line in fh:
                        cols = line.rstrip("\n").split(",")
                        if len(cols) > idx:
                            d = cols[idx].strip()
                            if d and os.path.isdir(d) and d not in found:
                                found.append(d)
        except OSError:
            pass
        break

    if found:
        return found

    for pat in (f"VS_GR_{pdbid}_*", f"VS_GR_{pdbid.lower()}_*"):
        for d in glob.glob(pat):
            if os.path.isdir(d) and d not in found:
                found.append(d)
    return found


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
    
    # 5. Carpetas VS_GR_* del caso (simulaciones MD): se mueven aparte, a una
    #    subcarpeta md_simulations/, y se identifican por sus nombres reales (no por
    #    el nombre del caso), porque VS_GR_* se nombra por el RECEPTOR
    #    (VS_GR_<pdbid>_<recchain>_...), no por la cadena de péptido del caso. Ver
    #    find_case_vs_folders() y la sección de movido más abajo.

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

    # Carpetas de simulación MD (VS_GR_*): identificarlas AHORA, antes de mover el
    # resultados_rmsd_md_<CASE>.csv del que dependemos para nombrarlas con precisión.
    print("Buscando carpetas de simulación MD (VS_GR_*)...")
    vs_folders_found = find_case_vs_folders(case_name)

    # Mostrar lo que se va a mover
    print(f"\nCarpetas encontradas: {len(folders_found)}")
    for folder in sorted(folders_found):
        size = get_dir_size(folder)
        stats['total_size'] += size
        print(f"  - {folder} ({format_size(size)})")

    print(f"\nCarpetas VS_GR_* (MD) encontradas: {len(vs_folders_found)}")
    for folder in sorted(vs_folders_found):
        size = get_dir_size(folder)
        stats['total_size'] += size
        print(f"  - {folder} ({format_size(size)}) -> {target_dir}/md_simulations/")

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

    # Mover carpetas de simulación MD (VS_GR_*) a la subcarpeta md_simulations/.
    # Mantiene la raíz del proyecto despejada y agrupa todo el caso en <CASE>_results/.
    # calculate_md_rmsd.find_vs_folders() también busca aquí, así que recomputar el
    # paso 13 tras organizar sigue encontrando las simulaciones.
    if vs_folders_found:
        md_dest_dir = os.path.join(target_dir, "md_simulations")
        os.makedirs(md_dest_dir, exist_ok=True)
        print(f"\nMoviendo {len(vs_folders_found)} carpetas VS_GR_* a {md_dest_dir}...")
        for folder in sorted(vs_folders_found):
            try:
                dest = os.path.join(md_dest_dir, os.path.basename(folder))
                if os.path.exists(dest):
                    print(f"  ⚠ {folder} -> Ya existe en destino, omitiendo")
                elif not os.path.isdir(folder):
                    print(f"  ⚠ {folder} -> ya no existe (movido por otro caso), omitiendo")
                else:
                    shutil.move(folder, dest)
                    print(f"  ✓ {folder} -> {dest}")
                    stats['folders_moved'] += 1
            except Exception as e:
                print(f"  ✗ Error moviendo {folder}: {e}")

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

