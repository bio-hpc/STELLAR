#!/usr/bin/env python3
"""
Compute MM/PBSA-style energies (van der Waals, electrostatic, total) from MD
trajectories and write one combined CSV.

Per combination folder:
1. Run grompp (in singularity/gr.simg) to build the .tpr
2. Run g_mmpbsa (local binary) to produce the MM .xvg
3. Parse the .xvg and average columns:
   - Protein-L01 VdW Energy
   - Protein-L01 Elec. Energy
   - Protein-L01 Total Energy

Output CSV columns:
combination,vs_folder,type,mmpbsa_vdw,mmpbsa_elec,mmpbsa_total,frames
"""

import argparse
import glob
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional


def run_cmd(cmd: List[str], cwd: Optional[str] = None, stdin_input: Optional[str] = None):
    """Run a command and return stdout; raise on failure."""
    print(" ".join(cmd))
    if stdin_input:
        print(f"  stdin: {stdin_input.strip()}")
    res = subprocess.run(
        cmd, 
        cwd=cwd, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, 
        text=True,
        input=stdin_input
    )
    if res.returncode != 0:
        print(res.stdout)
        print(res.stderr, file=sys.stderr)
        raise RuntimeError(f"Command failed with exit code {res.returncode}")
    return res.stdout


def find_first(pattern: str) -> Optional[str]:
    """First sorted glob match or None."""
    matches = sorted(glob.glob(pattern))
    return matches[0] if matches else None


def folder_has_mmpbsa_inputs(abs_folder: str) -> bool:
    """Return True if a VS_GR folder has the minimum inputs for MM/PBSA."""
    mdp = find_first(os.path.join(abs_folder, "grids", "*_md.mdp"))
    gro = find_first(os.path.join(abs_folder, "molecules", "*_complex_md.gro")) or find_first(
        os.path.join(abs_folder, "molecules", "*_complex*.gro")
    )
    top = find_first(os.path.join(abs_folder, "molecules", "*.top"))
    xtc = find_first(os.path.join(abs_folder, "molecules", "*_complex_md.xtc"))
    return all([mdp, gro, top, xtc])


def parse_mm_xvg(xvg_file: str) -> Dict[str, Optional[float]]:
    """
    Parse g_mmpbsa MM .xvg and average:
    s4: Protein-L01 VdW Energy
    s5: Protein-L01 Elec. Energy
    s6: Protein-L01 Total Energy
    """
    vdw_vals, elec_vals, total_vals = [], [], []
    with open(xvg_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("@"):
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            try:
                vdw_vals.append(float(parts[5]))   # s4
                elec_vals.append(float(parts[6]))  # s5
                total_vals.append(float(parts[7])) # s6
            except ValueError:
                continue
    frames = len(vdw_vals)
    def avg(lst):
        return round(sum(lst) / len(lst), 6) if lst else ""
    return {
        "vdw": avg(vdw_vals),
        "elec": avg(elec_vals),
        "total": avg(total_vals),
        "frames": frames,
    }


def extract_combo_id(folder: str) -> Optional[str]:
    """Parse combination index from folder name."""
    m = re.search(r"combination_(\d+)", folder)
    return m.group(1) if m else None


def detect_type(folder: str) -> str:
    """Infer GN/LF from folder name."""
    lower = folder.lower()
    if "_lf_" in lower or "lf_" in lower:
        return "LF"
    if "_gn_" in lower or "gn_" in lower:
        return "GN"
    return ""


# Index group names that are ions/solvent — not the L01 ligand
_NDX_ION_SOLVENT_NAMES = frozenset(
    {"NA", "CL", "K", "MG", "CA", "SOL", "WATER", "ION", "WATER_AND_IONS"}
)


def detect_ligand_group(ndx_file: str) -> Optional[int]:
    """
    Detect L01 ligand group index from an .ndx file.
    Returns 0-based index (g_mmpbsa: Group 0 = System, Group 13 = L01).
    Priority:
    1. Group name exactly "L01".
    2. Name contains L01 and is not ion/solvent (NA, CL, etc.).
    """
    if not Path(ndx_file).exists():
        return None

    def is_ion_or_solvent(name: str) -> bool:
        name_upper = name.upper().strip()
        if "WATER" in name_upper or "SOL" in name_upper or "ION" in name_upper:
            return True
        # Short ion names (NA, CL, ...)
        if name_upper in _NDX_ION_SOLVENT_NAMES:
            return True
        # Nombre que empieza por ion (NA_, NA , CL_, etc.) para no confundir "NA L01" con ligando
        first_token = name_upper.split()[0].split("_")[0] if name_upper else ""
        if first_token in _NDX_ION_SOLVENT_NAMES:
            return True
        # Skip Water_and_ions_L01, SOL_L01, etc.
        if "L01" in name_upper and ("WATER" in name_upper or "ION" in name_upper or "SOL" in name_upper):
            return True
        return False

    try:
        group_number = 0
        exact_l01_group = None  # grupo con nombre exacto "L01"
        first_l01_candidate = None  # primer grupo que contiene L01 y no es ion/solvente

        with open(ndx_file, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line.startswith("[") or not line.endswith("]"):
                    continue
                group_name = line[1:-1].strip()
                group_number += 1

                if is_ion_or_solvent(group_name):
                    continue

                # Exact "L01" ligand; g_mmpbsa uses 0-based indices
                if group_name.upper() == "L01":
                    exact_l01_group = group_number - 1
                    return exact_l01_group

                # Contiene L01 como palabra
                if "L01" in group_name.upper():
                    if first_l01_candidate is None:
                        first_l01_candidate = group_number - 1

        return first_l01_candidate
    except Exception as e:
        print(f"⚠ Error reading {ndx_file}: {e}")

    return None


def process_folder(folder: str, gr_simg: str, g_mmpbsa_bin: str, maxwarn: int, force: bool, base_dir: str = ".") -> Dict[str, str]:
    """
    Procesa una carpeta VS_GR_*:
    - Ejecuta grompp -> .tpr
    - Detecta grupos de proteína (1) y ligando (L01)
    - Ejecuta g_mmpbsa -> .xvg
    - Parsea y devuelve métricas
    """
    abs_folder = os.path.abspath(folder)
    combo_id = extract_combo_id(folder)
    combo_type = detect_type(folder)
    if not combo_id:
        raise RuntimeError(f"Could not parse combination_id from {folder}")

    # Main inputs
    mdp = find_first(os.path.join(abs_folder, "grids", "*_md.mdp"))
    # .gro may use _npt_8000 etc. instead of _md
    gro = find_first(os.path.join(abs_folder, "molecules", "*_complex_md.gro"))
    if not gro:
        gro = find_first(os.path.join(abs_folder, "molecules", "*_complex*.gro"))
    top = find_first(os.path.join(abs_folder, "molecules", "*.top"))
    xtc = find_first(os.path.join(abs_folder, "molecules", "*_complex_md.xtc"))
    ndx = find_first(os.path.join(abs_folder, "grids", "*_index.ndx")) or find_first(os.path.join(abs_folder, "soft_out", "*_index.ndx"))

    if not all([mdp, gro, top, xtc]):
        raise RuntimeError(f"Missing inputs for {folder}: mdp={mdp}, gro={gro}, top={top}, xtc={xtc}")

    # Base name for .tpr / mm_xvg from .gro
    gro_base = os.path.splitext(os.path.basename(gro))[0]
    # Normalize stem toward _complex_md
    if "_complex_md" in gro_base:
        base_no_ext = gro_base.replace("_complex_md", "_complex_md")
    elif "_complex" in gro_base:
        # Strip before _complex and add _complex_md
        base_no_ext = gro_base.split("_complex")[0] + "_complex_md"
    else:
        base_no_ext = gro_base
    
    tpr = os.path.join(abs_folder, "molecules", f"{base_no_ext}_pre.tpr")
    mm_xvg = os.path.join(abs_folder, "molecules", f"{base_no_ext}_mm_energy.xvg")

    # L01 ligand group
    protein_group = 1  # protein
    ligand_group = None
    
    if ndx:
        ligand_group = detect_ligand_group(ndx)
        if ligand_group:
            print(f"  L01 ligand group: {ligand_group}")
        else:
            print(f"  ⚠ Could not detect L01 group in {ndx}, using default group 13")
            ligand_group = 13
    else:
        print(f"  ⚠ No .ndx file; using default group 13 for L01")
        ligand_group = 13

    # 1) grompp (unless exists or --force)
    if force or not Path(tpr).exists():
        run_cmd([
            "singularity", "exec", gr_simg, "gmx", "grompp",
            "-f", mdp, "-c", gro, "-p", top,
            "-o", tpr, "-maxwarn", str(maxwarn)
        ])
    else:
        print(f"✓ skipping grompp (tpr exists): {tpr}")

    # 2) g_mmpbsa
    if force or not Path(mm_xvg).exists():
        # Resolve g_mmpbsa path to absolute if relative
        if not os.path.isabs(g_mmpbsa_bin):
            g_mmpbsa_bin_abs = os.path.abspath(os.path.join(base_dir, g_mmpbsa_bin))
        else:
            g_mmpbsa_bin_abs = g_mmpbsa_bin
        
        if not Path(g_mmpbsa_bin_abs).exists():
            raise RuntimeError(f"g_mmpbsa binary not found: {g_mmpbsa_bin_abs}")
        
        # Interactive group selection via stdin
        cmd = [
            g_mmpbsa_bin_abs,
            "-f", xtc,
            "-s", tpr,
            "-mm", mm_xvg
        ]
        if ndx:
            cmd.extend(["-n", ndx])
        
        # g_mmpbsa prompts for groups on stdin
        input_data = f"{protein_group}\n{ligand_group}\n"
        print(f"  Groups: protein={protein_group}, ligand L01={ligand_group}")
        
        run_cmd(cmd, cwd=abs_folder, stdin_input=input_data)
    else:
        print(f"✓ skipping g_mmpbsa (mm_energy exists): {mm_xvg}")

    # 3) Parse .xvg
    metrics = parse_mm_xvg(mm_xvg)
    return {
        "combination": combo_id,
        "vs_folder": os.path.basename(folder),
        "type": combo_type,
        "mmpbsa_vdw": metrics["vdw"],
        "mmpbsa_elec": metrics["elec"],
        "mmpbsa_total": metrics["total"],
        "frames": metrics["frames"],
    }


def main():
    parser = argparse.ArgumentParser(
        description="Compute MM/PBSA metrics for VS_GR_* folders and write a CSV."
    )
    parser.add_argument("--case", help="Case (e.g. 1CJR_A) to filter folders", default=None)
    parser.add_argument("--type", choices=["GN", "LF", "all"], default="all", help="Filter by GN/LF")
    parser.add_argument("--base-pattern", default="VS_GR_*", help="Glob for VS_GR folders")
    parser.add_argument("--gr-simg", default="singularity/gr.simg", help="Singularity image with GROMACS")
    parser.add_argument("--g-mmpbsa-bin",
                        default="./GROMACS/analyze_results/Simulation_gromacs/analyze_trajectory/extra/g_mmpbsa",
                        help="Path to g_mmpbsa binary")
    parser.add_argument("--maxwarn", type=int, default=5, help="grompp -maxwarn value")
    parser.add_argument("--output-csv", default="mmpbsa_results.csv", help="Output CSV path")
    parser.add_argument("--force", action="store_true", help="Recompute even if outputs exist")
    args = parser.parse_args()

    # Backward compatibility: accept legacy image at repository root.
    if not os.path.exists(args.gr_simg):
        legacy_gr = "gr.simg"
        if os.path.exists(legacy_gr):
            print(f"ℹ Using legacy GROMACS image path: {legacy_gr}")
            args.gr_simg = legacy_gr

    # Find VS_GR_* folders
    pattern = args.base_pattern
    folders = [d for d in glob.glob(pattern) if os.path.isdir(d)]

    # Case filter: PDB base (e.g. 1AQD from 1AQD_C), same as calculate_md_rmsd
    case_base = None
    if args.case:
        case_upper = args.case.strip().upper()
        case_base = case_upper.split("_")[0] if "_" in case_upper else case_upper

    def _case_base_from_folder(path: str) -> Optional[str]:
        m = re.search(r"VS_GR_([a-zA-Z0-9]+)_", os.path.basename(path), re.IGNORECASE)
        return m.group(1).upper() if m else None

    type_filter = args.type.lower() if args.type in ("GN", "LF") else None

    def folder_ok(path: str) -> bool:
        name = os.path.basename(path)
        low = name.lower()
        if case_base:
            if _case_base_from_folder(path) != case_base:
                return False
        if type_filter and type_filter not in low:
            return False
        return True

    folders = [d for d in folders if folder_ok(d)]
    if args.case and case_base:
        print(f"Case filter: {args.case.strip()} (base {case_base}) → {len(folders)} folders")

    print("=" * 70)
    print(f"Computing MM/PBSA for {len(folders)} folders")
    print("=" * 70)
    if not folders:
        print("No folders matched the pattern.")
        sys.exit(0)

    results = []
    errors = 0
    skipped_incomplete = 0
    base_dir = os.getcwd()
    for folder in sorted(folders):
        abs_folder = os.path.abspath(folder)
        if not folder_has_mmpbsa_inputs(abs_folder):
            skipped_incomplete += 1
            print(f"⊘ Skipping incomplete MD folder: {folder}")
            continue
        try:
            res = process_folder(folder, args.gr_simg, args.g_mmpbsa_bin, args.maxwarn, args.force, base_dir)
            results.append(res)
            print(f"✓ {folder} -> combination {res['combination']}")
        except Exception as e:
            errors += 1
            print(f"✗ Error in {folder}: {e}", file=sys.stderr)

    # Write CSV
    if results:
        fieldnames = ["combination", "vs_folder", "type", "mmpbsa_vdw", "mmpbsa_elec", "mmpbsa_total", "frames"]
        with open(args.output_csv, "w", newline="", encoding="utf-8") as f:
            import csv
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in results:
                writer.writerow(row)
        print(f"\n✓ Wrote CSV: {args.output_csv} ({len(results)} rows)")
    else:
        print("\n⚠ No CSV written (no successful results).")
    if skipped_incomplete:
        print(f"ℹ Skipped incomplete folders: {skipped_incomplete}")

    if errors:
        print(f"\n⚠ Errors in {errors} folder(s). Check stderr above.")


if __name__ == "__main__":
    main()
