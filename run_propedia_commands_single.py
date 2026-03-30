#!/usr/bin/env python3
"""
Utility to execute commands one by one from
Propedia_pdbqt_5Frag/[PDB_ID]/Commands/[PDB_ID]_commands.txt.

This script executes commands sequentially, automatically waiting for each
SLURM job (and its dependencies) to complete before starting the next one.
It can filter by a specific complex (PDB_ID) and allows specifying the
Propedia folder name via command-line parameter.
"""

import argparse
import getpass
import math
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterator, Optional, Tuple


REPO_ROOT = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Execute Propedia commands one by one, waiting for each to complete. "
            "Can filter by specific complex (PDB_ID) and specify Propedia folder."
        )
    )
    parser.add_argument(
        "--complex",
        "--pdb-id",
        dest="complex_id",
        metavar="PDB_ID",
        help="Process only commands for this specific complex (PDB_ID).",
    )
    parser.add_argument(
        "--propedia-folder",
        default="Propedia_pdbqt_5Frag",
        help="Name of the Propedia folder (default: Propedia_pdbqt_5Frag).",
    )
    parser.add_argument(
        "--user",
        default=getpass.getuser(),
        help="SLURM user to monitor with squeue (default: current user).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only print the commands that would be executed.",
    )
    parser.add_argument(
        "--max-commands",
        type=int,
        help="Stop after executing this many commands (useful for smoke tests).",
    )
    parser.add_argument(
        "--no-wait",
        action="store_true",
        help="Don't wait for SLURM jobs to complete (default: always wait for sbatch jobs).",
    )
    parser.add_argument(
        "--check-interval",
        type=int,
        default=30,
        help="Seconds between checks when waiting for jobs (default: 30).",
    )
    parser.add_argument(
        "--skip-if-exists",
        action="store_true",
        help=(
            "Skip docking for a complex if results already exist and are valid: "
            "for each fragment, either a .tar.gz with the same name as the folder exists, "
            "or the folder exists and contains results (e.g. .pse/.pdbqt in results/best_scores, or files in molecules/)."
        ),
    )
    return parser.parse_args()


def get_propedia_root(propedia_folder: str) -> Path:
    """Get the Propedia root path based on folder name."""
    return REPO_ROOT / propedia_folder


def iter_command_files(propedia_root: Path, complex_id: Optional[str]) -> Iterator[Path]:
    """Iterate over command files, optionally filtered by complex_id."""
    if not propedia_root.exists():
        raise SystemExit(f"Missing directory: {propedia_root}")

    command_files = sorted(propedia_root.glob("*/Commands/*_commands.txt"))
    
    if complex_id:
        complex_id_lower = complex_id.lower()
        command_files = [
            path
            for path in command_files
            if path.parent.parent.name.lower() == complex_id_lower
        ]
        if not command_files:
            raise SystemExit(
                f"No command files found for complex '{complex_id}' in {propedia_root}"
            )
    
    for path in command_files:
        yield path


def iter_commands(path: Path) -> Iterator[Tuple[int, str]]:
    """Iterate over commands in a command file."""
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            command = raw_line.strip()
            if not command or command.startswith("#"):
                continue
            yield line_number, command


def extract_job_id_from_sbatch(output: str) -> Optional[int]:
    """Extract SLURM job ID from sbatch output."""
    # sbatch typically outputs: "Submitted batch job 12345"
    match = re.search(r"Submitted batch job (\d+)", output)
    if match:
        return int(match.group(1))
    return None


def get_job_dependencies(job_id: int) -> list[int]:
    """Get all job IDs that this job depends on."""
    result = subprocess.run(
        ["scontrol", "show", "job", str(job_id)],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return []
    
    # Look for Dependency= line, e.g., "Dependency=afterok:12345:12346"
    deps = []
    for line in result.stdout.splitlines():
        if line.startswith("Dependency="):
            dep_str = line.split("=", 1)[1].strip()
            if dep_str and dep_str != "(null)":
                # Parse dependencies like "afterok:12345:12346" or "afterok:12345"
                for part in dep_str.split(":"):
                    if part and part != "afterok" and part != "afterany" and part != "afternotok":
                        try:
                            deps.append(int(part))
                        except ValueError:
                            pass
            break
    return deps


def is_job_running(job_id: int, user: str) -> bool:
    """Check if a SLURM job (or any of its dependencies) is still running."""
    # Check the job itself
    result = subprocess.run(
        ["squeue", "-h", "-j", str(job_id), "-u", user],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode == 0 and result.stdout.strip():
        return True
    
    # Also check with scontrol in case job is pending due to dependencies
    scontrol_result = subprocess.run(
        ["scontrol", "show", "job", str(job_id)],
        check=False,
        capture_output=True,
        text=True,
    )
    if scontrol_result.returncode == 0:
        # If scontrol can show the job, it exists (might be pending or running)
        for line in scontrol_result.stdout.splitlines():
            if line.startswith("JobState="):
                state = line.split("=", 1)[1].strip().split()[0]  # Get first word (e.g., "PENDING", "RUNNING")
                if state in ("PENDING", "RUNNING", "CONFIGURING"):
                    return True
                break
    
    # Check dependencies
    deps = get_job_dependencies(job_id)
    for dep_id in deps:
        dep_result = subprocess.run(
            ["squeue", "-h", "-j", str(dep_id), "-u", user],
            check=False,
            capture_output=True,
            text=True,
        )
        if dep_result.returncode == 0 and dep_result.stdout.strip():
            return True
        
        # Also check with scontrol for dependencies
        dep_scontrol = subprocess.run(
            ["scontrol", "show", "job", str(dep_id)],
            check=False,
            capture_output=True,
            text=True,
        )
        if dep_scontrol.returncode == 0:
            for line in dep_scontrol.stdout.splitlines():
                if line.startswith("JobState="):
                    state = line.split("=", 1)[1].strip().split()[0]
                    if state in ("PENDING", "RUNNING", "CONFIGURING"):
                        return True
                    break
    
    return False


def wait_for_job(job_id: int, user: str, check_interval: int) -> None:
    """Wait for a SLURM job and all its dependencies to complete."""
    deps = get_job_dependencies(job_id)
    if deps:
        print(f"Job {job_id} has dependencies: {deps}", flush=True)
    
    print(f"Waiting for job {job_id} (and dependencies) to complete...", flush=True)
    check_count = 0
    while is_job_running(job_id, user):
        check_count += 1
        time.sleep(check_interval)
        if check_count % 4 == 0:  # Print status every 4 checks
            print(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
                f"Job {job_id} still running, checking again in {check_interval}s...",
                flush=True,
            )
    print(f"Job {job_id} and all dependencies completed.", flush=True)


def get_current_job_ids(user: str) -> set[int]:
    """Get set of currently running/pending job IDs for the user."""
    result = subprocess.run(
        ["squeue", "-h", "-u", user, "-o", "%i"],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return set()
    
    job_ids = set()
    for line in result.stdout.splitlines():
        line = line.strip()
        if line:
            try:
                job_ids.add(int(line))
            except ValueError:
                pass
    return job_ids


def find_new_jobs(before_jobs: set[int], after_jobs: set[int]) -> list[int]:
    """Find job IDs that are new (in after_jobs but not in before_jobs)."""
    return sorted(list(after_jobs - before_jobs))


def wait_for_all_jobs(job_ids: list[int], user: str, check_interval: int) -> None:
    """Wait for all jobs and their dependencies to complete."""
    if not job_ids:
        return
    
    print(f"Waiting for {len(job_ids)} job(s) to complete: {job_ids}", flush=True)
    
    # Get all dependencies for all jobs
    all_jobs_to_wait = set(job_ids)
    for job_id in job_ids:
        deps = get_job_dependencies(job_id)
        all_jobs_to_wait.update(deps)
    
    if len(all_jobs_to_wait) > len(job_ids):
        print(f"Including dependencies: {sorted(all_jobs_to_wait)}", flush=True)
    
    check_count = 0
    while True:
        # Check if any job is still running
        still_running = False
        for job_id in all_jobs_to_wait:
            if is_job_running(job_id, user):
                still_running = True
                break
        
        if not still_running:
            break
        
        check_count += 1
        time.sleep(check_interval)
        if check_count % 4 == 0:  # Print status every 4 checks
            print(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
                f"Jobs still running, checking again in {check_interval}s...",
                flush=True,
            )
    
    print(f"All jobs completed: {sorted(all_jobs_to_wait)}", flush=True)


def is_sbatch_command(command: str) -> bool:
    """Check if command is an sbatch command."""
    stripped = command.strip()
    return stripped.startswith("sbatch") or " sbatch " in stripped


def run_command(command: str, wait_for_jobs: bool, user: str, check_interval: int) -> None:
    """Execute a command and wait for SLURM jobs to complete (unless --no-wait is used)."""
    print(f"Executing: {command}", flush=True)
    
    # Get current jobs before execution
    before_jobs = set()
    if wait_for_jobs:
        before_jobs = get_current_job_ids(user)
        time.sleep(2)  # Small delay to ensure consistency
    
    # Execute the command
    if is_sbatch_command(command) and wait_for_jobs:
        # For sbatch commands, capture output to get job ID
        completed = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            raise subprocess.CalledProcessError(
                returncode=completed.returncode,
                cmd=command,
                stderr=completed.stderr,
            )
        
        # Extract job ID and wait for it
        job_id = extract_job_id_from_sbatch(completed.stdout)
        if job_id:
            wait_for_job(job_id, user, check_interval)
        else:
            # Fallback: detect new jobs
            time.sleep(3)  # Give time for job to appear in queue
            after_jobs = get_current_job_ids(user)
            new_jobs = find_new_jobs(before_jobs, after_jobs)
            if new_jobs:
                wait_for_all_jobs(new_jobs, user, check_interval)
            else:
                print(
                    f"Warning: Could not extract job ID from sbatch output: {completed.stdout}",
                    file=sys.stderr,
                    flush=True,
                )
    else:
        # For other commands (like ./ms.sh), run and detect new jobs
        completed = subprocess.run(command, shell=True)
        if completed.returncode != 0:
            raise subprocess.CalledProcessError(
                returncode=completed.returncode,
                cmd=command,
            )
        
        # Wait for any new SLURM jobs that were launched
        if wait_for_jobs:
            time.sleep(3)  # Give time for jobs to appear in queue
            after_jobs = get_current_job_ids(user)
            new_jobs = find_new_jobs(before_jobs, after_jobs)
            if new_jobs:
                wait_for_all_jobs(new_jobs, user, check_interval)
            else:
                # No new jobs detected, command likely completed synchronously
                print("No SLURM jobs detected, command completed.", flush=True)


def resolve_receptor_path(target: str) -> str:
    """Resolve receptor path, handling case-insensitive matching."""
    target_path = (REPO_ROOT / target).resolve()
    if target_path.exists():
        return target

    folder = target_path.parent
    if not folder.exists():
        print(
            f"Warning: Folder for receptor not found ({folder}). "
            "Keeping original -t argument.",
            file=sys.stderr,
            flush=True,
        )
        return target

    candidates = list(folder.glob("*.pdbqt"))
    if not candidates:
        print(
            f"Warning: No .pdbqt files found in {folder}. "
            "Keeping original -t argument.",
            file=sys.stderr,
            flush=True,
        )
        return target

    expected_name = target_path.name.lower()
    matching = next(
        (cand for cand in candidates if cand.name.lower() == expected_name),
        None,
    )
    if not matching:
        if len(candidates) == 1:
            matching = candidates[0]
        else:
            print(
                f"Warning: Multiple .pdbqt files in {folder} "
                "and none match expected name. Keeping original -t argument.",
                file=sys.stderr,
                flush=True,
            )
            return target

    try:
        rel_path = matching.resolve().relative_to(REPO_ROOT)
    except ValueError:
        rel_path = matching.resolve()
    return str(rel_path)


def adjust_target_tokens(tokens: list[str]) -> None:
    """Adjust -t tokens to resolve receptor paths."""
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token == "-t":
            if i + 1 >= len(tokens):
                raise ValueError("Found -t without a following argument.")
            original_target = tokens[i + 1]
            adjusted_target = resolve_receptor_path(original_target)
            if adjusted_target != original_target:
                print(
                    f"Adjusted receptor path: {original_target} -> {adjusted_target}",
                    flush=True,
                )
                tokens[i + 1] = adjusted_target
            i += 2
            continue
        i += 1


def round_grid_tokens(tokens: list[str]) -> None:
    """Round grid tokens (-gx, -gy, -gz) to integers."""
    grid_flags = {"-gx", "-gy", "-gz"}
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token in grid_flags:
            if i + 1 >= len(tokens):
                raise ValueError(f"Found {token} without a following argument.")
            try:
                value = float(tokens[i + 1])
            except ValueError as exc:
                raise ValueError(f"Non-numeric value for {token}: {tokens[i + 1]}") from exc
            rounded = str(int(math.ceil(value)))
            if rounded != tokens[i + 1]:
                print(
                    f"Rounded {token} value: {tokens[i + 1]} -> {rounded}",
                    flush=True,
                )
                tokens[i + 1] = rounded
            i += 2
            continue
        i += 1


def extract_output_dir(tokens: list[str]) -> Optional[str]:
    """Extract output directory from command tokens."""
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token == "-d":
            if i + 1 >= len(tokens):
                raise ValueError("Found -d without a following argument.")
            return tokens[i + 1]
        i += 1
    return None


def normalize_command(raw_command: str) -> Tuple[str, Optional[str]]:
    """Normalize a command by adjusting paths and rounding grid values."""
    tokens = shlex.split(raw_command)
    adjust_target_tokens(tokens)
    round_grid_tokens(tokens)
    output_dir = extract_output_dir(tokens)
    return shlex.join(tokens), output_dir


def _get_fragment_numbers_from_command_file(command_file: Path) -> list[int]:
    """Extract fragment numbers from -d output dirs in the command file (e.g. VS_GN_1BXP_B_Frag3 -> 3)."""
    frag_numbers = []
    frag_pattern = re.compile(r"_Frag(\d+)$", re.IGNORECASE)
    for _line_number, raw_command in iter_commands(command_file):
        try:
            tokens = shlex.split(raw_command)
            output_dir = extract_output_dir(tokens)
            if output_dir:
                match = frag_pattern.search(output_dir.strip("/").replace("\\", "/").split("/")[-1])
                if match:
                    frag_numbers.append(int(match.group(1)))
        except (ValueError, IndexError):
            continue
    return sorted(frag_numbers)


def docking_complete_for_complex(
    complex_name: str,
    command_file: Path,
    repo_root: Path,
) -> bool:
    """
    Return True if docking for this complex is already valid and can be skipped.

    We only look in the folder where results are moved ({case}_GN/), not in the
    repo root, so we avoid repeating checks and only consider "already organized"
    results as valid.

    For each fragment we require either:
    - A .tar.gz file (not .tar) with the same name as the folder (VS_GN_{case}_Frag{i}.tar.gz)
      in {case}_GN/, or
    - The folder exists in {case}_GN/ and contains results: at least one file in results/best_scores/
      (.pse, .pdbqt, .json, etc.) or at least one file in molecules/ (docking output).
    """
    case_upper = complex_name.upper()
    frag_numbers = _get_fragment_numbers_from_command_file(command_file)
    if not frag_numbers:
        return False

    # Only search under the results folder ({case}_GN/), not project root
    search_roots = [
        repo_root / f"{case_upper}_GN",
        repo_root / f"{complex_name.lower()}_GN",
    ]
    search_roots = list(dict.fromkeys(search_roots))  # sin duplicados

    for frag_i in frag_numbers:
        folder_name = f"VS_GN_{case_upper}_Frag{frag_i}"
        tar_name = f"{folder_name}.tar.gz"

        # Check 1: .tar.gz con el mismo nombre que la carpeta en {case}_GN/
        tar_found = any((root / tar_name).exists() for root in search_roots)
        if tar_found:
            continue

        # Check 2: carpeta en {case}_GN/ con resultados (best_scores con ficheros o molecules/ con ficheros)
        folder_has_results = False
        for root in search_roots:
            folder_path = root / folder_name
            if not folder_path.is_dir():
                continue
            best_scores = folder_path / "results" / "best_scores"
            if best_scores.is_dir():
                # Cualquier fichero en best_scores (no solo .pse)
                if any(f.is_file() for f in best_scores.iterdir()):
                    folder_has_results = True
                    break
            molecules = folder_path / "molecules"
            if molecules.is_dir() and any(f.is_file() for f in molecules.iterdir()):
                folder_has_results = True
                break
        if folder_has_results:
            continue

        return False
    return True


def organize_complex_results(complex_name: str) -> None:
    """Organize generated files and directories for a complex into a single directory."""
    target_dir_name = f"{complex_name}_GN"
    target_dir = REPO_ROOT / target_dir_name
    
    # Create target directory if it doesn't exist
    target_dir.mkdir(exist_ok=True)
    print(f"\nOrganizing results for complex {complex_name} into {target_dir_name}/", flush=True)
    
    moved_items = []
    
    # Find and move VS_GN_* directories related to this complex
    complex_pattern = complex_name.replace("_", "_").upper()  # Normalize case
    for item in REPO_ROOT.iterdir():
        if item.is_dir() and item.name.startswith("VS_GN_") and complex_pattern in item.name.upper():
            try:
                dest = target_dir / item.name
                if dest.exists():
                    print(f"  Removing existing {dest.name} in target directory...", flush=True)
                    shutil.rmtree(dest)
                shutil.move(str(item), str(dest))
                moved_items.append(item.name)
                print(f"  Moved directory: {item.name} -> {target_dir_name}/{item.name}", flush=True)
            except Exception as exc:
                print(
                    f"  Warning: Could not move {item.name}: {exc}",
                    file=sys.stderr,
                    flush=True,
                )
        
        # Also move tar.gz files related to this complex
        elif item.is_file() and item.name.endswith(".tar.gz"):
            if complex_pattern in item.name.upper():
                try:
                    dest = target_dir / item.name
                    if dest.exists():
                        print(f"  Removing existing {dest.name} in target directory...", flush=True)
                        dest.unlink()
                    shutil.move(str(item), str(dest))
                    moved_items.append(item.name)
                    print(f"  Moved file: {item.name} -> {target_dir_name}/{item.name}", flush=True)
                except Exception as exc:
                    print(
                        f"  Warning: Could not move {item.name}: {exc}",
                        file=sys.stderr,
                        flush=True,
                    )
    
    if moved_items:
        print(f"  Organized {len(moved_items)} item(s) into {target_dir_name}/", flush=True)
    else:
        print(f"  No items found to organize for complex {complex_name}", flush=True)


def main() -> int:
    args = parse_args()
    propedia_root = get_propedia_root(args.propedia_folder)
    
    if not propedia_root.exists():
        print(f"Directory not found: {propedia_root}", file=sys.stderr)
        return 1

    total_files = 0
    executed_commands = 0

    for command_file in iter_command_files(propedia_root, args.complex_id):
        total_files += 1
        complex_name = command_file.parent.parent.name
        print(f"\nProcessing complex: {complex_name}", flush=True)

        if args.skip_if_exists and docking_complete_for_complex(
            complex_name, command_file, REPO_ROOT
        ):
            print(
                f"Docking already valid for {complex_name} (tar or results in best_scores/molecules). Skipping.",
                flush=True,
            )
            continue

        for line_number, raw_command in iter_commands(command_file):
            if args.max_commands is not None and executed_commands >= args.max_commands:
                break
            
            try:
                command, output_dir = normalize_command(raw_command)
            except ValueError as exc:
                print(
                    f"Error processing command in {command_file}:{line_number}: {exc}",
                    file=sys.stderr,
                )
                return 1

            if output_dir:
                output_path = REPO_ROOT / output_dir
                if output_path.exists():
                    print(
                        f"Output directory `{output_dir}` already exists. Removing it...",
                        flush=True,
                    )
                    try:
                        if output_path.is_dir():
                            shutil.rmtree(output_path)
                        else:
                            output_path.unlink()
                        print(f"Removed existing output directory: {output_dir}", flush=True)
                    except Exception as exc:
                        print(
                            f"Warning: Could not remove {output_dir}: {exc}",
                            file=sys.stderr,
                            flush=True,
                        )
                        print("Continuing anyway...", flush=True)

            executed_commands += 1
            if args.dry_run:
                print(f"[DRY-RUN] {command_file}:{line_number} -> {command}")
                continue
            
            try:
                run_command(
                    command,
                    not args.no_wait,  # Wait by default, unless --no-wait is specified
                    args.user,
                    args.check_interval,
                )
            except subprocess.CalledProcessError as exc:
                print(
                    f"Command failed (file {command_file}, line {line_number}): {exc}",
                    file=sys.stderr,
                )
                return exc.returncode

        if args.max_commands is not None and executed_commands >= args.max_commands:
            break
        
        # Organize results for this complex after processing all commands
        if not args.dry_run:
            organize_complex_results(complex_name)

    print(
        f"\nFinished processing {executed_commands} commands across {total_files} files.",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

