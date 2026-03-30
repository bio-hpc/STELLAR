#!/usr/bin/env python3
"""
Batch-launch molecular dynamics (GROMACS) simulations for valid combinations.

Runs:
./sm.sh -t {combination_dir}/1a1m_A.pdb -q {combination_dir}/query_combination_{num}/ -o VS -s GR ...
"""

import os
import glob
import subprocess
import sys
import argparse
import shutil
import re
import time
import getpass
import tempfile
from pathlib import Path


def find_vs_folders(combo_dir, prefix_type, debug=False):
    """
    Find VS_GR_* folders for a given combination.
    sm.sh creates these at the project root with pattern:
    VS_GR_{case}_{type}_query_combination_{num}_*

    Args:
        combo_dir: Combination directory (e.g. valid_GN_final/combination_102)
        prefix_type: 'GN' or 'LF'
        debug: Verbose folder-matching logs

    Returns:
        list of VS_GR_* directory paths
    """
    vs_folders = []
    combo_number = os.path.basename(combo_dir).replace("combination_", "")
    
    # VS_GR_* at project root; match combo id exactly (avoid 1 matching 10)
    project_root = os.getcwd()
    
    # Infer case id from parent dir: valid_{TYPE}_final_{CASE} or valid_{TYPE}_final
    base_dir = os.path.dirname(combo_dir)
    case_pattern = None
    if base_dir and '_' in base_dir:
        parts = base_dir.split('_')
        # Pattern valid_{TYPE}_final_{CASE}
        if len(parts) >= 4 and parts[0] == 'valid' and parts[2] == 'final':
            case_pattern = '_'.join(parts[3:]) if len(parts) > 3 else None
    
    # Else infer case from PDB name if present
    if not case_pattern:
        pdb_files = glob.glob(os.path.join(combo_dir, f"*_{prefix_type}.pdb"))
        if pdb_files:
            # First PDB basename → case tokens
            pdb_name = os.path.basename(pdb_files[0])
            parts = pdb_name.replace(f"_{prefix_type}.pdb", "").split('_')
            if len(parts) >= 2:
                potential_case = '_'.join(parts[:2])
                # Verificar si hay una carpeta VS_GR con este caso
                test_pattern = os.path.join(project_root, f"VS_GR_{potential_case.lower()}*")
                if glob.glob(test_pattern):
                    case_pattern = potential_case
                elif len(parts) >= 3:
                    potential_case = '_'.join(parts[:3])
                    test_pattern = os.path.join(project_root, f"VS_GR_{potential_case.lower()}*")
                    if glob.glob(test_pattern):
                        case_pattern = potential_case
    
    patterns = [
        os.path.join(project_root, f"VS_GR_*query_combination_{combo_number}_*"),
        os.path.join(project_root, f"VS_GR_*query_combination_{combo_number}/"),
    ]
    matches = []
    for pattern in patterns:
        matches.extend(glob.glob(pattern))
    
    # Manual check: glob cannot do lookahead for exact combo id
    all_vs_gr = glob.glob(os.path.join(project_root, "VS_GR_*"))
    for folder in all_vs_gr:
        folder_name = os.path.basename(folder)
        exact_pattern = rf"query_combination_{re.escape(combo_number)}(_|/|$)"
        if re.search(exact_pattern, folder_name) and folder not in matches:
            matches.append(folder)
    
    # Filter by GN/LF and case id
    for match in matches:
        if os.path.isdir(match):
            folder_name = os.path.basename(match)
            
            # Type must appear in folder name
            type_pattern = rf"_{prefix_type}_"
            if not re.search(type_pattern, folder_name):
                if debug:
                    print(f"    [DEBUG] Skipped (type mismatch): {folder_name}")
                continue
            
            if case_pattern:
                # El caso puede estar en diferentes formatos (4WVI_D, 4wvi_D, etc.)
                case_lower = case_pattern.lower()
                case_upper = case_pattern.upper()
                case_mixed = case_pattern
                if not (case_lower in folder_name.lower() or case_upper in folder_name or case_mixed in folder_name):
                    if debug:
                        print(f"    [DEBUG] Skipped (case mismatch): {folder_name}")
                    continue
            
            molecules_dir = os.path.join(match, "molecules")
            if os.path.exists(molecules_dir):
                vs_folders.append(match)
    
    if debug:
        print(f"    [DEBUG] Search root: {project_root}")
        print(f"    [DEBUG] Pattern: VS_GR_*query_combination_{combo_number}*")
        print(f"    [DEBUG] VS_GR_* folders matched: {len(vs_folders)}")
        for folder in vs_folders:
            print(f"    [DEBUG]   - {os.path.basename(folder)}")
    
    # Legacy: VS_GR under query_combination_X/
    query_dir = os.path.join(combo_dir, f"query_combination_{combo_number}")
    if os.path.exists(query_dir):
        query_patterns = [
            os.path.join(query_dir, "VS_GR_*"),
            os.path.join(query_dir, "VS_*"),
        ]
        for pattern in query_patterns:
            matches = glob.glob(pattern)
            for match in matches:
                if os.path.isdir(match):
                    molecules_dir = os.path.join(match, "molecules")
                    if os.path.exists(molecules_dir) and match not in vs_folders:
                        vs_folders.append(match)
    
    return vs_folders


def has_pdb_in_molecules(vs_folder):
    """
    Verifica si existe al menos un archivo PDB en molecules/ dentro de una carpeta VS_*.
    
    Args:
        vs_folder: Ruta a la carpeta VS_*
    
    Returns:
        bool: True si existe al menos un PDB, False en caso contrario
    """
    molecules_dir = os.path.join(vs_folder, "molecules")
    if not os.path.exists(molecules_dir):
        return False
    
    pdb_files = glob.glob(os.path.join(molecules_dir, "*.pdb"))
    return len(pdb_files) > 0


def check_existing_pdb(combo_dir, prefix_type, debug=False):
    """
    Whether this combination already has output PDBs under VS_*.

    Args:
        combo_dir: Combination directory
        prefix_type: 'GN' or 'LF'
        debug: Verbose logs

    Returns:
        tuple: (has_pdb, vs_folders_with_pdb, vs_folders_without_pdb, pdb_count)
    """
    vs_folders = find_vs_folders(combo_dir, prefix_type, debug)
    
    if not vs_folders:
        return False, [], [], 0
    
    total_pdbs = 0
    folders_with_pdb = []
    folders_without_pdb = []
    
    for vs_folder in vs_folders:
        if has_pdb_in_molecules(vs_folder):
            pdb_count = len(glob.glob(os.path.join(vs_folder, "molecules", "*.pdb")))
            total_pdbs += pdb_count
            folders_with_pdb.append(vs_folder)
        else:
            folders_without_pdb.append(vs_folder)
    
    return len(folders_with_pdb) > 0, folders_with_pdb, folders_without_pdb, total_pdbs


def check_topology_error_in_first_job(vs_folder):
    """
    Scan first MD step logs (jobs_out/0-1.err or 0-1.out) for topology/GROMACS failures.

    Args:
        vs_folder: VS_GR_* path

    Returns:
        tuple: (has_error, message)
    """
    jobs_out = os.path.join(vs_folder, "jobs_out")
    if not os.path.isdir(jobs_out):
        return False, None
    # GROMACS launcher: 0-1.err and 0-1.out
    err_file = os.path.join(jobs_out, "0-1.err")
    out_file = os.path.join(jobs_out, "0-1.out")
    text = ""
    for path in (err_file, out_file):
        if os.path.isfile(path):
            try:
                with open(path, "r", encoding="utf-8", errors="replace") as f:
                    text += f.read() + "\n"
            except OSError:
                pass
    if not text:
        return False, None
    text_lower = text.lower()
    err_ref = os.path.join(os.path.basename(vs_folder), "jobs_out", "0-1.err")
    msg = f"Topology/GROMACS failure (see {err_ref})"
    if any(x in text_lower for x in (
        "fatal error",
        "error in topology",
        "cannot resolve",
        "invalid order",
        "no such file",
        "command not found",
        "failed to",
        "back off",
        "exiting",
    )):
        return True, msg
    if "topology" in text_lower and ("error" in text_lower or "fatal" in text_lower or "cannot" in text_lower):
        return True, msg
    if ("gmx" in text_lower or "gromacs" in text_lower) and ("error" in text_lower or "fatal" in text_lower):
        return True, msg
    return False, None


def remove_vs_folders(vs_folders, dry_run=False):
    """
    Delete the given VS_* directories.

    Args:
        vs_folders: Paths to remove
        dry_run: List only; do not delete

    Returns:
        tuple: (success, deleted_count, error_message)
    """
    if not vs_folders:
        return True, 0, None
    
    deleted_count = 0
    errors = []
    
    for vs_folder in vs_folders:
        if dry_run:
            print(f"    [DRY RUN] Would remove: {vs_folder}")
            deleted_count += 1
        else:
            try:
                if os.path.exists(vs_folder):
                    shutil.rmtree(vs_folder)
                    deleted_count += 1
                    print(f"    🗑️  Removed: {os.path.basename(vs_folder)}")
            except Exception as e:
                errors.append(f"Error removing {vs_folder}: {e}")
    
    if errors:
        return False, deleted_count, "; ".join(errors)
    
    return True, deleted_count, None


def extract_job_id_from_sbatch(output: str):
    """Extract first SLURM job ID from sbatch output (legacy single-ID)."""
    ids = extract_job_ids_from_sbatch_output(output)
    return ids[0] if ids else None


def extract_job_ids_from_sbatch_output(output: str):
    """Extract all SLURM job IDs from sm.sh/sbatch output.
    Matches: "Submitted batch job 12345" (sbatch) and "+ JOB: 12345" (GROMACS execute_jobs.sh).
    Returns list of unique job IDs in order of appearance, so we only wait for jobs from this invocation.
    """
    ids = []
    seen = set()
    # Both patterns: sbatch direct output and GROMACS launcher's echo "+ JOB: <id>"
    for m in re.finditer(r"Submitted batch job (\d+)|(?:\+ JOB:\s*)(\d+)", output):
        jid = int(m.group(1) or m.group(2))
        if jid not in seen:
            seen.add(jid)
            ids.append(jid)
    return ids


def get_job_dependencies(job_id: int):
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


def cancel_job(job_id: int, cancel_dependencies: bool = False) -> bool:
    """Cancel a SLURM job (and optionally jobs that depend on it). Returns True if scancel succeeded."""
    to_cancel = [job_id]
    if cancel_dependencies:
        to_cancel.extend(get_job_dependencies(job_id))
    for jid in to_cancel:
        result = subprocess.run(
            ["scancel", str(jid)],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            print(f"  Job {jid} cancelled (scancel).", flush=True)
        else:
            # Job may already be finished
            if "Invalid job" not in (result.stderr or ""):
                print(f"  Could not cancel job {jid}: {result.stderr or result.stdout}", flush=True)
    return True


def wait_for_job(job_id: int, user: str, check_interval: int, timeout_seconds: int = None) -> bool:
    """
    Wait for a SLURM job and all its dependencies to complete.
    If timeout_seconds is set and the job runs longer, cancel it and return False.
    Returns True if the job completed normally, False if it was cancelled by timeout.
    """
    deps = get_job_dependencies(job_id)
    if deps:
        print(f"  Job {job_id} has dependencies: {deps}", flush=True)
    
    timeout_msg = f" (timeout {timeout_seconds}s)" if timeout_seconds else ""
    print(f"  Waiting for job {job_id} (and dependencies) to complete{timeout_msg}...", flush=True)
    check_count = 0
    start_time = time.monotonic()
    while is_job_running(job_id, user):
        elapsed = time.monotonic() - start_time
        if timeout_seconds is not None and elapsed >= timeout_seconds:
            print(
                f"  [{time.strftime('%Y-%m-%d %H:%M:%S')}] Job {job_id} exceeded time limit ({timeout_seconds}s). Cancelling; continuing.",
                flush=True,
            )
            cancel_job(job_id, cancel_dependencies=False)
            return False
        check_count += 1
        time.sleep(check_interval)
        if check_count % 4 == 0:  # Print status every 4 checks
            elapsed_str = f"{int(elapsed)}s" if timeout_seconds is None else f"{int(elapsed)}/{timeout_seconds}s"
            print(
                f"  [{time.strftime('%Y-%m-%d %H:%M:%S')}] Job {job_id} still running ({elapsed_str}), checking again in {check_interval}s...",
                flush=True,
            )
    print(f"  Job {job_id} and all dependencies completed.", flush=True)
    return True


def get_current_job_ids(user: str):
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


def find_new_jobs(before_jobs: set, after_jobs: set):
    """Find job IDs that are new (in after_jobs but not in before_jobs)."""
    return sorted(list(after_jobs - before_jobs))


def wait_for_all_jobs(job_ids: list, user: str, check_interval: int, timeout_seconds: int = None) -> bool:
    """
    Wait for all jobs and their dependencies to complete.
    If timeout_seconds is set and exceeded, cancel all job_ids and return False.
    Returns True if all completed normally, False if cancelled by timeout.
    """
    if not job_ids:
        return True
    
    timeout_msg = f" (timeout {timeout_seconds}s)" if timeout_seconds else ""
    print(f"  Waiting for {len(job_ids)} job(s) to complete: {job_ids}{timeout_msg}", flush=True)
    
    # Get all dependencies for all jobs
    all_jobs_to_wait = set(job_ids)
    for job_id in job_ids:
        deps = get_job_dependencies(job_id)
        all_jobs_to_wait.update(deps)
    
    if len(all_jobs_to_wait) > len(job_ids):
        print(f"  Including dependencies: {sorted(all_jobs_to_wait)}", flush=True)
    
    check_count = 0
    start_time = time.monotonic()
    while True:
        elapsed = time.monotonic() - start_time
        if timeout_seconds is not None and elapsed >= timeout_seconds:
            print(
                f"  [{time.strftime('%Y-%m-%d %H:%M:%S')}] Jobs exceeded time limit ({timeout_seconds}s). Cancelling; continuing.",
                flush=True,
            )
            for jid in job_ids:
                cancel_job(jid, cancel_dependencies=False)
            return False
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
        if check_count % 4 == 0:
            elapsed_str = f"{int(elapsed)}s" if timeout_seconds is None else f"{int(elapsed)}/{timeout_seconds}s"
            print(
                f"  [{time.strftime('%Y-%m-%d %H:%M:%S')}] Jobs still running ({elapsed_str}), checking again in {check_interval}s...",
                flush=True,
            )
    print(f"  All jobs completed: {sorted(all_jobs_to_wait)}", flush=True)
    return True


def run_simulation(combo_dir, prefix_type, sm_script="./sm.sh", dry_run=False, check_existing=False, debug=False, wait_for_completion=False, user=None, check_interval=30, job_timeout_seconds=None, **kwargs):
    """
    Launch MD for one combination via sm.sh.

    Args:
        combo_dir: e.g. valid_LF_final/combination_763
        prefix_type: 'GN' or 'LF'
        sm_script: Path to sm.sh (default ./sm.sh)
        dry_run: Print command only
        check_existing: Skip if PDB already under VS_*
        wait_for_completion: Wait for SLURM job(s) before returning
        user: SLURM user for squeue (default: current user)
        check_interval: Seconds between queue polls
        job_timeout_seconds: Cancel job if exceeded (None = no limit)
        **kwargs: Passed through to sm.sh

    Returns:
        tuple: (success, error_message, skipped)
    """
    combo_name = os.path.basename(combo_dir)
    
    # combination_763 -> 763
    combo_number = combo_name.replace("combination_", "")
    
    # Prefer *_{GN|LF}.pdb
    pdb_files = glob.glob(os.path.join(combo_dir, "*.pdb"))
    
    pdb_file = None
    for pdb_path in pdb_files:
        if pdb_path.endswith(f"_{prefix_type}.pdb"):
            pdb_file = pdb_path
            break
    
    # Fallback: any .pdb in folder
    if pdb_file is None:
        if pdb_files:
            pdb_file = pdb_files[0]
        else:
            pdb_file = os.path.join(combo_dir, f"*_{prefix_type}.pdb")
    
    query_dir = os.path.join(combo_dir, f"query_combination_{combo_number}")
    
    if check_existing:
        has_pdb, vs_folders_with_pdb, vs_folders_without_pdb, pdb_count = check_existing_pdb(combo_dir, prefix_type, debug)

        if has_pdb:
            return True, f"PDB already present ({pdb_count} file(s) in {len(vs_folders_with_pdb)} folder(s))", True

        if vs_folders_without_pdb:
            print(f"  🗑️  Removing {len(vs_folders_without_pdb)} empty VS_* folder(s) before launch...")
            success, deleted_count, error = remove_vs_folders(vs_folders_without_pdb, dry_run)
            if not success:
                return False, f"Error removing VS_* folders: {error}", False
            if not dry_run:
                print(f"    ✓ Removed {deleted_count} folder(s)")

    if not os.path.exists(pdb_file):
        return False, f"PDB not found: {pdb_file}", False

    if not os.path.exists(query_dir):
        return False, f"Query directory not found: {query_dir}", False

    if not os.path.exists(sm_script) and not dry_run:
        return False, f"Script not found: {sm_script}"

    current_dir = os.getcwd()
    
    pdb_rel = os.path.relpath(pdb_file, current_dir)
    query_rel = os.path.relpath(query_dir, current_dir)
    if not query_rel.endswith('/'):
        query_rel = query_rel + '/'
    
    cmd = [
        sm_script,
        "-t", pdb_rel,
        "-q", query_rel,
        "-o", kwargs.get("output", "VS"),
        "-s", kwargs.get("software", "GR"),
        "-j", str(kwargs.get("jobs", 1)),
        "-prp", kwargs.get("prp", "na"),
        "-prl", kwargs.get("prl", "na"),
        "-gp", str(kwargs.get("gp", 1)),
        "-co", str(kwargs.get("co", 24)),
        "-td", str(kwargs.get("td", 259200)),
        "-el", kwargs.get("el", "n"),
        "-step_min", str(kwargs.get("step_min", 1000)),
        "-write_data", str(kwargs.get("write_data", 3000)),
        "-step_npt", str(kwargs.get("step_npt", 2000)),
        "-step_nvt", str(kwargs.get("step_nvt", 2000)),
        "-step_md", str(kwargs.get("step_md", 2000)),
        "-solvent", kwargs.get("solvent", "tip3p"),
        "-force_field", kwargs.get("force_field", "amber99sb"),
        "-solvatation", kwargs.get("solvatation", "SOL"),
        "-temp", str(kwargs.get("temp", 300)),
        "-bt", kwargs.get("bt", "dodecahedron"),
        "-padding_grid", str(kwargs.get("padding_grid", 1.5)),
        "-seedg", str(kwargs.get("seedg", 2077)),
        "-prefix_gromacs", kwargs.get("prefix_gromacs", "/usr/local/gromacs/avx2_256/bin/gmx"),
    ]
    
    # Optional NPT pressure
    if "pressure_npt" in kwargs:
        cmd.extend(["-pressure_npt", str(kwargs["pressure_npt"])])
    
    # Per-dynamics walltime for sm.sh (-tj HH:MM:SS)
    cmd.extend(["-tj", kwargs.get("tj", "00:30:00")])
    
    if dry_run:
        print(f"  [DRY RUN] {' '.join(cmd)}")
        return True, None, False
    
    try:
        # Get current jobs before execution if we need to wait
        before_jobs = set()
        if wait_for_completion:
            if user is None:
                user = getpass.getuser()
            before_jobs = get_current_job_ids(user)
            time.sleep(2)  # Small delay to ensure consistency
        
        # Stream stdout to file to limit RAM (OOM on small nodes)
        stdout_path = None
        if wait_for_completion:
            fd, stdout_path = tempfile.mkstemp(suffix='.out', text=True)
            os.close(fd)
        fout = open(stdout_path, 'w') if stdout_path else None
        try:
            result = subprocess.run(
                cmd,
                stdout=fout if fout else subprocess.DEVNULL,
                stderr=subprocess.STDOUT if fout else subprocess.DEVNULL,
                text=True,
                check=False,
                cwd=current_dir
            )
        finally:
            if fout:
                fout.close()
        # Parse job IDs from this sm.sh run only
        job_ids_from_output = []
        err_lines = []
        if stdout_path and os.path.exists(stdout_path):
            try:
                with open(stdout_path, 'r') as f:
                    stdout_content = f.read()
                job_ids_from_output = extract_job_ids_from_sbatch_output(stdout_content)
                err_lines = stdout_content.splitlines()[-5:] if stdout_content else []
            finally:
                try:
                    os.unlink(stdout_path)
                except OSError:
                    pass
        
        if result.returncode == 0:
            # If we need to wait for completion, wait only for jobs from this sm.sh invocation
            if wait_for_completion:
                if user is None:
                    user = getpass.getuser()
                if job_ids_from_output:
                    completed_ok = wait_for_all_jobs(job_ids_from_output, user, check_interval, timeout_seconds=job_timeout_seconds)
                    if not completed_ok:
                        print(f"  ⏱️  Run cancelled by time limit; continuing.", flush=True)
                else:
                    time.sleep(3)  # Give time for job to appear in queue
                    after_jobs = get_current_job_ids(user)
                    new_jobs = find_new_jobs(before_jobs, after_jobs)
                    if new_jobs:
                        completed_ok = wait_for_all_jobs(new_jobs, user, check_interval, timeout_seconds=job_timeout_seconds)
                        if not completed_ok:
                            print(f"  ⏱️  Run cancelled by time limit; continuing.", flush=True)
                    else:
                        print("  No SLURM jobs detected, command completed.", flush=True)
                # After first MD step: topology error in jobs_out/0-1.err?
                vs_folders = find_vs_folders(combo_dir, prefix_type, debug=False)
                for vf in vs_folders:
                    has_err, err_msg = check_topology_error_in_first_job(vf)
                    if has_err:
                        return False, err_msg, False
            return True, None, False
        else:
            error_preview = '\n'.join(err_lines) if err_lines else f"Exit code {result.returncode}"
            return False, f"Execution error:\n{error_preview}", False

    except Exception as e:
        return False, f"Exception: {e}", False


def process_all_combinations(base_dir, prefix_type, sm_script="./sm.sh", dry_run=False, check_existing=False, debug=False, max_jobs=None, wait_for_completion=False, user=None, check_interval=30, job_timeout_seconds=None, **kwargs):
    """
    Run MD for every combination_* under base_dir.

    Args:
        base_dir: Root containing combination_* folders
        prefix_type: 'GN' or 'LF'
        sm_script: Path to sm.sh
        dry_run: Print only
        check_existing: Skip combinations that already have PDB output
        max_jobs: Cap on newly launched runs (None = no cap)
        wait_for_completion: Wait for each SLURM job before next launch
        user: SLURM user for monitoring
        check_interval: squeue poll interval
        **kwargs: sm.sh parameters

    Returns:
        tuple: (total_processed, total_skipped, total_errors, errors_list)
    """
    combo_dirs = sorted(glob.glob(os.path.join(base_dir, "combination_*")))
    
    if not combo_dirs:
        print(f"⚠ No combination_* folders in {base_dir}")
        return 0, 0, 0, []
    
    total_combinations = len(combo_dirs)
    if check_existing:
        print(f"Found {total_combinations} combinations")
        print("Checking which already have PDB output...")
        has_pdb_count = 0
        for combo_dir in combo_dirs:
            has_pdb, _, _, _ = check_existing_pdb(combo_dir, prefix_type)
            if has_pdb:
                has_pdb_count += 1
        to_process = total_combinations - has_pdb_count
        print(f"  ✓ {has_pdb_count} combinations already have PDB (will skip)")
        print(f"  → {to_process} combinations will run")
    else:
        print(f"Found {total_combinations} combinations")

    if max_jobs is not None:
        print(f"⚠ Job cap: {max_jobs} (at most {max_jobs} new launches)")

    if dry_run:
        print("⚠ DRY RUN — no commands will be executed")
    if check_existing:
        print("⚠ CHECK EXISTING — skip combinations with PDB already")
    print()
    
    total_processed = 0
    total_skipped = 0
    total_errors = 0
    errors_list = []
    
    for i, combo_dir in enumerate(combo_dirs, 1):
        combo_name = os.path.basename(combo_dir)
        print(f"[{i}/{len(combo_dirs)}] {combo_name}")
        
        success, error, skipped = run_simulation(combo_dir, prefix_type, sm_script, dry_run, check_existing, debug, wait_for_completion, user, check_interval, job_timeout_seconds, **kwargs)
        
        if skipped:
            total_skipped += 1
            print(f"  ⊘ Skipped: {error}")
        elif success:
            total_processed += 1
            if not dry_run:
                print(f"  ✓ Simulation launched")
        else:
            total_errors += 1
            errors_list.append((combo_name, error))
            print(f"  ✗ {error}")
        
        if max_jobs is not None and total_processed >= max_jobs:
            print(f"\n⚠ Reached job cap ({max_jobs}). Stopping.")
            break
        
        if (i % 10 == 0) and not dry_run:
            print(f"  Progress: {total_processed} run, {total_skipped} skipped, {total_errors} errors...")
    
    return total_processed, total_skipped, total_errors, errors_list


def main():
    parser = argparse.ArgumentParser(
        description="Launch MD simulations for all valid combination_* folders"
    )
    parser.add_argument(
        'type',
        choices=['GN', 'LF', 'all'],
        help='Combination type: GN, LF, or all (both)'
    )
    parser.add_argument(
        '--base-dir-gn',
        help='Base directory for GN (default: valid_GN_final)'
    )
    parser.add_argument(
        '--base-dir-lf',
        help='Base directory for LF (default: valid_LF_final)'
    )
    parser.add_argument(
        '--sm-script',
        default='./sm.sh',
        help='Path to sm.sh (default: ./sm.sh)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print commands only; do not run'
    )
    parser.add_argument(
        '--check-existing',
        action='store_true',
        help='Skip combinations that already have PDB output under VS_*'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Verbose VS_* folder discovery'
    )
    parser.add_argument(
        '--max-jobs',
        type=int,
        default=None,
        help='Max new launches (e.g. 10 for a controlled batch)'
    )
    parser.add_argument(
        '--sequential',
        '--wait',
        action='store_true',
        dest='wait_for_completion',
        help='Wait for each SLURM job to finish before launching the next'
    )
    parser.add_argument(
        '--user',
        default=getpass.getuser(),
        help='SLURM user for squeue (default: current user)'
    )
    parser.add_argument(
        '--check-interval',
        type=int,
        default=30,
        help='Seconds between queue polls when waiting (default: 30)'
    )
    parser.add_argument(
        '--job-timeout',
        type=int,
        default=0,
        metavar='MINUTES',
        help='When using --sequential: cancel a dynamics if walltime exceeds this many minutes; 0 = no cancel by this script (default: 0). Per-step GROMACS walltime is still set by -tj / sm.sh.'
    )

    # sm.sh passthrough defaults
    parser.add_argument('--output', '-o', default='VS', help='Output (default: VS)')
    parser.add_argument('--software', '-s', default='GR', help='Software (default: GR)')
    parser.add_argument('--jobs', '-j', type=int, default=1, help='Jobs (default: 1)')
    parser.add_argument('--prp', default='na', help='PRP (default: na)')
    parser.add_argument('--prl', default='na', help='PRL (default: na)')
    parser.add_argument('--gp', type=int, default=1, help='GP (default: 1)')
    parser.add_argument('--co', type=int, default=24, help='CO (default: 24)')
    parser.add_argument('--td', type=int, default=259200, help='TD (default: 259200)')
    parser.add_argument('--el', default='n', help='EL (default: n)')
    parser.add_argument('--step-min', type=int, default=1000, help='Step MIN (default: 1000)')
    parser.add_argument('--write-data', type=int, default=3000, help='Write data (default: 3000)')
    parser.add_argument('--step-npt', type=int, default=2000, help='Step NPT (default: 2000)')
    parser.add_argument('--step-nvt', type=int, default=2000, help='Step NVT (default: 2000)')
    parser.add_argument('--step-md', type=int, default=2000, help='Step MD (default: 2000)')
    parser.add_argument('--solvent', default='tip3p', help='Solvent (default: tip3p)')
    parser.add_argument('--force-field', default='amber99sb', help='Force field (default: amber99sb)')
    parser.add_argument('--solvatation', default='SOL', help='Solvatation (default: SOL)')
    parser.add_argument('--temp', type=int, default=300, help='Temperature (default: 300)')
    parser.add_argument('--bt', default='dodecahedron', help='Box type (default: dodecahedron)')
    parser.add_argument('--padding-grid', type=float, default=1.5, help='Padding grid (default: 1.5)')
    parser.add_argument('--seedg', type=int, default=2077, help='Seed G (default: 2077)')
    parser.add_argument('--prefix-gromacs', default='/usr/local/gromacs/avx2_256/bin/gmx',
                       help='Prefix GROMACS (default: /usr/local/gromacs/avx2_256/bin/gmx)')
    parser.add_argument('--pressure-npt', type=float, default=1.0, help='Pressure NPT (default: 1.0)')
    parser.add_argument('--tj', default='00:30:00', metavar='HH:MM:SS',
                        help='Per-dynamics job walltime passed to sm.sh -tj (default: 00:30:00)')

    args = parser.parse_args()

    if args.base_dir_gn is None:
        args.base_dir_gn = "valid_GN_final"
    if args.base_dir_lf is None:
        args.base_dir_lf = "valid_LF_final"
    
    # sm.sh keyword args
    kwargs = {
        "output": args.output,
        "software": args.software,
        "jobs": args.jobs,
        "prp": args.prp,
        "prl": args.prl,
        "gp": args.gp,
        "co": args.co,
        "td": args.td,
        "el": args.el,
        "step_min": args.step_min,
        "write_data": args.write_data,
        "step_npt": args.step_npt,
        "step_nvt": args.step_nvt,
        "step_md": args.step_md,
        "solvent": args.solvent,
        "force_field": args.force_field,
        "solvatation": args.solvatation,
        "temp": args.temp,
        "bt": args.bt,
        "padding_grid": args.padding_grid,
        "seedg": args.seedg,
        "prefix_gromacs": args.prefix_gromacs,
        "pressure_npt": args.pressure_npt,
        "tj": args.tj,
    }
    
    print("=" * 70)
    print("Molecular dynamics batch launcher")
    print("=" * 70)
    print(f"Script: {args.sm_script}")
    print(f"Args: -o {kwargs['output']} -s {kwargs['software']} -j {kwargs['jobs']} "
          f"-td {kwargs['td']} -temp {kwargs['temp']} -pressure_npt {kwargs['pressure_npt']} -tj {kwargs['tj']}")
    if args.wait_for_completion:
        print(f"⚠ SEQUENTIAL: wait for each job before the next")
        print(f"  User: {args.user}, poll interval: {args.check_interval}s")
        job_timeout_seconds = (args.job_timeout * 60) if args.job_timeout else None
        if job_timeout_seconds:
            print(f"  Per-dynamics cancel timeout: {args.job_timeout} min (then cancel and continue)")
        else:
            print(f"  Per-dynamics cancel timeout: none")
    else:
        job_timeout_seconds = None
    print()
    
    total_processed = 0
    total_skipped = 0
    total_errors = 0
    all_errors = []
    
    if args.type in ['GN', 'all']:
        if os.path.exists(args.base_dir_gn):
            print(f"\n📁 GN: {args.base_dir_gn}")
            _timeout = (args.job_timeout * 60) if args.job_timeout and args.wait_for_completion else None
            processed, skipped, errors, errors_list = process_all_combinations(
                args.base_dir_gn, 'GN', args.sm_script, args.dry_run, args.check_existing, args.debug, args.max_jobs, args.wait_for_completion, args.user, args.check_interval, _timeout, **kwargs
            )
            total_processed += processed
            total_skipped += skipped
            total_errors += errors
            all_errors.extend(errors_list)
        else:
            print(f"\n⚠ GN directory missing: {args.base_dir_gn}")

    if args.type in ['LF', 'all']:
        if os.path.exists(args.base_dir_lf):
            print(f"\n📁 LF: {args.base_dir_lf}")
            _timeout = (args.job_timeout * 60) if args.job_timeout and args.wait_for_completion else None
            processed, skipped, errors, errors_list = process_all_combinations(
                args.base_dir_lf, 'LF', args.sm_script, args.dry_run, args.check_existing, args.debug, args.max_jobs, args.wait_for_completion, args.user, args.check_interval, _timeout, **kwargs
            )
            total_processed += processed
            total_skipped += skipped
            total_errors += errors
            all_errors.extend(errors_list)
        else:
            print(f"\n⚠ LF directory missing: {args.base_dir_lf}")

    print("\n" + "=" * 70)
    if args.dry_run:
        print(f"Summary (DRY RUN): {total_processed} combinations would run")
    else:
        summary_parts = [f"{total_processed} simulations launched"]
        if total_skipped > 0:
            summary_parts.append(f"{total_skipped} skipped (PDB already present)")
        if total_errors > 0:
            summary_parts.append(f"{total_errors} errors")
        print(f"Summary: {', '.join(summary_parts)}")
    print("=" * 70)

    if all_errors:
        print(f"\nErrors ({len(all_errors)}):")
        for combo_name, error in all_errors[:10]:
            print(f"  - {combo_name}: {error}")
            if len(error.split('\n')) > 1:
                print()
        if len(all_errors) > 10:
            print(f"  ... and {len(all_errors) - 10} more")
    
    if total_errors > 0:
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

