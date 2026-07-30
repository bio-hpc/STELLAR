#!/usr/bin/env python3
"""Parallel AI-driven agentic loop for STELLAR using the Cursor SDK.

For each case it:
  1. Runs the STELLAR pipeline (run_case.py: stage -> steps 1..17 -> success).
  2. If it FAILS, spawns a Cursor agent (local runtime, Opus by default) that
     receives the future-cases playbook (CONTEXTO_FUTUROS_CASOS.md) as rules plus
     the failing log tail, diagnoses the root cause and applies a *targeted* fix
     respecting the project restrictions (no touching original case data, obrms
     -only rmsd_md, etc.).
  3. Re-runs the pipeline (resuming from the failing step when possible).
  4. Records the agent's reasoning + the actual code diff it produced.

Cases run in PARALLEL (thread pool). The AI *edit* phase is serialized with a
global lock so two agents never patch the shared STELLAR code at the same time;
the long pipeline/SLURM waits still overlap. (A fully-isolated alternative is one
git worktree per case -- see --help notes.)

At the end it writes a Markdown + PDF report of every modification the agent
considered and applied, under auto_fix_reports/RUN_<timestamp>/.

Usage
-----
  export CURSOR_API_KEY="crsr_3d186d74a388231f971e4d6787637efa2059c0f9d93472f54a2e3abcb5e99441"
  pip install cursor-sdk                       # into the 3.9 shim's interpreter
  agentic_loop_stellar/bin/python3 agentic_loop_stellar/auto_fix_loop.py \
      --workers 3 --max-fix-attempts 1 CASE1 CASE2 ...

  # all cases from a subset dir:
  CASES_DIR=/path/to/other_subset/cases \
  agentic_loop_stellar/bin/python3 agentic_loop_stellar/auto_fix_loop.py --workers 4

Model
-----
  --model opus-... forces a model id. If omitted the script queries
  Cursor.models.list() and auto-picks an Opus variant (Opus 4.8 preferred),
  falling back to "auto" if none is visible to the key.
"""

import argparse
import datetime as _dt
import os
import re
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CONTEXT_DOC = HERE / "CONTEXTO_FUTUROS_CASOS.md"
RUN_CASE = HERE / "run_case.py"
PY_SHIM = HERE / "bin" / "python3"

# Tracked source globs we attribute to the agent's edits (data dirs are untracked
# and never show up in `git diff`).
CODE_GLOBS = ["*.py", "*.sh", "*.cfg", "*.tex", "*.md", "*.mdp", "*.json"]

# Serializes the AI edit phase (snapshot + agent.prompt + diff capture) so
# concurrent agents don't clobber each other's changes to the shared codebase.
EDIT_LOCK = threading.Lock()
# Guards stdout so parallel workers don't interleave lines mid-write.
PRINT_LOCK = threading.Lock()


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def log(msg):
    with PRINT_LOCK:
        print(f"[{_dt.datetime.now():%H:%M:%S}] {msg}", flush=True)


def git(*args):
    """Run a git command at ROOT and return stripped stdout ('' on error)."""
    try:
        out = subprocess.run(
            ["git", "-C", str(ROOT), *args],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            universal_newlines=True,
        )
        return out.stdout if out.returncode == 0 else ""
    except Exception:  # noqa: BLE001
        return ""


def code_snapshot():
    """A commit-ish capturing tracked modifications *now*, without touching the
    working tree (so we can later diff the agent's edits against it)."""
    snap = git("stash", "create").strip()
    return snap or git("rev-parse", "HEAD").strip()


def code_diff_since(base):
    if not base:
        return ""
    return git("diff", base, "--", *CODE_GLOBS)


def tail_file(path, n=120):
    try:
        lines = Path(path).read_text(errors="replace").splitlines()
        return "\n".join(lines[-n:])
    except Exception:  # noqa: BLE001
        return "(no log)"


# --------------------------------------------------------------------------- #
# pipeline
# --------------------------------------------------------------------------- #
RESULT_RE = re.compile(r"RESULT\s+\S+:\s+(PASS|FAIL)\s+(.*)")
FAILSTEP_RE = re.compile(r"FAIL at step (\S+)")


def run_pipeline(case, knobs, from_step=None):
    """Invoke run_case.py for a single case. Returns (passed, detail, fail_step)."""
    cmd = [str(PY_SHIM), str(RUN_CASE), case,
           "--max-valid", str(knobs["max_valid"]),
           "--md-pool", str(knobs["md_pool"]),
           "--max-organize", str(knobs["max_organize"]),
           "--overlap-tol", str(knobs["overlap_tol"])]
    if from_step:
        cmd += ["--from-step", from_step]
    env = dict(os.environ)
    env.setdefault("CASES_DIR", str(ROOT / "improve5Frag_2" / "cases"))
    proc = subprocess.run(cmd, cwd=str(ROOT), env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          universal_newlines=True)
    out = proc.stdout or ""
    passed, detail, fail_step = False, "", None
    m = RESULT_RE.search(out)
    if m:
        passed = m.group(1) == "PASS"
        detail = m.group(2).strip()
    fs = FAILSTEP_RE.search(out)
    if fs:
        fail_step = fs.group(1)
    if not m:
        detail = f"driver rc={proc.returncode}"
    return passed, detail, fail_step


# --------------------------------------------------------------------------- #
# Cursor agent
# --------------------------------------------------------------------------- #
def pick_model(explicit, api_key):
    """Return a model id. Prefer explicit; else auto-detect Opus (4.8 first)."""
    if explicit:
        return explicit
    try:
        from cursor_sdk import Cursor  # noqa: WPS433
        models = Cursor.models.list()  # falls back to CURSOR_API_KEY
        ids = []
        for m in models:
            mid = getattr(m, "id", None) or (m.get("id") if isinstance(m, dict) else None) or str(m)
            ids.append(mid)
        opus = [i for i in ids if "opus" in i.lower()]
        # normalize dashes/dots/spaces so "claude-opus-4-8" matches "opus48"
        def _norm(s):
            return s.lower().replace("-", "").replace(".", "").replace(" ", "")
        # Opus 4.8 first (user preference), then any other Opus.
        for i in opus:
            if "opus48" in _norm(i):
                return i
        if opus:
            return sorted(opus)[-1]
        log(f"WARN: no Opus model visible to this key (saw: {ids}). Using 'auto'.")
    except Exception as e:  # noqa: BLE001
        log(f"WARN: could not list models ({e}); using 'auto'.")
    return "auto"


PROMPT_TEMPLATE = """You are fixing a STELLAR pipeline failure. Work in the repo at {root}.

## Project rules (authoritative playbook)
{context}

## Hard restrictions (do NOT violate)
- Do NOT modify original case data under improve5Frag_2/ (or $CASES_DIR). Staging is write-isolated.
- rmsd_md MUST come from obrms (robust backbone method). Never reintroduce a PyMOL `align` fallback.
- Prefer the smallest, most general fix. Fixes should generalize to future cases, not hardcode this one.
- Only edit tracked source (STELLAR/*.py, GROMACS/*, agentic_loop_stellar/run_case.py, etc.).

## Failing case
Case: {case}
Pipeline result: {detail}
Failing step: {fail_step}

### Tail of the master log (agentic_loop_stellar/logs/{case}/run.log)
{run_log}

### Tail of the failing step log
{step_log}

## Your task
1. Diagnose the ROOT cause of the failure for step {fail_step}.
2. Apply a targeted code fix consistent with the rules above.
3. Do NOT run the full pipeline yourself (it uses SLURM/MD and is slow); the harness re-runs it.

## Required output (end your reply with EXACTLY this block)
<<<SUMMARY
DIAGNOSIS: <1-3 sentences on the root cause>
CHANGES: <files touched and what you changed, one bullet per file>
GENERALIZES: <why this helps future cases, not just {case}>
RISK: <any risk / what to watch>
SUMMARY>>>
"""


def build_prompt(case, detail, fail_step, context):
    logdir = HERE / "logs" / case
    run_log = tail_file(logdir / "run.log", 80)
    step_log = ""
    if fail_step:
        step_log = tail_file(logdir / f"step_{fail_step}.log", 120)
    return PROMPT_TEMPLATE.format(
        root=ROOT, context=context, case=case, detail=detail,
        fail_step=fail_step or "?", run_log=run_log,
        step_log=step_log or "(no step log)")


def extract_summary(text):
    m = re.search(r"<<<SUMMARY\s*(.*?)\s*SUMMARY>>>", text or "", re.DOTALL)
    return (m.group(1).strip() if m else (text or "").strip())


def run_agent_fix(case, detail, fail_step, context, model, api_key, timeout_hint):
    """Run one Cursor agent to diagnose+fix. Returns dict with the record.

    The edit phase is serialized (EDIT_LOCK) so parallel agents don't collide on
    shared code. Returns even on failure (status recorded)."""
    from cursor_sdk import (Agent, AgentOptions, LocalAgentOptions,  # noqa: WPS433
                            LocalAgentStoreConfig, CursorAgentError)

    # The default local store uses node:sqlite (needs Node >= 22.13). The bundled
    # node we run on el7 is v20, so force a JSONL-backed store (no sqlite).
    store_dir = HERE / ".sdkenv" / "agent_store"
    store_dir.mkdir(parents=True, exist_ok=True)

    prompt = build_prompt(case, detail, fail_step, context)
    rec = {"case": case, "model": model, "status": "not_run",
           "agent_id": "", "run_id": "", "summary": "", "diff": "",
           "diff_stat": "", "error": ""}

    with EDIT_LOCK:
        base = code_snapshot()
        try:
            result = Agent.prompt(
                prompt,
                AgentOptions(
                    api_key=api_key,
                    model=model,
                    local=LocalAgentOptions(
                        cwd=str(ROOT),
                        store=LocalAgentStoreConfig(type="jsonl",
                                                    root_dir=str(store_dir)),
                    ),
                ),
            )
            rec["status"] = getattr(result, "status", "unknown")
            rec["run_id"] = getattr(result, "id", "") or ""
            rec["summary"] = extract_summary(getattr(result, "result", "") or "")
        except CursorAgentError as e:
            rec["status"] = "startup_error"
            rec["error"] = f"{getattr(e, 'message', e)} (retryable={getattr(e, 'is_retryable', '?')})"
            log(f"[{case}] agent startup error: {rec['error']}")
            return rec
        except Exception as e:  # noqa: BLE001
            rec["status"] = "exception"
            rec["error"] = str(e)
            log(f"[{case}] agent exception: {e}")
            return rec
        # capture what actually changed on disk
        rec["diff"] = code_diff_since(base)
        rec["diff_stat"] = git("diff", "--stat", base, "--", *CODE_GLOBS)
    return rec


# --------------------------------------------------------------------------- #
# per-case orchestration
# --------------------------------------------------------------------------- #
def process_case(case, knobs, context, model, api_key, max_fix_attempts, use_ai):
    record = {
        "case": case,
        "final": "FAIL",
        "detail": "",
        "attempts": [],   # list of agent fix records
        "history": [],     # list of (phase, passed, detail)
    }
    log(f"[{case}] pipeline run (initial)")
    passed, detail, fail_step = run_pipeline(case, knobs)
    record["history"].append(("initial", passed, detail))
    record["detail"] = detail
    if passed:
        record["final"] = "PASS"
        log(f"[{case}] PASS (no fix needed): {detail}")
        return record

    if not use_ai:
        log(f"[{case}] FAIL ({detail}); AI disabled (--no-ai).")
        return record

    for attempt in range(1, max_fix_attempts + 1):
        log(f"[{case}] FAIL at step {fail_step} ({detail}); agent fix attempt {attempt}/{max_fix_attempts}")
        fix = run_agent_fix(case, detail, fail_step, context, model, api_key,
                            timeout_hint=knobs)
        fix["attempt"] = attempt
        record["attempts"].append(fix)
        if fix["status"] in ("startup_error", "exception"):
            log(f"[{case}] agent could not run; stopping fix loop.")
            break
        # re-run: resume from the failing step when we have one (faster), else full
        resume = fail_step if fail_step and fail_step not in ("1", "2") else None
        log(f"[{case}] re-running pipeline (from_step={resume or 'start'})")
        passed, detail, fail_step = run_pipeline(case, knobs, from_step=resume)
        record["history"].append((f"after_fix_{attempt}", passed, detail))
        record["detail"] = detail
        if passed:
            record["final"] = "PASS"
            log(f"[{case}] PASS after fix attempt {attempt}: {detail}")
            return record

    log(f"[{case}] still FAIL after {max_fix_attempts} attempt(s): {detail}")
    return record


# --------------------------------------------------------------------------- #
# reports
# --------------------------------------------------------------------------- #
def _tex_escape(s):
    repl = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
            "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
            "^": r"\textasciicircum{}"}
    return "".join(repl.get(c, c) for c in s)


def _verbatim(s, cap=6000):
    s = (s or "").strip()
    if len(s) > cap:
        s = s[:cap] + "\n... (truncated)"
    return s


def write_reports(records, model, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    ts = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    npass = sum(1 for r in records if r["final"] == "PASS")
    nfail = len(records) - npass

    # ---- Markdown ----
    md = out_dir / "informe_auto_fix.md"
    with open(md, "w") as f:
        f.write(f"# STELLAR auto-fix loop — informe\n\n")
        f.write(f"- Fecha: {ts}\n- Modelo agente: `{model}`\n")
        f.write(f"- Casos: {len(records)} — PASS={npass}, FAIL={nfail}\n\n")
        f.write("## Resumen por caso\n\n")
        f.write("| Caso | Resultado | Intentos IA | Detalle |\n|---|---|---|---|\n")
        for r in records:
            f.write(f"| {r['case']} | {r['final']} | {len(r['attempts'])} | {r['detail']} |\n")
        f.write("\n## Modificaciones consideradas/aplicadas por el agente\n\n")
        for r in records:
            if not r["attempts"]:
                continue
            f.write(f"### {r['case']} ({r['final']})\n\n")
            for a in r["attempts"]:
                f.write(f"**Intento {a['attempt']}** — status `{a['status']}`"
                        f"{' — run `' + a['run_id'] + '`' if a['run_id'] else ''}\n\n")
                if a.get("error"):
                    f.write(f"> error: {a['error']}\n\n")
                if a.get("summary"):
                    f.write("_Razonamiento / cambios considerados:_\n\n")
                    f.write("```\n" + _verbatim(a["summary"]) + "\n```\n\n")
                if a.get("diff_stat"):
                    f.write("_Diff aplicado (stat):_\n\n")
                    f.write("```\n" + _verbatim(a["diff_stat"], 2000) + "\n```\n\n")
                if a.get("diff"):
                    patch = out_dir / f"{r['case']}_attempt{a['attempt']}.patch"
                    patch.write_text(a["diff"])
                    f.write(f"_Patch completo:_ `{patch.name}`\n\n")
    log(f"MD report -> {md}")

    # ---- LaTeX + PDF ----
    tex = out_dir / "informe_auto_fix.tex"
    rows = "\n".join(
        f"{_tex_escape(r['case'])} & {r['final']} & {len(r['attempts'])} & "
        f"{_tex_escape(r['detail'][:60])} \\\\"
        for r in records)
    sections = []
    for r in records:
        if not r["attempts"]:
            continue
        body = [f"\\subsection*{{{_tex_escape(r['case'])} ({r['final']})}}"]
        for a in r["attempts"]:
            body.append(f"\\textbf{{Intento {a['attempt']}}} — status "
                        f"\\texttt{{{_tex_escape(a['status'])}}}\\\\")
            if a.get("summary"):
                body.append("\\begin{lstlisting}")
                body.append(_verbatim(a["summary"], 4000))
                body.append("\\end{lstlisting}")
            if a.get("diff_stat"):
                body.append("\\begin{lstlisting}")
                body.append(_verbatim(a["diff_stat"], 1500))
                body.append("\\end{lstlisting}")
        sections.append("\n".join(body))
    tex_src = r"""\documentclass[10pt]{article}
\usepackage[a4paper,margin=2cm]{geometry}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{longtable}
\usepackage{listings}
\usepackage{xcolor}
\lstset{basicstyle=\ttfamily\footnotesize,breaklines=true,frame=single,
  columns=fullflexible,keepspaces=true,breakatwhitespace=false}
\title{STELLAR auto-fix loop --- informe}
\date{%(ts)s}
\begin{document}
\maketitle
\noindent Modelo agente: \texttt{%(model)s}\\
Casos: %(n)d --- PASS=%(npass)d, FAIL=%(nfail)d
\section*{Resumen por caso}
\begin{longtable}{llll}
\textbf{Caso} & \textbf{Resultado} & \textbf{Intentos IA} & \textbf{Detalle} \\
\hline
%(rows)s
\end{longtable}
\section*{Modificaciones consideradas/aplicadas}
%(sections)s
\end{document}
""" % {
        "ts": _tex_escape(ts), "model": _tex_escape(model), "n": len(records),
        "npass": npass, "nfail": nfail, "rows": rows,
        "sections": "\n\n".join(sections) if sections else "(sin intervenciones del agente)",
    }
    tex.write_text(tex_src)
    log(f"LaTeX -> {tex}")
    if _which("pdflatex"):
        for _ in range(2):
            subprocess.run(["pdflatex", "-interaction=nonstopmode",
                            "-halt-on-error", tex.name],
                           cwd=str(out_dir),
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        pdf = out_dir / "informe_auto_fix.pdf"
        if pdf.is_file():
            log(f"PDF report -> {pdf}")
        else:
            log("WARN: pdflatex did not produce a PDF (see .log).")
    else:
        log("WARN: pdflatex not found; skipping PDF.")


def _which(name):
    from shutil import which
    return which(name)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def discover_cases():
    cases_dir = Path(os.environ.get("CASES_DIR",
                                    str(ROOT / "improve5Frag_2" / "cases")))
    if not cases_dir.is_dir():
        return []
    return sorted(p.name for p in cases_dir.iterdir() if p.is_dir())


def main():
    ap = argparse.ArgumentParser(
        description="Parallel AI-driven STELLAR agentic loop (Cursor SDK).")
    ap.add_argument("cases", nargs="*",
                    help="cases to process (default: all under $CASES_DIR)")
    ap.add_argument("--workers", type=int, default=3,
                    help="parallel cases (default 3). SLURM/MD is the real limit.")
    ap.add_argument("--max-fix-attempts", type=int, default=1,
                    help="AI diagnose+fix+rerun attempts per failing case")
    ap.add_argument("--model", default=os.environ.get("STELLAR_FIX_MODEL"),
                    help="Cursor model id (default: auto-detect Opus)")
    ap.add_argument("--no-ai", action="store_true",
                    help="just run the pipeline in parallel, no agent fixes")
    # pilot knobs mirror run_loop.sh
    ap.add_argument("--max-valid", type=int, default=int(os.environ.get("MAX_VALID", 2)))
    ap.add_argument("--md-pool", type=int, default=int(os.environ.get("MD_POOL", 0)))
    ap.add_argument("--max-organize", type=int, default=int(os.environ.get("MAX_ORGANIZE", 400)))
    ap.add_argument("--overlap-tol", type=int, default=int(os.environ.get("OVERLAP_TOL", 450)))
    args = ap.parse_args()

    knobs = {"max_valid": args.max_valid, "md_pool": args.md_pool,
             "max_organize": args.max_organize, "overlap_tol": args.overlap_tol}

    cases = args.cases or discover_cases()
    if not cases:
        log("No cases found (set CASES_DIR or pass case names).")
        return 2

    api_key = os.environ.get("CURSOR_API_KEY", "")
    context = CONTEXT_DOC.read_text() if CONTEXT_DOC.is_file() else "(playbook missing)"

    model = "auto"
    if not args.no_ai:
        try:
            import cursor_sdk  # noqa: F401,WPS433
        except ImportError:
            log("ERROR: cursor-sdk not installed. `pip install cursor-sdk` into "
                f"{PY_SHIM}'s interpreter, or run with --no-ai.")
            return 2
        if not api_key:
            log("ERROR: CURSOR_API_KEY not set (needed for the agent). "
                "Export it or run with --no-ai.")
            return 2
        model = pick_model(args.model, api_key)
        log(f"Agent model: {model}")

    out_dir = HERE / "auto_fix_reports" / _dt.datetime.now().strftime("RUN_%Y%m%d_%H%M%S")
    log(f"Cases: {len(cases)} | workers: {args.workers} | AI: {not args.no_ai} "
        f"| report dir: {out_dir}")

    records = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(process_case, c, knobs, context, model, api_key,
                          args.max_fix_attempts, not args.no_ai): c
                for c in cases}
        for fut in as_completed(futs):
            c = futs[fut]
            try:
                records.append(fut.result())
            except Exception as e:  # noqa: BLE001
                log(f"[{c}] worker crashed: {e}")
                records.append({"case": c, "final": "FAIL",
                                "detail": f"worker crash: {e}",
                                "attempts": [], "history": []})

    records.sort(key=lambda r: r["case"])
    write_reports(records, model, out_dir)
    npass = sum(1 for r in records if r["final"] == "PASS")
    log(f"DONE. PASS={npass} FAIL={len(records) - npass}. Reports in {out_dir}")
    return 0 if npass == len(records) else 1


if __name__ == "__main__":
    sys.exit(main())
