"""
eurosem_run.py  --  run the headless EUROSEM executable for one case.

    from eurosem_io import Run, read_par, read_pcp
    from eurosem_run import run_case
    out = run_case(run, "/path/to/eurosem")

What it does
  1. validates the run (refuses to start if validate() reports errors)
  2. makes a fresh empty folder and writes model.par, model.pcp and the
     control file eurosem.ctl (stage 2 of the Fortran clean-up reads it
     and skips every keyboard question)
  3. runs the program with no keyboard, with a time limit
  4. reads eurosem.status ("OK" or "ERROR: message"), and if OK parses
     model.dyn / model.sta into the results structure used by the page

It returns a plain dict, ready to send to the browser as JSON:
  {"ok": True,  "results": {...}, "issues": [...], "seconds": 0.4}
  {"ok": False, "stage": "validation"|"run"|"output", "error": "...",
   "issues": [...], "log": "last lines of screen output",
   "aux": "last lines of the diagnostic file"}
"""
import copy
import os
import shutil
import subprocess
import tempfile
import time
from typing import Optional

import eurosem_io as E

TITLE_MAX = 80


def _tail(text: str, n: int = 25) -> str:
    return "\n".join(text.replace("\r", "").rstrip().split("\n")[-n:])


def run_case(run: E.Run, exe: str, timeout: float = 60.0,
             keep_dir: Optional[str] = None, accept: str = "Y",
             check_inputs: bool = True) -> dict:
    """Run one case.  `accept` answers the ROC / PAVE / SPLTEX questions
    (Y = accept the adjustment).  `check_inputs=False` skips validate()
    so the model's own checks can be tested.  `keep_dir`: copy the working folder
    there afterwards (for debugging)."""
    run = copy.deepcopy(run)
    # The results page needs the hydrograph of every element.
    for e in run.elements:
        if e.NPRINT <= 0:
            e.NPRINT = 1
    issues = E.validate(run)
    iss = [str(i) for i in issues]
    # check_inputs=False is for testing the model's own error handling.
    if check_inputs and E.has_errors(issues):
        return {"ok": False, "stage": "validation", "issues": iss,
                "error": "; ".join(s for s in iss if s.startswith("[ERROR]"))}

    work = tempfile.mkdtemp(prefix="eurosem_")
    try:
        title = " ".join(run.title.split())[:TITLE_MAX] or "EUROSEM run"
        files = {"model.par": E.write_par(run), "model.pcp": E.write_pcp(run),
                 "eurosem.ctl": "\n".join([title, "model.par", "model.pcp", "model.dyn",
                                           "model.aux", "model.sta", accept[:1]]) + "\n"}
        for name, text in files.items():
            with open(os.path.join(work, name), "w", newline="\n") as f:
                f.write(text)

        t0 = time.time()
        try:
            p = subprocess.run([os.path.abspath(exe)], cwd=work, stdin=subprocess.DEVNULL,
                               capture_output=True, text=True, timeout=timeout, errors="replace")
        except subprocess.TimeoutExpired:
            return {"ok": False, "stage": "run", "issues": iss,
                    "error": f"model did not finish within {timeout:.0f} s"}
        seconds = time.time() - t0

        def slurp(name):
            try:
                with open(os.path.join(work, name), errors="replace") as f:
                    return f.read()
            except OSError:
                return ""

        status = slurp("eurosem.status").strip()
        log = _tail(p.stdout + p.stderr)
        if status != "OK":
            msg = status[len("ERROR: "):] if status.startswith("ERROR: ") else (
                "model ended without reporting a status (exit code %d)" % p.returncode)
            if "excess rill fill" in msg.lower():
                msg += (". Sediment deposited in a rill is more than the rill can hold, so the model stops. "
                        "Make the rills larger (fewer, wider or deeper), or check cohesion, erodibility and "
                        "particle size, which set how much sediment reaches them.")
            return {"ok": False, "stage": "run", "issues": iss, "error": msg,
                    "log": log, "aux": _tail(slurp("model.aux"))}

        dyn, sta = slurp("model.dyn"), slurp("model.sta")
        if not dyn or not sta:
            return {"ok": False, "stage": "output", "issues": iss, "log": log,
                    "error": "model reported OK but an output file is missing or empty"}
        try:
            results = E.load_results(dyn, sta)
        except Exception as exc:            # unexpected layout
            return {"ok": False, "stage": "output", "issues": iss, "log": log,
                    "error": f"could not read the model output: {exc}"}
        if not results["elements"]:
            return {"ok": False, "stage": "output", "issues": iss, "log": log,
                    "error": "no element results found in the output files"}
        return {"ok": True, "results": results, "issues": iss, "seconds": round(seconds, 2)}
    finally:
        if keep_dir:
            shutil.copytree(work, keep_dir, dirs_exist_ok=True)
        shutil.rmtree(work, ignore_errors=True)
