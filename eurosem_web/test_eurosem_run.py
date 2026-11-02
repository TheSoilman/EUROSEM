"""Plumbing test for eurosem_run using a stand-in executable.
The stand-in is a small shell script that behaves like the headless
model: it needs eurosem.ctl, copies the reference output into the files
named there, and writes eurosem.status.  Run:  python3 test_eurosem_run.py"""
import os, sys, stat, tempfile, copy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import eurosem_io as E
from eurosem_run import run_case

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "reference")
fails = []
def check(name, cond, detail=""):
    print(("PASS  " if cond else "FAIL  ") + name + ("" if cond else "  " + str(detail)))
    if not cond: fails.append(name)

run = E.read_pcp(open(f"{REF}/test.pcp").read(), E.read_par(open(f"{REF}/testn.par").read()))
tmp = tempfile.mkdtemp()

def make_exe(name, body):
    p = os.path.join(tmp, name)
    open(p, "w").write("#!/bin/bash\n" + body)
    os.chmod(p, os.stat(p).st_mode | stat.S_IEXEC)
    return p

good = make_exe("good", f'''
[ -f eurosem.ctl ] || {{ echo "ERROR: no ctl" > eurosem.status; exit 2; }}
mapfile -t L < eurosem.ctl
[ "${{L[1]}}" = model.par ] && [ -s model.par ] && [ -s model.pcp ] || {{ echo "ERROR: inputs" > eurosem.status; exit 2; }}
cp {REF}/ref.dyn "${{L[3]}}"; cp {REF}/ref.aux "${{L[4]}}" 2>/dev/null; touch "${{L[4]}}"; cp {REF}/ref.sta "${{L[5]}}"
echo OK > eurosem.status
''')
bad = make_exe("bad", 'echo "ERROR: Porosity less than THMX on element 1" > eurosem.status; echo " FATAL ERROR" ; exit 2')
silent = make_exe("silent", 'echo hello; exit 0')
slow = make_exe("slow", 'sleep 5')

r = run_case(run, good)
check("good run ok", r["ok"], r)
if r["ok"]:
    el = r["results"]["elements"]
    check("results parsed", abs(el["2"]["outflow_m3"] - 295.2) < 1.0 and el["2"]["sediment_total_kg"] == 7878.425)
    check("global error", r["results"]["global"]["vol_error_pct"] == 0.0648)
check("exit status error is reported", (lambda r: (not r["ok"]) and r["stage"] == "run" and "Porosity" in r["error"])(run_case(run, bad)))
r = run_case(run, silent); check("missing status is an error", not r["ok"] and "without reporting" in r["error"], r)
r = run_case(run, slow, timeout=1); check("timeout is reported", not r["ok"] and "did not finish" in r["error"], r)
broken = copy.deepcopy(run); broken.element(1).THI = 0.9
r = run_case(broken, good); check("validation blocks the run", not r["ok"] and r["stage"] == "validation" and "THI" in r["error"], r)
nop = copy.deepcopy(run); nop.element(2).NPRINT = 0
seen = {}
r = run_case(nop, good, keep_dir=os.path.join(tmp, "kept"))
par = open(os.path.join(tmp, "kept", "model.par")).read()
check("NPRINT forced to 1 for every element", E.read_par(par).element(2).NPRINT == 1)
check("control file written", open(os.path.join(tmp, "kept", "eurosem.ctl")).read().split("\n")[1:7] ==
      ["model.par", "model.pcp", "model.dyn", "model.aux", "model.sta", "Y"])
print("\n%d failure(s)" % len(fails)); sys.exit(1 if fails else 0)
