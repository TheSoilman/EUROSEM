"""End-to-end test with the REAL compiled model.
    cd eurosem_web && python3 test_real_model.py
Compiles ../headless/<stage> (default stage4_rval; pass another name as an
argument) with gfortran, then:
  1. runs the reference case through eurosem_run and compares every number
  2. breaks parameters that only the Fortran checks catch and confirms each
     gives a clear error, not a crash or a hang
  3. runs a case with ROC > 0 and PAVE > 0 (the questions answered by line 7
     of eurosem.ctl) and checks it completes"""
import os, sys, copy, subprocess, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import eurosem_io as E
from eurosem_run import run_case

STAGE = sys.argv[1] if len(sys.argv) > 1 else "stage4_rval"
stage = os.path.join(HERE, "..", "headless", STAGE)
print("Testing", STAGE)
exe = os.path.join(tempfile.mkdtemp(), "eurosem")
srcs = sorted(f for f in os.listdir(stage) if f.endswith(".for"))
r = subprocess.run(["gfortran", "-std=legacy", "-fno-automatic", "-fcheck=bounds", "-g", "-w", "-o", exe]
                   + [os.path.join(stage, f) for f in srcs], capture_output=True, text=True)
if r.returncode:
    sys.exit("COMPILE FAILED\n" + r.stderr[:1500])

REF = os.path.join(HERE, "reference")
base = E.read_pcp(open(f"{REF}/test.pcp").read(), E.read_par(open(f"{REF}/testn.par").read()))
fails = []
def check(name, cond, detail=""):
    print(("PASS  " if cond else "FAIL  ") + name + ("" if cond else "   " + str(detail)[:300]))
    if not cond: fails.append(name)

# 1. reference through the full Python path, compared with the frozen reference
out = run_case(base, exe)
check("reference case runs", out["ok"], out)
if out["ok"]:
    ref = E.load_results(open(f"{REF}/ref.dyn").read(), open(f"{REF}/ref.sta").read())
    new = out["results"]
    same = all(abs(new["elements"][j][k] - ref["elements"][j][k]) <= 1e-5 * max(1, abs(ref["elements"][j][k]))
               for j in ref["elements"] for k in ("sediment_total_kg", "peak_flow_mmh", "outflow_m3"))
    check("results match the frozen reference", same)
    check("global runoff 295.1976 m3", new["global"]["runoff_m3"] == 295.1976, new["global"])
    check("run time under 5 s", out["seconds"] < 5, out["seconds"])

# 2. errors the model's own checks must report (validator bypassed)
def broken(label, mutate, expect):
    r = copy.deepcopy(base); mutate(r)
    o = run_case(r, exe, check_inputs=False, timeout=30)
    check(label, (not o["ok"]) and o["stage"] == "run" and expect.lower() in o["error"].lower(), o)
broken("THI > THMX is caught by the model", lambda r: setattr(r.element(1), "THI", 0.35), "THI greater than THMX")
broken("porosity < THMX is caught", lambda r: setattr(r.element(1), "POR", 0.25), "Porosity")
broken("ROC out of range is caught", lambda r: setattr(r.element(1), "ROC", 1.5), "ROC")
broken("XL = 0 is caught", lambda r: setattr(r.element(2), "XL", 0), "XL")
broken("too many time steps is caught", lambda r: setattr(r, "DELT", 0.01), "TFIN/DELT")

# 3. cases that trigger the questions answered by eurosem.ctl
q = copy.deepcopy(base); q.element(1).ROC = 0.2; q.element(1).PAVE = 0.3
o = run_case(q, exe, accept="Y")
check("ROC/PAVE questions answered Y: runs", o["ok"], o)
o2 = run_case(q, exe, accept="N")
check("ROC/PAVE questions answered N: runs", o2["ok"], o2)
if o["ok"] and o2["ok"]:
    check("answer Y vs N changes the water results", o["results"]["global"]["runoff_m3"] != o2["results"]["global"]["runoff_m3"])

# 4. ROC on the FIRST element only must change the water results (needs stage 4)
r4 = copy.deepcopy(base); r4.element(1).ROC = 0.3
ya = run_case(r4, exe, accept="Y"); na = run_case(r4, exe, accept="N")
check("ROC on element 1 only: both answers run", ya["ok"] and na["ok"], (ya.get("error"), na.get("error")))
if ya["ok"] and na["ok"]:
    qy, qn = ya["results"]["global"]["runoff_m3"], na["results"]["global"]["runoff_m3"]
    check("ROC=0.3 on element 1 changes runoff (per-element ROC)", qy != qn, f"Y {qy}  N {qn}")
    check("accepting ROC reduces infiltration", ya["results"]["global"]["infiltration_planes_mm"] <
          na["results"]["global"]["infiltration_planes_mm"])
    print(f"      runoff Y {qy} m3, N {qn} m3")

# 5. the web form's default values give the reference result
from eurosem_form import build_run, DEFAULT_SPEC
frun, facc = build_run(DEFAULT_SPEC)
fo = run_case(frun, exe, accept=facc)
check("form defaults run", fo["ok"], fo)
if fo["ok"]:
    fe, re_ = fo["results"], E.load_results(open(f"{REF}/ref.dyn").read(), open(f"{REF}/ref.sta").read())
    check("form defaults reproduce the reference", all(
        abs(fe["elements"][j][k] - re_["elements"][j][k]) <= 1e-5 * max(1, abs(re_["elements"][j][k]))
        for j in re_["elements"] for k in ("sediment_total_kg", "peak_flow_mmh", "outflow_m3")) and
        fe["global"]["runoff_m3"] == re_["global"]["runoff_m3"], fe["global"])

print("\n%d failure(s)" % len(fails)); sys.exit(1 if fails else 0)
