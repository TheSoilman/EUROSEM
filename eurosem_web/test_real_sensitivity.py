"""Sensitivity test against the REAL model.  Run:  python3 test_real_sensitivity.py [path/to/eurosem]
Uses the reference two-plane case (values unchanged) with the texture 'Loam' chosen on both planes,
so the ranges come from the guide but the base case is still the reference run."""
import os, sys, copy
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from eurosem_form import DEFAULT_SPEC, build_run
from eurosem_run import run_case
import eurosem_sens as S
exe = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "eurosem")
fails = []
def ck(n, c, d=""):
    print(("PASS  " if c else "FAIL  ") + n + ("" if c else "   " + str(d)[:300]))
    if not c: fails.append(n)
def run_one(sp):
    run, acc = build_run(sp)
    o = run_case(run, exe, accept=acc)
    return (True, o["results"]) if o["ok"] else (False, o.get("error"))
spec = copy.deepcopy(DEFAULT_SPEC)
for p in spec["planes"]: p.update(texture="Loam", condition="compacted")
out = S.run_sensitivity(spec, ["fmin", "g", "porosity", "erodibility", "cohesion"], {}, run_one)
ck("sensitivity ran", out["ok"], out)
if out["ok"]:
    rows = {r["key"]: r for r in out["rows"]}
    ck("all five parameters ran, no failed run", len(rows) == 5 and all(r[c]["error"] is None for r in rows.values() for c in ("low", "high")),
       [(r["key"], c, r[c]["error"]) for r in rows.values() for c in ("low", "high") if r[c]["error"]])
    ck("base case is the reference run (runoff 26.84 mm, 295.2 m3 outlet)", abs(out["base"]["runoff_mm"] - 26.84) < 0.05 and abs(out["base"]["sediment_kg"] - 7878.4) < 5, out["base"])
    f = rows["fmin"]; ck("lower conductivity gives more runoff than higher", f["low"]["metrics"]["runoff_mm"] > f["high"]["metrics"]["runoff_mm"], (f["low"]["metrics"], f["high"]["metrics"]))
    e = rows["erodibility"]; ck("higher erodibility gives no less sediment", e["high"]["metrics"]["sediment_kg"] >= e["low"]["metrics"]["sediment_kg"] - 1e-6, (e["low"]["metrics"]["sediment_kg"], e["high"]["metrics"]["sediment_kg"]))
    p = rows["porosity"]; ck("porosity runs valid (moisture limits adjusted if needed)", p["low"]["metrics"] is not None and p["high"]["metrics"] is not None, p)
    print("\n  change from base, %:")
    b = out["base"]
    for k, r in rows.items():
        c = lambda case, m: (r[case]["metrics"][m] - b[m]) / b[m] * 100
        print("   %-12s runoff %+7.1f/%+7.1f   peak %+7.1f/%+7.1f   sediment %+7.1f/%+7.1f" % (k, c("low","runoff_mm"), c("high","runoff_mm"), c("low","peak_flow_mmh"), c("high","peak_flow_mmh"), c("low","sediment_kg"), c("high","sediment_kg")))
print("\n%d failure(s)" % len(fails)); sys.exit(1 if fails else 0)
