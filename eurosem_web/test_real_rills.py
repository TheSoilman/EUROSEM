"""Rills against the REAL model.  Run:  python3 test_real_rills.py [path/to/eurosem]
Runs the reference two-plane case with rills off (must equal the reference) and then with rills
on plane 1 at several settings, printing what changes.  The checks are loose on purpose: they
confirm the run completes, the water balance closes and the rill erosion line is read; the
printed numbers are for you to judge against what you know the model should do."""
import os, sys, copy
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from eurosem_form import DEFAULT_SPEC, build_run
from eurosem_run import run_case
exe = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "eurosem")
fails = []
def ck(n, c, d=""):
    print(("PASS  " if c else "FAIL  ") + n + ("" if c else "   " + str(d)[:400]))
    if not c: fails.append(n)
def go(**kw):
    sp = copy.deepcopy(DEFAULT_SPEC); sp["planes"][0].update(kw)
    run, acc = build_run(sp)
    return run_case(run, exe, accept=acc)
def show(tag, o):
    r = o["results"]; e1 = r["elements"]["1"] if "1" in r["elements"] else r["elements"][1]
    e2 = r["elements"]["2"] if "2" in r["elements"] else r["elements"][2]
    er = e1.get("erosion", {})
    print("   %-34s runoff %.3f mm  outlet sed %9.1f kg  plane1: rill %s kg, interrill %s kg, net %s kg  vol err %s %%" % (
        tag, r["global"]["runoff_mm"], e2.get("sediment_total_kg", float("nan")), er.get("gross_rill_kg", "-"),
        er.get("gross_interrill_kg", "-"), er.get("net_kg", "-"), e1.get("vol_error_pct")))
off = go()
ck("rills off runs", off["ok"], off.get("error"))
if off["ok"]:
    ck("rills off equals the reference (runoff 26.84 mm, sediment 7878 kg)", abs(off["results"]["global"]["runoff_mm"] - 26.84) < 0.05, off["results"]["global"])
    show("off (reference)", off)
cases = [("20 rills, 0.3 m wide, 0.15 m deep (default)", dict(rills=1)),
         ("default, depth grows downslope", dict(rills=1, rill_uniform=0)),
         ("10 rills, 0.5 m wide, 0.25 m deep", dict(rills=1, rill_number=10, rill_width=0.5, rill_depth=0.25)),
         ("default, rill Manning n 0.05", dict(rills=1, rill_manning=0.05)),
         ("default, interrill slope 0.2", dict(rills=1, interrill_slope=0.2))]
for tag, kw in cases:
    o = go(**kw)
    ck("runs: " + tag, o["ok"], o.get("error"))
    if not o["ok"]: continue
    show(tag, o)
    r = o["results"]; e1 = r["elements"].get("1", r["elements"].get(1))
    ck("  water balance closes (< 1 %)", abs(r["global"].get("vol_error_pct", 0) or 0) < 1, r["global"])
    ck("  rill erosion line was read", "gross_rill_kg" in e1.get("erosion", {}), e1.get("erosion"))
    ck("  rill geometry was echoed by the model", e1["input"].get("DEPNO") is not None, e1["input"])
print("\n%d failure(s)" % len(fails)); sys.exit(1 if fails else 0)
