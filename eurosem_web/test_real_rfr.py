"""RFR against the REAL model.  Run:  python3 test_real_rfr.py [path/to/eurosem]
RFR 0 must reproduce the reference.  For RFR > 0 the model's echoed surface storage should equal
0.001 x exp(-6.66 + 0.27 RFR) (the model prints it in metres although the label says mm), and the printed results show what depression storage and the Manning n
addition do to runoff and sediment."""
import os, sys, copy, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from eurosem_form import DEFAULT_SPEC, build_run
from eurosem_run import run_case
exe = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "eurosem")
fails = []
def ck(n, c, d=""):
    print(("PASS  " if c else "FAIL  ") + n + ("" if c else "   " + str(d)[:300]))
    if not c: fails.append(n)
def go(rfr):
    sp = copy.deepcopy(DEFAULT_SPEC)
    for p in sp["planes"]: p["rfr"] = rfr
    run, acc = build_run(sp); return run_case(run, exe, accept=acc)
for rfr in (0, 10, 20, 31.5):
    o = go(rfr)
    ck("RFR %g runs" % rfr, o["ok"], o.get("error"))
    if not o["ok"]: continue
    r = o["results"]; els = r["elements"]; e1 = els.get("1", els.get(1)); e2 = els.get("2", els.get(2))
    if rfr == 0:
        ck("RFR 0 equals the reference", abs(r["global"]["runoff_mm"] - 26.84) < 0.05 and abs(e2["sediment_total_kg"] - 7878.4) < 5)
    ss = e1["input"].get("SurfStor_mm"); want = 0.001 * math.exp(-6.66 + 0.27 * rfr)   # the model prints metres under an "mm" label
    ck("  echoed SurfStor %s (metres, printed with 6 decimals) matches 0.001*exp(-6.66+0.27 RFR) = %.6f" % (ss, want), ss is not None and abs(ss - want) < 1.5e-6, ss)
    print("   RFR %-5g runoff %.3f mm  outlet sediment %.1f kg  plane 1: time to runoff %s min, peak %s mm/h" % (
        rfr, r["global"]["runoff_mm"], e2["sediment_total_kg"], e1["hydrology"].get("time_to_runoff_min"), e1["hydrology"].get("peak_flow_mmh")))
print("\n%d failure(s)" % len(fails)); sys.exit(1 if fails else 0)
