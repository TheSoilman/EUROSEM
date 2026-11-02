"""Compare a new EUROSEM run with the frozen reference.
Usage:  python3 compare_runs.py new.dyn new.sta [ref.dyn ref.sta]
Defaults to reference/ref.dyn and reference/ref.sta.
Exit code 0 = same within tolerance, 1 = something moved.

Tolerances: the model prints 4-6 significant figures, so a rewrite that
changes nothing should reproduce every printed number exactly.  Small
allowances cover compiler/optimisation differences."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import eurosem_io as E

here = os.path.dirname(os.path.abspath(__file__))
a = sys.argv[1:]
if len(a) not in (2, 4):
    sys.exit(__doc__)
ref = (a[2], a[3]) if len(a) == 4 else (os.path.join(here, "reference/ref.dyn"), os.path.join(here, "reference/ref.sta"))
new = E.load_results(open(a[0], errors="replace").read(), open(a[1], errors="replace").read())
old = E.load_results(open(ref[0], errors="replace").read(), open(ref[1], errors="replace").read())

REL = 1e-5     # relative tolerance on printed values
ABS = 1e-6
bad = []

def cmp(label, x, y):
    if x is None or y is None:
        if x != y: bad.append(f"{label}: {x} vs reference {y}")
        return
    if abs(x - y) > ABS + REL * max(abs(x), abs(y)):
        bad.append(f"{label}: {x} vs reference {y}")

if set(new["elements"]) != set(old["elements"]):
    bad.append(f"elements differ: {sorted(new['elements'])} vs {sorted(old['elements'])}")
for j in sorted(set(old["elements"]) & set(new["elements"])):
    o, n = old["elements"][j], new["elements"][j]
    for k in ("vol_error_pct", "sediment_total_kg", "peak_flow_mmh", "time_to_peak_min", "contributing_area_ha"):
        cmp(f"el{j} {k}", n.get(k), o.get(k))
    for grp in ("erosion", "hydrology"):
        for k, v in o.get(grp, {}).items():
            cmp(f"el{j} {grp}.{k}", n.get(grp, {}).get(k), v)
    ho, hn = o.get("hydrograph"), n.get("hydrograph")
    if (ho is None) != (hn is None):
        bad.append(f"el{j}: hydrograph present in only one run")
    elif ho:
        if len(ho["time_min"]) != len(hn["time_min"]):
            bad.append(f"el{j}: {len(hn['time_min'])} rows vs reference {len(ho['time_min'])}")
        else:
            for col in ("rain_mmh", "q_m3min", "q_mmh", "conc", "qs_kgmin"):
                worst = max(((abs(x - y) - ABS - REL * max(abs(x), abs(y)), i) for i, (x, y) in enumerate(zip(hn[col], ho[col]))))
                if worst[0] > 0:
                    i = worst[1]
                    bad.append(f"el{j} hydrograph {col} at t={ho['time_min'][i]}: {hn[col][i]} vs reference {ho[col][i]}")
for k, v in old["global"].items():
    cmp(f"global {k}", new["global"].get(k), v)

if bad:
    print(f"DIFFERENT: {len(bad)} value(s) moved")
    for b in bad[:25]: print("  ", b)
    sys.exit(1)
print("SAME: every printed number matches the reference")
