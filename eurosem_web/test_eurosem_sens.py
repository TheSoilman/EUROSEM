"""Tests for eurosem_sens with a fake model.  Run: python3 test_eurosem_sens.py"""
import os, sys, csv, io, copy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import eurosem_io as E, eurosem_sens as S, eurosem_guide as G
from eurosem_form import DEFAULT_SPEC
REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference")
base_res = E.load_results(open(f"{REF}/ref.dyn").read(), open(f"{REF}/ref.sta").read())
f = []
def ck(n, c, d=""):
    print(("PASS  " if c else "FAIL  ") + n + ("" if c else "   " + str(d)[:300]))
    if not c: f.append(n)

calls = []
def fake(spec):
    """pretend runoff falls as conductivity rises and sediment rises with erodibility"""
    calls.append(copy.deepcopy(spec))
    pl = spec["planes"]
    for p in pl:                      # same rules as the real validator
        if not (0 <= float(p["theta_i"]) <= float(p["theta_max"]) <= float(p["porosity"])):
            return False, "need 0 <= THI <= THMX <= POR"
    r = copy.deepcopy(base_res)
    k = 5.0 / float(pl[0]["fmin"])
    r["global"]["runoff_mm"] *= k
    o = r["elements"]["2"]; o["peak_flow_mmh"] *= k
    o["sediment_total_kg"] *= float(pl[0]["erodibility"]) / 0.9
    return True, r

spec = copy.deepcopy(DEFAULT_SPEC)
tex = lambda t, c="compacted": dict(texture=t, condition=c, **G.guide_fill(t, c))
spec["planes"][0].update(tex("Loam")); spec["planes"][1].update(tex("Loam"))
pl = S.plan(spec)
fm = next(x for x in pl if x["key"] == "fmin")
ck("plan: conductivity range from A4.2", fm["guide"][0]["low"] == 2 and fm["guide"][0]["high"] == 65)
ck("plan: base values listed per plane", fm["base"] == [spec["planes"][0]["fmin"], spec["planes"][1]["fmin"]] and fm["base"][0] == 13, fm["base"])
ck("plan: manning has no guide range", not next(x for x in pl if x["key"] == "manning")["has_guide"])
ck("plan: cohesion uses compacted values", next(x for x in pl if x["key"] == "cohesion")["guide"][0]["low"] == 7)
sp2 = copy.deepcopy(spec); sp2["planes"][0]["condition"] = "loose_dry"
ck("plan: cohesion uncompacted when loose", S.guide_range(sp2["planes"][0], "cohesion")[:2] == (2, 4))
ck("plan: anomalous cohesion row gives no range", S.guide_range({"texture": "Silty clay"}, "cohesion") is None)
ck("plan: no texture gives no range", S.guide_range({}, "fmin") is None)

out = S.run_sensitivity(spec, ["fmin", "erodibility", "manning", "porosity"], {}, fake)
ck("ran ok", out["ok"], out)
rows = {r["key"]: r for r in out["rows"]}
ck("manning skipped with a reason", [s["key"] for s in out["skipped"]] == ["manning"] and "enter" in out["skipped"][0]["reason"], out["skipped"])
ck("runs = 1 base + 2 per parameter", out["runs"] == 7 and len(calls) == 7, (out["runs"], len(calls)))
ck("fmin low/high applied to both planes", calls[1]["planes"][0]["fmin"] == 2 and calls[1]["planes"][1]["fmin"] == 2 and calls[2]["planes"][1]["fmin"] == 65)
ck("other parameters untouched in a run", calls[1]["planes"][0]["g"] == spec["planes"][0]["g"] and calls[1]["planes"][0]["erodibility"] == spec["planes"][0]["erodibility"])
ck("base spec not modified", spec["planes"][0]["fmin"] == 13)
lo, hi = rows["fmin"]["low"]["metrics"]["runoff_mm"], rows["fmin"]["high"]["metrics"]["runoff_mm"]
ck("lower conductivity gives more runoff (fake model)", lo > out["base"]["runoff_mm"] > hi or lo > hi)
ck("erodibility range is loose-dry to loose-moist (1.0 to 2.7)", rows["erodibility"]["low"]["values"] == [1.0, 1.0] and rows["erodibility"]["high"]["values"] == [2.7, 2.7])
por = rows["porosity"]
ck("low porosity run is valid: moisture lowered and noted", por["low"]["error"] is None and por["low"]["notes"] and por["low"]["values"][0] == 0.38, por["low"])
lowpor = next(c for c in calls if c["planes"][0]["porosity"] == 0.38)
ck("capped moisture really is <= porosity", float(lowpor["planes"][0]["theta_max"]) <= 0.38 and float(lowpor["planes"][0]["theta_i"]) <= float(lowpor["planes"][0]["theta_max"]), lowpor["planes"][0])
# user override replaces guide ranges and works for parameters with none
o2 = S.run_sensitivity(spec, ["manning", "fmin"], {"manning": {"low": 0.1, "high": 0.3}, "fmin": {"low": 1, "high": 5}}, fake)
ck("override accepted for manning", any(r["key"] == "manning" and r["low"]["values"] == [0.1, 0.1] and r["source"] == "entered by you" for r in o2["rows"]), o2["skipped"])
ck("override beats guide", next(r for r in o2["rows"] if r["key"] == "fmin")["high"]["values"] == [5, 5])
o3 = S.run_sensitivity(spec, ["manning"], {"manning": {"low": 0.3, "high": 0.1}}, fake)
ck("reversed low/high are swapped", o3["rows"][0]["low"]["values"] == [0.1, 0.1])
# a plane with no texture keeps its value
sp3 = copy.deepcopy(spec); sp3["planes"][1].pop("texture"); calls.clear()
S.run_sensitivity(sp3, ["fmin"], {}, fake)
ck("plane without a range is left alone", calls[1]["planes"][1]["fmin"] == sp3["planes"][1]["fmin"] and calls[1]["planes"][0]["fmin"] == 2)
# failing run is recorded, not fatal
def flaky(s):
    return (False, "model stopped") if float(s["planes"][0]["fmin"]) == 65 else fake(s)
o4 = S.run_sensitivity(spec, ["fmin"], {}, flaky)
ck("a failed run is reported", o4["rows"][0]["high"]["error"] == "model stopped" and o4["rows"][0]["high"]["metrics"] is None)
ck("failed base run stops everything", not S.run_sensitivity(spec, ["fmin"], {}, lambda s: (False, "no"))["ok"])
many = S.run_sensitivity(spec, [p[0] for p in S.PARAMS], {k: {"low": 0.01, "high": 0.02} for k, *_ in S.PARAMS}, fake)
ck("run cap respected", many["runs"] <= S.MAX_RUNS, many["runs"])
text = S.sensitivity_csv(spec, out)
rows_ = list(csv.reader(io.StringIO(text)))
k = next(i for i, r in enumerate(rows_) if r and r[0] == "parameter")
ck("csv: header rows are comments", all(r[0].startswith("# ") for r in rows_[:k] if r), rows_[:2])
tab = rows_[k:]
ck("csv: 2 rows per parameter + skipped", len(tab) == 1 + 6 + 1, len(tab))
ck("csv: change from base computed", any(r[2] == "low" and r[0] == "Conductivity (Ks)" and float(r[10]) > 0 for r in tab[1:]))
ck("csv: records the guide choices", "Soil texture chosen (guide)" in text)
try:
    import pandas as pd; df = pd.read_csv(io.StringIO(text), comment="#"); ck("csv: pandas reads it", df.shape[0] == 7, df.shape)
except ImportError: pass
print("\n%d failure(s)" % len(f)); sys.exit(1 if f else 0)
