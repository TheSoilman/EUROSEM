"""Tests for eurosem_export.  Run: python3 test_eurosem_export.py"""
import os, sys, csv, io
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import eurosem_io as E
from eurosem_form import DEFAULT_SPEC
from eurosem_export import make_results_csv, make_settings_csv
REF = os.path.join(HERE, "reference")
fails = []
def check(name, cond, detail=""):
    print(("PASS  " if cond else "FAIL  ") + name + ("" if cond else "   " + str(detail)[:300]))
    if not cond: fails.append(name)
res = E.load_results(open(f"{REF}/ref.dyn").read(), open(f"{REF}/ref.sta").read())
text = make_results_csv(DEFAULT_SPEC, res)
raw = text.split("\n")
kt = next(i for i, l in enumerate(raw) if l.startswith("time_min"))
check("every header line begins with '#' in the raw text", all(l.startswith("#") for l in raw[:kt]), [l[:40] for l in raw[:kt] if not l.startswith("#")])
check("no byte-order mark in the results file", not text.startswith("\ufeff"))
rows = list(csv.reader(io.StringIO(text)))
head = [r for r in rows if r and r[0].startswith("# ")]
k = next(i for i, r in enumerate(rows) if r and r[0] == "time_min")
table = rows[k:]
check("all header rows come before the table", all(not (r and r[0].startswith("# ")) for r in table))
check("header rows are all comments", len(head) == k, (len(head), k))
check("table: 177 time rows + header", len(table) == 178, len(table))
check("columns: 2 + 5 per plane", len(table[0]) == 12 and table[0][:3] == ["time_min", "rain_mmh", "plane1_q_mmh"], table[0])
check("every row has the same width", all(len(r) == 12 for r in table))
col = {n: i for i, n in enumerate(table[0])}
last = table[-1]
check("last time is 88", float(last[0]) == 88.0)
peak = max(table[1:], key=lambda r: float(r[col["plane2_q_mmh"]]))
check("plane 2 peak 77.287 mm/h at 30 min", float(peak[col["plane2_q_mmh"]]) == 77.287 and float(peak[0]) == 30.0)
check("plane 1 peak sediment 773.4 kg/min", abs(max(float(r[col["plane1_sed_kgmin"]]) for r in table[1:]) - 773.4) < 0.1)
check("cumulative sediment ends near 7878 kg (plane 2)", abs(float(last[col["plane2_cum_sed_kg"]]) - 7878.4) < 40, last[col["plane2_cum_sed_kg"]])
flat = "\n".join(",".join(r) for r in head)
check("header records the parameters", "Manning n" in flat and "0.095" in flat and "Plant height" in flat and "40" in flat)
check("header records the rainfall table", "Rainfall: time (min)" in flat and "50" in flat)
check("header records summary", "Sediment out" in flat and "7878.425" in flat and "295.1976" in flat)
check("header names the model", "LOG10" in flat)
s = list(csv.reader(io.StringIO(make_settings_csv(DEFAULT_SPEC).lstrip("\ufeff"))))
check("settings file: plane columns", ["Plane parameter", "unit", "plane 1", "plane 2"] in s)
check("settings file: slope row", ["Slope", "m/m", "0.095", "0.05"] in s, [r for r in s if r and r[0] == "Slope"])
# a title with commas and quotes survives
sp = dict(DEFAULT_SPEC, title='Field "A", north')
s2 = list(csv.reader(io.StringIO(make_settings_csv(sp).lstrip("\ufeff"))))
check("title with comma and quotes is kept", s2[0] == ["Run title", "", 'Field "A", north'], s2[0])

# --- guide choices are recorded ---
import eurosem_guide as G
sp = dict(DEFAULT_SPEC); sp["planes"] = [dict(pl) for pl in DEFAULT_SPEC["planes"]]
sp["planes"][0].update(texture="Loam", condition="loose_dry", **G.guide_fill("Loam", "loose_dry"))
sp["planes"][1].update(texture="Loam", condition="compacted", **G.guide_fill("Loam", "compacted")); sp["planes"][1]["fmin"] = 99
g = list(csv.reader(io.StringIO(make_settings_csv(sp).lstrip("\ufeff"))))
row = lambda name: next(r for r in g if r and r[0] == name)
check("texture recorded per plane", row("Soil texture chosen (guide)")[2:] == ["Loam", "Loam"], row("Soil texture chosen (guide)"))
check("surface state recorded", row("Surface state chosen")[2:] == ["Loose and dry", "Sealed or compacted"])
cmp_ = row("Values compared with the guide fill")[2:]
check("untouched plane is 'as filled'", cmp_[0] == "as filled", cmp_)
check("edited value is named", cmp_[1] == "edited after filling: conductivity", cmp_)
check("hand-entered planes say so", row("Soil texture chosen (guide)")[2:] == ["Loam", "Loam"] and list(csv.reader(io.StringIO(make_settings_csv(DEFAULT_SPEC).lstrip("\ufeff"))))[-1][0] != "Guide source")
bad = dict(DEFAULT_SPEC); bad["planes"] = [dict(DEFAULT_SPEC["planes"][0], texture="=cmd()"), dict(DEFAULT_SPEC["planes"][1])]
gb = list(csv.reader(io.StringIO(make_settings_csv(bad).lstrip("\ufeff"))))
check("unknown texture text is never written", "=cmd()" not in make_settings_csv(bad) and next(r for r in gb if r and r[0] == "Soil texture chosen (guide)")[2] == "none")
gsp = list(csv.reader(io.StringIO(make_results_csv(sp, res))))
check("results header also carries the guide rows", any(r and r[0] == "# Soil texture chosen (guide)" for r in gsp))
check("results table still intact", next(i for i, r in enumerate(gsp) if r and r[0] == "time_min") > 0 and len([r for r in gsp if r and r[0] == "time_min"]) == 1)
try:
    import pandas as pd
    df = pd.read_csv(io.StringIO(text), comment="#")
    check("pandas reads it with comment='#'", df.shape == (177, 12) and df["plane2_q_mmh"].max() == 77.287, df.shape)
except ImportError:
    print("      (pandas not installed: skipped the pandas check)")
print("\n%d failure(s)" % len(fails)); sys.exit(1 if fails else 0)
