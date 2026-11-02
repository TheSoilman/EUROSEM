"""Tests for eurosem_io against the baseline run.
Run:  python3 test_eurosem_io.py [folder-with-baseline-files]
Baseline files expected: testn.par, test.pcp, test26.dyn, test26.sta"""
import sys, os, copy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import eurosem_io as E

D = sys.argv[1] if len(sys.argv) > 1 else "."
rd = lambda n: open(os.path.join(D, n), errors="replace").read()
fails = []


def check(name, cond, detail=""):
    print(("PASS  " if cond else "FAIL  ") + name + (("  " + detail) if detail and not cond else ""))
    if not cond:
        fails.append(name)


def close(a, b, tol):
    return abs(a - b) <= tol


# ---- 1. read the baseline inputs ------------------------------------
run = E.read_par(rd("testn.par"))
run = E.read_pcp(rd("test.pcp"), run)
check("par: 2 elements, order 1,2", run.NELE == 2 and run.order == [1, 2])
check("par: system record", (run.CLEN, run.TFIN, run.DELT, run.THETA, run.TEMP) == (110, 88, 0.5, 0.7, 20))
check("par: NTIME/NEROS", (run.NTIME, run.NEROS) == (2, 2))
e1, e2 = run.element(1), run.element(2)
check("par: element 1 values", (e1.XL, e1.W, e1.S, e1.SIR, e1.FMIN, e1.G, e1.D50, e1.NPRINT) ==
      (100, 100, 0.085, 0.095, 5, 400, 25, 7))
check("par: element 2 values", (e2.NU, e2.XL, e2.S, e2.SIR, e2.NPRINT) == (1, 10, 0.05, 0.05, 0))
check("pcp: 1 gage, weights", len(run.gages) == 1 and run.weights == {1: (1, 1.0), 2: (1, 1.0)})
check("pcp: only first ND=2 pairs kept", run.gages[0].times == [0, 30] and run.gages[0].depths == [0, 50])
check("pcp: gage id", run.gages[0].ident == "eurosem test data", run.gages[0].ident)
check("effective slope: no rills -> SIR", e1.effective_slope == 0.095)

# ---- 2. round trip: write -> read gives identical objects -------------
run2 = E.read_pcp(E.write_pcp(run), E.read_par(E.write_par(run)))
check("round trip par+pcp", run2.elements == run.elements and run2.order == run.order and
      run2.gages == run.gages and run2.weights == run.weights and
      (run2.TFIN, run2.DELT, run2.NTIME) == (run.TFIN, run.DELT, run.NTIME))

# ---- 3. written files have the same record positions as baseline -----
def record_lines(text):
    """line numbers (1-based) of lines that are purely numeric"""
    out = []
    for n, ln in enumerate(text.replace("\r", "").split("\n"), 1):
        t = ln.replace(",", " ").split()
        if t:
            try:
                [float(x) for x in t]; out.append(n)
            except ValueError:
                pass
    return out

check("par: numeric records on same lines as baseline",
      record_lines(E.write_par(run)) == record_lines(rd("testn.par")))
# pcp baseline has 3 pairs listed with ND=2 -> compare up to the 2nd pair
base_pcp = record_lines(rd("test.pcp"))
new_pcp = record_lines(E.write_pcp(run))
check("pcp: numeric records on same lines as baseline (first 2 pairs)", new_pcp == base_pcp[:len(new_pcp)] and len(new_pcp) == len(base_pcp) - 1,
      f"{new_pcp} vs {base_pcp}")

# ---- 4. validation ------------------------------------------------------
iss = E.validate(run)
errs = [i for i in iss if i.level == "error"]
check("baseline passes validation (no errors)", not errs, "; ".join(map(str, errs)))
check("baseline: note that rain stops at 30 min", any(i.level=="note" and "before TFIN" in i.message for i in iss))
check("baseline: SIR!=S note on element 1",
      any(i.element == 1 and i.field == "SIR" and i.level == "note" for i in iss))
check("baseline: NPRINT 7 warning on element 1",
      any(i.element == 1 and i.field == "NPRINT" and i.level == "warning" for i in iss))
check("baseline: NPRINT 0 note on element 2",
      any(i.element == 2 and i.field == "NPRINT" and i.level == "note" for i in iss))

def bad(mutate, field, level="error"):
    r = copy.deepcopy(run); mutate(r)
    return any(i.field.startswith(field) and i.level == level for i in E.validate(r))

check("note (not error) for rain ending before TFIN", bad(lambda r: setattr(r.gages[0], "times", [0, 60]), "gage 1", "note"))
check("reject >99 pairs", bad(lambda r: (setattr(r.gages[0], "times", list(range(100))), setattr(r.gages[0], "depths", list(range(100)))), "gage 1"))
check("reject channel (W=0)", bad(lambda r: setattr(r.element(1), "W", 0), "W"))
check("reject channel inflow NC1", bad(lambda r: setattr(r.element(2), "NC1", 1), "NR/NL"))
check("reject missing upstream", bad(lambda r: setattr(r.element(2), "NU", 9), "NU"))
check("reject wrong order", bad(lambda r: setattr(r, "order", [2, 1]), "NU"))
check("reject THI > THMX", bad(lambda r: setattr(r.element(1), "THI", 0.35), "THI"))
check("reject too many steps", bad(lambda r: setattr(r, "DELT", 0.01), "DELT"))
check("reject NTIME 3", bad(lambda r: setattr(r, "NTIME", 3), "NTIME"))
check("reject SIR=0", bad(lambda r: setattr(r.element(1), "SIR", 0), "SIR"))
check("warn steep slope", bad(lambda r: setattr(r.element(1), "SIR", 0.7), "S", "warning") or
      bad(lambda r: (setattr(r.element(1), "DEPNO", 5), setattr(r.element(1), "RILLW", .1), setattr(r.element(1), "S", 0.7)), "S", "warning"))
check("reject decreasing cumulative rain", bad(lambda r: setattr(r.gages[0], "depths", [10, 5]), "gage 1"))

# ---- 5. parse outputs, compare with the numbers from the baseline run ----
res = E.load_results(rd("test26.dyn"), rd("test26.sta"))
g = res["global"]; a, b = res["elements"]["1"], res["elements"]["2"]
check("title", res["title"].strip() == "test 5", repr(res["title"]))
check("hydrograph rows: 177 per element", len(a["hydrograph"]["time_min"]) == 177 and len(b["hydrograph"]["time_min"]) == 177,
      f'{len(a["hydrograph"]["time_min"])} {len(b["hydrograph"]["time_min"])}')
check("time axis 0..88 step 0.5", a["hydrograph"]["time_min"][0] == 0 and a["hydrograph"]["time_min"][-1] == 88)
check("el1 sediment total 16530.156 kg", close(a["sediment_total_kg"], 16530.156, 1e-3))
check("el2 sediment total 7878.426 kg", close(b["sediment_total_kg"], 7878.426, 1e-3))
check("el1 vol error 0.0712 %", close(a["vol_error_pct"], 0.0712, 1e-4))
check("el2 vol error -0.0236 %", close(b["vol_error_pct"], -0.0236, 1e-4))
check("areas 1 ha / 1.1 ha", close(a["contributing_area_ha"], 1.0, 1e-6) and close(b["contributing_area_ha"], 1.1, 1e-6))
check("peak 77.555 mm/h at 30 min (el 1)", a["peak_flow_mmh"] == 77.555 and a["time_to_peak_min"] == 30.0)
check("peak 77.287 mm/h (el 2)", b["peak_flow_mmh"] == 77.287)
check("sta: el1 peak sediment 773.39 kg/min", a["hydrology"]["peak_sed_kgmin"] == 773.39)
check("sta: el2 peak sediment 409.42 kg/min", b["hydrology"]["peak_sed_kgmin"] == 409.42)
check("sta: runoff durations 53.5 / 56.5", a["hydrology"]["runoff_duration_min"] == 53.5 and b["hydrology"]["runoff_duration_min"] == 56.5)
check("sta: time to runoff 2.5", a["hydrology"]["time_to_runoff_min"] == 2.5)
check("sta: net rainfall 49.865, peak rain 99.73", a["hydrology"]["net_rainfall_mm"] == 49.865 and a["hydrology"]["peak_rain_mmh"] == 99.73)
check("sta: net erosion el2 = -8651.723 (deposition)", b["erosion"]["net_kg"] == -8651.723 and b["erosion"]["net_tha"] == -86.517)
check("sta: gross interrill el1 16530.055", a["erosion"]["gross_interrill_kg"] == 16530.055)
check("sta: echoed inputs", a["input"]["XL"] == 100 and a["input"]["SIR"] == 0.095 and a["input"]["D50"] == 25 and b["input"]["NU"] == 1,
      str(a["input"]))
check("global: runoff 295.1976 m3 / 26.83614 mm", g["runoff_m3"] == 295.1976 and g["runoff_mm"] == 26.83614)
check("global: infiltration 22.99585, error 0.0648", g["infiltration_planes_mm"] == 22.99585 and g["vol_error_pct"] == 0.0648)
check("global: basin area 11000", g["basin_area_m2"] == 11000)

# integrals of the printed hydrograph agree with the model's own totals
check("el1 hydrograph volume ~271.1 m3", close(a["outflow_m3"], 271.1, 1.0), f'{a["outflow_m3"]:.2f}')
check("el2 hydrograph volume ~295.2 m3 (global runoff)", close(b["outflow_m3"], g["runoff_m3"], 1.5), f'{b["outflow_m3"]:.2f}')
check("el1 hydrograph sediment ~16530 kg", close(a["outflow_sediment_kg"], 16530.156, 80), f'{a["outflow_sediment_kg"]:.1f}')
check("el2 hydrograph sediment ~7878 kg", close(b["outflow_sediment_kg"], 7878.426, 40), f'{b["outflow_sediment_kg"]:.1f}')
check("peak of listed q equals reported peak", close(max(a["hydrograph"]["q_mmh"]), 77.555, 1e-3))

print("\n%d failure(s)" % len(fails))
sys.exit(1 if fails else 0)
