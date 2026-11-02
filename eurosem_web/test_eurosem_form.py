"""Tests for eurosem_form.  Run: python3 test_eurosem_form.py"""
import os, sys, copy
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import eurosem_io as E
from eurosem_form import build_run, DEFAULT_SPEC, FormError, PLANE_DEFAULTS
REF = os.path.join(HERE, "reference")
fails = []
def check(name, cond, detail=""):
    print(("PASS  " if cond else "FAIL  ") + name + ("" if cond else "   " + str(detail)[:300]))
    if not cond: fails.append(name)

base = E.read_pcp(open(f"{REF}/test.pcp").read(), E.read_par(open(f"{REF}/testn.par").read()))
run, acc = build_run(DEFAULT_SPEC)
check("default form validates with no errors", not E.has_errors(E.validate(run)), [str(i) for i in E.validate(run)])
check("order is 1,2", run.order == [1, 2])
# every parameter equals the reference file, except what the form sets deliberately
b = copy.deepcopy(base)
for e in b.elements:
    e.NPRINT = 1; e.S = e.SIR
ok = all(a == c for a, c in zip(run.elements, b.elements))
check("default planes equal the reference planes", ok,
      [(x.J, {k: (getattr(x, k), getattr(y, k)) for k in vars(x) if getattr(x, k) != getattr(y, k)}) for x, y in zip(run.elements, b.elements)])
check("system record equals reference", (run.CLEN, run.TFIN, run.DELT, run.THETA, run.TEMP, run.NTIME, run.NEROS) ==
      (base.CLEN, base.TFIN, base.DELT, base.THETA, base.TEMP, base.NTIME, base.NEROS))
check("rainfall equals reference", run.gages[0].times == base.gages[0].times and run.gages[0].depths == base.gages[0].depths)
check("written files re-read identically", E.read_par(E.write_par(run)).elements == run.elements)

def bad(mut, text):
    s = copy.deepcopy(DEFAULT_SPEC); mut(s)
    try: build_run(s)
    except FormError as e: return text.lower() in str(e).lower()
    return False
check("no planes", bad(lambda s: s.update(planes=[]), "at least one plane"))
check("upstream loop", bad(lambda s: (s["planes"][0].update(upstream=2)), "loop"))
check("plane feeding two planes", bad(lambda s: s["planes"].append(dict(s["planes"][0], upstream=1)), "more than one"))
check("two outlets", bad(lambda s: s["planes"][1].update(upstream=0), "single outlet"))
check("non-numeric value", bad(lambda s: s["planes"][0].update(length="abc"), "not a number"))
check("too many planes", bad(lambda s: s.update(planes=[dict(DEFAULT_SPEC["planes"][0]) for _ in range(9)]), "at most"))
check("rain needs rows", bad(lambda s: s.update(rain=[[0, 0]]), "two rows"))
# downstream listed first still orders upstream-first
s = copy.deepcopy(DEFAULT_SPEC); s["planes"][0]["upstream"] = 2; s["planes"][1]["upstream"] = 0
r, _ = build_run(s); check("order follows upstream links", r.order == [2, 1], r.order)
# plant shape, stem angle and stone indicator now come from the form (defaults are the reference values)
r0, _ = build_run(DEFAULT_SPEC); e0 = r0.element(1)
check("shape/angle/stone defaults match the reference file", (e0.SHAPE, e0.PLANGLE, e0.STONE) == (2.0, 55.0, 1.0), (e0.SHAPE, e0.PLANGLE, e0.STONE))
s2 = copy.deepcopy(DEFAULT_SPEC); s2["planes"][0].update(plant_shape=1, plant_angle=85, stone=-1)
e2 = build_run(s2)[0].element(1); check("shape/angle/stone are taken from the form", (e2.SHAPE, e2.PLANGLE, e2.STONE) == (1.0, 85.0, -1.0))

# --- rills ---------------------------------------------------------------
def _rp(**kw):
    sp = {"title": "r", "duration": 60, "timestep": 0.5, "temperature": 20, "rain": [[0, 0], [30, 30]],
          "planes": [dict(PLANE_DEFAULTS, rills=1, slope=0.08, **kw)]}
    return sp
def _try(sp):
    try:
        return build_run(sp)[0], None
    except FormError as ex:
        return None, " | ".join(ex.problems)
r, err = _try(_rp())
e = r.elements[0] if r else None
check("rills on: geometry reaches the element", e and (e.DEPNO, e.RILLW, e.RILLD, e.ZLR, e.RS) == (20.0, 0.3, 0.15, 1.0, 1.0) and e.has_rills, err)
check("rills on: S is the slope, SIR defaults to 1.4 x S", e and abs(e.S - 0.08) < 1e-12 and abs(e.SIR - 0.112) < 1e-12, e and (e.S, e.SIR))
r, err = _try(_rp(interrill_slope=0.2, rill_manning=0.05, rill_uniform=0))
e = r.elements[0]
check("explicit interrill slope, rill Manning n and growing depth", (e.SIR, e.MANN, e.MANN_IR, e.RS) == (0.2, 0.05, 0.188, 0.0), (e.SIR, e.MANN, e.MANN_IR, e.RS))
check("interrill slope not steeper than S refused", _try(_rp(interrill_slope=0.08))[1] is not None)
check("rills that do not fit the width refused", "does not fit" in (_try(_rp(rill_number=80, rill_width=1.5))[1] or ""))
check("zero depth refused", _try(_rp(rill_depth=0))[1] is not None)
r, err = _try(dict(_rp(), planes=[dict(PLANE_DEFAULTS, rills=0, rill_number=99, rill_depth=0)]))
e = r.elements[0]
check("rills off: fields ignored, no rills", not e.has_rills and e.DEPNO == 0 and e.SIR == e.S, err)

# --- RFR ---------------------------------------------------------------
_s = copy.deepcopy(DEFAULT_SPEC); _s["planes"][0]["rfr"] = 25
check("RFR default 0 and taken from the form", build_run(DEFAULT_SPEC)[0].elements[0].RFR == 0 and build_run(_s)[0].elements[0].RFR == 25)
_s["planes"][0]["rfr"] = 90
check("RFR out of range refused", _try(_s)[1] is not None)
print("\n%d failure(s)" % len(fails)); sys.exit(1 if fails else 0)
