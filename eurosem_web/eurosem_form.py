"""
eurosem_form.py  --  turns the web form's values into a model run.

The form deals in plain quantities (a plane is a length, a width, a slope,
soil and cover values; rain is a table of time and cumulative depth).
This module builds the full Run that eurosem_io writes to the model's
files, and refuses combinations the model cannot handle.

Planes only, in a single cascade:
  * planes are numbered 1, 2, 3 ... in the order the form lists them
  * each plane may receive runoff from ONE other plane ("upstream")
  * a plane may feed at most one other plane; exactly one plane (the
    outlet) feeds none
  * the computation order (upstream first) is worked out here
  * one slope per plane: it is written to both S and SIR, because for a
    plane without rills the model uses SIR and ignores S.  With the rills
    switch on, the slope is the slope along the rills and an interrill
    slope (steeper) is added.
"""
from typing import List, Tuple

import eurosem_io as E

MAX_PLANES = 8          # the form's limit (the model allows 60)
MAX_RAIN_ROWS = 99      # the model's limit per gage

# --- default form: the two-plane reference case ------------------------
PLANE_DEFAULTS = dict(
    length=100.0, width=100.0, slope=0.095, manning=0.188,
    fmin=5.0, g=400.0, porosity=0.4, theta_i=0.15, theta_max=0.3, roc=0.0,
    interception=0.9, storage=1.0, rfr=0.0,
    cover=0.15, plant_height=40.0, plant_base=0.05,
    d50=25.0, erodibility=0.9, cohesion=1.0, pave=0.0,
    plant_shape=2.0, plant_angle=55.0, stone=1.0,
    splash_texture=2.0, particle_density=2.6, upstream=0,
    # rills (off by default; the model sees no rills when rills is 0)
    rills=0, rill_number=20.0, rill_width=0.3, rill_depth=0.15, rill_side=1.0,
    rill_uniform=1, interrill_slope=0.0, rill_manning=0.0)

DEFAULT_SPEC = {
    "title": "test",
    "duration": 88.0,
    "timestep": 0.5,
    "temperature": 20.0,
    "rain": [[0.0, 0.0], [30.0, 50.0]],
    "planes": [
        dict(PLANE_DEFAULTS),
        dict(PLANE_DEFAULTS, length=10.0, slope=0.05, upstream=1),
    ],
}


class FormError(ValueError):
    """The form values cannot be turned into a run.  .problems is a list
    of plain-language messages."""
    def __init__(self, problems: List[str]):
        super().__init__("; ".join(problems))
        self.problems = problems


def _num(d: dict, key: str, label: str, problems: List[str], default=None) -> float:
    v = d.get(key, default)
    try:
        f = float(v)
        if f != f or f in (float("inf"), float("-inf")):
            raise ValueError
        return f
    except (TypeError, ValueError):
        problems.append(f"{label} is not a number")
        return 0.0


def _rills(i, pl, e, g, p):
    """Rill geometry for plane i.  The form's slope is the slope ALONG the
    rills (S); the interrill slope SIR must be steeper (blank or 0 = the
    model's own default, 1.4 x S).  RILLD must be given (> 0)."""
    n, w = g("rill_number", "number of rills"), g("rill_width", "rill width")
    d, z = g("rill_depth", "rill depth"), g("rill_side", "rill side slope")
    sir, mr = g("interrill_slope", "interrill slope"), g("rill_manning", "rill Manning n")
    uni = int(g("rill_uniform", "rill depth pattern") or 0)
    if n < 1 or w <= 0 or d <= 0 or z < 0:
        p.append(f"Plane {i}: rills need a number of 1 or more, a width and depth above 0, and a side slope of 0 or more")
        return
    if e.S <= 0:
        p.append(f"Plane {i}: with rills the slope along the rills must be above 0")
        return
    if sir <= 0:
        sir = 1.4 * e.S
    elif sir <= e.S:
        p.append(f"Plane {i}: the interrill slope ({sir:g}) must be steeper than the slope along the rills ({e.S:g})")
        return
    top = w + 2.0 * z * d
    if n * top >= e.W:
        p.append(f"Plane {i}: {n:g} rills with a top width of {top:g} m take up {n * top:g} m, "
                 f"which does not fit in the plane width of {e.W:g} m")
        return
    e.DEPNO, e.RILLW, e.RILLD, e.ZLR = n, w, d, z
    e.RS = 1.0 if uni else 0.0
    e.SIR = sir
    if mr > 0:
        e.MANN = mr


def build_run(spec: dict) -> Tuple[E.Run, str]:
    """Return (run, accept).  `accept` is the Y/N answer for the model's
    ROC / PAVE questions: Y, because a value entered in the form is meant
    to be used.  Raises FormError."""
    p: List[str] = []
    planes = spec.get("planes")
    if not isinstance(planes, list) or not planes:
        raise FormError(["add at least one plane"])
    if len(planes) > MAX_PLANES:
        raise FormError([f"at most {MAX_PLANES} planes are allowed in the form"])
    rain = spec.get("rain")
    if not isinstance(rain, list) or len(rain) < 2:
        raise FormError(["the rainfall table needs at least two rows"])
    if len(rain) > MAX_RAIN_ROWS:
        raise FormError([f"the rainfall table can have at most {MAX_RAIN_ROWS} rows"])

    run = E.Run()
    run.title = " ".join(str(spec.get("title", "EUROSEM run")).split())[:80] or "EUROSEM run"
    run.TFIN = _num(spec, "duration", "Duration", p, 60)
    run.DELT = _num(spec, "timestep", "Time step", p, 0.5)
    run.TEMP = _num(spec, "temperature", "Temperature", p, 20)
    run.NTIME, run.NEROS, run.NPART, run.THETA = 2, 2, 0, 0.7

    n = len(planes)
    ups = []
    for i, pl in enumerate(planes, start=1):
        u = int(_num(pl, "upstream", f"Plane {i} upstream", p, 0))
        if u < 0 or u > n or u == i:
            p.append(f"Plane {i}: upstream must be none or another plane")
            u = 0
        ups.append(u)
    feeds = [u for u in ups if u]
    for u in sorted(set(feeds)):
        if feeds.count(u) > 1:
            p.append(f"Plane {u} feeds more than one plane; a plane can feed only one other plane")
    outlets = [i for i in range(1, n + 1) if i not in feeds]
    if len(outlets) != 1 and len(feeds) < n:
        p.append("the planes must drain to a single outlet plane "
                 f"(found {len(outlets)} planes that feed no other plane)")

    # computation order: upstream before downstream
    order, placed = [], set()
    for _ in range(n):
        for i in range(1, n + 1):
            if i not in placed and (ups[i - 1] == 0 or ups[i - 1] in placed):
                order.append(i); placed.add(i)
                break
        else:
            break
    if len(order) != n:
        p.append("the planes feed each other in a loop")

    for i, pl in enumerate(planes, start=1):
        g = lambda k, lab: _num(pl, k, f"Plane {i} {lab}", p, PLANE_DEFAULTS[k])
        e = E.Element(J=i, NU=ups[i - 1], NPRINT=1)
        e.XL, e.W = g("length", "length"), g("width", "width")
        e.S = e.SIR = g("slope", "slope")
        e.MANN = e.MANN_IR = g("manning", "Manning n")
        if int(g("rills", "rills switch") or 0):
            _rills(i, pl, e, g, p)
        e.FMIN, e.G = g("fmin", "saturated conductivity"), g("g", "capillary drive")
        e.POR, e.THI, e.THMX = g("porosity", "porosity"), g("theta_i", "initial moisture"), g("theta_max", "maximum moisture")
        e.RFR = g("rfr", "roughness ratio RFR")
        if e.RFR < 0 or e.RFR > 60:
            p.append(f"Plane {i}: roughness ratio RFR must be between 0 and 60 (guide values are 13 to 35)")
        e.ROC, e.DINT, e.RECS = g("roc", "rock fraction"), g("interception", "interception"), g("storage", "surface storage")
        e.COVER, e.PLANTH, e.PLANTBASE = g("cover", "cover"), g("plant_height", "plant height"), g("plant_base", "plant base")
        e.D50, e.EROD, e.COH, e.PAVE = g("d50", "D50"), g("erodibility", "erodibility"), g("cohesion", "cohesion"), g("pave", "pavement")
        e.SPLTEX, e.RHOS = g("splash_texture", "splash texture"), g("particle_density", "particle density")
        e.SHAPE, e.PLANGLE, e.STONE = g("plant_shape", "plant shape"), g("plant_angle", "plant stem angle"), g("stone", "stone indicator")
        run.elements.append(e)
    run.order = order

    times, depths = [], []
    for k, row in enumerate(rain, start=1):
        try:
            times.append(float(row[0])); depths.append(float(row[1]))
        except (TypeError, ValueError, IndexError):
            p.append(f"Rainfall row {k} needs a time and a depth")
    run.gages = [E.Gage(ident="form rainfall", times=times, depths=depths)]
    run.weights = {e.J: (1, 1.0) for e in run.elements}

    if p:
        raise FormError(p)
    # CLEN (characteristic length) is reported by the model; the reference
    # case uses 110 = 100 + 10, the sum of the plane lengths.
    run.CLEN = sum(e.XL for e in run.elements)
    return run, "Y"
