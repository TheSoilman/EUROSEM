"""
eurosem_sens.py -- one-at-a-time sensitivity to parameter uncertainty.

For each chosen parameter the model is run twice more, with the parameter at
a LOW and a HIGH value while everything else stays as in the form (the base
case).  The parameter is changed on every plane that has a range for it.

Ranges
  * guide ranges come from the plane's chosen soil texture (eurosem_guide.py):
        conductivity  Table A4.2      capillary drive, porosity  Table A4.1
        erodibility   Table A9.1      cohesion  Table A9.2 (uncompacted or
                                      compacted, following the surface state)
    Rows the guide prints inconsistently (A9.2 anomalies) get no range.
  * the user can type a low and high for any parameter; that replaces the
    guide range on all planes.
  * porosity / moisture: if a low porosity would fall below a plane's maximum
    moisture (or initial moisture), those are lowered to match so the run is
    valid, and the row says so.

This shows how much the output moves for plausible changes in one input at a
time.  It is not a confidence interval and it ignores interactions.
"""
import copy
import datetime
import time
from typing import Callable, Dict, List, Optional, Tuple

import eurosem_guide as G
import eurosem_export as X

# (key, label, unit, where a guide range comes from or None)
PARAMS = [
    ("manning", "Manning n", "", None),
    ("slope", "Slope", "m/m", None),
    ("fmin", "Conductivity (Ks)", "mm/h", "Table A4.2"),
    ("g", "Capillary drive", "mm", "Table A4.1"),
    ("porosity", "Porosity", "", "Table A4.1"),
    ("theta_i", "Initial moisture", "", None),
    ("cover", "Canopy cover", "fraction", None),
    ("plant_height", "Plant height", "cm", None),
    ("interception", "Interception", "mm", None),
    ("d50", "Median particle size", "um", None),
    ("erodibility", "Erodibility", "g/J", "Table A9.1"),
    ("cohesion", "Cohesion", "kPa", "Table A9.2"),
]
LABEL = {k: (l, u) for k, l, u, _ in PARAMS}
MAX_RUNS = 40


def guide_range(plane: dict, key: str) -> Optional[Tuple[float, float, str]]:
    """(low, high, source) from the plane's chosen texture, or None."""
    tex = plane.get("texture")
    if tex not in G.A41:
        return None
    if key == "fmin":
        lo, _, hi = G.A42[tex]; return lo, hi, "Table A4.2"
    if key == "g":
        p = G.A41[tex]; return p[5], p[7], "Table A4.1"
    if key == "porosity":
        p = G.A41[tex]; return p[0], p[2], "Table A4.1"
    if key == "erodibility" and tex in G.A91:
        lo, _, hi = G.A91[tex]; return lo, hi, "Table A9.1"
    if key == "cohesion" and tex in G.A92 and tex not in G.ANOMALIES["A9.2"]:
        c = G.A92[tex]
        if plane.get("condition") == "compacted" or not plane.get("condition"):
            # default surface state in the form is "sealed or compacted"
            if c[3] is not None:
                return c[3], c[5], "Table A9.2 (compacted)"
            return None
        return c[0], c[2], "Table A9.2 (uncompacted)"
    return None


def plan(spec: dict) -> List[dict]:
    """What each parameter's range would be, per plane, before any run."""
    out = []
    for key, label, unit, src in PARAMS:
        per = []
        for i, pl in enumerate(spec["planes"], start=1):
            r = guide_range(pl, key)
            per.append(None if r is None else {"plane": i, "low": r[0], "high": r[1], "source": r[2]})
        out.append({"key": key, "label": label, "unit": unit, "guide": per,
                    "has_guide": any(per), "base": [pl.get(key) for pl in spec["planes"]]})
    return out


def _f(x):
    return float(x)


def _outlet(spec: dict) -> int:
    feeds = {int(float(p.get("upstream", 0) or 0)) for p in spec["planes"]} - {0}
    for i in range(1, len(spec["planes"]) + 1):
        if i not in feeds:
            return i
    return len(spec["planes"])


def metrics(results: dict, spec: dict) -> dict:
    o = results["elements"][str(_outlet(spec))]
    g = results["global"]
    return {"runoff_mm": g["runoff_mm"], "infiltration_mm": g["infiltration_planes_mm"],
            "peak_flow_mmh": o["peak_flow_mmh"], "time_to_peak_min": o["time_to_peak_min"],
            "sediment_kg": o["sediment_total_kg"]}


def _apply(spec: dict, key: str, values: List[Optional[float]]) -> Tuple[dict, List[str]]:
    s = copy.deepcopy(spec)
    notes = []
    for pl, v in zip(s["planes"], values):
        if v is None:
            continue
        pl[key] = v
        if key in ("porosity", "theta_i"):
            por, tm, ti = _f(pl["porosity"]), _f(pl["theta_max"]), _f(pl["theta_i"])
            if tm > por:
                pl["theta_max"] = por; tm = por
                notes.append("maximum moisture lowered to equal porosity")
            if ti > tm:
                pl["theta_i"] = tm
                notes.append("initial moisture lowered to stay within maximum moisture")
    return s, sorted(set(notes))


def run_sensitivity(spec: dict, keys: List[str], overrides: Dict[str, dict],
                    run_one: Callable[[dict], Tuple[bool, object]],
                    max_seconds: float = 150.0) -> dict:
    """run_one(spec) -> (True, results) or (False, message)."""
    t0 = time.time()
    overrides = overrides or {}
    n = len(spec["planes"])
    ok, base = run_one(spec)
    if not ok:
        return {"ok": False, "error": "the base case did not run: %s" % base}
    out = {"ok": True, "outlet": _outlet(spec), "base": metrics(base, spec),
           "rows": [], "skipped": [], "runs": 1}
    wanted = [p for p in PARAMS if p[0] in set(keys)]
    for key, label, unit, _ in wanted:
        ov = overrides.get(key) or {}
        if ov.get("low") not in (None, "") and ov.get("high") not in (None, ""):
            lo, hi = _f(ov["low"]), _f(ov["high"])
            if lo > hi:
                lo, hi = hi, lo
            lows, highs = [lo] * n, [hi] * n
            source = "entered by you"
        else:
            rs = [guide_range(pl, key) for pl in spec["planes"]]
            if not any(rs):
                out["skipped"].append({"key": key, "label": label,
                    "reason": "no range: choose a soil texture with a guide value, or enter a low and high"})
                continue
            lows = [None if r is None else r[0] for r in rs]
            highs = [None if r is None else r[1] for r in rs]
            source = "; ".join(sorted({r[2] for r in rs if r}))
        if out["runs"] + 2 > MAX_RUNS or time.time() - t0 > max_seconds:
            out["skipped"].append({"key": key, "label": label, "reason": "run limit reached"})
            continue
        row = {"key": key, "label": label, "unit": unit, "source": source}
        for case, vals in (("low", lows), ("high", highs)):
            s, notes = _apply(spec, key, vals)
            ok, res = run_one(s)
            out["runs"] += 1
            row[case] = {"values": vals, "notes": notes,
                         "metrics": metrics(res, s) if ok else None,
                         "error": None if ok else str(res)}
        out["rows"].append(row)
    out["seconds"] = round(time.time() - t0, 1)
    return out


def sensitivity_csv(spec: dict, out: dict) -> str:
    """Rows per parameter and case, with the base settings recorded above."""
    head = [["EUROSEM sensitivity: one parameter at a time, guide or user-entered low and high values"],
            ["Model", X.MODEL],
            ["Created (UTC)", datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M")],
            ["Reading it", "base = the form's values. Each row changes ONE parameter (on every plane that has a range) "
                           "and shows the outlet results and the change from base. Not a confidence interval; interactions are ignored."],
            ["Base case", "", ] + [f"{k}={round(v, 6)}" for k, v in out["base"].items()],
            []]
    head += X._rows_settings(spec)
    head.append(["End of header: the table below starts at the row beginning parameter"])
    cols = ["parameter", "unit", "case", "values_by_plane", "range_source", "adjustments",
            "runoff_mm", "peak_flow_mmh", "time_to_peak_min", "sediment_kg",
            "runoff_change_pct", "peak_flow_change_pct", "sediment_change_pct", "error"]
    body = [cols]
    b = out["base"]
    pct = lambda v, k: "" if v is None or not b[k] else round((v - b[k]) / b[k] * 100, 3)
    for r in out["rows"]:
        for case in ("low", "high"):
            c = r[case]; m = c["metrics"]
            vals = "; ".join(f"plane{i}={v}" for i, v in enumerate(c["values"], 1) if v is not None)
            body.append([r["label"], r["unit"], case, vals, r["source"], "; ".join(c["notes"]),
                         "" if not m else round(m["runoff_mm"], 5), "" if not m else round(m["peak_flow_mmh"], 4),
                         "" if not m else m["time_to_peak_min"], "" if not m else round(m["sediment_kg"], 3),
                         pct(m and m["runoff_mm"], "runoff_mm"), pct(m and m["peak_flow_mmh"], "peak_flow_mmh"),
                         pct(m and m["sediment_kg"], "sediment_kg"), c["error"] or ""])
    for s in out["skipped"]:
        body.append([s["label"], "", "skipped", "", "", "", "", "", "", "", "", "", "", s["reason"]])
    return X._write(head, prefix="# ") + X._write(body)
