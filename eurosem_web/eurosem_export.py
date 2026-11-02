"""
eurosem_export.py  --  CSV files for plotting and record keeping.

make_results_csv(spec, results)   one wide table, one row per time step:
      time_min, rain_mmh, then for every plane
      plane1_q_mmh, plane1_q_m3min, plane1_conc, plane1_sed_kgmin,
      plane1_cum_sed_kg, plane2_...
    Above the table are comment rows (first cell starts with "# ") that
    record the run: title, date, model version, event settings, rainfall
    table, every plane's parameter values and the summary results.
    Reading it:  pandas.read_csv(path, comment="#")  or in Excel the
    comment rows simply appear above the table.  (No byte-order mark, so
    that every comment row begins exactly with "#".)
make_settings_csv(spec)           the settings alone, as an ordinary table.
"""
import csv
import datetime
import io
from typing import List

import eurosem_guide as G

MODEL = "EUROSEM v3.2L (headless build, LOG10 kinetic-energy correction)"

# (key, label, unit) in the order they appear in the form
PLANE_PARAMS = [
    ("upstream", "Receives runoff from plane (0 = none)", ""),
    ("length", "Length", "m"), ("width", "Width", "m"), ("slope", "Slope", "m/m"),
    ("manning", "Manning n", ""),
    ("fmin", "Saturated conductivity (Ks)", "mm/h"), ("g", "Capillary drive", "mm"),
    ("porosity", "Porosity", ""), ("theta_i", "Initial moisture", ""),
    ("theta_max", "Maximum moisture", ""), ("roc", "Rock fraction", ""),
    ("cover", "Canopy cover", "fraction"), ("plant_height", "Plant height", "cm"),
    ("interception", "Interception", "mm"), ("plant_base", "Plant base fraction", ""),
    ("storage", "Microrelief height RECS (infiltration recession)", "mm"),
    ("rfr", "Roughness ratio RFR (sets depression storage)", "cm/m"),
    ("plant_shape", "Plant shape (0 none, 1 bladed or needle, 2 broad)", ""),
    ("plant_angle", "Plant stem angle", "deg"), ("stone", "Stone indicator (+1 or -1)", ""),
    ("d50", "Median particle size", "um"), ("erodibility", "Erodibility", "g/J"),
    ("cohesion", "Cohesion", "kPa"), ("pave", "Pavement fraction", ""),
    ("rills", "Rills (1 = on)", ""), ("rill_number", "Number of rills across width", ""),
    ("rill_width", "Rill bottom width", "m"), ("rill_depth", "Rill depth", "m"),
    ("rill_side", "Rill side slope (1 : x)", ""), ("rill_uniform", "Rill depth uniform along plane (1) or grows downslope (0)", ""),
    ("interrill_slope", "Interrill slope (0 = 1.4 x slope)", "m/m"), ("rill_manning", "Rill Manning n (0 = same as plane)", ""),
    ("splash_texture", "Splash texture", ""), ("particle_density", "Particle density", "g/cm3"),
]
SERIES = [("q_mmh", "q_mmh"), ("q_m3min", "q_m3min"), ("conc", "conc"),
          ("qs_kgmin", "sed_kgmin"), ("cum_sediment_kg", "cum_sed_kg")]


def _num(v):
    try:
        f = float(v)
        return int(f) if f == int(f) and abs(f) < 1e12 else f
    except (TypeError, ValueError):
        return v


def _rows_settings(spec: dict) -> List[list]:
    planes = spec["planes"]
    rows = [["Run title", "", spec.get("title", "")],
            ["Duration", "min", _num(spec.get("duration"))],
            ["Time step", "min", _num(spec.get("timestep"))],
            ["Water temperature", "C", _num(spec.get("temperature"))],
            ["Number of planes", "", len(planes)]]
    rows.append([])
    rows.append(["Rainfall: time (min)", "cumulative depth (mm)"])
    for t, d in spec["rain"]:
        rows.append([_num(t), _num(d)])
    rows.append([])
    rows.append(["Plane parameter", "unit"] + [f"plane {i}" for i in range(1, len(planes) + 1)])
    for key, label, unit in PLANE_PARAMS:
        rows.append([label, unit] + [_num(p.get(key, "")) for p in planes])
    rows.append([])
    rows.append(["Guide choices", ""] + [f"plane {i}" for i in range(1, len(planes) + 1)])
    info = [_guide_info(p) for p in planes]
    rows.append(["Soil texture chosen (guide)", ""] + [x[0] for x in info])
    rows.append(["Surface state chosen", ""] + [x[1] for x in info])
    rows.append(["Values compared with the guide fill", ""] + [x[2] for x in info])
    picks = [("Roughness cover chosen (guide A3.1)", "manning_ref", {G._m_label(l, q) for l, q, *_ in G.A31}),
             ("Vegetation type chosen (guide A7.2)", "veg_ref", set(G.A72)),
             ("Interception type chosen (guide A7.1)", "dint_ref", set(G.A71)),
             ("Basal area class chosen (guide A7.3)", "pbase_ref", {G._pb_label(l, c) for l, c, _ in G.A73}),
             ("Root type chosen (guide A9.3)", "root_ref", set(G.A93)),
             ("Tillage chosen for RFR (guide A6.1)", "till_ref", set(G.A61))]
    for label, key, allowed in picks:
        vals = [p.get(key) if p.get(key) in allowed else "" for p in planes]
        if any(vals):                       # only names the guide module knows are ever written
            rows.append([label, ""] + [v or "none" for v in vals])
    if any(x[0] != "none" for x in info):
        rows.append(["Guide source", "", "EUROSEM user guide (Morgan, Quinton et al. 1998): Tables A4.1, A4.2 "
                     "(soil water), A9.1 (erodibility), A9.2 (cohesion). Mean values are used for soil water."])
    return rows


def _guide_info(p: dict):
    """(texture, surface state, comparison) for one plane, as text.  Only a
    texture the guide module knows is written (anything else is ignored)."""
    tex = p.get("texture")
    if tex not in G.A41:
        return ("none", "", "values entered by hand")
    cond = p.get("condition") if p.get("condition") in G.CONDITIONS else "compacted"
    fill = G.guide_fill(tex, cond)
    changed = []
    for k, v in fill.items():
        try:
            if abs(float(p.get(k)) - float(v)) > 1e-9:
                changed.append(G.FILL_LABEL[k])
        except (TypeError, ValueError):
            changed.append(G.FILL_LABEL[k])
    absent = [G.FILL_LABEL[k] for k in ("erodibility", "cohesion") if k not in fill]
    text = "as filled" if not changed else "edited after filling: " + "; ".join(changed)
    if absent:
        text += " (no guide value for: " + "; ".join(absent) + ")"
    return (tex, G.CONDITIONS[cond][0], text)


def _write(rows, prefix="") -> str:
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    for r in rows:
        r = list(r)
        if prefix:
            # the first cell must stay unquoted so the row begins with "#"
            first = str(r[0]) if r else ""
            first = first.replace(",", ";").replace('"', "'").replace("\n", " ")
            r = [prefix + first] + r[1:]
        w.writerow(r)
    return buf.getvalue()


def make_settings_csv(spec: dict) -> str:
    return "\ufeff" + _write(_rows_settings(spec))


def make_results_csv(spec: dict, results: dict) -> str:
    ids = sorted(results["elements"], key=int)
    g = results["global"]
    head = [["EUROSEM results: wide format, one row per time step"],
            ["Model", MODEL],
            ["Created (UTC)", datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M")],
            ["Columns", "plane<N>_q_mmh = discharge (mm/h); _q_m3min = discharge (m3/min); "
                        "_conc = sediment concentration as printed by the model (CONC.); "
                        "_sed_kgmin = sediment discharge (kg/min); _cum_sed_kg = cumulative sediment out (kg)"],
            ["Discharge and sediment are the flow leaving each plane, including inflow from any plane above it"],
            []]
    head += _rows_settings(spec)
    head.append([])
    head.append(["Summary", "unit"] + [f"plane {i}" for i in ids])
    def per(key, sub=None):
        out = []
        for i in ids:
            e = results["elements"][i]
            v = e.get(sub, {}).get(key) if sub else e.get(key)
            out.append("" if v is None else round(v, 6))
        return out
    for lab, unit, key, sub in [
        ("Contributing area", "ha", "contributing_area_ha", None),
        ("Time to runoff", "min", "time_to_runoff_min", "hydrology"),
        ("Duration of runoff", "min", "runoff_duration_min", "hydrology"),
        ("Peak flow", "mm/h", "peak_flow_mmh", "hydrology"),
        ("Time to peak flow", "min", "time_to_peak_flow_min", "hydrology"),
        ("Peak sediment discharge", "kg/min", "peak_sed_kgmin", "hydrology"),
        ("Water out", "m3", "outflow_m3", None),
        ("Sediment out", "kg", "sediment_total_kg", None),
        ("Gross rill erosion (rilled planes only)", "t/ha", "gross_rill_tha", "erosion"),
        ("Gross interrill erosion", "t/ha", "gross_interrill_tha", "erosion"),
        ("Net erosion (negative = deposition)", "t/ha", "net_tha", "erosion"),
        ("Volume error", "%", "vol_error_pct", None)]:
        head.append([lab, unit] + per(key, sub))
    head.append(["Basin runoff", "mm", round(g["runoff_mm"], 5)])
    head.append(["Basin runoff", "m3", round(g["runoff_m3"], 4)])
    head.append(["Global volume error", "%", g["vol_error_pct"]])
    head.append(["End of header: the table below starts at the row beginning time_min"])

    first = results["elements"][ids[0]]["hydrograph"]
    cols = ["time_min", "rain_mmh"]
    for i in ids:
        cols += [f"plane{i}_{name}" for _, name in SERIES]
    body = [cols]
    for k, t in enumerate(first["time_min"]):
        row = [t, first["rain_mmh"][k]]
        for i in ids:
            h = results["elements"][i]["hydrograph"]
            row += [h[key][k] for key, _ in SERIES]
        body.append(row)
    return _write(head, prefix="# ") + _write(body)
