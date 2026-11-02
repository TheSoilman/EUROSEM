"""
eurosem_guide.py -- guide values transcribed from the EUROSEM user guide
(Morgan, Quinton et al. 1998, Cranfield/Silsoe).  Every number below was read
from the page images, not from the garbled text layer.

Each table gives (low, mean, high) unless stated.  Where the PRINTED guide is
internally inconsistent the value is kept exactly as printed and the row is
listed in ANOMALIES so the interface can say so instead of hiding it.

Units follow the guide's tables; conversion to the model's input units is the
caller's job (see NOTES): the guide quotes POR/THR/THMAX as v/v fractions,
vegetation height in METRES (the model input PLANTH is in CM).
"""

TEXTURES = ["Sand", "Loamy sand", "Sandy loam", "Loam", "Silt loam",
            "Sandy clay loam", "Clay loam", "Silty clay loam",
            "Sandy clay", "Silty clay", "Clay"]

# Table A4.1 (guide p.10 of Appendix 4) -- Woolhiser et al. (1990) KINEROS values.
# Printed rows are shifted one line below "Sandy clay loam" (label wraps); the
# mapping below was checked against standard porosity values.
# texture: (POR low, mean, high, THR mean, THMAX mean, G low, mean, high mm)
A41 = {
    "Sand":            (0.37, 0.44, 0.50, 0.020, 0.42,  22, 101,  207),
    "Loamy sand":      (0.37, 0.44, 0.51, 0.035, 0.41,  41, 147,  323),
    "Sandy loam":      (0.35, 0.45, 0.56, 0.040, 0.41,  98, 248,  526),
    "Loam":            (0.38, 0.46, 0.55, 0.030, 0.43, 185, 375,  937),
    "Silt loam":       (0.42, 0.50, 0.58, 0.015, 0.47, 220, 485, 1043),
    "Sandy clay loam": (0.33, 0.40, 0.46, 0.070, 0.33, 220, 617, 1070),
    "Clay loam":       (0.41, 0.46, 0.52, 0.070, 0.39, 250, 533, 1174),
    "Silty clay loam": (0.42, 0.47, 0.52, 0.380, 0.43, 370, 720, 1470),
    "Sandy clay":      (0.37, 0.43, 0.49, 0.110, 0.32, 373, 768, 1730),
    "Silty clay":      (0.43, 0.48, 0.53, 0.060, 0.42, 430, 812, 1700),
    "Clay":            (0.43, 0.48, 0.53, 0.090, 0.39, 460, 890, 1830),
}

# Table A4.2: saturated hydraulic conductivity, mm/h (low, mean, high)
A42 = {
    "Sand": (170, 210, 600), "Loamy sand": (18, 61, 800),
    "Sandy loam": (7, 26, 190), "Loam": (2, 13, 65), "Silt loam": (3, 7, 25),
    "Sandy clay loam": (1, 4, 50), "Clay loam": (0.4, 2, 38),
    "Silty clay loam": (0.6, 1.5, 12), "Sandy clay": (0.6, 1.2, 25),
    "Silty clay": (0.5, 0.9, 5), "Clay": (0.1, 0.6, 12),
}

# Table A9.1: splash detachability EROD, g/J (low, mean, high).  Different
# texture list from A4 (no sandy clay loam / silty clay etc.).
# low = loose dry soil; high = loose moist; mean = sealed or compacted topsoil.
A91 = {
    "Clay": (1.7, 2.0, 2.4), "Clay loam": (1.4, 1.7, 1.9),
    "Silt": (0.8, 1.2, 1.6), "Silt loam": (0.8, 1.5, 2.3),
    "Loam": (1.0, 2.0, 2.7), "Sandy loam": (1.7, 2.6, 3.1),
    "Loamy sand": (1.9, 3.0, 4.0), "Fine sand": (2.0, 3.5, 6.0),
    "Sand": (1.0, 1.9, 3.0),
}

# Table A9.2: cohesion COH, kPa, saturated: (uncompacted low, mean, high,
#                                            compacted low, mean, high) or None
A92 = {
    "Clay": (10, 12, 14, 29, 33, 44), "Clay loam": (9, 10, 14, None, None, None),
    "Silty clay": (9, 15, 11, None, None, None),
    "Silty clay loam": (10, 9, 26, None, None, None),
    "Sandy clay loam": (8, 3, 10, None, None, None),
    "Silt loam": (2, 3, 5, 6, 9, 17), "Loam": (2, 3, 4, 7, 7, 8),
    "Fine sandy loam": (2, 3, 3, 5, 6, 8), "Sandy loam": (2, 2, 4, 4, 7, 10),
    "Loamy sand": (2, 2, 3, 6, 8, 9), "Sand": (2, 2, 3, 8, 8, 9),
}

# Table A3.1: Manning's n  (land use, qualifier, low, mean, high); None = blank
A31 = [
    ("Bare soil", "roughness depth < 25 mm", 0.010, 0.020, 0.030),
    ("Bare soil", "roughness depth 25-50 mm", 0.014, 0.025, 0.033),
    ("Bare soil", "roughness depth 50-100 mm", 0.023, 0.030, 0.038),
    ("Bare soil", "roughness depth > 100 mm", 0.045, 0.047, 0.049),
    ("Bermuda grass, sparse to good cover", "very short grass > 50 mm", 0.015, 0.023, 0.040),
    ("Bermuda grass, sparse to good cover", "short grass 50-100 mm", 0.030, 0.046, 0.060),
    ("Bermuda grass, sparse to good cover", "medium grass 150-200 mm", 0.030, 0.074, 0.085),
    ("Bermuda grass, sparse to good cover", "long grass 250-600 mm", 0.040, 0.100, 0.150),
    ("Bermuda grass, sparse to good cover", "very long grass > 600 mm", 0.060, 0.150, 0.200),
    ("Bermuda grass, dense cover", "", 0.300, 0.410, 0.480),
    ("Other dense sod-forming grasses", "", 0.390, 0.450, 0.630),
    ("Dense bunch grasses", "", None, 0.150, None),
    ("Annual grasses (e.g. Sudan grass)", "", None, 0.200, None),
    ("Kudzu", "", 0.070, 0.150, 0.230),
    ("Lespedeza (legumes)", "", None, 0.100, None),
    ("Natural rangeland", "", 0.100, 0.130, 0.320),
    ("Clipped range", "", 0.020, 0.150, 0.240),
    ("Wheat straw mulch", "2.5 t/ha", 0.050, 0.055, 0.080),
    ("Wheat straw mulch", "5.0 t/ha", 0.075, 0.100, 0.150),
    ("Wheat straw mulch", "7.5 t/ha", 0.100, 0.150, 0.200),
    ("Wheat straw mulch", "10.0 t/ha", 0.130, 0.180, 0.250),
    ("Chopped maize stalks", "2.5 t/ha", 0.012, 0.020, 0.050),
    ("Chopped maize stalks", "5.0 t/ha", 0.020, 0.040, 0.075),
    ("Chopped maize stalks", "10.0 t/ha", 0.023, 0.070, 0.130),
    ("Cotton", "", 0.070, 0.080, 0.090),
    ("Wheat", "", 0.100, 0.125, 0.300),
    ("Sorghum", "", 0.040, 0.090, 0.110),
    ("Mouldboard plough", "", 0.020, 0.060, 0.100),
    ("Chisel plough", "residue < 0.6 t/ha", 0.010, 0.070, 0.170),
    ("Chisel plough", "residue 0.6-2.5 t/ha", 0.070, 0.180, 0.340),
    ("Chisel plough", "residue 2.5-7.5 t/ha", 0.190, 0.300, 0.470),
    ("Chisel plough", "residue > 7.5 t/ha", 0.340, 0.400, 0.460),
    ("Disc/harrow", "residue < 0.6 t/ha", 0.010, 0.080, 0.410),
    ("Disc/harrow", "residue 0.6-2.5 t/ha", 0.100, 0.160, 0.250),
    ("Disc/harrow", "residue 2.5-7.5 t/ha", 0.140, 0.250, 0.530),
    ("Disc/harrow", "residue > 7.5 t/ha", None, 0.300, None),
    ("No tillage", "residue < 0.6 t/ha", 0.030, 0.040, 0.070),
    ("No tillage", "residue 0.6-2.5 t/ha", 0.010, 0.070, 0.130),
    ("No tillage", "residue 2.5-7.5 t/ha", 0.160, 0.300, 0.470),
    ("Coulter", "", 0.050, 0.100, 0.130),
]

# Table A7.1: maximum interception storage of mature plants, DINTR mm
A71 = {
    "Fescue grass": 1.2, "Molinia": 0.2, "Rye grass": 2.5,
    "Meadow grass, clover": 2.3, "Blue stem grass": 2.3, "Heather": 1.5,
    "Bracken": 1.3, "Tropical rain forest": 2.5,
    "Temperate deciduous woodland: winter": 1.0,
    "Temperate deciduous woodland: summer": 2.5,
    "Needleleaf forest: pines": 1.0, "Needleleaf forest: spruce, firs": 1.5,
    "Evergreen hardwood forest": 0.8, "Apple": 0.5, "Soya beans": 0.7,
    "Potatoes": 0.9, "Cabbage": 0.5, "Brussels sprouts": 1.0,
    "Sugar beet": 0.6, "Millet": 0.3, "Spring wheat": 1.8, "Winter wheat": 3.0,
    "Barley, rye, oats": 1.2, "Maize": 0.8, "Tobacco": 1.8, "Alfalfa": 2.8,
}

# Table A7.2 (first page only so far): plant height in METRES (range),
# stem angle degrees (range), SHAPE factor 1 (thin blade) or 2 (broad leaf).
# Model input PLANTH is in CM: multiply the height by 100.
A72 = {
    "Temperate deciduous forest": ((20, 40), (10, 80), 2),
    "Coniferous forest: pine": ((30, 40), (10, 80), 1),
    "Coniferous forest: spruce, fir": ((50, 60), (10, 80), 1),
    "Apple": ((10, 15), (10, 20), 2), "Peach, nectarine": ((6, 7), (40, 60), 2),
    "Citrus": ((6, 12), (10, 80), 2), "Olive": ((12, 15), (30, 40), 2),
    "Banana": ((2, 5), (20, 80), 2), "Grape": ((0.8, 1), (40, 80), 2),
    "Fescue grass": ((0.05, 0.06), (60, 90), 1), "Molinia": ((0.02, 1.2), (75, 80), 1),
    "Rye grass": ((0.1, 0.9), (45, 60), 1), "Timothy grass": ((0.5, 1), (70, 75), 1),
    "Oat grass": ((0.5, 1.5), (20, 90), 1), "Bermuda grass": ((0.3, 0.6), (50, 60), 1),
    "Kikuyu grass": ((0.2, 0.2), (40, 70), 1), "Guinea grass": ((2, 3), (20, 60), 1),
    "Napier grass": ((2, 6), (70, 90), 2), "Rhodes grass": ((0.5, 2), (50, 80), 1),
    "Vetiver grass": ((1, 3), (60, 80), 2), "Prairie grass": ((0.8, 1), (40, 80), 1),
    "Buffel grass": ((0.1, 1), (50, 80), 1), "Elymus": ((0.3, 0.5), (50, 90), 1),
    "Bent grass": ((0.4, 0.5), (60, 80), 1), "Clover": ((0.3, 0.6), (10, 60), 2),
    "Alfalfa, lucerne": ((0.3, 0.9), (50, 70), 2), "Heather": ((0.5, 0.6), (0, 90), 2),
    "Beans (Phaseolus, Vicia)": ((1, 3), (60, 80), 2),
    "Mung bean, Black gram": ((0.3, 1), (60, 80), 2), "Soya bean": ((0.1, 0.2), (20, 40), 2),
    "Pigeon pea, Red gram": ((3, 4), (20, 40), 2), "Chick pea": ((0.4, 0.5), (60, 70), 2),
    "Cotton": ((0.8, 1.2), (0, 20), 2), "Groundnut": ((0.2, 0.6), (40, 80), 2),
    "Hops": ((5, 6), (15, 70), 2), "Maize": ((2, 3), (50, 80), 2),
    "Millet, sorghum": ((1, 2), (50, 80), 2),
    # second page of the table
    "Oilseed rape": ((1, 1.4), (25, 60), 2), "Linseed": ((0.8, 1.6), (60, 90), 2),
    "Pineapple": ((0.5, 1), (70, 90), 2), "Potato": ((0.6, 1), (30, 50), 2),
    "Cassava": ((2.5, 3), (70, 90), 2), "Rice": ((0.5, 1), (70, 80), 1),
    "Sugar beet": ((0.8, 1), (70, 80), 2), "Sugar cane": ((3, 6), (70, 90), 2),
    "Tobacco": ((1.5, 2), (10, 60), 2), "Wheat, barley, oats": ((0.5, 1.5), (80, 90), 1),
    "Rye": ((1, 2), (80, 90), 1), "Rubber": ((18, 30), (20, 80), 2),
    "Oil palm": ((9, 10), (0, 90), 2), "Coffee": ((4, 4.5), (40, 80), 2),
    "Tea": ((1, 1.5), (60, 80), 2), "Cocoa": ((4.5, 7), (60, 80), 2),
    "Coconut": ((18, 30), (0, 90), 2),
}

# Table A7.3: proportional basal area PBASE by land use and cover condition
A73 = [
    ("Fallow: after row crops", "", 0.1), ("Fallow: after sod", "", 0.3),
    ("Row crops", "poor", 0.1), ("Row crops", "good", 0.2),
    ("Small grain", "poor", 0.2), ("Small grain", "good", 0.3),
    ("Hay - legume", "poor", 0.2), ("Hay - legume", "good", 0.4),
    ("Hay - sod", "poor", 0.4), ("Hay - sod", "good", 0.6), ("Hay - sod", "excellent", 0.8),
    ("Pasture or range (bunch grass)", "poor", 0.2), ("Pasture or range (bunch grass)", "fair", 0.3),
    ("Pasture or range (bunch grass)", "good", 0.4),
    ("Temporary pasture - sod", "poor", 0.4), ("Temporary pasture - sod", "fair", 0.5),
    ("Temporary pasture - sod", "good", 0.6),
    ("Permanent pasture or meadow", "poor", 0.8), ("Permanent pasture or meadow", "good", 1.0),
    ("Woods and forest", "", 1.0),
]

# Table A6.1: roughness ratio RFR (cm/m, i.e. per cent) just after tillage (low, high)
A61 = {
    "Mouldboard plough": (30, 33), "Chisel plough": (24, 27), "Cultivator": (15, 23),
    "Tandem disc": (25, 28), "Offset disc": (32, 35), "Paraplow": (32, 35),
    "Spike-tooth harrow": (17, 23), "Spring-tooth harrow": (25, 25), "Rotary hoe": (21, 22),
    "Rototiller": (23, 23), "Drill": (20, 21), "Row planter": (13, 22),
}

# Table A9.3: increase in cohesion from root reinforcement, kPa (low, high), mature plants
A93 = {
    "Barley": (0.2, 0.6), "Grass": (1, 8), "Marram grass": (1.5, 15),
    "Chaparral, matorral": (0.3, 3), "Alfalfa": (10, 10), "Alder": (2, 12),
    "Sitka spruce": (4, 12), "Hemlock": (1, 8), "Willow": (6, 6), "Poplar": (2, 2),
    "Maple": (4, 6), "Pines": (4, 10), "Coniferous forest": (1, 17.5),
    "Candlenut": (15, 35), "Acacia": (1, 5),
}

# Rows printed inconsistently in the guide (kept as printed).
ANOMALIES = {
    "A4.1": {"Silty clay loam": "residual saturation THR 0.380 is far above its neighbours "
                                "(0.015-0.11); probably a misprint (0.038?). Check before relying on it."},
    "A9.2": {"Silty clay": "mean (15) lies above high (11)",
             "Silty clay loam": "mean (9) lies below low (10)",
             "Sandy clay loam": "mean (3) lies below low (8)"},
}

NOTES = {
    "PLANTH": "model input is in cm; Table A7.2 heights are in metres (x100). "
              "Use the height that intercepts raindrops, not the plant top, if ground "
              "cover or litter is dense (guide A7).",
    "A4.1 header": "the column headed 'Maximum saturation (THR)' is THMAX (guide typo).",
    "EROD": "low = loose dry soil, high = loose moist soil, mean = sealed or compacted topsoil (A9.1).",
    "COH": "table values are for bare saturated soil; add a root-reinforcement value from "
           "Table A9.3 when vegetated (A9.3 not yet transcribed).",
}


def soil_preset(texture):
    """Return form-ready values for an A4.1/A4.2 texture, or None."""
    if texture not in A41:
        return None
    p = A41[texture]
    ks = A42[texture]
    return {"porosity": p[1], "theta_max": p[4], "fmin": ks[1], "g": p[6],
            "ranges": {"porosity": (p[0], p[2]), "fmin": (ks[0], ks[2]), "g": (p[5], p[7])},
            "thr": p[3]}


def _m_label(land, qual):
    return land + (": " + qual if qual else "")


def _pb_label(land, cond):
    return land + (" (" + cond + ")" if cond else "")


def guide_json():
    """Everything the form's pickers need, as plain JSON data."""
    out = []
    for t in TEXTURES:
        p, ks = A41[t], A42[t]
        out.append({
            "name": t,
            "porosity": [p[0], p[1], p[2]],
            "theta_max": p[4],
            "thr": p[3],
            "thr_note": ANOMALIES["A4.1"].get(t),
            "fmin": list(ks),
            "g": [p[5], p[6], p[7]],
            "erod": list(A91[t]) if t in A91 else None,
            "coh": list(A92[t]) if t in A92 else None,
            "coh_note": ANOMALIES["A9.2"].get(t),
        })
    return {
        "textures": out,
        "manning": [{"label": _m_label(l, q), "low": lo, "mean": m, "high": hi} for l, q, lo, m, hi in A31],
        "dint": [{"name": k, "value": v} for k, v in A71.items()],
        "veg": [{"name": k, "height": list(v[0]), "angle": list(v[1]), "shape": v[2]} for k, v in A72.items()],
        "pbase": [{"label": _pb_label(l, c), "value": v} for l, c, v in A73],
        "root": [{"name": k, "range": list(v)} for k, v in A93.items()],
        "tillage": [{"name": k, "range": list(v)} for k, v in A61.items()],
    }


CONDITIONS = {"loose_dry": ("Loose and dry", 0), "compacted": ("Sealed or compacted", 1),
              "loose_moist": ("Loose and moist", 2)}
FILL_LABEL = {"fmin": "conductivity", "g": "capillary drive", "porosity": "porosity",
              "theta_max": "maximum moisture", "erodibility": "erodibility", "cohesion": "cohesion"}


def guide_fill(texture, condition="compacted"):
    """The values the form's texture picker puts in a plane (the same rules as
    applyGuide() in app.html).  Returns {field: value}; a field is absent when
    the guide has no usable value for it.  None for an unknown texture."""
    if texture not in A41:
        return None
    p, ks = A41[texture], A42[texture]
    out = {"fmin": ks[1], "g": p[6], "porosity": p[1], "theta_max": p[4]}
    ix = CONDITIONS.get(condition, CONDITIONS["compacted"])[1]
    if texture in A91:
        out["erodibility"] = A91[texture][ix]
    coh = A92.get(texture)
    if coh and texture not in ANOMALIES["A9.2"]:
        if condition == "compacted":
            if coh[3] is not None:
                out["cohesion"] = coh[4]
        else:
            out["cohesion"] = coh[1]
    return out
