"""
eurosem_io.py  --  file interface for EUROSEM v3.2L (planes only)
=================================================================

Purpose
-------
The Fortran model is a file-in / file-out program.  This module is the
Python side of that contract, so a web page (or any script) never has to
touch the fixed-format text files directly:

    inputs  : Run  -> write_par() / write_pcp()  -> .par and .pcp text
              text -> read_par()  / read_pcp()   -> Run   (round trip)
    checks  : validate(run)                       -> list of Issue
    outputs : parse_dyn() / parse_sta()           -> plain dicts / lists
              load_results()                      -> one JSON-ready dict

The .par / .pcp files are "template" files: the Fortran reader skips a
fixed NUMBER of lines between records (READTM) and never looks at what
the skipped lines say.  The writer here therefore reproduces the exact
line counts of the baseline files (testn.par / test.pcp).  The skipped
lines carry the human-readable headings, which we keep so the files stay
readable by a person.

Scope
-----
Planes only.  Channels (W = 0), ponds and lateral inflows are rejected
by validate() with a clear message; they can be added later.

Units follow the model: metres, millimetres, minutes (NTIME = 2).

Requires only the Python standard library (3.8+).
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional

MAX_ELEMENTS = 60          # array limit in the Fortran COMMON blocks
MAX_TIME_STEPS = 1000      # model assumes fewer than ~1000 output steps


# ---------------------------------------------------------------------
#  Data model
# ---------------------------------------------------------------------
@dataclass
class Element:
    """One plane.  Field names follow the Fortran variable names."""
    # --- topology -----------------------------------------------------
    J: int = 1            # element number
    NU: int = 0           # upstream element feeding this one (0 = none)
    NR: int = 0           # lateral inflow (right)  - channels only
    NL: int = 0           # lateral inflow (left)   - channels only
    NC1: int = 0          # channel inflow 1        - channels only
    NC2: int = 0          # channel inflow 2        - channels only
    NPRINT: int = 1       # >0 writes this element's hydrograph (7 = very verbose)
    # --- geometry -----------------------------------------------------
    XL: float = 100.0     # slope length (m)
    W: float = 100.0      # width (m)
    S: float = 0.05       # slope (m/m)  [see note on SIR below]
    ZR: float = 0.0       # channel side slope   (unused for planes)
    ZL: float = 0.0       # channel side slope   (unused for planes)
    BW: float = 0.0       # channel bottom width (unused for planes)
    MANN: float = 0.1     # Manning n, rill flow
    MANN_IR: float = 0.1  # Manning n, interrill flow
    # --- infiltration (Smith-Parlange) --------------------------------
    FMIN: float = 5.0     # saturated conductivity (mm/h)
    G: float = 400.0      # capillary drive (mm)
    POR: float = 0.4      # porosity
    THI: float = 0.15     # initial moisture content
    THMX: float = 0.3     # maximum moisture content
    ROC: float = 0.0      # rock fraction
    RECS: float = 1.0     # random roughness / surface storage (mm)
    DINT: float = 0.9     # interception depth (mm)
    # --- rills --------------------------------------------------------
    DEPNO: float = 0.0    # number of rills
    RILLW: float = 0.0    # rill width (m)
    RILLD: float = 0.0    # rill depth (m)
    ZLR: float = 0.0      # rill side slope
    RS: float = 0.0       # rill spacing / roughness switch (see source)
    RFR: float = 0.0      # rill flow ratio (%)
    SIR: float = 0.05     # interrill slope (m/m)
    # --- cover --------------------------------------------------------
    COVER: float = 0.0    # canopy cover fraction
    SHAPE: float = 2.0    # plant shape code
    PLANGLE: float = 55.0  # leaf drip angle (deg)
    PLANTBASE: float = 0.05  # plant base fraction
    PLANTH: float = 3.0   # plant height (cm)
    DERO: float = 5.0     # effective drop-fall height / erosion depth term
    STONE: float = 1.0    # stone code
    # --- soil / erosion ----------------------------------------------
    D50: float = 25.0     # median particle size (micron)
    EROD: float = 0.9     # soil erodibility (g/J)
    SPLTEX: float = 2.0   # splash texture class
    COH: float = 1.0      # cohesion (kPa)
    RHOS: float = 2.6     # particle density (g/cm3)
    PAVE: float = 0.0     # pavement fraction
    SIGMAS: float = 0.0   # std. dev. of particle size
    MCODE: int = 1        # method code

    @property
    def has_rills(self) -> bool:
        """PLANE treats the element as rilled only if both are > 0."""
        return self.DEPNO > 0 and self.RILLW > 0

    @property
    def effective_slope(self) -> float:
        """Slope PLANE actually uses.  Without rills it overwrites S by
        SIR (eplanel.for: ``s(j) = sir(j)``), so S is ignored."""
        return self.S if self.has_rills else self.SIR


@dataclass
class Gage:
    ident: str = "gage 1"
    times: List[float] = field(default_factory=lambda: [0.0, 60.0])   # min
    depths: List[float] = field(default_factory=lambda: [0.0, 30.0])  # mm, cumulative


@dataclass
class Run:
    title: str = "EUROSEM run"
    # system
    NPART: int = 0
    CLEN: float = 100.0   # characteristic length (m)
    TFIN: float = 60.0    # finish time (min)
    DELT: float = 0.5     # time step (min)
    THETA: float = 0.7    # implicit weighting
    TEMP: float = 20.0    # water temperature (deg C)
    # options
    NTIME: int = 2        # 1 = seconds, 2 = minutes
    NEROS: int = 2        # >=1 switches sediment on
    # topology
    order: List[int] = field(default_factory=list)   # element numbers, computation order
    elements: List[Element] = field(default_factory=list)
    # rainfall
    gages: List[Gage] = field(default_factory=list)
    weights: Dict[int, tuple] = field(default_factory=dict)  # element -> (gage, weight)
    MAXND: Optional[int] = None

    @property
    def NELE(self) -> int:
        return len(self.elements)

    def element(self, j: int) -> Element:
        for e in self.elements:
            if e.J == j:
                return e
        raise KeyError(j)


# ---------------------------------------------------------------------
#  Number formatting
# ---------------------------------------------------------------------
def _f(x: float) -> str:
    """Compact, Fortran list-directed-safe real (no exponent surprises)."""
    s = repr(float(x))
    return s[:-2] if s.endswith(".0") else s   # '5.0' -> '5'; list-directed reads either


def _row(*vals) -> str:
    out = []
    for v in vals:
        out.append(str(int(v)) if isinstance(v, int) and not isinstance(v, bool) else _f(v))
    return "    " + "     ".join(out)


# ---------------------------------------------------------------------
#  Writers
# ---------------------------------------------------------------------
def write_par(run: Run) -> str:
    """Return the text of the parameter file.  Line counts between
    records match testn.par exactly (the Fortran reader counts lines)."""
    L: List[str] = []
    A = L.append
    A(f" MEUROSEM V. 3/93 Parameter Input File {run.title}"[:79])
    A("#")
    A("**********************************")
    A("***********  S Y S T E M  ********")
    A("**********************************")
    A("*  NELE   NPART  CLEN(M)  TFIN(min) DELT(min)  THETA   TEMP")
    A(_row(run.NELE, run.NPART, run.CLEN, run.TFIN, run.DELT, run.THETA, run.TEMP))
    A("#")
    A("***********************************")
    A("********* O P T I O N S ***********")
    A("***********************************")
    A("  NTIME   NEROS")
    A(_row(run.NTIME, run.NEROS))
    A("#")
    A("***********************************************")
    A("****   C O M P U T A T I O N   O R D E R   ****")
    A("***********************************************")
    A("  There must be NELE elements in the list. NLOG")
    A("  must be sequential. ELEMENT NUM. need not be.")
    A("#")
    A("  COMP. ORDER    ELEMENT")
    A("     (NLOG)      NUM. (J)")
    A("     ------      --------")
    for n, j in enumerate(run.order, start=1):
        A(_row(n, j))
    A("#")
    A("***********************************************")
    A("******   E L E M E N T - W I S E   I N F O  ***")
    A("***********************************************")
    A("    There must be NELE sets of the ELEMENT-WISE  prompts and data")
    A("    records; duplicate records from * to * for each element.  The")
    A("    elements may be entered in any order.")
    dash = "-" * 67
    for e in run.elements:
        A(f"*      element {e.J}")
        A(" J     NU        NR       NL      NC1      NC2     NPRINT")
        A(_row(e.J, e.NU, e.NR, e.NL, e.NC1, e.NC2, e.NPRINT))
        A(dash)
        A("     XL(M)    W(M)     S     ZR     ZL    BW(M)  MANN(Rill) Mann(IR)")
        A(_row(e.XL, e.W, e.S, e.ZR, e.ZL, e.BW, e.MANN, e.MANN_IR))
        A(dash)
        A("     FMIN(mm/h)  G(mm)  POR  ThI   ThMX    ROC   RECS(mm) DINT(mm) ")
        A(_row(e.FMIN, e.G, e.POR, e.THI, e.THMX, e.ROC, e.RECS, e.DINT))
        A(dash)
        A("     DEPNO    RILLW(m)   RILLD(m)   ZLR    RS      RFR   SIR")
        A(_row(e.DEPNO, e.RILLW, e.RILLD, e.ZLR, e.RS, e.RFR, e.SIR))
        A(dash)
        A("     COVER     SHAPE   PLANGLE   PLANTBASE   PLANTH(cm)  Dero  stone")
        A(_row(e.COVER, e.SHAPE, e.PLANGLE, e.PLANTBASE, e.PLANTH, e.DERO, e.STONE))
        A(dash)
        A("     D50(u)    EROD   SPLTEX   COH    RHOS   PAVE   SIGMAS  MCODE")
        A(_row(e.D50, e.EROD, e.SPLTEX, e.COH, e.RHOS, e.PAVE, e.SIGMAS, e.MCODE))
        A("*")
    return "\n".join(L) + "\n"


def write_pcp(run: Run) -> str:
    """Return the text of the rainfall file.

    The Fortran reader keeps only the first ND pairs of each gage.  If
    the last time is before TFIN it adds a closing point (TFIN + 0.05
    min, same cumulative depth), so rain simply stops.  MAXND must be
    >= the largest ND."""
    maxnd = max(len(g.times) for g in run.gages)
    if run.MAXND is not None:
        maxnd = max(maxnd, run.MAXND)
    L: List[str] = []
    A = L.append
    A("   EUROSEM   Rainfall Input Data V3.5 " + run.title)
    A("#")
    A("****************************")
    A("      Gage Network Data")
    A("****************************")
    A("#")
    A("  NUM. OF RAINGAGES    MAX. NUM. OF TIME-DEPTH DATA PAIRS FOR ALL GAGES")
    A("      (NGAGES)            (MAXND)")
    A("      --------            -------")
    A(_row(len(run.gages), maxnd))
    A("#")
    A(" There must be NELE  pairs of (GAGE  WEIGHT) data")
    A("*")
    A(" ELE. NUM. (J)     RAINGAGE     WEIGHT")
    A(" -------------     --------     ------")
    for e in run.elements:
        g, w = run.weights.get(e.J, (1, 1.0))
        A(_row(e.J, g, w))
    A("#")
    A("*****************************")
    A("       Rainfall Data")
    A("*****************************")
    A(" There must be NGAGES sets of rainfall data. Repeat lines from * to *")
    A(" for each gage inserting a variable number of TIME-DEPTH data pairs")
    A(" (see example in User Manual).")
    A("#")
    for n, g in enumerate(run.gages, start=1):
        A(f"* ALPHA-NUMERIC GAGE ID: {g.ident}")
        A("#")
        A("  GAGE NUM.     NUM. OF DATA PAIRS (ND)")
        A("  ---------     -----------------------")
        A(_row(n, len(g.times)))
        A("#")
        A(" There must be ND pairs of time-depth (T D) data: NOTE: The last time")
        A(" must be greater than TFIN (the total computational time).")
        A("#")
        A("  TIME(min)  ACCUM. DEPTH(mm)")
        A("  ------    ------------")
        for t, d in zip(g.times, g.depths):
            A(_row(t, d))
        A("*")
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------
#  Readers (mirror the Fortran line-counting logic)
# ---------------------------------------------------------------------
class _Lines:
    """Cursor that skips N lines then returns the numbers on the next."""
    def __init__(self, text: str):
        self.lines = text.replace("\r", "").replace("\x1a", "").split("\n")
        self.i = 0
        self.last_skipped: List[str] = []

    def skip(self, n: int):
        self.last_skipped = self.lines[self.i:self.i + n]
        self.i += n

    def record(self, skip: int, expect: int = 0) -> List[float]:
        self.skip(skip)
        if self.i >= len(self.lines):
            raise ValueError("file ended early (line %d)" % self.i)
        raw = self.lines[self.i]
        self.i += 1
        toks = raw.replace(",", " ").split()
        try:
            vals = [float(t.replace("D", "E").replace("d", "e")) for t in toks]
        except ValueError:
            raise ValueError("line %d is not numeric: %r" % (self.i, raw))
        if expect and len(vals) < expect:
            raise ValueError("line %d has %d values, expected %d: %r"
                             % (self.i, len(vals), expect, raw))
        return vals


def read_par(text: str, run: Optional[Run] = None) -> Run:
    """Parse a parameter file exactly as READER/PAREAD would."""
    run = run or Run()
    c = _Lines(text)
    run.title = c.lines[0].split("File", 1)[-1].strip() if "File" in c.lines[0] else c.lines[0].strip()
    v = c.record(6, 7)
    nele = int(v[0]); run.NPART = int(v[1])
    run.CLEN, run.TFIN, run.DELT, run.THETA, run.TEMP = v[2:7]
    v = c.record(5, 2); run.NTIME, run.NEROS = int(v[0]), int(v[1])
    run.order = []
    c.skip(10 - 1)           # record() below skips the last of the 10
    for k in range(nele):
        v = c.record(1 if k == 0 else 0, 2)
        run.order.append(int(v[1]))
    c.skip(7)
    run.elements = []
    for _ in range(nele):
        e = Element()
        v = c.record(2, 7); (e.J, e.NU, e.NR, e.NL, e.NC1, e.NC2, e.NPRINT) = [int(x) for x in v[:7]]
        v = c.record(2, 8); (e.XL, e.W, e.S, e.ZR, e.ZL, e.BW, e.MANN, e.MANN_IR) = v[:8]
        v = c.record(2, 8); (e.FMIN, e.G, e.POR, e.THI, e.THMX, e.ROC, e.RECS, e.DINT) = v[:8]
        v = c.record(2, 7); (e.DEPNO, e.RILLW, e.RILLD, e.ZLR, e.RS, e.RFR, e.SIR) = v[:7]
        v = c.record(2, 7); (e.COVER, e.SHAPE, e.PLANGLE, e.PLANTBASE, e.PLANTH, e.DERO, e.STONE) = v[:7]
        v = c.record(2, 8); (e.D50, e.EROD, e.SPLTEX, e.COH, e.RHOS, e.PAVE, e.SIGMAS) = v[:7]; e.MCODE = int(v[7])
        c.skip(1)             # the terminating '*'
        run.elements.append(e)
    return run


def read_pcp(text: str, run: Optional[Run] = None) -> Run:
    """Parse a rainfall file.  Like the Fortran, keeps only the first ND
    pairs of each gage (extra pairs listed in the file are ignored)."""
    run = run or Run()
    c = _Lines(text)
    v = c.record(9, 2); ngages, run.MAXND = int(v[0]), int(v[1])
    nele = run.NELE if run.elements else 0
    run.weights = {}
    n_ele_records = nele or 0
    if n_ele_records == 0:
        raise ValueError("read the .par first so NELE is known")
    for k in range(n_ele_records):
        v = c.record(5 if k == 0 else 0, 3)
        run.weights[int(v[0])] = (int(v[1]), v[2])
    c.skip(8)
    run.gages = []
    for _ in range(ngages):
        v = c.record(4, 2); nd = int(v[1])
        ident = c.last_skipped[0].split(":", 1)[-1].strip() if c.last_skipped else ""
        g = Gage(ident=ident, times=[], depths=[])
        for k in range(nd):
            p = c.record(6 if k == 0 else 0, 2)
            g.times.append(p[0]); g.depths.append(p[1])
        c.skip(1)
        run.gages.append(g)
    return run


# ---------------------------------------------------------------------
#  Validation
# ---------------------------------------------------------------------
@dataclass
class Issue:
    level: str            # 'error' (do not run) | 'warning' | 'note'
    element: Optional[int]
    field: str
    message: str

    def __str__(self):
        w = f"element {self.element}: " if self.element is not None else ""
        return f"[{self.level.upper()}] {w}{self.field}: {self.message}"


def validate(run: Run) -> List[Issue]:
    """Check a run before it is sent to the model.  Any 'error' means
    the run must not be started."""
    out: List[Issue] = []
    err = lambda j, f, m: out.append(Issue("error", j, f, m))
    warn = lambda j, f, m: out.append(Issue("warning", j, f, m))
    note = lambda j, f, m: out.append(Issue("note", j, f, m))

    # ---- system ------------------------------------------------------
    n = run.NELE
    if n < 1: err(None, "NELE", "at least one element is needed")
    if n > MAX_ELEMENTS: err(None, "NELE", f"model limit is {MAX_ELEMENTS} elements")
    if run.TFIN <= 0: err(None, "TFIN", "must be > 0")
    if run.DELT <= 0: err(None, "DELT", "must be > 0")
    if run.DELT > 0 and run.TFIN / run.DELT >= MAX_TIME_STEPS:
        err(None, "DELT", f"TFIN/DELT = {run.TFIN/run.DELT:.0f} steps; model limit is < {MAX_TIME_STEPS}")
    if not (0.5 <= run.THETA <= 1.0):
        warn(None, "THETA", "implicit weighting is normally 0.5-1.0 (0.7 typical)")
    if run.NTIME not in (1, 2):
        err(None, "NTIME", "must be 1 (seconds) or 2 (minutes); the model rejects 3")
    if run.NEROS < 1:
        note(None, "NEROS", "sediment is switched off (NEROS < 1): no erosion output")
    if run.NPART != 0:
        warn(None, "NPART", "particle-size classes are untested in this interface; use 0")

    # ---- topology ----------------------------------------------------
    ids = [e.J for e in run.elements]
    if len(set(ids)) != len(ids): err(None, "J", "element numbers must be unique")
    if sorted(run.order) != sorted(ids):
        err(None, "order", "computation order must list every element exactly once")
    pos = {j: k for k, j in enumerate(run.order)}
    for e in run.elements:
        j = e.J
        if any((e.NR, e.NL, e.NC1, e.NC2)):
            err(j, "NR/NL/NC1/NC2", "lateral and channel inflows are not supported yet (planes only)")
        if e.NU != 0:
            if e.NU not in ids:
                err(j, "NU", f"upstream element {e.NU} does not exist")
            elif e.NU == j:
                err(j, "NU", "an element cannot feed itself")
            elif pos.get(e.NU, -1) > pos.get(j, -1):
                err(j, "NU", f"upstream element {e.NU} must come EARLIER in the computation order")
        if e.NPRINT <= 0:
            note(j, "NPRINT", "no hydrograph will be written for this element; use 1")
        elif e.NPRINT > 1:
            warn(j, "NPRINT", f"NPRINT={e.NPRINT} writes very verbose diagnostics; 1 is recommended")
    feeders = [e.NU for e in run.elements if e.NU]
    dup = {u for u in feeders if feeders.count(u) > 1}
    for u in sorted(dup):
        warn(u, "NU", "feeds more than one element; only the last element should be the outlet - check layout")

    # ---- per-plane physical checks -------------------------------------
    for e in run.elements:
        j = e.J
        if e.W <= 0:
            err(j, "W", "width must be > 0 (W = 0 means a channel, which is not supported yet)")
        if e.BW or e.ZR or e.ZL:
            err(j, "BW/ZR/ZL", "channel geometry present; planes only")
        if e.XL <= 0: err(j, "XL", "slope length must be > 0")
        if e.S <= 0 and e.has_rills: err(j, "S", "slope must be > 0 when rills are present")
        if e.S <= 0: warn(j, "S", "slope is <= 0 (ignored without rills, but check the value)")
        if e.S > 0.5: warn(j, "S", f"slope {e.S} is very steep (> 0.5); kinematic wave may be unreliable")
        if e.SIR <= 0: err(j, "SIR", "interrill slope must be > 0")
        if not e.has_rills:
            if abs(e.SIR - e.S) > 1e-9:
                note(j, "SIR", f"no rills: the model uses SIR = {e.SIR} as the plane slope and "
                               f"ignores S = {e.S}. Set them equal to avoid confusion.")
        else:
            if e.SIR <= e.S:
                warn(j, "SIR", "with rills SIR must exceed S; the model will reset SIR = 1.4*S")
            if e.MANN <= 0: err(j, "MANN", "rill Manning n must be > 0")
        if not (0.005 <= e.MANN <= 1.0): warn(j, "MANN", "Manning n outside 0.005-1.0")
        if not (0.005 <= e.MANN_IR <= 1.0): warn(j, "Mann(IR)", "Manning n outside 0.005-1.0")
        if e.FMIN < 0: err(j, "FMIN", "must be >= 0")
        if e.G <= 0: err(j, "G", "capillary drive must be > 0")
        if not (0 < e.POR < 1): err(j, "POR", "porosity must be between 0 and 1")
        if not (0 <= e.THI <= e.THMX <= e.POR):
            err(j, "THI/THMX", "need 0 <= THI <= THMX <= POR")
        if not (0 <= e.ROC < 1): err(j, "ROC", "rock fraction must be in [0, 1)")
        if e.RECS < 0: err(j, "RECS", "must be >= 0")
        if e.DINT < 0: err(j, "DINT", "must be >= 0")
        if not (0 <= e.COVER <= 1): err(j, "COVER", "cover must be in [0, 1]")
        if not (0 <= e.PLANGLE <= 90): err(j, "PLANGLE", "leaf angle must be 0-90 degrees")
        if e.PLANTH < 0: err(j, "PLANTH", "must be >= 0")
        if e.DEPNO < 0 or e.RILLW < 0 or e.RILLD < 0:
            err(j, "DEPNO/RILLW/RILLD", "must be >= 0")
        if (e.DEPNO > 0) != (e.RILLW > 0):
            warn(j, "DEPNO/RILLW", "rills need both a number and a width; otherwise the plane is run without rills")
        if run.NEROS >= 1:
            if e.D50 <= 0: err(j, "D50", "must be > 0 micron")
            if e.EROD < 0: err(j, "EROD", "must be >= 0")
            if e.COH < 0: err(j, "COH", "must be >= 0")
            if not (2.0 <= e.RHOS <= 3.0): warn(j, "RHOS", "particle density outside 2.0-3.0")
            if not (0 <= e.PAVE <= 1): err(j, "PAVE", "pavement fraction must be in [0, 1]")

    # ---- rainfall ----------------------------------------------------
    if not run.gages:
        err(None, "gages", "at least one raingage is needed")
    if len(run.gages) > 20:
        err(None, "gages", "model limit is 20 raingages")
    for gi, g in enumerate(run.gages, start=1):
        if len(g.times) != len(g.depths) or len(g.times) < 2:
            err(None, f"gage {gi}", "need at least two (time, depth) pairs of equal length")
            continue
        if g.times[0] != 0: warn(None, f"gage {gi}", "first time should be 0")
        if any(b <= a for a, b in zip(g.times, g.times[1:])):
            err(None, f"gage {gi}", "times must be strictly increasing")
        if any(b < a for a, b in zip(g.depths, g.depths[1:])):
            err(None, f"gage {gi}", "cumulative depth must not decrease")
        if g.times[-1] < run.TFIN:
            note(None, f"gage {gi}", f"rain data end at {g.times[-1]} min, before TFIN = {run.TFIN}: the model "
                                     "holds the depth constant, i.e. no rain after the last point")
        if len(g.times) > 99:
            err(None, f"gage {gi}", "model limit is 99 time-depth pairs per gage")
        if g.depths[-1] <= 0: warn(None, f"gage {gi}", "no rain: there will be no runoff")
    for e in run.elements:
        gw = run.weights.get(e.J)
        if gw is None: err(e.J, "weights", "no raingage assigned")
        elif not (1 <= gw[0] <= len(run.gages)): err(e.J, "raingage", f"gage {gw[0]} does not exist")
        elif gw[1] <= 0: err(e.J, "weight", "gage weight must be > 0")
    return out


def has_errors(issues: List[Issue]) -> bool:
    return any(i.level == "error" for i in issues)


# ---------------------------------------------------------------------
#  Output parsers
# ---------------------------------------------------------------------
_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?"


def _num(s: str) -> float:
    return float(s)


def parse_dyn(text: str) -> dict:
    """Parse the dynamic (.dyn) output file.

    Returns {'title', 'par_file', 'pcp_file', 'elements': {j: {...}},
             'global': {...}}.  Each element holds the summary-table
    values and (if NPRINT > 0) the hydrograph as column lists:
    time_min, rain_mmh, q_m3min, q_mmh, conc, qs_kgmin."""
    L = text.replace("\r", "").split("\n")
    res: dict = {"elements": {}, "global": {}}
    for ln in L[:5]:
        m = re.search(r"INPUT PARAMETER FILE:\s*(\S+)", ln)
        if m: res["par_file"] = m.group(1)
        m = re.search(r"INPUT RAINFALL FILE:\s*(\S+)", ln)
        if m: res["pcp_file"] = m.group(1)
    for k, ln in enumerate(L[:12]):
        if "DESCRIPTIVE RUN TITLE" in ln:
            res["title"] = L[k + 1].strip()
            break

    # per-element balance rows:  "  1   PLANE   0.712E-01   16530.156"
    row = re.compile(r"^\s*(\d+)\s+(PLANE|CHANNEL|POND\w*|\w+)\s+(%s)\s+(%s)\s*$" % (_NUM, _NUM))
    for ln in L:
        m = row.match(ln)
        if m:
            j = int(m.group(1))
            res["elements"].setdefault(j, {}).update(
                type=m.group(2), vol_error_pct=_num(m.group(3)), sediment_total_kg=_num(m.group(4)))

    i = 0
    while i < len(L):
        m = re.search(r"HYDROGRAPH FOR ELEMENT\s+(\d+)", L[i])
        if not m:
            i += 1
            continue
        j = int(m.group(1))
        el = res["elements"].setdefault(j, {})
        cols = {k: [] for k in ("time_min", "rain_mmh", "q_m3min", "q_mmh", "conc", "qs_kgmin")}
        i += 1
        m = re.search(r"CONTRIBUTING AREA =\s*(%s) SQ\. METER OR\s*(%s)" % (_NUM, _NUM), L[i])
        if m:
            el["contributing_area_m2"] = _num(m.group(1))
            el["contributing_area_ha"] = _num(m.group(2))
        while i < len(L) and "TIME TO PEAK FLOW RATE" not in L[i]:
            toks = L[i].split()
            if len(toks) == 6:
                try:
                    vals = [float(t) for t in toks]
                except ValueError:
                    vals = None
                if vals:
                    for k, v in zip(cols, vals):
                        cols[k].append(v)
            i += 1
        el["hydrograph"] = cols
        while i < len(L) and not re.match(r"\s*PEAK FLOW RATE =", L[i]):
            m = re.search(r"TIME TO PEAK FLOW RATE =\s*(%s)" % _NUM, L[i])
            if m: el["time_to_peak_min"] = _num(m.group(1))
            i += 1
        if i < len(L):
            m = re.search(r"PEAK FLOW RATE =\s*(%s)" % _NUM, L[i])
            if m: el["peak_flow_mmh"] = _num(m.group(1))
        i += 1

    g = res["global"]
    def grab(key, pat):
        for ln in L:
            m = re.search(pat, ln)
            if m:
                g[key] = [_num(x) for x in m.groups()]
                if len(g[key]) == 1: g[key] = g[key][0]
                return
    grab("total_rainfall_mm", r"TOTAL RAINFALL DEPTH =\s*(%s)" % _NUM)
    grab("basin_area_m2", r"BASIN AREA =\s*(%s)" % _NUM)
    grab("storage_planes_mm", r"STORAGE REMAINING ON ALL PLANES\s*=\s*(%s)" % _NUM)
    grab("storage_channels_mm", r"STORAGE REMAINING IN CHANNELS\+CONDUITS\s*=\s*(%s)" % _NUM)
    grab("storage_ponds_mm", r"STORAGE REMAINING IN PONDS\s*=\s*(%s)" % _NUM)
    grab("infiltration_planes_mm", r"TOTAL INFILTRATION FROM ALL PLANES\s*=\s*(%s)" % _NUM)
    grab("infiltration_channels_mm", r"TOTAL INFILTRATION FROM ALL CHANNELS\s*=\s*(%s)" % _NUM)
    grab("runoff_mm_m3", r"TOTAL BASIN RUNOFF\s*=\s*(%s) \(MM\)\s*(%s)" % (_NUM, _NUM))
    grab("sum_terms_mm", r"TOTAL OF STOR\., INFIL\. AND RUNOFF TERMS\s*=\s*(%s)" % _NUM)
    grab("vol_error_pct", r"GLOBAL VOL\. ERROR\s*=\s*(%s)" % _NUM)
    if isinstance(g.get("runoff_mm_m3"), list):
        g["runoff_mm"], g["runoff_m3"] = g.pop("runoff_mm_m3")
    return res


def parse_sta(text: str) -> dict:
    """Parse the static summary (.sta) file: echoed input parameters,
    erosion and hydrology summary per element, global balance."""
    L = text.replace("\r", "").split("\n")
    res: dict = {"elements": {}, "global": {}}
    for ln in L:
        m = re.search(r"RUN TITLE:\s*(.*)$", ln)
        if m:
            res["title"] = m.group(1).strip()
            break

    # split into element sections
    starts = [(k, int(re.search(r"ELEMENT\s+(\d+)", ln).group(1)))
              for k, ln in enumerate(L) if "INPUT DATA FOR ELEMENT" in ln]
    gstart = next((k for k, ln in enumerate(L) if ln.strip().startswith("GLOBAL VOLUME BALANCE")), len(L))
    for n, (k, j) in enumerate(starts):
        end = starts[n + 1][0] if n + 1 < len(starts) else gstart
        block = L[k:end]
        el = {"input": {}, "erosion": {}, "hydrology": {}}
        # echoed inputs: "LABEL:  value" pairs; units follow the number
        for ln in block:
            if "Derived parameters" in ln or "EROSION SUMMARY" in ln:
                pass
            for lab, val in re.findall(r"([A-Za-z0-9()]+):\s*(%s)" % _NUM, ln):
                el["input"].setdefault(lab, _num(val))
            m = re.search(r"SurfStor:\s*(%s)" % _NUM, ln)
            if m: el["input"]["SurfStor_mm"] = _num(m.group(1))
        def f(key, pat):
            for ln in block:
                m = re.match(r"\s*" + pat, ln)
                if m:
                    return m.groups()
            return None
        m = f("rill", r"GROSS RILL EROSION\s+(%s) kg\s+(%s) t/ha" % (_NUM, _NUM))
        if m: el["erosion"].update(gross_rill_kg=_num(m[0]), gross_rill_tha=_num(m[1]))
        m = f("gross", r"GROSS INTERRILL EROSION\s+(%s) kg\s+(%s) t/ha" % (_NUM, _NUM))
        if m: el["erosion"].update(gross_interrill_kg=_num(m[0]), gross_interrill_tha=_num(m[1]))
        m = f("net", r"NET EROSION/DEPOSITION\s+(%s) kg\s+(%s) t/ha" % (_NUM, _NUM))
        if m: el["erosion"].update(net_kg=_num(m[0]), net_tha=_num(m[1]))
        for key, pat in (
            ("net_rainfall_mm", r"NET RAINFALL\s*=\s*(%s)"),
            ("peak_rain_mmh", r"PEAK RAINFALL RATE\s*=\s*(%s)"),
            ("time_to_runoff_min", r"TIME TO RUNOFF\s*=\s*(%s)"),
            ("runoff_duration_min", r"DURATION OF RUNOFF\s*=\s*(%s)"),
            ("time_to_peak_flow_min", r"TIME TO PEAK FLOW RATE\s*=\s*(%s)"),
            ("peak_flow_mmh", r"PEAK FLOW RATE\s*=\s*(%s)"),
            ("time_to_peak_sed_min", r"TIME TO PEAK SEDIMENT DISCHARGE\s*=\s*(%s)"),
            ("peak_sed_kgmin", r"PEAK SEDIMENT DISCHARGE\s*=\s*(%s)"),
        ):
            m = f(key, pat % _NUM)
            if m: el["hydrology"][key] = _num(m[0])
        res["elements"][j] = el

    gb = parse_dyn("\n".join(L[gstart:]))["global"]
    res["global"] = gb
    return res


def load_results(dyn_text: str, sta_text: str) -> dict:
    """Combine .dyn and .sta into one JSON-ready dict for the web page.
    Adds a few derived numbers: outflow volume (m3), area, and a
    cumulative sediment curve for each element."""
    dyn, sta = parse_dyn(dyn_text), parse_sta(sta_text)
    out = {"title": sta.get("title") or dyn.get("title", ""),
           "elements": {}, "global": {**dyn["global"], **sta["global"]}}
    for j in sorted(set(dyn["elements"]) | set(sta["elements"])):
        d, s = dyn["elements"].get(j, {}), sta["elements"].get(j, {})
        el = {**{k: v for k, v in d.items() if k != "hydrograph"}, **s}
        h = d.get("hydrograph")
        if h and len(h["time_min"]) > 1:
            t, q, qs = h["time_min"], h["q_m3min"], h["qs_kgmin"]
            vol = sed = 0.0
            cum_sed = [0.0]
            for k in range(1, len(t)):               # trapezoid rule
                dt = t[k] - t[k - 1]
                vol += 0.5 * (q[k] + q[k - 1]) * dt
                sed += 0.5 * (qs[k] + qs[k - 1]) * dt
                cum_sed.append(sed)
            el["outflow_m3"] = vol
            el["outflow_sediment_kg"] = sed
            h = {**h, "cum_sediment_kg": cum_sed}
            el["hydrograph"] = h
        out["elements"][str(j)] = el
    return out


# ---------------------------------------------------------------------
if __name__ == "__main__":      # quick command-line use:  python eurosem_io.py file.dyn file.sta
    import sys
    if len(sys.argv) == 3:
        r = load_results(open(sys.argv[1]).read(), open(sys.argv[2]).read())
        print(json.dumps(r, indent=1)[:4000])
    else:
        print(__doc__)
