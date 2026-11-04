# EUROSEM planes: local web app

Runs EUROSEM v3.2L from a browser form and shows the results as charts.
Run locally, everything stays on your computer (see "Public (hosted) copy" for the web version). Planes only (a cascade of up to 8),
each with optional user-defined rills; channels and ponds are not in the form yet.

**Disclaimer.** Results are model estimates and depend on the parameters entered. They are provided without warranty and should not be the sole basis for engineering, safety, regulatory or financial decisions. You are responsible for checking the inputs and for any use of the results.

## First time

    cd ~/Downloads/eurosem/eurosem_web
    bash build_exe.sh                  # compiles the model (stage 4) to ./eurosem; needs gfortran
    python3 test_real_model.py         # checks the form's defaults against the reference run

## Start the app

    python3 server.py
    # then open http://127.0.0.1:8080 in a browser. Ctrl-C stops it.

Options: `--port 9000`, `--exe /path/to/eurosem`, `--timeout 60`.
Without `--host` the service answers only on this computer. After changing
any file, stop the server (Ctrl-C), start it again and hard-refresh the page
(Safari: Cmd-Option-R; Chrome: Cmd-Shift-R).

If the page looks like plain text with a "Summary" box, Safari has switched
to Reader view: click the page icon at the left of the address bar.

## Using it

1. **Event and rainfall.** Title, duration, time step, water temperature, and
   a table of cumulative rain depth against time (or fill it from a constant
   rate).
2. **Planes.** One card per plane. "Receives runoff from" links the cascade;
   the last plane is the outlet. Every field has a "?" button with the
   guide's definition and the table it comes from.
3. **Run model.** Results: outlet summary tiles, four linked charts
   (rainfall, discharge, sediment discharge, cumulative sediment), a table per
   plane and the global water balance.
4. **Download.** "Results (CSV, wide)" has comment rows at the top recording
   every setting used, then one row per time step. "Settings only" is the
   same record without results.

## Saving and reloading settings

"Save settings to a file" (next to Run model) downloads every setting on the page (event, rainfall
table, all planes, and the guide choices such as soil texture) as a small `.json` file.
"Load settings from a file" fills the whole form from such a file, so a run can be repeated or
adjusted without retyping. This happens in the browser only: nothing is uploaded to the server or
stored there, so it works the same on the hosted copy. A file is checked before it is used: it must
be an EUROSEM settings file, every number must be a number, names must be ones the form knows,
planes must link up, and missing fields come from the example case. If anything is wrong the form
is left unchanged and the problems are listed. Unknown fields are ignored with a note. After loading,
check the values and press Run model. The CSV downloads are for records and plotting and cannot be
loaded back; use the settings file for that. Sensitivity ranges are not saved.

## Public (hosted) copy

`python3 server.py --public` (or `EUROSEM_PUBLIC=1`) is the setting for a hosted copy: sensitivity is
off, runs are rate-limited per visitor and overall, and it listens on all addresses on `$PORT`.
The Dockerfile and hosting steps are in `../DEPLOY.md`.

## Rills

Each plane has a **Rills** switch (under Geometry). Off, the whole plane is one interrill
surface, as before. On, you give the rill geometry and the model routes interrill flow into
the rills and reports rill and interrill erosion separately.

* **Rills across the width** (DEPNO), **bottom width** and **depth** (m), and **side slope**
  (1 : x). Spacing is width / number. The rills' top width (bottom width + 2 x side slope
  x depth) times the number must be less than the plane width, or the model stops; the
  form says so as you type.
* **Slope** in Geometry becomes the slope along the rills. **Interrill slope** must be
  steeper; 0 means the model's own default, 1.4 x the rill slope.
* **Rill Manning n**: 0 uses the plane's Manning n for both.
* **Depth along the plane**: "same depth" or "deepens downslope" (the depth you enter is then the
  value at the bottom; the model applies this only to a plane that receives no runoff from
  another plane). The guide's text describes this switch the opposite way round; the form
  follows the model source and Chapter 3.
* The CSV records every rill setting, and the CSV summary has gross rill and interrill erosion rows
  (rill erosion is blank for planes without rills). Sensitivity runs work with rills on (slope then changes the rill slope and,
  if left at 0, the interrill slope with it).
* Reference numbers from real-model runs (reference two-plane case, rills on plane 1 only;
  `python3 test_real_rills.py` reproduces them):

  | Case | Runoff mm | Outlet sediment kg | Plane 1 rill / interrill kg |
  |---|---|---|---|
  | No rills (reference) | 26.84 | 7,878 | n/a / 16,530 |
  | 20 rills, 0.3 m wide, 0.15 m deep (default) | 30.15 | 14,326 | 38,680 / 22,817 |
  | same, depth deepens downslope | 30.16 | 14,338 | 38,751 / 22,818 |
  | 10 rills, 0.5 m wide, 0.25 m deep | 29.97 | 14,089 | 33,746 / 25,703 |
  | default, rill Manning n 0.05 | 30.60 | 18,014 | 62,263 / 22,818 |
  | default, interrill slope 0.2 | 30.22 | 15,973 | 26,913 / 47,244 |

  The default interrill slope (1.4 x the rill slope, 0.133 here) is steeper than the no-rills
  plane (0.095), so the first step down the table mixes the rill effect with a steeper interrill
  surface. For a like-for-like comparison set the interrill slope just above the rill slope.
  Rill depth pattern hardly matters; rill Manning n and interrill slope matter most here.
  Volume errors were under 0.25 %. Tiny rills (100 of 0.1 x 0.05 m) stopped the model.
* If deposition fills a rill completely the model stops ("excess rill fill"). In the first real test 100 rills of 0.1 m x 0.05 m did this, while 20 rills of 0.3 m x 0.15 m ran. The default is the second.
* Checked against the real model (numbers above).

## Surface roughness (RFR) and depression storage

"Roughness ratio (RFR)", in cm/m, is the downslope roughness of the surface (guide Appendix 6).
The model turns it into depression storage D = exp(-6.66 + 0.27 RFR) mm, which has to fill before
runoff starts, and, above about 2.9, it also raises the interrill Manning n (a roughness height
of 60 RFR/29 - 6 mm is combined with n; this comes from the source, not the guide). The form
shows both effects live. Warning: if the Manning n you chose already allows for a rough tilled
surface, the roughness is then counted twice. The default is 0 (smooth), as in the reference
run, which gives about 0.001 mm of storage. A tillage picker fills Table A6.1 mid-values; those
are for freshly tilled soil, and equations A6.1 (clay) and A6.2 (rain since tillage) are given in
the "?" help but not applied for you. The older "Surface storage" field is really RECS, the
microrelief height that drives infiltration after the rain stops, and is now labelled that way.
The .sta file's "SurfStor ... mm" line is actually in metres (a label slip in the model). `python3 test_real_rfr.py` checks it against the formula.

## Guide support (values from the 1998 EUROSEM user guide)

All numbers were transcribed from the page images, not the garbled text
layer, and are in `eurosem_guide.py` with tests.

* **Soil texture** (11 USDA classes) fills conductivity, capillary drive,
  porosity and maximum moisture with the guide means (Tables A4.1, A4.2).
  **Surface state** (loose dry / sealed or compacted / loose moist) fills
  erodibility (A9.1) and cohesion (A9.2). Changing the surface state does not
  overwrite soil values you have edited.
* **Pickers** fill Manning n (A3.1), vegetation height, stem angle and shape
  (A7.2, height converted from metres to cm), interception (A7.1) and basal
  area (A7.3). A roots picker (A9.3) shows how much roots add to cohesion; it
  does not fill anything.
* **Warnings, never blocks:** values outside the guide range are flagged;
  maximum moisture above porosity, and initial moisture above maximum or below
  the residual moisture (THR), are flagged; a live line shows the model's
  cohesion factor.
* **Misprints in the guide**, kept as printed and flagged in the form:
  Table A4.1 silty clay loam THR 0.380 (neighbours 0.015 to 0.11); Table A9.2
  rows for silty clay, silty clay loam and sandy clay loam have a mean outside
  low to high. Where a fill would rely on such a row, the field is left alone
  and the form says why. Table A4.1's printed rows are shifted one line below
  "Sandy clay loam"; the mapping used was checked against porosity values.
* Guide units: the guide quotes porosity and moisture as % v/v; this form
  uses fractions (0.46). Table A7.2 heights are in metres; the model's input
  is cm.

## Sensitivity to parameter uncertainty

"Set up sensitivity" lists each parameter with Low, Base (the form's value)
and High. Low and High come from the guide range for each plane's chosen
texture, or you type your own. "Run sensitivity" runs the base case and then
each ticked parameter at its low and at its high value, one at a time, on every
plane that has a range (up to 40 runs). Results: a ranked chart of the percent
change from the base case for basin runoff, peak outflow or sediment leaving the
outlet, a table of the values and results, and a CSV. If a low porosity would
fall below a plane's maximum moisture, the moisture values are lowered to match
and the row says so. This is not a confidence interval and interactions
between parameters are ignored. "Reset" clears results and ranges.

What the first results showed (reference case, 100 mm/h for 30 min):
conductivity and capillary drive dominate runoff; porosity has no effect
(infiltration uses maximum minus initial moisture); erodibility has no visible
effect on sediment leaving the outlet, while cohesion above a few kPa cuts it
sharply. In the source, erodibility only scales splash detachment, and cohesion
scales the flow-detachment exchange rate through `0.79*exp(-0.85*COH)` (0.335
at or below 1 kPa). One event only: other conditions (short, weak flow, steep
early-event) may differ.

## Licence and credits

Results are model estimates and depend on the parameters entered. They are provided without warranty and should not be the sole basis for engineering, safety, regulatory or financial decisions. You are responsible for checking the inputs and for any use of the results. (The same notice is at the top of the page. The GPL's own warranty disclaimer and limitation of liability are sections 15 and 16 of the licence.)

The web interface, the Docker and hosting files, and the changes made to the model source for this
version are released under the **GNU General Public License v3** (see the `LICENSE` file in the
repository root; add it on GitHub with Add file > Create new file > `LICENSE` > Choose a license
template > GNU GPLv3). You may use, study, change and share the code; if you distribute a changed
version, you must do so under the same licence and make its source available.

Credits, to be kept with any copy:

* EUROSEM: Morgan, R.P.C., Quinton, J.N., Smith, R.E., Govers, G., Poesen, J.W.A., Auerswald, K.,
  Chisci, G., Torri, D., Styczen, M.E., Folly, A.J.V. (1998). *The European soil erosion model
  (EUROSEM): documentation and user guide.* Silsoe College, Cranfield University.
* Morgan, R.P.C., Quinton, J.N., Smith, R.E., Govers, G., Poesen, J.W.A., Auerswald, K., Chisci, G.,
  Torri, D. and Styczen, M.E. (1998). The European Soil Erosion Model (EUROSEM): a dynamic approach for
  predicting sediment transport from fields and small catchments. *Earth Surface Processes and
  Landforms*, 23, 527-544. (The journal paper; used here for the equations and for checking the
  guide and the source against each other.)
* The model's overland-flow and infiltration framework derives from KINEROS (USDA Agricultural
  Research Service).
* The guide tables built into the form (A3.1, A4.1, A4.2, A6.1, A7.1 to A7.3, A9.1 to A9.3) were
  transcribed from that user guide; misprints are flagged in the form.

## Discrepancies between the guide, the paper and the source

Found while building this version; each was checked against the Fortran source where it matters.
"Guide" is the 1998 user guide, "paper" is Morgan et al. (1998) in Earth Surface Processes and Landforms.

| Item | What the documents say | What the source does | Handling here |
|---|---|---|---|
| **Plant height units** | Guide Chapter 3, the worked example and Appendix 7 text: cm (`PLANTH`, 15.0 in the example file, echoed as 0.15 m). Guide Table A7.2: metres. Paper eq. 23: `PH` in metres. | Reads cm, divides by 100 (`emainl`, record 5), uses metres in `KINET`: `KE(LD) = 15.8 sqrt(PH) - 5.87`. Negative values are clamped to zero, which happens below about 14 cm. | Not a bug in the model; the documents mix the input unit (cm) with the internal unit (m). The form takes cm, the picker converts Table A7.2 (x 100), and a live note warns when a value of 2 or less (probably metres) is entered. |
| **Canopy cover units** | Guide Chapter 3 and Appendix 7: "percentage cover", while the same Chapter 3 entry goes on to say proportion 0 to 1. Paper: a fraction. | Used as a fraction (`COVER * RAIN`, `COVER * DICMAX`); no division by 100. | Form takes a fraction (0 to 1). |
| **Rill depth switch `RS`** | Guide example text describes it the opposite way round from Chapter 3. | `RS <= 0` on a plane with no upstream plane: depth scales with the square root of distance, so the entered depth is the bottom value. | Form follows the source and Chapter 3. |
| **`SurfStor` in the .sta file** | Labelled "mm". | Prints metres (`0.001 exp(-6.66 + 0.27 RFR)`). | Parser and notes treat it as metres; the real-model RFR test checks it. |
| **Table A4.1, silty clay loam THR** | Printed 0.380; neighbouring rows 0.015 to 0.11. | n/a | Kept as printed and flagged in the form. |
| **Table A9.2, three rows** | Silty clay, silty clay loam and sandy clay loam: mean outside the low to high range. | n/a | Kept as printed; the form does not fill cohesion from them and says why. |
| **Porosity, maximum and initial moisture** | Guide quotes porosity and moisture as % v/v. | Fractions. | Form takes fractions. |
| **"Surface storage"** | Presented as storage. | `RECS`, the microrelief height that drives infiltration after rain stops. Depression storage comes from `RFR` instead. | Field relabelled "Microrelief height (RECS)"; RFR is exposed. |
| **Erodibility sensitivity** | Guide Chapter 5 calls the model only moderately sensitive to erodibility. | In the reference case, erodibility 0 changed outlet sediment very little because flow detachment dominates; cohesion above a few kPa cut it sharply. | Recorded in the sensitivity notes; one storm and one soil only, so not a general conclusion. |
| **`FMINEW(I)` in the channel branch; stale COMMON blocks in `Euroa.for`** | n/a | Probable bugs; they do not affect planes. | Not fixed; channels are not in the form. |
| **Rill fill** | n/a | If deposition fills a rill the model stops ("excess rill fill") instead of continuing with a full rill. | The form explains the message in plain words. |

Anything else you find, send it with the line of source or the page of the guide or paper and it goes in this table.

## What is where

| File | Job |
|---|---|
| `app.html` | the form, help, pickers, sensitivity panel and results page |
| `server.py` | local web service (`/api/defaults`, `/api/guide`, `/api/run`, `/api/sensitivity`, `/api/sensitivity/plan`) |
| `eurosem_form.py` | form values to a model run; planes ordered upstream first |
| `eurosem_run.py` | writes the model's input files, runs it, reads the results |
| `eurosem_io.py` | .par / .pcp writers and readers, validation, .dyn / .sta parser |
| `eurosem_guide.py` | the guide's tables, with the misprints flagged; `guide_fill`, `guide_json` |
| `eurosem_sens.py` | one-at-a-time sensitivity and its CSV |
| `eurosem_export.py` | results and settings CSV files |
| `compare_runs.py` | compares a run with the frozen reference |
| `reference/` | the 40 cm, LOG10 reference run |
| `build_exe.sh` | compiles the model |

## Tests

Run each from this folder; each should end with `0 failure(s)`.

| Test | Needs |
|---|---|
| `test_eurosem_io.py` | none (see its header for the baseline folder argument) |
| `test_eurosem_form.py`, `test_eurosem_run.py` | none (the latter uses a stand-in model) |
| `test_eurosem_guide.py`, `test_eurosem_export.py`, `test_eurosem_sens.py` | none (`test_eurosem_export.py` skips one check if pandas is missing; that is fine) |
| `test_real_model.py` | gfortran and the compiled model: form defaults reproduce the reference; ROC checks |
| `test_real_sensitivity.py` | the compiled model: prints percent changes per parameter |
| `test_eurosem_server.py` | none (uses a stand-in model): public mode, run limits, health check |
| `test_real_rfr.py` | the compiled model: RFR 0 equals the reference; echoed depression storage matches the formula; prints effect of RFR |
| `test_real_rills.py` | the compiled model: rills off equals the reference; rills on completes, balance closes, rill erosion read; prints the comparison |

The Fortran regression harness is in `../headless` (see its README):
`./run_regression.sh stage4_rval` must print `SAME`.

## Rules the form enforces

* Planes only, in one cascade: each plane receives runoff from at most one
  other plane and feeds at most one; the last plane is the outlet.
* One slope per plane (the model uses the interrill slope for a plane
  without rills, so the form writes the slope to both).
* Rock fraction (ROC) and pavement fraction (PAVE) are per plane and are
  used as entered.
* Limits: 8 planes in the form (the model allows 60), 99 rainfall rows,
  duration / time step below 1000. Splash texture is fixed at 2.0 (the
  guide's constant). If the rain ends before the run does, no rain falls after
  the last row.

## Changes made to the model source

Commented sources are in the parent folder; compiled stages in `../headless`:
`FATAL`/`ENDOK` and a status file instead of `STOP`; a control file
`eurosem.ctl` replaces the keyboard questions; `MB` / `NI` set per element for
upstream inflow; ROC held per element. The LOG10 kinetic-energy correction in
`eurseml` is your own change and is kept. All stages reproduce the reference run
exactly.

## Known issues and not done yet

* Channels and ponds are not in the form. Rills are in and checked against the real model (see Rills).
* `FMINEW(I)` in the channel branch and stale COMMON blocks in `Euroa.for`
  are not fixed (they do not affect planes).
* Local use only: no sign-in, no limits on use, two model runs at a time.
  A public version needs a separate server for the Fortran (WordPress.com
  cannot run it) and run limits; the page can then be linked from WordPress.
* Sensitivity: one parameter at a time; no net-erosion-per-plane metric yet;
  maximum and initial moisture have no guide range (type your own).
