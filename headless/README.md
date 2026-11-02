# EUROSEM headless clean-up: staged versions

Each folder is a complete set of sources (6 commented files + `ehless.for`).
Compile all the `.for` files in the folder together. Each stage includes the
previous one, so if a stage fails you know which change to look at.

| Stage | Change | Results should be |
|---|---|---|
| `stage1_fatal` | Every `STOP` becomes `CALL FATAL(msg)` (message to screen, aux file and `eurosem.status`; exit code 2). Normal end calls `ENDOK` (status `OK`). Still asks its questions. | identical to reference |
| `stage2_headless` | If `eurosem.ctl` exists, all keyboard questions are skipped and file names (up to 80 characters) come from it. Without it, behaves as before. | identical to reference |
| `stage3_mbfix` | `PLANE` sets `MB` and `NI` itself for an element with upstream inflow (previously left over from the previous element). | identical to reference |
| `stage4_rval` | ROC (rock content) now works per element. Before, one value was shared and reset for every element, so only the last element's ROC counted. | identical to reference (ROC = 0 there) |

## Testing a stage

    cd headless
    ./run_regression.sh stage1_fatal
    ./run_regression.sh stage2_headless
    ./run_regression.sh stage3_mbfix
    ./run_regression.sh stage4_rval

Each run compiles, runs the reference case (`eurosem_web/reference/`, your
40 cm LOG10 run), and prints `SAME` or `DIFFERENT`. If `DIFFERENT`, the
first values that moved are listed. If compilation fails, send me the
first few error lines.

## eurosem.ctl (stage 2 onwards)

Seven lines, in the folder the program runs in:

    run title (up to 80 characters)
    parameter file name
    rainfall file name
    dynamic output file name
    auxiliary output file name
    static output file name
    Y        <- answer to the ROC / PAVE / SPLTEX questions (Y = accept)

After a run, `eurosem.status` contains `OK` or `ERROR: <message>`.

## What was changed, exactly

`python3 codediff.py <folder A> <folder B>` lists every code line that
differs between two sets of sources (comments ignored). Totals: stage 1
changes 57 code lines (all `STOP` statements) plus the new file; stage 2
adds the `eurosem.ctl` handling (83 lines); stage 3 adds 3 lines.
Rebuild all three stages from the commented sources with
`python3 build_headless.py`.

## Not changed

Channel and pond code (apart from the `STOP` statements), the
`FMINEW(I)` index in the channel branch of `IMPLCT`, and the stale
COMMON blocks in `Euroa.for`. A plane-only run does not use them.
