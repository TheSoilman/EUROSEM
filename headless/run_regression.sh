#!/bin/bash
# Compile one stage, run the reference case, compare with the reference.
#   ./run_regression.sh stage1_fatal
#   ./run_regression.sh stage2_headless
#   ./run_regression.sh stage3_mbfix
# Needs gfortran and python3.  Works in a temporary folder.
STAGE="$1"
[ -d "$(dirname "$0")/$STAGE" ] || { echo "usage: $0 stage1_fatal|stage2_headless|stage3_mbfix"; exit 1; }
HERE="$(cd "$(dirname "$0")" && pwd)"
REF="$HERE/../eurosem_web/reference"
W="$(mktemp -d)"
cp "$REF/testn.par" "$REF/test.pcp" "$W/"
echo "Compiling $STAGE ..."
gfortran -std=legacy -fno-automatic -fcheck=bounds -g -w -o "$W/eurosem" "$HERE/$STAGE"/*.for || { echo "COMPILE FAILED"; exit 1; }
cd "$W"
if [ "$STAGE" = stage1_fatal ]; then
  # stage 1 is still interactive: answer the questions on standard input
  printf '\ntest\nn\ntestn.par\ntest.pcp\nnew.sta\nnew.dyn\nnew.aux\nn\n' | ./eurosem > run.log 2>&1
else
  # stages 2 and 3: no keyboard, names come from eurosem.ctl
  printf 'test\ntestn.par\ntest.pcp\nnew.dyn\nnew.aux\nnew.sta\nY\n' > eurosem.ctl
  ./eurosem < /dev/null > run.log 2>&1
fi
echo "exit code: $?   (0 = normal, 2 = stopped by FATAL)"
echo "status file: $(cat eurosem.status 2>/dev/null || echo 'none')"
echo "--- comparison with reference ---"
python3 "$HERE/../eurosem_web/compare_runs.py" new.dyn new.sta
echo "Files are in $W (run.log has the screen output)"
