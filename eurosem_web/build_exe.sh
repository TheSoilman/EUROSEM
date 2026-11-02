#!/bin/bash
# Compile the model (stage 4) into ./eurosem next to this script.
HERE="$(cd "$(dirname "$0")" && pwd)"
STAGE="${1:-stage4_rval}"
gfortran -std=legacy -fno-automatic -fcheck=bounds -g -w -o "$HERE/eurosem" "$HERE/../headless/$STAGE"/*.for \
  && echo "Built $HERE/eurosem from $STAGE" || { echo "COMPILE FAILED"; exit 1; }
