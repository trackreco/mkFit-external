#!/bin/bash
# Driver for the pure-geometry seeding study (standalone/SeedGeom.cc).
#
#   <mkFit-external checkout>/mkfit-standalone-seedgeom/seedgeom-run.sh <n_events> <name> [sg_* cmds...]
#
# Run from the standalone build directory.  Output lands in
# standalone/test-seedgeom/<name>.{txt,root} unless <name> is a path.
#
# e.g.  ./seedgeom-run.sh 20 sg-l012 'sg_ptmin(0.9)' 'sg_d0max(0.1)' 'sg_layers(0,1,2)'
#
# Run from the standalone build directory.
set -e
N=${1:-10}
PREFIX=${2:-seedgeom}
case "$PREFIX" in */*) ;; *) mkdir -p test-seedgeom; PREFIX=test-seedgeom/$PREFIX ;; esac
shift 2 || true

# overridable: SAMPLE=... GEOM=... XOPTS='--read-sim-hit-states' seedgeom-run.sh ...
SAMPLE=${SAMPLE:-/foo/matevz/mic-dev/trackingNtuple_HLT_2026_March.bin}
GEOM=${GEOM:-CMS-phase2}
# the loader next to THIS script, so the study is compiled from this checkout
LOADER=$(dirname "$(realpath "$0")")/seedgeom-load.C

unset DISPLAY
export LD_LIBRARY_PATH=.

CMD=(./mkFit --geom "$GEOM" --seed-input cmssw --input-file "$SAMPLE" $XOPTS
     --num-events "$N" --num-thr 1 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".x $LOADER\")")
for c in "$@"; do CMD+=(--shell-command "$c"); done
CMD+=(--shell-command 'sg_reset()')
for ((i=1;i<=N;i++)); do
  CMD+=(--shell-command "s.GoToEvent($i)" --shell-command 'sg_ev(s.event())')
done
CMD+=(--shell-command 'sg_report()' --shell-command "sg_write(\"$PREFIX\")")

echo .q | "${CMD[@]}" 2>&1 | tee "$PREFIX.txt"
