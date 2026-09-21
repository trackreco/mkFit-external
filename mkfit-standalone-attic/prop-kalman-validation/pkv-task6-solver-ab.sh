#!/bin/bash
# Task 6, three solver configurations over the SAME tracks in ONE process.
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.

SAMPLE=/foo/matevz/mic-dev/TTbar_noPU_100events_HLT_Phase2_singleTrackingIteration_CAExtSeedsUpToThreeOTBarrelLayer_into_LSTSeedsIncludingT5s.bin
N=${1:-20}
MAXHPL=${2:-3}
OUT=${3:-pkv-task6-solver-ab}
OUTDIR=../RecoTracker/MkFitCore/standalone/test

CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L ../RecoTracker/MkFitCore/standalone/test/prop-kalman-validation.C\")")

add_config () {   # $1 = hermite, $2 = stable_root, $3 = label
  CMD+=(--shell-command "pkv_solver($1,$2)")
  CMD+=(--shell-command 'pkv_task6_reset()')
  CMD+=(--shell-command "pkv_task6_maxhpl($MAXHPL)")
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)" --shell-command 'pkv_task6_ev(s.event())')
  done
  CMD+=(--shell-command 'pkv_task6_report()')
  CMD+=(--shell-command "pkv_task6_stash(\"$3\")")
}

add_config 0 0 A_orig_getS
add_config 0 1 B_getS_step2
add_config 1 1 C_hermite_step3

CMD+=(--shell-command "pkv_task6_root(\"${OUTDIR}/${OUT}.root\")")
CMD+=(--shell-command "pkv_task6_badcases(\"../RecoTracker/MkFitCore/standalone/test/${OUT}-badcases.txt\", \"$SAMPLE\")")

echo .q | "${CMD[@]}" 2>&1 | tee "${OUTDIR}/${OUT}.txt"
