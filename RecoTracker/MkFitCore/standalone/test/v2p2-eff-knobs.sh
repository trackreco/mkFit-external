#!/bin/bash
# Several finder configurations in one process, same events and seeds, each
# given as NAME=cmd;cmd;... (shell commands applied before its events). A
# configuration's commands are applied on top of the previous one's, so each
# must set every knob it cares about. The first is the reference.
#
#   usage:  v2p2-eff-knobs.sh n_events sample out_prefix 'NAME=cmds' ...
#   ASSOC=mtv, ITCONF=cmssw, FLAGS="..." as in v2p2-eff-v2-vs-v2p2.sh.
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=$1; SAMPLE=$2; OUT=$3; shift 3
T=../RecoTracker/MkFitCore/standalone/test
REF=${1%%=*}
CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --read-cmssw-tracks --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 ${FLAGS:-} --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command "val_assoc_mtv($([ "${ASSOC:-}" = mtv ] && echo true || echo false))"
     --shell-command 'val_eff_reset()'
     --shell-command "val_eff_ref(\"$REF\")")
if [ "${ITCONF:-}" = cmssw ]; then
  CMD+=(--shell-command 'val_max_cands(6)'
        --shell-command 'val_track_scorer("phase1:default")'
        --shell-command 'val_bkw_pickups(5, 15, 7, 15, 5)'
        --shell-command 's.SetBkwSearchV2p2(true)')
fi
for cfg in "$@"; do
  name=${cfg%%=*}; cmds=${cfg#*=}
  IFS=';' read -ra cl <<< "$cmds"
  for c in "${cl[@]}"; do CMD+=(--shell-command "$c"); done
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)" --shell-command 's.ProcessEventStd()'
          --shell-command "val_eff_ev(s.event(), \"$name\")")
  done
done
CMD+=(--shell-command "val_eff_report(\"$OUT\")")
echo .q | "${CMD[@]}"
