#!/bin/bash
# Several finder configurations in one process, same events and seeds, each
# given as NAME=cmd;cmd;... (shell commands applied before its events). A
# configuration's commands are applied on top of the previous one's, so each
# must set every knob it cares about. The first is the reference.
#
#   usage:  v2p2-eff-knobs.sh n_events sample out_prefix 'NAME=cmds' ...
#   Run it from a build directory; out_prefix is relative to it.
#   ASSOC=mtv, ITCONF=cmssw, FLAGS="..." as in v2p2-eff-v2-vs-v2p2.sh, and the
#   sample must carry the CMSSW tracks as there.
set -e
if [ ! -x ./mkFit ] || [ ! -f Makefile ]; then
  echo "v2p2-eff-knobs.sh: run it from a build directory, one with ./mkFit and the Makefile of configure." >&2
  exit 1
fi
if [ $# -lt 4 ]; then
  echo "usage: v2p2-eff-knobs.sh n_events sample out_prefix 'NAME=cmds' ..." >&2
  exit 1
fi
unset DISPLAY
export LD_LIBRARY_PATH=.
N=$1; SAMPLE=$2; OUT=$3; shift 3
if [ ! -f "$SAMPLE" ]; then
  echo "v2p2-eff-knobs.sh: sample '$SAMPLE' not found." >&2
  exit 1
fi
P=$(sed -n 's/^export SACMS *:= *//p' Makefile)/Probe
REF=${1%%=*}
CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --read-cmssw-tracks --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 ${FLAGS:-} --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $P/probe.C\")"
     --shell-command "probe_assoc_mtv($([ "${ASSOC:-}" = mtv ] && echo true || echo false))"
     --shell-command 'probe_eff_reset()'
     --shell-command "probe_eff_ref(\"$REF\")")
if [ "${ITCONF:-}" = cmssw ]; then
  CMD+=(--shell-command 's.SetMaxCands(6)'
        --shell-command 's.SetTrackScorer("phase1:default")'
        --shell-command 's.SetBkwSearchPickups(41, 7, 7, 7, 19)'
        --shell-command 's.SetBkwSearchV2p2(true)')
fi
for cfg in "$@"; do
  name=${cfg%%=*}; cmds=${cfg#*=}
  IFS=';' read -ra cl <<< "$cmds"
  for c in "${cl[@]}"; do CMD+=(--shell-command "$c"); done
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)" --shell-command 's.ProcessEventStd()'
          --shell-command "probe_eff_ev(s.event(), \"$name\")")
  done
done
CMD+=(--shell-command "probe_eff_report(\"$OUT\")")
echo .q | "${CMD[@]}"
