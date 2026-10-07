#!/bin/bash
# V2 against v2p2 in one process, same events and seeds: the standalone
# counterpart of the CMSSW MTV comparison (standard V2 vs the v2p2 branch).
# Production defaults otherwise; the switch is Shell::SetUseV2p2(). Also
# reports val_hitorder: whether the hits the CMSSW refit would get, after
# MkFitOutputConverter's sort, are in propagation order.
#
#   usage:  v2p2-eff-v2-vs-v2p2.sh [n_events] [sample] [out_prefix]
#   ASSOC=mtv associates by MTV's rule (> 75 % of all hits true) instead of
#   quality-val's (>= 50 % of the non-seed hits).
#   ITCONF=cmssw sets the things in which the CMSSW release's
#   mkfit-phase2-initialStep.json (data-RecoTracker-MkFit V00-20-00) differs
#   from the plugin: maxCandsPerSeed 6 (plugin 3), the phase1:default track
#   scorer (the plugin's default too since 2026-10-04; it was phase2:LstIntoPix
#   before), backward-search pickups at plan index
#   5 15 7 15 5 (plugin 7 27 10 27 7). Found by diffing --json-save-iterations
#   against the release file. It also runs the backward search with v2p2 for
#   the v2p2 configuration (Shell::SetBkwSearchV2p2), as run_OneIteration()
#   does in CMSSW; the shell's default is the clone engine. (--json-load parses a file and DROPS it, and
#   --json-patch does not know these keys, hence the shell setters.)
#   FLAGS="..." extra mkFit options, e.g. "--backward-fit --remove-dup", which
#   CMSSW runs (MkFitProducer defaults) and a bare driver does not.
#   out_prefix is relative to the standalone build dir; pass an absolute path.
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-30}
SAMPLE=${2:-/foo/matevz/mic-dev/ttbar-PU200-D121-C22-100ev-rt.bin}
OUT=${3:-eff-v2-vs-v2p2}
T=../RecoTracker/MkFitCore/standalone/test

CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --read-cmssw-tracks --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2
     ${FLAGS:-} --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command "val_assoc_mtv($([ "${ASSOC:-}" = mtv ] && echo true || echo false))"
     --shell-command 'val_eff_reset()'
     --shell-command 'val_eff_ref("v2")')

if [ "${ITCONF:-}" = cmssw ]; then
  CMD+=(--shell-command 'val_max_cands(6)'
        --shell-command 'val_track_scorer("phase1:default")'
        --shell-command 'val_bkw_pickups(5, 15, 7, 15, 5)'
        --shell-command 's.SetBkwSearchV2p2(true)')
fi

for c in v2 v2p2; do
  CMD+=(--shell-command "s.SetUseV2p2($([ $c = v2p2 ] && echo true || echo false))")
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)"
          --shell-command 's.ProcessEventStd()'
          --shell-command "val_eff_ev(s.event(), \"$c\")"
          --shell-command "val_ho_ev(s.event(), \"$c\")")
  done
done
CMD+=(--shell-command "val_eff_report(\"$OUT\")" --shell-command 'val_ho_report()')
echo .q | "${CMD[@]}"
