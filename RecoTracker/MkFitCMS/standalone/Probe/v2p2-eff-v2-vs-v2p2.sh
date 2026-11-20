#!/bin/bash
# V2 against v2p2 in one process, same events and seeds: the standalone
# counterpart of the CMSSW MTV comparison (standard V2 vs the v2p2 branch).
# Production defaults otherwise; the switch is Shell::SetUseV2p2().
#
#   usage:  v2p2-eff-v2-vs-v2p2.sh <n_events> <sample> [out_prefix]
#   Run it from a build directory. The sample must carry the CMSSW tracks (a .bin
#   written with --write-rec-tracks), which --read-cmssw-tracks requires.
#   ASSOC=mtv associates by MTV's rule (> 75 % of all hits true) instead of
#   quality-val's (>= 50 % of the non-seed hits).
#   ITCONF=cmssw sets the things in which the CMSSW release's
#   mkfit-phase2-initialStep.json (data-RecoTracker-MkFit V00-20-00) differs
#   from the plugin: maxCandsPerSeed 6 (the plugin's too since 2026-10-27; it
#   was 3 before), the phase1:default track scorer (the plugin's default too
#   since 2026-10-04; it was phase2:LstIntoPix before), backward-search pickups
#   at plan index
#   5 15 7 15 5 (plugin 7 27 10 27 7). Found by diffing --json-save-iterations
#   against the release file. The script names them by layer, 41 7 7 7 19, which
#   are the layers at those indices of the plugin's plans (its own pickups are
#   at layers 43 53 10 31 21). It also runs the backward search with v2p2 for
#   the v2p2 configuration (Shell::SetBkwSearchV2p2), as run_OneIteration()
#   does in CMSSW; the shell's default is the clone engine. (--json-load parses a file and DROPS it, and
#   --json-patch does not know these keys, hence the shell setters.)
#   FLAGS="..." extra mkFit options, e.g. "--backward-fit --remove-dup", which
#   CMSSW runs (MkFitProducer defaults) and a bare driver does not.
#   out_prefix is relative to the build directory.
set -e
if [ ! -x ./mkFit ] || [ ! -f Makefile ]; then
  echo "v2p2-eff-v2-vs-v2p2.sh: run it from a build directory, one with ./mkFit and the Makefile of configure." >&2
  exit 1
fi
if [ $# -lt 2 ]; then
  echo "usage: v2p2-eff-v2-vs-v2p2.sh <n_events> <sample> [out_prefix]" >&2
  exit 1
fi
unset DISPLAY
export LD_LIBRARY_PATH=.
N=$1
SAMPLE=$2
OUT=${3:-eff-v2-vs-v2p2}
if [ ! -f "$SAMPLE" ]; then
  echo "v2p2-eff-v2-vs-v2p2.sh: sample '$SAMPLE' not found." >&2
  exit 1
fi
P=$(sed -n 's/^export SACMS *:= *//p' Makefile)/Probe

CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --read-cmssw-tracks --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2
     ${FLAGS:-} --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $P/probe.C\")"
     --shell-command "probe_assoc_mtv($([ "${ASSOC:-}" = mtv ] && echo true || echo false))"
     --shell-command 'probe_eff_reset()'
     --shell-command 'probe_eff_ref("v2")')

if [ "${ITCONF:-}" = cmssw ]; then
  CMD+=(--shell-command 's.SetMaxCands(6)'
        --shell-command 's.SetTrackScorer("phase1:default")'
        --shell-command 's.SetBkwSearchPickups(41, 7, 7, 7, 19)'
        --shell-command 's.SetBkwSearchV2p2(true)')
fi

for c in v2 v2p2; do
  CMD+=(--shell-command "s.SetUseV2p2($([ $c = v2p2 ] && echo true || echo false))")
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)"
          --shell-command 's.ProcessEventStd()'
          --shell-command "probe_eff_ev(s.event(), \"$c\")")
  done
done
CMD+=(--shell-command "probe_eff_report(\"$OUT\")")
echo .q | "${CMD[@]}"
