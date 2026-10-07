#!/bin/bash
# Single-threaded build timing: V1/CloneEngine against v2p2 as a single greedy
# path (cap 1, max_hits_per_layer 1) and at two caps.
#
# TWO NUMBERS PER CONFIGURATION, and they mean different things:
#   "Iteration 0 build time"  -- runBtpCe_MultiIter brackets ONLY
#                               (builder.*FindTracks)(IT_FwdSearch), i.e. the
#                               forward search. Not seed loading, not the
#                               backward fit, not duplicate removal. This is the
#                               inner path.
#   "Total event loop time"   -- read_in + loadHitsAndBeamSpot + the build +
#                               export, wall clock over the whole loop. The file
#                               is ~9.7 GB, so this is partly I/O and only
#                               comparable between configurations run back to
#                               back on a warm page cache.
# Both are printed per run; the report takes the MIN over repetitions.
#
# THE BINARY MUST BE BUILT WITH EVENT_RDF_TRACE OFF. With it on, Makefile.config
# also sets TBB_DEBUG, which turns TBB_PARALLEL_FOR into a serial loop, and the
# search writes a trace record per scanned hit. Timing that measures the
# instrument. The script refuses to run against a trace build.
#
# --num-iters-cmssw 1 is REQUIRED: without it runBtpCe_MultiIter walks
# Config::ItrInfo past what the phase-2 plugin filled and segfaults.
#
#   usage:  v2p2-timing.sh [n_events] [n_reps] [sample]
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-100}
REPS=${2:-3}
SAMPLE=${3:-/foo/matevz/mic-dev/ttbar-PU200-D121-C22-100ev.bin}
OUT=${TMPDIR:-/tmp}/mkfit-timing-$$

if ./mkFit --help 2>&1 | head -40 | grep -q 'shell'; then :; fi
if ldd ./mkFit | grep -qi libCore; then
  echo "REFUSING: mkFit is linked against ROOT, so this is very likely a trace build." >&2
  echo "Rebuild with EVENT_RDF_TRACE commented out in Makefile.config." >&2
  exit 1
fi

BASE=(--geom CMS-phase2 --seed-input cmssw --input-file "$SAMPLE"
      --num-events "$N" --num-thr 1 --num-thr-ev 1 --num-iters-cmssw 1 --silent)

#     label            extra flags
CFG=("V1_ce_c6|--build-mimi --max-cands-per-seed 6"
     "greedy_c1|--build-mimi --build-mimi-v2p2 --max-cands-per-seed 1 --v2p2-max-hits-per-layer 1"
     "comb_c3|--build-mimi --build-mimi-v2p2 --max-cands-per-seed 3"
     "comb_c4|--build-mimi --build-mimi-v2p2 --max-cands-per-seed 4")

mkdir -p "$OUT"
echo "# mkFit single-thread timing, $N events, $REPS repetitions, $(date -Is)"
echo "# sample: $SAMPLE"
echo "# host:   $(hostname), $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2-)"
echo

for c in "${CFG[@]}"; do
  label=${c%%|*}
  flags=${c##*|}
  for ((r=1;r<=REPS;r++)); do
    ./mkFit "${BASE[@]}" $flags > "$OUT/$label.$r.log" 2>&1
    inner=$(grep -m1 "Iteration 0 build time =" "$OUT/$label.$r.log" | awk '{print $6}')
    total=$(grep -m1 "Total event loop time" "$OUT/$label.$r.log" | awk '{print $5}')
    echo "$label rep$r inner=$inner total=$total"
  done
done

echo
echo "# MIN over $REPS reps, seconds for $N events"
printf "%-14s %12s %12s %12s %12s\n" config inner_s inner_ms_ev total_s total_ms_ev
for c in "${CFG[@]}"; do
  label=${c%%|*}
  grep -h "Iteration 0 build time =" "$OUT/$label".*.log | awk '{print $6}' | sort -g | head -1 > "$OUT/$label.inner"
  grep -h "Total event loop time"    "$OUT/$label".*.log | awk '{print $5}' | sort -g | head -1 > "$OUT/$label.total"
  i=$(cat "$OUT/$label.inner"); t=$(cat "$OUT/$label.total")
  printf "%-14s %12.3f %12.1f %12.3f %12.1f\n" "$label" "$i" "$(echo "$i*1000/$N" | bc -l)" "$t" "$(echo "$t*1000/$N" | bc -l)"
done
echo
echo "# logs in $OUT"
