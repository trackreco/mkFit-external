#!/bin/bash
# The chosen seedsurf configuration (2026-09-27): the feed-forward chain,
# starting up to two crossed layers late, no other holes, cleaning N = 3,
# with the D121 window tables in windows-D121/.
#
#   seedsurf-chain.sh [seedsurf options ...]
#
# Environment: B (the build directory, default the isolated seeding build),
# S (the sample), GEOM, BIND (truth binding in cm; needs SimHitStates in the
# sample, empty to turn it off). Anything on the command line is appended, e.g.
#
#   seedsurf-chain.sh --first-event 40 --num-events 60 --truth truth.txt
#
# gives the chain row of the README (1431.5 found tracks / ev, 27.5k quads, fake 0.343).

set -e
D=$(cd "$(dirname "$0")" && pwd)
B=${B:-/foo/matevz/mic-dev/CMSSW_14_1_0_pre0-p2p/src-seeding/standalone}
S=${S:-/foo/matevz/mic-dev/ttbar-PU200-D121-C22-100ev.bin}
GEOM=${GEOM:-CMS-phase2-Run4D121}
BIND=${BIND-0.05}

WIN=$(cat "$D"/windows-D121/pixel.txt "$D"/windows-D121/ot1p.txt "$D"/windows-D121/skip.txt)

cd "$B"
# $WIN is deliberately unquoted: one option word per token
LD_LIBRARY_PATH=. exec test-seedgeom/bin/seedsurf --input-file "$S" --geom "$GEOM" \
  --pt-min 0.9 --marg-b 0.001 ${BIND:+--bind $BIND} \
  $WIN \
  --chain 0 --chain-start-holes 2 --chain-lead-only --dedup 3 \
  "$@"
