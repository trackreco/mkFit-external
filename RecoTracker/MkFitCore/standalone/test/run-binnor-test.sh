#!/bin/bash
# Build and run the binnor correctness suite. No ROOT, no geometry, no TBB --
# it compiles in a couple of seconds and is meant to be run before and after any
# change to binnor.h, to a bin-range call site, or to the seed maker's use of
# either. Exit 0 means every REQUIRED check passed; documented defects are
# reported as "known" and do not fail the suite.
set -e
cd "$(dirname "$0")"
SRC=/foo/matevz/mic-dev/current/src
BIN=${TMPDIR:-/tmp}/binnor_test.$$
c++ -o "$BIN" -O2 -std=c++20 -I"$SRC" binnor_test.cxx "$SRC/RecoTracker/MkFitCore/src/radix_sort.cc"
"$BIN"; rc=$?
rm -f "$BIN"
exit $rc
