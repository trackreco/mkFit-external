# Backward-fit hook + its analysis — retired 2026-09-21

Superseded by `TrBkFitUpdate` in the trace graph
(`MkFitCore/standalone/DataFormats/RntStructs.h`), written from
`MkFinder::bkFitFitTracksProp2Plane` and analysed by
`val_bkfit_trace_{reset,event,report}` in `RdfTrace/ValProp.cc`.

## Why it went

The hook was a `void (*)(const BkFitHitRec &)`, null by default, called per hit
in the backward fit, with the ROOT-writing collector on the standalone side. It
predated the trace port and was the only subsystem instrumented that way: the
search writes into the trace directly from `MkFinderV2p2.cc`, inside the same
nested TBB structure, with no hook. Keeping one for the fit alone was an
inconsistency with nothing behind it.

## What replaced it, and the evidence it is equivalent

Both paths were run on the SAME 20 events in one process — hook armed and trace
accumulating together, **49765 records each** — and every data line of
`an-val-bkfit.C` matches the trace report, checked with `diff`, not by eye:

    per track   pure 0.659   gf>=0.8 0.934   gf<0.8 8.91
    per hit     pure 0.867   gf>=0.8 1.15    gf<0.8 4.61
    truth split pure 23980 matched / 0 other      (0 is correct: gf==1.0)
                gf>=0.8 10379 @ 1.08 / 1421 @ 2.07
                gf<0.8   4219 @ 2.89 / 2891 @ 10

A separate element-wise comparator (`val_bkfit_compare_trace`, also retired)
reported 0 mismatches in layer, step, fail, chi2 and all three module-frame
residuals, plus 0 chain breaks and 0 wrong-stage states.

## Two defects the port exposed, both since fixed in the live code

- **`mc_match` as a bool cannot be computed in `MkFinder`.** It needs the
  track's SIM label, and the only label available there is the seed's, which
  after `relabelSeedTracksSequentially()` is the seed's INDEX — a different
  namespace. `TrBkFitUpdate` stores `mc_track_id` and resolves the comparison in
  analysis via `TrCandMeta::global_seed`.
- **`fail` is 0 on every record**, because Leonardo's `getS` carries no bracket
  and so detects no failures. Not a bug, and not evidence the fit did not fail.

## Running the old macro, if it is ever wanted

It needs the **dev** ROOT — `root.exe` on `PATH` is the system build, and
loading a dictionary built against dev ROOT into it segfaults in
`TGenericClassInfo::~TGenericClassInfo` at exit, after appearing to run and
writing no output file:

    cd <standalone>
    export LD_LIBRARY_PATH=/baz/matevz/root-dev/dev-1-bld/lib:.
    export ROOT_INCLUDE_PATH=<src root>
    /baz/matevz/root-dev/dev-1-bld/bin/root.exe -l -b -q \
      '../RecoTracker/MkFitCore/standalone/test/an-val-bkfit.C("in.root","prefix")'

It also needs `ValFitHit` (removed from `ValStructs.h`) and its `LinkDef`
entries restored.
