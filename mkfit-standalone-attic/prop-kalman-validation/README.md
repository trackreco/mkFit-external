# Synthetic propagation / Kalman / fit validation suite — atticked 2026-09-21

3651 lines, seven synthetic tasks plus a real-sim-track one. Removed from the
build (`MkFitCore/standalone/Makefile` takes `$(wildcard ${SADIR}/*.cc)`, so
moving the file out is sufficient — nothing else had to change).

## Why it went here rather than staying

It is not run live for algorithm tuning, and it was **27 % of the branch diff
against `treco/master`**. Keeping it in-tree also meant keeping it compiling
against five PRIVATE `MkFitCore/src/` headers — `Matrix.h`,
`PropagationMPlex.h`, `KalmanUtilsMPlex.h`, `MiniPropagators.h`, `MkBins.h` —
while three changes to exactly those signatures are in flight: the `pea`
removal, the mini-propagator `MPF` parametric-B overload, and the
covariance-growth estimator that replaces the propagated window covariance.
Atticking decouples that work from this suite instead of making it drag the
suite along.

Note the boundary: the other mkFit code in this repo (`mkfit-geom-cms-phase2`)
includes **interface-only**. This suite cannot, which is the second reason it
belongs in an attic rather than as live external code.

## What it established, so nobody re-runs it to find out

All of this is written up in `RecoTracker/CLAUDE.md`; the raw output is beside
this README.

- **Step 2** — `getS()`'s quadratic root was solved in the cancelling form.
  Fixing it moved the plane miss p90 from 54 um to 2.7e-6 cm. Superseded by
  Leonardo's PR #181, which does the same thing better and is now merged here.
- **Step 3** — the residual failure is the solver landing on the WRONG CROSSING
  of the module plane, not a convergence failure; `Config::nSStepsInProp2Plane`
  2 -> 8 changes the closure not at all.
- **The three-way solver A/B** — Leonardo's `getS` beats both the original and
  our Hermite on speed and at fit level. Do NOT merge the Hermite.
- **Finding 11** — the transported covariance is CURVILINEAR, blind along the
  momentum, with the surface term supplied downstream by `jacCurv2Loc`'s `cosz`.
- **Task 6** — pT resolution from a forward+backward fit of REAL sim tracks,
  0.4-1.6 % in the barrel, and it is what proved
  `/foo/matevz/mic-dev/ntuple_ttbar-p2.bin` geometry-incompatible (hit-to-its-
  own-module-plane distance median 0.26 cm against 2.9e-9 for the LST sample).

## THE ONE PART THAT IS STILL WANTED: Task 6

`pkv_task6_*`, ~1111 lines, fits sim tracks through their own rec hits. That is
the stated basis for **S13c** (forward-fit reference for the handover state).
Restoring it means copying it back and fixing it against whatever the propagator
API looks like then — budget for that in S13c rather than assuming it still
compiles.

Tasks 0-5 are synthetic and their conclusions are recorded; they are unlikely to
be needed again.

## Two traps if it is ever restored

- **`std::isfinite()` is unusable in this build.** `-Ofast` implies
  `-ffinite-math-only`, so a runtime 1/0 gives inf and `std::isfinite` returns
  **true**. Use `mkfit::isFinite` (`cms_common_macros.h`), which is bit-based.
  This silently corrupted the suite's own first pass.
- `Config::usePtMultScat` is false standalone and true in CMSSW, and phase-2
  sets neither. With it false, multiple scattering never reaches `1/pT`.
