# Pure-geometry pixel seeding — first measurement, 2026-09-21

Sample `trackingNtuple_HLT_2026_March.bin`, phase-2, PU200.

**Nothing of this lives in the CMSSW source tree.** Code is here in
`mkFit-external/mkfit-standalone-seedgeom/`; it is compiled on demand by
ACLiC, which also generates the dictionary that makes the entry points
callable from cling. Artefacts go to the build dir, `standalone/test-seedgeom/`.

| file | what |
|---|---|
| `SeedGeom.{h,cc}` | the study; short global-scope `sg_*` wrappers at the end of the .cc |
| `seedgeom-load.C` | ACLiC loader — include paths, defines, build dir, lib symlinks |
| `seedgeom-run.sh` | driver; generates the per-event `--shell-command` loop |
| `README.md` | this file: the durable record. The raw logs are regenerable |

## The flags are NOT hardcoded

`seedgeom-load.C` gets them from **`make echo-aclic`**, a target added to
`MkFitCore/standalone/Makefile` (and to the build-dir Makefile that `configure`
writes -- `configure` is bash, not perl). A stale copy in the loader would be
silent corruption rather than a build error: `standalone/Event.h` gains members
under `MKFIT_TRACE`, so an ACLiC build without the identical `-D` set does not
fail to link, it reads every member at the wrong offset. If the target cannot
be reached the loader **refuses to load** rather than guessing.

    ACLIC_DEFS=-DMPLEX_USE_INTRINSICS -DMKFIT_STANDALONE -DMKFIT_TRACE ...
    ACLIC_OPT=-g0 -Ofast -mavx -fPIC -ftree-vectorize ... -mrecip=none

Two things the target deliberately does NOT emit, and one thing the loader
must do:

- **No ROOT `-I` or `-std`.** ACLiC supplies those from the RUNNING ROOT, which
  is the ROOT the library was built against. Letting the Makefile supply them
  would mean compiling against whichever ROOT the sub-shell's `root-config`
  finds.
- **The loader runs make under `gROOT->GetRootSys()`.** Without it the
  sub-shell picks up system ROOT: measured, `-I/usr/include/root` and
  `-std=c++17` instead of `-I/baz/matevz/root-dev/dev-1-bld/include` and
  `c++20`. The `-D` set survives that (it comes from `EVENT_RDF_TRACE`, not
  from root-config) but `CXX_STD` does not, so the flags would describe a
  different build from the one being linked against.
- **`ACLIC_OPT` is not cosmetic.** `-Ofast` implies fast-math; without it the
  results differ in the last digit and at the 1e-7 level in counts
  (`ref-aclic.txt`) and timings are not comparable. Loaded with `.L ...cc+O`.

**Standing check:** `ref-lib.txt` (the same file compiled into libMkFitCore)
vs `ref-mk2.txt` (ACLiC, flags from the Makefile) are **bit-identical**. Re-run
that pair if the build configuration changes.

**Question.** Can a good starting sample be extracted from layer-to-layer
deltas alone — no fit, no covariance, no vertex as an input — with the physics
entering only as two scalars that widen one integer window?

    phi   w(a->b) = (r_b - r_a)/(2 R_min)  +  D0_max (1/r_a - 1/r_b)
    q     z at the target predicted from the doublet, linear in r

Three space points are 6 measurements against 5 helix parameters, so there is
exactly ONE constraint and it is longitudinal. It involves neither the beamline
nor D0, so displaced tracks pass it unchanged.

## A. Truth vertices (20 events)

The `.bin` has no vertex collection. Sim-track parameters are at the production
point, so grouping recovers them — but NOT by an exact bit-pattern group-by:
alongside each PU vertex sit a few hundred low-multiplicity points at the same
z (decays inside the beam width), so a **z cluster** of on-beamline points is
the vertex. `ProdType::Signal` labels every particle of the signal event
INCLUDING its secondaries out at r = 130 cm, so it does not by itself identify
the vertex.

| | |
|---|---|
| sim tracks | 45924 /ev |
| **primary vertices** (dr < 200 um, z-gap 50 um) | **211.7 /ev** |
| ... with >= 2 tracks of pT > 0.9 | 125.3 /ev |
| vertex z spread (robust sigma) | 4.39 cm (beamspot sigmaZ 4.26) |
| vertex multiplicity, median | 46 tracks — but only **3** above pT 0.9 |
| nearest-neighbour dz | p10 116 um, median **506 um**, p90 2.27 mm |
| off-beamline production points | 15767 /ev, median r 27 cm — a free displaced sample |

## B. Geometric quadruplets, pixel barrel layers 0,1,2 + 3

pT_min 0.9, D0_max 0.1 cm, qwin 0.035 / 0.025 cm, phi window at layer 4
0.002 rad. 20 events. (`qA.txt`)

| stage | per event | per previous stage |
|---|---|---|
| hits, layers 0/1/2/3 | 9123 / 7326 / 5836 / — | |
| doublets (phi window only) | **837 323** | 91.8 per layer-0 hit |
| layer-c hits touched | 1.89 M | 2.25 per doublet |
| triplets after the q prediction | 43 815 | **28.9x rejection** |
| layer-d hits touched | 14 861 | **0.34 per triplet** |
| **quadruplets** | **812** | **58.3x rejection** |

- **purity 73.7 %**, efficiency **92.9 %** against the findable denominator
  (sim tracks with hits in all four layers AND pT > pT_min AND production
  radius < D0_max — quoting against every sim track with four hits would be
  quoting against a ceiling the window cannot reach).
- efficiency vs pT: 0.891 / 0.938 / 0.969 / 0.976 / 0.982 / 0.976 over
  0.9-1.2 / 1.2-1.6 / 1.6-2.2 / 2.2-3.2 / 3.2-5 / >5 GeV. The low-pT loss is
  consistent with the q window being ~1.5 sigma there — the windows should be
  pT-dependent, they are not yet.
- **z_v (beamline intercept) per quadruplet: 76 um**, and it scales as 1/p
  exactly as multiple scattering says: 135 / 103 / 70 / 43 / 28 um over
  pT < 0.7 / 0.7-1 / 1-2 / 2-5 / > 5 GeV.
- three-point curvature: (pT_fit - pT_true)/pT_true robust sigma **4.6 %**.
- residuals at the 4th layer, TRUE quads: **dphi 0.63 mrad, dz 83 um**.
- **0.30 s/event marginal**, single thread, scalar, including the event read.

The 4th layer costs 0.34 hits touched per triplet and buys 58x — the helix is
determined by then, so it is a lookup, not a search. Only the first pair is
O(N^2).

## C. Displacement costs exactly what the formula says

`s1.txt` vs `s2.txt` (10 ev): D0_max 0.1 -> 1.0 cm on layers 0,1,2+3.

| | D0 0.1 cm | D0 1.0 cm |
|---|---|---|
| doublets /ev | 801 769 | 3 938 795 (**4.7x**) |
| findable tracks /ev | 363 | 371 (**+8**) |
| quad purity | 0.751 | 0.412 |

Predicted 4.9x from `(0.0198 + 0.166 + 0.002)/(0.0198 + 0.0166 + 0.002)`.
Measured 4.7x. And it buys 8 tracks/ev — at pixel radii there is very little
displaced production to find, so the 1 cm allowance is poor value HERE.

The D0 term goes as `1/r_inner`, the curvature term as `Dr`, so **starting
further out is much cheaper for displaced tracks**: layers 1,2,3 at D0 1.0 cm
give 1 268 808 doublets against 3 938 795 (`s3.txt`) — **3.1x fewer**.

## D. Vertices from the votes alone — works for the busy ones, not the median

Greedy peak extraction on the quadruplet z_v values, no tracking anywhere.
Match radius 500 um.

| pT_min | doublets /ev | quads /ev | quad purity | found /ev | **eff, all 211** | **eff, the 125 with >= 2 tk** | fake frac | vertex sigma_z |
|---|---|---|---|---|---|---|---|---|
| 0.9 | 837 k | 812 | 0.737 | 69 | **0.320** | **0.492** | 0.024 | **54 um** |
| 0.5 | 1129 k | 2220 | 0.495 | 204 | 0.631 | 0.764 | 0.350 | 96 um |
| 0.3 | 1619 k | 5681 | 0.264 | 631 | **0.866** | **0.902** | 0.712 | 137 um |

Efficiency vs the vertex's own n(pT > 0.9), at pT_min 0.9 / >= 3 votes:
0.056 / 0.094 / 0.118 / 0.234 / **0.708** for n = 0 / 1 / 2-3 / 4-7 / >= 8.

**This is the population statement, not a machinery failure.** The median PU
vertex has 3 tracks above 0.9 GeV, so at pT_min 0.9 it can produce at most a
handful of votes. The hard-scatter and the busy vertices are found at 71-88 %.
Reaching all 211 needs pT_min ~ 0.3, which costs only **1.9x** in doublets
(the window is partly floored by the D0 and margin terms) but drops quad purity
to 26 %, and the fake-vertex rate then reflects that. The fake numbers are from
a deliberately crude finder — no vote weighting, no quality cut — so treat them
as an upper bound, not a limit.

## E. Extending outward into TBPS-P (layer 4) — works, needs pT-dependent windows

`qC.txt`, layers 1,2,3 + 4, 10 events. Efficiency **74.3 %**, purity 41.9 %,
vertex sigma_z 97 um, doublets 475 620/ev (43 % fewer than 0,1,2+3, since layer
0 is the busiest). The transverse prediction over the 10 cm gap is still
excellent — **dphi 1.03 mrad** — but the true dz spread is **0.62 mm** (MS),
against a 1.2 mm window, i.e. 2 sigma, and that is where the efficiency goes:
0.627 at pT 0.9-1.2 rising to 0.983 above 5 GeV. Not a physics limit, a
window-shape limit.

## Known limitations of this measurement

- Windows are **flat in pT**; every efficiency loss above is at the low-pT end
  and is consistent with the window being 1.5-2 sigma there.
- The layer's ~1 cm radial thickness forces the q bin-range to be widened by
  `|cot theta| * dr_half`, which is why 200 cells are visited per doublet.
  `ModuleInfo::pos` would bound r per MODULE instead of per layer and cut that;
  it is a cost win, not an accuracy one (the per-hit cut already uses the hit's
  own r). It becomes necessary for the tilted TBPS layers.
- Barrel only; nothing here touches the forward discs.
- The vertex finder is a greedy peak extractor with no vote weighting.
- Duplicates are not removed: 812 quads/ev against 349 found tracks.

## Reproduce

Run from the standalone BUILD directory; output lands in `test-seedgeom/`:

    cd /foo/matevz/mic-dev/current/src/standalone
    mkFit-external/mkfit-standalone-seedgeom/seedgeom-run.sh 20 myrun \
        'sg_layers(0,1,2)' 'sg_layer4(3)' 'sg_qwin(0.035)' 'sg_qwind(0.025)' \
        'sg_phiwind(0.002)' 'sg_vf(0.02,3,0.05)'

Interactively, inside an mkFit `--shell`:

    .x mkFit-external/mkfit-standalone-seedgeom/seedgeom-load.C
    sg_layers(0,1,2); sg_layer4(3); sg_reset();
    s.GoToEvent(1); sg_ev(s.event());
    sg_report(); sg_write("test-seedgeom/myrun");

Knobs: `sg_ptmin sg_d0max sg_layers sg_layer4 sg_qwin sg_qwind sg_phiwind
sg_phimargin sg_nr sg_grid sg_light sg_vf sg_vtxgap sg_blrcut sg_partb
sg_verbose`, plus `sg_dg(ev)` (production points) and `sg_rdg(ev)` (per-layer
radial structure).

---

# Addendum, same day: what "200 cells per doublet" was, and the radial axis

## The 200 was an artifact of MY grid, not a property of the method

`HitGrid` is 1024 phi x 512 z, i.e. **4x finer in phi and 25x finer in q than
production**: `LayerOfHits` uses 256 phi N-bins and `q_bin` 2.0 cm for PixB.
Per doublet that rectangle was 13 phi x ~15 z = 200 cells holding **2.25 hits**
-- 0.011 hits/cell, a scan of an empty array. At production N-granularity the
same query is ~6 cells. Do not read the 200 as a finding.

## The real over-scan is binning-independent, and there are TWO of them

    q bin-range half-width  =  qwin  +  |cot theta| * dr_half
                               0.035  +   ~1.0      *  0.53     cm

1. **Bin granularity.** PixB `q_bin` = 2.0 cm against a physics window of
   0.07 cm: ~28x more hits pulled than wanted, purely from bin size. This one
   is free to fix -- finer bins for seeding, collapse afterwards.
2. **Layer radial thickness.** Even with infinitely fine q bins the box is
   `2 qwin + |cot| dr` = 0.07 + 1.12 cm, i.e. **16x too wide**. No bin size
   removes it. Only making r a bin axis does.

## Layer radial structure IS discrete: two ladder shells

Measured hit-r distribution, event 1 (`sg_rdg(s.event())`):

| layer | shells | edges [cm] | gap |
|---|---|---|---|
| 0 | 2 | 2.750-2.889, 3.258-3.819 | exactly empty |
| 1 | 2 | 5.900-5.978, ... | exactly empty |

Each shell's internal width is the planar-module `R/cos(u)` tail -- a sharp
peak at the module centre radius with a tail outward -- so it is computable
from `ModuleInfo`, not stochastic.

## Making r a bin axis: measured, and it goes as 1/n_r

`sg_nr(n)` adds n radial sub-bins per layer, key `(iphi*nr + ir)*NZ + iz`, each
sub-bin getting its own z window. 5 events, layers 0,1,2+3:

| n_r | layer-c cells/doublet | **layer-c hits/doublet** | layer-d hits/triplet | quads | purity | eff |
|---|---|---|---|---|---|---|
| 1 | 186.1 | **2.037** | 0.331 | 823 | 0.7447 | 0.9339 |
| 2 | 199.8 | **1.138** | 0.194 | 823 | 0.7448 | 0.9339 |
| 4 | 227.3 | **0.657** | 0.123 | 823 | 0.7448 | 0.9339 |
| 8 | 282.4 | **0.412** | 0.087 | 823 | 0.7448 | 0.9339 |
| 16 | 392.6 | **0.287** | 0.069 | 823 | 0.7448 | 0.9339 |
| 32 | 613.0 | **0.225** | 0.060 | 823 | 0.7448 | 0.9339 |

Quads, purity and efficiency are **bit-identical** throughout -- same hit set,
less scanning. Fitting `hits = A + B/n_r` gives A = 0.166, B = 1.87, which
predicts 0.633 at n_r = 4 (measured 0.657) and 0.400 at n_r = 8 (0.412). So the
scanned area is `2 qwin dr + |cot| dr^2 / n_r` as the algebra says, the floor
is ~0.17 hits/doublet, and **n_r = 8 buys 4.9x of the available 12x**.

Cells go the other way (186 -> 282) because each sub-bin needs at least one z
bin. With production's coarse N-bins that term is small; with my over-fine grid
it dominates. The two costs trade, and the trade is what the binnor's M/N split
exists to manage.

## Concrete proposal

- `qbar` **already holds the hit's r** for barrel layers
  (`HitStructures.cc:109`), so the exact per-hit test costs nothing extra
  today. What is missing is r as a *bin/sort* axis.
- `binnor` is templated on exactly two axes, so a third is a template change.
  The cheaper route: the q axis carries **M=16 / N=8 bits, i.e. 256 fine bins
  inside each 2.0 cm N-bin (78 um)** -- far finer than the physics needs.
  Feeding `register_m_bins()` a q M-index whose high bits are a 3-bit radial
  sub-index would give r as a within-cell sort axis for **no new memory and no
  template change**, and hits in a cell would then be grouped by radial shell.
  Worth checking what else relies on the current q ordering inside a cell
  (`m_cons` is kept for preselection) before doing it.
- Independently: the window should be pT-dependent. Every efficiency loss
  measured above is at the low-pT end and is the flat window being ~1.5 sigma
  there.

---

# Addendum 2: what "efficiency" and "purity" actually mean here

## Definitions, from the code

Truth label of a hit is `MCHitInfo::mcTrackID_` reached through `Hit::mcHitID()`
— the **rec -> sim** direction, i.e. the one `bestTkIdx` arbitration writes.

- **findable (denominator)**: a sim track with (a) at least one hit carrying its
  label in **each** required layer, (b) `pT > pT_min`, (c) production point
  within `D0_max` of the beamspot. (b) and (c) matter: without them the
  denominator contains 0.2 GeV curlers the phi stencil can never accept, and
  the efficiency is quoted against an unreachable ceiling.
- **found (numerator)**: **at least one** built candidate has all of its hits
  carrying that same label.

So efficiency is **per TRACK, not per candidate**, and it is satisfied by any
one of a track's candidates. 804 quads/ev contain 600 true ones but only 341
distinct found tracks — ~1.8 true quads per found track. Per-candidate goodness
is `purity`, a different number.

It does not require the candidate to be unique or best, and it tests label
agreement, not that the hits are the *right* ones. The `mc_match` similarity
trap is mostly disarmed here by construction — the four hits must also satisfy
the phi stencil and a 250 um q prediction, so a geometrically impossible match
cannot pass — but a track's own overlap hit in the same layer does count, which
is correct.

## Purity has three outcomes, not two — 26 % of it is UNDECIDABLE

Only **70-77 %** of pixel-barrel hits carry a truth link at all (per layer:
0.7675 / 0.7466 / 0.7298 / 0.7021). Most of the rest belong to no sim track —
noise, sub-threshold secondaries — not to arbitration loss. 10 events,
layers 0,1,2+3, n_r = 8:

| outcome | /ev | fraction |
|---|---|---|
| **true** — all four hits one label | 600 | **0.7463** |
| **genuine fake** — four valid labels that disagree | 56 | **0.0692** |
| **undecidable** — >= 1 hit has no truth link | 148 | **0.1846** |
| (>= 3 of 4 hits share a label) | 680 | (0.8459) |

- **Decidable purity = 600/(600+56) = 0.9152.**
- **84.6 % of quads have at least 3 of their 4 hits from one sim track** — the
  most robust single statement, since it does not depend on the fourth hit
  having a link.

**Do NOT "correct" the purity by dividing by `prod(link fraction)` = 0.294.**
That was tried and gives 2.54: the model is wrong, because an unlabelled hit
is usually one that could never have belonged to a true quad in the first
place. The decidable ratio is the honest correction.

So the fair reading of the earlier tables: **purity 0.75 is a lower bound;
0.92 is the value among candidates the truth can adjudicate; the truth is
somewhere between and this sample cannot say where.** Efficiency (0.93-0.94)
is affected far less, because a track needs only one of its ~1.8 quads to be
fully labelled.

---

# Addendum 3: SCOPE, and why a single efficiency number is wrong

## Scope

**Pixel barrel only** — mkFit layers 0,1,2,3. **No pixel disks anywhere.** The
only exception is `qC` / `s3` / `s4` / `eta-ot`, which use layers 1,2,3 +
TBPS-P (layer 4).

**And the denominator caps the acceptance.** Requiring a hit in all four
barrel layers limits |eta| to the OUTERMOST layer's `asinh(z_half/R)`, which
for IT4 is **1.12**. Measured: 2702 / 2643 / 1962 / 203 findable tracks per
20 events in |eta| 0-0.4 / .4-.8 / .8-1.2 / 1.2-1.6, and **zero above 1.6**.
So this whole study covers the central ~40 % of the eta range, with 71 % of
its own weight below |eta| 0.8. Forward coverage needs mixed barrel/disc
layer triples, which the code does not do.

## Efficiency, pixel barrel 0,1,2+3, 20 events (`eta-pixb.txt`)

| pT \ \|eta\| | 0-0.4 | .4-.8 | .8-1.2 | 1.2-1.6 | all |
|---|---|---|---|---|---|
| 0.9-1.2 | 0.947 | 0.935 | **0.833** | **0.707** | 0.905 |
| 1.2-1.6 | 0.976 | 0.958 | 0.907 | 0.796 | 0.946 |
| 1.6-2.2 | 0.979 | 0.975 | 0.974 | - | 0.976 |
| 2.2-3.2 | 0.986 | 0.980 | 0.989 | - | 0.981 |
| 3.2-5 | 1.000 | 1.000 | 0.974 | - | 0.993 |
| >5 | 1.000 | 1.000 | 0.963 | - | 0.992 |
| **all pT** | **0.967** | **0.955** | **0.897** | **0.783** | 0.940 |

The corner that matters runs **0.947 -> 0.707**; high pT is 1.0 everywhere.
Flat-in-pT windows again.

## Purity vs |eta|, from the candidate's own cot theta (no truth needed)

| \|eta\| | quads/ev | true | genuine fake | undecidable | decidable |
|---|---|---|---|---|---|
| 0-0.4 | 294 | 0.786 | 0.063 | 0.152 | 0.926 |
| .4-.8 | 284 | 0.774 | 0.066 | 0.160 | 0.922 |
| .8-1.2 | 204 | 0.727 | 0.072 | 0.201 | 0.910 |
| 1.2-1.6 | 41 | 0.306 | 0.208 | 0.486 | 0.595 |

The last row is the acceptance edge showing up as junk.

## Layers 1,2,3 + TBPS-P: the eta trend REVERSES at low pT (`eta-ot.txt`)

| pT \ \|eta\| | 0-0.4 | .4-.8 | .8-1.2 | all |
|---|---|---|---|---|
| 0.9-1.2 | **0.524** | 0.667 | **0.733** | 0.635 |
| 1.6-2.2 | 0.809 | 0.916 | 0.936 | 0.879 |
| 2.2-3.2 | 0.966 | 0.972 | 0.965 | 0.965 |
| >5 | 1.000 | 1.000 | 0.962 | 0.992 |

Purity: true 0.50, undecidable 0.39, decidable 0.81. **No explanation offered
for the reversal** — the layer-4 dz window is 2 sigma at low pT, so this
configuration is window-limited before it is anything else, and the trend
should not be interpreted until the windows are pT-dependent.

## Correction to addenda 1-2 and to the main body

Every "flat in pT" statement there was flat because it was averaged over
|eta| < 1.2 with most of the weight below 0.8. Quote the table, not the
number: **0.95-0.98 in the central barrel above pT 1.6, degrading toward both
low pT and the eta edge, and entirely unmeasured forward.**

---

# Addendum 4: where the time actually goes (and what it is NOT)

Ryzen 7 2700 (Zen+, AVX2 only, no AVX-512), **single thread, scalar**, box
under other load. Marginal cost from (t[11 events] - t[1 event])/10.

**Repeatability first: 190 / 246 / 213 ms on an identical config** -- +-15 %.
Nothing below that is resolvable, and several things tried are not.

| | ms/ev |
|---|---|
| Part A alone (event read + sim-track truth map) | **98** |
| Part A + Part B, layers 0,1,2+3, light, grid 256x32, nr=1 | **216 +- 25** |
| **=> Part B proper** | **~118 +- 25**, i.e. **~140 ns/doublet** |

So a third to a half of the earlier "0.30 s/event" was event I/O and my truth
bookkeeping, not the algorithm.

## Hypotheses tested, and two of three were WRONG

| hypothesis | result |
|---|---|
| cells/doublet is the cost | **NO.** 186 -> 5.8 cells/doublet (grid 1024x512 -> 256x32) moved 396 -> 393 ms |
| L2 vs L3 residency | **partly** -- 493 KB (fits the 512 KB/core L2) vs 8.5 MB gives 226 vs 267 ms, **1.2x** |
| scalar `double` divisions in the window test | **unresolved** -- precomputing per-hit `1/r` and a branch-free `wrap_pi` landed inside the +-15 % noise |

**Quadruplet count and purity are IDENTICAL (808/809, 0.914) across every grid
granularity and every `n_r`** -- the grid is purely a lookup structure, which
is the validation that all of this is behaviour-preserving.

Also measured: **`n_r > 1` costs more than it saves at these settings.** It cuts
hits/doublet 4.06 -> 1.25 but the grid build (`start` array grows x n_r) and the
extra cell lookups outweigh it: 226 -> 394 ms. The n_r win is real in *hit
tests* and not yet real in *wall clock*; it would need the build cost amortised
(the binnor is built once per event for the whole finding, not per study).

## What is left, and what it means for a CPU-vs-GPU argument

At ~140 ns/doublet the inner loop is roughly 4 loads + 6 flops + 2 compares per
candidate. On an 8-wide AVX2 core the arithmetic floor is single-digit ns, so
there are one to two orders of margin **in principle** -- but none of it is
demonstrated here, and the one optimisation attempted was below the noise
floor. **Do not quote a projected speedup from this work.**

The structural facts that ARE established:

- **The working set is small: 493 KB for four pixel layers**, i.e. L2-resident
  per core. This is latency/cache work, not bandwidth streaming -- the regime
  where a GPU's main advantage does not apply.
- **The doublet count is the hard floor**: 837 k/ev at pT_min 0.9, 1.6 M at 0.3,
  x4.7 at D0_max 1 cm. Every optimisation is per-doublet cost; none of it
  reduces this.
- **Trivially parallel in phi**, no communication.

No Patatrack timing was measured on this sample, and none should be quoted
until it is.

---

# Addendum 5: extending outward — bias vs width, and where propagation belongs

## The long path does NOT break the geometry. It widens it, anisotropically.

TRUE quads, 20 events (`bias-pixb.txt`, `bias-ot.txt`):

| prediction | dphi median | dphi sigma | dz median | dz sigma |
|---|---|---|---|---|
| -> layer 3, 4 cm hop inside IT | -7e-6 rad = **-0.011 sigma** | 0.632 mrad | -1.0 um = **-0.012 sigma** | 84 um |
| -> layer 4 TBPS-P, 10 cm across the IT/OT gap | -7e-6 rad = **-0.006 sigma** | 1.093 mrad | -9.6 um = **-0.015 sigma** | **630 um** |

**Bias below 1.5 % of sigma in both directions, even across the gap.** Uniform
B, no material, no energy loss, closed-form circle -- unbiased at r = 25 cm.
What degrades is width, and it is anisotropic: **sigma_q grows 7.5x while
sigma_phi grows 1.7x.** That is the phi-discriminates / q-contains asymmetry,
getting worse with path length.

## Consequence: extending with a GEOMETRIC q window costs purity

Same events, layers 1,2,3 + TBPS-P, q window 0.15 cm (2.4 sigma):

| | 0,1,2+3 | 1,2,3+4 |
|---|---|---|
| quads /ev | 832 | 1658 |
| 4th-layer rejection | 58.2x | 15.1x |
| purity | **0.735** | **0.286** |
| efficiency | 0.940 | 0.777 |

So the naive extension costs more than it buys. The fix is not a tighter flat
window -- 0.15 cm is already only 2.4 sigma -- it is a window that is *per
track*. sigma_q out there is MS-dominated, and MS is knowledge only a
propagated covariance carries: a flat window must be sized for the worst pT in
a steeply falling spectrum.

**CORRECTED in addendum 7 (2026-09-24): the 630 um is NOT multiple scattering.**
It is flat in pT, and a truth split shows it is the P sensor's own q extent: a
1.5 mm macro-pixel, which is uniform, so the IQR-based sigma reads 1.28x its
standard deviation. The prediction's own error is the smaller term (282 um at
pT 0.9-1.5, falling as 1/p). A propagated covariance therefore cannot recover
the purity. A per-hit containment window can.

## Where the handover to propagation belongs: the item count decides

| stage | items /ev | cost/item affordable in ~100 ms |
|---|---|---|
| doublets | 837 000 | ~120 ns -- geometry only, no state exists |
| triplets | 66 000 | ~1.5 us |
| **quadruplets** | **800** | **~100 us** -- a full propagation is free |

Each stage is 10-80x fewer items, so the budget stays flat if each is ~10x
dearer. **After the 4th point there are 800 objects/event**, and a real
`propagateHelixToPlaneMPlex` with material and parametric B costs nothing
measurable there. The triplet stage is still 66 k items -- 13x too many.

So: **pure geometry where there is no state and 10^6 items; propagation from
the 4th point outward.** Not an arbitrary line.

## Dead / noisy modules: the redundancy is worth 6 points of acceptance

With per-layer coverage 98.4 % (`an-layercov`, pixel barrel):

- require 4 of 4: `0.984^4` = **93.7 %**
- require any 4 of 6 (add TBPS-P and a disc): **99.9 %**

A floor no window tuning can reach, and larger once real dead/noisy modules are
included. `LayerOfHits::suckInDeads()` / `isBinDead(pi,qi)` already exist.
Structurally free: "any 4 of 6" is not C(6,4) searches, it is the
extend-and-confirm chain with a hole counter.

## Anchoring improves before the 4th point is even added

3-point curvature `(pT_fit - pT_true)/pT_true` robust sigma: **0.044** from
layers 0,1,2 (lever 3 -> 10.5 cm) vs **0.035** from 1,2,3 (lever 6.2 -> 15),
consistent with 1/L^2 sagitta scaling. A 4-point fit through TBPS-P at r = 25
nearly doubles the lever again -- and curvature at handover is what sets
mkFit's downstream window sizes.

## Reading

The three motivations for extending outward -- acceptance, fake/duplicate
reduction, anchoring -- are **not independent**. Acceptance and anchoring
already pay; purity currently does NOT (0.735 -> 0.286); and the fix for purity
is the same propagated covariance that makes the anchoring real. One change,
three payoffs.

Caveat: pixel barrel, |eta| < 1.2, one sample, and the TBPS-P configuration is
still window-limited at low pT, so its 0.777 efficiency is a floor rather than
a measurement of the method.

---

# Addendum 6 (2026-09-24): layer-c phi window from a linear-in-r extrapolation

Suggested by the MkFinderV2p2 session (recotracker-02), by analogy with its
line pre-cut. `sg_philin(mode, margin)`, default **off**, so the standing check
still holds (`ref-lib.txt` bit-identical, re-verified).

The layer-c window was centred on `phi_b` with the generic half-width
`(rc-rb)/(2 R_min) + D0_max (1/rb - 1/rc) + margin`. With the switch on it is
centred on `phi_b + slope (r_c - r_b)`, where `slope = (phi_b - phi_a)/(r_b - r_a)`.
The per-hit test uses the hit's own `r_c`. Extrapolating linearly in r absorbs
the curvature, since `phi(r) ~ phi0 - r/2R`. The D0 term goes as `D0/r` and is
not linear, so the tolerance is
`D0_max |1/rc - 1/rb - (rc-rb)/(rb-ra) (1/rb - 1/ra)| + margin`. That is
0.0145 rad at radii 3.3 / 6.4 / 10.8 cm with D0_max = 1 mm. Mode 1 applies the
linear test alone. Mode 2 applies it AND the generic band.

**Mode 1 alone drops the pT_min bound at layer c.** The a->b doublet window
bounds pT only weakly over its 3 cm lever arm. So mode 1 admits sub-threshold
tracks: true quads go 612 -> 891 /ev, while the findable (pT > 0.9)
efficiency does not move. Use mode 2 for a like-for-like comparison.

20 events, layers 0,1,2 + 3, grid 256 x 32, `sg_light(true)`, two reps each.
The physics output is identical between reps. Timing comes from the new
in-code timer around the search loop only, so event I/O and grid build are
excluded. Load was 1.5-3.4 on 16 threads.

| | off | mode 1, m 3 mrad | **mode 2, m 2 mrad** | mode 2, m 3 | mode 2, m 5 |
|---|---|---|---|---|---|
| layer-c cells / doublet | 5.82 | 5.19 | **4.98** | 5.19 | 5.49 |
| layer-c hits touched / doublet | 4.21 | 3.75 | **3.60** | 3.75 | 3.97 |
| triplets on phi alone / doublet | 2.46 | 1.36 | **0.87** | 0.91 | 1.00 |
| triplets after q / ev | 48458 | 27763 | **17585** | 18447 | 20154 |
| quads / ev | 833 | 1090 | **736** | 739 | 745 |
| true quads / ev | 612 | 891 | **605** | 605 | 606 |
| raw purity | 0.735 | 0.818 | **0.822** | 0.819 | 0.814 |
| decidable purity | 0.908 | 0.941 | **0.946** | 0.945 | 0.942 |
| genuine fakes / ev | 62 | 56 | **35** | 35 | 37 |
| findable efficiency (of 7510) | 0.9395 | 0.9389 | **0.9389** | 0.9389 | 0.9390 |
| search loop, ms/ev, rep 1 / rep 2 | 171.9 / 175.3 | 154.5 / 154.1 | **147.8 / 147.3** | 153.1 / 154.2 | 160.3 / 162.7 |

Readings:

- **Mode 2 at 2 mrad is 14 % faster in the search loop** (173.6 -> 147.6
  ms/ev, averaged over the two reps). The rep-to-rep spread is 2 %. It also
  cuts genuine fakes 62 -> 35 and raises decidable purity 0.908 -> 0.946.
- **It cuts phi-only triplets 2.8x but hits touched only 15 %.** At 256 phi
  bins one bin is 24.5 mrad. Narrowing the window from ~0.036 to ~0.02 rad
  therefore goes from about 5 bins to about 3, and the z window then decides
  how many hits are touched. The saving comes from fewer triplets reaching the
  q test and the 4th layer, not from fewer cells.
- **Efficiency costs 0.06 %, about 4.5 of 7510 findable tracks, all in
  pT 0.9-1.2.** The loss is the same at margins 2, 3 and 5 mrad, so it is not
  a multiple-scattering margin. Not investigated.
- A larger margin only costs. 2 mrad is the best point measured, and nothing
  below 2 mrad was tried.

**A defect of mine found while timing this, recorded because it cost 40x.** The
first version of the tolerance lambda captured `[=]` and read
`gb.invr[ib]` inside, so it copied the whole layer-b `HitGrid` on every
doublet: 5.7 us/doublet instead of 0.2. The physics output was unaffected,
which is why only the timer showed it. The lambda now captures named scalars.

**Open lead from the same exchange, not tested.** Addendum 5 attributed the
7.5x growth of sigma_q at the TBPS-P 4th point (84 -> 630 um) to multiple
scattering. recotracker-02 found that on tilted TBPS modules a hit's `qbar` is
its centroid radius. The true crossing lies anywhere along the sensor, which
has a radial component, so a z prediction at the hit's own r shifts by
`|cot theta| * hl_fac * sigma_r(hit)`. Without that term their pre-cut lost
44 % of TBPS pairs. For a P macro-pixel (half-length 0.075 cm, tilt up to 72
deg) the shift is of order a few hundred um at |eta| ~ 1. That is the same
order as the width attributed to MS. `LayerOfHits::hit_qbar_half_extent(i)`
now caches the term. **Bin the TBPS-P q residual by `|cot theta| * qbar
half-extent` before running the propagation test**, because the propagation
test assumes the width is MS.

---

# Addendum 7 (2026-09-24): the TBPS-P q width is the sensor, not multiple scattering

Addendum 5 attributed the 630 um q width at the TBPS-P 4th point to multiple
scattering and made a propagated covariance "the next measurement". This
addendum tests that attribution, and it fails.

Configuration: layers 1,2,3 + 4 (TBPS-P), 20 events, 4th-layer q window 0.5 cm
so the distribution is not clipped. Output `tbps-q050.txt` (April HLT sample)
and `tbps-d121-truth.txt` (D121 ttbar PU200, `--read-sim-hit-states`,
`--geom CMS-phase2-Run4D121`). The two samples give the same widths to within
7 %, except one sparse bin (flat modules at |cot| 0.75-1.0, n = 110).

## 1. The width is flat in pT, so it is not multiple scattering

Robust sigma of the residual of true quads [um], April sample, binned by
`|cot theta| * sigma_r` of the 4th hit:

| \|cot\| sigma_r | pT 0.9-1.5 | 1.5-3 | 3-10 | > 10 |
|---|---|---|---|---|
| < 20 um (flat modules) | 536 | 538 | 523 | 529 |
| 100-200 um | 609 | 594 | 491 | -- |
| 200-400 um | 772 | 721 | 734 | 798 |
| > 400 um | 911 | 919 | 804 | -- |

Multiple scattering would shrink the width 5-10x between the first and last
column. It moves by at most 12 %.

## 2. The truth split: the hit's own error dominates, the prediction's is MS

With a truth state per sim hit, the residual splits at the true crossing
radius into `e_hit = hit - truth` and `e_pred = prediction - truth`. D121,
9660 true quads, all with a state:

| \|cot\| | flat: e_hit | e_pred | tilted: e_hit | e_pred |
|---|---|---|---|---|
| 0-0.25 | 563 | 161 | -- | -- |
| 0.25-0.5 | 542 | 180 | 557 | 219 |
| 0.5-0.75 | 531 | 167 | 586 | 207 |
| 0.75-1.0 | 506 | 144 | 730 | 270 |
| 1.0-1.4 | -- | -- | 791 | 353 |

| pT | 0.9-1.5 | 1.5-3 | 3-10 | > 10 |
|---|---|---|---|---|
| e_hit | 638 | 645 | 590 | 642 |
| e_pred | **282** | **153** | **78** | **45** |

`e_pred` falls as 1/p, which is multiple scattering, and it is the smaller term
at every pT. `e_hit` is flat in pT and carries the width.

## 3. The hit's error agrees with its own covariance, once the estimator is right

The robust sigma used throughout this study is IQR/1.349, correct for a
Gaussian. A single macro-pixel is UNIFORM over its length. For a uniform
distribution on [-L, L] the IQR sigma is 0.741 L while the standard deviation
is 0.577 L, a ratio of **1.284**. For the P macro-pixel (L = 750 um) that is
556 um against 433 um.

The prediction from the hit's own covariance is exact when it includes the
signed r-z term: `var(z - cot r) = var_z - 2 cot cov_rz + cot^2 var_r`. On a
tilted sensor facing the IP, a displacement along the sensor moves r and z
together, and the terms ADD. A first version that used magnitudes only had them
cancelling, and was wrong. Scaled by 1.284:

| | covariance predicts x 1.284 | measured e_hit |
|---|---|---|
| flat, all \|cot\| | 543 | 506-563 |
| tilted, \|cot\| 0.25-0.5 | 548 | 557 |
| tilted, 0.5-0.75 | 630 | 586 |
| tilted, 0.75-1.0 | 715 | 730 |
| tilted, 1.0-1.4 | 827 | 791 |

**Agreement is within 7 % everywhere, so the covariance is correct.** The
excess over the covariance that the first pass of this addendum could not
explain, about 330 um at central eta, was the estimator.

## Consequences

- **The recorded "next measurement" is withdrawn.** A propagated covariance at
  the 4th point would carry the MS term, which is the minor one. It cannot
  bring the 0.286 purity back toward 0.735.
- **What can is a per-hit containment window**: the hit's own q extent from its
  covariance, projected with the signed r-z term, plus a 1/p prediction term.
  The flat 0.15 cm window is twice the half-length of a flat-module
  macro-pixel. In the tilted rings at |cot| 1.0-1.4 the hit's projected
  half-extent is about 0.11 cm (sqrt(3) x 644 um), before the prediction term
  is added. So one number is too loose in one place and too tight in the other.
- **recotracker-02's `|cot| * qbar_half_extent` term in the pre-cut is a valid
  upper bound**, by the triangle inequality on the exact expression. It is the
  cot^2 var_r term without the correlation.
- **Every IQR sigma in this README is 1.28x high for a uniform distribution.**
  It is right for the pixel-barrel residuals, which are close to Gaussian. It is
  not right for macro-pixels or strips. Report a standard deviation, or a
  half-extent, for those.
