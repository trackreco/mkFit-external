# mkfit-seeding: a geometric quadruplet seeder on the mkFit binnor

Milestone 1 of the seeding work that started as the study in
`../mkfit-standalone-seedgeom/`: barrel quadruplets on layers 0,1,2 + 3,
reproducing that prototype's output exactly, then made fast.

## Files

| file | what |
|---|---|
| `SeedLayer.h` | a layer's hits in bin order, struct-of-arrays, registered as `LayerOfHits::suckInHits()` does, with the axis types as template parameters |
| `SeedFinder.h` | `find_quads()`, the scalar port. Per-hit cuts copied from the prototype expression for expression; only the fetch is new |
| `SeedFinderStaged.h` | the same search as a pipeline of flat loops in blocks of a-hits, optionally with candidate generation fused into the triplet test. `finish_triplets()` is the helix and 4th-layer stage all variants share |
| `SeedFinderBMajor.h` | the b-hit as the outer loop; per b-hit a z-bucketed list of the c-hits that pass the phi and r tests, shared by every doublet through it |
| `SeedMath.h` | the per-triplet and per-quad arithmetic, templated on a policy: `ArithRef` (double, the prototype's expressions), `ArithFast` (float + vdt, circle-centre form), `ArithFastK` (float + vdt, curvature form). Every cut also returns its signed margin |
| `SeedMargins.h` | `eval_quad<A>()`: every cut of one quad in arithmetic A, pass flag and margin. It calls the functions `finish_triplets<A>()` runs, so it cannot drift from the finder |
| `SeedFinderTile.h` | step B in float: per b-hit kernels K1 (doublets and slopes), K3 (the c list from a phi-only layer c), K4 (the r-z slope pre-filter in slope buckets, then the exact z test). `--tile`; `--tile-brute` keeps the dense-tile form for comparison |
| `SeedStats.h` | `seedfind --stats`: trip counts and value ranges of the per-pair stages, walking the b-major finder's windows with its cuts (step B0). Its doublet and triplet counts must equal the finder's |
| `seedfind.cc` | standalone driver: geometry plugin and events as `mkFit.cc` loads them, timers around the fill and the search only |
| `SeedOT1.h` | `seedfind --ot1`: a compatible hit in OT1 (P or S) for every quad, with truth |
| `SeedSeq.h` | `seedfind --sequences[-all]`: the layer sequence every sim track crosses, rec -> sim and sim -> rec |
| `SeedSurf.h` | the surface finder: quads on any four layers, barrel or disc, double precision; `SurfChain`, the feed-forward chain over the layers in crossing order; `SurfOwnership`, the phase-space ownership test |
| `seedsurf.cc` | its driver, with truth by pattern and for the union, residual dump, cleaning, ownership, `--bind`, `--chain` |
| `windows-D121/` | the c and d windows per layer combination, as `seedsurf` options, one combination per line: `pixel.txt` (11 pixel patterns), `ot1p.txt` (3 triplet + OT1-P), `skip.txt` (29 combinations with a skipped layer) |
| `seedsurf-chain.sh` | runs `seedsurf` in the chosen configuration, the chain with a late start |
| `Makefile` | flags from the build's `make echo-aclic`, re-read on every build; ROOT only if `libMkFitCore.so` links it |

## Build and run

The shared build in `src/standalone/` is rebuilt, and switched between the
trace build and a ROOT-off build, by other work. So this runs against an
**isolated build**:

- `CMSSW_14_1_0_pre0-p2p/src-seeding/`: a detached, sparse worktree of the
  CMSSW repo at `59b10cfe7db`, with the committed `Makefile.config` (ROOT-off,
  no trace).
- `src-seeding/standalone/`: its build directory. `mkFit-external` there is a
  symlink to `src/standalone/mkFit-external-seeding/`, the worktree of this
  branch, and the geometry `.bin` files are symlinks into the shared build.

```
B=/foo/matevz/mic-dev/CMSSW_14_1_0_pre0-p2p/src-seeding/standalone
make BLD=$B                                  # -> $B/test-seedgeom/bin/seedfind
make BLD=$B ARCH=-march=native NAME=seedfind-native
cd $B && LD_LIBRARY_PATH=. test-seedgeom/bin/seedfind \
    --input-file /foo/matevz/mic-dev/trackingNtuple_HLT_2026_March.bin \
    --num-events 5 --reps 3 --qbin-c 2.0 --qbin-d 0.5 --bmajor --lbin 0.25 \
    --dump test-seedgeom/quads/x.txt
# step A: float arithmetic, and the difference tool against a reference list
... --arith fastk --margins test-seedgeom/quads/A-ref-50ev.txt
```

`--arith ref|fast|fastk` selects the per-triplet arithmetic (b-major, staged
and fused finders; the scalar port is always `ref`). `--margins REF` runs the
difference tool, below. `--eps-cm` / `--eps-rad` set its epsilon (default
1e-4 cm = 1 um and 1e-6 rad), `--margins-print N` how many differing quads it
lists.

**Defaults changed 2026-09-26 (maintainer): the doublet-slope phi cut is on
(`phi_lin` 2, 2 mrad) and the three fixed windows are 1.4x wider (3rd-hit z
0.049 cm, 4th-hit z 0.035 cm, 4th-hit phi 2.8 mrad).** pT > 0.9, 20 events:
efficiency 0.972, fake 0.21 by the consistent rule, against 0.940 and 0.19
before. **Every number and reference list in this file dated before that was
taken with the old values; `seedfind --first-look` restores them** (put any
window option after it). With `--first-look` the 50-event list is identical to
the one before the change. `--stats` does not apply the phi cut, so with it on
its triplet count no longer equals the finder's.

**Residuals of the fixed windows on true quads, `seedfind --residuals OUT.txt`
(2026-09-26).** One row per findable track, found or not: truth pT, p and eta,
the triplet circle's R (pT_est = 0.0114 R [GeV, cm]) and cot theta, and the
|residual| of c_z, d_phi and d_z on the track's own hits. 100 events at pT_min
0.5, 91 145 tracks. pT_est / pT has median 0.996 and 16-84 % range 0.951-1.044.
- **q95 x pT_est is constant from 0.5 to ~2 GeV**: about 270 um GeV for c_z,
  1.9 mrad GeV for d_phi and 255 um GeV for d_z. That is the multiple-scattering
  1/pT scaling of the f x pT result.
- **It grows with |cot theta|.** For d_phi it follows sin(theta)^-0.5 (x1.38 at
  |cot| 1.6-2.5 against x1.40 predicted). For c_z and d_z it rises somewhat more
  slowly than sin^-1.5 (x2.3-2.6 against x2.7).
- **Above pT_est ~3 GeV the tails grow again. That is migration, not high-pT
  tracks.** By TRUE pT, tracks above 3 GeV are found 99 %. By pT_est, 24 % of
  the tracks above 10 GeV have pT_est / pT outside [0.8, 1.25]: a circle made
  too straight by a kink or a bad hit. A window that shrinks with pT_est
  therefore needs a constant floor.

**4th-hit windows from the curvature, `--win-scaled F PT_KEEP` (2026-09-26):
built, measured, and NOT better than widening the fixed windows.** Per triplet,
w = F sqrt(a^2 + (b x)^2), x = g(theta) / pT_est, with a, b the 95 % contour
fitted on true quads whose curvature estimate is within 20 % (d_phi a 0.2 mrad,
b 1.6 mrad GeV; d_z a 63 um, b 177 um GeV). 1 / pT_est is capped at 1 / pT_min.
Below PT_KEEP the window is the wider of the scaled and the fixed one, or with
`--win-low-fixed` the fixed one exactly. `--win-floor A_PHI A_Z` sets a. With
the option off, every finder's 50-event list is bit-identical to before. With
it on, the tile, b-major and brute lists agree, and why-missed finds no track
that passes every cut. 20 events, pT > 0.9, fakes by the consistent rule:

| 4th-hit windows | efficiency | 0.9-1.2 / 1.2-2 / 2-5 / >5 GeV | fake |
|---|---|---|---|
| fixed x1.0 | 0.9443 | 0.912 / 0.959 / 0.986 / 0.992 | 0.129 |
| fixed x1.4 (default) | 0.9722 | 0.957 / 0.979 / 0.992 / 0.992 | 0.211 |
| fixed x1.8 | 0.9787 | 0.969 / 0.983 / 0.993 / 0.992 | 0.288 |
| fixed x2.2 | 0.9828 | 0.976 / 0.985 / 0.994 / 0.992 | 0.364 |
| fixed x2.6 | 0.9842 | 0.979 / 0.986 / 0.994 / 0.992 | 0.432 |
| scaled F 2.0, no cap | 0.9748 | 0.973 / 0.976 / 0.980 / 0.959 | 0.289 |
| scaled F 2.5, cap | 0.9806 | 0.977 / 0.982 / 0.986 / 0.984 | 0.340 |
| scaled F 2.5, cap, floor 0.8 mrad / 100 um | 0.9826 | 0.978 / 0.984 / 0.990 / 0.992 | 0.376 |
| scaled F 2.0 above 1.2 GeV, fixed below | 0.9722 | 0.958 / 0.979 / 0.989 / 0.992 | 0.195 |

Every scaled point lies on or just below the fixed curve. The best one keeps the
default's efficiency with 8 % fewer fakes, and it is 0.3 points lower at 2-5 GeV
(1.2 sigma). Without the cap, combinatorial triplets with a small pT_est opened
the widest windows. A floor under ~0.8 mrad / 100 um lost high-pT tracks with a
kink, which the curvature estimate places at too high a pT. Not measured: the
fake rate binned in pT_est, which would show where the fakes are. Raw logs in the
working report's prep/res-2026-09-26/s-*.log.

**Seed purity against pT_est, `seedfind --truth` (2026-09-26).** The truth
report now bins the quads in pT_est, from the circle through a, b, c, and splits
the consistent-rule fakes by whether a, b, c carry one label. 100 events, current
defaults, fake among decidable: 0.13 at 0.7-0.9 GeV, 0.16-0.18 at 0.9-2, 0.23 at
2-3, 0.37 at 3-5, **0.58 at 5-10 and 0.72 above 10 GeV**. Only 3.2 % of fakes are
a true triplet with a wrong d (2.2 % with the true quad also in the list), so
choosing the best d per triplet would remove about 2 % of fakes. The scaled
4th-hit windows with fixed windows below 1.2 GeV cut the fakes above 5 GeV:
f 2.0, floor 0.8 mrad / 100 um gives 0.35 / 0.46 at unchanged efficiency;
f 1.5, floor 0.3 mrad / 50 um gives 0.12 / 0.10, but costs 2-4 points of
efficiency above 1.2 GeV true pT (20 events: 0.961 / 0.963 / 0.951 at 1.2-2 /
2-5 / >5 against 0.979 / 0.992 / 0.992). Logs in the working report's
prep/res-2026-09-26/pe-*.log.

**A compatible hit in OT1, `seedfind --ot1 OUT.txt [--ot1-layer L] [--ot1-fetch
RAD CM]` (`SeedOT1.h`, 2026-09-26).** A truth study; the finder is unchanged. Each
quad's circle through a, c, d (tangent at d) and its r-z line are extended to
layer L, and every candidate is predicted at its own radius, since the OT1
modules are tilted. The hit's error on each residual comes from its full
covariance, including the r term. The analysis is the working report's
`prep/ot1-analysis.py` (window, match rates), `ot1-classes.py` (three classes)
and `ot1-zfloor.py`, with raw output and tables in `prep/ot1-2026-09-26/`. 100
events, pT_min 0.9, current defaults.

The window is w_phi = sqrt((3 s_phi)^2 + (b_phi / pT_est)^2) and w_z = sqrt(3)
s_z + b_z / pT_est, each with the hit's own s. b is the envelope over pT_est bins
that puts at least 95 % of the true hits of good-curvature quads (pT_est within
20 % of the truth) inside the window in every bin. For OT1-P: **b_phi 3.20 mrad
GeV, b_z 316 um GeV**. The per-bin products are 2.7-3.2 mrad GeV in phi, so the
phi prediction scales as 1/pT; in z the hit's extent alone covers above ~2 GeV.
A match is |dphi| < f w_phi and |dz| < f w_z.

**The first fit had a 486 um z floor, and it was the fit, not the physics.** It
added the q95 of the residual and the bin's median hit term in quadrature.
Deconvolving the true-hit dz as U(-h, h) + N(0, sigma) with each hit's own h
(`ot1-zfloor.py`) gives sigma = 0, on q68 and on q95 alike, above 3 GeV and in
every |cot| bin: there the residual is the macro-pixel alone. h is not one
number: 734 um on the flat modules, 940-1110 um at |cot| > 0.7 from the r-z
projection on the tilted ones, so a median h under-states the pool's q95. With
that fit the z window was 15-20 % wider than the corrected one (880 against 754
um above 10 GeV, 1330 against 1130 um at 0.8 GeV). At f = 2 the corrected window
moves the fake match rate by under one point in every pT_est bin above 0.7 GeV
(0.477 -> 0.472 at 0.7-0.9, 0.095 -> 0.094 at 5-10) and 2.3 points at 0.5-0.7
GeV; true quads with their own hit stay at 0.988. Phi does most of the
rejection: above 5 GeV, 12.5 % of fakes pass phi alone, 37 % pass z alone, and
10.1 % pass both.

- **OT1-P (layer 4):** 93.3 % of true quads have a hit of their own track there
  (acceptance, the P sensor missing and module inefficiency together).
- **OT1-S (layer 5) does not discriminate.** Its z window is the strip. P or S
  keeps 98.5 % of true quads and 62 % of fakes, against 96.2 % and 32 % for P.
- **Three classes at f = 2, by pT_est** (fractions of true / fake quads):

| pT_est | true /ev | P | S only | none | fake /ev | P | S only | none |
|---|---|---|---|---|---|---|---|---|
| 0.5-1 | 287.3 | 0.964 | 0.024 | 0.012 | 55.4 | 0.442 | 0.301 | 0.257 |
| 1-2 | 253.1 | 0.962 | 0.021 | 0.017 | 55.4 | 0.323 | 0.305 | 0.372 |
| 2-5 | 69.7 | 0.956 | 0.025 | 0.018 | 27.1 | 0.180 | 0.300 | 0.521 |
| > 5 | 9.3 | 0.952 | 0.022 | 0.027 | 17.5 | 0.101 | 0.275 | 0.623 |

  Above 5 GeV a P match keeps 95 % of true quads and 10 % of fakes. The P-matched
  set there is 8.9 true against 1.8 fake quads per event (83 % of the decidable),
  against 35 % in the whole list. The 5 % of true quads without a P match (0.45
  per event) cannot be told from the 16 unmatched fakes at the seed level, so a
  P match can rank or flag seeds but not veto them.
- **It already depends on |eta| inside the barrel.** Above 5 GeV the fake P-match
  rate is 0.115 / 0.139 / 0.064 at |eta| < 0.5 / 0.5-0.9 / > 0.9. These are
  pixel-barrel quads, |eta| up to ~1.1; OT1's tilted rings at higher |eta| need
  their own measurement.

## Tried and not kept

Every change that was built and measured, and then dropped or kept off by
default because it bought nothing or cost more than it saved. Each row points
at the section below that has the full numbers. Times are search ms per event
on black unless stated. The quad list was identical in every speed row.

| what | measured | verdict | section |
|---|---|---|---|
| branch-free binary search in the per-b c list | no gain over the branching one; ~48 ns per doublet, latency bound | replaced by z buckets | Measurements |
| dense slope tile (`--tile-brute`) | 45.2 ns per doublet against 29.5 with slope buckets; volume bound | kept as an option for comparison | Step B, first half |
| K4 masks at fixed per-doublet slots | K4 5.9 against 5.5 ms | dropped | K4: three restructurings |
| K4 bucket against bucket (doublets sorted too) | 8.4 against 5.5 | dropped | same |
| K4 4-wide SSE step | 6.4 against 5.5 | dropped | same |
| int16 per-pair arithmetic (`mkfit-seeding-int16`) | -1.5 to -5 % on four machines; the test arithmetic is ~7 % of the time | branch kept, not merged | After the float kernels; Four machines |
| forced 512-bit vectors on Skylake-SP | helix -18 %, every other stage +4-5 %, net loss | GCC's 256-bit preference kept | Four machines |
| `__builtin_prefetch` in the 4th-layer fetch | stage 6 1.95 -> 2.47 ms; IPC 2.2, 2.2 % L1 misses | reverted | Before and after on four machines |
| 2, 4 or 8 counter sets in K3's sort | K3 4.33-4.70 against 4.44 with cursors, 4.23 / 4.30 against 3.93 with ranks | dropped | K3's counting sort with ranks |
| drop the per-entry slope tolerance t_c | K3 -0.25-0.4 ms, K4 +0.35-0.8 ms (32 % more pre-filter hits) | dropped | same |
| 4th-hit windows from the triplet curvature (`--win-scaled`) | every configuration on or below the fixed-window curve; best: same efficiency, 8 % fewer fakes | built, off by default | 4th-hit windows from the curvature |
| tight scaled windows, f 1.5, floor 0.3 mrad / 50 um | fakes above 5 GeV 0.58 -> 0.12, efficiency -2 to -4 points above 1.2 GeV | not taken | Seed purity against pT_est |
| phi-cut tolerance M, 1-10 mrad | efficiency 0.9386-0.9393 over the whole range; fake 0.094-0.117 | 2 mrad chosen, the value does not matter | The doublet-slope phi cut |

Not built, sized from the truth tables only:
- choosing the best 4th hit per triplet would remove about 2 % of fakes (3.2 %
  of fakes are a true triplet with a wrong d);
- Matriplex for the per-triplet stages could buy at most ~15 % of the search
  time, the share of those stages after the float port.

In the prototype study (`../mkfit-standalone-seedgeom/`): radial sub-bins
(n_r = 8) touched 4.9x fewer hits and made the search slower (226 -> 394 ms),
since the start table grows with n_r.

## Acceptance: the quad list, not the physics numbers

The reference is the prototype in cover mode (`sg_cover(true)`,
`sg_dump()`): every bin range there provably covers its per-hit cut, so the
list does not depend on the grid. Checked on grids 64x8 to 1024x512 with 1-4
radial sub-bins. 5 events of the April HLT sample give **4115 quads**. Every
variant and setting below reproduces that list exactly, compared as sorted
(event, ia, ib, ic, id) tuples. The doublet and triplet counts also agree to
the unit: 835863 and 47282 per event.

**Write each dump to a new file name.** A rejected command-line option leaves
the previous dump in place, and comparing it again reads as a pass. That
happened twice before the driver runs were changed.

## The difference tool: `--margins REF`

Once the arithmetic changes, "the same list" can no longer be the test. The
test becomes: **every quad that differs is within epsilon of some cut in the
reference arithmetic.** Per event, the tool takes the symmetric difference
with REF and evaluates each differing quad in both arithmetics. The cut whose
pass flag differs is the one that flipped. Each quad is classed by the
reference margin of the closest flipped cut:

- `edge`: within epsilon. That is rounding.
- `FAR`: beyond epsilon. That is a bug, or a precision defect.
- `NOFLIP`: no cut flips, so the evaluator and the finder disagree, or a fetch
  missed a quad its own cut accepts.

It also prints two things that make the verdict meaningful:

- **the baseline**: how many reference quads have any cut within epsilon at
  all. It is 67 of 41226 (0.16 %). So a difference that lands at an edge is
  not a coincidence.
- **the precision table**: `|margin_ref - margin_new|` per cut over every
  reference quad, with the worst quad's pT and cot. This measures how far the
  new arithmetic moves each cut, whether or not any quad flipped. It is the
  number to read first.

Two self-checks are built in. The evaluator in the run's own arithmetic must
accept every quad that run found. `--arith ref` must give zero differences.
Both hold.

`test-seedgeom/quads/A-ref-50ev.txt` is the 50-event reference, made with
`--arith ref`. Its first 5 events equal the prototype's list exactly, and the
arithmetic is the prototype's, so it extends the reference rather than
replacing it.

## Measurements

Single thread, Ryzen 7 2700 (Zen+), minimum of 3 repetitions per event, 5
events, layers 0,1,2 + 3, the prototype's default window (`sg_philin` off).
"ns/doublet" is the search time over the doublet count. The fill of the four
layers is 1.2-1.6 ms per event on top.

| variant | ms/ev | ns/doublet |
|---|---|---|
| prototype (ACLiC, own CSR grid, 256x32) | -- | 176.6 |
| scalar port, q bins 2.0 cm | 140 | 167.9 |
| scalar port, q bins 0.5 cm | 114 | 136.0 |
| staged, 0.5 cm | 104 | 124.7 |
| fused 3+4, 1.0 cm | 86 | 102.8 |
| b-major, binary search | 75 | 89.4 |
| b-major, z buckets 0.25 cm | 53.5 | 64.0 |
| + x, y precomputed at fill | 51.8 | 61.9 |
| same, `-march=native` | 50.5 | 60.4 |
| *step A, re-timed at load ~4:* | | |
| templated, `--arith ref` | 54.3 | 65.0 |
| `--arith fast` (float, vdt, circle centre) | 43.6 | 52.2 |
| **`--arith fastk`** (float, vdt, curvature) | **42.8** | **51.2** |

| *step B, float kernels, re-timed at load ~2-3:* | | |
| b-major, fastk (same run) | 44.2 | 52.9 |
| `--tile-brute` (dense slope tile, 8-wide AVX) | 37.8 | 45.2 |
| **`--tile`** (slope buckets, one masked 8-wide step per doublet) | **24.6** | **29.5** |
| same, `-march=native` (AVX2, FMA) | 24.4 | 29.1 |

The step-A rows were run back to back, twice, on a box at load average ~4;
the two passes agree to 0.3 %. Read them against each other, not against the
rows above, which were taken on a quieter box. `fastk` against `ref`: helix
stage 13.0 -> 5.2 ms per event, 4th-layer test 4.4 -> 2.0.

Last row by stage, ms per event: doublets 3.5, per-doublet window 6.1,
per-b-hit lists 8.6, lookup and z test 16.4, helix 11.3, 4th-layer candidates
2.0, 4th-layer test 4.2.

What the profile said at each step, since that is what drove the order:

- **Scalar port:** ~600 instructions and 5.5 branch mispredicts per doublet,
  with the data in L2 (LL miss rate 0.1 %). Cutting layer-c hits touched 2.4x
  bought only 15 %: the cost was per-doublet bookkeeping, not per hit.
- **Staged:** the candidate test cost 12 ns per candidate. The reason was the
  gather chain candidate -> doublet -> hit arrays. Fusing it into candidate
  generation removed that.
- **b-major:** the binary search in the per-b list cost ~48 ns per doublet,
  latency bound (~12 dependent loads). A branch-free version did not help.
  Buckets turned it into two O(1) lookups.
- **`-march=native`:** 3.5 %. Expected while the loops are scalar; AVX2 pays
  once they vectorise.

## Where it goes next: the plan (agreed with the maintainer, 2026-09-24)

**Bit-equality with the prototype ends here, deliberately.** Everything above
keeps the prototype's arithmetic exactly, doubles, `hypot`, `asin` and all,
because that made the quad list a complete check of every restructuring. A
vectorised seeder cannot keep it, so the acceptance test changes from "the
same list" to **"differences only at a cut edge, counted and bounded"**.

### Step A: float, -Ofast, vdt in the per-triplet stages

**Status 2026-09-24: done for float + vdt. Matriplex is not done.**

- The build was already `-Ofast`, so this step was float + vdt.
- **The naive port fails, and the difference tool is what showed it.** The
  circle-centre form in float (`--arith fast`) produces 3 FAR differences in
  50 events. Two of them are at pT ~2 TeV, where the centre sits ~10^5 cm away
  and `dC^2 - R^2` cancels. There `d_phi` moves by 1.5 mrad. The third is at
  pT 19 GeV: `d_phi` flipped 2.5 urad from its edge after the float path
  moved it by 6 urad. On 5 events the same port gave the identical list, so a
  5-event null would have passed it.
- **The curvature form fixes it** (`--arith fastk`). The circle is given by
  the signed Menger curvature k and by the point and tangent at the third hit.
  The crossing with |X| = r is the line X.(k P0 + n) = k (r^2 + r0^2)/2 + P0.n,
  which is the circle equation multiplied through by k. All its terms stay
  O(r) as k -> 0. The arc length is L asin(h)/h, exact at k = 0.
- Result on 50 events: **0 differing quads out of 41226.** Largest margin
  shift per cut: `d_phi` 0.48 urad, `d_z` 83 nm (at |cot| ~3), `reach` and
  `d_cross` 10 nm. All are under the 1 um / 1 urad epsilon, so any future flip
  must be an edge flip.
- The degeneracy cut `c3` of the reference (collinear points, |G| < 1e-12)
  has no counterpart in the curvature form, since k = 0 is simply a straight
  line.

Original plan text:

- The helix and 4th-layer stages (`finish_triplets()`) go to float. `circle3`,
  `circle_cross_r` and `arc` in float, `hypot` -> `sqrt`, `asin` and `atan2`
  from vdt (`vdt::fast_asinf`, `vdt::fast_atan2f`), which the build already
  carries in `mkFit-external/vdt`. That is 11.3 + 4.2 ms of the 51.8 per event
  today, at ~240 ns per triplet.
- **The difference tool first, before any change.** `seedfind --margins REF`:
  for every quad in the symmetric difference with a reference list, evaluate
  every cut in BOTH arithmetics and print the value, the tolerance and the
  distance to it. Acceptance: every differing quad has at least one cut within
  a stated epsilon of its threshold (say 1 um in z, 1 urad in phi). Report the
  count as well. A difference far from every cut is a bug, not rounding.
- The exact list stays the reference: `test-seedgeom/quads/cov-g256x32r1.txt`
  in the isolated build, 4115 quads over 5 events.
- Then Matriplex with cdt for the per-triplet stages. There is enough
  arithmetic per item there for it to pay (~47 k triplets per event).
  *After the float port those stages are 7.2 of 42.8 ms per event, so the
  most this can buy is ~15 %. Step B addresses the ~35 ms in stages 1-4.*

### Step B0: the numbers step B is designed around (done 2026-09-24)

`seedfind --stats`, 5 events, fastk. Its doublet and triplet counts equal the
finder's: 835863 and 47282 per event. Per b-hit, 7318 b-hits per event:

| quantity | mean | p50 | p99 | max |
|---|---|---|---|---|
| a-hits fetched | 165 | 163 | 232 | 290 |
| contiguous a runs | 1.01 | 1 | 2 | 2 |
| doublets (69 % of fetched) | 114 | 114 | 163 | 190 |
| c-hits fetched | 100 | 99 | 147 | 171 |
| **contiguous c runs** | **20.9** | 21 | 33 | 40 |
| c list (67 % of fetched) | 67 | 66 | 97 | 108 |
| dense tile, doublets x list | 7702 | 7446 | 13460 | 15390 |
| z-bucket candidates | 203 | 187 | 601 | 860 |

Per doublet: z-bucket candidates 1.78 (p99 11), triplets 0.057.

Five findings:

1. **Stages 1 and 3 are dense loops over contiguous hits against one
   broadcast b-hit.** That is the ideal SIMD shape. Stage 1 is one run of
   ~165 hits. Stage 3 is fragmented into 21 runs only because layer c is
   binned in q, while the list build asks for all of q. A phi-only copy of
   layer c makes it one run.
2. **The z test is exactly a slope test relative to b.** Dividing the
   determinant by dr_a dr_c > 0 gives |s_c - s_a| < q_win / dr_c, with
   s = dz/dr from the b-hit. The slope form counts 236417 triplets over 5
   events against 236412, and the 5 extra are edge rounding. Candidates within
   the list's LARGEST tolerance number **0.058 per doublet, against 1.78 in z
   buckets**, so the slope pre-filter is 98 % pure. The layer's radial spread,
   which widens every z window, drops out.
3. **So the dense tile is the right shape after all.** Its 38x more pair
   tests than the buckets are each a subtract, abs and compare, 16 per AVX2
   instruction in int16. The current scalar stage 4 costs ~11 ns per bucket
   candidate. Tile passes are rare (0.058 per doublet), so the movemask is
   almost always zero.
4. **No multiply is needed in the default mode.** Every per-pair test becomes
   subtract, abs, compare against a per-b constant plus a per-hit term. The
   phi band is w = [r_b/2R - D0/r_b + marg] + [D0/r_a - r_a/2R], per-b plus
   per-a, and the same for c. The pmaddwd determinant is needed only for the
   linear-phi mode.
5. **The bit budget, from the data.**
   - Slopes: |s_a| <= 18.5 and |s_c| <= 11.3. int16 over +-16 with saturation
     is conservative, since a saturated pair can only pass the pre-filter, not
     fail it. The LSB is 4.9e-4, 7.5 % of the tightest tolerance (0.0065).
   - Phi: |dphi| <= 0.069 rad fetched. The uint32 difference >> 16 gives
     int16 with a 96 urad LSB over +-pi and no saturation, ~0.5 % of a typical
     window.
   - Both are pre-filter precisions. Survivors get the exact step-A test, so
     the output list stays identical.

Profile of stages 1-4 as they stand (cachegrind, 1 event): ~265 instructions
and ~1.4 branch mispredicts per doublet. Stage 4 is 76 of them and 0.68 of
the mispredicts. The bucket lambdas take another 37, the per-doublet window
of stage 2 51, and `wrap_pi` ~30. Stage 2 disappears in the slope form except
for one division per doublet.

**Step B as it stands after B0**, per b-hit:

- **K1**: the a-window, one run, test phi (and r order), then s_a.
- **K3**: the c-window from a phi-only copy of layer c, test phi (and r
  order), then s_c and t_c = q_win/dr_c plus the rounding slack.
- **K4**: the tile |s_a - s_c| <= t_c. Its rare hits get the exact c_z test
  and go to `finish_triplets<ArithFastK>`.

Order: first in float with 8 lanes, which the current `-mavx` build already
has for float, to validate the shape with an identical list. Then int16 with
16 lanes, which needs AVX2. On NEON the one missing piece is movemask
(emulated with a narrowing shift).

### Step B, first half: the float kernels (done 2026-09-24)

`--tile`. On 50 events it gives the identical list of 41226 quads and the
same doublet and triplet counts. The pre-filter passes 0.0585 candidates per
doublet into the exact z test. Per event, min of 5 reps:

| stage | b-major | tile |
|---|---|---|
| doublets (stage 1 / K1) | 4.0 ms | 4.6 |
| c lists (stage 3 / K3) | 8.2 | 5.4 (with the counting sort) |
| z test (stages 2 + 4 / K4) | 23.7 | 5.7 |
| helix + 4th layer (5-7) | 9.5 | 9.6 |

Branch misses fell from ~1.4 to ~0.3 per doublet (perf).

Three things that had to be found by measuring:

- **The dense tile was volume-bound** at ~9 cycles per 8 pairs. Slope
  buckets reduce it to one masked step per doublet: 19 -> 5.7 ms.
- **A count captured by reference into a lambda is kept in memory.** The
  compaction loops reloaded it on every element (`movl -0x124(%rbp)` hot in
  perf annotate). Explicit run loops with a local count fixed it: K1
  7.3 -> 4.6 ms.
- **AVX2 buys ~1 % on this Zen+ box.** With `-march=native` GCC still
  vectorises at 16 bytes. Its tuning for Zen 1/Zen+ treats 128-bit vectors as
  optimal, since the units are 128 bits wide. Time on a Zen 2+ or Intel box
  before judging AVX2.

**The same on uaf-4** (AMD EPYC 7662, Zen 2, which has full 256-bit vector
units), 2026-09-24:

- Setup under `/ceph/users/matevz/seeding-bench/`: `git archive` of the same
  two commits, GCC 15.3.1 and TBB 2022.3 from CVMFS `el9_amd64_gcc15`.
- The 50-event `--arith ref` reference list is byte-identical to black's (same
  md5). `--tile` with `-mavx` and with `-march=native` both reproduce it: 0
  differing quads of 41226.
- Timings, min of 5 reps, two passes within 0.5 %, load ~0.5:

| ns / doublet | `-mavx` | `-march=native` (znver2) |
|---|---|---|
| b-major fastk | 46.7 | 45.3 |
| tile brute | 41.4 | 40.9 |
| **tile** | **25.8** | **24.4** |

  With `-march=native`, tile per event is K1 3.5, K3 4.2, K4 4.7, helix 4.6,
  4th layer 3.7 ms. With `-mavx`: 4.1 / 4.6 / 4.7 / 4.7 / 3.5.

- **AVX2 is worth 5-6 % on Zen 2.** GCC vectorises K1 and K3 at 32 bytes
  there, not 16 as on Zen+, and K1 gains 15 %. K4 does not move, since its
  intrinsics are 256-bit AVX in both builds. The rest of the tile finder is
  still scalar, so this is the ceiling of auto-vectorisation, not of AVX2.
- Per core, uaf-4 runs the tile finder 13-17 % faster than black (25.8 against
  29.5 ns with the same flags).

### After the float kernels: bookkeeping, then vectorising stages 5 and 7 (2026-09-24)

perf on the tile finder put the test arithmetic of K1 and K3 at only ~7 % of
the samples. That is what int16 would speed up, so **int16 was deferred**.
Three changes went in first, each checked against the identical 50-event
list:

- `8162325`: K1 drops the constant b-hit stream, K3 sums prefixes over the
  occupied slope buckets only, and K4 runs in two phases. 30.0 -> 28.1 ns per
  doublet.
- `e1ec231`: the curvature-form helix is branch-free (`helix_k`) and runs in
  a simd loop over struct-of-arrays. Three obstacles had to go: calls that
  were not inlined (fixed with `gnu::flatten`), vdt's union in `fast_asinf`
  (replaced by `vdtv::fast_asinf`, which gives the same values through
  `copysign`), and a store through a member's address. Helix 5.6 -> 2.4 ms per
  event.
- `4b66a05`: the 4th-layer test is done the same way (`quad_k`). 2.0 -> 0.93 ms
  per event.

The margin evaluator calls the same `helix_k` / `quad_k`. Its per-cut shift
table is unchanged to every digit, so the vectorised stages reproduce the
scalar float results.

| ns / doublet | black, `-mavx` | uaf-4, `-mavx` | uaf-4, `-march=native` |
|---|---|---|---|
| tile, `d4dd6f3` | 29.5 | 25.8 | 24.4 |
| **tile, `4b66a05`** | **23.2** | **21.1** | **18.75** |

Against the prototype's 176.6 ns per doublet, that is **7.6x on black and
9.4x on uaf-4 with AVX2**, for the same quad list.

uaf-4 native per event: K1 3.3, K3 3.9, **K4 5.1**, helix 1.2, 4th-layer
candidates 1.7, 4th-layer test 0.6 ms. K4 is now a third of the time.

**Four machines, 2026-09-24.** `/ceph` is shared between uaf-4, uaf-9 and
phi3, so all three ran the same sources: `git archive` of `59b10cfe7db` and
`mkfit-seeding-int16` at `bb0b0d0`, and the same sample. `fastk` gave the
identical 50-event list everywhere. `fastk16` (branch `mkfit-seeding-int16`)
gave the same 128 edge-flip quads everywhere. Min of 5 reps, 5 events, ns per
doublet:

| machine | compiler | flags | fastk | fastk16 |
|---|---|---|---|---|
| black, Ryzen 7 2700 (Zen+) | GCC 15.2 | `-mavx` | 23.2 | 22.6-23.2 |
| uaf-9, Xeon E5-2670 v3 (Haswell) | GCC 15.3 | `-mavx` | 29.5 | 28.0 |
| uaf-9 | GCC 15.3 | `-march=haswell` | 25.7 | 24.6 |
| uaf-4, EPYC 7662 (Zen 2) | GCC 15.3 | `-march=native` | 18.7 | 18.05 |
| phi3, Xeon Gold 6130 (Skylake-SP) | GCC 14.3 | `-march=skylake` (AVX2) | 17.9 | 17.3 |
| phi3 | GCC 14.3 | `-march=skylake-avx512` | 17.1 | 16.9 |
| phi3 | GCC 14.3 | `... -mprefer-vector-width=512` | 17.5 | 17.1 |

- AVX2 over `-mavx` is worth 13 % on Haswell, 5-6 % on Zen 2 and ~1 % on
  Zen+.
- AVX-512 at GCC's default 256-bit preference is worth another 4.5 % on
  Skylake-SP.
- Forcing 512-bit vectors speeds up the helix by 18 % (1.05 -> 0.86 ms per
  event) but slows every other stage by 4-5 %, so it is a net loss. This is
  consistent with the known 512-bit clock reduction on that part; the
  frequency was not measured.
- int16 is -1.5 to -5 % everywhere, including at 32 lanes, because the test
  arithmetic is a small share of the time.
- phi3 has GCC 14.3 from CVMFS `el8_amd64_gcc14` (Alma 8 has no gcc15
  build), so its rows carry a compiler difference too.

**Reference timing set (2026-09-24 22:40-22:45).** One session, all four
machines, the first **20** events (837323 doublets and 833 quads per event),
fastest of 5 reps per event, fastest of 3 interleaved passes; search time per
event [ms], single thread. The tables above are on 5 events and at various
loads; these supersede them for quoting.

| machine, build | load | fastk | fastk16 |
|---|---|---|---|
| black, `-mavx` | 1.3-1.7 | 19.76 | 19.14 |
| black, `-march=native` | 1.3-1.7 | 19.80 | -- |
| uaf-9, `-mavx` | 0.0-0.4 | 25.43 | 23.91 |
| uaf-9, `-march=haswell` | 0.2-0.5 | 21.98 | 20.90 |
| uaf-4, `-mavx` | 4.0-4.2 | 17.68 | 17.19 |
| uaf-4, `-march=native` | 4.0 | 15.89 | 15.43 |
| phi3, `-march=skylake` | 0.0-0.4 | 15.06 | 14.70 |
| phi3, `-march=skylake-avx512` | 0.1-0.4 | 14.59 | 14.30 |
| phi3, `... -mprefer-vector-width=512` | 0.1-0.4 | 14.96 | 14.38 |

black, every version, `-mavx`: scalar port q 2 cm 142.30, q 0.5 cm 120.11,
staged 89.66, fused 81.05, b-major ref 53.08, fast 44.13, fastk 39.31, tile
brute 27.65, tile 19.76 ms/ev. The prototype was timed on 5 events only
(176.6 ns of search time per doublet); x 837323 that is ~148 ms/ev, an
estimate. Hence 7.5x on black and 10.3x at best (phi3, AVX-512, int16).
AVX2 over `-mavx`: Haswell 14 %, Zen 2 10 %, Zen+ 0 %. AVX-512 at 256-bit
preference over AVX2: 3 %.

**K4: three restructurings measured, none kept** (black, interleaved against
the lookup form at 5.5 ms per event; the list was identical in all three):

| variant | K4 ms / event |
|---|---|
| lookup form, one 8-wide masked step per doublet (kept) | 5.5 |
| masks stored at fixed per-doublet slots, to break a suspected loop-carried dependency | 5.9 |
| bucket against bucket: doublets counting-sorted too, list buckets loaded once | 8.4 |
| 4-wide SSE step, since a range holds ~3 entries | 6.4 |
| *lookup form with the exact confirmation skipped, as a cost split* | *5.0* |

The exact confirmation costs only 0.45 ms. The rest is the pre-filter step,
~19 cycles per doublet at IPC ~1. The join removed the two table lookups
per doublet and added a second counting sort, and it was slower. So the
lookups are not what K4 pays for. The hot store that suggested a dependency
was most likely perf sampling skid.

### Seed quality from truth: `seedfind --truth` (2026-09-24)

`SeedTruth.h` matches the finder's quads to simulation with the study's
definitions (`../mkfit-standalone-seedgeom/SeedGeom.cc`): a hit's label is
`mcHitID -> simHitsInfo_ -> mcTrackID`; findable is a hit with the label in all
four layers, pT > pT_min and the production point within D0_max of the beam
spot in xy; found is one quad of four hits with the label; a quad is true,
fake (four valid labels that disagree) or undecidable (a label is -1). New
options `--pt-min`, `--d0-max`, `--truth OUT.txt`. `truth2root.py` turns the
tables into TGraphAsymmErrors and TMultiGraphs versus |eta| for JSROOT.

At pT_min 0.9 it reproduces the study: efficiency **0.9395** on the same 20
events (the study: 0.9395), decidable purity 0.910. The tile finder was also
checked with the difference tool against `--bmajor --arith ref` at the other
windows, 5 events: 0 differences at pT_min 0.5 and 2.0. At 0.2 one quad is
extra and 2 of 12903 are rejected by the finder's own evaluator, by 3e-8 rad
(d_phi) and 15 nm (d_z): the vectorised helix/quad loops and the scalar
evaluator differ in the last rounding. So they are not bit-identical in
general, although the margin table at 0.9 is.

20 events, tile finder, fastk, D0_max 1 mm. Only the phi band follows pT_min;
the z windows and the 4th-hit phi window stay at their 0.9 GeV values. Time:
black, fastest of 5 reps, fastest of 3 interleaved passes, load 0.4-1.4:

| pT_min [GeV] | findable /ev | efficiency | fake / decidable | undecidable | extra true quads per found track | doublets /ev | ms /ev |
|---|---|---|---|---|---|---|---|
| 0.2 | 1870 | 0.522 | 0.292 | 0.380 | 0.108 | 2.33 M | 111.8 |
| 0.5 | 940 | 0.795 | 0.119 | 0.211 | 0.107 | 1.18 M | 33.1 |
| 0.9 | 376 | 0.940 | 0.090 | 0.184 | 0.117 | 0.84 M | 20.0 |
| 2.0 | 72 | 0.984 | 0.091 | 0.192 | 0.150 | 0.60 M | 12.8 |

- **Duplicates are 0.11-0.15 per found track, flat in |eta|.** The study's
  "~1.8 true quads per found track" divided all true quads by found tracks.
  At pT_min 0.9, 218 of the 612 true quads per event belong to non-findable
  tracks: 203 below pT_min, 15 produced beyond D0_max.
- **pT_min is a lower bound, not a cut.** Only the curvature term of the phi
  band depends on it; the D0 term and the margin dominate at high pT_min. With
  the windows at 2 GeV, 171 true quads per event belong to tracks below 2 GeV,
  against 72 findable above.
- **Low-pT efficiency is the fixed windows.** The z windows and the 2 mrad at
  the 4th hit are sized for 0.9 GeV.
- The fake rate at pT_min 0.9 is flat at 6-9 % up to |eta| 1.0 and rises past
  the four-layer acceptance edge (1.12): 23 % at 1.2-1.3, > 80 % beyond 1.4.

**3 of 4 as a match**, same runs (columns 12-20 of the tables). The study's
undecidable class is lenient: a quad whose valid labels disagree cannot be
true whatever its unlabelled hits are, so the consistent rule (fake if max
label count + unlabelled < 4) is the fair 4-of-4 baseline:

| pT_min | eff 4/4 | eff 3/4 | fake 4/4 study rule | 4/4 consistent | 3/4 | dup 4/4 | 3/4 |
|---|---|---|---|---|---|---|---|
| 0.2 | 0.522 | 0.533 | 0.292 | 0.529 | 0.420 | 0.108 | 0.146 |
| 0.5 | 0.795 | 0.801 | 0.119 | 0.246 | 0.163 | 0.107 | 0.146 |
| 0.9 | 0.940 | 0.942 | 0.090 | 0.190 | 0.118 | 0.117 | 0.158 |
| 2.0 | 0.984 | 0.985 | 0.091 | 0.191 | 0.116 | 0.150 | 0.208 |

Over tracks with a labelled hit in >= 3 of the 4 layers (550 /ev at 0.9) the
3-of-4 efficiency is 0.70.

**Window scale.** `--qwin`, `--qwin-d`, `--phiwin-d` (defaults 0.035 cm, 0.025
cm, 0.002 rad), all three x f, 20 events, study-rule fakes:

| pT_min | f | eff | fake | quads /ev | ms /ev |
|---|---|---|---|---|---|
| 0.5 | 1 / 1.4 / 1.8 / 2.5 / 3.5 | 0.795 / 0.903 / 0.945 / 0.971 / 0.981 | 0.12 / 0.23 / 0.37 / 0.58 / 0.78 | 1390 / 2179 / 3289 / 6560 / 15320 | 33.7 / 41.0 / 45.7 / 58.0 / 70.5 |
| 0.9 | 1 / 1.4 / 1.8 / 2.5 | 0.940 / 0.973 / 0.982 / 0.990 | 0.09 / 0.20 / 0.33 / 0.55 | 833 / 1161 / 1642 / 3101 | 19.8 / 22.4 / 24.6 / 28.9 |
| 0.2 | 1 / 2.5 / 4.5 | 0.522 / 0.847 / 0.940 | 0.29 / 0.79 / 0.95 | 2657 / 23881 / 128087 | 113 / 221 / 378 |

f = 0.9/0.5 = 1.8 restores the 0.9 GeV efficiency at 0.5 GeV, at 5.4x the fake
quads. Times: fastest of 3 reps, load 2.4-3.3. The tile finder at pT_min 0.5,
f 1.8, against `--bmajor --arith ref`: 2 edge flips, self-check clean.

**Why a findable track is missed: `seedfind --why-missed`** (`SeedMissed.h`).
Every cut is evaluated on the missed track's own hits with `eval_quad<A>`, the
best combination kept when a layer holds several, and the first failing cut
in pipeline order is the reason. pT_min 0.9, 20 events: 22.7 of 375.5
findable tracks per event are missed (6.05 %). First failing cut: d_z 8.05 /ev
(median 52 um beyond 250), d_phi 7.25 (0.51 mrad beyond 2), c_z 6.80 (141 um
beyond 350), ab_phi 0.35, bc_phi 0.25; **no missed track passes every cut**,
so none is a finder or fetch loss. By pT: 15.1 / 6.5 / 1.1 / 0.05 per event
below 1.2 / 1.2-2 / 2-5 / above 5 GeV: multiple scattering at the lowest
momenta, typically 1.2-1.4x outside the fixed windows.

**Before and after on four machines, 2026-09-26.** db9ae41 (before) and 803b561
(after, K4 in its own function, the left-pack, the sort function) built on each
machine with the reference set's compiler and flags, 20 events, fastest of 5
reps per event and of 3 interleaved passes, one session per machine. Search time
per event [ms], before / after / after with `--phi-lin 2 0.002`:

| machine, build | before | after | + phi cut |
|---|---|---|---|
| black, Zen+, AVX (GCC 15.2) | 19.59 | 16.09 | 13.22 |
| uaf-4, Zen 2, AVX | 17.75 | 13.74 | 11.33 |
| uaf-4, Zen 2, native | 15.85 | 11.95 | 10.24 |
| uaf-9, Haswell, AVX | 25.38 | 20.64 | 17.71 |
| uaf-9, Haswell, AVX2 | 22.20 | 17.56 | 15.55 |
| phi3, Skylake-SP, AVX2 (GCC 14.3) | 14.98 | 12.14 | 10.25 |
| phi3, Skylake-SP, AVX-512 | 14.91 | 11.86 | 10.02 |

The 50-event quad lists of before and after are identical on every machine and
build. "Before" reproduces the reference set within 2.2 %. Raw output in the
working report's prep/machines-2026-09-26/, bench sources in
/ceph/users/matevz/seeding-bench/seedsrc/.

**Tried and not kept:** a two-pass fourth-layer fetch with `__builtin_prefetch`
of the layer-d start-table entries. Stage 6 went 1.95 -> 2.47 ms per event on
black, 2.02 without the prefetch. On black the whole run has 2.21 instructions
per cycle and 2.2 % of level-1 loads missing, so the finder is not memory bound.

**Vector left-pack in K1 and K3, 2026-09-26.** Both kernels compute their
pass mask vectorised, then compacted the passing entries with a scalar loop of
four or five stores per fetched hit, and Zen+ retires one store per cycle. The
compaction is now `detail::left_pack<NS>`: per 4 lanes one lookup of a pshufb
control, one byte shuffle and one 16-byte store per 32-bit stream. To make every
stream 32 bits wide, a doublet's two 16-bit bucket indices are one word,
`d_bb = blo | bhi << 16`. The 50-event list is identical in the default mode,
with the phi cut and in brute mode. Black, one session, fastest of 3 passes: K1
4.31 -> 2.91 ms per event, K3 4.91 -> 4.40, search 18.11 -> 16.06.

**K3's counting sort with ranks, 2026-09-26.** Per event the sort sees 7325
b-hits with lists of 68 entries on average (all between 32 and 127) spread over
19 of the 64 slope buckets, and it took ~350 ns per list. It is already a radix
sort with one 6-bit digit, so more digits would only add passes. In the profile
the samples sat on the scatter's running cursor per bucket. The histogram pass
now records each entry's rank within its bucket, and the scatter computes its
slot as bucket start + rank, from loads alone. Black, one session, fastest of
3 passes, load 2.4-2.5: K3 4.45 -> 3.94 ms per event, search 16.15 -> 15.56.
The 50-event list is identical, in the same order, by default and with the phi
cut.

Tried and not kept, same session and machine:
- Two, four or eight counter sets, entry j using set j mod S, to break the
  chain of increments on one counter: no gain with cursors (K3 4.33-4.70
  against 4.44), and a loss with ranks (4.23 / 4.30 against 3.93). The time is
  not the chained increment itself. Most likely the core held each cursor load
  until the older scattered stores had their addresses.
- Dropping the per-entry slope tolerance t_c, i.e. one stream less in K3 and
  the sort, with K4's pre-filter using the layer-wide t_max instead: the list is
  identical, K3 falls 0.25-0.4 ms, but 32 % more pre-filter hits (64330 against
  48870 per event) reach the exact test and K4 rises 0.35-0.8 ms.

**The doublet-slope phi cut in the tile finder, 2026-09-26.** `--phi-lin 2 M`
now works with `--tile`. It cuts the layer-c hit on its phi predicted linearly in
r from the doublet's own phi slope, |phi_c - phi_b - s_phi (r_c - r_b)| <=
D0_max |1/r_c - 1/r_b - (r_c - r_b)(1/r_b - 1/r_a)/(r_b - r_a)| + M, in addition
to the generic band. It runs in K4's confirm step next to the exact c_z test, in
the b-major finder's float expression: the 50-event lists of the two finders are
identical at M = 2 and 3 mrad. Mode 1 (linear only) is not supported, because the
c-list is fetched with the generic band.

pT_min 0.9, 20 events, strict 4-of-4, fakes by the consistent rule:

| M | efficiency | triplets / ev | quads / ev | fake | undecidable |
|---|---|---|---|---|---|
| off | 0.9395 | 48458 | 823.7 | 0.190 | 0.184 |
| 1 mrad | 0.9386 | 16689 | 729.6 | 0.094 | 0.127 |
| 2 mrad | 0.9389 | 17585 | 732.8 | 0.097 | 0.128 |
| 3 mrad | 0.9389 | 18447 | 735.6 | 0.101 | 0.130 |
| 5 mrad | 0.9390 | 20154 | 741.0 | 0.106 | 0.133 |
| 10 mrad | 0.9393 | 24215 | 752.3 | 0.117 | 0.140 |

Duplicates stay at 0.117 per found track. Timing on black, one session, fastest
of 5 reps over 3 interleaved passes: the search went 17.93 -> 15.14 ms per event
at 2 mrad. K4 costs 0.51 ms more, the confirm step now evaluating the phi cut on
49000 candidates per event. The circle and fourth-layer stages save 3.29 ms,
since 64 % fewer triplets reach them. The default stays off.

Spent on wider windows (all three fixed windows x f), the saving buys
efficiency. At pT_min 0.9 with the cut at 2 mrad: f = 1.4 gives efficiency 0.972
at fake rate 0.21 in 16.5 ms per event, against 0.940 / 0.19 / 17.9 today;
f = 1.8 gives 0.981 / 0.34 / 18.0. Without the cut f = 1.4 costs 20.6 ms at
fake rate 0.37. At pT_min 0.5 with windows x 1.8: efficiency 0.9449 -> 0.9430,
fake rate 0.594 -> 0.318, triplets per event 195703 -> 50158, search 43.9 ->
29.0 ms per event.

**K4's 19 cycles per doublet, explained, 2026-09-26.** The disassembly of
phase 1 (the per-doublet slope step) showed about 40 instructions per doublet,
with the loop counter kept on the stack (`addl $1,-0xc0(%rbp)` every iteration)
and the array pointers reloaded from `TileWork` every iteration. The byte store
to `h_m` may alias `TileWork`'s members, and the in-loop `StagedWork::fit` calls
of the long-range branch added register pressure. Phase 1 is now
`detail::k4_phase1<Brute>`, with every array a `__restrict` argument, the output
arrays sized once per b-hit, and the long-range steps in a cold helper. In
bucket mode the first step needs no lane mask, because the slots past the range
end hold higher slope buckets and cannot pass. The entries are stored only when
the mask is nonzero. Black, 20 events, fastest of 5 reps, 3 interleaved passes,
same session, K4 stage / whole search in ms per event:

| phase 1 form | K4 | search |
|---|---|---|
| before (inlined, spilled, stores every doublet) | 5.47 | 19.70 |
| own function, registers only, stores every doublet | 4.30 | 18.51 |
| ... and no lane mask (26 instructions per doublet) | 4.23 | 18.50 |
| **... and stores only for a nonzero mask** | **3.64** | **17.90** |
| diagnostic: the test kept, nothing stored (wrong results) | 2.32 | |

Removing 14 instructions bought ~1 %, so the loop was not instruction bound.
Storing at slot `nh` for every doublet made each store address wait for the
previous doublet's mask, the end of the longest chain in the loop. That cost
~0.6 ms, and the branch on the mask is cheaper. The quad list is identical on
50 events (41226 quads), in bucket and brute mode. At 3.2 GHz the phase-1 step
went from ~19 to ~12 cycles per doublet; the test alone is ~7 (two dependent
index loads, then two unaligned 8-wide loads). The reference timing set
predates this and was not re-measured.

**`--why-missed` extended, 2026-09-25.** The report now also gives the missed
fraction per pT bin with its own findable denominator, and for the three
fixed-window cuts the |residual| / window of the first failing cut (median,
90th percentile, maximum, fractions beyond 2x and 5x) per pT bin. It ends with
`effpt lo hi findable missed` lines in fine pT bins, 0.5-50 GeV, for plotting.
pT_min 0.9, windows x1, 20 events: missed 9.5 / 4.5 / 1.7 / 0.8 % at pT
0.9-1.2 / 1.2-2 / 2-5 / above 5 GeV (the last is 1 of 123 tracks). d_z failures
sit at a median 1.21x the window with 0.6 % beyond 5x; d_phi at 1.25x with 4 %
beyond 5x, up to 14x; c_z at 1.41x with 8 % beyond 5x, up to 30x. Tracks more
than 5x outside are 0.9 per event, 0.24 % of findable.

Scanning the three fixed windows by f = 1 / 1.4 / 1.8 / 2.5 at pT_min 0.5 and
0.9 reproduced the efficiencies of the earlier window scan exactly. The missed
fraction depends on f and pT only through f x pT: 12-13 % at f x pT ~ 1 GeV for
every combination, ~2 % at 2 GeV, ~1.2 % at 3 GeV. That is the signature of
multiple scattering against fixed windows. Above f x pT ~ 3.5 GeV it flattens
at 0.3-0.9 %, about 0.7 %, and widening does not recover it. At pT_min 0.9 and
f = 2.5 the 1.0 % left is 0.6 tracks per event failing the azimuth bands and
3.3 per event outside the widened windows by a median 1.5-1.9x of them. The
plots are slides 145-146 of the deck.

Also fixed: the `K4 pre-filter hits` counter summed over repetitions. It is
0.0573 per doublet (47859 per event against 47282 triplets, 98.8 % pure).

### Step B, second half: fixed-point integers (the original plan text)

The per-pair tests run 10^6 times per event and are pure geometry. Nothing in
them needs an exponent: every quantity is bounded by the detector. What they
need is resolution, and fixed point gives it directly.

- **phi as `uint32`, 2 pi = 2^32.** The LSB is 1.5 nrad, and a difference wraps
  by integer overflow, so every `wrap_pi` disappears.
- **Divide-free consistency tests, as integer determinants.** The layer-c z test
  `|z_c - z_a - cot (r_c - r_a)| < q_win` with `cot = (z_b - z_a)/(r_b - r_a)`
  becomes
  `|(z_c - z_a)(r_b - r_a) - (z_b - z_a)(r_c - r_a)| < q_win (r_b - r_a)`,
  and the linear phi extrapolation
  `|(phi_c - phi_b)(r_b - r_a) - (phi_b - phi_a)(r_c - r_b)| < w (r_b - r_a)`.
  This is the "equality of dz/dphi via dr" the whole study started from, done
  exactly.
- **16-bit window-local deltas, 32-bit products.** A 2x2 determinant
  `a d - b c` is one `pmaddwd` (`_mm256_madd_epi16`) on (a, -b) and (d, c):
  8 determinants per instruction from 16 int16 lanes.
- **The bit budget has to be written down per expression.** Absolute z at pixel
  precision needs ~20 bits, so coordinates are stored as int32 and only DELTAS
  are narrowed. Numbers for the pixel barrel, to be confirmed on the data:
  - z test: `|dz| <= 40 cm`. A 10 um LSB gives 40000, which does NOT fit int16;
    16 um gives 25000, which does, and is 4.6 % of the 350 um tolerance.
    `dr <= 15 cm` at 10 um is 15000. The product stays under 2^30.
  - phi test: `|dphi|` within a window is <= ~0.05 rad, so +-0.05 rad in int16
    is a 1.5 urad LSB. `dr` as above. The product stays under 2^30.
  - The outer tracker and the discs need their own per-layer-pair scales.
- **Where it stays float:** the per-triplet curvature and pT cut (degree 4-6 in
  the coordinates) and everything after the quad. Integers for the geometry,
  floats from the triplet on, the physics after the handover.
- Order: 32-bit scalar integers first, to validate the determinants against
  step A with the difference tool, then 16-bit window-local SIMD.

### Portability: AVX2 now, NEON in mind

- 256-bit integer SIMD needs **AVX2**. The build's `-mavx` gives only 128-bit
  integer operations. AVX2 dates from 2013 (Haswell), and the target for
  phase-2 HLT nodes is worth confirming. Measured so far, `-march=native` on
  this Zen+ box is 3.5 %, as expected while the code is scalar.
- **NEON (AArch64)** is 128-bit: 8 int16 lanes. The determinant maps onto
  `vmull_s16` + `vmlsl_s16` (widening multiply and multiply-subtract into
  int32x4), plus the `_high` forms for the upper half. So the int16-delta /
  int32-product DESIGN is portable. Only the instruction selection differs.
- Therefore: keep the algorithm in terms of 16-bit deltas and 32-bit products,
  and keep the intrinsics behind one thin layer, or first try whether the
  compiler vectorises plain loops over int16 arrays (with -Ofast and
  -fopenmp-simd it often does). Test on uaf-4 (a newer x86 box) as well as
  here, and on an ARM box when one is at hand.

### Later milestones (unchanged)

2. Pure-disc triplets and quads (TFPX): phi is linear in z for D0 = 0, so
   this is the easy forward case, and a check that predictions are written in
   the target layer's own (q, qbar).
3. Mixed barrel/disc combinations in the transition region.
4. Iterations and a larger D0_max, where starting further out is cheaper.

## Beyond the pixel barrel: `seedsurf` (2026-09-26/27)

`SeedSurf.h` is a scalar, double-precision reference finder for quads on any
four surfaces, barrel layers or discs, and `seedsurf.cc` is its driver
(`make BLD=... seedsurf`). It is written to be right, not fast; nothing of the
tile finder's arithmetic is in it. Every hit is predicted at its own qbar in the
target layer's (q, qbar): a barrel layer is parametrised by r and measures z, a
disc by z and measures r.

- b: the barrel finder's geometric phi bound, and the a-b line reaching r = 0
  within `--zv` (25 cm) of the beam spot. The bound is loose everywhere: true
  doublets sit 0.4-28 mrad inside it at q99.
- c: the a-b line in (qbar, q) and (qbar, phi). For D0 = 0 a point's azimuth is
  phi0 + s/2R and z is linear in s, so phi is exactly linear in z on a disc.
  Tolerance: D0_max times the departure of 1/r from linear in qbar, plus a margin.
- d: the circle through a, b, c in its curvature form (signed k, point and
  tangent at c), z linear in arc length. On a disc the arc runs to the hit's z;
  on a barrel layer Newton solves |P(s)| = r_hit.

Windows are set per pattern from true-quad residuals (`--resid`, analysed by
`prep/surf-resid.py` in the working report): fixed at the q99 of pT 0.9-1 GeV
for pixel-only patterns, and a + b/pT_est for patterns ending in OT1-P, fitted
as an envelope of the q95 per true-pT bin below 3 GeV. `--pattern-win` and
`--pattern-dwin` take them; `--win-scale` scales all.

**Which patterns, from the census** (`seedfind --sequences[-all]`, the pixel-layer
crossing order of every sim track, 100 events, pT > 0.9): B1 B2 B3 B4 below
|eta| 1.0, B1 B2 B3 F1 at 1.0-1.4, B1 B2 F1 F2 at 1.4-1.8, B1 F1 F2 F3 at
1.8-2.6, sliding Fk..Fk+3 at 2.6-3.8, F6 F7 F8 E1 and F7 F8 E1 E2 above. At
|eta| 1.0-1.4, 38 % of tracks with >= 3 pixel layers have exactly three (B1 B2 B3
or B1 B2 F1): the gap between B4's end (|eta| 1.12) and F1's inner edge (1.24).
90 % of them have an OT1-P hit, hence pixel triplet + OT1-P (0 1 2 4, 0 1 16 4,
0 16 17 4). The phase-2 HLT covers this population with its HighPtTripletStep.

**Overlaps are configurable.** A track leaves ~1.7 labelled hits per disc, so one
disc pattern holds ~8 true quads per track and sliding patterns find it 3-4 times.
- `--dedup N` drops a quad sharing >= N hits with a better kept one, over all
  patterns; better = fewer outer-tracker layers in the pattern, then the smaller
  sum of squared c and d residuals over the windows. Without the tier, fake OT1-P
  quads outranked true B1 B2 B3 B4 quads.
- `--own DELTA --own-skip K --own-skip-ot K`: a pattern keeps a quad only if its
  layers can be the first four the quad's a-d line crosses, each layer tested over
  its whole qbar slab; K definite layers may be skipped. Testing the mid-radius
  alone rejected true hits at the z edge of B3 and B4.

**Truth.** On the D121 PU200 sample (bestTkIdx arbitration fixed) true-quad
residuals have a q99 ten times the q95: labels on off-trajectory hits. `--bind
0.05` keeps a label only if the hit's SimHitState is within 500 um (2.6 cm in the
outer tracker); the residuals then match the April sample's.

**Result, D121 PU200, --bind 0.05, windows from events 0-39, measured on 40-99,
pT > 0.9, cleaning N = 3** (`prep/endcap-2026-09-26/cmp-D121-E2.txt`):

| configuration | found tracks / ev | quads / ev | fake among decidable |
|---|---|---|---|
| 11 pixel patterns (+ mirrors) | 1382 | 24.7k | 0.25 |
| + 3 OT1-P patterns, no ownership | 1452 | 54.5k | 0.65 |
| + OT1-P, own 0.2 skip 1, OT skip 0 | 1423 | 26.2k | 0.38 |
| + OT1-P, own 0.2 skip 1, OT skip 1 | 1442 | 41.3k | 0.61 |
| + OT1-P, own 0.2 skip 2, OT skip 0 | 1428 | 28.1k | 0.36 |

At |eta| 1.0-1.2 and 1.2-1.4 the OT1-P patterns take found tracks from 47.8 and
48.0 per event to 65.3 and 67.3 (last row). With OT skip 1 they also rescue
barrel tracks that lost a pixel hit, +19 tracks per event for +15k quads. The
OT1-P windows are wide at low pT (B1 B2 F1 + OT1-P: 1.8 mm + 3.6 mm GeV / pT in
z), from multiple scattering over the step; that is the open transition problem.
Timing is not measured: the finder is the scalar reference.

### The feed-forward chain (`SurfChain`, `seedsurf --chain`)

Maintainer's proposal: instead of a list of four-layer patterns, one pass over
the layers of each z side in crossing order (B1..B4, the discs by |z|, then
OT1-P). A doublet's layer sequence is determined by its own line: it is queued
at the next layer its line crosses, and after a hit at the next one after that.
A candidate that misses its target can be forwarded to the following crossed
layer with a hole (a hole is charged only if the target was a definite crossing).
Start pairs are derived at setup by scanning lines from the beam region, not
listed. The per-pattern window tables are looked up by layer combination; a
combination without one is not searched (`--chain-any` searches it with the
default windows, and that produced 2k fake quads per event in single
combinations). OT1-P is entered only without a charged hole (`--chain-holes-ot 0`).

Measured on D121 PU200, --bind 0.05, events 40-99, cleaning N = 3
(`prep/endcap-2026-09-26/cmp-chain*.txt`):

| | found / ev | quads / ev | fake among decidable |
|---|---|---|---|
| pattern list, own 0.2 skip 2 | 1428.3 | 28.1k | 0.36 |
| chain, no holes | 1403.6 | 22.1k | 0.48 |
| chain, holes when extending (1 or 2) | 1405.4 | 36.0k | 0.68 |
| chain, holes also in the start doublet | 1432.7 | 48.1k | 0.59 |
| **chain, start up to 2 crossed layers late, no other holes** | **1430.8** | **29.0k** | **0.38** |

Two readings. Holes during extension buy +2 tracks per event for +14k quads: a
fake doublet always fails at its third layer, so forwarded candidates are
almost all fakes getting a second chance. The recoverable tracks are lost at
the START (no hit on the first or second crossed layer); a late start
(`--chain-start-holes 2 --chain-lead-only`) recovers them for 8 % more doublets,
and it needs neither the ownership test nor the sliding pattern list.

**Decision (maintainer, 2026-09-27): the chain is the method from here.** The
late-start chain finds as many tracks as the best pattern list, at the same quad
count and fake fraction, with one mechanism: no ownership test, and no list of
sliding patterns, since each doublet's own line decides its layer sequence. The
pattern options stay, for two jobs: they carry the window tables (the chain
looks its windows up by layer combination), and they define the truth
denominator. Ownership and the plain pattern loop are kept as the comparison and
are not developed further.

`seedsurf-chain.sh` runs that configuration: `--chain 0 --chain-start-holes 2
--chain-lead-only --dedup 3`, `--pt-min 0.9 --marg-b 0.001 --bind 0.05`, and the
43 window lines in `windows-D121/`. These are the options of the measured
chain row, with the tables copied from the working report's
`prep/endcap-2026-09-26/patterns-D121-{px,ot-q0.95,skip}.txt`. Rerunning the
script on events 40-99, with the OT1-P windows as they were then, gave a
truth report identical to the recorded one (`truth-CHL2.txt`) in every line
except the timings. The OT1-P windows have since been refitted (next section).

**Chain counters for that run, per event:** 4.64 M doublets, 236 k triplets,
82 k quads before cleaning, 29.0 k after. The scalar double-precision
reference takes 2.07 s per event, which says nothing about a real
implementation.

### The multiple-scattering term of the d windows (2026-09-27)

The d windows are `a + b / max(pT_est, pt_min)`, and the b term is multiple
scattering. Scattering itself cannot be reduced; what can be got right is the
variable b is fitted in and the population it is fitted on. Residuals: D121
PU200, events 0-39, `--bind 0.05`, pT down to 0.5 GeV; `prep/ms-eta.py` fits
`a + b / pT` to the q95 in true-pT bins, per pattern and per |eta| bin. Raw
output in the working report's `prep/ms-eta-2026-09-27/`.

**b does not scale as 1/p.** If it did, b cosh(eta) would be flat. It is not,
in any pattern:
- pixel patterns, phi: b is roughly flat in |eta|, 1.3-2.1 mrad GeV;
- pixel discs, r: b falls, but more slowly than 1/cosh(eta) (F1 F2 F3 F4:
  157 -> 49 um GeV while cosh(eta) goes 3.7 -> 8.2);
- barrel z: b rises (B1 B2 B3 B4: 170 -> 347 um GeV over |eta| 0-1.2), because
  the projection onto a barrel layer grows faster than the angle falls with p;
- pixel triplet + OT1-P: b rises steeply, e.g. B1 B2 F1 + OT1-P b_z 1985 ->
  8574 um GeV and b_phi 6.4 -> 19.8 mrad GeV from |eta| 1.05-1.31 to 1.83-2.09.

A per-|eta| (a, b) instead of one per pattern would shrink the window area at
95 % containment by 5-10 % in the pixel patterns and by 16 % and 20 % in
B1 B2 B3 + OT1-P and B1 B2 F1 + OT1-P. Not implemented.

**The OT1-P windows were fitted on tracks the chain does not build.** `--resid`
now reports, per combination, the definite crossings of other layers on its a-d
line before hit a (lead) and between a and d (inner). The chain, with no holes
after the start, builds a pattern only if inner = 0. For B1 B2 B3 + OT1-P that is
11 k of 83 k true combinations: the rest cross B4 on the way and the chain builds
B1 B2 B3 B4 for them. Refitted with `prep/surf-resid.py` (q95 envelope) on the
inner = 0 rows, the d windows move both ways, since each pattern is used only in
part of its |eta| range:

| pattern | before: a_phi, b_phi, a_z, b_z | inner = 0 |
|---|---|---|
| B1 B2 B3 + OT1-P | 2.2 mrad, 6.1 mrad GeV, 1.24 mm, 0.78 mm GeV | 2.9, 7.9, 1.35, 1.03 |
| B1 B2 F1 + OT1-P | 2.7, 12.0, 1.80, 3.64 | 2.9, 4.7, 1.71, 1.85 |
| B1 F1 F2 + OT1-P | 4.6, 12.2, 2.46, 5.98 | 1.9, 6.7, 1.86, 7.10 |

The c windows are kept (the refit of B1 B2 F1's c window came out at 6.8 mrad,
from a q99 in a sparse pT bin). `windows-D121/ot1p.txt` holds the refitted d
windows. The chain with them, events 40-99 (not the events the windows were
fitted on), `prep/ms-eta-2026-09-27/cmp-chainfit.txt`:

| OT1-P d windows | found / ev | kept quads / ev | fake among decidable |
|---|---|---|---|
| fitted on all combinations | 1430.8 | 28 975 | 0.378 |
| **fitted on inner = 0** | **1431.5** | **27 468** | **0.343** |

Found tracks are unchanged (+0.4 at |eta| 1.0-1.2, +0.2 at 1.2-1.4). The quads
fall in the transition, 1273 -> 917 per event at |eta| 1.0-1.2 and 1585 -> 1159
at 1.2-1.4, and are unchanged above 1.8. Any window table for the chain should
be fitted on the inner = 0 rows.

**The material scan (geo-stuff-9d, `/baz/matevz/root-dev/geo-stuff/seedscan.txt`).**
Straight lines from (0, 0, z0), z0 uniform in +-5 cm, phi uniform, 400 per |eta|
point. Each line is integrated from the third layer C (its envelope included) to
the OT1-P sensor mid-plane, or to ITLayer4 as the reference. M0 is the
path-integrated x/X0 and M2 = integral of (s_d - s)^2 dX. The scan was
cross-checked line by line against geo-stuff-7a's independent implementation.
At |eta| 2.0, TFPX1 -> OT1-P crosses TFPX2-4 and their service sheets unused, and
those carry 0.10 of an M0 of 0.31.

**Highland on the scan reproduces the |eta| shape of b_z.** The prediction is
sigma_perp = 0.0136 GeV / p sqrt(<M2>) (1 + 0.038 ln M0), with p = pT cosh(eta). On
a barrel target at fixed r, dz = sigma_perp cosh(eta) and r dphi = sigma_perp,
and q95 = 1.96 sigma. That gives b_z = 1.96 * 0.0136 sqrt(<M2>) H and b_phi =
b_z / (cosh(eta) r_d) per GeV of pT. `prep/ms-highland.py` compares them with b
fitted in 0.1 bins of |eta|; the output is `prep/ms-eta-2026-09-27/ms-highland.txt`:
- measured / predicted b_z is 1.2-1.7 and flat in |eta| for B1 B2 B3 + OT1-P
  (|eta| 0.1-1.6), B1 B2 F1 + OT1-P (1.0-2.0), B1 F1 F2 + OT1-P (1.3-2.4), and
  the pixel-only B1 B2 B3 B4 against IT3 -> IT4. Over those ranges the prediction
  itself grows by a factor of 3-6. So the material explains the rise, up to a
  factor of ~1.45 that is the same everywhere. The low-statistics edges
  (|eta| > 2 for B1 B2 F1, n < 700) run to 2-2.5.
- measured / predicted b_phi is 2-3.5. The azimuthal residual also carries the
  triplet's curvature error from scattering before C, which the step from C to
  the target does not include.

**phi: the material has a 40 degree structure, and it shows in one pattern.**
The scan's M0 at |eta| 0.2 on IT3 -> OT1-P varies by a factor 2.2 in phi. Every
significant harmonic is a multiple of 9: a 40 degree period with a sharp comb,
strongest at 13.33 degrees (n = 27, 20 % of the mean) and n = 9 at 12 %, phase
+10.2 degrees. A first check in 24 absolute bins of 15 degrees saw nothing, and
that was aliasing: the bin width nearly equals the 13.33 degree period. `--resid`
rows now carry the d hit's and the track's azimuth. `prep/ms-phifold.py` folds
residual / (a + b/pT) modulo the period in 9 bins, pT < 1.5, and compares the
first harmonic of the per-bin q95 with 200 azimuth shuffles. The azimuth where
the track crosses radius r is interpolated between the vertex momentum phi and
the d hit's phi. Output: `prep/ms-eta-2026-09-27/ms-phifold.txt`.
- B1 B2 B3 + OT1-P, |eta| < 1, modulo 40: amplitude 0.075-0.093 in z and
  0.047-0.095 in phi (depending on the radius phi is taken at), against a
  shuffled q95 of 0.04-0.046. The phase is +9 to +11 degrees, matching the
  scan's +10.2.
- modulo 13.33: nothing above the shuffles, plausibly because the 1-3 degree
  uncertainty in where a curving pT < 1.5 track crosses the material washes the
  comb out.
- B1 B2 B3 B4 and the forward OT1-P patterns: nothing. So the source is not in
  IT3 or IT4. It lies between the IT4 envelope (r 15.6) and the OT1-P sensor.
  geo-stuff-9d has been asked to localize it in r.

**The source is the OT1 flat rods, and it is the lever arm, not material.**
geo-stuff-7a pointed out that the OT1 flat rods alternate between r = 22.11 and
24.49 cm every 20 degrees. `--resid` rows now carry the d hit's radius. For
B1 B2 B3 + OT1-P at |eta| < 0.5 the hits split into inner rods (median r 22.30)
and outer rods (24.55), and phi mod 40 separates the two classes. Within each
class, the q95 of residual / window over 5 degree bins has max/min 1.03-1.16,
against a shuffled q95 of 1.16-1.25: the phi structure vanishes. Between the
classes the ratio steps from 0.90 to 1.06-1.10. The window fitted per class at
pT 1 is 1.21x larger on the outer rods, and the lever-arm ratio from IT3 is
(24.55 - 10.52) / (22.30 - 10.52) = 1.19. geo-stuff-9d's radial slices agree:
M0 has no 9-fold structure below r = 21.5. The structure in M2 at smaller radii
comes only from the (s_d - s)^2 weight, i.e. from where the endpoint sits
(`prep/ms-eta-2026-09-27/ms-rodclass.txt`).

**Scaling b with each candidate's own lever arm: built, measured, a null, off by
default.** `--pattern-sref A B C D S` scales a listed pattern's b terms by s_cd /
S. s_cd is the 3D path length on the helix from hit c to the candidate, which the
prediction already solves for. The fetch uses the longer end of the layer's
slab. `--resid` rows carry s_cd. `prep/surf-resid-lever.py` fits a + b (s / s_ref)
/ pT with the same envelope as `surf-resid.py`, s_ref being the median s_cd;
without `--lever` it reproduces `windows-D121/ot1p.txt` to the last digit.
Chain, events 40-99 (`prep/ms-eta-2026-09-27/cmp-lever.txt`):

| OT1-P d windows | found / ev | kept quads / ev | fake among decidable |
|---|---|---|---|
| a + b / pT (committed) | 1431.5 | 27 468 | 0.343 |
| a + b (s / s_ref) / pT | 1431.5 | 27 451 | 0.343 |

Found tracks are identical in every |eta| bin. The quads move from short paths
(-100 to -134 per event at |eta| 0.8-1.2) to long ones (+110 to +154 at 1.4-1.8).
They do not fall. Without `--pattern-sref` the chain's truth report is identical
to before in every line but the timings.

**Per-|eta| d windows: built, measured, a null, off by default.**
`--pattern-dwin-eta A B C D LO HI APHI BPHI AQ BQ` adds an |eta| slice to a listed
pattern's d windows. The triplet's helix |eta| selects the slice; outside every
slice the pattern's own windows apply. `prep/surf-resid-eta.py` fits four slices
of equal count per pattern, on the inner = 0 rows, with the surf-resid envelope.
A slice needs at least 3 pT bins, otherwise the pattern keeps its single window,
which is what happens for B1 F1 F2 + OT1-P (911 tracks). B1 B2 B3 + OT1-P comes out
monotonic in |eta| (b_z 0.058 / 0.087 / 0.109 / 0.157 cm GeV). B1 B2 F1 + OT1-P
has one slice at 0.34 between neighbours at 0.08 and 0.20: an envelope that takes
the max over pT bins of ~200 tracks is noise-limited per slice. Chain, events
40-99 (`prep/ms-eta-2026-09-27/cmp-eta.txt`):

| OT1-P d windows | found / ev | kept quads / ev | fake among decidable |
|---|---|---|---|
| one a, b per pattern (committed) | 1431.5 | 27 468 | 0.343 |
| per-\|eta\| slices | 1431.5 | 27 705 | 0.349 |

Quads fall 4-5 % at |eta| 0.6-1.0 (the clean B1 B2 B3 slices) and rise at 1.2-1.8
(the noisy B1 B2 F1 slice). No bin moves in found tracks. The 5-20 % area gain
estimated at fixed containment above does not survive the per-slice fit
statistics of 40 events. Together with the lever-arm null, this says the d
windows fitted on the chain's rows are as good as a window table gets here. The
scattering term is left as a + b / pT, with pt_min as the knob.

### Making the chain fast (2026-09-27)

**Timing rule for this section.** `SurfChain::run`'s own time (the `CHAIN ...
ms/ev` line: the search only, no fill, no cleaning, no truth), D121 PU200
events 40-59, single thread, black, the chosen configuration without `--bind`;
the fastest of 3 interleaved passes. Raw output in the working report's
`prep/chain-kernels-2026-09-27/`.

**Acceptance: `seedsurf --margins REF`.** REF is a `--dump` of a run with the
same pattern options and events. Per event the tool compares the chain's quads
before cleaning with REF, as (pattern, hits). Each quad in the symmetric
difference is evaluated in double, with the windows the chain uses for its
layer combination: b (r order, phi bound, z0), c (phi, q), d (phi, q). The tool
reports the smallest |margin| / window over the cuts, so a flip at a threshold
can be told from a bug. A quad with no cut near a threshold is labelled
structural: a crossing or queue decision. The reference list for this section
is `ref-quads.txt`, from the scalar chain at 3d10608: 1516437 quads in 20
events. Run against itself, the tool
reports 0 differences.

**Profile of the scalar chain** (perf, 6 events): `SurfChain::run` 34 % self,
`atan2` 15 % from inside it, and `hypot` 25 % from `main` (the cleaning's
residuals, outside the chain's timer). Every use of a hit recomputed its r and
phi from x, y in double.

**K1, bit-identical: 1866.9 -> 708.7 ms per event, x2.63, the same 1516437
quads.**
- `surf::P3` computes r and phi once at construction. `SurfLayer` caches them per
  hit in double, from the layer's float x, y, with the same expressions, so
  `surf_p3` returns identical numbers (-42 %).
- The chain keeps each candidate's r-z line, which the crossing test takes. The
  line is rebuilt only when a hit is added, not on every test.
- The chain resolves the window parameters once per chain-position combination,
  at its first run, into a table indexed by position. Before, each queued
  candidate did a `std::map` lookup and a `SurfParams` copy.
- The per-layer queues persist across events, and a candidate is copied only
  when it is forwarded.

- A queued candidate holds its hits as (chain position, index in the layer),
  not as four points. The points come back from the layer's cache, identical to
  storing them, and the candidate shrinks from ~200 to ~56 bytes. It is copied
  into a queue about 2.3 M times per event: 708.7 -> 661.7 ms (interleaved with
  the step below).

Per doublet that is 162 ns, against 17-23 ns for the barrel tile finder. After
K1 the profile is spread out: stage b about 15 % (building the P3, the phi bound
with two divisions, z0 with one), the crossing tests about 9 %, the binnor run
loop 4 %, vector construction 4 %. That is where float kernels on structure of
arrays come in (K2), and they change the arithmetic: accept them by
`--margins` against `ref-quads.txt`, not by list identity.

**K2a, stage b in float (`--chain-fast`): 661.7 -> 611.3 ms per event, 0 of
1516437 quads differ.** `surf_stage_b_fast()` uses the same fetch as
`surf_stage_b()` (now `surf_b_fetch()`, shared). Per candidate b hit, it applies
the same three cuts in float on the SeedLayer arrays, with no division: the phi
bound from the cached 1/r, and z0 in [zlo, zhi] multiplied through by dr = r_b -
r_a > 0. It runs in masked chunks of 64 over the binnor's contiguous runs. The
doublet count is identical. Forwarded and dropped candidates move by 0.1 and 0.6
per event, and none of those becomes a quad. Fastest of 3 interleaved passes:

| | ms / ev | vs scalar |
|---|---|---|
| scalar chain (3d10608) | 1867.2 | -- |
| K1 (P3 cache, line, parameter table, queues, slim candidate) | 661.7 | x2.82 |
| + K2a (stage b in float) | 611.3 | x3.05 |

Stage b was 15 % of the profile after K1 and the float loop takes 7.6 % off, so
what remains is per candidate: the crossing tests on each doublet's line
(before b for the holes, after it for the next target), the queueing, and stage
c for ~2.2 M doublets per event that mostly find nothing.

**Where the time went after K2a** (`--chain-phases`: the start doublets by wall
clock, stage c and d by rdtsc, 20 events): start doublets 262 ms, forward pass
367 ms, of which stage c 40 % over 1.43 M doublets and stage d 60 % over 216 k
triplets. Stage d cost ~1 us per triplet: each triplet made ~6 exact predictions
(the fetch's two edges and one per candidate, 4.2 candidates per triplet), each
a Newton solve with a sincos per step on a barrel target.

**K2c, stage d as pre-filter plus exact confirm (`--chain-fast`): 615.0 ->
519.8 ms per event, 0 of 1516437 quads differ.**
- `Helix::cross_r_cf()`: the barrel crossing in closed form. The helix circle
  against |P| = r by the radical line, with d^2 - rho^2 taken as |c|^2 + 2 (c.n)/k
  so it does not cancel for a stiff track, then the arc to each point.
  `--chain-fast-check` measures it against Newton over 3.1 M node predictions:
  max 9.2e-9 cm in z and 1.4e-7 rad in phi, and never succeeds where Newton fails.
- `surf_stage_d_fast()` predicts at the layer's two qbar edges (also the fetch)
  and its middle, and takes q, phi and the path length as quadratics in qbar.
  Over 4.3 M candidates the quadratic's worst error is 2.0 % of the window in q
  and 0.42 % in phi. Candidates are tested against it in float with the window
  x1.05. Survivors are confirmed with `surf_stage_d`'s own prediction and double
  cut, so the accepted quads are decided as in the reference.
- Without the confirm, the quadratic alone flipped 52 of 1.5 M quads (26 each
  way). All were in OT1-P patterns, at the d q cut, within 5.4e-3 of the window:
  its interpolation error over the 8.4 cm OT1 slab.

| | ms / ev | vs scalar |
|---|---|---|
| scalar chain (3d10608) | 1884.7 | -- |
| + K1, K2a | 615.0 | x3.06 |
| + K2c | 519.8 | x3.63 |

(Fastest of 3 interleaved passes; the load rose from 1.75 to 3.68 during them.)

**This is still the scalar double prototype**, one candidate at a time: only the
stage b and stage d candidate tests run in float, and every decision that
reaches a quad is made in double. The profile after K2c is flat. At the top are
the crossing tests on each doublet's line (~11 %), building P3s from five
arrays (5 %) and the line construction (2 %). What remains is structural.

Next, in order (agreed 2026-09-27):
1. **Done: see "The batched float finder" below.** A batched float finder in its own file, next to `SurfChain`, as
   `SeedFinderTile.h` sits next to the scalar barrel port. It keeps the same
   chain logic, with candidates as structure-of-arrays batches per target
   layer. `SurfChain` stays, in double, as the reference. Accepted by `--margins`
   against `ref-quads.txt`.
2. **Measured and not taken: see "Stage d: one point and a direct solve" below.** Its stage d on MkFitCore's `mini_propagators` (maintainer's suggestion).
   - From NN triplets at a time, set up an `InitialStatePlex` at hit c (position,
     direction from the helix tangent, k).
   - `propagate_to_r` / `propagate_to_z` to the target layer's two edges. This is
     the same closed form as `cross_r_cf`, which then goes.
   - Build a two-point `Hermite3D` across the slab, and test each candidate
     against the cubic, or with `Hermite3DOnPlane` on its own module plane (the
     right comparison for the tilted OT1-P modules).
   - Expected error: R_c dalpha^4 / 384, ~0.3 um across OT1 at 0.9 GeV
     (estimated), against up to ~10 um for the quadratic. Measure it with
     `--chain-fast-check`, then decide whether the exact confirm can go.
   - Guard the degenerate span near the turning radius (m_Hderfac -> 0 when both
     edges clamp) with the fallback `surf_stage_d_fast` already has.
   - The mini-propagators use uniform B, as the helix does.
   - The one-point mode is 16x worse and not enough for the 15-30 cm c->d steps.
3. **Done in the batched finder.** Crossing tests vectorized over the batch. Add the crossing margin (the
   line's distance to the layer envelope edge) to `--margins` first, so a flipped
   crossing is classified instead of labelled structural.
4. Stage c with the fetch shared per b hit (b-major, as the barrel's K3/K4).
5. Iterations and larger D0.

### The batched float finder (2026-09-27)

`SeedSurfBatch.h`, run with `seedsurf --chain-batch`. `SurfChainBatch` takes a
set-up `SurfChain` for everything that is configuration: the chain order, the
start pairs, the window tables and the crossing envelopes. `SurfChain` stays the
double reference. What the batch finder changes:
- A queued candidate is 28 bytes: three hit indices, three chain positions, its
  holes, the target's crossing state, and its r-z line (z0, cot) in float.
- Each stage takes the candidates of one target layer in blocks of 256. The
  crossing tests run over the block as lanes, one layer at a time. They cover
  the holes before b, the start layers' own states, and the next crossed layer
  after b or c. Each line is kept as (z0, cot) for a barrel crossing and as
  (1/cot, -z0/cot) for a disc crossing, so the lane loops have no division and
  no branch. The work arrays are float and int32, which GCC vectorizes with
  `-mavx`. The same loops over `unsigned char` arrays stayed scalar.
- Stage b is `surf_stage_b_fast`'s float loop. The side test and the doublet's
  line are computed inside it, and the survivors are compacted without a branch.
  A phi-major copy of the b layers, which makes each fetch one run instead of
  one per q bin, was tried and dropped: it touched 46.9 M hits per event against
  12.6 M.
- Stage c is `surf_stage_c` in float. The phi residual is taken relative to
  hit a, `wrap(wrap(phi - phi_a) - dphi_ab t)`, so no phi of order pi enters a
  difference at the 0.5 mrad scale of `phi_c`.
- Stage d predicts at each hit's own qbar from one point, the helix at hit c
  (curvature k, point c, tangent at c), with no double confirm:
  - barrel, |P| = r: the radical line of |P| = r and the helix circle,
    multiplied through by k: P.(k c + n) = k (r^2 + |c|^2) / 2 + c.n, with n the
    left normal of the tangent. Every term is O(1) for any k, and the helix
    centre, 1/k away, never appears. Of the two points, the one ahead of c with
    the shorter chord L is taken. The arc is L asin(|k| L / 2) / (|k| L / 2),
    and z = z_c + cot * arc.
  - disc, z = u: the arc is (u - z_c) / cot, and the point is c plus the chord
    along the tangent turned by half the turning angle.
  - The phi cut |dphi| < w is taken as dot > 0 and cross^2 < sin^2(w) |P|^2
    |h|^2, on the hit's own x, y, so there is no atan2 per hit.
- `--margins` now also reports a crossing margin. It is the smallest distance
  of any envelope crossing of the a-b or a-c line to a crossing-state
  threshold, over the 0.2 cm crossing margin. A flipped crossing is then
  labelled `crossing` instead of structural.

**Accepted: 6 of 1516437 quads differ from `ref-quads.txt`**, 4 only in the
reference and 2 only in the batch run. All 6 are at a cut threshold: 5 at d_q
within 5.5e-5 of the window, 1 at c_phi at 4.0e-6. The doublet, triplet,
forwarded and dropped counts are those of the reference. With truth (events
40-99, `--bind 0.05`) the union row is unchanged to the printed digit: 1431.5
found tracks per event, 27468 quads, fake 0.3433, for K2c and for the batch
finder.

**Stage d precision** (`--chain-fast-check`, 20 events, 18.4 M fetched hits):
the float single-point prediction against surf::Helix::predict in double gives
max |dq| 2.9e-5 cm and max |dphi| 3.0e-7 rad. That is at most 2.1e-4 of the q
window and 7.6e-5 of the phi window. There are 0 hits where only one of the two
predictions succeeds.

| | ms / ev | vs scalar |
|---|---|---|
| scalar chain (3d10608) | 1884.7 | -- |
| K2c, `--chain-fast` | 520.7 | x3.62 |
| batch, `--chain-batch` | 385.5 | x4.89 |

(Fastest of 3 interleaved passes, K2c and batch only, load 1.3-1.8;
`prep/chain-kernels-2026-09-27/ab-b1.txt`.) Phases of the batch run
(`--chain-phases`): start doublets 180 ms, forward pass 210 ms, of which stage c
71 % over 1.43 M doublets and stage d 29 % over 216 k triplets.

Where the start doublets go (6 events, cycle counters since removed):
- stage b runs and mask, 1.36 M runs of ~9 hits: ~300 Mcyc per event;
- the flush of 3.3 M lanes: lines, holes and states 63 Mcyc, compaction 37,
  routing 62, queueing 105 (a 28-byte store per candidate into queues larger
  than the caches).

The three barrel start pairs (B1 B2, B2 B3, B3 B4) cost 102 Mcyc per side and
are fetched once for each side, so fetching them once for both sides would save
~30 ms per event.

### Stage d: one point and a direct solve, not the two-point Hermite (2026-09-27)

Plan step 2 was to put the batch finder's stage d on MkFitCore's
`mini_propagators` with a two-point `Hermite3D` across the target slab. The
direct single-point solve above was measured against the two `Hermite3D` modes
instead, as options of the batch finder (`--chain-batch-d N`). Modes 1 and 2
reproduce `Hermite3D`'s cubics in scalar float, not its Matriplex code:
- 0 (default): direct, from the helix at hit c, as described above;
- 1: `Hermite3D`'s one-point mode, the Taylor cubic about c in the transverse
  arc s, P(s) = c + t (s - k^2 s^3 / 6) + n k s^2 / 2. Barrel: |P(s)| = r by
  three Newton steps from the straight line. Disc: s from z, which is exact;
- 2: `Hermite3D`'s two-point mode. The exact helix is taken at the slab's two
  qbar edges, and the cubic in t runs through both points with the tangents
  scaled by the transverse arc between them. z and the arc are linear in t.
  Barrel: |H(t)| = r by three Newton steps from t linear in r. Disc: t from z.
  It falls back to mode 0 where an edge is not reached.

The same 20 events (40-59), each prediction against surf::Helix::predict in double
at every fetched hit (`--chain-fast-check`, 18.4 M hits), and `--margins` against
`ref-quads.txt`:

| stage d | max \|dq\| | max \|dphi\| | worst error / window, q and phi | quads differing |
|---|---|---|---|---|
| 0 direct | 0.29 um | 3.0e-7 rad | 2.1e-4, 7.6e-5 | 6 |
| 1 one-point cubic | 256 um | 10 mrad | 0.10, 1.15 | 2454 |
| 2 two-point Hermite | 5.1 um | 1.5e-5 rad | 2.0e-3, 1.9e-3 | 10 |

Stage d time, from `--chain-phases`, fastest of 3 interleaved passes (load
2.0-2.4): mode 0 ~66 ms per event, mode 1 ~76, mode 2 ~85. Mode 2 needs two edge
solves per candidate and a Newton solve per barrel hit, where mode 0 solves each
hit in closed form.

**So the second point is not needed.** From one point, the exact helix solved
directly at each hit's qbar is both ~10x more precise than the two-point cubic
and cheaper. The one-point cubic is too coarse for the 15-30 cm c-d steps: its
worst phi error is larger than the window. Modes 1 and 2 stay as options for the
record; the batch finder uses mode 0.

What stage d costs is per candidate, not the prediction. The 66 ms is ~300 ns
for each of 216 k triplets, of which the per-hit predictions for ~4.3 fetched
hits are a small part. The rest is the gathers of nine hit coordinates from
three layers, the helix set-up, the edge predictions, the fetch and the
forwarding of misses.

### Fakes: where they are and what removes them (2026-09-27)

`seedsurf --quad-dump F` writes one row per quad kept after the cleaning. Each
row carries:
- its truth class, |eta| (a-d line) and pT_est;
- the transverse distance of the (a, b, c) and (a, c, d) circles from the beam
  spot;
- the stage c and d residuals over their windows;
- the four hits' cluster spans;
- the label composition;
- the best hit on the next outer P layer for the helix through b, c, d, with
  its residuals and label: OT1-P after a pixel d, OT2-P after an OT1-P d.

Analysis scripts are in the working report, `prep/fakes-2026-09-27/`. The chain
configuration is `seedsurf-chain.sh --chain-batch`, D121 PU200 events 40-69.
Tables are fitted on the first 15 events and applied to the last 15. "Tracks
lost" is findable tracks that lose every true quad, over all findable tracks.

Fake / decidable by the quad's |eta|: 0.47 (0-0.8), **0.81 (0.8-1.6)**, 0.55
(1.6-2.4), 0.12 (2.4-4). The fakes are not near-misses. In the barrel and
transition 91-95 % have every labelled hit from a different particle, and 20-23 %
have two unlabelled hits. At |eta| 0.8-1.6 the quads with OT1-P as the fourth
hit, which bridge the pixel gap, are 93 % fakes, 64.7 k of them in 30 events.
Quads with a pixel fourth hit are 61 % fakes.

What separates them:
- **Residual score**, the sum of the squared stage c and d residuals over their
  windows: median 0.09-0.18 for true quads, 1.1-1.6 for fakes. A cut at 1
  takes the regions to 0.17 / 0.57 / 0.28 / 0.06 for 0.37 % of tracks.
- **Cluster length along z in the barrel pixels.** `spanCols()` grows with
  |cot theta|: median 1 column at |cot| < 0.3, 8 at 4-12. A hit of another track
  has its own track's length. Keeping 99.5 % of true hits per layer and |cot|
  bin (0.1) takes the regions to 0.32 / 0.66 / 0.45 / 0.12 for 0.62 %. On the
  discs both spans are 1-2 and carry nothing.
- **OT2-P after an OT1-P fourth hit.** The window was fitted on true quads (q97,
  a + b / pT): 6.2 / pT mrad in phi and 2.5 + 1.5 / pT mm in z. At x1.5, 93 %
  of true quads in OT2 acceptance pass and 13 % of fakes.
- Weak or not usable:
  - OT1-P after a pixel fourth hit: the window is 1.4 + 11.3 / pT mrad and
    2.0 + 3.6 / pT mm; 93 % of true pass and 43 % of fakes. It stays
    information on the seed, not a cut.
  - The circle's d0: stage b already bounds D0, so barrel fakes have small d0
    too.

| cuts, last 15 events | tracks lost | quads / ev | fake: 0-0.8 | 0.8-1.6 | 1.6-2.4 | 2.4-4 | all |
|---|---|---|---|---|---|---|---|
| none | -- | 28746 | 0.47 | 0.81 | 0.55 | 0.12 | 0.353 |
| score < 1.5, OT2-P, shape | 0.95 % | 21901 | 0.17 | 0.32 | 0.33 | 0.08 | 0.157 |
| score < 1, OT2-P, shape | 1.17 % | 19998 | 0.10 | 0.21 | 0.22 | 0.05 | 0.098 |

Caveats:
- The cluster-shape tables come from simulated clusters.
- The score cut is a tighter window in disguise: an ellipse instead of the box
  of the q95 windows.
- The track loss is small partly because most found tracks have more than one
  true quad.

### Fakes in the finder (2026-09-27)

The three cuts of the section above now run inside the batch finder
(`--chain-batch`). Each is off by default:
- `--fk-score S`: the residual score, (dq_c/w)^2 + (dphi_c/w)^2 + (dphi_d/w)^2 +
  (dq_d/w)^2, each residual over its own window, below S. Stage c drops a
  triplet whose c part alone reaches S and stores the c part in the candidate,
  which grows from 28 to 32 bytes. Stage d adds the d part. The d phi term is
  sin^2(dphi) / sin^2(w), from the cross and dot products the phi cut already
  has.
- `--fk-shape`: the cluster length along z (`Hit::spanCols()`, now cached per
  hit in `SurfLayer::span_`) of every hit on a barrel pixel layer, within the
  band of true hits for |cot theta|. Stage b tests hits a and b on the a-b
  line. Stage c tests hit c on the a-b line, and stage d tests hit d on the a-c
  line. The bands are `windows-D121/shape.txt` (`--shape-win`, read by
  `seedsurf-chain.sh`): the central 99.5 % of true hits per layer and 0.1 bin of
  |cot|, with no cut in a bin of fewer than 50 hits.
- `--fk-ot2 F`: a quad whose d is on OT1-P needs an OT2-P hit for the helix
  through b, c, d. `SurfChainBatch::next_hit()` takes the best hit by the score
  of the offline study, among the hits within F x the window of the prediction.
  That best hit must lie within F x (a + b / pT) in phi and z. A helix that
  misses OT2-P, or reaches it outside |z| < zmax - 2 cm, passes. The window is
  `--ot2-win`, default -0.81 + 5.92 / pT mrad and 3.08 + 0.91 / pT mm.
- `--attach-ot1 F`: after the cleaning, each kept quad with a pixel d in OT1-P
  acceptance gets its best OT1-P hit, if the hit is within F x (2.16 + 9.99 / pT
  mrad, 2.32 + 3.57 / pT mm) (`--ot1-win`). This is information for the seed
  and not a cut.

The tables and windows were fitted on events 0-39 (q97 of true quads for the
windows; the working report's `prep/fakecuts-2026-09-27/fit.py`). They were
applied to events 40-99, which the numbers below come from. Union after the
cleaning, `--bind 0.05`, with OT2-P at x1.5 throughout:

| cuts | found / ev | tracks lost | quads / ev | fake: all | 0-0.8 | 0.8-1.6 | 1.6-2.4 | 2.4-4 |
|---|---|---|---|---|---|---|---|---|
| none | 1431.5 | -- | 27468 | 0.343 | 0.45 | 0.80 | 0.53 | 0.12 |
| score < 0.75 | 1423.2 | 0.55 % | 19338 | 0.115 | 0.15 | 0.46 | 0.22 | 0.04 |
| score < 0.75, OT2-P | 1421.2 | 0.69 % | 18890 | 0.094 | 0.13 | 0.28 | 0.21 | 0.04 |
| score < 1, OT2-P, shape | 1421.7 | 0.66 % | 20024 | 0.121 | 0.14 | 0.33 | 0.26 | 0.06 |
| **score < 0.75, OT2-P, shape** | 1417.3 | 0.95 % | 18653 | **0.083** | 0.09 | 0.23 | 0.18 | 0.04 |
| score < 0.5, OT2-P, shape | 1407.3 | 1.62 % | 16933 | 0.050 | 0.05 | 0.14 | 0.11 | 0.03 |

(Tracks lost: findable tracks no longer found, over all findable tracks.)
Changing the OT2-P factor between x1 and x2 moves the transition by +-0.02 and
the lost tracks by < 0.1 %.

Inside the finder the cuts act before the cleaning, and the offline study
applied them after it. At the same cuts (score < 1, OT2-P, shape) the finder
keeps more tracks and more fakes: fake 0.121 for 0.66 % of tracks lost,
against 0.098 for 1.17 % offline. A quad that the cleaning used to drop can
survive once the better quad sharing its hits is cut. At score < 0.75 the
finder does better than the offline cut at 1: fake 0.083 for 0.95 %.

The OT1-P attach (x1.5): of the true quads with a pixel d that are in OT1-P
acceptance, 92.9 % get a hit, and 89.2 % of those hits are the track's own. Of
the ~19 k kept quads with a pixel d per event, ~4.2 k are in OT1-P
acceptance. It costs ~12 ms per event (measured while other jobs ran).

(Measured with the nominal layer extents. With the extents of 29ca69d below,
the same rows read 1433.9 / 0.346 (none), 1424.0 / 0.122 (score < 1), 1419.7 /
0.084 (score < 0.75) and 1409.7 / 0.050 (score < 0.5), all with OT2-P and
shape.)

Time (the chain's own, events 40-59, fastest of 3 interleaved passes, black):
389 ms per event before, 394 with the cuts off, and **350 with score < 0.75,
OT2-P and shape**. That is less than with no cuts, because the stage c score
cut removes 43 % of the stage d candidates (216 k to 123 k per event). All
three tests run on the survivors of the existing box cuts, outside the vector
loops. Written into the loops, the same tests cost +30 ms per event with the
cuts off. With the cuts off, `--margins` against `ref-quads.txt` gives the same
6 differing quads as before.

### Speed after the fake cuts (2026-09-27)

Timing rule as in "Making the chain fast": the chain's own time, events 40-59,
black, fastest of 3 interleaved passes, each row measured against the row
before it in one session (the same binary reads within ~1 % across sessions).
"Cuts on" is `--fk-score 0.75 --fk-ot2 1.5 --fk-shape`. Acceptance: `--margins`
against the reference list, and for the cleaning the union truth row on events
40-99.

| commit | change | cuts off | cuts on |
|---|---|---|---|
| d7a3c1e | before this section | 389.5 | -- |
| 0b24bc4 | fake cuts in the finder (off: +5 ms) | 394.2 | 350.3 |
| 5bcee31 | a cot range per start pair | 367.6 | 321.6 |
| 1aef0d4 | barrel start pairs once for both sides; disc phi fetch | 343.2 | 299.8 |
| 87ec7f0 | stage d prediction as lanes | 339.8 | 297.6 |
| d7e144a | stage c and d fetch ranges in float | 333.4 | -- |
| 61837b8 | stage c one hit at a time, q first | 318.5 | -- |
| 29ca69d | layer extents over the hits, 2048 phi bins | 295.5 | 258.6 |
| f3440c2 | stage d q pre-filter, exact prediction for its survivors | **286.5** | **253.8** |

- **A cot range per start pair.** At set-up, `scan_starts()` scans r-z lines
  over the beam region with the flush's and `route()`'s own crossing tests. Each
  start pair gets the cot range of the lines it can use, widened by 0.02 in eta.
  A pair no line can use is skipped: the pairs ending on the last disc (23 27,
  26 27 and mirrors), which OT1-P closes to doublets. Before, they cost 19 Mcyc
  per side per event and pushed nothing. Stage b tests each doublet's cot against
  the range and narrows the fetch on B to the q range the cot range allows.
- **Shared barrel pairs.** `run_both()` fetches B1 B2, B2 B3 and B3 B4 once and
  sends each doublet to the side its line goes to. On a disc B, the phi fetch
  takes the geometric bound at the largest r the cot range reaches.
- **Stage c one hit at a time.** A stage c fetch returns 4.8 hits per candidate
  in 1.6 runs. The chunked mask and scan loops cost ~236 cycles per candidate
  (cycle counters), more than the arithmetic. A plain loop that tests q first,
  where ~95 % of the hits fail, is 15 ms per event faster.
- **The layer extents.** `LayerInfo::rin()` of the pixel barrel lies inside the
  inner radial shell: in layer 0 it is 2.868 cm, with hits from 2.750 cm. 54 / 48
  / 50 / 49 % of the hits of layers 0-3 have r below it (one event); the discs
  are exact. The fetches of stages b, c and d were built on the nominal extent,
  the double reference's too. The 256 phi bins hid it: with 1024 bins the batch
  finder lost 15 reference quads. `SurfLayer::fill()` now widens the extents to
  the event's hits. With it the list no longer depends on the phi binning, and
  2048 bins (3.07 mrad) take 20 ms off stages c and d. Physics, events 40-99:
  1431.5 -> 1433.9 found tracks per event and fake 0.343 -> 0.346 with the cuts
  off; 1417.3 -> 1419.7 and 0.083 -> 0.084 with them on.
- **Stage d q pre-filter.** Stage d spent ~490 cycles per triplet in the fetch
  and mask for ~4.3 hits. It now also predicts q at the target's middle qbar and
  takes q as a quadratic in qbar through the three predictions, which K2c
  measured within 2 % of the window. A hit off it by more than 1.1 x the window
  + 20 um is skipped; the others get the exact prediction and cut, one at a
  time. The truth rows are unchanged.
- **The cleaning** (`seedsurf.cc`, outside the chain's time) computed each
  quad's tier inside the sort comparator: 6 % of the job's branch mispredicts.
  Computed once per quad, the job takes 1.0 s less per 60 events (c102f12).
- **The reference list** is regenerated with the double chain as
  `ref-quads-2.txt` in the working report (`prep/chain-kernels-2026-09-27/`),
  1518676 quads. The double chain at 256 and at 2048 bins differs there by 1
  quad, 0.17 % inside the stage c phi window. That is an edge of the double
  chain's own stage c fetch. Against it the batch finder differs by 4 + 3
  quads, all within 1e-4 of a cut.

Tried and not kept:
- Gathering the block's a and b hit coordinates into arrays before the stage c
  loop, so their loads overlap: 5 ms per event slower.
- The fake-cut tests inside the vector loops: +30 ms per event with the cuts
  off. They run on the survivors of the box cuts.
- Stage c's q test over a run without a branch, with the passing hits compacted
  first: 17 ms per event slower than the plain q-first loop.

Where the time is now (cuts off): start doublets ~116 ms, of which the flush
(holes, states, routing, queueing of 2.6 M lanes per event, ~63 cycles each) is
~41 %; stage c ~120 ms over 1.43 M candidates (~290 cycles each); stage d ~51 ms
over 216 k triplets. Each is roughly proportional to the number of doublets the
configuration makes. Next: b-major stage c, a c-hit list per b hit sorted by
slope, as the barrel finder's K3/K4. Measured for it: 14.6 stage c candidates
per (b hit, target) on average, but the queues are ordered by start pair and a
hit, not by b. So it needs either a b-major stage b or a sort of the stage c
queues, ~1.4 M candidates per event. The chain's crossing envelopes (`SurfOwnership`)
also take the barrel extents from `LayerInfo`, and are not changed yet.

### The rest of the seeding: layer fill and cleaning (2026-09-27)

The chain's own time was the only number measured so far. With timers around
every phase of a `seedsurf-chain.sh --chain-batch` job (events 40-59, cuts off,
ms per event): reading the event 50, layer fill 22, the chain 288, the per-quad
loop after it 72 (truth labels, and the cleaning score recomputed in double with
`surf_eval` for all 76k quads), the cleaning 58. Of these the layer fill, the
cleaning score and the cleaning belong to the seeding, and ~150 ms per event of
it had no timer. The job now prints the cleaning time. Fastest of 3 interleaved
passes, the baseline being f3440c2 with the same timers:

| | chain | layer fill | cleaning |
|---|---|---|---|
| cuts off, before | 289.5 | 21.4 | 55.7 |
| cuts off, after (1ce7f8e) | 293.9 | 14.5 | 9.4 |
| cuts on, before | 258.3 | 22.9 | 39.3 |
| cuts on, after | 262.5 | 15.7 | 7.3 |

- **The score from the finder** (0b083e6). The batch finder computes the
  cleaning score at `take()` in float, with the pattern's own windows, and hands
  it out with the quad. That costs the chain ~4 ms per event and takes the
  double `surf_eval` off the batch path (it stays for the double chain).
- **The cleaning** (03d1feb): a stable radix sort on one 32-bit key, per-hit
  lists in flat arrays, and only the 5 - N shortest lists walked for "shares
  >= N hits". The steps are in the commit message.
- **No double cache** (1ce7f8e): `SurfLayer::fill()` builds the double r, phi
  of every hit only for the double chain.

The union truth row on events 40-99 is identical to the printed precision, cuts
off and on, and the OT1-P attach row too. The truth files differ by one
undecidable quad in a few eta bins: the float score reorders near-ties.
`--margins` against `ref-quads-2` is unchanged.

**vdt, measured null.** Replacing the finder's float `atan2`, `sin`, `asin`
and `sincos` with vdt, forced inline, left the margins identical and the time
unchanged to +5 ms. A microbenchmark (4096 floats, ns per element, `-mavx`)
says why: `sin` is 1.98 in a loop GCC vectorises through glibc's libmvec
against 1.86 for vdt; `atan2` is 2.66 through libmvec but 10.25 for vdt,
whose `fast_atan2f` GCC does not vectorise (its branches and swaps). The
finder's calls are one per candidate, not in loops, where vdt `atan2` is ~5 ns
faster than glibc (10.5 against 16.0), a few ms per event at most. Not kept.

**`LayerInfo::rin()`, the cause.** `MkFitGeometryESProducer` takes it as the
smallest radius over each module's 8 corners, and a flat module's closest point
to the beam is the foot of the perpendicular. On D121 it is high by 0.06-0.14 cm
in the pixel barrel, 0.23-0.53 cm in TBPS and 0.10-0.16 cm in TB2S; layer 0's
true minimum is 2.7425 cm, the hits' 2.750 less the 75 um half-thickness.
`rout` and the z extents are right.

**Stage c by b hit, measured slower, not kept (2026-09-28).** The doublets
queued at a target, as they come, share their b hit in runs of 1.01; sorted
within blocks of 256 or 4096 in groups of 3.4 or 4.6, over the whole queue in
groups of 14.6 (14.1 per start pair and b hit, which a b-major stage b would
give without a sort). With the queues radix-sorted by b hit, which costs 35 ms
per event, three grouped kernels were built behind a flag. Each fetches the
target over the union of the group's windows once. The members then test
(a) all hits of the union as lanes, (b) a range of the hits sorted by their r-z
slope from b, or (c) the slope buckets of a counting sort. The union holds 43
hits per member, the sum of the members' own fetches (14.6 x 2.9): the
windows of doublets that share b are almost disjoint in q on the target,
because their a-b slopes differ. Of the 43 hits, 0.85 pass a member's q cut and
7.2 its phi cut. So a group shares the b hit and nothing else, and the old
per-candidate fetch was not the expensive part. Stage c, the queue sort not
counted, events 40-59: 118 ms per event one at a time; (a) ~450, (b) ~190,
(c) ~148. Prefetching the a and b hits of the candidate 8 ahead in the old loop:
the forward pass +6 ms per event. The margins of (a)-(c) against `ref-quads-2`
were the 4 + 2 quads plus one new quad 0.17 % inside the c phi window, a hit
the old per-candidate fetch missed at its edge.
