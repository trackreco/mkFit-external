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

20 events, tile finder, fastk, one repetition, D0_max 1 mm:

| pT_min [GeV] | findable /ev | efficiency | fake / decidable | undecidable | extra true quads per found track | doublets /ev | ms /ev |
|---|---|---|---|---|---|---|---|
| 0.2 | 1870 | 0.522 | 0.292 | 0.380 | 0.108 | 2.33 M | 116 |
| 0.5 | 940 | 0.795 | 0.119 | 0.211 | 0.107 | 1.18 M | 36 |
| 0.9 | 376 | 0.940 | 0.090 | 0.184 | 0.117 | 0.84 M | 22 |
| 2.0 | 72 | 0.984 | 0.091 | 0.192 | 0.150 | 0.60 M | 13 |

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
