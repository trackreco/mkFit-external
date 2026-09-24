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

### Step B: fixed-point integers in the per-pair tests

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
