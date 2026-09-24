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
```

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

## Where it goes next

Everything above keeps the prototype's arithmetic exactly, doubles, `hypot`,
`asin` and all. The remaining large items are the helix (240 ns per triplet in
double precision) and the per-doublet lookup. Both want float arithmetic and
vectorised maths, and that ends bit-equality with the prototype. The
acceptance test would then become "differences only at a cut edge", counted
and bounded.
