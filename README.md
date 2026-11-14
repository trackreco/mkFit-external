# mkFit-external

External header-only packages required for standalone build of mkFit from
CMSSW packages RecoTracker/MkFitCore and RecoTracker/MkFitCMS.

1. nlohman/json

2. Extract of SMatrix package from ROOT

## Setting up standalone build

The sources of `RecoTracker/MkFitCore`, `RecoTracker/MkFitCMS` and `FWCore/Utilities` come from a CMSSW src directory. Everything else comes from this repository: the standalone trees `RecoTracker/MkFitCore/standalone` and `RecoTracker/MkFitCMS/standalone`, at the same paths as they had in CMSSW, the external headers and the geometry plugins. Some data files (geometry binary dumps) are automatically retrived from the web (from `phi1.t2.ucsd.edu`).

Run `configure` from a build directory outside this checkout. It writes the top-level `Makefile`, creates the `objs-*` directories and copies the build options, `Makefile.config`, into the build directory; edit that copy. With an existing CMSSW environment:
```
git clone git@github.com:trackreco/mkFit-external.git
# the build also needs FWCore/Utilities in the src directory:
(cd $CMSSW_BASE/src && git cms-addpkg FWCore/Utilities)
mkdir standalone
cd standalone
../mkFit-external/configure $CMSSW_BASE/src
export -n INTEL_LICENSE_FILE
make -j 32
# or, for build with ROOT enabled
make WITH_ROOT=1 -j 32
```

`configure` refuses a CMSSW src directory that still has its own `RecoTracker/MkFitCore/standalone`, i.e. one from before the standalone trees moved here. Such a tree builds with its own `RecoTracker/MkFitCore/standalone/configure`.

Without CMSSW environment you still need to run `configure` as above for initial setup. `tbb` can be taken from the system ... but will likely require some sym-link magick as CMSSW includes pick it up from `oneapi/tbb` directory which is not the standard tbb deployment on GNU/Linux.

Notes on running standalone:
```
export LD_LIBRARY_PATH=.
# or
export LD_LIBRARY_PATH=.:${LD_LIBRARY_PATH}

# phase1 multi-iteration test
./mkFit --cmssw-n2seeds --input-file /data2/slava77/samples/2021/11834.0_TTbar_14TeV+2021/AVE_50_BX01_25ns/memoryFile.fv6.default.211008-c6b7c67.bin --read-cmssw-tracks --build-mimi --num-iters-cmssw 9 --backward-fit --num-events 100 --quality-val --num-thr 32 --num-thr-ev 8

# phase2 initialStep test
./mkFit --input-file /data2/slava77/analysis/CMSSW_12_4_0-mkFit/work/18f4388/39634.0_TTbar_14TeV+2026D88PU/memoryFile.fv7.clean.writeAll.recT.allSeeds.220712-703c46d.bin  --geom CMS-phase2 --seed-input cmssw --seed-cleaning n2 --build-mimi --num-iters-cmssw 1 --num-events 100 --quality-val

# binary file maker
 ./writeMemoryFile  --input ~/mnt/phi3-slash/ceph/cms/store/user/slava77/CMSSW_11_2_0_mkFit/work/24feee2/2021/10muPt0p2to1HS/trackingNtuple.root --output 10mu-lowpT-1000.bin  --maxevt 1000 --write-rec-tracks

# JSON default dumper
./mkFit --json-save-iterations mkfit-phase1-%s.json --geom CMS-phase1 --num-events 0
./mkFit --json-save-iterations mkfit-phase2-%s.json --geom CMS-phase2 --num-events 0

# Validation in standalone/
SA_PATH=<path-to>/RecoTracker/MkFitCore/standalone
ln -s ${SA_PATH}/plotting .
ln -s ${SA_PATH}/val_scripts .
ln -s ${SA_PATH}/web .
ln -s ${SA_PATH}/xeon_scripts .
# Note the last argument, this is for phase2
val_scripts/validation-cmssw-benchmarks-multiiter.sh forPR --mtv-like-val TTbar_phase2
web/collectBenchmarks-multi.sh <dir-name> forPR
```
Up to date samples for running standalone are listed below (tracking ntuples and step1/step2 files can be found in the same folder):

```      
/ceph/cms/store/user/legianni/generate-phase2/10mu-newgeom/ntuple_pt0p1to1_*.bin  
/ceph/cms/store/user/legianni/generate-phase2/10mu-newgeom/ntuple_pt1to1000_*.bin  
/ceph/cms/store/user/legianni/generate-phase2/ttbar-newgeom/ntuple_ttbar.bin
/ceph/cms/store/user/legianni/generate-phase2/ttbar-newgeom/ntuple_ttbar_PU.bin

/ceph/cms/store/user/legianni/generate-phase1/10mu-phase1-newgeom/ntuple_pt0p1to1_*.bin   
/ceph/cms/store/user/legianni/generate-phase1/10mu-phase1-newgeom/ntuple_pt1to1000_*.bin
/ceph/cms/store/user/legianni/generate-phase1/ttbar-phase1-newgeom/ntuple_ttbar_PU.bin
```

## Probes, knobs and tracing

Two READMEs map the tools for measuring the v2p2 finder from the standalone shell
(builds with `WITH_ROOT`):
- [RecoTracker/MkFitCMS/standalone/Probe/README.md](RecoTracker/MkFitCMS/standalone/Probe/README.md):
  the probes, which read after the event and change nothing, the knobs and the
  taps, which are `Shell` methods, and the driver scripts;
- [RecoTracker/MkFitCore/standalone/RdfTrace/README.md](RecoTracker/MkFitCore/standalone/RdfTrace/README.md):
  the trace of the v2p2 search (`EVENT_RDF_TRACE`) and its RDataFrame analysis.

## RntDumpers -- ROOT NTuple Dumpers

This is the latest attempt to have a somewhat standard way of dumping tunbing / debugging information directly into ROOT RNTuples or TTrees. Before we had dumpers via ROOT text tree files (slurped into ROOT via `TTree::ReadFile()`) -- this was slow and produced HUGE intermediate files. This resides in `MkFitCore/standalone/RntDumpers`.

So far there is a single implementation of such dumper, implemented in `MkFinder_selectHitIndices.icc` with call-out scattered through `Mkfinder.cc` as if-defed pieces of code, `#ifdef RNT_DUMP_MkF_SelHitIdcs`. To enable it, build WITH_ROOT and uncomment `CPPFLAGS += -DRNT_DUMP_MkF_SelHitIdcs` in `Makefile.config`, around line 98.

This will dump one tree entry per layer loop in `MkBuilder::find_tracks_in_layers()`. It contains a `std::vector<CandInfo> ci` structs -- detailed inforamtion about each candidate processed in this loop, including the search windows and propagation dcetails as well as detailes about considered hits. All the structrues used in dumping are defined in `RntStruct.h` and "standard" conversions from mkFit types to them are in `RntConversions.h`. Further, each `CandInfo` object contains a `std::vector<HitMatchInfo> hmi` -- this contains information about how each hit was matched in `selectHitIndicesV2()` and also about how it was scored (if it passed all the intermediate checks) in `findCandidatesCloneEngine()`.

Notes:
1. A lot of score information is missing because the hit chi2 is too large or because it is not passing the compatibility tests. One might need to rework score dumping logick to proceed anyway (or loosen the cuts).
2. Invalid hits that are added to canidadets are clearly not dumped.
3. One might need enlarge bin-search window scale-up factors in `selectHitIndicesV2()` for some very specific cases ... they are already plenty large for most use-cases (i.e., they work OK for pixelLess iteration).
