# RdfTrace: tracing the v2p2 search

The v2p2 finder records what each candidate did, layer by layer, into the
`Event`. RDataFrame reads it back, one entry per event.

## Turn it on

In the build directory's `Makefile.config`:

    EVENT_RDF_TRACE := 1

This sets `MKFIT_TRACE`, `WITH_ROOT`, `WITH_REVE` and `TBB_DEBUG`. With
`TBB_DEBUG` the TBB loops run serially, so an event is traced on one thread.
Parallelism comes from several events in flight.

Only the v2p2 path records anything: run with `--build-mimi --build-mimi-v2p2`.
V1 and V2 run fine in a trace build, but their trace is empty.

## Where it is written

| core file | records |
|---|---|
| `MkBuilder.cc`, `trace_begin_v2p2_search()` | at the start of each v2p2 search: a meta per seed (once), a stage, the root state |
| `MkFinderV2p2.cc` | `TrLayerSearch`, `TrHitMatch`, the new `TrCandState`s after each layer |
| `MkFinderV2p2Structures.cc` | `TrKalmanUpdate` |
| `MkFinderV2p2-BkFit.cc` | the backward fit's stage and states, `TrBkFitUpdate` |
| `TrackStructures.h` | `CombCandidate::m_trace_meta_id`, which ties a seed's candidates to their meta across stages |

All under `#ifdef MKFIT_TRACE`. Core calls the `trace_*()` methods of `Event`
(`Event.h`); `DataFormats/RntConversions.h` fills records from core state
(`fill_layer()`, `fill_window()`, `fill_bins()`, `make_hitmatch()`).

## Where it lands

`Event` members, cleared per event in `Event::reset()`:

| member | one per | links |
|---|---|---|
| `trCandMetas_` | seed's candidate family | `seed` (into `trSeeds_`), `global_seed` (the seed's label, see Traps), `cand` (into `candidateTracks_`, final only), `stage_ids[3]` |
| `trCandStages_` | stage of a family | `stage` 0 forward search, 1 backward fit, 2 backward search; `root_state_id`, `final_state_id`, `parent_stage_id` |
| `trCandStates_` | candidate state at a layer | `parent_id` (the tree), `meta_id`, `stage_id`, `layer` |
| `trLayerSearches_` | search of one layer from one state | `state_id`; windows, bin ranges, mini-propagator ends |
| `trHitMatches_` | hit scanned | `state_id`, `search_id` |
| `trKalmanUpdates_` | Kalman update tried | `hit_match_id`, `state_id_in`, `state_id_out` (-1 if not accepted) |
| `trBkFitUpdates_` | hit of the backward fit | `state_id_in`, `state_id_out`; a chain, not a tree |

Accessors `tr_candmeta(i)`, `tr_candstate(i)`, ... in `Event.h`. Record types
in `DataFormats/RntStructs.h`, dictionary in `DataFormats/RootDFs_LinkDef.h`.

## Reading it

| file | what |
|---|---|
| `RdfSources.{h,cc}` | `MakeEventDF()`: one RDF entry per `Event` |
| `AnRunEvMap.h` | projects an `Event`'s vectors into `RVec` columns |
| `AnRun.{h,cc}` | the canned analyses: `RunBasicSeedCandCheck()`, duplicate checks, T5s into pixels, stage-2 root-state checks |
| `AnRun-SimEff.cc` | `RunSimTrackEfficiency()` |
| `CanvasGroup.{h,cc}` | groups of histograms drawn and written together |
| `RdfVectorSource.{h,cc}` | one RDF entry per vector element. Not used by the current path; kept as a working example of a custom `RDataSource` |

Probes that read these records event by event: `RecoTracker/MkFitCMS/standalone/Probe/`.

## Drivers, in the shell (`--shell`)

| call | does |
|---|---|
| `s.TraceFwdSearch(n_ev, prefix, n_thr)` | outward forward search on the iteration's seeds, then `AnRun`. The events stay reachable as `ar.m_ev_vec` |
| `s.ProcessEventStd()` | one event, standard search, traced |
| `s.LoopNEventsHlt(n)`, `s.ProcessEventHlt()`, `s.RunLSTintoPix()` | inward search of LST T5 seeds into the pixels (`Shell-LST.cc`) |
| `s.TraceVectorSource()`, `s.TraceEventSource()` | smoke tests of the two RDF sources |

Example:

    ./mkFit --input-file <sample.bin> --geom CMS-phase2 --seed-input cmssw \
      --build-mimi --build-mimi-v2p2 --num-iters-cmssw 1 --num-events 10 --shell \
      --shell-command 's.TraceFwdSearch(10, "fwd", 4)'

## Traps

- Ids can be -1, and the `tr_*()` accessors do not check. Test before indexing.
- `search_id` is per event: key a map on `(event, search_id)`.
- `TrHitMatch::layer` is the layer searched; the `TrCandState` it points to sits
  on the previous layer.
- `trCandStates_` holds all stages: filter on the stage.
- Iterations are not separated: `TrCandStage::iteration_idx` is never set, and
  everything goes into the same per-event vectors.
- Seeds: `TrCandMeta::seed` indexes `trSeeds_`, the trace's own copy of the
  iteration's seeds, made by `Event::build_trace_maps_etc()`, because the
  driver's per-iteration seed vector may be temporary. With several iterations
  it holds the last one. Seed truth per meta is in `trSIFHforSeedByMeta_`.
  `TrCandMeta::global_seed` holds the seed's label, not an index: it equals the
  index in `seedTracks_` only after `relabelSeedTracksSequentially()`, and the
  per-iteration copy keeps no other link to `seedTracks_`.
- `fail` is 0 on every backward-fit record, because `getS` has no bracket. It
  does not mean the fit did not fail.

Not to be confused with `RntDumper/`, the older TTree dumper for V2's hit
selection.
