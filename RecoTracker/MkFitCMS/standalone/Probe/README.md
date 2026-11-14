# Probe: measuring the standalone mkFit from its shell

Four kinds of tool, by what they do to a run:

- **probe**: reads after the event, or a synthetic sample, and changes nothing.
- **knob**: changes what the algorithm does. A `Shell` method, in the section
  "Knobs" of `../Shell.h`, bodies in `../Shell-Knobs.cc`: `s.SetMaxCands(6)`,
  `s.SetWindowDq(...)`, `s.SetTrackScorer(...)`, ...
- **tap**: in-core instrumentation, switched by a flag, that collects during the
  run: `s.SetTapFinalBeam(bool)`, `s.SetTapScoreStats(bool)`. A probe reads it
  afterwards.
- **driver**: runs configurations over events and calls the above: `probe.C`
  and the three scripts here.

## Probe groups

| group | file | reads |
|---|---|---|
| closure | `Probe-Closure.cc` | synthetic tracks with an exact answer, no event |
| checks | `Probe-Checks.cc`, `probe-shs.C` | the input, with a verdict |
| perf | `Probe-Perf.cc` | the final tracks against truth |
| trace | `Probe-Trace.cc` | the `Tr*` records of a trace build, see `RecoTracker/MkFitCore/standalone/RdfTrace/` |
| tap readers | `Probe-Taps.cc` | what the taps collected |

`Probe.h` declares them. They are built into `libMkFitCMS.so` with `WITH_ROOT`.
`Probe-Trace.cc` needs the trace build (`EVENT_RDF_TRACE`); without it the file
compiles to nothing and its probes are not declared. `probe_kinfo_*` and the
`cov` tree of `probe_search_*` also need `MKFIT_TRACE_KALMAN_DEBUG` in
`Makefile.config`. The records the probes write are `Probe*` structs,
`RecoTracker/MkFitCore/standalone/DataFormats/ProbeStructs.h`, with a dictionary
in `libMkFitRootDataFormats.so`.

## Driving them

From a build directory, load `probe.C` in the shell. `<SACMS>` is the directory
its `Makefile` exports as `SACMS`:

    ./mkFit --geom CMS-phase2 --seed-input cmssw --input-file <sample.bin> \
      --num-events 20 --num-thr 1 --build-mimi --build-mimi-v2p2 --shell \
      --shell-command 'gROOT->ProcessLine(".L <SACMS>/Probe/probe.C")'

`probe.C` says `using namespace mkfit`, so the probes are called by their names,
and adds short forms (`probe_eff_ev()`, `probe_sr_ev()`, ...). An event's probes
take `s.event()` after `s.ProcessEventStd()`:

    probe_assoc_mtv(true); probe_eff_reset();
    s.GoToEvent(1); s.ProcessEventStd(); probe_eff_ev(s.event(), "v2p2");
    ...
    probe_eff_report("eff");

The trace probes take a seed's truth from the trace's own copy of the seeds,
which `Event::build_trace_maps_etc()` makes. The collecting drivers
(`s.TraceFwdSearch()`, ...) call it; `s.ProcessEventStd()` does not, so call
`s.event()->build_trace_maps_etc()` after it. Without it the trace probes say
so and record no seed truth.

The scripts run from a build directory and load `probe.C` themselves:

| script | does |
|---|---|
| `v2p2-eff-v2-vs-v2p2.sh <n_events> <sample> [out_prefix]` | V2 against v2p2 in one process, `probe_eff` |
| `v2p2-eff-knobs.sh <n_events> <sample> <out_prefix> 'NAME=cmd;cmd' ...` | several knob settings in one process, `probe_eff` |
| `v2p2-timing.sh <n_events> <n_reps> <sample>` | single-thread timing; refuses a build linked with ROOT |

The two `v2p2-eff-*` scripts read `--read-cmssw-tracks`, so the sample must carry
the CMSSW tracks.

## Tap mode

The taps accumulate into unsynchronised globals, so they are correct only with
one thread inside the event and one event in flight. Switching any tap on sets
tap mode: `collect_events_std()` and `collect_events_hlt()` then run one event in
flight, and in a build without the trace (whose in-event loops are parallel) a
`tbb::global_control` limits TBB to one thread. Switching all taps off clears it.
`s.TapMode()` says whether it is on.

## The live families

| family (short forms in `probe.C`) | group | writes |
|---|---|---|
| `probe_gen_closure` (`probe_gen`) | closure | TTrees `closure` (`ProbeClosure`) and `cfgs` (`ProbeCfgInfo`) |
| `probe_cov_transport` (`probe_covx`) | closure | TTree `covxport` (`ProbeCovXport`) |
| `probe_dalpha_*` (`probe_dalpha_ev`) | closure | TTree `steps` (`ProbeStep`): turn angles of real tracks, reads events |
| `probe_geom_check_*` (`probe_geom_*`) | checks | printout: is the sample's geometry this build's |
| `probe_sample_info` (`probe_info`) | checks | printout: what is in the sample |
| `probe_shs_*` (`probe-shs.C`) | checks | printout: the per-sim-hit truth states against their hits |
| `probe_assoc_mtv` | perf | sets the association rule of the perf probes and of quality-val |
| `probe_eff_*` (`probe_eff_ev`, `probe_eff_cmssw_ev`) | perf | `<prefix>.txt`, `<prefix>.root`: efficiency, fakes, duplicates, length, d(pT)/pT, paired across configurations |
| `probe_track_eff_*` (`probe_te_*`) | perf | printout: true pixel and disc hits on the final tracks |
| `probe_trkdump_*` (`probe_trkdump_ev`) | perf | TTree `trk`: the best-matched reco track of each sim track |
| `probe_search_*` (`probe_sr_*`) | trace | TTrees `search` (`ProbeSearchHit`), `miss` (`ProbeSearchMiss`), `cov` (`ProbeCovStep`) |
| `probe_bkfit_trace_*` (`probe_bt_*`) | trace | `<prefix>.txt`: chi2 through the backward fit |
| `probe_kinfo_*` (`probe_kinfo_ev`) | trace | printout: the exact gain test of the Kalman update |
| `probe_rank_*` (`probe_rank_ev`) | trace | printout: the worst end-of-layer rank of each final candidate |
| `probe_final_beam_report` | tap reader | printout: the end-of-search beam purity |
| `probe_score_term_stats_report` | tap reader | printout: the mean likelihood terms, the constants of `s.SetScoreTerms()` |

The records of the studies of 2026-09 and 2026-10 that used the removed probes
are kept outside this repository.
