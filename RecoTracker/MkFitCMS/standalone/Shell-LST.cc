#include "RecoTracker/MkFitCMS/standalone/Shell.h"

#include "RecoTracker/MkFitCore/src/Debug.h"

#include "RecoTracker/MkFitCMS/interface/runFunctions.h"

#include "RecoTracker/MkFitCore/interface/HitStructures.h"
#include "RecoTracker/MkFitCore/interface/MkBuilder.h"
#include "RecoTracker/MkFitCore/src/MkFitter.h"
#include "RecoTracker/MkFitCMS/interface/MkStdSeqs.h"
#include "RecoTracker/MkFitCMS/standalone/MkStandaloneSeqs.h"

#include "RecoTracker/MkFitCore/interface/Config.h"
#include "RecoTracker/MkFitCore/interface/IterationConfig.h"
#include "RecoTracker/MkFitCore/standalone/ConfigStandalone.h"

#include "RecoTracker/MkFitCore/standalone/Event.h"

#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"

#ifdef MKFIT_TRACE
#include "RecoTracker/MkFitCore/standalone/RdfTrace/AnRun.h"
#endif

#include "TROOT.h"

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"
#include "oneapi/tbb/task_arena.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <mutex>

namespace mkfit {

  // Runtime selectors for the HLT-seed driver; see Shell.h.
  Shell::HltSeedKind_e Shell::s_hlt_seed_kind = Shell::HSK_T5;
  bool                 Shell::s_hlt_chop_pixels = false;
  std::map<int, std::vector<HitOnTrack>> Shell::s_chopped_hits;

  namespace {
    // Track scorer "phase2:LstIntoPix", for LST T5 tracks extended into the pixels:
    // a found hit is worth more than the largest chi2 (30), and tail holes cost as
    // much as inside holes. Select it with val_track_scorer("phase2:LstIntoPix").
    float trackScoreLstIntoPixels(const TrackScoreInput &in) {
      const int ntailholes = in.penalize_tail_holes ? in.n_tail_holes : 0;
      float bonus = 30;
      float penalty = Config::missingHitPenalty_;
      float tailPenalty = Config::missingHitPenalty_;
      float score = bonus * in.n_found_hits - penalty * in.n_inside_holes - tailPenalty * ntailholes - in.chi2;
      return score;
    }

    struct register_lst_track_scorer {
      register_lst_track_scorer() {
        IterationConfig::register_track_scorer("phase2:LstIntoPix", trackScoreLstIntoPixels);
      }
    } rlts_instance;
  }  // namespace


  //===========================================================================
  #pragma region Event Loops
  //===========================================================================

  void Shell::LoopNEvents(int N_events) {
    N_events = std::min(N_events, m_evs_in_file);
    for (int ev = 1; ev <= N_events; ++ev) {

      printf("\n##### BEG Event %d #####\n\n", ev);

      GoToEvent(ev);

      // find max pT5 index (seeds are: pT5, T5, p)
      int min_seed = 0; // for now always 0
      int max_seed = 0;
      auto const &seeds = m_ctx.ev->seedTracks_;
      int ns = seeds.size();
      for (int si = 0; si < ns; ++si) {
        auto const &s = seeds[si];
        auto hot = s.getHitOnTrack(0);
        auto &li = tracker_info()->layer(hot.layer);
        if ( ! li.is_pixel() || s.algoint() != 4) {
          max_seed = si - 1;
          break;
        }
      }
      int num_seeds = max_seed - min_seed + 1;
      printf("Selected %d LST pT5 seeds [%d, %d].\n", num_seeds, min_seed, max_seed);

      if (num_seeds > 0) {
        RunLSTintoPix(m_ctx, mkfit::Shell::SS_IndexPreCleaning, min_seed, num_seeds);
      }

      printf("\n##### END Event %d #####\n", ev);
    }
  }

  //===========================================================================
  // HLT seed pre-selection: which seeds are the pT5s / T5s / pixel-only ones
  //===========================================================================

  const char* Shell::hlt_seed_kind_name(HltSeedKind_e k) {
    switch (k) {
      case HSK_pT5: return "LST pT5s";
      case HSK_T5:  return "LST T5s";
      case HSK_pix: return "pixel tracks";
    }
    return "<unknown>";
  }

  Shell::HltSeedBlocks Shell::hlt_seed_preselect(const Event &ev, const TrackerInfo &ti,
                                                 int algo, HltSeedKind_e want)
  {
    HltSeedBlocks r;
    const auto &seeds = ev.seedTracks_;
    const int ns = (int) seeds.size();

    bool in_algo = false;
    int last_match = -2;   // index of the previous seed of the wanted kind

    for (int si = 0; si < ns; ++si) {
      const Track &s = seeds[si];

      if (s.algoint() != algo) {
        if (in_algo) break;   // seeds are grouped by algo, we are past our block
        continue;
      }
      in_algo = true;

      const bool f_is_pix  =   ti.layer(s.getHitOnTrack(0).layer).is_pixel();
      const bool l_is_strp = ! ti.layer(s.getLastHitLyr()).is_pixel();

      HltSeedKind_e k;
      if      ( f_is_pix &&  l_is_strp) k = HSK_pT5;
      else if (!f_is_pix &&  l_is_strp) k = HSK_T5;
      else if ( f_is_pix && !l_is_strp) k = HSK_pix;
      else continue;   // starts in strips and ends in pixels -- not one of ours

      switch (k) {
        case HSK_pT5: ++r.n_pT5; break;
        case HSK_T5:  ++r.n_T5;  break;
        case HSK_pix: ++r.n_pix; break;
      }

      if (k == want) {
        if (r.first < 0)
          r.first = si;
        else if (si != last_match + 1)
          r.contiguous = false;   // a gap: the one-solid-block assumption broke
        ++r.count;
        last_match = si;
      }
    }
    return r;
  }

  // Per-event HLT/LST processing: seed classification by hit pattern, then
  // RunLSTintoPix() on the selected block, quality validation and the bad-seed
  // report. Drivers below own the event loop and the quality accumulator.
  void Shell::ProcessEventHlt(EvCtx &ctx, const int wanted_algo) {
    printf("\n##### BEG Event %d ##### HLT seeds\n\n", ctx.ev->evtID());

    ctx.ev->filterOutMislabeledHitsInSimTracks();
    ctx.ev->relabelSeedTracksSequentially();

    // Which kind of HLT seeds to run over, and whether to chop a pT5's pixel
    // hits off before searching. Runtime now, so the pT5-chopped control can be
    // run against plain T5 in one binary: a pT5's pixel hits were FOUND by the
    // upstream reconstruction, so chopping them gives a denominator of hits
    // known to be findable -- unlike a plain T5, where they were not found.
    const HltSeedKind_e wanted_kind = s_hlt_seed_kind;

    const HltSeedBlocks blk =
        hlt_seed_preselect(*ctx.ev, *tracker_info(), wanted_algo, wanted_kind);
    const int min_seed = blk.first;
    const int num_seeds = blk.count;

    printf("HLT seed kinds for algo %d: pT5 = %d, T5 = %d, pix = %d\n",
           wanted_algo, blk.n_pT5, blk.n_T5, blk.n_pix);
    if ( ! blk.contiguous)
      printf("*** WARNING *** %s do not form one contiguous block -- everything below\n"
             "                takes [first, first+count), so the selection is wrong.\n",
             hlt_seed_kind_name(wanted_kind));

  #ifdef MKFIT_TRACE
    ctx.ev->seedVecInsp_.n_pTNs = blk.n_pT5;
    ctx.ev->seedVecInsp_.n_TNs  = blk.n_T5;
    ctx.ev->seedVecInsp_.n_ps   = blk.n_pix;
  #endif

    printf("Selected %d %s as seeds [%d, %d].\n", num_seeds, hlt_seed_kind_name(wanted_kind), min_seed, min_seed + num_seeds);

    const bool print_seed_summary = false;
    const bool print_seed_details = false;
    if (print_seed_summary) {
      for (int si = min_seed; si < min_seed + num_seeds; ++si) {
        const Track &t = ctx.ev->seedTracks_[si];
        auto sifh = ctx.ev->simInfoForTrack(t);
        if (sifh.is_set()) {
          printf("  si=%d n_hits=%d n_match=%d n_pix=%d n_pix_match=%d frac=%f lab=%d\n",
                si, sifh.n_hits, sifh.n_match, sifh.n_pix, sifh.n_pix_match, sifh.good_frac(), sifh.label);
        } else {
          printf("  si=%d no label\n", si);
        }
        if (print_seed_details) {
          mkfit::print("  SEED", si, t, *ctx.ev);
        }
      }
    }

    // HACK -- seed validation -- comment out the RunLSTintoPix() call below.
    // ctx.tracks = ctx.ev->seedTracks_;

    if (num_seeds > 0) {
      RunLSTintoPix(ctx, mkfit::Shell::SS_IndexPreCleaning, min_seed, num_seeds);
    }

    {
      // quality_val() ends in add_to_quality_sum(), which writes the static
      // Quality::s_quality_sum -- guard it so several events can be in flight.
      // Also keeps each event's quality printout in one piece.
      static std::mutex s_qual_mutex;
      std::lock_guard<std::mutex> qlock(s_qual_mutex);
      StdSeq::Quality qval;
      qval.quality_val(ctx.ev);
    }

    // "BAD" seeds -- leading to tracks with certain badness, see selection below.
    // Currently: |eta| < 1, no pixel hits added
    const bool print_bad_reco_seed_vector = true;
    std::vector<int> bad_seeds;

    int NT = ctx.tracks.size();
    for (int i = 0; i < NT; ++i) {
      const Track &t = ctx.tracks[i];
      auto sifh = ctx.ev->simInfoForTrack(t);
      if (!sifh.is_set()) // only take tracks with sim match
        continue;

      if (print_bad_reco_seed_vector) {
        if (sifh.n_pix_match == 0 && std::abs(t.momEta()) < 1.0f) { // catching cases where we add no pixel hits
          printf("bad cand %3d, seed %3d  pt=%6.3f, eta=% 5.3f, phi=% 5.3f -- ", i, t.label(), t.pT(), t.momEta(), t.momPhi());
          print("", sifh);
          bad_seeds.push_back(t.label());
        }
      }

    }

    if (print_bad_reco_seed_vector) {
      printf("s.SetSeedsFromIdcs({ ");
      int nbs = bad_seeds.size();
      if (nbs > 0) printf("%d", bad_seeds[0]);
      for (int bsi = 1; bsi < nbs; ++bsi)
        printf(", %d", bad_seeds[bsi]);
      printf(" })\n");
      printf("// NOTE: these are labels, sequntial ids before cleaning (as seeds are relabeled at the start).\n");
    }

    printf("\n##### END Event %d ##### HLT seeds\n", ctx.ev->evtID());
  }

  void Shell::LoopNEventsHlt(int N_events, const int wanted_algo) {
    StdSeq::Quality::s_quality_sum.quality_reset();

    int ev_first, ev_last;
    resolve_event_range(N_events, ev_first, ev_last);
    printf("Shell::LoopNEventsHlt over events [%d, %d], algo %d\n", ev_first, ev_last, wanted_algo);

    for (int ev = ev_first; ev <= ev_last; ++ev) {
      GoToEvent(ev);
      ProcessEventHlt(m_ctx, wanted_algo);
    }

    StdSeq::Quality::s_quality_sum.quality_print();
  }

  #pragma endregion Event Loops

  //===========================================================================
  #pragma region RunLSTintoPix
  //===========================================================================

  void Shell::RunLSTintoPix(EvCtx &ctx, SeedSelect_e seed_select, int selected_seed, int count) {
    const IterationConfig &itconf = Config::ItrInfo[m_it_index];
    IterationMaskIfc mask_ifc;
    ctx.ev->fill_hitmask_bool_vectors(itconf.m_track_algorithm, mask_ifc.m_mask_vector);

    const TrackerInfo &trackerInfo = Config::TrkInfo;

    ctx.tracks.clear();

    select_seeds(*ctx.ev, itconf.m_track_algorithm, seed_select, selected_seed, count, ctx.seeds);

    printf("Shell::RunLSTintoPix running over %d seeds\n", (int) ctx.seeds.size());

    // The track scorer decides the final pick of each seed, the score after the
    // backward fit and the duplicate cleaner's ranking. For this search it is
    // whatever the iteration config names; none of them was made for it.
    static std::once_flag s_scorer_warning;
    std::call_once(s_scorer_warning, [&]() {
      const SteeringParams &sp = itconf.m_steering_params[0];
      const std::string &name =
          sp.m_track_scorer_name.empty() ? itconf.m_default_track_scorer_name : sp.m_track_scorer_name;
      printf("WARNING Shell::RunLSTintoPix: track scorer '%s' -- the track score for the LST-into-pixels "
             "search needs work.\n",
             name.c_str());
    });

    for (int is = 0; auto &s : ctx.seeds) {
      // Chomp off pixels, pretending we have T5s
      // s.sortHitsByLayer();

      // Spring 2026, trouble with pickup of T5s (and especially T4 - 1 hit).
      // This was both in barrel and transition region. Also needed some for FWD (PIX only).
      // Eventually the following helps:
      // - increas inward plan pickup layer range (see CMS-phase2.cc)
      // - enable sorting of hits by layer below (to handle "swapped PS layers")
      // Some things are expected to improve.
      // - dropping of some hits from TXs (Slava is on to it)
      // - eventually finder-v2p2 will handle this better (or at least more flexibly).
      s.sortHitsByLayer();

      // print("seed-post-sort", is, s, *ctx.ev);

      // TODO (running pT5 + T5 together): this is the "chomp off the pixels,
      // pretend we have a T5" switch from the comment at the top of the loop,
      // and it is most likely how pT5s would be fed in alongside real T5s --
      // strip a pT5's pixel hits and let the backward search find them again.
      // Currently a dead const bool.
      //
      // It can stay one global switch: a T5 starts in the strips and so has no
      // pixel hits to chop, making this a no-op for them. What actually blocks
      // the combined run is elsewhere -- hlt_seed_preselect() returns a single
      // kind and everything downstream takes [first, first+count).
      //
      // Note it interacts with the hit sorting above; those comments are about
      // exactly this pickup being fragile in barrel and transition.
      const bool clear_out_pixel_hits = s_hlt_chop_pixels;
      if (clear_out_pixel_hits) {
        std::vector<HitOnTrack> ohits;
        s.swapOutAndResetHits(ohits);
        auto &rec = s_chopped_hits[s.label()];
        for (auto &oh : ohits) {
          if (trackerInfo[oh.layer].is_pixel()) {
            if (oh.index >= 0) rec.push_back(oh);
            continue;
          }
          s.addHitIdx(oh, 0.0f);
        }
      }

      // QQQQ test hacking seed charge if it is zero
      if (s.charge() == 0) {
        auto sifh = ctx.ev->simInfoForTrack(s);
        if (sifh.is_set()) {
          int chg = ctx.ev->simTracks_[sifh.label].charge();
          printf("ZERO-CHG: setting charge to %d for seed %d\n", chg, is);
          s.setCharge(chg);
        } else {
          printf("ZERO-CHG: failed getting label from hits for seed %d\n", is);
        }
      }
      // print("seed-post-pix-removal", is, s, *ctx.ev);

      ++is;
    }

    {
      const EventOfHits &eoh = *ctx.eoh;
      const IterationMaskIfcBase &it_mask_ifc = mask_ifc;
      MkBuilder &builder = *ctx.bld;
      TrackVec &seeds = ctx.seeds;
      TrackVec &out_tracks = ctx.tracks;

      MkJob job({trackerInfo, itconf, eoh, eoh.refBeamSpot(), &it_mask_ifc});

      builder.begin_event(&job, ctx.ev, __func__);

      // Check nans in seeds.
      builder.seed_post_cleaning(seeds);

      ctx.ev->setCurrentSeedTracks(seeds);

      builder.find_tracks_load_seeds(seeds, false); // false - seeds not sorted

// ********** NOTE: FORWARD SEARCH COMMENTED OUT **********
      // builder.findTracksStandardv2p2();

      job.switch_to_backward();

      builder.compactifyHitStorageForBestCand(false, 99); // do not remove anything
      builder.backwardFit(); // prop_to_plane depends on Config::usePropToPlane

      builder.beginBkwSearch();

      // QQQQQ [ April 2026 - this is mostly obsolete now, to be traced down for V2p2 ]
      // On a barrel track, index 1, pt 9.8
      // - selectHitIndicesV1 & V2
      //   . no scale misses the hit on layer 3 (as overlap)
      //   . scale 10 gets it right
      //   . scale 100 picks up an extra hit in layer 3
      // - selectHitIndicesV2
      //   . no scale misses the hit on layer 3 (as overlap)
      //   . scale 10 gets it right
      //   . scale 100 picks up an extra hit in layer 3
      // On a track going through 3 pixel disks both fail, index 0, pt 0.8
      // builder.ref_eocc_nc().scaleErrors(100.0f);

      // print("post-bkfit-n-scale state", builder.ref_eocc()[0][0].state());

      // builder.findTracksCloneEngine(SteeringParams::IT_BkwSearch);
      builder.findTracksStandardv2p2(SteeringParams::IT_BkwSearch);

      // April 2026: For tracing we really want to keep the final track.
      // Problem: it might not be the same one as the selected one.
    #ifndef MKFIT_TRACE
      filter_candidates_func post_filter;
      post_filter = StdSeq::qfilter_nan_n_silly<TrackCand>;
      // post_filter is always at least doing nan_n_silly filter.
      builder.filter_comb_cands(post_filter, true);
    #else
      for (int i = 0; i < builder.ref_eocc_nc().size(); ++i) {
        if (builder.ref_eocc_nc().cands_in_backward_rep())
          builder.ref_eocc_nc()[i].repackCandPostBkwSearch(0);
      }
    #endif

      builder.endBkwSearch();

      builder.export_best_comb_cands(out_tracks, false /*true*/); // do not remove missing hits

      ctx.ev->candidateTracks_ = ctx.tracks; // For quality-val and tracing

      // Final fit
      builder.import_tracks(out_tracks);
      builder.fittracks();
      builder.export_tracks(ctx.ev->fitTracks_);

    #ifndef MKFIT_TRACE
      ctx.ev->resetCurrentSeedTracks();
    #endif

      builder.end_event();
    }

    // print("END LST into Pix", ctx.tracks, *ctx.ev);
  }

  #pragma endregion RunLSTintoPix

  //===========================================================================
  #pragma region Tests etc
  //===========================================================================

#ifdef MKFIT_TRACE
  namespace {
    void export_AnRun(AnRun *ar) {
      char buf[256];
      sprintf(buf, "AnRun &ar = * (AnRun*) %p;", ar);
      gROOT->ProcessLine(buf);
      printf("AnRun &ar variable is set: ");
      gROOT->ProcessLine("ar");

      gROOT->ProcessLine("#define EV ar.CTX.ev");
      gROOT->ProcessLine("#define TI ar.CTX.trk_info");
      printf("EV and TI macros are set.\n");
    }
  }
#endif

  void Shell::Test(int Nevents) {
    printf("Shell::Test running basic test via LoopNEventsHlt(), algo 4\n");
    LoopNEventsHlt(Nevents, 4);
  }

#ifdef MKFIT_TRACE
  void Shell::TestVectorSource() {

    LoopNEventsHlt(1, 4);
    event()->build_trace_maps_etc();

    AnRun *ar = new AnRun(event(), *tracker_info());
    ar->RunOldVecBased();
    export_AnRun(ar);
  }
#endif

  //===========================================================================
  // collect_events_hlt -- find tracks over an event range, hand the Events out
  //===========================================================================

  void Shell::collect_events_hlt(int ev_first, int ev_last, int n_thr, int wanted_algo,
                                 std::vector<const Event *> &out)
  {
    const int n_ev = ev_last - ev_first + 1;
    out.assign(n_ev, nullptr);

    StdSeq::Quality::s_quality_sum.quality_reset();

    if (n_thr <= 1) {
      for (int i = 0; i < n_ev; ++i) {
        GoToEvent(ev_first + i);
        ProcessEventHlt(m_ctx, wanted_algo);
#ifdef MKFIT_TRACE
        m_ctx.ev->build_trace_maps_etc();
#endif
        m_ctx.ev->printMemUsage();
        out[i] = RelinquishEvent();
      }
      StdSeq::Quality::s_quality_sum.quality_print();
      return;
    }

    // Several events in flight. Per slot: Event, EventOfHits, MkBuilder and its
    // own read handle. Shared and already thread-safe: the g_exe_ctx finder
    // pools (Pool::makeOrGet() creates on demand, so it cannot starve) and the
    // DataFile event cursor (advancePosToNextEvent() is mutex-guarded).
    //
    // NOTE: tbb::parallel_for, not TBB_PARALLEL_FOR -- the macro degrades to a
    // serial loop under TBB_DEBUG, which EVENT_RDF_TRACE turns on. That is what
    // we want *inside* an event (vectorization only), not for this loop.
    tbb::task_arena arena(n_thr);
    const int n_slots = arena.max_concurrency();

    std::vector<EvCtx> ctxs(n_slots);
    for (auto &c : ctxs) {
      c.ev  = new Event(0, Config::TrkInfo.n_layers());
      c.eoh = new EventOfHits(Config::TrkInfo);
      c.bld = new MkBuilder(Config::silent);
      c.fp  = fopen(m_in_file.c_str(), "r");
      if (c.fp == nullptr)
        throw std::runtime_error("collect_events_hlt: could not open input file for a worker slot");
    }

    // Park the shared cursor on the first wanted event; workers then claim
    // events off it one at a time, each reading through its own handle.
    m_data_file->rewind();
    if (ev_first > 1)
      m_data_file->skipNEvents(ev_first - 1);

    printf("\nShell::collect_events_hlt: %d events over %d slots (--num-thr-ev %d)\n\n",
           n_ev, n_slots, Config::numThreadsEvents);

    std::mutex read_mutex;
    int claimed = 0;
    std::atomic<double> t_in_read{0.0};   // seconds spent holding read_mutex
    const auto t_coll_beg = std::chrono::steady_clock::now();

    arena.execute([&]() {
      tbb::parallel_for(tbb::blocked_range<int>(0, n_ev, 1),
        [&](const tbb::blocked_range<int> &br) {
          const int slot = tbb::this_task_arena::current_thread_index();
          if (slot < 0 || slot >= n_slots)
            throw std::runtime_error("collect_events_hlt: task arena slot index out of range");
          EvCtx &ctx = ctxs[slot];

          for (int k = br.begin(); k != br.end(); ++k) {
            int idx;
            {
              // Claim an event and read it under one lock, so that the id we
              // assign really is the event we got. The body read is cheap
              // next to the finding, which stays fully parallel.
              const auto tb = std::chrono::steady_clock::now();
              std::lock_guard<std::mutex> rlock(read_mutex);
              idx = claimed++;
              ctx.ev->reset(ev_first + idx);
              ctx.ev->read_in(*m_data_file, ctx.fp);
              const double dt = std::chrono::duration<double>(
                  std::chrono::steady_clock::now() - tb).count();
              for (double e = t_in_read.load(); !t_in_read.compare_exchange_weak(e, e + dt); ) {}
            }

            StdSeq::loadHitsAndBeamSpot(*ctx.ev, *ctx.eoh);
            if (Config::useDeadModules)
              StdSeq::loadDeads(*ctx.eoh, m_deadvectors);

            ProcessEventHlt(ctx, wanted_algo);
#ifdef MKFIT_TRACE
            ctx.ev->build_trace_maps_etc();
#endif

            // Hand the Event over and take a fresh one for the next round.
            out[idx] = ctx.ev;
            ctx.ev = new Event(0, Config::TrkInfo.n_layers());
          }
        });
    });

    for (auto &c : ctxs) {
      delete c.ev;
      delete c.bld;
      delete c.eoh;
      if (c.fp) fclose(c.fp);
    }

    // Events are claimed in file order and stored by claim index, so `out` is
    // already ordered. Sort anyway: it is cheap and keeps the guarantee local.
    std::sort(out.begin(), out.end(),
              [](const Event *a, const Event *b) { return a->evtID() < b->evtID(); });

    const double t_coll = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_coll_beg).count();
    printf("\nShell::collect_events_hlt: %d events in %.2f s over %d slots; "
           "%.2f s cumulative inside the read lock (%.0f%% of one slot's share)\n\n",
           n_ev, t_coll, n_slots, t_in_read.load(), 100.0 * t_in_read.load() / (t_coll * n_slots));

    StdSeq::Quality::s_quality_sum.quality_print();
  }

#ifdef MKFIT_TRACE
  void Shell::TestEventSource(int Nevents, const char *prefix, int n_thr) {
    std::vector<const Event*> ev_vec;

    int ev_first, ev_last;
    resolve_event_range(Nevents, ev_first, ev_last);
    if (n_thr < 0)
      n_thr = Config::numThreadsEvents; // --num-thr-ev
    printf("\n######### Shell::TestEventSource() -- running over events [%d, %d], n_thr = %d.\n\n",
           ev_first, ev_last, n_thr);

    // Track finding, optionally several events in flight. The RDF analysis
    // below runs afterwards and is deliberately left SINGLE THREADED.
    //
    // ROOT implicit MT works now (it used to deadlock -- the data sources were
    // not chaining to RDataSource::SetNSlots(), so RSlotStack got zero slots and
    // spun forever), and Shell::EnableRdfMT() turns it on. It is simply not
    // worth it here: one RDF entry = one Event, so the loop is ~1 s over ~100
    // entries, and IMT made it slower (1.02 s -> 1.56 s at 8 threads). The run
    // is dominated by fixed cost instead -- ~11 s ROOT/geometry startup plus
    // ~16 s of cling JIT for the string Defines/Filters.
    //
    // Revisit when the analysis gets heavier or interactive: in a live session
    // the fixed cost is paid once and the event loop is what repeats. Then also
    // consider handing back finer entry ranges than slots, since EventSource
    // currently returns exactly nSlots ranges and TBB has nothing to steal if
    // some events are much more expensive than others.
    collect_events_hlt(ev_first, ev_last, n_thr, 4, ev_vec);

    AnRun *ar = new AnRun(*tracker_info());
    ar->SetPrefix(prefix); // output goes to <prefix>.root and <prefix>.txt
    ar->SetupRdfEvent(ev_vec); // ev_vec is swapped, AnRun assumes ownership

    // Figuring out a good structure for basic columns, derived rdfs, etc.

    ar->RunBasicSeedCandCheck();
    // ar->RunMetaVsSeedDuplicateCheck();
    // ar->Run_T5_vs_pT5_AsSeeds_DuplicateCount();
    // ar->Run_T5s_into_Pix();
    ar->Run_Stage2_RootState_QualityCheck();
    ar->Run_Stage2_RootState_Covariance_Check();

    ar->DrawCanvasGroups();
    ar->WriteCanvasGroupsToFile();

    export_AnRun(ar);
  }
#endif

  #pragma endregion Tests etc

  //===========================================================================
  #pragma region Standard forward search, traced
  //===========================================================================

  // Per-event processing for the ordinary OUTWARD forward search, as opposed to
  // the inward T5-into-pixels search the *Hlt functions above drive.
  //
  // ProcessEvent() already does all of it and already traces -- on the v2p2
  // path findTracksStandardv2p2() creates the TrCandMeta / TrCandStage /
  // TrCandState nodes -- so this is just the per-event wrapper the collector
  // needs, mirroring ProcessEventHlt().
  //
  // Deliberately a parallel implementation rather than a refactor of the Hlt
  // path. The differences that are real, not incidental:
  //  - no seed-kind classification: we want every seed of the iteration, so
  //    SS_UseAll, not a contiguous [first, count) block picked by index;
  //  - hence no relabelSeedTracksSequentially() -- nothing indexes seeds here;
  //  - seed cleaning, backward fit and backward search follow the
  //    IterationConfig as in a real job, instead of the Hlt path's hand-picked
  //    selection.
  //===========================================================================
  // Sim-seeded forward search: start at layer 0 and go OUTWARD through the
  // pixels.
  //
  // Why this exists: the normal forward search is seeded in the pixels, so it
  // never scans the pixel barrel -- yet that is exactly where the inward search
  // shows the q covariance 22x too small. Without an outward pass through the
  // same layers there is no like-for-like comparison. Seeding from the sim
  // tracks gives one: the state is true by construction, so anything the
  // covariance does wrong afterwards is the code's, not the seed's.
  //
  // The seeds are "completely fake, well, completely true": state copied from
  // the sim track, hits taken as that track's innermost n_seed_hits.
  //===========================================================================

  int Shell::MakeSimSeeds(EvCtx &ctx, int n_seed_hits, float pt_min) {
    ctx.seeds.clear();
    for (int is = 0; is < (int) ctx.ev->simTracks_.size(); ++is) {
      const Track &st = ctx.ev->simTracks_[is];
      if (st.pT() < pt_min) continue;
      // innermost hits, by layer number -- the sim hit list is not guaranteed
      // ordered, and for an outward search the seed must sit at the INSIDE.
      std::vector<HitOnTrack> hots;
      for (int i = 0; i < st.nTotalHits(); ++i) {
        const HitOnTrack hot = st.getHitOnTrack(i);
        if (hot.index >= 0 && hot.layer >= 0) hots.push_back(hot);
      }
      if ((int) hots.size() < n_seed_hits) continue;
      std::sort(hots.begin(), hots.end(),
                [](const HitOnTrack &a, const HitOnTrack &b) { return a.layer < b.layer; });
      // require the seed to actually start in the innermost pixel layer, or the
      // outward pass does not cover what it is meant to cover
      if (hots[0].layer != 0) continue;
      Track seed(st.state(), 0.0f, is, 0, nullptr);
      // The sim track's own covariance is a PLACEHOLDER -- err(i,i) = value^2,
      // i.e. a 100 % relative error, and singular at the origin.
      // Left as is it gives sigma_q ~ 4.5 cm at layer 1 and the first two
      // layers measure nothing but the placeholder washing out. Replace it with
      // a small, non-degenerate, deliberately TIGHT covariance: the state is
      // true, so the only honest prior is a narrow one, and the filter forgets
      // a prior within 3-4 hits anyway (measured: the bkfit error-scale scan is
      // flat over 1 .. 1e4).
      {
        SMatrixSym66 &e = seed.errors_nc();
        for (int a = 0; a < 6; ++a) for (int b = 0; b <= a; ++b) e(a,b) = 0.0f;
        const float sp = 1.0e-3f;             // 10 um on each position
        const float sa = 1.0e-4f;             // 0.1 mrad on each angle
        e(0,0) = e(1,1) = e(2,2) = sp*sp;
        e(3,3) = 1.0e-4f * st.invpT() * st.invpT();   // 1 % on 1/pT
        e(4,4) = e(5,5) = sa*sa;
      }
      for (int i = 0; i < n_seed_hits; ++i) seed.addHitIdx(hots[i], 0.0f);
      seed.setLabel(is);
      ctx.seeds.push_back(seed);
    }
    printf("Shell::MakeSimSeeds: %d sim-seeded tracks (pt > %.2f, starting at layer 0, "
           "%d seed hits)\n", (int) ctx.seeds.size(), pt_min, n_seed_hits);
    return (int) ctx.seeds.size();
  }

  void Shell::ProcessEventSimSeeded(EvCtx &ctx, int n_seed_hits, float pt_min) {
    if (MakeSimSeeds(ctx, n_seed_hits, pt_min) == 0) return;
    // The forward pickup is a FIXED plan index (make_iterator does
    // `m_cur_index = m_fwd_search_pickup`), NOT something derived from the
    // seed's last hit. The barrel default is 2, so a seed ending at layer 0
    // would have the search begin at layer 2 and silently skip layer 1. Move
    // the pickup to just past the seed's last layer.
    // A candidate is picked up at the plan entry matching its seed's LAST hit
    // layer, and the pickup entry is itself search-free (is_pickup_only() is
    // `m_cur_index == m_fwd_search_pickup`). Our seeds end at layer 0, so the
    // pickup must be plan index 0 -- with anything else nothing ever matches
    // and the search produces literally zero layer-searches.
    // Every region's plan starts with fill_plan(0, ...), so index 0 is layer 0
    // in all of them; set it everywhere, not just the barrel, since a sim seed
    // can land in any region.
    IterationConfig &ic = Config::ItrInfo[m_it_index];
    int saved[TrackerInfo::Reg_Count];
    for (int r = 0; r < TrackerInfo::Reg_Count; ++r) {
      saved[r] = ic.m_steering_params[r].m_fwd_search_pickup;
      ic.m_steering_params[r].m_fwd_search_pickup = n_seed_hits - 1;
    }
    ctx.tracks.clear();
    ProcessEvent(ctx, SS_PreSet);
    for (int r = 0; r < TrackerInfo::Reg_Count; ++r)
      ic.m_steering_params[r].m_fwd_search_pickup = saved[r];
  }

  void Shell::ProcessEventStd(EvCtx &ctx) {
    printf("\n##### BEG Event %d ##### standard forward search\n\n", ctx.ev->evtID());

    // Sim-truth cleanup, as in the Hlt path. The trace's mc_match comes off
    // these labels, so it matters for any pull or efficiency study.
    ctx.ev->filterOutMislabeledHitsInSimTracks();

    ProcessEvent(ctx, SS_UseAll);

    // Publish the found tracks into the Event, as RunLSTintoPix() does. Both
    // quality_val() and the tracing need it: TrCandMeta::cand is an index into
    // candidateTracks_, assigned in MkBuilder::export_best_comb_cands(), and
    // ProcessEvent() exports into the caller's vector (ctx.tracks) only.
    ctx.ev->candidateTracks_ = ctx.tracks;

    {
      // quality_val() ends in add_to_quality_sum(), which writes the static
      // Quality::s_quality_sum -- guard it so several events can be in flight.
      // Also keeps each event's quality printout in one piece.
      static std::mutex s_qual_mutex;
      std::lock_guard<std::mutex> qlock(s_qual_mutex);
      StdSeq::Quality qval;
      qval.quality_val(ctx.ev);
    }

    printf("\n##### END Event %d ##### standard forward search\n\n", ctx.ev->evtID());
  }

  // Near-copy of collect_events_hlt() with ProcessEventStd() in place of
  // ProcessEventHlt(). See the comments there for why tbb::parallel_for and not
  // TBB_PARALLEL_FOR, and for what is and is not shared between slots.
  void Shell::collect_events_std(int ev_first, int ev_last, int n_thr,
                                 std::vector<const Event *> &out)
  {
    const int n_ev = ev_last - ev_first + 1;
    out.assign(n_ev, nullptr);

    StdSeq::Quality::s_quality_sum.quality_reset();

    if (n_thr <= 1) {
      for (int i = 0; i < n_ev; ++i) {
        GoToEvent(ev_first + i);
        ProcessEventStd(m_ctx);
#ifdef MKFIT_TRACE
        m_ctx.ev->build_trace_maps_etc();
#endif
        out[i] = RelinquishEvent();
      }
      StdSeq::Quality::s_quality_sum.quality_print();
      return;
    }

    tbb::task_arena arena(n_thr);
    const int n_slots = arena.max_concurrency();

    std::vector<EvCtx> ctxs(n_slots);
    for (auto &c : ctxs) {
      c.ev  = new Event(0, Config::TrkInfo.n_layers());
      c.eoh = new EventOfHits(Config::TrkInfo);
      c.bld = new MkBuilder(Config::silent);
      c.fp  = fopen(m_in_file.c_str(), "r");
      if (c.fp == nullptr)
        throw std::runtime_error("collect_events_std: could not open input file for a worker slot");
    }

    m_data_file->rewind();
    if (ev_first > 1)
      m_data_file->skipNEvents(ev_first - 1);

    printf("\nShell::collect_events_std: %d events over %d slots (--num-thr-ev %d)\n\n",
           n_ev, n_slots, Config::numThreadsEvents);

    std::mutex read_mutex;
    int claimed = 0;
    const auto t_coll_beg = std::chrono::steady_clock::now();

    arena.execute([&]() {
      tbb::parallel_for(tbb::blocked_range<int>(0, n_ev, 1),
        [&](const tbb::blocked_range<int> &br) {
          const int slot = tbb::this_task_arena::current_thread_index();
          if (slot < 0 || slot >= n_slots)
            throw std::runtime_error("collect_events_std: task arena slot index out of range");
          EvCtx &ctx = ctxs[slot];

          for (int k = br.begin(); k != br.end(); ++k) {
            int idx;
            {
              std::lock_guard<std::mutex> rlock(read_mutex);
              idx = claimed++;
              ctx.ev->reset(ev_first + idx);
              ctx.ev->read_in(*m_data_file, ctx.fp);
            }

            StdSeq::loadHitsAndBeamSpot(*ctx.ev, *ctx.eoh);
            if (Config::useDeadModules)
              StdSeq::loadDeads(*ctx.eoh, m_deadvectors);

            ProcessEventStd(ctx);
#ifdef MKFIT_TRACE
            ctx.ev->build_trace_maps_etc();
#endif

            out[idx] = ctx.ev;
            ctx.ev = new Event(0, Config::TrkInfo.n_layers());
          }
        });
    });

    for (auto &c : ctxs) {
      delete c.ev;
      delete c.bld;
      delete c.eoh;
      if (c.fp) fclose(c.fp);
    }

    std::sort(out.begin(), out.end(),
              [](const Event *a, const Event *b) { return a->evtID() < b->evtID(); });

    const double t_coll = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_coll_beg).count();
    printf("\nShell::collect_events_std: %d events in %.2f s over %d slots\n\n",
           n_ev, t_coll, n_slots);

    StdSeq::Quality::s_quality_sum.quality_print();
  }

  // Driver: traced forward search over the configured event range, Events then
  // handed to AnRun.
  //
  // Intended for a sample with ordinary CMSSW iteration seeds, e.g.
  //   ./mkFit --input-file .../ntuple_ttbar_PU-p2.bin --geom CMS-phase2
  //     --seed-input cmssw --seed-cleaning n2 --build-mimi --build-mimi-v2p2
  //     --num-iters-cmssw 1 --num-events 10 --shell
  //     --shell-command "gROOT->SetBatch(kTRUE)"
  //     --shell-command 's.TraceFwdSearch(10, "fwd", 4)'
  // (one line in the shell; the breaks above are only for readability)
  //
  // The canned AnRun analyses used by the Hlt driver are mostly about the
  // backward-search stage (stage 2) and are not run here. For the forward-search
  // pull work the Events are reachable from the interpreter as ar.m_ev_vec,
  // which is what the exploratory macros use.
#ifdef MKFIT_TRACE
  void Shell::TraceFwdSearch(int Nevents, const char *prefix, int n_thr) {
    std::vector<const Event*> ev_vec;

    int ev_first, ev_last;
    resolve_event_range(Nevents, ev_first, ev_last);
    if (n_thr < 0)
      n_thr = Config::numThreadsEvents;
    printf("\n######### Shell::TraceFwdSearch() -- events [%d, %d], n_thr = %d.\n\n",
           ev_first, ev_last, n_thr);

    collect_events_std(ev_first, ev_last, n_thr, ev_vec);

    AnRun *ar = new AnRun(*tracker_info());
    ar->SetPrefix(prefix);
    ar->SetupRdfEvent(ev_vec);

    ar->RunBasicSeedCandCheck();
    ar->RunSimTrackEfficiency();

    ar->DrawCanvasGroups();
    ar->WriteCanvasGroupsToFile();

    export_AnRun(ar);
  }
#endif

  #pragma endregion Standard forward search, traced

} // end namespace mkfit
