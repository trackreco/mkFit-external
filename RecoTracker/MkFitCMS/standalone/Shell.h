#ifndef RecoTracker_MkFitCMS_standalone_Shell_h
#define RecoTracker_MkFitCMS_standalone_Shell_h

#include "RecoTracker/MkFitCore/interface/Hit.h"
#include "RecoTracker/MkFitCore/interface/Track.h"
#include "RecoTracker/MkFitCore/interface/DeadRegion.h"

#include <map>

class TTree;

#ifdef MKFIT_WITH_REVE
namespace ROOT::Experimental {
  class REveManager;
  class REveTrackPropagator;
}
#endif

namespace mkfit {

  class DataFile;
  class Event;
  class EventOfHits;
  class MkBuilder;
  class TrackerInfo;

  // Per-event processing context: everything the track-finding pipeline needs
  // that cannot be shared between concurrently processed events.
  //
  // Shared, and deliberately NOT in here: the MkFinder / CandCloner / MkFitter
  // pools in g_exe_ctx (thread-safe, create-on-demand), the DataFile (its event
  // cursor is mutex-guarded) and Config.
  //
  // Shell holds one of these for interactive use; a multi-threaded driver holds
  // one per event in flight.
  struct EvCtx {
    Event       *ev  = nullptr;
    EventOfHits *eoh = nullptr;
    MkBuilder   *bld = nullptr;
    FILE        *fp  = nullptr;  // own read handle, only needed when running MT

    TrackVec seeds;   // working copies; published into ev during processing
    TrackVec tracks;
  };

  class Shell {
  public:
    enum SeedSelect_e { SS_UseAll = 0, SS_Label, SS_IndexPreCleaning, SS_IndexPostCleaning, SS_PreSet };

    // HLT seeds arrive as one iteration, three kinds back to back, told apart
    // only by hit pattern: does the seed start in the pixels, does it end in
    // the strips. In file order: pT5 (T5 matched to pixel tracks), T5 (the
    // leftover, unmatched ones), pix (leftover pixel-only).
    enum HltSeedKind_e { HSK_pT5 = 0, HSK_T5, HSK_pix };

    struct HltSeedBlocks {
      int first = -1;   // first seed index of the wanted kind, -1 if there are none
      int count = 0;    // how many of them
      int n_pT5 = 0, n_T5 = 0, n_pix = 0;   // totals, all three kinds
      bool contiguous = true;  // false if the wanted kind was not one solid block
    };

    static const char* hlt_seed_kind_name(HltSeedKind_e k);

    // Runtime selectors for ProcessEventHlt(), so the pT5-chopped control can
    // be run against plain T5 in one binary. A pT5's pixel hits were FOUND by
    // the upstream reconstruction, so chopping them off gives a denominator of
    // hits known to be findable -- which a plain T5 does not.
    static HltSeedKind_e s_hlt_seed_kind;
    static bool          s_hlt_chop_pixels;

    // What the chop actually removed, keyed by the seed's (sequential) label.
    // This is the exact denominator for "did the inward search put them back":
    // these hits were found by the upstream reconstruction, so they are known
    // findable, and no truth matching is involved in the comparison.
    static std::map<int, std::vector<HitOnTrack>> s_chopped_hits;

    // Copy the seeds of one algo into `out`, honouring the SS_* selector.
    // Seeds are taken to be grouped by algo, so the scan stops once it leaves
    // that block. SS_PreSet leaves `out` untouched.
    static void select_seeds(const Event &ev, int algo, SeedSelect_e seed_select,
                             int selected_seed, int count, TrackVec &out);

    // Classify the seeds of one algo by hit pattern and return the block of the
    // wanted kind. Everything downstream takes [first, first+count), so a false
    // `contiguous` means the file layout changed and that range is wrong.
    static HltSeedBlocks hlt_seed_preselect(const Event &ev, const TrackerInfo &ti,
                                            int algo, HltSeedKind_e want);

    Shell();
    Shell(std::vector<DeadVec> &dv, const std::string &in_file, int start_ev);
    ~Shell();
    void RunShell(const std::vector<std::string> &commands);

    void Status();

    void GoToEvent(int eid);
    void NextEvent(int skip = 1);

    // Event range used by the Loop / Test drivers below.
    // 1-based and inclusive: event 1 is the first event in the file, as on the
    // command line. Defaults come from --start-event and --num-events; a
    // positive count argument to a driver overrides the length for that call.
    void SetEventRangeFirstLast(int first, int last);
    void SetEventRangeBegCnt(int beg, int count);
    int EventRangeFirst() const { return m_ev_first; }
    int EventRangeLast() const { return m_ev_last; }
    int EventsInFile() const { return m_evs_in_file; }
    // Pipeline functions take the context explicitly so several events can be
    // in flight at once. The no-context overloads run on the shell's own
    // context -- that is what you want from the ROOT prompt.
    void ProcessEvent(EvCtx &ctx, SeedSelect_e seed_select = SS_UseAll, int selected_seed = -1, int count = 1);
    void ProcessEvent(SeedSelect_e seed_select = SS_UseAll, int selected_seed = -1, int count = 1)
      { ProcessEvent(m_ctx, seed_select, selected_seed, count); }
    Event* RelinquishEvent();

    void SelectIterationIndex(int itidx);
    void SelectIterationAlgo(int algo);
    void PrintIterations();

    // ROOT implicit MT for the RDF analysis. n_thr <= 0 disables it.
    // NOTE: mkFit.cc caps total TBB parallelism at Config::numThreadsFinder via
    // tbb::global_control, so run with --num-thr >= n_thr or this starves.
    void EnableRdfMT(int n_thr);

    bool GetDebug() const;
    void SetDebug(bool b);
    void SetCleanSeeds(bool b);
    void SetBackwardFit(bool b);
    void SetBackwardSearch(bool b);
    void SetRemoveDuplicates(bool b);
    void SetUseDeadModules(bool b);
    void SetUseV2p2(bool b);
    // Backward search with findTracksStandardv2p2 when v2p2 is on, as CMSSW's
    // run_OneIteration() does; the default keeps findTracksCloneEngine.
    void SetBkwSearchV2p2(bool b) { m_bkw_search_v2p2 = b; }

    Event *event() { return m_ctx.ev; }
    EventOfHits *eoh() { return m_ctx.eoh; }
    MkBuilder *builder() { return m_ctx.bld; }
    TrackerInfo *tracker_info();

    EvCtx &ctx() { return m_ctx; }

    const TrackVec &seeds() const { return m_ctx.seeds; }
    const TrackVec &tracks() const { return m_ctx.tracks; }

    void SetSeedsFromIdcs(std::vector<int> idcs);

    // Build "fake but true" seeds from the SIM tracks: the seed state is the
    // sim track's own state and its hits are that track's innermost ones. Run
    // with SS_PreSet, this lets the FORWARD search start at layer 0 and go
    // outward through the pixel barrel -- the one region the normal forward
    // search never scans, because that is where its real seeds already are.
    // n_seed_hits: how many innermost hits to hand the seed (>=1).
    // Returns the number of seeds built.
    int MakeSimSeeds(EvCtx &ctx, int n_seed_hits = 1, float pt_min = 0.5f);
    int MakeSimSeeds(int n_seed_hits = 1, float pt_min = 0.5f)
      { return MakeSimSeeds(m_ctx, n_seed_hits, pt_min); }
    void ProcessEventSimSeeded(EvCtx &ctx, int n_seed_hits = 1, float pt_min = 0.5f);
    void ProcessEventSimSeeded(int n_seed_hits = 1, float pt_min = 0.5f)
      { ProcessEventSimSeeded(m_ctx, n_seed_hits, pt_min); }

    // --------------------------------------------------------
    // Analysis helpers

    int LabelFromHits(Track &t, bool replace, float good_frac);
    void FillByLabelMaps_CkfBase();

    bool CheckMkFitLayerPlanVsReferenceHits(const Track &mkft, const Track &reft, const std::string &name);

    // --------------------------------------------------------
    // Analysis drivers / main functions / Comparators

    void Compare();

    // --------------------------------------------------------
    // Seed study prototype
    using seed_selector_cf = bool(const Track &);
    using seed_selector_func = std::function<seed_selector_cf>;

    void StudySimAndSeeds(bool report_lost_seeds=true);
    void PreSelectSeeds(int iter_idx, seed_selector_func selector = [](const Track&) {return true;});

    void FindInterestingSimTracks();

    void WriteSimTree();
    void ReadSimTree();

    // --------------------------------------------------------
    // Low-level checks
    TTree* CheckHitVsModulePosition();

    // --------------------------------------------------------
    // Visualization stuff
#ifdef MKFIT_WITH_REVE
    void ReveInit();
    void ShowTracker(int lay_first, int lay_last);
    void ShowSimTrack(int sim_idx);

    ROOT::Experimental::REveManager& EveMgr() { return *m_reve_mgr; }
#endif

    // --------------------------------------------------------
    // Experimental phase2 / LST stuff, in Shell-LST.cc
    void RunLSTintoPix(EvCtx &ctx, SeedSelect_e seed_select = SS_UseAll, int selected_seed = -1, int count = 1);
    void RunLSTintoPix(SeedSelect_e seed_select = SS_UseAll, int selected_seed = -1, int count = 1)
      { RunLSTintoPix(m_ctx, seed_select, selected_seed, count); }

    void LoopNEvents(int N_events);

    // One event: seed classification + RunLSTintoPix + quality-val + reports.
    // Assumes ctx is already loaded with the wanted event.
    void ProcessEventHlt(EvCtx &ctx, const int wanted_algo = 4);
    void ProcessEventHlt(const int wanted_algo = 4) { ProcessEventHlt(m_ctx, wanted_algo); }
    void LoopNEventsHlt(int N_events = -1, const int wanted_algo = 4);

    // ---- Standard (non-HLT) forward search, traced ------------------------
    // The Hlt drivers above run the inward T5-into-pixels search on LST seeds.
    // These run the ordinary outward forward search on the iteration's own
    // seeds. ProcessEvent() already produces a full trace on the v2p2 path, so
    // all that was missing was a driver. Written as a deliberate parallel of
    // the *Hlt versions rather than factored with them -- commonalities to be
    // pulled out once we know which parts really are common.
    // Motivation: pre-selection pulls are only the visible end. Kalman, chi2
    // and scoring are where a bad covariance actually costs tracks, and that
    // wants the outward direction as well as the inward one.
    void ProcessEventStd(EvCtx &ctx);
    void ProcessEventStd() { ProcessEventStd(m_ctx); }
#ifdef MKFIT_TRACE
    void TraceFwdSearch(int Nevents = -1, const char *prefix = "mkfit-fwd", int n_thr = 1);
#endif

    // Current default processing. N_events <= 0 means "use the configured
    // event range", see SetEventRangeFirstLast() / SetEventRangeBegCnt().
    // NOTE on naming: Test() really is a test; TestVectorSource() and
    // TestEventSource() are tracing drivers and want to be Trace*.
    void Test(int Nevents = -1);
#ifdef MKFIT_TRACE
    void TestVectorSource();
    // n_thr > 1 runs that many events in flight, each in its own EvCtx; the RDF
    // analysis afterwards is unchanged. n_thr < 0 takes --num-thr-ev.
    void TestEventSource(int Nevents = -1, const char *prefix = "mkfit", int n_thr = 1);
#endif

  protected:
    int select_seeds_for_algo(int algo, TrackVec &seeds);

    // Find tracks over [ev_first, ev_last] and hand the processed Events out,
    // in event order. Serial for n_thr <= 1.
    void collect_events_hlt(int ev_first, int ev_last, int n_thr, int wanted_algo,
                            std::vector<const Event *> &out);
    void collect_events_std(int ev_first, int ev_last, int n_thr,
                            std::vector<const Event *> &out);

    // Effective range for one driver call: count > 0 overrides the length of
    // the configured range, starting from its first event.
    void resolve_event_range(int count, int &first, int &last) const;

  private:
    std::vector<DeadVec> &m_deadvectors;
    std::string m_in_file;             // kept so MT slots can open their own handles
    DataFile *m_data_file = nullptr;   // shared: event cursor is mutex-guarded
    EvCtx m_ctx;                       // the interactive / single-threaded context
    int m_evs_in_file = -1;
    int m_ev_first = 1;
    int m_ev_last = 1;
    int m_it_index = 0;
    bool m_clean_seeds = true;
    bool m_backward_fit = true;
    // Gates the backward SEARCH only; the backward FIT is m_backward_fit above.
    // Initialised from Config::backwardSearch (--no-backward-search), which this
    // path used to ignore -- see ProcessEvent().
    bool m_backward_search = true;
    bool m_bkw_search_v2p2 = false;
    bool m_remove_duplicates = true;

    using map_t = std::map<int, Track *>;
    using map_i = map_t::iterator;

    std::map<int, Track *> m_ckf_map, m_sim_map, m_seed_map, m_mkf_map;

#ifdef MKFIT_WITH_REVE
    ROOT::Experimental::REveManager *m_reve_mgr = nullptr;
    ROOT::Experimental::REveTrackPropagator *m_reve_track_prop = nullptr;
#endif

  };

}  // namespace mkfit

#endif
