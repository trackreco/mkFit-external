#ifndef RecoTracker_MkFitCore_standalone_V2p2Diag_h
#define RecoTracker_MkFitCore_standalone_V2p2Diag_h

// Standalone diagnostics of the v2p2 finder. The core sources keep only the calls,
// under MKFIT_STANDALONE; the definitions are in V2p2Diag.cc.

#include <atomic>

namespace mkfit {

  class Event;
  class EventOfCombCandidates;

  // Did the per-layer policy actually fire? quality-val cannot answer that: a cut
  // that never fires and a cut that fires on candidates which were doomed anyway
  // both show up as "no change". Atomic because the search runs one finder per
  // thread; summed over threads and events, and printed by mkFit.cc with the
  // quality-val summary.
  //
  // The finder increments them through V2P2_COUNT() and V2P2_COUNT_ADD(), which
  // MkFinderV2p2.h defines to nothing in the CMSSW build.
  struct V2p2PolicyCounters {
    std::atomic<long> n_quadrant_skip{0};  // pull-in: coarse rz check said the layer is behind us
    std::atomic<long> n_stop_minpt{0};     // pull-in: pT below minPtCut
    std::atomic<long> n_stop_looper{0};    // pull-in: past pi/2 - 0.2 to the radial direction
    std::atomic<long> n_hole{0};           // kHitMissIdx -- a real hole, counted against the limits
    std::atomic<long> n_hot_edge{0};       // kHitEdgeIdx -- not counted
    std::atomic<long> n_hot_gap{0};        // kHitInGapIdx -- not counted
    std::atomic<long> n_stop_holes{0};     // kHitStopIdx -- out of hole budget
    std::atomic<long> n_ccand_retired{0};  // every TrackCand under it stopped
    std::atomic<long> n_nodes{0};          // in-layer tree nodes built
    std::atomic<long> n_nodes_deep{0};     // ... of those, at depth 2 or more
    std::atomic<long> n_path_taken{0};     // layers where a path was registered
    std::atomic<long> n_extra_hits{0};     // hits taken beyond the first in one layer
    std::atomic<long> n_sel_entries{0};    // competitors entering the end-of-layer selection
    std::atomic<long> n_sel_kept{0};       // ... and surviving it
    std::atomic<long> n_selections{0};     // end-of-layer selections run
    std::atomic<long> n_same_module_vetoed{0};  // extensions refused for sharing a module
    std::atomic<long> n_hole_slot_reserved{0};  // beam slots given to an outranked decliner
    std::atomic<long> n_best_short_offered{0};  // stopped candidates removed from the beam
    std::atomic<long> n_best_short_taken{0};    // ... and that became the seed's best short
    // Matriplex lane occupancy of the Kalman batches: is the expansion still
    // vectorised, or is it running ragged tails? mean = lanes / (calls * NN).
    std::atomic<long> n_kalman_calls{0};
    std::atomic<long> n_kalman_lanes{0};
    std::atomic<long> n_kalman_calls_d0{0};   // depth 0 only
    std::atomic<long> n_kalman_lanes_d0{0};
    std::atomic<long> n_late_max_cands{0};          // CombCandidate activations at InLayer::late_max_cands
    std::atomic<long> n_late_max_cands_dropped{0};  // ... TrackCands dropped by it at activation
    std::atomic<long> n_late_max_cands_long{0};     // activations kept at full width by a long step
    std::atomic<long> n_sister_hole_dropped{0};     // holes dropped for a sister-sensor sibling

    void reset();
    void print(const char *tag) const;
  };
  extern V2p2PolicyCounters g_v2p2_policy_counters;

#define V2P2_COUNT(field) (++g_v2p2_policy_counters.field)
#define V2P2_COUNT_ADD(field, n) (g_v2p2_policy_counters.field += (n))

  // Diag::final_beam_purity: each seed's beam at the end of the v2p2 forward search,
  // against truth. Called from MkFinderV2p2 (bg_record) and from
  // MkBuilder::findTracksStandardv2p2() (diag), single-threaded accumulation.
  // layer_sec: the second sub-layer of the plan entry, or -1.
  void v2p2_final_beam_bg_record(int seed, int layer, int layer_sec, float n_bg);
  class CombCandidate;
  void v2p2_final_pick_record(const CombCandidate &cc);
  void v2p2_final_beam_diag(const Event *ev, const EventOfCombCandidates &eoccs);
  void v2p2_final_beam_diag_reset();
  void v2p2_final_beam_diag_report();

}  // namespace mkfit

#endif
