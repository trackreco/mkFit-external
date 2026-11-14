// Probe-Taps.cc -- the readers of the in-core taps. Shell::SetTapFinalBeam() and
// Shell::SetTapScoreStats() switch the taps and set tap mode, Shell.h.
#include "RecoTracker/MkFitCMS/standalone/Probe/Probe.h"
#include "RecoTracker/MkFitCore/standalone/V2p2Diag.h"

#include "RecoTracker/MkFitCore/src/V2p2Config.h"

#include <cstdio>

namespace mkfit {
  namespace V2p2 = Config::V2p2;

  // The end-of-search beam purity collected while Shell::SetTapFinalBeam(true) is on.
  void probe_final_beam_report() { v2p2_final_beam_diag_report(); }

  // Mean ln(rho) and mean ln(det V) per hit taken, so the constants for
  // Shell::SetScoreTerms() are measured. The sums are collected while
  // Shell::SetTapScoreStats(true) is on.
  void probe_score_term_stats_report() {
    const long n = V2p2::ScoreStats::n_hits;
    printf("probe_score_term_stats_report: %ld hits scored; mean ln(rho) = %.4f, "
           "mean ln(det V) = %.4f\n", n,
           n ? V2p2::ScoreStats::sum_log_rho / n : 0.0,
           n ? V2p2::ScoreStats::sum_log_detv / n : 0.0);
  }

}  // namespace mkfit
