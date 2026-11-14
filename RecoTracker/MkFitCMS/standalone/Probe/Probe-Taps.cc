#include "RecoTracker/MkFitCMS/standalone/Probe/Probe.h"
#include "RecoTracker/MkFitCore/standalone/V2p2Diag.h"

#include "RecoTracker/MkFitCore/src/V2p2Config.h"

#include <cstdio>

namespace mkfit {
  namespace V2p2 = Config::V2p2;

  void probe_final_beam(bool on) {
    V2p2::Diag::final_beam_purity = on;
    v2p2_final_beam_diag_reset();
    printf("probe_final_beam: %d (end-of-search beam purity, reset)\n", (int) on);
  }
  void probe_final_beam_report() { v2p2_final_beam_diag_report(); }

  // Mean ln(rho) and mean ln(det V) per hit taken, so the constants for
  // val_score_terms() are measured. Accumulation is not thread safe; a trace build serialises the
  // in-event loops, which is where this is meant to run.
  void probe_score_term_stats(bool on) {
    if (on) {
      V2p2::ScoreStats::n_hits = 0;
      V2p2::ScoreStats::sum_log_rho = V2p2::ScoreStats::sum_log_detv = 0.0;
      V2p2::ScoreStats::accum = true;
      printf("probe_score_term_stats: accumulating.\n");
      return;
    }
    V2p2::ScoreStats::accum = false;
    const long n = V2p2::ScoreStats::n_hits;
    printf("probe_score_term_stats: %ld hits scored; mean ln(rho) = %.4f, "
           "mean ln(det V) = %.4f\n", n,
           n ? V2p2::ScoreStats::sum_log_rho / n : 0.0,
           n ? V2p2::ScoreStats::sum_log_detv / n : 0.0);
  }

}  // namespace mkfit
