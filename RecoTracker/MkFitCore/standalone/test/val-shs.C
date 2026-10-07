// Validate the per-sim-hit truth states (Event::simHitStates_, --write-sim-hit-states).
//
// The check is the distance from each truth state to the rec hit it belongs to.
// In the PIXELS that must come out at about one pixel pitch, which happens only
// if the two really are the same hit; in the STRIPS it is legitimately larger,
// because the sim hit is the true crossing while the rec hit is the cluster
// centroid, so they differ along the strip by up to its half length.
//
//   val_shs_reset();  val_shs_event(s.event()); ... ; val_shs_report();

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {
  struct ShsAcc {
    std::vector<double> d_pix, d_strb, d_stre;
    long n_hits = 0, n_no_mcid = 0, n_invalid = 0, n_valid = 0;
    long n_states = 0, n_hitinfo = 0;
  };
  ShsAcc g_shs;

  double shs_q(std::vector<double> &v, double f) {
    if (v.empty()) return 0;
    size_t i = std::min(v.size() - 1, (size_t)(f * v.size()));
    std::nth_element(v.begin(), v.begin() + i, v.end());
    return v[i];
  }
  void shs_row(const char *name, std::vector<double> &v) {
    if (v.empty()) { printf("  %-18s      (none)\n", name); return; }
    printf("  %-18s n=%-8zu p50 %8.4f  p90 %8.4f  p99 %8.4f  max %8.4f  [cm]\n",
           name, v.size(), shs_q(v, 0.50), shs_q(v, 0.90), shs_q(v, 0.99), shs_q(v, 1.0));
  }
}  // namespace

void val_shs_reset() { g_shs = ShsAcc(); }

void val_shs_event(mkfit::Event *ev) {
  if (ev == nullptr) { printf("val_shs_event: null event\n"); return; }
  const mkfit::TrackerInfo &ti = mkfit::Config::TrkInfo;
  g_shs.n_states += (long)ev->simHitStates_.size();
  g_shs.n_hitinfo += (long)ev->simHitsInfo_.size();
  if (ev->simHitStates_.empty()) return;

  for (int l = 0; l < (int)ev->layerHits_.size(); ++l) {
    const bool is_pix = ti[l].is_pixel();
    const bool is_brl = ti[l].is_barrel();
    for (const auto &h : ev->layerHits_[l]) {
      ++g_shs.n_hits;
      const int mcid = h.mcHitID();
      if (mcid < 0 || mcid >= (int)ev->simHitStates_.size()) { ++g_shs.n_no_mcid; continue; }
      const mkfit::SimHitState &s = ev->simHitStates_[mcid];
      if (!s.is_valid()) { ++g_shs.n_invalid; continue; }
      ++g_shs.n_valid;
      const double d = std::sqrt((s.x() - h.x()) * (s.x() - h.x()) + (s.y() - h.y()) * (s.y() - h.y()) +
                                 (s.z() - h.z()) * (s.z() - h.z()));
      if (is_pix) g_shs.d_pix.push_back(d);
      else if (is_brl) g_shs.d_strb.push_back(d);
      else g_shs.d_stre.push_back(d);
    }
  }
}

void val_shs_report() {
  printf("\n=== sim-hit truth states ===\n");
  printf("  states %ld, simHitsInfo %ld  -- %s\n", g_shs.n_states, g_shs.n_hitinfo,
         g_shs.n_states == g_shs.n_hitinfo ? "parallel, as required" : "*** SIZES DISAGREE ***");
  printf("  rec hits %ld: no mcHitID %ld, state invalid (zero momentum) %ld, usable %ld (%.1f %%)\n",
         g_shs.n_hits, g_shs.n_no_mcid, g_shs.n_invalid, g_shs.n_valid,
         100.0 * g_shs.n_valid / std::max(1L, g_shs.n_hits));
  printf("  distance from the truth state to its own rec hit:\n");
  shs_row("pixel", g_shs.d_pix);
  shs_row("strip barrel", g_shs.d_strb);
  shs_row("strip endcap", g_shs.d_stre);
  printf("  (pixel p50 should be about one pixel pitch, 15-20 um; the strip rows are\n"
         "   larger by construction -- true crossing against cluster centroid.)\n");
}
