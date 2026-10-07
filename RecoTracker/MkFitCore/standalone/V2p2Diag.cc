// Standalone diagnostics of the v2p2 finder, moved out of the core sources: the
// policy counters, the dominance count of the end-of-layer selection, the MkBins
// surface reference of dq_track, and the end-of-search beam purity. See V2p2Diag.h.

#include "RecoTracker/MkFitCore/standalone/V2p2Diag.h"

#include "RecoTracker/MkFitCore/interface/IterationConfig.h"
#include "RecoTracker/MkFitCore/interface/TrackStructures.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/MkFitCore/src/MkBins.h"
#include "RecoTracker/MkFitCore/src/MkFinderV2p2.h"
#include "RecoTracker/MkFitCore/src/V2p2Config.h"
#include "RecoTracker/MkFitCore/standalone/ConfigStandalone.h"
#include "RecoTracker/MkFitCore/standalone/Event.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <mutex>
#include <vector>

namespace mkfit {

  using namespace Config::V2p2;

  // Per-layer policy counters, see MkFinderV2p2.h.
  V2p2PolicyCounters g_v2p2_policy_counters;

  void V2p2PolicyCounters::reset() {
    n_quadrant_skip = 0; n_stop_minpt = 0; n_stop_looper = 0;
    n_wsr_inside = 0; n_wsr_edge = 0; n_wsr_outside = 0; n_wsr_in_gap = 0;
    n_layer_skipped = 0;
    n_hole = 0; n_hot_edge = 0; n_hot_gap = 0; n_stop_holes = 0;
    n_ccand_retired = 0;
    n_nodes = 0; n_nodes_deep = 0; n_path_taken = 0; n_extra_hits = 0;
    n_sel_entries = 0; n_sel_kept = 0; n_selections = 0;
    n_same_module_vetoed = 0;
    n_hole_slot_reserved = 0; n_best_short_offered = 0; n_best_short_taken = 0;
    n_kalman_calls = 0; n_kalman_lanes = 0; n_kalman_calls_d0 = 0; n_kalman_lanes_d0 = 0;
    n_arena_layers = 0; n_arena_hw_sum = 0; n_arena_hw_max = 0; n_early_selections = 0;
    n_late_max_cands = 0; n_late_max_cands_dropped = 0; n_late_max_cands_long = 0;
    n_sister_hole_dropped = 0;
    for (int l = 0; l < k_dom_layers; ++l) {
      n_dom_sel[l] = 0; n_dom_kept[l] = 0; n_dom_sel_full[l] = 0; n_dom_kept_full[l] = 0;
      n_dom_sub[l] = 0; n_dom_hole[l] = 0; n_dom_sister[l] = 0;
      n_dom_sub_full[l] = 0; n_dom_hole_full[l] = 0; n_dom_sister_full[l] = 0;
    }
  }

  void V2p2PolicyCounters::print_dominance() const {
    // Phase-2 layer groups; the OT pairs split into their first (even) and second
    // (odd) sub-layer, the one where a sister hit follows.
    struct G { const char *name; int lo, hi, parity; };  // parity -1: all layers
    const G groups[] = {{"PixB 0-3", 0, 3, -1},           {"TBPS 1st (4,6,8)", 4, 9, 0},
                        {"TBPS 2nd (5,7,9)", 4, 9, 1},    {"TB2S 1st", 10, 15, 0},
                        {"TB2S 2nd", 10, 15, 1},          {"FPix +-", 16, 27, -1},
                        {"TEDD 1st", 28, 37, 0},          {"TEDD 2nd", 28, 37, 1}};
    printf("  dominated beam slots after the end-of-layer selection (all / in full selections):\n"
           "    %-18s %9s %10s %7s %7s %7s %10s %10s %7s %7s %7s\n", "layers", "sel", "kept", "sub%",
           "hole%", "sist%", "sel full", "kept full", "sub%", "hole%", "sist%");
    for (const G &g : groups) {
      long a[10] = {0};
      for (int l = 0; l < k_dom_layers; ++l) {
        int lm = l;
        if (g.lo == 16) lm = (l >= 38 && l <= 49) ? l - 22 : l;   // FPix- onto FPix+
        if (g.lo == 28) lm = (l >= 50 && l <= 59) ? l - 22 : l;   // TEDD- onto TEDD+
        if (lm < g.lo || lm > g.hi || (g.parity >= 0 && (lm & 1) != g.parity))
          continue;
        a[0] += n_dom_sel[l]; a[1] += n_dom_kept[l]; a[2] += n_dom_sub[l]; a[3] += n_dom_hole[l];
        a[4] += n_dom_sister[l]; a[5] += n_dom_sel_full[l]; a[6] += n_dom_kept_full[l];
        a[7] += n_dom_sub_full[l]; a[8] += n_dom_hole_full[l]; a[9] += n_dom_sister_full[l];
      }
      auto pc = [](long x, long n) { return n > 0 ? 100.0 * x / n : 0.0; };
      printf("    %-18s %9ld %10ld %7.1f %7.1f %7.1f %10ld %10ld %7.1f %7.1f %7.1f\n", g.name, a[0], a[1],
             pc(a[2], a[1]), pc(a[3], a[1]), pc(a[4], a[1]), a[5], a[6], pc(a[7], a[6]), pc(a[8], a[6]),
             pc(a[9], a[6]));
    }
  }

  void V2p2PolicyCounters::print(const char *tag) const {
    const long n_wsr = n_wsr_inside + n_wsr_edge + n_wsr_outside;
    // Denominator is layer searches that got as far as a WSR verdict, i.e. after
    // the pull-in skips. Those are reported separately because they are a
    // different population: a candidate skipped at pull-in never became one.
    const double f = n_wsr > 0 ? 100.0 / n_wsr : 0.0;
    printf("MkFinderV2p2 layer policy (%s):\n"
           "  pull-in    : rz-quadrant skip %ld, stop-minPt %ld, stop-looper %ld\n"
           "  WSR of %ld : inside %ld (%.1f%%), edge %ld (%.1f%%), outside %ld (%.1f%%), in-gap %ld\n"
           "  no-hit HoT : hole %ld, edge %ld, gap %ld, stop-out-of-holes %ld; layers skipped %ld\n"
           "  retired    : %ld CombCandidates\n"
           "  in-layer   : %ld tree nodes (%ld at depth >= 2), %ld paths taken, %ld extra hits\n"
           "  selection  : %ld of them, %ld competitors -> %ld kept (%.2f -> %.2f per seed)\n"
           "  extra hits : %ld same-module extensions vetoed\n"
           "  hole slots : %ld reserved for an outranked decliner\n"
           "  best-short : %ld stopped cands left the beam, %ld became the seed's best short\n"
           "  Mplex lanes: depth 0 %.2f of %d over %ld calls; deeper %.2f over %ld calls\n"
           "  arena      : high water per layer %.1f nodes mean, %ld max; %ld of %ld selections early\n"
           "  late max_cands: %ld activations at it, %ld TrackCands dropped, %ld kept wide by a long step\n",
           tag,
           n_quadrant_skip.load(), n_stop_minpt.load(), n_stop_looper.load(),
           n_wsr, n_wsr_inside.load(), f * n_wsr_inside, n_wsr_edge.load(), f * n_wsr_edge,
           n_wsr_outside.load(), f * n_wsr_outside, n_wsr_in_gap.load(),
           n_hole.load(), n_hot_edge.load(), n_hot_gap.load(), n_stop_holes.load(),
           n_layer_skipped.load(),
           n_ccand_retired.load(),
           n_nodes.load(), n_nodes_deep.load(), n_path_taken.load(), n_extra_hits.load(),
           n_selections.load(), n_sel_entries.load(), n_sel_kept.load(),
           n_selections > 0 ? (double) n_sel_entries / n_selections : 0.0,
           n_selections > 0 ? (double) n_sel_kept / n_selections : 0.0,
           n_same_module_vetoed.load(), n_hole_slot_reserved.load(),
           n_best_short_offered.load(), n_best_short_taken.load(),
           n_kalman_calls_d0 > 0 ? (double) n_kalman_lanes_d0 / n_kalman_calls_d0 : 0.0, NN,
           n_kalman_calls_d0.load(),
           (n_kalman_calls - n_kalman_calls_d0) > 0 ?
             (double)(n_kalman_lanes - n_kalman_lanes_d0) / (n_kalman_calls - n_kalman_calls_d0) : 0.0,
           (long)(n_kalman_calls - n_kalman_calls_d0),
           n_arena_layers > 0 ? (double) n_arena_hw_sum / n_arena_layers : 0.0, n_arena_hw_max.load(),
           n_early_selections.load(), n_selections.load(),
           n_late_max_cands.load(), n_late_max_cands_dropped.load(), n_late_max_cands_long.load());
    printf("  sister hole: %ld holes dropped for a sibling on the sister sensor\n", n_sister_hole_dropped.load());
    print_dominance();
  }

  // Counts the dominated entries among the kept m_sel[0, n_keep), see
  // V2p2PolicyCounters::n_dom_*. Measurement only: nothing is changed.
  void MkFinderV2p2::count_dominated_kept(const CombCandidate &ccand, int n_keep) const {
    auto &C = g_v2p2_policy_counters;
    const int lay = m_rz_limits.layer_info_1().layer_id();
    if (lay < 0 || lay >= V2p2PolicyCounters::k_dom_layers)
      return;
    const bool full = (int) m_sel.size() > n_keep;
    ++C.n_dom_sel[lay];
    C.n_dom_kept[lay] += n_keep;
    if (full) {
      ++C.n_dom_sel_full[lay];
      C.n_dom_kept_full[lay] += n_keep;
    }
    auto detid_of = [this](int layer, int idx) {
      const LayerOfHits &L = mp_job->m_event_of_hits[layer];
      return L.layer_info().module_info(L.refHit(idx).detIDinLayer()).detid;
    };
    for (int k = 0; k < n_keep; ++k) {
      const SelEntry &e = m_sel[k];
      bool sub = false, hole = false, sister = false;
      if (e.node_idx >= 0) {
        for (int j = 0; j < n_keep && !sub; ++j) {
          if (j == k || m_sel[j].node_idx < 0 || m_sel[j].tcand_idx != e.tcand_idx)
            continue;
          for (int ci = m_node_arena[m_sel[j].node_idx].m_parent_idx; ci >= 0; ci = m_node_arena[ci].m_parent_idx)
            if (ci == e.node_idx) {
              sub = true;
              break;
            }
        }
      } else if (e.add_fake) {
        for (int j = 0; j < n_keep && !sister; ++j) {
          if (j == k || m_sel[j].node_idx < 0 || m_sel[j].tcand_idx != e.tcand_idx)
            continue;
          hole = true;
          const TrackCand &tc = ccand[e.tcand_idx];
          const int pl = tc.getLastHitLyr(), pi = tc.getLastHitIdx();
          if (pi < 0 || pl != lay - 1 || (lay & 1) == 0 || lay < 4 || (lay > 15 && lay < 28))
            continue;
          int root = m_sel[j].node_idx;
          while (m_node_arena[root].m_parent_idx >= 0)
            root = m_node_arena[root].m_parent_idx;
          const auto &hot = m_node_arena[root].m_hot;
          sister = detid_of(hot.layer, hot.index) == detid_of(pl, pi) + 1;
        }
      }
      if (sub) {
        ++C.n_dom_sub[lay];
        if (full) ++C.n_dom_sub_full[lay];
      }
      if (hole) {
        ++C.n_dom_hole[lay];
        if (full) ++C.n_dom_hole_full[lay];
      }
      if (sister) {
        ++C.n_dom_sister[lay];
        if (full) ++C.n_dom_sister_full[lay];
      }
    }
  }

  //----------------------------------------------------------------------------
  // surface_reference_dq() -- dq_track referenced to the layer surface (radial
  // normal in the barrel, z in the endcap), evaluated at m_sp2 where cov_ex
  // lives. Off by default (Diag::mkbins_surface_q); the per-hit version with the
  // module normal is MkFinderV2p2::surface_referenced_dq(). See
  // doc/MkFinderV2p2-DesignNotes.md, "Search window".
  //----------------------------------------------------------------------------

  void MkBins::surface_reference_dq(const MkBinTrackCovExtract &cov_ex) {
    // Clamp on the amplification. g = cot(theta) for a radial barrel track, so 20
    // is |eta| ~ 3.7, beyond the tracker.
    constexpr float kMaxSlope = 20.0f;

    for (int i = 0; i < m_n_proc; ++i) {
      const float x = m_sp2.x[i], y = m_sp2.y[i];
      const float r2 = x * x + y * y;
      if (r2 <= 0.0f)
        continue;
      const float rinv = 1.0f / std::sqrt(r2);
      const float nx = x * rinv, ny = y * rinv;   // radial unit vector

      const float pr = nx * m_sp2.px[i] + ny * m_sp2.py[i];  // p . n_radial
      const float pz = m_sp2.pz[i];

      const float c00 = cov_ex.m_cov_0_0[i], c01 = cov_ex.m_cov_0_1[i];
      const float c11 = cov_ex.m_cov_1_1[i], c22 = cov_ex.m_cov_2_2[i];
      const float c02 = cov_ex.m_cov_0_2[i], c12 = cov_ex.m_cov_1_2[i];

      float var;
      if (m_is_barrel) {
        // v = e_z - (p_z / (p.n)) * n ; note |p| cancels out of the ratio.
        if (pr == 0.0f)
          continue;
        float g = pz / pr;
        g = std::clamp(g, -kMaxSlope, kMaxSlope);
        const float v0 = -g * nx, v1 = -g * ny;
        var = v0 * v0 * c00 + v1 * v1 * c11 + c22 + 2.0f * (v0 * v1 * c01 + v0 * c02 + v1 * c12);
      } else {
        // w = r^ - ((r^.p) / p_z) * e_z
        if (pz == 0.0f)
          continue;
        float ginv = pr / pz;
        ginv = std::clamp(ginv, -kMaxSlope, kMaxSlope);
        var = nx * nx * c00 + ny * ny * c11 + ginv * ginv * c22 +
              2.0f * (nx * ny * c01 - ginv * nx * c02 - ginv * ny * c12);
      }

      if (var > 0.0f)
        m_dq_track[i] = 3.0f * std::sqrt(var);
    }
  }

  //==============================================================================
  // Diag::final_beam_purity: each seed's beam at the end of the v2p2 forward search,
  // after the final sort, against truth. The particle is the majority sim track of
  // the seed's hits; a candidate is "clean" when more than 3/4 of its found hits are
  // that particle's (the MTV rule). Per population: the pick (rank 0) is clean; it is
  // not but a clean candidate is in the beam (its rank, hits, OT hits); or no clean
  // candidate is left. Accumulated across events, single-threaded.
  //==============================================================================
  namespace {
    struct FinalBeamDiag {
      long n = 0, pick_clean = 0, other_clean = 0, none_clean = 0;
      long sum_rank = 0, sum_pick_hits = 0, sum_clean_hits = 0, sum_pick_ot = 0, sum_clean_ot = 0;
      long clean_more_ot = 0, beam_size = 0;
    };
    FinalBeamDiag g_fbd[3];  // 0 all, 1 seed 1.7 <= |eta| < 2.7 pT < 0.9, 2 the same pT >= 0.9

    // Truth-free features of a candidate, for "could anything at the final pick tell the clean
    // candidate from the dirty one". Per-hit quantities over the non-seed found hits.
    struct FbdFeat {
      float nfound, chi2, chi2_per_hit, max_chi2, ot_chi2, nholes, score, n_ot;
    };
    constexpr int FBD_NF = 8;
    const char *fbd_feat_name[FBD_NF] = {"found hits", "chi2", "chi2 / hit", "max hit chi2",
                                         "mean OT hit chi2", "inside holes", "track_score_func", "OT hits"};
    const bool fbd_feat_higher_better[FBD_NF] = {true, false, false, false, false, false, true, true};
    float fbd_feat_get(const FbdFeat &f, int k) {
      const float v[FBD_NF] = {f.nfound, f.chi2, f.chi2_per_hit, f.max_chi2, f.ot_chi2, f.nholes, f.score, f.n_ot};
      return v[k];
    }
    // Pairwise: [population][A = dirty pick vs the first clean, B = clean pick vs the first dirty]
    // per feature, how often it prefers the alternative (ties count half).
    struct FbdPair {
      long n = 0;
      double prefer_alt[FBD_NF] = {};
    };
    FbdPair g_fbp[3][2];

    // Expected background per (seed origin index, layer) of the forward search: the hits in the
    // density region, averaged over the seed's candidates that searched the layer.
    std::map<std::pair<int, int>, std::pair<double, int>> g_fbd_bg;
    // Stop-rule emulation: stop a pick at the first layer whose expected background exceeds X.
    // [population][X][0 dirty picks made clean, 1 dirty picks, 2 clean picks made dirty,
    //  3 clean picks stopped, 4 clean picks, 5 found hits lost on clean picks, 6 OT hits lost]
    constexpr int FBD_NX = 7;
    const float fbd_x[FBD_NX] = {0.5f, 1.f, 2.f, 3.f, 5.f, 8.f, 12.f};
    long g_fbs[3][FBD_NX][7];
    // [population][background bin] first wrong hit of a dirty pick, and every OT hit of clean picks
    constexpr int FBD_NB = 7;
    const float fbd_bin[FBD_NB] = {0.5f, 1.f, 2.f, 3.f, 5.f, 8.f, 1e9f};
    long g_fbh_wrong[3][FBD_NB], g_fbh_clean[3][FBD_NB];
    // Seed purity: [population][0 seed all L, 1 seed with one wrong or unlinked hit, 2 worse]
    // x [pick clean, pick dirty]; and where a dirty pick's first wrong hit sits: [0 seed, 1 pixel
    // after the seed, 2 OT]
    long g_fbseed[3][3][2], g_fbfirst[3][3];
    int fbd_bin_of(float v) { int k = 0; while (v > fbd_bin[k]) ++k; return k; }
    // Margins: [population][A/B][margin feature] the pick-minus-alternative difference, for a
    // threshold scan (switch only when the pick is worse by more than t).
    constexpr int FBD_NM = 4;
    const char *fbd_margin_name[FBD_NM] = {"chi2", "chi2 / hit", "mean OT hit chi2", "chi2 per extra hit"};
    std::vector<float> g_fbm[3][2][FBD_NM];

    bool fbd_is_ot(int l) { return (l >= 4 && l < 16) || (l >= 28 && l < 38) || l >= 50; }
  }  // namespace

  void v2p2_final_beam_bg_record(int seed, int layer, int layer_sec, float n_bg) {
    for (int l : {layer, layer_sec}) {
      if (l < 0)
        continue;
      auto &e = g_fbd_bg[{seed, l}];
      e.first += n_bg;
      ++e.second;
    }
  }

  //==============================================================================
  // Final-pick scorers. v2p2_final_pick_record() keeps each seed's candidates as they
  // stand before mergeCandsAndBestShortOne(), the best short one included, while
  // score() still holds the summed v2p2 layer-step likelihood. v2p2_final_beam_diag()
  // then asks, per seed, which candidate each scorer of a family would pick and
  // whether it is clean (the MTV 3/4 rule over its found hits):
  //   0  phase1:default, as the final sort uses it
  //   1  the likelihood as accumulated in the search
  //   2+ the likelihood with every non-seed hit re-weighted to one hit efficiency
  //      eps, minus h per inside hole and t per tail hole.
  // The merge itself can drop a candidate, so the set is taken before it.
  //==============================================================================
  namespace {
    std::mutex g_fpr_mutex;
    std::map<int, std::vector<TrackCand>> g_fpr;  // seed_origin_index -> candidates
    // eps < 0: the search's own per-group efficiencies, no re-weighting.
    constexpr int FPR_NE = 6, FPR_NH = 5, FPR_NT = 4;
    const float fpr_eps[FPR_NE] = {-1.f, 0.03f, 0.05f, 0.10f, 0.20f, 0.30f};
    const float fpr_hole[FPR_NH] = {0.f, 4.f, 8.f, 12.f, 16.f};
    const float fpr_tail[FPR_NT] = {0.f, 3.f, 6.f, 9.f};
    constexpr int FPR_NS = 2 + FPR_NE * FPR_NH * FPR_NT;
    // [population][scorer]; population as g_fbd: 0 all, 1 / 2 the forward low / high pT.
    long g_fpr_n[3] = {}, g_fpr_any[3] = {}, g_fpr_p1_same[3] = {};
    long g_fpr_clean[3][FPR_NS] = {}, g_fpr_gain[3][FPR_NS] = {}, g_fpr_loss[3][FPR_NS] = {};

    int fpr_group(int layer) {
      const LayerInfo &li = Config::TrkInfo[layer];
      if (li.is_pixel())
        return li.is_barrel() ? 0 : 1;
      if (li.is_barrel())
        return layer < 10 ? 2 : 3;
      return 4;
    }
    float fpr_logit(float e) { return std::log(e / (1.0f - e)); }
    void fpr_grid(int s, float &eps, float &h, float &t) {
      const int g = s - 2;
      eps = fpr_eps[g / (FPR_NH * FPR_NT)];
      h = fpr_hole[(g / FPR_NT) % FPR_NH];
      t = fpr_tail[g % FPR_NT];
    }
  }  // namespace

  void v2p2_final_pick_record(const CombCandidate &cc) {
    std::lock_guard<std::mutex> lk(g_fpr_mutex);
    auto &v = g_fpr[cc.seed_origin_index()];
    v.clear();
    for (int j = 0; j < (int)cc.size(); ++j)
      v.push_back(cc[j]);
    if (cc.refBestShortCand().combCandidate())
      v.push_back(cc.refBestShortCand());
  }

  void v2p2_final_beam_diag(const Event *ev, const EventOfCombCandidates &eoccs) {
    auto label_of = [&](int lyr, int idx) {
      if (idx < 0) return -1;
      const int id = ev->layerHits_[lyr][idx].mcHitID();
      return id >= 0 ? ev->simHitsInfo_[id].mcTrackID() : -1;
    };
    for (int i = 0; i < eoccs.size(); ++i) {
      const CombCandidate &cc = eoccs[i];
      if (cc.empty()) continue;
      const Track &seed = ev->currentSeed(cc.seed_origin_index());
      std::map<int, int> cnt;
      for (int h = 0; h < seed.nTotalHits(); ++h) {
        const int l = label_of(seed.getHitLyr(h), seed.getHitIdx(h));
        if (l >= 0) ++cnt[l];
      }
      if (cnt.empty()) continue;
      int L = -1, nL = 0;
      for (auto &[l, c] : cnt) if (c > nL) { L = l; nL = c; }
      const float aeta = std::abs(seed.momEta()), pt = seed.pT();
      int pop = -1;
      if (aeta >= 1.7f && aeta < 2.7f) pop = pt < 0.9f ? 1 : 2;

      auto features = [&](const TrackCand &tc) {
        FbdFeat f{};
        f.nfound = tc.nFoundHits();
        f.chi2 = tc.chi2();
        f.nholes = tc.nInsideMinusOneHits();
        f.score = tc.score();
        std::vector<const HoTNode *> nodes;
        int nh = tc.nTotalHits(), ch = tc.lastCcIndex();
        while (--nh >= 0 && ch >= 0) {
          const HoTNode &hn = cc.hot_node(ch);
          nodes.push_back(&hn);
          ch = hn.m_prev_idx;
        }
        // nodes run from the last hit back to the first; the seed's hits are the last ones.
        const int n_seed = tc.getNSeedHits();
        int nfs = 0, n_ot = 0;
        float sum = 0, mx = 0, sum_ot = 0;
        for (int j = 0; j < (int)nodes.size() - n_seed; ++j) {
          const HoTNode &hn = *nodes[j];
          if (hn.m_hot.index < 0) continue;
          ++nfs; sum += hn.m_chi2; mx = std::max(mx, hn.m_chi2);
          if (fbd_is_ot(hn.m_hot.layer)) { ++n_ot; sum_ot += hn.m_chi2; }
        }
        f.chi2_per_hit = nfs > 0 ? sum / nfs : 0.f;
        f.max_chi2 = mx;
        f.ot_chi2 = n_ot > 0 ? sum_ot / n_ot : 0.f;
        f.n_ot = n_ot;
        return f;
      };
      auto pair_fill = [&](int pop_k, int ab, const FbdFeat &pick, const FbdFeat &alt) {
        FbdPair &p = g_fbp[pop_k][ab];
        ++p.n;
        const float dh = pick.nfound - alt.nfound;
        const float m[FBD_NM] = {pick.chi2 - alt.chi2, pick.chi2_per_hit - alt.chi2_per_hit,
                                 pick.ot_chi2 - alt.ot_chi2,
                                 dh > 0 ? (pick.chi2 - alt.chi2) / dh : pick.chi2 - alt.chi2};
        for (int k = 0; k < FBD_NM; ++k) g_fbm[pop_k][ab][k].push_back(m[k]);
        for (int k = 0; k < FBD_NF; ++k) {
          const float a = fbd_feat_get(alt, k), b = fbd_feat_get(pick, k);
          if (a == b) p.prefer_alt[k] += 0.5;
          else if ((a > b) == fbd_feat_higher_better[k]) p.prefer_alt[k] += 1.0;
        }
      };

      int first_clean = -1, first_dirty = -1, pick_hits = 0, pick_ot = 0, clean_hits = 0, clean_ot = 0;
      for (int r = 0; r < (int)cc.size(); ++r) {
        const Track t = cc[r].exportTrack(true);
        int nf = 0, nl = 0, not_l = 0;
        for (int h = 0; h < t.nTotalHits(); ++h) {
          const int idx = t.getHitIdx(h), lyr = t.getHitLyr(h);
          if (idx < 0) continue;
          ++nf;
          if (label_of(lyr, idx) == L) { ++nl; if (fbd_is_ot(lyr)) ++not_l; }
        }
        const bool clean = 4 * nl > 3 * nf;
        if (r == 0) { pick_hits = nf; pick_ot = not_l; }
        if (clean && first_clean < 0) { first_clean = r; clean_hits = nf; clean_ot = not_l; }
        if (!clean && first_dirty < 0) first_dirty = r;
      }
      // Stop-rule emulation on the pick, and the background where hits were taken.
      {
        const TrackCand &pk = cc[0];
        const int seed_idx = cc.seed_origin_index();
        std::vector<const HoTNode *> nodes;
        int nh = pk.nTotalHits(), ch = pk.lastCcIndex();
        while (--nh >= 0 && ch >= 0) {
          const HoTNode &hn = cc.hot_node(ch);
          nodes.push_back(&hn);
          ch = hn.m_prev_idx;
        }
        std::reverse(nodes.begin(), nodes.end());  // first hit first
        const int n_seed = pk.getNSeedHits();
        auto bg_at = [&](int layer) {
          auto it = g_fbd_bg.find({seed_idx, layer});
          return it == g_fbd_bg.end() ? -1.f : float(it->second.first / it->second.second);
        };
        const bool pick_clean = first_clean == 0;
        {
          int ns = 0, nsl = 0;
          for (int j = 0; j < n_seed && j < (int)nodes.size(); ++j) {
            const HoTNode &hn = *nodes[j];
            if (hn.m_hot.index < 0) continue;
            ++ns;
            if (label_of(hn.m_hot.layer, hn.m_hot.index) == L) ++nsl;
          }
          const int sc = nsl == ns ? 0 : (nsl == ns - 1 ? 1 : 2);
          int fw = -1;
          for (int j = 0; j < (int)nodes.size() && fw < 0; ++j) {
            const HoTNode &hn = *nodes[j];
            if (hn.m_hot.index < 0) continue;
            if (label_of(hn.m_hot.layer, hn.m_hot.index) != L)
              fw = j < n_seed ? 0 : (fbd_is_ot(hn.m_hot.layer) ? 2 : 1);
          }
          for (int k : {0, pop}) {
            if (k < 0) continue;
            ++g_fbseed[k][sc][pick_clean ? 0 : 1];
            if (!pick_clean && fw >= 0) ++g_fbfirst[k][fw];
          }
        }
        bool seen_wrong = false;
        for (int j = n_seed; j < (int)nodes.size(); ++j) {
          const HoTNode &hn = *nodes[j];
          if (hn.m_hot.index < 0) continue;
          const float bg = bg_at(hn.m_hot.layer);
          if (bg < 0) continue;
          const bool right = label_of(hn.m_hot.layer, hn.m_hot.index) == L;
          for (int k : {0, pop}) {
            if (k < 0) continue;
            if (!pick_clean && !right && !seen_wrong) ++g_fbh_wrong[k][fbd_bin_of(bg)];
            if (pick_clean && fbd_is_ot(hn.m_hot.layer)) ++g_fbh_clean[k][fbd_bin_of(bg)];
          }
          if (!right) seen_wrong = true;
        }
        for (int xi = 0; xi < FBD_NX; ++xi) {
          int stop = (int)nodes.size();
          for (int j = n_seed; j < (int)nodes.size(); ++j) {
            const float bg = bg_at(nodes[j]->m_hot.layer);
            if (bg > fbd_x[xi]) { stop = j; break; }
          }
          int nf = 0, nl = 0, lost = 0, lost_ot = 0;
          for (int j = 0; j < (int)nodes.size(); ++j) {
            const HoTNode &hn = *nodes[j];
            if (hn.m_hot.index < 0) continue;
            if (j < stop) { ++nf; if (label_of(hn.m_hot.layer, hn.m_hot.index) == L) ++nl; }
            else { ++lost; if (fbd_is_ot(hn.m_hot.layer)) ++lost_ot; }
          }
          const bool clean_after = 4 * nl > 3 * nf;
          for (int k : {0, pop}) {
            if (k < 0) continue;
            long *c = g_fbs[k][xi];
            if (!pick_clean) { ++c[1]; if (clean_after) ++c[0]; }
            else { ++c[4]; if (!clean_after) ++c[2]; if (stop < (int)nodes.size()) ++c[3]; c[5] += lost; c[6] += lost_ot; }
          }
        }
      }
      if (first_clean > 0 || (first_clean == 0 && first_dirty > 0)) {
        const int alt = first_clean > 0 ? first_clean : first_dirty;
        const FbdFeat fp = features(cc[0]), fa = features(cc[alt]);
        for (int k : {0, pop})
          if (k >= 0) pair_fill(k, first_clean > 0 ? 0 : 1, fp, fa);
      }
      // Final-pick scorers on the candidates as they stood before the merge.
      if (auto it = g_fpr.find(cc.seed_origin_index()); it != g_fpr.end() && !it->second.empty()) {
        static const track_score_func p1 = IterationConfig::get_track_scorer("phase1:default");
        const auto &cands = it->second;
        const int nc = cands.size();
        std::vector<char> cl(nc);
        std::vector<float> s_p1(nc), s_l(nc);
        std::vector<int> nin(nc), ntail(nc);
        std::vector<std::array<int, 5>> ng(nc);
        bool any = false;
        for (int j = 0; j < nc; ++j) {
          const TrackCand &tc = cands[j];
          const Track t = tc.exportTrack(true);
          int nf = 0, nl = 0;
          for (int h = 0; h < t.nTotalHits(); ++h) {
            const int idx = t.getHitIdx(h), lyr = t.getHitLyr(h);
            if (idx < 0) continue;
            ++nf;
            if (label_of(lyr, idx) == L) ++nl;
          }
          cl[j] = 4 * nl > 3 * nf;
          any |= cl[j];
          s_p1[j] = getScoreCand(p1, tc);
          s_l[j] = tc.score();
          nin[j] = tc.nInsideMinusOneHits();
          ntail[j] = tc.nTailMinusOneHits();
          ng[j] = {0, 0, 0, 0, 0};
          int nh = tc.nTotalHits(), ch = tc.lastCcIndex();
          int n_walk = nh - tc.getNSeedHits();  // the chain runs from the last hit back; seed hits last
          while (--nh >= 0 && ch >= 0 && n_walk-- > 0) {
            const HoTNode &hn = tc.combCandidate()->hot_node(ch);
            if (hn.m_hot.index >= 0)
              ++ng[j][fpr_group(hn.m_hot.layer)];
            ch = hn.m_prev_idx;
          }
        }
        const auto &fp = Config::ItrInfo[0].m_v2p2.score_outward;
        float le[5];
        for (int g = 0; g < 5; ++g)
          le[g] = fpr_logit(fp.hit_eff_grp[g] >= 0.f ? fp.hit_eff_grp[g] : fp.hit_eff);
        auto score_of = [&](int s, int j) -> float {
          if (s == 0) return s_p1[j];
          if (s == 1) return s_l[j];
          float eps, h, t;
          fpr_grid(s, eps, h, t);
          const float l = eps > 0.f ? fpr_logit(eps) : 0.f;
          float v = s_l[j] - h * nin[j] - t * ntail[j];
          if (eps > 0.f)
            for (int g = 0; g < 5; ++g) v += ng[j][g] * (l - le[g]);
          return v;
        };
        auto pick_of = [&](int s) {
          int b = 0;
          for (int j = 1; j < nc; ++j)
            if (score_of(s, j) > score_of(s, b)) b = j;
          return b;
        };
        const int p1_pick = pick_of(0);
        const bool p1_same = cands[p1_pick].lastCcIndex() == cc[0].lastCcIndex();
        for (int k : {0, pop}) {
          if (k < 0) continue;
          ++g_fpr_n[k];
          if (any) ++g_fpr_any[k];
          if (p1_same) ++g_fpr_p1_same[k];
        }
        for (int s = 0; s < FPR_NS; ++s) {
          const bool c = cl[pick_of(s)], c1 = cl[p1_pick];
          for (int k : {0, pop}) {
            if (k < 0) continue;
            if (c) ++g_fpr_clean[k][s];
            if (c && !c1) ++g_fpr_gain[k][s];
            if (!c && c1) ++g_fpr_loss[k][s];
          }
        }
      }

      for (int k : {0, pop}) {
        if (k < 0) continue;
        FinalBeamDiag &d = g_fbd[k];
        ++d.n; d.beam_size += cc.size();
        if (first_clean == 0) ++d.pick_clean;
        else if (first_clean > 0) {
          ++d.other_clean; d.sum_rank += first_clean;
          d.sum_pick_hits += pick_hits; d.sum_clean_hits += clean_hits;
          d.sum_pick_ot += pick_ot; d.sum_clean_ot += clean_ot;
          if (clean_ot > pick_ot) ++d.clean_more_ot;
        } else ++d.none_clean;
      }
    }
    g_fbd_bg.clear();
    g_fpr.clear();
  }

  void v2p2_final_beam_diag_reset() {
    for (auto &d : g_fbd) d = FinalBeamDiag();
    for (auto &pp : g_fbp) for (auto &p : pp) p = FbdPair();
    for (auto &a : g_fbm) for (auto &b : a) for (auto &v : b) v.clear();
    g_fbd_bg.clear();
    std::memset(g_fbs, 0, sizeof(g_fbs));
    std::memset(g_fbh_wrong, 0, sizeof(g_fbh_wrong));
    std::memset(g_fbh_clean, 0, sizeof(g_fbh_clean));
    std::memset(g_fbseed, 0, sizeof(g_fbseed));
    std::memset(g_fbfirst, 0, sizeof(g_fbfirst));
    g_fpr.clear();
    std::memset(g_fpr_n, 0, sizeof(g_fpr_n));
    std::memset(g_fpr_any, 0, sizeof(g_fpr_any));
    std::memset(g_fpr_p1_same, 0, sizeof(g_fpr_p1_same));
    std::memset(g_fpr_clean, 0, sizeof(g_fpr_clean));
    std::memset(g_fpr_gain, 0, sizeof(g_fpr_gain));
    std::memset(g_fpr_loss, 0, sizeof(g_fpr_loss));
  }

  void v2p2_final_beam_diag_report() {
    const char *pn[3] = {"all seeds", "seed 1.7-2.7 pT < 0.9", "seed 1.7-2.7 pT >= 0.9"};
    printf("\nFINAL BEAM PURITY, v2p2 forward search, after the final sort (clean = > 3/4 of the found\n"
           "hits from the seed's majority particle)\n");
    printf("%-24s %8s %6s %10s %10s %10s | %6s %7s %7s %7s %7s %10s\n", "population", "seeds", "beam",
           "pick clean", "other cln", "none clean", "rank", "hits pk", "hits cl", "OT pk", "OT cl", "cl more OT");
    for (int k = 0; k < 3; ++k) {
      const FinalBeamDiag &d = g_fbd[k];
      const double n = std::max(1L, d.n), o = std::max(1L, d.other_clean);
      printf("%-24s %8ld %6.2f %9.1f%% %9.1f%% %9.1f%% | %6.2f %7.2f %7.2f %7.2f %7.2f %9.1f%%\n", pn[k], d.n,
             d.beam_size / n, 100. * d.pick_clean / n, 100. * d.other_clean / n, 100. * d.none_clean / n,
             d.sum_rank / o, d.sum_pick_hits / o, d.sum_clean_hits / o, d.sum_pick_ot / o, d.sum_clean_ot / o,
             100. * d.clean_more_ot / o);
    }
    printf("\nCOULD A FEATURE TELL THEM APART? Pairwise, pick against one alternative in the beam.\n"
           "  A: pick dirty, alternative = the best-ranked clean candidate -- want 'prefers alt' HIGH.\n"
           "  B: pick clean, alternative = the best-ranked dirty candidate -- want it LOW.\n"
           "  net = n_A * A - n_B * B: seeds a pairwise switch on this feature alone would fix minus break.\n");
    for (int k = 0; k < 3; ++k) {
      const FbdPair &A = g_fbp[k][0], &B = g_fbp[k][1];
      printf("  %s: n_A %ld, n_B %ld\n", pn[k], A.n, B.n);
      printf("    %-18s %9s %9s %10s\n", "feature", "A alt", "B alt", "net");
      for (int f = 0; f < FBD_NF; ++f) {
        const double a = A.n ? A.prefer_alt[f] / A.n : 0, b = B.n ? B.prefer_alt[f] / B.n : 0;
        printf("    %-18s %8.1f%% %8.1f%% %10.0f\n", fbd_feat_name[f], 100. * a, 100. * b,
               A.prefer_alt[f] - B.prefer_alt[f]);
      }
    }
    printf("\nEXPECTED BACKGROUND (hits in the density region of the layer step, seed average)\n"
           "  where a DIRTY pick took its first wrong hit, and where CLEAN picks took their OT hits\n");
    for (int k = 0; k < 3; ++k) {
      long tw = 0, tc = 0;
      for (int b = 0; b < FBD_NB; ++b) { tw += g_fbh_wrong[k][b]; tc += g_fbh_clean[k][b]; }
      printf("  %-24s %-14s", pn[k], "bg <=");
      for (int b = 0; b < FBD_NB; ++b) { if (b + 1 < FBD_NB) printf(" %6g", fbd_bin[b]); else printf("   more"); }
      printf("\n  %-24s %-14s", "", "first wrong");
      for (int b = 0; b < FBD_NB; ++b) printf(" %5.1f%%", 100. * g_fbh_wrong[k][b] / std::max(1L, tw));
      printf("   (%ld)\n  %-24s %-14s", tw, "", "clean OT hits");
      for (int b = 0; b < FBD_NB; ++b) printf(" %5.1f%%", 100. * g_fbh_clean[k][b] / std::max(1L, tc));
      printf("   (%ld)\n", tc);
    }
    printf("\nSEED PURITY (seed hits of the seed's majority particle) against the pick, and where a\n"
           "dirty pick's FIRST wrong hit is\n");
    for (int k = 0; k < 3; ++k) {
      const long (*q)[2] = g_fbseed[k];
      const long d = q[0][1] + q[1][1] + q[2][1], c = q[0][0] + q[1][0] + q[2][0];
      const long fw = g_fbfirst[k][0] + g_fbfirst[k][1] + g_fbfirst[k][2];
      printf("  %-24s clean picks %7ld: seed pure %5.1f%%, one off %5.1f%%, worse %5.1f%% | dirty picks %6ld: "
             "seed pure %5.1f%%, one off %5.1f%%, worse %5.1f%% | first wrong in seed %5.1f%%, pixel %5.1f%%, OT %5.1f%%\n",
             pn[k], c, 100. * q[0][0] / std::max(1L, c), 100. * q[1][0] / std::max(1L, c), 100. * q[2][0] / std::max(1L, c),
             d, 100. * q[0][1] / std::max(1L, d), 100. * q[1][1] / std::max(1L, d), 100. * q[2][1] / std::max(1L, d),
             100. * g_fbfirst[k][0] / std::max(1L, fw), 100. * g_fbfirst[k][1] / std::max(1L, fw),
             100. * g_fbfirst[k][2] / std::max(1L, fw));
    }
    printf("\nSTOP-RULE EMULATION on the pick: stop at the first layer after the seed whose expected\n"
           "background exceeds X. dirty->clean of dirty picks; clean->dirty and stopped of clean\n"
           "picks; found and OT hits lost per clean pick.\n");
    for (int k = 0; k < 3; ++k) {
      printf("  %s\n", pn[k]);
      for (int xi = 0; xi < FBD_NX; ++xi) {
        const long *c = g_fbs[k][xi];
        const double nc = std::max(1L, c[4]);
        printf("    X %4g: dirty->clean %5.1f%% (%5ld of %5ld)  clean->dirty %4.1f%%  clean stopped %5.1f%%"
               "  lost/clean: hits %5.2f OT %5.2f\n",
               fbd_x[xi], 100. * c[0] / std::max(1L, c[1]), c[0], c[1], 100. * c[2] / nc, 100. * c[3] / nc,
               c[5] / nc, c[6] / nc);
      }
    }
    printf("\nMARGIN SCAN: switch to the alternative only when pick - alt > t. Fraction above t in A (fix)\n"
           "and B (break), and net = n_A * fA - n_B * fB.\n");
    const float ts[] = {0.f, 2.f, 5.f, 10.f, 20.f, 40.f};
    for (int k = 0; k < 3; ++k) {
      printf("  %s\n", pn[k]);
      for (int mi = 0; mi < FBD_NM; ++mi) {
        const auto &va = g_fbm[k][0][mi], &vb = g_fbm[k][1][mi];
        printf("    %-20s", fbd_margin_name[mi]);
        for (float t : ts) {
          long na = 0, nb = 0;
          for (float v : va) na += v > t;
          for (float v : vb) nb += v > t;
          printf("  t=%-3g %5.1f%%/%4.1f%% %+6ld", t, 100. * na / std::max<size_t>(1, va.size()),
                 100. * nb / std::max<size_t>(1, vb.size()), na - nb);
        }
        printf("\n");
      }
    }

    printf("\nFINAL-PICK SCORERS on each seed's candidates before the final sort: clean picks (MTV 3/4),\n"
           "and against phase1:default the seeds a scorer turns clean (gain) or dirty (loss). Grid:\n"
           "the likelihood with every non-seed hit at one hit efficiency eps, minus h per inside hole\n"
           "and t per tail hole. The top grid scorers by net.\n");
    const char *fpn[3] = {"all seeds", "seed 1.7-2.7 pT < 0.9", "seed 1.7-2.7 pT >= 0.9"};
    for (int k = 0; k < 3; ++k) {
      const double n = std::max(1L, g_fpr_n[k]);
      printf("  %-24s seeds %7ld, some candidate clean %5.1f%%, phase1 re-pick = actual pick %5.1f%%\n", fpn[k],
             g_fpr_n[k], 100. * g_fpr_any[k] / n, 100. * g_fpr_p1_same[k] / n);
      auto line = [&](int s, const char *name) {
        printf("    %-34s clean %6.2f%%  gain %6ld  loss %6ld  net %+6ld\n", name, 100. * g_fpr_clean[k][s] / n,
               g_fpr_gain[k][s], g_fpr_loss[k][s], g_fpr_gain[k][s] - g_fpr_loss[k][s]);
      };
      line(0, "phase1:default");
      line(1, "likelihood as in the search");
      line(2 + 2 * FPR_NT + 1, "llh eps search, h 8, t 3 (always)");
      std::vector<int> ord;
      for (int s = 2; s < FPR_NS; ++s) ord.push_back(s);
      std::sort(ord.begin(), ord.end(), [&](int a, int b) {
        return g_fpr_gain[k][a] - g_fpr_loss[k][a] > g_fpr_gain[k][b] - g_fpr_loss[k][b];
      });
      for (int i = 0; i < 10 && i < (int)ord.size(); ++i) {
        float eps, h, t;
        fpr_grid(ord[i], eps, h, t);
        char nm[64];
        if (eps > 0.f)
          snprintf(nm, sizeof(nm), "llh eps %.2f, h %g, t %g", eps, h, t);
        else
          snprintf(nm, sizeof(nm), "llh eps search, h %g, t %g", h, t);
        line(ord[i], nm);
      }
    }
  }

  //==============================================================================
  // Diag::dupclean_pairs: the pixel-seed duplicate cleaner's decisions against truth
  //==============================================================================

  // Each pair the cleaner marks: the track it kept and the one it dropped. A track is
  // clean when more than 3/4 of its found hits are one particle's (the MTV rule, all
  // hits counted). When exactly one of the two is clean, the rules below say which one
  // they would have kept, from what the cleaner can see on a Track.
  namespace {
    struct DcTrk {
      std::vector<std::pair<int, int>> hits;  // (layer, index) of found hits
      int nfound, nholes, ntail;
      float chi2, pt, eta, score;
    };
    struct DcPair {
      DcTrk kept, dropped;
      int rule;
    };
    std::vector<DcPair> g_dc_pending;

    constexpr int DC_NP = 3;  // populations by the kept track: all, 1.7-2.7 pT < 0.9, 1.7-2.7 pT >= 0.9
    constexpr int DC_NC = 5;  // same particle, two particles, kept clean, kept dirty, both dirty
    constexpr int DC_NR = 6;  // alternative rules, see dc_rule_names
    const char *dc_class_names[DC_NC] = {
        "both clean, same particle", "both clean, two particles", "kept clean, dropped dirty",
        "kept dirty, dropped clean", "both dirty"};
    const char *dc_rule_names[DC_NR] = {"score (the cleaner)",  "more found hits",
                                        "fewer inside holes",   "lower chi2",
                                        "lower chi2 / hit",     "more hits, then chi2"};
    long g_dc_n[DC_NP][2][DC_NC] = {};      // [pop][rule dR / shared][class]
    double g_dc_keep_clean[DC_NP][DC_NR] = {};  // on mixed pairs, ties count 0.5
    long g_dc_mixed[DC_NP] = {};
    long g_dc_events = 0;

    // +1: the rule keeps a, -1: b, 0: tie.
    int dc_rule(int r, const DcTrk &a, const DcTrk &b) {
      auto cmp = [](float x, float y) { return x > y ? 1 : (x < y ? -1 : 0); };
      const float ca = a.nfound > 0 ? a.chi2 / a.nfound : 0.f, cb = b.nfound > 0 ? b.chi2 / b.nfound : 0.f;
      switch (r) {
        case 0: return cmp(a.score, b.score);
        case 1: return cmp(a.nfound, b.nfound);
        case 2: return cmp(b.nholes, a.nholes);
        case 3: return cmp(b.chi2, a.chi2);
        case 4: return cmp(cb, ca);
        default: {
          const int c = cmp(a.nfound, b.nfound);
          return c != 0 ? c : cmp(b.chi2, a.chi2);
        }
      }
    }

    DcTrk dc_trk(const Track &t) {
      DcTrk d;
      for (int h = 0; h < t.nTotalHits(); ++h)
        if (t.getHitIdx(h) >= 0) d.hits.emplace_back(t.getHitLyr(h), t.getHitIdx(h));
      d.nfound = t.nFoundHits();
      d.nholes = t.nInsideMinusOneHits();
      d.ntail = t.nTailMinusOneHits();
      d.chi2 = t.chi2();
      d.pt = t.pT();
      d.eta = t.momEta();
      d.score = t.score();
      return d;
    }
  }  // namespace

  void v2p2_dupclean_record(const Track &kept, const Track &dropped, int rule) {
    g_dc_pending.push_back({dc_trk(kept), dc_trk(dropped), rule});
  }

  void v2p2_dupclean_diag(const Event *ev) {
    auto label_of = [&](int lyr, int idx) {
      const int id = ev->layerHits_[lyr][idx].mcHitID();
      return id >= 0 ? ev->simHitsInfo_[id].mcTrackID() : -1;
    };
    // Majority particle and whether it holds more than 3/4 of the found hits.
    auto truth = [&](const DcTrk &t, int &L) {
      std::map<int, int> cnt;
      for (auto &[l, i] : t.hits) {
        const int lab = label_of(l, i);
        if (lab >= 0) ++cnt[lab];
      }
      int nL = 0;
      L = -1;
      for (auto &[l, c] : cnt)
        if (c > nL) { L = l; nL = c; }
      return 4 * nL > 3 * (int)t.hits.size();
    };
    ++g_dc_events;
    for (const DcPair &p : g_dc_pending) {
      int lk, ld;
      const bool ck = truth(p.kept, lk), cd = truth(p.dropped, ld);
      const int cls = ck && cd ? (lk == ld ? 0 : 1) : ck ? 2 : cd ? 3 : 4;
      const float aeta = std::abs(p.kept.eta);
      const int pops[2] = {0, aeta >= 1.7f && aeta < 2.7f ? (p.kept.pt < 0.9f ? 1 : 2) : -1};
      for (int k : pops) {
        if (k < 0) continue;
        ++g_dc_n[k][p.rule][cls];
        if (cls == 2 || cls == 3) {
          ++g_dc_mixed[k];
          const DcTrk &cl = cls == 2 ? p.kept : p.dropped, &di = cls == 2 ? p.dropped : p.kept;
          for (int r = 0; r < DC_NR; ++r) {
            const int c = dc_rule(r, cl, di);
            g_dc_keep_clean[k][r] += c > 0 ? 1.0 : (c == 0 ? 0.5 : 0.0);
          }
        }
      }
    }
    g_dc_pending.clear();
  }

  void v2p2_dupclean_diag_reset() {
    g_dc_pending.clear();
    for (auto &a : g_dc_n) for (auto &b : a) for (auto &c : b) c = 0;
    for (auto &a : g_dc_keep_clean) for (auto &b : a) b = 0;
    for (auto &a : g_dc_mixed) a = 0;
    g_dc_events = 0;
  }

  void v2p2_dupclean_diag_report() {
    const char *pn[DC_NP] = {"all pairs", "kept 1.7-2.7 pT < 0.9", "kept 1.7-2.7 pT >= 0.9"};
    printf("\nDUPLICATE CLEANER DECISIONS, clean_duplicates_sharedhits_pixelseed, %ld events (clean = > 3/4\n"
           "of the found hits from one particle). Each marked pair, by the dR rule and the shared-hit rule;\n"
           "on mixed pairs (one clean), how often each rule keeps the clean one (ties 0.5).\n",
           g_dc_events);
    for (int k = 0; k < DC_NP; ++k) {
      long tot[2] = {0, 0};
      for (int r = 0; r < 2; ++r) for (int c = 0; c < DC_NC; ++c) tot[r] += g_dc_n[k][r][c];
      printf("  %-24s  pairs dR %ld, shared %ld\n", pn[k], tot[0], tot[1]);
      for (int c = 0; c < DC_NC; ++c)
        printf("    %-28s dR %7ld  shared %7ld\n", dc_class_names[c], g_dc_n[k][0][c], g_dc_n[k][1][c]);
      printf("    mixed pairs %ld, kept the clean one:\n", g_dc_mixed[k]);
      for (int r = 0; r < DC_NR; ++r)
        printf("      %-24s %6.2f %%\n", dc_rule_names[r],
               g_dc_mixed[k] > 0 ? 100.0 * g_dc_keep_clean[k][r] / g_dc_mixed[k] : 0.0);
    }
  }

}  // namespace mkfit
