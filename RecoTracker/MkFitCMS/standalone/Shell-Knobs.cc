// Shell knobs and tap switches, the section "Knobs" of Shell.h.
//
// A knob changes what the algorithm does in the next event processed. The
// per-iteration ones write every IterationConfig of Config::ItrInfo, the next
// run reads them through MkJob; the others write the standalone knobs of
// Config::V2p2 (src/V2p2Config.h) or the PropagationConfig.

#include "RecoTracker/MkFitCMS/standalone/Shell.h"

#include "RecoTracker/MkFitCore/interface/Config.h"
#include "RecoTracker/MkFitCore/interface/IterationConfig.h"
#include "RecoTracker/MkFitCore/interface/PropagationConfig.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/MkFitCore/interface/V2p2Params.h"
#include "RecoTracker/MkFitCore/src/V2p2Config.h"
#include "RecoTracker/MkFitCore/standalone/ConfigStandalone.h"
#include "RecoTracker/MkFitCore/standalone/V2p2Diag.h"

#include "oneapi/tbb/global_control.h"

#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>

namespace mkfit {

  namespace V2p2 = Config::V2p2;

  namespace {
    // The V2p2Params setters change every iteration config, as SetMaxCands()
    // does for maxCandsPerSeed: the next run reads them through MkJob::v2p2().
    template <typename F>
    void for_each_v2p2_params(F &&set) {
      const int ni = Config::ItrInfo.size();
      for (int i = 0; i < ni; ++i)
        set(Config::ItrInfo[i].m_v2p2);
    }

    // Held while tap mode is on, in a build whose in-event loops are parallel.
    std::unique_ptr<tbb::global_control> s_tap_thread_limit;
  }  // namespace

  //===========================================================================
  // Beam width and in-layer search
  //===========================================================================

  // Both parameter sets, not just the forward one: the beam width a candidate
  // actually gets is CombCandidate::capacity(), reserved from
  // MkJob::max_max_cands() = max(params(), params_bks()), so setting the forward
  // one alone leaves the capacity at whichever is larger and the knob does
  // nothing. That is not hypothetical -- it is what an earlier version of this
  // setter did, and cap 3 and cap 6 then came out bit-identical.
  void Shell::SetMaxCands(int fwd, int bkw) {
    if (bkw < 0)
      bkw = fwd;
    const int ni = Config::ItrInfo.size();
    for (int i = 0; i < ni; ++i) {
      Config::ItrInfo[i].m_params.maxCandsPerSeed = fwd;
      Config::ItrInfo[i].m_backward_params.maxCandsPerSeed = bkw;
    }
    printf("Shell::SetMaxCands: maxCandsPerSeed = %d fwd, %d bkw for %d iteration configs\n", fwd, bkw, ni);
  }

  // Most hits one path may take in a layer. 1 with SetMaxCands(1) is a single
  // greedy path, the closest the search gets to one best hit per layer.
  void Shell::SetMaxHitsPerLayer(int n) {
    if (n < 1 || n > V2p2Params::max_hits_per_layer_limit) {
      printf("Shell::SetMaxHitsPerLayer: %d out of range [1, %d], ignored\n", n,
             V2p2Params::max_hits_per_layer_limit);
      return;
    }
    for_each_v2p2_params([n](V2p2Params &p) { p.max_hits_per_layer = n; });
    printf("Shell::SetMaxHitsPerLayer: %d\n", n);
  }

  void Shell::SetLayerPolicy(bool wsr, bool hole_limits, bool stop_cuts) {
    V2p2::Policy::use_wsr = wsr;
    V2p2::Policy::use_hole_limits = hole_limits;
    V2p2::Policy::use_stop_cuts = stop_cuts;
    printf("Shell::SetLayerPolicy: wsr=%d hole_limits=%d stop_cuts=%d\n", (int)wsr, (int)hole_limits,
           (int)stop_cuts);
  }

  void Shell::SetLateMaxCands(int max_cands, int hits, float step_cm) {
    V2p2::InLayer::late_max_cands = max_cands;
    V2p2::InLayer::late_max_cands_hits = hits;
    V2p2::InLayer::late_max_cands_step_cm = step_cm;
    printf("Shell::SetLateMaxCands: V2p2::InLayer::late_max_cands = %d above %d found hits, except steps over "
           "%.1f cm (0 = off)\n",
           max_cands, hits, step_cm);
  }

  // The sister-hole dominance, V2p2::InLayer::sister_hole_drop and sister_hole_chi2.
  void Shell::SetSisterHole(bool on, float chi2) {
    V2p2::InLayer::sister_hole_drop = on;
    V2p2::InLayer::sister_hole_chi2 = chi2;
    printf("Shell::SetSisterHole: V2p2::InLayer::sister_hole_drop = %d, sister chi2 below %.1f (0 = any accepted)\n",
           (int)on, chi2);
  }

  // Keep one beam slot for a continuation that declined the layer. A beam policy,
  // not a score: declining is the only move that does not shrink the covariance,
  // so it is the only branch left open to an earlier hit having been wrong.
  void Shell::SetReserveHoleSlot(bool on) {
    V2p2::InLayer::reserve_hole_slot = on;
    printf("Shell::SetReserveHoleSlot: V2p2::InLayer::reserve_hole_slot = %d\n", (int)on);
  }

  // Backward-search pickup per region, named by LAYER and set through
  // SteeringParams::set_bkw_search_pickup_at_layer(), which finds the plan entry
  // covering the layer (either member of a paired entry). A layer, unlike a plan
  // index, does not move when the plan is rebuilt. -1 keeps a region's pickup.
  // Every iteration config is checked before any is changed; a layer missing
  // from a plan throws.
  void Shell::SetBkwSearchPickups(int l0, int l1, int l2, int l3, int l4) {
    const int lay[5] = {l0, l1, l2, l3, l4};
    const int ni = Config::ItrInfo.size();
    for (int i = 0; i < ni; ++i) {
      const auto &sps = Config::ItrInfo[i].m_steering_params;
      for (int r = 0; r < 5 && r < (int)sps.size(); ++r) {
        if (lay[r] >= 0 && sps[r].plan_index_of_layer(lay[r]) == -1) {
          printf("Shell::SetBkwSearchPickups: layer %d is not in the plan of iteration config %d, region %d\n",
                 lay[r], i, r);
          throw std::runtime_error("Shell::SetBkwSearchPickups: layer not in the layer plan");
        }
      }
    }
    for (int i = 0; i < ni; ++i) {
      auto &sps = Config::ItrInfo[i].m_steering_params;
      for (int r = 0; r < 5 && r < (int)sps.size(); ++r) {
        if (lay[r] >= 0)
          sps[r].set_bkw_search_pickup_at_layer(lay[r]);
      }
    }
    printf("Shell::SetBkwSearchPickups: layers %d %d %d %d %d for %d iteration configs; plan indices of config 0:",
           l0, l1, l2, l3, l4, ni);
    if (ni > 0)
      for (const auto &sp : Config::ItrInfo[0].m_steering_params)
        printf(" %d", sp.m_bkw_search_pickup);
    printf("\n");
  }

  //===========================================================================
  // Search window and hit acceptance
  //===========================================================================

  // The dq cut, split.  trk_fac multiplies dq_track (itself 3 sigma); hit_fac
  // multiplies the hit's own half-extent and has a FLOOR OF EXACTLY 1.0, below
  // which the window stops reaching the strip it is containing. extra_bins is
  // the fetch margin beyond the cut, in whole q bins.
  void Shell::SetWindowDq(float trk_fac, float hit_fac, int extra_bins) {
    V2p2::Window::dq_trk_fac = trk_fac;
    V2p2::Window::dq_hit_fac = hit_fac;
    V2p2::Window::q_extra_bins = extra_bins;
    printf("Shell::SetWindowDq: trk_fac = %.3f  hit_fac = %.3f%s  extra_bins = %d\n", trk_fac, hit_fac,
           hit_fac < 1.0f ? "  *** BELOW THE GEOMETRIC FLOOR OF 1.0 ***" : "", extra_bins);
  }

  // Phi pre-selection. The fetch range is DERIVED from the cut, so there is no
  // way to ask for a cut wider than the fetch -- the failure mode that made the
  // first version of this scan flatline for a non-physics reason.
  //   trk_fac    factor on dphi_track, applied to cut AND fetch
  //   flat_rad   the flat per-hit tolerance, RADIANS (one phi bin = 0.024544);
  //              used only with the per-hit extent off, see SetWindowPhiPerHit()
  //   extra_bins fetch safety margin beyond the cut, in WHOLE bins
  void Shell::SetWindowDphi(float trk_fac, float flat_rad, int extra_bins) {
    V2p2::Window::dphi_trk_fac = trk_fac;
    V2p2::Window::dphi_flat_rad = flat_rad;
    V2p2::Window::phi_extra_bins = extra_bins;
    printf("Shell::SetWindowDphi: trk_fac = %.3f  flat_rad = %.5f rad (%.3f bins)  extra_bins = %d\n", trk_fac,
           flat_rad, flat_rad / (2.0f * float(M_PI) / 256.0f), extra_bins);
  }

  // Phi cut from the hit's OWN covariance instead of the flat half-bin constant.
  // fac means what SetWindowDq()'s hit_fac means: 1.0 exactly contains the hit's
  // extent. Expect a large tightening -- ~380x in TB2S -- so scan it, do not
  // assume it.
  void Shell::SetWindowPhiPerHit(bool on, float fac) {
    V2p2::Window::phi_per_hit = on;
    V2p2::Window::dphi_hit_fac = fac;
    printf("Shell::SetWindowPhiPerHit: %s  fac = %.2f%s\n", on ? "ON (per-hit covariance)" : "off (flat)", fac,
           (on && fac < 1.0f) ? "  *** BELOW THE CONTAINMENT FLOOR ***" : "");
  }

  // Per-hit surface reference, using the HIT'S OWN MODULE NORMAL. This is the
  // real fix; MkBins::surface_reference_dq (SetMkBinsSurfaceQ) is the
  // layer-cylinder scaffold that proved the mechanism and over-widens tilted
  // TBPS ~8x.
  void Shell::SetWindowSurfaceQ(bool on) {
    V2p2::Window::surface_q = on;
    printf("Shell::SetWindowSurfaceQ: V2p2::Window::surface_q = %d\n", (int)on);
  }

  // Reference the pre-selection q window to the layer SURFACE rather than to the
  // fixed path length errPropFromPathL_impl transports the covariance to. See
  // MkBins::surface_reference_dq(). Changes dq_track, which the fetch, the WSR,
  // the line pre-cut and the rho region read -- the Kalman update already
  // carries the term, via jacCurv2Loc's cosz.
  void Shell::SetMkBinsSurfaceQ(bool on) {
    V2p2::Diag::mkbins_surface_q = on;
    printf("Shell::SetMkBinsSurfaceQ: V2p2::Diag::mkbins_surface_q = %d\n", (int)on);
  }

  void Shell::SetMiniPropParamB(int mode) {
    V2p2::Window::mini_prop_param_b = mode;
    printf("Shell::SetMiniPropParamB: %d (0 uniform Config::Bfield, 1 bFieldFromZR at each state, 2 as the "
           "Kalman)\n",
           mode);
  }

  void Shell::SetPreCut(bool q, bool phi) {
    V2p2::PreCut::q = q;
    V2p2::PreCut::phi = phi;
    printf("Shell::SetPreCut: q %s, phi %s\n", q ? "on" : "off", phi ? "on" : "off");
  }

  void Shell::SetHitChi2(float cut, float trk_fac) {
    for_each_v2p2_params([cut, trk_fac](V2p2Params &p) {
      p.hit_chi2_cut = cut;
      p.chi2_trk_fac = trk_fac;
    });
    printf("Shell::SetHitChi2: cut %g, trk_fac %g\n", cut, trk_fac);
  }

  //===========================================================================
  // Score and track scorers
  //===========================================================================

  // The three hit_eff setters each set one thing: SetScoreMode() the mode only,
  // SetHitEff() the base hit_eff of each direction, SetHitEffGroups() the forward
  // per-group values, of which a negative one takes the base.
  void Shell::SetScoreMode(int mode) {
    V2p2::Score::mode = mode;
    printf("Shell::SetScoreMode: V2p2::Score::mode = %d (0 linear, 1 loglh)\n", mode);
  }

  void Shell::SetHitEff(float fwd, float bkw) {
    for_each_v2p2_params([fwd, bkw](V2p2Params &p) {
      p.score_outward.hit_eff = fwd;
      p.score_inward.hit_eff = bkw;
    });
    printf("Shell::SetHitEff: hit_eff fwd %g, bkw %g (per-group values unchanged)\n", fwd, bkw);
  }

  void Shell::SetHitEffGroups(float pix_brl, float pix_ec, float strip_brl_ps, float strip_brl_2s, float strip_ec) {
    const float g[5] = {pix_brl, pix_ec, strip_brl_ps, strip_brl_2s, strip_ec};
    for_each_v2p2_params([&g](V2p2Params &p) {
      for (int i = 0; i < 5; ++i)
        p.score_outward.hit_eff_grp[i] = g[i];
    });
    printf("Shell::SetHitEffGroups: fwd hit_eff pixB %g pixE %g TBPS %g TB2S %g TEDD %g (negative = hit_eff %g)\n",
           g[0], g[1], g[2], g[3], g[4], Config::ItrInfo[0].m_v2p2.score_outward.hit_eff);
  }

  // Fix one term of the log-likelihood at a constant, at fixed eps. Passing the term's own
  // MEAN as the constant removes its variation and nothing else -- see the note
  // in V2p2Score.h for why that is the only comparison that can attribute a
  // regional effect to rho rather than to eps.
  void Shell::SetScoreTerms(bool use_rho, float rho_const, bool use_detv, float detv_const) {
    V2p2::Score::use_rho = use_rho;
    V2p2::Score::rho_const = rho_const;
    V2p2::Score::use_detv = use_detv;
    V2p2::Score::detv_const = detv_const;
    printf("Shell::SetScoreTerms: rho %s (const %.4f), det V %s (const %.4f)\n", use_rho ? "PER STEP" : "FLAT",
           rho_const, use_detv ? "PER HIT" : "FLAT", detv_const);
  }

  // The score of Score::mode 0. miss_fwd / miss_bkw are the head-body
  // asymmetry: outward a trailing hole is at large radius and cheap, inward it
  // is at small radius and is the most expensive hole there is.
  void Shell::SetLinearScore(float hit_bonus, float chi2_weight, float miss_fwd, float miss_bkw) {
    for_each_v2p2_params([=](V2p2Params &p) {
      p.score_outward.hit_bonus = p.score_inward.hit_bonus = hit_bonus;
      p.score_outward.chi2_weight = p.score_inward.chi2_weight = chi2_weight;
      p.score_outward.miss_penalty = miss_fwd;
      p.score_inward.miss_penalty = miss_bkw;
    });
    printf("Shell::SetLinearScore: hit_bonus=%.2f chi2_weight=%.2f miss_fwd=%.2f miss_bkw=%.2f "
           "(hole beats a hit above chi2 = %.1f inward)\n",
           hit_bonus, chi2_weight, miss_fwd, miss_bkw, chi2_weight > 0 ? (hit_bonus + miss_bkw) / chi2_weight : 0.0f);
  }

  // Default track scorer of every iteration config, by registered name (e.g.
  // "phase1:default", which the CMSSW phase-2 initialStep JSON and, since
  // 2026-10-04, the phase-2 plugin name, or "phase2:LstIntoPix", the plugin's
  // earlier default). Re-resolves the function pointers.
  void Shell::SetTrackScorer(const char *name) {
    const int ni = Config::ItrInfo.size();
    for (int i = 0; i < ni; ++i) {
      Config::ItrInfo[i].m_default_track_scorer_name = name;
      Config::ItrInfo[i].setupStandardFunctionsFromNames();
    }
    printf("Shell::SetTrackScorer: %s for %d iteration configs\n", name, ni);
  }

  // Per-task track scorer of every iteration config, by registered name; task is
  // "final_pick", "post_bkfit" or "duplicate_cleaner", and "" restores the fallback
  // (the search scorer, or no re-scoring before the duplicate cleaner).
  void Shell::SetTrackScorerTask(const char *task, const char *name) {
    const std::string t(task);
    const int ni = Config::ItrInfo.size();
    for (int i = 0; i < ni; ++i) {
      auto &ic = Config::ItrInfo[i];
      if (t == "final_pick")
        ic.m_final_pick_track_scorer_name = name;
      else if (t == "post_bkfit")
        ic.m_post_bkfit_track_scorer_name = name;
      else if (t == "duplicate_cleaner")
        ic.m_duplicate_cleaner_track_scorer_name = name;
      else {
        printf("Shell::SetTrackScorerTask: unknown task '%s'\n", task);
        return;
      }
      ic.setupStandardFunctionsFromNames();
    }
    printf("Shell::SetTrackScorerTask: %s = '%s' for %d iteration configs\n", task, name, ni);
  }

  // Final pick by "v2p2:llh" (on) or by the search scorer (off), with its two
  // penalties. The penalties are process-wide knobs, the scorer name is set per
  // iteration config.
  void Shell::SetFinalPickLlh(bool on, float hole_penalty, float tail_penalty) {
    V2p2::Policy::final_pick_hole_penalty = hole_penalty;
    V2p2::Policy::final_pick_tail_penalty = tail_penalty;
    SetTrackScorerTask("final_pick", on ? "v2p2:llh" : "");
    printf("Shell::SetFinalPickLlh: %d, layer-step score - %g per inside hole - %g per tail hole\n", (int)on,
           hole_penalty, tail_penalty);
  }

  //===========================================================================
  // Material, backward fit, duplicate cleaner
  //===========================================================================

  // Toggle material on the SEARCH's inter-layer propagation, the Kalman step to
  // the hit's module plane.
  void Shell::SetSearchMaterial(bool on) {
    auto &pc = const_cast<PropagationConfig &>(Config::TrkInfo.prop_config());
    pc.finding_inter_layer_pflags.apply_material = on;
    printf("Shell::SetSearchMaterial: finding_inter_layer_pflags.apply_material = %d\n", (int)on);
  }

  // Toggle PF_apply_material on the backward-fit propagation only. Diagnostic:
  // material inflates the propagated covariance, so turning it off must RAISE
  // chi2; how much says how big the material term is.
  void Shell::SetBkFitMaterial(bool on) {
    auto &pc = const_cast<PropagationConfig &>(Config::TrkInfo.prop_config());
    pc.backward_fit_pflags.apply_material = on;
    printf("Shell::SetBkFitMaterial: backward_fit_pflags.apply_material = %d\n", (int)on);
  }

  // Outlier rejection in the backward fit (v2p2's and MkFinder's), m_backward_fit_outlier_*
  // of every iteration config (chi2 0 = off).
  void Shell::SetBkFitOutliers(float chi2, int max_outliers, float min_pt) {
    const int ni = Config::ItrInfo.size();
    for (int i = 0; i < ni; ++i) {
      Config::ItrInfo[i].m_backward_fit_outlier_chi2 = chi2;
      Config::ItrInfo[i].m_backward_fit_max_outliers = max_outliers;
      Config::ItrInfo[i].m_backward_fit_outlier_min_pt = min_pt;
    }
    printf("Shell::SetBkFitOutliers: chi2 > %.1f, at most %d per track, pT > %.2f (chi2 0 = off) for %d iteration "
           "configs\n",
           chi2, max_outliers, min_pt, ni);
  }

  // The duplicate cleaner's dR rule thresholds of every iteration config (dc_drth_central /
  // _obarrel / _forward, as set_dupl_params()); 0 switches the rule off. The shared-hit rule
  // is not touched.
  void Shell::SetDuplDrth(float central, float obarrel, float forward) {
    const int ni = Config::ItrInfo.size();
    for (int i = 0; i < ni; ++i) {
      Config::ItrInfo[i].dc_drth_central = central;
      Config::ItrInfo[i].dc_drth_obarrel = obarrel;
      Config::ItrInfo[i].dc_drth_forward = forward;
    }
    printf("Shell::SetDuplDrth: dR rule thresholds %g / %g / %g for %d iteration configs\n", central, obarrel,
           forward, ni);
  }

  //===========================================================================
  // HLT seeds
  //===========================================================================

  void Shell::SetHltSeeds(int kind, bool chop_pixels) {
    s_hlt_seed_kind = (HltSeedKind_e)kind;
    s_hlt_chop_pixels = chop_pixels;
    printf("Shell::SetHltSeeds: kind=%s chop_pixels=%d\n", hlt_seed_kind_name(s_hlt_seed_kind), (int)chop_pixels);
  }

  //===========================================================================
  // Taps and tap mode
  //===========================================================================

  // End-of-search beam purity, Diag::final_beam_purity. Switching it on clears
  // what it collected before; probe_final_beam_report() prints it.
  void Shell::SetTapFinalBeam(bool on) {
    V2p2::Diag::final_beam_purity = on;
    if (on)
      v2p2_final_beam_diag_reset();
    printf("Shell::SetTapFinalBeam: %d (end-of-search beam purity%s)\n", (int)on, on ? ", reset" : "");
    m_tap_final_beam = on;
    update_tap_mode();
  }

  // Running sums of the likelihood terms, Config::V2p2::ScoreStats, from which
  // probe_score_term_stats_report() gives the constants for SetScoreTerms().
  // Switching it on zeroes the sums.
  void Shell::SetTapScoreStats(bool on) {
    if (on) {
      V2p2::ScoreStats::n_hits = 0;
      V2p2::ScoreStats::sum_log_rho = V2p2::ScoreStats::sum_log_detv = 0.0;
    }
    V2p2::ScoreStats::accum = on;
    printf("Shell::SetTapScoreStats: %d%s\n", (int)on, on ? " (sums zeroed)" : "");
    m_tap_score_stats = on;
    update_tap_mode();
  }

  // The taps accumulate into unsynchronised globals keyed by seed, without the
  // event, so they are correct only with one thread inside the event and one
  // event in flight. Tap mode guarantees both: collect_events_std() and
  // collect_events_hlt() clamp the events in flight to one, and in a build whose
  // in-event loops are parallel (no TBB_DEBUG, i.e. no EVENT_RDF_TRACE) a
  // tbb::global_control limits TBB to one thread. The shell's own limit is the
  // one mkFit.cc sets from --num-thr.
  void Shell::update_tap_mode() {
    const bool mode = m_tap_final_beam || m_tap_score_stats;
    if (mode == m_tap_mode)
      return;
    m_tap_mode = mode;
#ifdef TBB_DEBUG
    printf("Shell: tap mode %s; the in-event loops are serial in this build (TBB_DEBUG).\n", mode ? "on" : "off");
#else
    if (mode) {
      s_tap_thread_limit = std::make_unique<tbb::global_control>(tbb::global_control::max_allowed_parallelism, 1);
      printf("Shell: tap mode on; TBB limited to one thread, so the search inside an event runs on one thread.\n");
    } else {
      s_tap_thread_limit.reset();
      printf("Shell: tap mode off; the TBB thread limit is back to --num-thr (%d).\n", Config::numThreadsFinder);
    }
#endif
  }

}  // namespace mkfit
