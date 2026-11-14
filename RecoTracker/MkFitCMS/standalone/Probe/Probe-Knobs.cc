#include "RecoTracker/MkFitCMS/standalone/Probe/Probe.h"

#include "RecoTracker/MkFitCore/interface/Config.h"
#include "RecoTracker/MkFitCore/interface/IterationConfig.h"
#include "RecoTracker/MkFitCore/interface/PropagationConfig.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/MkFitCore/interface/V2p2Params.h"
#include "RecoTracker/MkFitCore/src/V2p2Config.h"
#include "RecoTracker/MkFitCore/standalone/ConfigStandalone.h"

#include <cstdio>
#include <stdexcept>
#include <string>

namespace mkfit {
  namespace V2p2 = Config::V2p2;

  // Both parameter sets, not just the forward one: the beam width a candidate
  // actually gets is CombCandidate::capacity(), reserved from
  // MkJob::max_max_cands() = max(params(), params_bks()), so setting the forward
  // one alone leaves the capacity at whichever is larger and the knob does
  // nothing. That is not hypothetical -- it is what an earlier version of this
  // function did, and cap 3 and cap 6 then came out bit-identical.
  void val_max_cands(int n) {
    const int ni = Config::ItrInfo.size();
    for (int i = 0; i < ni; ++i) {
      Config::ItrInfo[i].m_params.maxCandsPerSeed = n;
      Config::ItrInfo[i].m_backward_params.maxCandsPerSeed = n;
    }
    printf("val_max_cands: maxCandsPerSeed = %d (fwd and bkw) for %d iteration configs\n", n, ni);
  }

  // Outlier rejection in the backward fit (v2p2's and MkFinder's), m_backward_fit_outlier_*
  // of every iteration config (chi2 0 = off).
  void val_bkfit_outliers(float chi2, int max_outliers, float min_pt) {
    const int ni = Config::ItrInfo.size();
    for (int i = 0; i < ni; ++i) {
      Config::ItrInfo[i].m_backward_fit_outlier_chi2 = chi2;
      Config::ItrInfo[i].m_backward_fit_max_outliers = max_outliers;
      Config::ItrInfo[i].m_backward_fit_outlier_min_pt = min_pt;
    }
    printf("val_bkfit_outliers: chi2 > %.1f, at most %d per track, pT > %.2f (chi2 0 = off) for %d iteration configs\n",
           chi2, max_outliers, min_pt, ni);
  }

  // Forward and backward beam widths separately (val_max_cands sets both).
  void val_max_cands_fb(int n_fwd, int n_bkw) {
    const int ni = Config::ItrInfo.size();
    for (int i = 0; i < ni; ++i) {
      Config::ItrInfo[i].m_params.maxCandsPerSeed = n_fwd;
      Config::ItrInfo[i].m_backward_params.maxCandsPerSeed = n_bkw;
    }
    printf("val_max_cands_fb: maxCandsPerSeed = %d fwd, %d bkw for %d iteration configs\n", n_fwd, n_bkw, ni);
  }

  // Default track scorer of every iteration config, by registered name (e.g.
  // "phase1:default", which the CMSSW phase-2 initialStep JSON and, since
  // 2026-10-04, the phase-2 plugin name, or "phase2:LstIntoPix", the plugin's
  // earlier default). Re-resolves the function pointers.
  void val_track_scorer(const char *name) {
    const int ni = Config::ItrInfo.size();
    for (int i = 0; i < ni; ++i) {
      Config::ItrInfo[i].m_default_track_scorer_name = name;
      Config::ItrInfo[i].setupStandardFunctionsFromNames();
    }
    printf("val_track_scorer: %s for %d iteration configs\n", name, ni);
  }

  // Per-task track scorer of every iteration config, by registered name; task is
  // "final_pick", "post_bkfit" or "duplicate_cleaner", and "" restores the fallback
  // (the search scorer, or no re-scoring before the duplicate cleaner).
  // The duplicate cleaner's dR rule thresholds of every iteration config (dc_drth_central /
  // _obarrel / _forward, as set_dupl_params()); 0 switches the rule off. The shared-hit rule
  // is not touched.
  void val_dupl_drth(float central, float obarrel, float forward) {
    const int ni = Config::ItrInfo.size();
    for (int i = 0; i < ni; ++i) {
      Config::ItrInfo[i].dc_drth_central = central;
      Config::ItrInfo[i].dc_drth_obarrel = obarrel;
      Config::ItrInfo[i].dc_drth_forward = forward;
    }
    printf("val_dupl_drth: dR rule thresholds %g / %g / %g for %d iteration configs\n", central, obarrel, forward,
           ni);
  }

  void val_track_scorer_task(const char *task, const char *name) {
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
        printf("val_track_scorer_task: unknown task '%s'\n", task);
        return;
      }
      ic.setupStandardFunctionsFromNames();
    }
    printf("val_track_scorer_task: %s = '%s' for %d iteration configs\n", task, name, ni);
  }

  // Backward-search pickup of iteration 0, as a PLAN INDEX per region (the
  // form the JSON stores; the plugin derives it from a layer). -1 keeps it.
  void val_bkw_pickups(int r0, int r1, int r2, int r3, int r4) {
    auto &sps = Config::ItrInfo[0].m_steering_params;
    const int v[5] = {r0, r1, r2, r3, r4};
    for (int r = 0; r < 5 && r < (int)sps.size(); ++r) {
      if (v[r] < 0)
        continue;
      if (v[r] >= (int)sps[r].m_layer_plan.size())
        throw std::runtime_error("val_bkw_pickups: plan index out of range");
      sps[r].m_bkw_search_pickup = v[r];
    }
    printf("val_bkw_pickups: %d %d %d %d %d\n", r0, r1, r2, r3, r4);
  }

  // Toggle material on the SEARCH's inter-layer propagation, the Kalman step to
  // the hit's module plane.
  void val_search_material(bool on) {
    auto &pc = const_cast<PropagationConfig&>(Config::TrkInfo.prop_config());
    pc.finding_inter_layer_pflags.apply_material = on;
    printf("val_search_material: finding_inter_layer_pflags.apply_material = %d\n", (int) on);
  }

  // Reference the pre-selection q window to the layer SURFACE rather than to the
  // fixed path length errPropFromPathL_impl transports the covariance to. See
  // MkBins::surface_reference_dq(). Affects the WINDOW AND THE dq CUT ONLY --
  // the Kalman update already carries the term, via jacCurv2Loc's cosz.
  void val_mkbins_surface_q(bool on) {
    V2p2::Diag::mkbins_surface_q = on;
    printf("val_mkbins_surface_q: V2p2::Diag::mkbins_surface_q = %d\n", (int) on);
  }

  // The dq cut, split.  trk_fac multiplies dq_track (itself 3 sigma); hit_fac
  // multiplies the hit's own half-extent and has a FLOOR OF EXACTLY 1.0, below
  // which the window stops reaching the strip it is containing. extra_bins is
  // the fetch margin beyond the cut, in whole q bins.
  void val_dq(float trk_fac, float hit_fac, int extra_bins) {
    V2p2::Window::dq_trk_fac   = trk_fac;
    V2p2::Window::dq_hit_fac   = hit_fac;
    V2p2::Window::q_extra_bins = extra_bins;
    printf("val_dq: trk_fac = %.3f  hit_fac = %.3f%s  extra_bins = %d\n",
           trk_fac, hit_fac,
           hit_fac < 1.0f ? "  *** BELOW THE GEOMETRIC FLOOR OF 1.0 ***" : "", extra_bins);
  }

  // Phi cut from the hit's OWN covariance instead of the flat half-bin constant.
  // fac means what val_dq's hit_fac means: 1.0 exactly contains the hit's extent.
  // Expect a large tightening -- ~380x in TB2S -- so scan it, do not assume it.
  void val_phi_per_hit(bool on, float fac) {
    V2p2::Window::phi_per_hit  = on;
    V2p2::Window::dphi_hit_fac = fac;
    printf("val_phi_per_hit: %s  fac = %.2f%s\n", on ? "ON (per-hit covariance)" : "off (flat)",
           fac, (on && fac < 1.0f) ? "  *** BELOW THE CONTAINMENT FLOOR ***" : "");
  }

  // The V2p2Params setters below change every iteration config, as val_max_cands()
  // does for maxCandsPerSeed: the next run reads them through MkJob::v2p2().
  template <typename F>
  void for_each_v2p2_params(F &&set) {
    const int ni = Config::ItrInfo.size();
    for (int i = 0; i < ni; ++i)
      set(Config::ItrInfo[i].m_v2p2);
  }

  // Phi pre-selection. The fetch range is DERIVED from the cut, so there is no
  // way to ask for a cut wider than the fetch -- the failure mode that made the
  // first version of this scan flatline for a non-physics reason.
  //   trk_fac    factor on dphi_track, applied to cut AND fetch
  //   hit_rad    the flat per-hit tolerance, RADIANS (one phi bin = 0.024544);
  //              used only with the per-hit extent off, see val_phi_per_hit()
  //   extra_bins fetch safety margin beyond the cut, in WHOLE bins
  void val_hit_chi2(float cut, float trk_fac) {
    for_each_v2p2_params([cut, trk_fac](V2p2Params &p) {
      p.hit_chi2_cut = cut;
      p.chi2_trk_fac = trk_fac;
    });
    printf("val_hit_chi2: cut %g, trk_fac %g\n", cut, trk_fac);
  }

  void val_precut(bool q, bool phi) {
    V2p2::PreCut::q = q;
    V2p2::PreCut::phi = phi;
    printf("val_precut: q %s, phi %s\n", q ? "on" : "off", phi ? "on" : "off");
  }

  void val_dphi(float trk_fac, float hit_rad, int extra_bins) {
    V2p2::Window::dphi_trk_fac   = trk_fac;
    V2p2::Window::dphi_flat_rad  = hit_rad;
    V2p2::Window::phi_extra_bins = extra_bins;
    printf("val_dphi: trk_fac = %.3f  hit_rad = %.5f rad (%.3f bins)  extra_bins = %d\n",
           trk_fac, hit_rad, hit_rad / (2.0f * float(M_PI) / 256.0f), extra_bins);
  }

  // Per-hit surface reference, using the HIT'S OWN MODULE NORMAL. This is the
  // real fix; MkBins::surface_reference_dq (val_surf_q) is the layer-cylinder
  // scaffold that proved the mechanism and over-widens tilted TBPS ~8x.
  void val_layer_policy(bool wsr, bool hole_limits, bool stop_cuts) {
    V2p2::Policy::use_wsr = wsr;
    V2p2::Policy::use_hole_limits = hole_limits;
    V2p2::Policy::use_stop_cuts = stop_cuts;
    printf("val_layer_policy: wsr=%d hole_limits=%d stop_cuts=%d\n",
           (int) wsr, (int) hole_limits, (int) stop_cuts);
  }

  // Keep one beam slot for a continuation that declined the layer. A beam policy,
  // not a score: declining is the only move that does not shrink the covariance,
  // so it is the only branch left open to an earlier hit having been wrong.
  void val_late_max_cands(int max_cands, int hits, float step_cm) {
    V2p2::InLayer::late_max_cands = max_cands;
    V2p2::InLayer::late_max_cands_hits = hits;
    V2p2::InLayer::late_max_cands_step_cm = step_cm;
    printf("val_late_max_cands: V2p2::InLayer::late_max_cands = %d above %d found hits, except steps over %.1f cm"
           " (0 = off)\n", max_cands, hits, step_cm);
  }

  // The sister-hole dominance, V2p2::InLayer::sister_hole_drop and sister_hole_chi2.
  void val_sister_hole(bool on, float chi2) {
    V2p2::InLayer::sister_hole_drop = on;
    V2p2::InLayer::sister_hole_chi2 = chi2;
    printf("val_sister_hole: V2p2::InLayer::sister_hole_drop = %d, sister chi2 below %.1f (0 = any accepted)\n",
           (int) on, chi2);
  }

  void val_reserve_hole_slot(bool on) {
    V2p2::InLayer::reserve_hole_slot = on;
    printf("val_reserve_hole_slot: V2p2::InLayer::reserve_hole_slot = %d\n", (int) on);
  }

  void val_hit_eff(float fwd, float bkw) {
    for_each_v2p2_params([fwd, bkw](V2p2Params &p) {
      p.score_outward.hit_eff = fwd;
      p.score_inward.hit_eff = bkw;
    });
    printf("val_hit_eff: hit_eff fwd %g, bkw %g\n", fwd, bkw);
  }

  void val_hit_eff_grp(float pix_brl, float pix_ec, float strip_brl_ps, float strip_brl_2s, float strip_ec) {
    const float g[5] = {pix_brl, pix_ec, strip_brl_ps, strip_brl_2s, strip_ec};
    for_each_v2p2_params([&g](V2p2Params &p) {
      for (int i = 0; i < 5; ++i)
        p.score_outward.hit_eff_grp[i] = g[i];
    });
    printf("val_hit_eff_grp: fwd hit_eff pixB %g pixE %g TBPS %g TB2S %g TEDD %g (negative = hit_eff %g)\n",
           g[0], g[1], g[2], g[3], g[4], Config::ItrInfo[0].m_v2p2.score_outward.hit_eff);
  }

  // Final pick by "v2p2:llh" (on) or by the search scorer (off), with its two penalties.
  void val_final_pick_llh(bool on, float hole_penalty, float tail_penalty) {
    V2p2::Policy::final_pick_hole_penalty = hole_penalty;
    V2p2::Policy::final_pick_tail_penalty = tail_penalty;
    val_track_scorer_task("final_pick", on ? "v2p2:llh" : "");
    printf("val_final_pick_llh: %d, layer-step score - %g per inside hole - %g per tail hole\n", (int) on,
           hole_penalty, tail_penalty);
  }

  void val_mini_prop_param_b(int on) {
    V2p2::Window::mini_prop_param_b = on;
    printf("val_mini_prop_param_b: %d (0 uniform Config::Bfield, 1 bFieldFromZR at each state, 2 as the Kalman)\n", on);
  }


  // Fix one term of the log-likelihood at a constant, at fixed eps. Passing the term's own
  // MEAN as the constant removes its variation and nothing else -- see the note
  // in V2p2Score.h for why that is the only comparison that can attribute a
  // regional effect to rho rather than to eps.
  void val_score_terms(bool use_rho, float rho_const, bool use_detv, float detv_const) {
    V2p2::Score::use_rho = use_rho;      V2p2::Score::rho_const = rho_const;
    V2p2::Score::use_detv = use_detv;    V2p2::Score::detv_const = detv_const;
    printf("val_score_terms: rho %s (const %.4f), det V %s (const %.4f)\n",
           use_rho ? "PER STEP" : "FLAT", rho_const,
           use_detv ? "PER HIT" : "FLAT", detv_const);
  }

  void val_score_mode(int mode, float hit_eff) {
    V2p2::Score::mode = mode;
    for_each_v2p2_params([hit_eff](V2p2Params &p) { p.score_outward.hit_eff = p.score_inward.hit_eff = hit_eff; });
    printf("val_score_mode: V2p2::Score::mode = %d (0 linear, 1 loglh), hit_eff = %g\n",
           mode, hit_eff);
  }

  // Most hits one path may take in a layer. 1 with val_max_cands(1) is a single
  // greedy path, the closest the search gets to one best hit per layer.
  void val_max_hits_per_layer(int n) {
    if (n < 1 || n > V2p2Params::max_hits_per_layer_limit) {
      printf("val_max_hits_per_layer: %d out of range [1, %d], ignored\n", n, V2p2Params::max_hits_per_layer_limit);
      return;
    }
    for_each_v2p2_params([n](V2p2Params &p) { p.max_hits_per_layer = n; });
    printf("val_max_hits_per_layer: %d\n", n);
  }

  // miss_fwd / miss_bkw are the head-body asymmetry: outward a trailing hole is
  // at large radius and cheap, inward it is at small radius and is the most
  // expensive hole there is.
  void val_score(float hit_bonus, float chi2_weight, float miss_fwd, float miss_bkw) {
    for_each_v2p2_params([=](V2p2Params &p) {
      p.score_outward.hit_bonus = p.score_inward.hit_bonus = hit_bonus;
      p.score_outward.chi2_weight = p.score_inward.chi2_weight = chi2_weight;
      p.score_outward.miss_penalty = miss_fwd;
      p.score_inward.miss_penalty = miss_bkw;
    });
    printf("val_score: hit_bonus=%.2f chi2_weight=%.2f miss_fwd=%.2f miss_bkw=%.2f "
           "(hole beats a hit above chi2 = %.1f inward)\n",
           hit_bonus, chi2_weight, miss_fwd, miss_bkw,
           chi2_weight > 0 ? (hit_bonus + miss_bkw) / chi2_weight : 0.0f);
  }

  void val_surf_q_hit(bool on) {
    V2p2::Window::surface_q = on;
    printf("val_surf_q_hit: V2p2::Window::surface_q = %d\n", (int) on);
  }

  // Toggle PF_apply_material on the backward-fit propagation only. Returns the
  // previous state. Diagnostic: material inflates the propagated covariance, so
  // turning it off must RAISE chi2; how much says how big the material term is.
  void val_bkfit_material(bool on) {
    auto &pc = const_cast<PropagationConfig&>(Config::TrkInfo.prop_config());
    pc.backward_fit_pflags.apply_material = on;
    printf("val_bkfit_material: backward_fit_pflags.apply_material = %d\n", (int) on);
  }

}  // namespace mkfit
