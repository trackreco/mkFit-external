// Compiled only in a trace build (EVENT_RDF_TRACE, which defines MKFIT_TRACE): every probe here
// reads the Tr* records the v2p2 search writes into the Event.
#ifdef MKFIT_TRACE

#include "RecoTracker/MkFitCMS/standalone/Probe/Probe.h"
#include "RecoTracker/MkFitCore/standalone/DataFormats/ProbeStructs.h"
#include "RecoTracker/MkFitCore/standalone/DataFormats/RntConversions.h"

#include "RecoTracker/MkFitCore/interface/Config.h"
#include "RecoTracker/MkFitCore/interface/cms_common_macros.h"
#include "RecoTracker/MkFitCore/interface/Hit.h"
#include "RecoTracker/MkFitCore/interface/IterationConfig.h"
#include "RecoTracker/MkFitCore/interface/Track.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/MkFitCore/src/V2p2Config.h"
#include "RecoTracker/MkFitCore/standalone/ConfigStandalone.h"
#include "RecoTracker/MkFitCore/standalone/Event.h"
#include "RecoTracker/MkFitCore/standalone/TrackExtra.h"

#include "TFile.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TString.h"
#include "TTree.h"
#include "TVectorD.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <functional>
#include <map>
#include <vector>

namespace mkfit {

  // Lite trace: keep only pixel-barrel hits and skip the exact double-precision
  // helix-plane solve. That solve plus the sheer volume (~142k records/event) is
  // what makes probe_search_event ~45 s/event; with this it is seconds, which is
  // what makes a covariance SCAN affordable.
  static int g_probe_search_lite = 0;   // bit0 skip exact solve, bit1 pixel barrel only

  // bit 0 = skip the exact double-precision helix-plane solve (the expensive
  // part); bit 1 = keep only pixel-barrel hits. 3 = both, which is what a
  // pixel-barrel covariance scan wants; 1 = whole detector but still fast,
  // which is what a phi-everywhere measurement wants.
  void probe_search_lite(int mode) {
    g_probe_search_lite = mode;
    printf("probe_search_lite: mode %d (skip_solve=%d pixb_only=%d)\n",
           mode, mode & 1, (mode >> 1) & 1);
  }

  // ==========================================================================
  // Flatten the MKFIT_TRACE graph of ONE event's hit matches into a TTree.
  //
  // The trace is a graph of parent-linked records; everything below is a join,
  // not a measurement: TrHitMatch -> TrLayerSearch (for the track's own sigmas)
  // and TrHitMatch -> TrKalmanUpdate (for chi2 and the module-frame residual).
  // Both ids can legitimately be -1 on an untraced path, and every tr_*()
  // accessor is an UNCHECKED vector index, so each is guarded.
  namespace {
    std::vector<ProbeSearchHit> g_sh;
    std::vector<ProbeSearchMiss> g_sm;
    std::vector<ProbeCovStep>    g_cs;

    // Project a 6x6 CCS covariance onto the search's own coordinates. Reading
    // err(2,2) as "sigma_q" is only right in the barrel; in the endcap q is r
    // and needs the full transverse block.
    [[maybe_unused]] void project_sigma(const mkfit::TrackState &ts, bool barrel,
                                        float &sig_q, float &sig_phi) {
      const double x = ts.parameters[0], y = ts.parameters[1];
      const double r2 = x*x + y*y;
      const double exx = ts.errors.At(0,0), eyy = ts.errors.At(1,1), exy = ts.errors.At(0,1);
      sig_q = barrel ? (float) std::sqrt(std::max(0.0, (double) ts.errors.At(2,2)))
                     : (r2 > 0 ? (float) std::sqrt(std::max(0.0,
                         (x*x*exx + 2*x*y*exy + y*y*eyy) / r2)) : -1.f);
      sig_phi = r2 > 0 ? (float) std::sqrt(std::max(0.0,
                  (y*y*exx - 2*x*y*exy + x*x*eyy) / (r2*r2))) : -1.f;
    }
    int g_sh_evt = -1;
  }


  namespace {
    // Exact crossing of a uniform-B helix with a plane, in double.
    // Substituting the helix into n.(x(a) - P) = 0 gives
    //   f(a) = A sin a + B (1 - cos a) + C a + D
    // which is transcendental unless n_z p_z = 0, hence scan-and-bisect rather
    // than a closed form. Returns the number of roots found in [-lim, lim].
    struct PlaneRoots {
      double a[8]; int n = 0;
    };
    PlaneRoots exact_plane_roots(const double x0[3], const double p0[3], double k,
                                 const double nrm[3], const double pnt[3], double lim) {
      const double A =  k * (nrm[0]*p0[0] + nrm[1]*p0[1]);
      const double B =  k * (nrm[1]*p0[0] - nrm[0]*p0[1]);
      const double C =  k * nrm[2]*p0[2];
      const double D =  nrm[0]*(x0[0]-pnt[0]) + nrm[1]*(x0[1]-pnt[1]) + nrm[2]*(x0[2]-pnt[2]);
      auto fv = [&](double a){ return A*std::sin(a) + B*(1.0-std::cos(a)) + C*a + D; };
      PlaneRoots R;
      const int N = 8000;
      double ap = -lim, fp = fv(-lim);
      for (int i = 1; i <= N && R.n < 8; ++i) {
        const double aa = -lim + 2.0*lim*i/N, ff = fv(aa);
        if ((fp < 0.0) != (ff < 0.0)) {
          double lo = ap, hi = aa, flo = fp;
          for (int it = 0; it < 100; ++it) {
            const double m = 0.5*(lo+hi), fm = fv(m);
            if ((flo < 0.0) == (fm < 0.0)) { lo = m; flo = fm; } else hi = m;
          }
          R.a[R.n++] = 0.5*(lo+hi);
        }
        ap = aa; fp = ff;
      }
      return R;
    }
  }

  namespace {
    // Exact crossing of the helix with a cylinder r = R (barrel) or a plane
    // z = Z (endcap), in double, by the same scan-and-bisect. Returns the root
    // of smallest |alpha|, or -999 if the surface is not reached.
    double exact_surface_alpha(const double x0[3], const double p0[3], double k,
                               double target, bool barrel, double lim) {
      auto fv = [&](double a) {
        const double sa = std::sin(a), ca = std::cos(a);
        if (!barrel) return x0[2] + k*p0[2]*a - target;
        const double X = x0[0] + k*(p0[0]*sa - p0[1]*(1.0-ca));
        const double Y = x0[1] + k*(p0[1]*sa + p0[0]*(1.0-ca));
        return std::hypot(X, Y) - target;
      };
      const int N = 4000;
      double best = -999.0, ap = -lim, fp = fv(-lim);
      for (int i = 1; i <= N; ++i) {
        const double aa = -lim + 2.0*lim*i/N, ff = fv(aa);
        if ((fp < 0.0) != (ff < 0.0)) {
          double lo = ap, hi = aa, flo = fp;
          for (int it = 0; it < 100; ++it) {
            const double m = 0.5*(lo+hi), fm = fv(m);
            if ((flo < 0.0) == (fm < 0.0)) { lo = m; flo = fm; } else hi = m;
          }
          const double r = 0.5*(lo+hi);
          if (best < -900 || std::fabs(r) < std::fabs(best)) best = r;
        }
        ap = aa; fp = ff;
      }
      return best;
    }
  }

  void probe_search_reset() { g_sh.clear(); g_sm.clear(); g_cs.clear(); g_sh_evt = -1; }

  void probe_search_event(const Event *ev, int event_idx) {
    if (ev == nullptr) return;
    // Seed truth, keyed by TrCandMeta::seed -- the index into the seed vector the
    // builder was given, which Event::currentSeedTracks() is.
    //
    // NOT by TrCandMeta::global_seed. MkBuilder.cc:425 sets that to
    // currentSeed(cm.seed).label(), a seed's LABEL, and a label equals its index
    // only after relabelSeedTracksSequentially() -- which ProcessEventHlt() calls
    // and ProcessEventStd deliberately does not. On the forward driver a seed's
    // label is its file label, i.e. the SIM index of the rec track that used it
    // (WriteMemoryFile's seedSimIdx) and -1 for seeds that produced none, so
    // indexing the seed array with it is silently WRONG rather than empty. This
    // key is correct on both paths and retires the relabel dependency, the same
    // trap val_bkfit_seed_truth() was removed for.
    const TrackVec *cseeds = nullptr;
    try { cseeds = &ev->currentSeedTracks(); } catch (...) { cseeds = nullptr; }
    const size_t n_seed = cseeds ? cseeds->size() : 0;
    std::vector<int> seed_sim(n_seed, -1);
    std::vector<float> seed_gf(n_seed, -1.f);
    std::vector<int> seed_nv(n_seed, -1), seed_nm(n_seed, -1);
    for (size_t i = 0; i < n_seed; ++i) {
      auto si = ev->simInfoForCurrentSeed((int) i);
      if (si.is_set()) {
        seed_sim[i] = si.label;
        seed_gf[i]  = (float) si.good_frac();
        seed_nv[i]  = si.n_valid;
        seed_nm[i]  = si.n_match;
      }
    }
    g_sh_evt = event_idx;
    const int n_ls = (int) ev->trLayerSearches_.size();
    const int n_ku = (int) ev->trKalmanUpdates_.size();
    const int n_cs = (int) ev->trCandStates_.size();
    int n_add = 0;
    for (const TrHitMatch &hm : ev->trHitMatches_) {
      if ((g_probe_search_lite & 2) && (hm.layer < 0 || hm.layer >= 4))
        continue;                                   // pixel barrel only
      ProbeSearchHit v;
      v.event = event_idx;
      v.layer = hm.layer;  v.hit = hm.hit;  v.mc_match = hm.mc_match;
      v.dphi = hm.dphi;  v.dq = hm.dq;  v.hit_q_half_len = hm.hit_q_half_len;
      v.passed_preselect = hm.passed_preselect;
      v.passed_pqueue = hm.passed_pqueue;  v.rank = hm.sub_rank;
      v.t_hermite = hm.t_hermite;  v.d_plane_h3 = hm.d_plane_h3;
      v.search_id = hm.search_id;
      // signed residuals, from the predicted point on the module plane
      if (hm.layer >= 0 && hm.layer < (int) ev->layerHits_.size() && hm.hit >= 0 &&
          hm.hit < (int) ev->layerHits_[hm.layer].size()) {
        const Hit &h = ev->layerHits_[hm.layer][hm.hit];
        const LayerInfo &li_t = Config::TrkInfo[hm.layer];
        const unsigned int sid_t = h.detIDinLayer();
        if ((int) sid_t < li_t.n_modules())
          v.mod_tilt = (float) std::fabs(li_t.module_info(sid_t).zdir[2]);
        // hit's own phi sigma, from its covariance: sigma_phi^2 =
        // (y^2 exx - 2xy exy + x^2 eyy) / r^4 -- same projection project_sigma()
        // uses for the track, applied to the hit.
        {
          const double hx = h.x(), hy = h.y(), hr2 = hx*hx + hy*hy;
          if (hr2 > 0) {
            const double pv = (hy*hy*h.exx() - 2*hx*hy*h.exy() + hx*hx*h.eyy()) / (hr2*hr2);
            v.hit_phi_sigma = (float) std::sqrt(std::max(0.0, pv));
          }
        }
        const auto &kp = hm.kine_on_plane.pos;
        const bool bar = Config::TrkInfo[hm.layer].is_barrel();
        const double pq = bar ? kp[2] : std::hypot(kp[0], kp[1]);
        const double hq = bar ? h.z()  : std::hypot(h.x(), h.y());
        v.dq_s = (float)(pq - hq);
        double dp = std::atan2((double)kp[1], (double)kp[0]) - std::atan2(h.y(), h.x());
        while (dp >  M_PI) dp -= 2.0*M_PI;
        while (dp <= -M_PI) dp += 2.0*M_PI;
        v.dphi_s = (float) dp;
      }

      if (hm.search_id >= 0 && hm.search_id < n_ls) {
        const TrLayerSearch &ls = ev->trLayerSearches_[hm.search_id];
        v.is_barrel = ls.is_barrel;  v.is_outward = ls.is_outward;
        // dphi_track / dq_track are stated as 3 sigma in the struct comment.
        v.sigma_q_trk   = ls.dq_track   / 3.0f;
        v.sigma_phi_trk = ls.dphi_track / 3.0f;
        v.cov_xx = ls.cov_xx;  v.cov_xy = ls.cov_xy;
        v.cov_yy = ls.cov_yy;  v.cov_zz = ls.cov_zz;
        if (ls.state_id >= 0 && ls.state_id < n_cs) {
          const TrCandState &cs = ev->trCandStates_[ls.state_id];
          v.pt = cs.state.pT();  v.eta = cs.state.momEta();  v.step = cs.step;
          v.theta = cs.state.theta();
          if (cs.meta_id >= 0 && cs.meta_id < (int) ev->trCandMetas_.size()) {
            const TrCandMeta &cm = ev->trCandMetas_[cs.meta_id];
            v.seed = cm.seed;  v.global_seed = cm.global_seed;  v.sim = cm.sim;
            if (cm.seed >= 0 && cm.seed < (int) seed_sim.size()) {
              v.sim_label = seed_sim[cm.seed];
              v.seed_good_frac = seed_gf[cm.seed];
              v.seed_n_valid   = seed_nv[cm.seed];
              v.seed_n_match   = seed_nm[cm.seed];
              if (v.sim_label >= 0 && hm.layer >= 0)
                v.n_sim_hits_in_layer = ev->countSimHitsInLayer(v.sim_label, hm.layer);
            }
          }
        }
      }
      if (hm.kalman_id >= 0 && hm.kalman_id < n_ku) {
        const TrKalmanUpdate &ku = ev->trKalmanUpdates_[hm.kalman_id];
        v.had_kalman = true;  v.accepted = ku.accepted;  v.chi2 = ku.chi2;
        v.residual_x = hm.residual_x;  v.residual_y = hm.residual_y;
        v.residual_z = hm.residual_z;
      }
      // --- exact helix-plane crossing from the candidate's own state
      if (!(g_probe_search_lite & 1) &&
          hm.search_id >= 0 && hm.search_id < n_ls && hm.hit >= 0 && hm.layer >= 0) {
        const TrLayerSearch &ls = ev->trLayerSearches_[hm.search_id];
        if (ls.state_id >= 0 && ls.state_id < n_cs) {
          const mkfit::TrackState &st = ev->trCandStates_[ls.state_id].state;
          const auto &lhv = ev->layerHits_[hm.layer];
          if (hm.hit < (int) lhv.size()) {
            const Hit &hit = lhv[hm.hit];
            const LayerInfo &li = Config::TrkInfo[hm.layer];
            const unsigned int sid = hit.detIDinLayer();
            if ((int) sid < li.n_modules()) {
              const auto &mi = li.module_info(sid);
              const double nrm[3] = {mi.zdir[0], mi.zdir[1], mi.zdir[2]};
              const double pnt[3] = {mi.pos[0],  mi.pos[1],  mi.pos[2]};
              const double x0[3] = {st.parameters[0], st.parameters[1], st.parameters[2]};
              const double pt = 1.0 / st.invpT();
              const double p0[3] = {pt*std::cos(st.momPhi()), pt*std::sin(st.momPhi()),
                                    pt / std::tan((double) st.theta())};
              const double inv_k = ((st.charge < 0) ? 0.01 : -0.01) * (double)Const::sol * Config::Bfield;
              const double k = 1.0 / inv_k;
              // the Hermite's own turn angle, from the momentum it reports
              const double ph0 = std::atan2(p0[1], p0[0]);
              const double ph1 = std::atan2((double)hm.kine_on_plane.mom[1],
                                            (double)hm.kine_on_plane.mom[0]);
              double dah = ph1 - ph0;
              while (dah >  M_PI) dah -= 2.0*M_PI;
              while (dah <= -M_PI) dah += 2.0*M_PI;
              v.dalpha_hermite = (float) dah;
              // what sp1 / sp2 should have been
              v.alpha_in  = (float) exact_surface_alpha(x0, p0, k,
                              li.is_barrel() ? li.rin()  : li.zmin(), li.is_barrel(), 1.5);
              v.alpha_out = (float) exact_surface_alpha(x0, p0, k,
                              li.is_barrel() ? li.rout() : li.zmax(), li.is_barrel(), 1.5);
              PlaneRoots R = exact_plane_roots(x0, p0, k, nrm, pnt, 1.5);
              v.n_roots = R.n;
              if (R.n > 0) {
                int i_near = 0, i_small = 0;
                for (int i = 1; i < R.n; ++i) {
                  if (std::fabs(R.a[i] - dah) < std::fabs(R.a[i_near] - dah)) i_near = i;
                  if (std::fabs(R.a[i]) < std::fabs(R.a[i_small])) i_small = i;
                }
                const double ae = R.a[i_near];
                v.dalpha_exact = (float) ae;
                v.wrong_crossing = (i_near != i_small);
                v.alpha_small = (float) R.a[i_small];
                const double sa = std::sin(ae), ca = std::cos(ae);
                const double xe[3] = { x0[0] + k*(p0[0]*sa - p0[1]*(1.0-ca)),
                                       x0[1] + k*(p0[1]*sa + p0[0]*(1.0-ca)),
                                       x0[2] + k*p0[2]*ae };
                v.d_hermite_exact = (float) std::sqrt(
                    std::pow(xe[0]-(double)hm.kine_on_plane.pos[0], 2) +
                    std::pow(xe[1]-(double)hm.kine_on_plane.pos[1], 2) +
                    std::pow(xe[2]-(double)hm.kine_on_plane.pos[2], 2));
                if (std::fabs(ae) > 1e-9) v.rel_ds = (float)((dah - ae) / std::fabs(ae));
              }
            }
          }
        }
      }

      g_sh.push_back(v);
      ++n_add;
    }
    // ---- one row per KALMAN UPDATE: the covariance before and after it
    // (needs MKFIT_TRACE_KALMAN_DEBUG for the two states; otherwise no rows)
#ifdef MKFIT_TRACE_KALMAN_DEBUG
    for (const TrKalmanUpdate &ku : ev->trKalmanUpdates_) {
      if (ku.hit_match_id < 0 || ku.hit_match_id >= (int) ev->trHitMatches_.size()) continue;
      const TrHitMatch &hm = ev->trHitMatches_[ku.hit_match_id];
      ProbeCovStep c;
      c.event = event_idx;  c.layer = hm.layer;  c.mc_match = hm.mc_match;
      c.accepted = ku.accepted;  c.chi2 = ku.chi2;
      c.state_in = ku.state_id_in;  c.state_out = ku.state_id_out;
      if (hm.search_id >= 0 && hm.search_id < n_ls) {
        const TrLayerSearch &ls = ev->trLayerSearches_[hm.search_id];
        c.is_barrel = ls.is_barrel;
        if (ls.state_id >= 0 && ls.state_id < n_cs) {
          const TrCandState &cst = ev->trCandStates_[ls.state_id];
          c.step = cst.step;  c.pt = cst.state.pT();  c.eta = cst.state.momEta();
        }
      }
      project_sigma(ku.propagated_state, c.is_barrel, c.sig_q_prop, c.sig_phi_prop);
      project_sigma(ku.updated_state,    c.is_barrel, c.sig_q_upd,  c.sig_phi_upd);
      g_cs.push_back(c);
    }
#endif

    // ---- one row per LAYER-SEARCH: where was the sim track's own hit?
    // Roll up what the search actually did, keyed by search id.
    struct Roll { int nsc = 0; bool mc_sc = false, mc_pre = false, mc_kal = false, mc_pass = false, mc_kept = false;
                  int npass = 0, npass_wrong = 0, nkept_wrong = 0;
                  float mc_c2 = 1e30f, mc_c2a = 1e30f, best_c2 = 1e30f, wrong_c2 = 1e30f, wrong_c2a = 1e30f; };
    // One roll-up per SUB-LAYER of a search: a paired plan entry searches both
    // sub-layers at once, and each gets its own row below (sub 1 = layer_sec).
    std::vector<Roll> roll(2 * n_ls);
    struct KuBest { float c2 = 1e30f, c2a = 1e30f; bool pass = false, kept = false; };
    std::vector<KuBest> ku_best(ev->trHitMatches_.size());
    for (const TrKalmanUpdate &ku : ev->trKalmanUpdates_) {
      if (ku.hit_match_id < 0 || ku.hit_match_id >= (int) ku_best.size()) continue;
      KuBest &kb = ku_best[ku.hit_match_id];
      const bool pass = ku.chi2_acc > -900.f ? ku.chi2_acc < Config::ItrInfo[0].m_v2p2.hit_chi2_cut : ku.chi2 < 30.f;
      kb.pass |= pass;
      kb.kept |= ku.accepted;
      if (ku.chi2 < kb.c2) { kb.c2 = ku.chi2; kb.c2a = ku.chi2_acc; }
    }
    for (const TrHitMatch &hm : ev->trHitMatches_) {
      if (hm.search_id < 0 || hm.search_id >= n_ls) continue;
      const TrLayerSearch &hls = ev->trLayerSearches_[hm.search_id];
      const int sub = (hls.layer_sec >= 0 && hm.is_sec_layer) ? 1 : 0;
      Roll &r = roll[2 * hm.search_id + sub];
      ++r.nsc;
      // A hit can be evaluated more than once: at depth 0 and as the child of
      // tree nodes, and hm.kalman_id names only the LAST evaluation. Take the
      // best of all of them (ku_best, built above from TrKalmanUpdate::hit_match_id).
      float c2 = -999.f, c2a = -999.f;
      bool pass = false, kept = false;
      const int hm_idx = (int)(&hm - &ev->trHitMatches_[0]);
      if (hm.kalman_id >= 0 && hm.kalman_id < n_ku) {
        const KuBest &kb = ku_best[hm_idx];
        c2 = kb.c2;  c2a = kb.c2a;  kept = kb.kept;  pass = kb.pass;
      }
      if (pass) {
        ++r.npass;
        if (c2 < r.best_c2) r.best_c2 = c2;
        if (!hm.mc_match) { ++r.npass_wrong;  if (kept) ++r.nkept_wrong;
          if (c2 < r.wrong_c2) r.wrong_c2 = c2;
          if (c2a < r.wrong_c2a) r.wrong_c2a = c2a; }
      }
      if (hm.mc_match) { r.mc_sc = true;
        if (hm.passed_preselect) r.mc_pre = true;
        if (hm.kalman_id >= 0) { r.mc_kal = true;
          if (pass) r.mc_pass = true;
          if (kept) r.mc_kept = true;
          if (c2 > -900.f && c2 < r.mc_c2) { r.mc_c2 = c2; r.mc_c2a = c2a; } } }
    }
    for (int isub = 0; isub < 2 * n_ls; ++isub) {
      const int is = isub / 2, sub = isub % 2;
      const TrLayerSearch &ls = ev->trLayerSearches_[is];
      if (ls.layer < 0 || ls.state_id < 0 || ls.state_id >= n_cs) continue;
      if (sub == 1 && ls.layer_sec < 0) continue;
      const int lay = sub ? ls.layer_sec : ls.layer;
      ProbeSearchMiss m;
      m.event = event_idx;  m.search_id = is;  m.layer = lay;
      m.is_barrel = ls.is_barrel;
      m.wsr = ls.wsr;  m.wsr_in_gap = ls.wsr_in_gap;
      const TrCandState &cs = ev->trCandStates_[ls.state_id];
      m.step = cs.step;  m.pt = cs.state.pT();  m.eta = cs.state.momEta();
      m.n_scanned = roll[isub].nsc;  m.mc_scanned = roll[isub].mc_sc;
      m.mc_preselect = roll[isub].mc_pre;  m.mc_kalman = roll[isub].mc_kal;
      if (roll[isub].mc_c2   < 1e29f) m.mc_chi2   = roll[isub].mc_c2;
      if (roll[isub].best_c2 < 1e29f) m.best_chi2 = roll[isub].best_c2;
      if (roll[isub].mc_c2a < 1e29f) m.mc_chi2_acc = roll[isub].mc_c2a;
      m.n_pass = roll[isub].npass;  m.n_pass_wrong = roll[isub].npass_wrong;  m.n_kept_wrong = roll[isub].nkept_wrong;
      if (roll[isub].wrong_c2 < 1e29f) m.best_wrong_chi2 = roll[isub].wrong_c2;
      if (roll[isub].wrong_c2a < 1e29f) m.best_wrong_chi2_acc = roll[isub].wrong_c2a;
      if (cs.meta_id >= 0 && cs.meta_id < (int) ev->trCandMetas_.size()) {
        const TrCandMeta &cm = ev->trCandMetas_[cs.meta_id];
        if (cm.seed >= 0 && cm.seed < (int) seed_sim.size())
          m.sim_label = seed_sim[cm.seed];
      }
      if (m.sim_label < 0) { g_sm.push_back(m); continue; }

      // Reconstruct the window centre in 3-D. q_center() is z (barrel) or r
      // (endcap); the other coordinate is taken at the layer's mid-surface.
      const LayerInfo &li = Config::TrkInfo[lay];
      const double rc = m.is_barrel ? 0.5*(li.rin() + li.rout()) : ls.q_center();
      const double zc = m.is_barrel ? ls.q_center() : 0.5*(li.zmin() + li.zmax());
      const double cx = rc*std::cos(ls.phi_center), cy = rc*std::sin(ls.phi_center);

      const Track &st = ev->simTracks_[m.sim_label];
      double best_here = 1e30, best_any = 1e30;
      for (int i = 0; i < st.nTotalHits(); ++i) {
        const HitOnTrack hot = st.getHitOnTrack(i);
        if (hot.index < 0 || hot.layer < 0) continue;
        if (hot.layer >= (int) ev->layerHits_.size()) continue;
        const auto &lv = ev->layerHits_[hot.layer];
        if (hot.index >= (int) lv.size()) continue;
        const Hit &h = lv[hot.index];
        const double d3 = std::sqrt(std::pow(h.x()-cx,2) + std::pow(h.y()-cy,2)
                                  + std::pow(h.z()-zc,2));
        if (d3 < best_any) { best_any = d3; m.near_layer = hot.layer; m.near_d3d = (float) d3; }
        if (hot.layer != lay) continue;
        ++m.n_sim_in_layer;
        double dphi = std::atan2(h.y(), h.x()) - ls.phi_center;
        while (dphi >  M_PI) dphi -= 2.0*M_PI;
        while (dphi <= -M_PI) dphi += 2.0*M_PI;
        if (std::fabs(dphi) < best_here) {
          best_here = std::fabs(dphi);
          m.has_sim_here = true;
          m.sim_dphi = (float) dphi;
          const double q = m.is_barrel ? h.z() : std::hypot(h.x(), h.y());
          m.sim_dq = (float) (q - ls.q_center());
          m.sim_d3d = (float) d3;
          m.sim_dphi_norm = ls.phi_delta > 0 ? (float)(std::fabs(dphi)/ls.phi_delta) : -999.f;
          const double half = 0.5*(ls.q_max - ls.q_min), mid = 0.5*(ls.q_max + ls.q_min);
          m.sim_dq_norm = half > 0 ? (float)((q - mid)/half) : -999.f;
        }
      }
      if (!m.has_sim_here)            m.verdict = 0;
      else if (m.mc_kalman) {
        if (m.mc_chi2 < -900.f)              m.verdict = 5;   // no usable chi2
        else if (!roll[isub].mc_pass)          m.verdict = 5;   // killed by the acceptance cut
        else if (!roll[isub].mc_kept)          m.verdict = 6;   // passed, on no surviving candidate
        else                                 m.verdict = 7;   // kept
      }
      else if (m.mc_preselect)        m.verdict = 4;
      else if (m.mc_scanned)          m.verdict = 3;
      else if (m.sim_dphi_norm > 1.f || std::fabs(m.sim_dq_norm) > 1.f) m.verdict = 1;
      else                            m.verdict = 2;
      g_sm.push_back(m);
    }

    printf("probe_search_event: ev %d -- %d hit matches, %d searches, %d kalman updates\n",
           event_idx, n_add, n_ls, n_ku);
  }

  void probe_search_write(const char *out_file) {
    TFile f(out_file, "RECREATE");
    TTree *t = new TTree("search", "backward-search hit matches");
    ProbeSearchHit v;  t->Branch("h", &v);
    for (const auto &x : g_sh) { v = x; t->Fill(); }
    TTree *t3 = new TTree("cov", "one row per Kalman update: covariance before/after");
    ProbeCovStep c;  t3->Branch("c", &c);
    for (const auto &x : g_cs) { c = x; t3->Fill(); }

    TTree *t2 = new TTree("miss", "one row per layer-search: where was the sim hit");
    ProbeSearchMiss m;  t2->Branch("m", &m);
    for (const auto &x : g_sm) { m = x; t2->Fill(); }
    t->Write();  t2->Write();  t3->Write();  f.Close();
    printf("probe_search_write: %zu hit matches, %zu searches, %zu kalman updates -> '%s'\n",
           g_sh.size(), g_sm.size(), g_cs.size(), out_file);
  }

  // ==========================================================================
  // chi2 through the backward fit, hit by hit -- the TRACE port of
  // test/an-val-bkfit.C. Reads trBkFitUpdates_ instead of a ValFitHit TTree.
  //
  // Two things come for free that the TTree version had to build by hand:
  // the per-track grouping (records carry state_id_in -> TrCandState.meta_id,
  // so a track IS a meta) and the seed truth (TrCandMeta::global_seed indexes
  // seedTracks_ directly, so no label->truth map and no relabel trap).
  //
  // Step 0 is EXCLUDED throughout: the seed state already sits on the outermost
  // hit's plane, so that update re-uses the hit that produced the state and its
  // chi2 measures nothing. That is also the zero-length propagation flagged as
  // TrBkFitUpdate::degenerate_step.
  namespace {
    struct BTAcc {
      std::vector<double> red;                 // chi2/(2n) per track, all cats
      std::vector<double> red_cat[3], per_hit[3], per_hit_mc[3][2];
      long nfail[3] = {0,0,0}, ntot[3] = {0,0,0}, ntrk[3] = {0,0,0};
      long n_nogf = 0, n_short = 0, n_hits = 0, n_trk = 0, n_ev = 0;
      long n_degen = 0, n_fail_nondegen = 0;
      void clear() { *this = BTAcc(); }
    };
    BTAcc g_bt;
  }

  void probe_bkfit_trace_reset() { g_bt.clear(); }

  void probe_bkfit_trace_event(const Event *ev) {
    if (ev == nullptr) return;
    ++g_bt.n_ev;
    // group records by meta -- i.e. by track
    std::map<int, std::vector<const TrBkFitUpdate*>> by_meta;
    for (const auto &t : ev->trBkFitUpdates_) {
      ++g_bt.n_hits;
      if (t.step == 0) ++g_bt.n_degen;     // zero-length by construction
      else if (t.fail) ++g_bt.n_fail_nondegen;
      if (t.state_id_in < 0 || t.state_id_in >= (int) ev->trCandStates_.size()) continue;
      by_meta[ev->trCandStates_[t.state_id_in].meta_id].push_back(&t);
    }
    for (auto &kv : by_meta) {
      const int meta_id = kv.first;
      if (meta_id < 0 || meta_id >= (int) ev->trCandMetas_.size()) continue;
      auto &v = kv.second;
      ++g_bt.n_trk;
      std::sort(v.begin(), v.end(),
                [](const TrBkFitUpdate *a, const TrBkFitUpdate *b){ return a->step < b->step; });

      double c2 = 0.0; int n = 0;
      for (const auto *x : v)
        if (x->step > 0 && isFinite(x->chi2)) { c2 += x->chi2; ++n; }
      if (n < 2) { ++g_bt.n_short; continue; }
      const double r = c2 / (2.0 * n);
      g_bt.red.push_back(r);

      // Seed truth, straight off the meta. No label->truth map, no relabel.
      const TrCandMeta &cm = ev->trCandMetas_[meta_id];
      if (cm.global_seed < 0 || cm.global_seed >= (int) ev->seedTracks_.size()) {
        ++g_bt.n_nogf; continue;
      }
      auto si = ev->simInfoForTrack(ev->seedTracks_[cm.global_seed]);
      if (!si.is_set()) { ++g_bt.n_nogf; continue; }
      const float gf = si.good_frac();
      const int sim_label = si.label;   // THE namespace for mc_track_id
      const int cat = gf >= 0.9999f ? 0 : (gf >= 0.8f ? 1 : 2);

      g_bt.red_cat[cat].push_back(r);
      ++g_bt.ntrk[cat];
      for (const auto *x : v) {
        if (x->step == 0 || !isFinite(x->chi2)) continue;
        g_bt.per_hit[cat].push_back(x->chi2);
        ++g_bt.ntot[cat];
        if (x->fail) ++g_bt.nfail[cat];
        g_bt.per_hit_mc[cat][(x->mc_track_id == sim_label) ? 1 : 0].push_back(x->chi2);
      }
    }
  }

  void probe_bkfit_trace_report(const char *prefix) {
    FILE *log = fopen(Form("%s.txt", prefix), "w");
    auto P = [&](const char *fmt, ...) {
      va_list a; va_start(a, fmt); vprintf(fmt, a); va_end(a);
      if (log) { va_start(a, fmt); vfprintf(log, fmt, a); va_end(a); }
    };
    auto q = [](std::vector<double> &u, double p) {
      return u.empty() ? 0.0 : u[(size_t)(0.01 * p * (u.size() - 1))];
    };

    P("### %s -- chi2 through the backward fit, FROM THE TRACE GRAPH\n###\n", prefix);
    P("### source     : Event::trBkFitUpdates_, grouped by TrCandState.meta_id.\n");
    P("### definition : per-hit chi2 from kalmanOperationPlaneLocal. Each hit is\n");
    P("###              a 2-D measurement in the module plane, so chi2/hit has\n");
    P("###              expectation 2, median of chi2_2 = 1.386, and the track\n");
    P("###              total ~2*n_hits.\n");
    P("### population : %ld hit records, %ld tracks, %ld events.\n",
      g_bt.n_hits, g_bt.n_trk, g_bt.n_ev);
    P("### step 0     : EXCLUDED -- the seed state already lies on the outermost\n");
    P("###              hit's plane, so that update re-uses the hit that made it.\n");
    P("### step-0 hits : %ld (one per track, the zero-length propagation).\n", g_bt.n_degen);
    P("### prop fails : %ld beyond step 0. NOTE Leonardo's getS carries no\n", g_bt.n_fail_nondegen);
    P("###              bracket, so it DETECTS no failures -- 0 here does not\n");
    P("###              mean the fit did not fail.\n###\n");

    std::sort(g_bt.red.begin(), g_bt.red.end());
    P("tracks with >=2 usable hits: %zu  (%ld too short)\n", g_bt.red.size(), g_bt.n_short);
    P("chi2 / (2 * n_hits), per track:\n");
    P("  p10 %.3g   p25 %.3g   median %.3g   p75 %.3g   p90 %.3g   p99 %.3g   max %.3g\n",
      q(g_bt.red,10), q(g_bt.red,25), q(g_bt.red,50), q(g_bt.red,75),
      q(g_bt.red,90), q(g_bt.red,99), g_bt.red.empty() ? 0.0 : g_bt.red.back());
    P("  expectation ~1.0 (median of chi2_k/k is 1 - 2/(3k), so ~0.96 at 15 dof)\n\n");

    // Categories by the SEED's truth purity, not by chi2 itself: splitting on
    // chi2 would slice the symptom by itself, purity tests a cause.
    const char *cname[3] = {"pure  (gf = 1.0)", "mostly (gf >= 0.8)", "dirty (gf < 0.8)"};
    P("tracks with no truth match at all: %ld\n\n", g_bt.n_nogf);
    P("PER TRACK, chi2/(2*n_hits) -- expectation 1.0\n");
    P("  %-20s %7s %9s %9s %9s %9s %9s\n", "category", "tracks", "p25", "median", "p75", "p90", "max");
    for (int c = 0; c < 3; ++c) {
      auto &u = g_bt.red_cat[c]; std::sort(u.begin(), u.end());
      P("  %-20s %7ld %9.3g %9.3g %9.3g %9.3g %9.3g\n", cname[c], g_bt.ntrk[c],
        q(u,25), q(u,50), q(u,75), q(u,90), u.empty() ? 0.0 : u.back());
    }
    P("\nPER HIT, chi2 -- expectation 2, median of chi2_2 = 1.386\n");
    P("  %-20s %7s %9s %9s %9s %9s %9s\n", "category", "hits", "p25", "median", "p75", "p90", "max");
    for (int c = 0; c < 3; ++c) {
      auto &u = g_bt.per_hit[c]; std::sort(u.begin(), u.end());
      P("  %-20s %7zu %9.3g %9.3g %9.3g %9.3g %9.3g\n", cname[c], u.size(),
        q(u,25), q(u,50), q(u,75), q(u,90), u.empty() ? 0.0 : u.back());
    }
    P("\nPER HIT, split by whether THAT HIT's mcTrackID matches the track's sim label\n");
    P("  %-20s %8s %9s %9s | %8s %9s %9s\n", "category",
      "n match", "med", "p90", "n other", "med", "p90");
    for (int c = 0; c < 3; ++c) {
      auto &m = g_bt.per_hit_mc[c][1]; std::sort(m.begin(), m.end());
      auto &o = g_bt.per_hit_mc[c][0]; std::sort(o.begin(), o.end());
      P("  %-20s %8zu %9.3g %9.3g | %8zu %9.3g %9.3g\n", cname[c],
        m.size(), q(m,50), q(m,90), o.size(), q(o,50), q(o,90));
    }
    P("\n  propagation failures: pure %ld/%ld, mostly %ld/%ld, dirty %ld/%ld\n",
      g_bt.nfail[0], g_bt.ntot[0], g_bt.nfail[1], g_bt.ntot[1], g_bt.nfail[2], g_bt.ntot[2]);
    if (log) fclose(log);
  }

  // ==========================================================================
  // THE EXACT GAIN TEST -- information gain, no jacobian reimplementation.
  //
  // For a true Kalman update with measurement matrix H and measurement
  // covariance R,
  //        P'^-1 - P^-1  =  H^T R^-1 H
  // exactly. The right-hand side is PSD and has RANK EQUAL TO THE NUMBER OF
  // MEASURED COORDINATES -- 2 here, the module plane. So from the two stored
  // covariances alone:
  //   * rank > 2  -> the update is injecting information it has no right to;
  //   * a negative eigenvalue -> it is REMOVING information, which no update can;
  //   * the two positive eigenvalues are 1/sigma^2 of the effective measurement,
  //     so they say directly whether the hit error used is too small.
  //
  // This needs no model of jacCurv2Loc and makes no 1-D approximation, which is
  // what the sigma-ratio version could not avoid.
  // ==========================================================================
  namespace {
    struct KIAcc {
      std::vector<double> l1, l2, l3, neg;   // top 3 eigenvalues, most negative
      std::vector<double> cond;              // condition number of whitened P
      std::vector<int> layer;
      long n = 0, n_singular = 0;
    };
    KIAcc g_ki;
  }

  void probe_kinfo_reset() { g_ki = KIAcc(); }

  void probe_kinfo_event(const Event *ev) {
    if (ev == nullptr) return;
#ifdef MKFIT_TRACE_KALMAN_DEBUG
    const int n_hm = (int) ev->trHitMatches_.size();
    for (const TrKalmanUpdate &ku : ev->trKalmanUpdates_) {
      if (ku.hit_match_id < 0 || ku.hit_match_id >= n_hm) continue;
      const TrHitMatch &hm = ev->trHitMatches_[ku.hit_match_id];
      if (!hm.mc_match || hm.layer < 0 || hm.layer >= 4) continue;   // pixel barrel
      // P is SINGULAR by construction -- the transported covariance is rank 5
      // with a null direction along the momentum (S11). So P^-1 does not exist
      // and the information form P'^-1 - P^-1 is ill-posed; measured, it gave a
      // negative eigenvalue in 54 % of updates, which is that null direction and
      // not a defect of the update. Use the form that needs NO inverse of P:
      //      P - P' = P H^T (H P H^T + R)^-1 H P
      // which is PSD of rank <= 2 for a 2-D measurement. Same content, no H, no
      // P^-1, and well conditioned. Whitened by P's own diagonal so the
      // eigenvalues are dimensionless -- a congruence, so rank and inertia,
      // which are the only things read here, are unchanged.
      double sc[6];
      bool okp = true;
      for (int i = 0; i < 6; ++i) {
        const double d = ku.propagated_state.errors.At(i,i);
        if (!(d > 0)) { okp = false; break; }
        sc[i] = 1.0 / std::sqrt(d);
      }
      if (!okp) { ++g_ki.n_singular; continue; }
      TMatrixDSym D(6);
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j <= i; ++j) {
          const double v = ((double) ku.propagated_state.errors.At(i,j) -
                            (double) ku.updated_state.errors.At(i,j)) * sc[i] * sc[j];
          D(i,j) = v; D(j,i) = v;
        }
      g_ki.cond.push_back(1.0);   // not used in this form
      TMatrixDSymEigen eig(D);
      TVectorD ev6 = eig.GetEigenValues();              // descending
      std::vector<double> e(6);
      for (int i = 0; i < 6; ++i) e[i] = ev6[i];
      std::sort(e.begin(), e.end(), std::greater<double>());
      g_ki.l1.push_back(e[0]); g_ki.l2.push_back(e[1]); g_ki.l3.push_back(e[2]);
      g_ki.neg.push_back(e[5]);
      g_ki.layer.push_back(hm.layer);
      ++g_ki.n;
    }
#else
    printf("probe_kinfo_event: needs MKFIT_TRACE_KALMAN_DEBUG\n");
#endif
  }

  void probe_kinfo_report() {
    auto q = [](std::vector<double> u, double p) {
      if (u.empty()) return 0.0;
      std::sort(u.begin(), u.end());
      return u[(size_t)(p/100.*(u.size()-1))]; };
    printf("\n===== EXACT GAIN TEST: P - P' must be PSD of RANK <= 2 =====\n");
    printf("Pixel barrel, MC-matched. n = %ld (singular %ld)\n\n", g_ki.n, g_ki.n_singular);
    if (g_ki.n < 40) { printf("  too few\n"); return; }
    printf("  eigenvalues of the information gain (descending), quantiles:\n");
    printf("    %-14s %12s %12s %12s\n", "", "p25", "median", "p75");
    printf("    %-14s %12.4g %12.4g %12.4g\n", "lambda_1", q(g_ki.l1,25), q(g_ki.l1,50), q(g_ki.l1,75));
    printf("    %-14s %12.4g %12.4g %12.4g\n", "lambda_2", q(g_ki.l2,25), q(g_ki.l2,50), q(g_ki.l2,75));
    printf("    %-14s %12.4g %12.4g %12.4g   <- must be ~0 (rank<=2)\n", "lambda_3", q(g_ki.l3,25), q(g_ki.l3,50), q(g_ki.l3,75));
    printf("    %-14s %12.4g %12.4g %12.4g   <- must be >= 0 (PSD)\n", "lambda_min", q(g_ki.neg,25), q(g_ki.neg,50), q(g_ki.neg,75));
    // rank: how often is lambda_3 non-negligible against lambda_2
    long bad = 0, negs = 0;
    for (size_t i = 0; i < g_ki.l1.size(); ++i) {
      if (g_ki.l2[i] > 0 && g_ki.l3[i] > 1e-2 * g_ki.l2[i]) ++bad;
      if (g_ki.neg[i] < -1e-2 * std::fabs(g_ki.l1[i])) ++negs;
    }
    printf("\n    rank > 2 (lambda_3 > 1e-2 lambda_2): %ld of %ld = %.2f %%\n",
           bad, g_ki.n, 100.*bad/g_ki.n);
    printf("    a NEGATIVE eigenvalue:               %ld of %ld = %.2f %%\n",
           negs, g_ki.n, 100.*negs/g_ki.n);
    printf("\n  THRESHOLDS: lambda_3 and lambda_min sit at ~1e-3 against lambda_1 ~ 2.4,\n");
    printf("  a relative 6e-4 -- float32 round-off through the whitening and the\n");
    printf("  eigen-decomposition. They are noise, not rank or indefiniteness, so the\n");
    printf("  flags above are cut at 1e-2 relative. A real rank-3 component or a real\n");
    printf("  negative direction would be orders of magnitude larger.\n");
    printf("\n  lambda_1, lambda_2 are the variance REMOVED by the hit in the two\n");
    printf("  measured directions, in units of P's own sigma. They are NOT bounded by\n");
    printf("  1: whitened, P - P' is bounded by the CORRELATION matrix of P, whose\n");
    printf("  eigenvalues reach 6 when the parameters are correlated. So lambda_1 ~ 2.4\n");
    printf("  is allowed and says the hit constrains a strongly correlated combination.\n");
    printf("\n  WHAT THIS TESTS, AND WHAT IT DOES NOT: it proves the update is a valid\n");
    printf("  rank-2 Kalman update -- the gain injects no spurious information and\n");
    printf("  removes none outside the measured plane. It says NOTHING about whether R\n");
    printf("  (the hit covariance in the module frame) has the right MAGNITUDE; that\n");
    printf("  needs H, i.e. the jacCurv2Loc chain.\n");
    printf("\n");
  }

  //============================================================================
  // probe_rank -- the worst end-of-layer selection rank along the ancestry of each
  // seed's final forward-search candidate (TrCandMeta::fwd_rank_*, MKFIT_TRACE,
  // MkFinderV2p2). The layer-step score is additive, so a candidate whose
  // ancestors all ranked <= R at beam width W > R also survives at width R: its
  // competitors at width R are a subset of those at W. So in a run at width 10,
  // the final tracks with fwd_rank_max > 6 are the ones width 6 cannot make.
  // Track-to-meta by TrCandMeta::cand, the export index: run WITHOUT
  // --remove-dup, which erases tracks and shifts it. Truth as probe_eff, by the
  // association mode in force (probe_assoc_mtv).
  //============================================================================

  namespace {
    struct RkCfg {
      std::string name;
      // [class][rank bin]: class 0 good, 1 fake, 2 other; rank bin 0: <= rank_over,
      // 1..4: rank_over + 1 .. + 4 and above.
      long n[3][5] = {};
      // first layer above rank_over, by group, per class
      long over_grp[3][9] = {};
      // per |eta| region of the reco track: tracks, of them above rank_over
      long reg[3][3][2] = {};
      long n_nometa = 0;
      // per event (in call order), seed index -> (class, hash of the hit list,
      // worst rank), for the per-seed comparison against the first config.
      struct Fin { int cls; size_t hash; int rank; unsigned long long fwd_hash; };
      std::vector<std::map<int, Fin>> per_ev;
    };
    std::vector<RkCfg> g_rk;

    int rk_group(int l) {
      // 0 PixB, 1 TBPS 1st, 2 TBPS 2nd, 3 TB2S 1st, 4 TB2S 2nd, 5 FPix, 6 TEDD 1st, 7 TEDD 2nd
      if (l < 0) return 8;
      if (l >= 38 && l <= 59) l -= 22;
      if (l <= 3) return 0;
      if (l <= 9) return 1 + (l & 1);
      if (l <= 15) return 3 + (l & 1);
      if (l <= 27) return 5;
      return 6 + (l & 1);
    }
  }

  void probe_rank_reset() { g_rk.clear(); }

  void probe_rank_event(const Event *ev, const char *cfg) {
    if (ev == nullptr) return;
    RkCfg *C = nullptr;
    for (auto &c : g_rk) if (c.name == cfg) C = &c;
    if (!C) { g_rk.push_back(RkCfg()); g_rk.back().name = cfg; C = &g_rk.back(); }
    const int R = Config::V2p2::Diag::rank_over;
    C->per_ev.emplace_back();
    auto &fin = C->per_ev.back();
    std::vector<int> meta_of(ev->candidateTracks_.size(), -1);
    for (int m = 0; m < (int) ev->trCandMetas_.size(); ++m) {
      const int ci = ev->trCandMetas_[m].cand;
      if (ci >= 0 && ci < (int) meta_of.size()) meta_of[ci] = m;
    }
    for (int ci = 0; ci < (int) ev->candidateTracks_.size(); ++ci) {
      if (meta_of[ci] < 0) { ++C->n_nometa; continue; }
      const TrCandMeta &cm = ev->trCandMetas_[meta_of[ci]];
      const Track &c = ev->candidateTracks_[ci];
      TrackExtra extra(c.label());
      extra.setMCTrackIDInfo(c, ev->layerHits_, ev->simHitsInfo_, ev->simTracks_, false, false);
      const int mc = extra.mcTrackID();
      const bool fake = Config::mtvLikeValidation ? (mc == -1 || (mc <= -5 && mc >= -9))
                                                  : (mc < 0 || mc >= (int) ev->simTracks_.size());
      const int cls = fake ? 1 : (mc >= 0 && mc < (int) ev->simTracks_.size()) ? 0 : 2;
      const int rb = cm.fwd_rank_max <= R ? 0 : std::min(4, cm.fwd_rank_max - R);
      ++C->n[cls][rb];
      if (rb > 0) ++C->over_grp[cls][rk_group(cm.fwd_rank_over_layer)];
      const float ae = std::abs(c.momEta());
      const int reg = ae < 0.9f ? 0 : ae < 1.7f ? 1 : 2;
      ++C->reg[cls][reg][0];
      if (rb > 0) ++C->reg[cls][reg][1];
      fin[cm.seed] = {cls, hits_hash(c), cm.fwd_rank_max, cm.fwd_hits_hash};
    }
  }

  void probe_rank_report() {
    const int R = Config::V2p2::Diag::rank_over;
    const char *cls_name[3] = {"good", "fake", "other"};
    const char *grp_name[9] = {"PixB", "TBPS1", "TBPS2", "TB2S1", "TB2S2", "FPix", "TEDD1", "TEDD2", "none"};
    printf("\n--- probe_rank: worst end-of-layer rank along the final candidate's ancestry (rank_over %d) ---\n", R);
    for (auto &C : g_rk) {
      printf("config %s (tracks without a trace meta: %ld)\n", C.name.c_str(), C.n_nometa);
      printf("  %-6s %9s %9s   rank %d / %d / %d / >=%d\n", "class", "tracks", "<=", R + 1, R + 2, R + 3, R + 4);
      for (int k = 0; k < 3; ++k) {
        long t = 0;
        for (int b = 0; b < 5; ++b) t += C.n[k][b];
        printf("  %-6s %9ld %9ld   %6ld %6ld %6ld %6ld\n", cls_name[k], t, C.n[k][0], C.n[k][1], C.n[k][2],
               C.n[k][3], C.n[k][4]);
      }
      printf("  first layer above rank %d, by group:\n  %-6s", R, "class");
      for (int g = 0; g < 9; ++g) printf(" %6s", grp_name[g]);
      printf("\n");
      for (int k = 0; k < 3; ++k) {
        printf("  %-6s", cls_name[k]);
        for (int g = 0; g < 9; ++g) printf(" %6ld", C.over_grp[k][g]);
        printf("\n");
      }
      printf("  by reco |eta|, tracks / above rank %d:   barrel        transition    endcap\n", R);
      for (int k = 0; k < 3; ++k)
        printf("  %-6s %32s %6ld / %-5ld %6ld / %-5ld %6ld / %-5ld\n", cls_name[k], "", C.reg[k][0][0],
               C.reg[k][0][1], C.reg[k][1][0], C.reg[k][1][1], C.reg[k][2][0], C.reg[k][2][1]);
    }
    // Per seed against the first configuration: the same seeds in every config.
    if (g_rk.size() < 2) return;
    const RkCfg &A = g_rk[0];
    for (size_t ic = 1; ic < g_rk.size(); ++ic) {
      const RkCfg &B = g_rk[ic];
      // [class A][class B][same hits][B rank above]; tf: the same, split by the
      // forward search's final candidate instead of the final track.
      long t[3][3][2][2] = {}, tf[3][3][2][2] = {};
      long only_a = 0, only_b = 0;
      const size_t ne = std::min(A.per_ev.size(), B.per_ev.size());
      for (size_t e = 0; e < ne; ++e) {
        for (auto &[s, fa] : A.per_ev[e]) {
          auto it = B.per_ev[e].find(s);
          if (it == B.per_ev[e].end()) { ++only_a; continue; }
          const auto &fb = it->second;
          ++t[fa.cls][fb.cls][fa.hash == fb.hash][fb.rank > R];
          ++tf[fa.cls][fb.cls][fa.fwd_hash == fb.fwd_hash][fb.rank > R];
        }
        for (auto &[s, fb] : B.per_ev[e])
          if (A.per_ev[e].find(s) == A.per_ev[e].end()) ++only_b;
      }
      printf("per seed, %s -> %s (%zu events; seeds only in %s %ld, only in %s %ld):\n", A.name.c_str(),
             B.name.c_str(), ne, A.name.c_str(), only_a, B.name.c_str(), only_b);
      printf("  %-14s %10s %10s %12s %12s\n", "", "same hits", "", "other hits", "");
      printf("  %-14s %10s %10s %12s %12s\n", "", "rank<=", "rank>", "rank<=", "rank>");
      for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
          printf("  %-5s -> %-5s  %10ld %10ld %12ld %12ld\n", cls_name[a], cls_name[b], t[a][b][1][0],
                 t[a][b][1][1], t[a][b][0][0], t[a][b][0][1]);
      printf("  split by the FORWARD search's final candidate (same / other forward hits):\n");
      for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
          printf("  %-5s -> %-5s  %10ld %10ld %12ld %12ld\n", cls_name[a], cls_name[b], tf[a][b][1][0],
                 tf[a][b][1][1], tf[a][b][0][0], tf[a][b][0][1]);
    }
  }

}  // namespace mkfit

#endif  // MKFIT_TRACE
