#ifndef RecoTracker_MkFitCore_standalone_DataFormats_RntConversions_h
#define RecoTracker_MkFitCore_standalone_DataFormats_RntConversions_h

#include "RecoTracker/MkFitCore/standalone/DataFormats/RntStructs.h"
#include "RecoTracker/MkFitCore/standalone/Event.h"

#include "RecoTracker/MkFitCore/interface/HitStructures.h"
#include "RecoTracker/MkFitCore/interface/SteeringParams.h"
#include "RecoTracker/MkFitCore/interface/Track.h"
#include "RecoTracker/MkFitCore/src/Matrix.h"
#include "RecoTracker/MkFitCore/src/MiniPropagators.h"
#include "RecoTracker/MkFitCore/src/MkBins.h"
#include "RecoTracker/MkFitCore/src/MkRZLimits.h"

namespace mkfit {
  namespace miprops = mkfit::mini_propagators;

  inline EVec3 state2pos(const miprops::State &s) { return {s.x, s.y, s.z}; }
  inline EVec3 state2mom(const miprops::State &s) { return {s.px, s.py, s.pz}; }
  inline EBiVec3 state2bivec3(const miprops::State &s) { return {state2pos(s), state2mom(s)}; }

  inline EVec3 statep2pos(const miprops::StatePlex &s, int i) { return {s.x[i], s.y[i], s.z[i]}; }
  inline EVec3 statep2mom(const miprops::StatePlex &s, int i) { return {s.px[i], s.py[i], s.pz[i]}; }
  inline EBiVec3 statep2bivec3(const miprops::StatePlex &s, int i) { return {statep2pos(s, i), statep2mom(s, i)}; }
  inline PropInfo statep2propinfo(const miprops::StatePlex &s, int i) {
    return {statep2bivec3(s, i), s.dalpha[i], s.fail_flag[i]};
  }

  inline EVec3 hit2pos(const Hit &h) { return {h.x(), h.y(), h.z()}; }
  inline EVec3 track2pos(const TrackBase &s) { return {s.x(), s.y(), s.z()}; }
  inline EVec3 track2mom(const TrackBase &s) { return {s.px(), s.py(), s.pz()}; }
  inline EBiVec3 track2bivec3(const TrackBase &s) { return {track2pos(s), track2mom(s)}; }
  inline EBiVec3 state2bivec3(const TrackState &s) { return {{s.x(), s.y(), s.z()}, {s.px(), s.py(), s.pz()}}; }

  // A hash of a track's hits, to tell whether two runs ended a seed on the same hits.
  // FNV-1a (64-bit offset basis and prime) over one key per hit, layer * 1000003 +
  // index; 1000003 is a prime above any hit index in a layer.
  inline unsigned long long hits_hash(const Track &t) {
    constexpr unsigned long long kFnvOffsetBasis = 1469598103934665603ull;
    constexpr unsigned long long kFnvPrime = 1099511628211ull;
    unsigned long long h = kFnvOffsetBasis;
    for (int ih = 0; ih < t.nTotalHits(); ++ih) {
      const HitOnTrack hot = t.getHitOnTrack(ih);
      h = (h ^ (unsigned long long)(hot.layer * 1000003 + hot.index)) * kFnvPrime;
    }
    return h;
  }

  // TrLayerSearch from MkFinderV2p2::determine_search_windows(), lane i.
  inline void fill_layer(TrLayerSearch &ls, const LayerControl &lc, const MkRZLimits &rz) {
    ls.layer = lc.m_layer;
    ls.layer_sec = lc.m_layer_sec;
    ls.is_barrel = rz.m_is_barrel;
    ls.is_outward = rz.m_is_outward;
  }
  inline void fill_window(TrLayerSearch &ls, const MkBins &B, const MkBinTrackCovExtract &TCE, int i) {
    ls.prop_entry = statep2propinfo(B.m_sp1, i);
    ls.prop_exit = statep2propinfo(B.m_sp2, i);
    ls.phi_center = B.m_phi_center[i];
    ls.phi_delta = B.m_phi_delta[i];
    ls.q_min = B.m_q_min[i];
    ls.q_max = B.m_q_max[i];
    ls.dphi_track = B.m_dphi_track[i];
    ls.dq_track = B.m_dq_track[i];
    ls.cov_xx = TCE.m_cov_0_0[i];
    ls.cov_xy = TCE.m_cov_0_1[i];
    ls.cov_yy = TCE.m_cov_1_1[i];
    ls.cov_zz = TCE.m_cov_2_2[i];
  }
  // bl_s is null for a single layer.
  inline void fill_bins(TrLayerSearch &ls, const MkBinLimits &bl_p, const MkBinLimits *bl_s, int i) {
    ls.p1 = bl_p.p1[i];
    ls.p2 = bl_p.p2[i];
    ls.q1 = bl_p.q1[i];
    ls.q2 = bl_p.q2[i];
    if (bl_s) {
      ls.p1_sec = bl_s->p1[i];
      ls.p2_sec = bl_s->p2[i];
      ls.q1_sec = bl_s->q1[i];
      ls.q2_sec = bl_s->q2[i];
    }
  }

  // The pre-selection part of a TrHitMatch, from MkFinderV2p2::preselect_hit_batch():
  // lane h of s is the state on the hit's module plane the cut was taken on, t the
  // Hermite parameter there and d_plane the distance to the plane left after the
  // solve. hit_orig_idx is the hit's index in the layer's hit vector. The residuals
  // are against the hit, in the module frame.
  inline void fill_hitmatch_presel(TrHitMatch &hm,
                                   const miprops::StatePlex &s,
                                   int h,
                                   float t,
                                   float d_plane,
                                   const LayerOfHits &L,
                                   int hit_orig_idx,
                                   float ddphi,
                                   float ddq,
                                   float q_half_len,
                                   bool passed) {
    const Hit &hit = L.refHit(hit_orig_idx);
    const ModuleInfo &mi = L.layer_info().module_info(hit.detIDinLayer());
    hm.kine_on_plane = statep2bivec3(s, h);
    hm.t_hermite = t;
    hm.d_plane_h3 = d_plane;
    hm.dphi = ddphi;
    hm.dq = ddq;
    hm.hit_q_half_len = q_half_len;
    hm.passed_preselect = passed;
    const SVector3 ydir = mi.calc_ydir();
    const EVec3 res = statep2pos(s, h) - hit2pos(hit);
    hm.residual_x = res.Dot(EVec3(mi.xdir[0], mi.xdir[1], mi.xdir[2]));
    hm.residual_y = res.Dot(EVec3(ydir[0], ydir[1], ydir[2]));
    hm.residual_z = res.Dot(EVec3(mi.zdir[0], mi.zdir[1], mi.zdir[2]));
  }

  // TrHitMatch for a pre-selected hit, without the state and search ids. mc_match:
  // the hit was made by the seed's sim track. layer is L.layer_id(), not the plan
  // entry's layer: with the secondary sub-layer scanned the two differ.
  inline TrHitMatch make_hitmatch(const Event *ev, int seed_idx, const LayerOfHits &L, int hit_idx, bool is_sec) {
    const int sim_lbl = ev->simInfoForCurrentSeed(seed_idx).label;
    const int hit_lbl = ev->simHitsInfo_[L.refHit(hit_idx).mcHitID()].mcTrackID();
    TrHitMatch hm{-1, -1, -1, (int)L.layer_id(), hit_idx, sim_lbl == hit_lbl};
    hm.is_sec_layer = is_sec;
    return hm;
  }

  inline SimSeedInfo evsi2ssinfo(const Event *ev, int seed_idx) {
    SimSeedInfo ssi;
    Event::SimInfoFromHits sifh = ev->simInfoForCurrentSeed(seed_idx);
    if (sifh.is_set()) {
      ssi.s_sim = track2bivec3(ev->simTracks_[sifh.label]);
      ssi.sim_lbl = sifh.label;
      ssi.n_hits = sifh.n_hits;
      ssi.n_match = sifh.n_match;
      ssi.has_sim = true;
    }
    auto seed = ev->currentSeed(seed_idx);
    ssi.s_seed = track2bivec3(seed);
    ssi.seed_lbl = seed.label();
    ssi.seed_idx = seed_idx;
    return ssi;
  }
}  // namespace mkfit

#endif
