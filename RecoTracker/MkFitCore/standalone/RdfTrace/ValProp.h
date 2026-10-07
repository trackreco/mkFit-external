#ifndef RecoTracker_MkFitCore_standalone_RdfTrace_ValProp_h
#define RecoTracker_MkFitCore_standalone_RdfTrace_ValProp_h

// ===========================================================================
// Synthetic validation of propagation, Kalman operations and fitting.
//
// GENERATION lives here and needs the mkFit build (the propagator). ANALYSIS
// does not: the output is a plain TTree of dictionary-ed Val* structs, so a
// bare root.exe session can read it after
//     gSystem->Load("libMkFitRootDataFormats.so")
// and run RDataFrame over it with no mkFit at all.
//
// No detector geometry is used. Closure is a property of the PROPAGATOR, not
// of the detector, so the planes here are constructed analytically at chosen
// radii. That also removes a whole class of contamination -- a synthetic plane
// cannot be geometry-incompatible with a sample. Real modules earn their way
// back in at a higher rung of the ladder.
//
// Declaration-only on purpose: this header is parsed by cling through
// test/val-prop.C, so it must pull in no Matriplex header.
//
//   cd /foo/matevz/mic-dev/current/src/standalone
//   unset DISPLAY; export LD_LIBRARY_PATH=.
//   echo .q | ./mkFit --geom CMS-phase2 --input-file <any .bin> --num-events 1
//     --num-thr 1 --shell
//     --shell-command 'gROOT->SetBatch(kTRUE)'
//     --shell-command 'gROOT->ProcessLine(".L ../RecoTracker/MkFitCore/standalone/test/val-prop.C")'
//     --shell-command 'val_gen("val-closure.root")'
// ===========================================================================

namespace mkfit {

  class EventOfHits;

  class Event;


  // Generates a jittered stratified sample of truth helices and, for each of
  // the built-in plane configurations, runs the three-leg round trip
  //     O -> A (conditioning) -> B (the step under test) -> A'
  // writing one ValClosure per trial into TTree "closure", plus the
  // configuration legend into TTree "cfgs".
  //
  // The sampling is a grid with per-cell jitter: every cell is populated, so
  // per-bin statistics stay exact, but no two tracks share a pT, an eta or a
  // phi -- which is what keeps the scan from aliasing against plane radii or
  // against the float grid itself.
  //
  // pT is flat in log10 over [0.05, 200] GeV, eta flat over [-4, +4]. The range
  // reaches below the physics band on purpose: the pixel-scale plane pairs can
  // only approach tangency there. Cut at pT >= 0.3 for physics statements.
  // d0 / z0 displace the production vertex (0 = from the origin).
  void val_gen_closure(const char *out_file = "val-closure.root",
                       int n_pt = 18,
                       int n_eta = 32,
                       int n_phi = 8,
                       unsigned seed = 20260911u,
                       float d0 = 0.0f,
                       float z0 = 0.0f);

  // What turn angles does a REAL track present between consecutive hits?
  //
  // Accumulates across events; call val_dalpha_reset() first, then
  // val_dalpha_add_event(ev) per event, then val_dalpha_write(file).
  // `which`: 0 = cmsswTracks_, 1 = simTracks_, 2 = candidateTracks_.
  void val_dalpha_reset();
  void val_dalpha_add_event(const Event *ev, int which = 0, float pt_min = 0.5f);
  void val_dalpha_write(const char *out_file = "val-dalpha.root");

  // Is a sample's geometry the one this build has?
  //
  // Hit::detIDinLayer() is a SHORT ID assigned when the ntuple was written,
  // by tkinfo[ilay].short_id(detId). It is only meaningful against the
  // geometry the file was written with. A hit lies on its own module's plane
  // by construction, so |n.(hit - plPnt)| must be ~0; anything else is proof
  // that the module lookup is wrong, and then every module-frame quantity on
  // that sample -- residuals, the Kalman projection, chi2 -- is meaningless.
  //
  // Reference values measured with this test: 2.9e-9 cm median on a compatible
  // sample, 0.26 cm median (max 220 cm) on an incompatible one.
  void val_geom_check_reset();
  void val_geom_check_event(const Event *ev);
  void val_geom_check_report(const char *sample_name = "");

  // What is actually in a sample: collection sizes, seed algorithms, and the
  // kinematic reach of each. Decides what can be asked of it before anything
  // is measured on it.
  void val_sample_info(const Event *ev);

  // Where are a sim track's parameters given? At the production vertex, or at
  // its first hit? Everything built on the truth helix depends on the answer.
  void val_simtrack_anchor(const Event *ev, int n_print = 5);

  // Seed diagnostics, before any fit: is the input covariance a fit result or a
  // canned matrix, is it positive definite, and where is the state anchored?
  // `kind`: 0 = pT5 (pix->strip), 1 = T5 (strip->strip), 2 = pix, -1 = all.
  void val_seed_cov(const Event *ev, int kind = 1);

  void val_cov_transport(const char *out_file, int n_samp = 10000, unsigned seed = 20260913);

  void val_max_cands(int n);
  void val_search_material(bool on);
  void val_mkbins_surface_q(bool on);
  void val_extra_dq(float f);
  void val_dq(float trk_fac, float hit_fac, int extra_bins);
  void val_phi_per_hit(bool on, float fac);
  void val_hit_extents(const EventOfHits *eoh, int lay_beg, int lay_end);
  void val_dphi(float trk_fac, float hit_rad, int extra_bins);
  void val_precut(bool q, bool phi);
  void val_hit_chi2(float cut, float trk_fac);
  void val_assoc_mtv(bool on);
  void val_hitorder_reset();
  void val_hitorder_event(const Event *ev, const char *cfg);
  void val_hitorder_report();
  void val_track_scorer(const char *name);
  void val_track_scorer_task(const char *task, const char *name);
  void val_dupl_drth(float central, float obarrel, float forward);
  void val_bkw_pickups(int r0, int r1, int r2, int r3, int r4);
  void val_qprecut_reset();
  void val_qprecut_event(const Event *ev);
  void val_qprecut_report(const char *tag);
  void val_surf_q_hit(bool on);
  // MkFinderV2p2 in-layer combinatorial search and its score, so a scan costs one
  // build and every configuration sees identical events and seeds in one process.
  void val_max_hits_per_layer(int n);
  void val_layer_policy(bool wsr, bool hole_limits, bool stop_cuts);
  void val_late_max_cands(int max_cands, int hits, float step_cm = 0.0f);
  void val_sister_hole(bool on, float chi2 = 0.0f);
  void val_max_cands_fb(int n_fwd, int n_bkw);
  void val_bkfit_outliers(float chi2, int max_outliers, float min_pt);
  void val_dup_drth(float central, float obarrel, float forward);
  void val_reserve_hole_slot(bool on);
  void val_score_mode(int mode, float hit_eff);
  void val_eff_signal_only(bool on);
  void val_final_beam(bool on);
  void val_final_beam_report();
  void val_dupclean(bool on);
  void val_dupclean_event(const Event *ev);
  void val_dupclean_report();
  void val_hit_eff(float fwd, float bkw);
  void val_mini_prop_param_b(int on);
  void val_hit_eff_grp(float pix_brl, float pix_ec, float strip_brl_ps, float strip_brl_2s, float strip_ec);
  void val_final_pick_llh(bool on, float hole_penalty, float tail_penalty);
  void val_score_terms(bool use_rho, float rho_const, bool use_detv, float detv_const);
  void val_score_term_stats(bool on);
  void val_score(float hit_bonus, float chi2_weight, float miss_fwd, float miss_bkw);
  void val_search_lite(int mode);
  void val_cluster_sizes(const Event *ev);
  void val_material_profile(float z);
  void val_kgain_reset();
  void val_kgain_event(const Event *ev);
  void val_kgain_report();
  void val_kinfo_reset();
  void val_kinfo_event(const Event *ev);
  void val_kinfo_report();
  void val_chop_recovery_event(const Event *ev);
  void val_chop_recovery_report(const char *tag);
  void val_track_eff_reset();
  void val_track_eff_event(const Event *ev);
  void val_track_eff_report(const char *tag);

  void val_search_reset();
  void val_search_event(const Event *ev, int event_idx);
  void val_search_write(const char *out_file);

  void val_bkfit_trace_reset();
  void val_bkfit_trace_event(const Event *ev);
  void val_bkfit_trace_report(const char *prefix);
  void val_bkfit_material(bool on);

  // Are cmsswTracks_ independent of the seeds in the same file, or are they
  // built FROM them? If each track contains its seed's hits, the collection is
  // downstream of the seeds and cannot serve as an independent reference.
  void val_seed_vs_cmssw(const Event *ev);

  // How often a hit has a sister in the partner sub-layer of the same physical
  // layer, and how far away. Decides whether in-layer processing may anchor on
  // one sensor. See ValProp.cc.
  // Per-sim-track efficiency, resolved in |eta|, pT and hits-per-layer, paired
  // across configurations. val_eff_ref() names the configuration everything else
  // is differenced against; without it the first one accumulated is used.
  void val_hitpur_reset();
  void val_hitpur_event(const Event *ev, const char *cfg);
  void val_hitpur_report();
  void val_rank_reset();
  void val_rank_event(const Event *ev, const char *cfg);
  void val_rank_report();
  void val_trkdump_reset();
  void val_trkdump_event(const Event *ev, const char *cfg);
  void val_trkdump_write(const char *fn);
  void val_eff_reset();
  void val_eff_ref(const char *cfg);
  void val_eff_ref2(const char *cfg);
  void val_eff_event(const Event *ev, const char *cfg);
  void val_eff_cmssw_event(const Event *ev, const char *cfg);
  void val_eff_report(const char *prefix);

  // pT5 pixel-chop recovery, resolved and paired. Truth-free.
  void val_chopres_reset();
  void val_chopres_ref(const char *cfg);
  void val_chopres_event(const Event *ev, const char *cfg);
  void val_chopres_report(const char *prefix);

  void val_qbins_reset();
  void val_qbins_event(const Event *ev);
  void val_qbins_report();

  void val_sister_reset();
  void val_sister_mcfilter(bool on);
  void val_sister_event(const Event *ev, float pt_min = 0.5f);
  void val_sister_report();

}  // namespace mkfit

#endif
