#ifndef RecoTracker_MkFitCMS_standalone_Probe_Probe_h
#define RecoTracker_MkFitCMS_standalone_Probe_Probe_h

// ===========================================================================
// Probes: they read an event after it was processed, or a synthetic sample,
// and change nothing. By what they read:
//   closure  Probe-Closure.cc  synthetic tracks with an exact answer, no event
//   checks   Probe-Checks.cc   the input, with a verdict
//   perf     Probe-Perf.cc     the final tracks against truth
//   trace    Probe-Trace.cc    the Tr* records of a trace build (MKFIT_TRACE)
//   taps     Probe-Taps.cc     the readers of the in-core taps
// The knobs and the tap switches are Shell methods, in the section "Knobs" of
// Shell.h.
//
// The records the probes write are Probe* structs (DataFormats/ProbeStructs.h)
// with a dictionary in libMkFitRootDataFormats.so, so a bare root.exe session
// can read them after gSystem->Load("libMkFitRootDataFormats.so").
//
// Declaration-only on purpose: this header is parsed by cling through
// probe.C, so it must pull in no Matriplex header. From a build directory,
// with <SACMS> the directory its Makefile exports as SACMS:
//
//   unset DISPLAY; export LD_LIBRARY_PATH=.
//   echo .q | ./mkFit --geom CMS-phase2 --input-file <any .bin> --num-events 1
//     --num-thr 1 --shell
//     --shell-command 'gROOT->SetBatch(kTRUE)'
//     --shell-command 'gROOT->ProcessLine(".L <SACMS>/Probe/probe.C")'
//     --shell-command 'probe_gen("probe-closure.root")'
// ===========================================================================

namespace mkfit {

  class Event;

  //--------------------------------------------------------------------------
  // Closure: Probe-Closure.cc

  // Generates a jittered stratified sample of truth helices and, for each of
  // the built-in plane configurations, runs the three-leg round trip
  //     O -> A (conditioning) -> B (the step under test) -> A'
  // writing one ProbeClosure per trial into TTree "closure", plus the
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
  void probe_gen_closure(const char *out_file = "probe-closure.root",
                         int n_pt = 18,
                         int n_eta = 32,
                         int n_phi = 8,
                         unsigned seed = 20260911u,
                         float d0 = 0.0f,
                         float z0 = 0.0f);

  // What turn angles does a REAL track present between consecutive hits?
  //
  // Accumulates across events; call probe_dalpha_reset() first, then
  // probe_dalpha_add_event(ev) per event, then probe_dalpha_write(file).
  // `which`: 0 = cmsswTracks_, 1 = simTracks_, 2 = candidateTracks_.
  void probe_dalpha_reset();
  void probe_dalpha_add_event(const Event *ev, int which = 0, float pt_min = 0.5f);
  void probe_dalpha_write(const char *out_file = "probe-dalpha.root");

  void probe_cov_transport(const char *out_file, int n_samp = 10000, unsigned seed = 20260913);

  //--------------------------------------------------------------------------
  // Checks: Probe-Checks.cc

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
  void probe_geom_check_reset();
  void probe_geom_check_event(const Event *ev);
  void probe_geom_check_report(const char *sample_name = "");

  // What is actually in a sample: collection sizes, seed algorithms, and the
  // kinematic reach of each. Decides what can be asked of it before anything
  // is measured on it.
  void probe_sample_info(const Event *ev);

  //--------------------------------------------------------------------------
  // Perf: Probe-Perf.cc

  void probe_assoc_mtv(bool on);

  void probe_track_eff_reset();
  void probe_track_eff_event(const Event *ev);
  void probe_track_eff_report(const char *tag);

  // Per-sim-track efficiency, resolved in |eta|, pT and hits-per-layer, paired
  // across configurations. probe_eff_ref() names the configuration everything else
  // is differenced against; without it the first one accumulated is used.
  void probe_eff_signal_only(bool on);
  void probe_eff_reset();
  void probe_eff_ref(const char *cfg);
  void probe_eff_ref2(const char *cfg);
  void probe_eff_event(const Event *ev, const char *cfg);
  void probe_eff_cmssw_event(const Event *ev, const char *cfg);
  void probe_eff_report(const char *prefix);

  void probe_trkdump_reset();
  void probe_trkdump_event(const Event *ev, const char *cfg);
  void probe_trkdump_write(const char *fn);

  //--------------------------------------------------------------------------
  // Trace: Probe-Trace.cc, built only with MKFIT_TRACE

#ifdef MKFIT_TRACE
  void probe_search_lite(int mode);
  void probe_search_reset();
  void probe_search_event(const Event *ev, int event_idx);
  void probe_search_write(const char *out_file);

  void probe_bkfit_trace_reset();
  void probe_bkfit_trace_event(const Event *ev);
  void probe_bkfit_trace_report(const char *prefix);

  void probe_kinfo_reset();
  void probe_kinfo_event(const Event *ev);
  void probe_kinfo_report();

  void probe_rank_reset();
  void probe_rank_event(const Event *ev, const char *cfg);
  void probe_rank_report();
#endif

  //--------------------------------------------------------------------------
  // Tap readers: Probe-Taps.cc. The taps are switched by Shell::SetTapFinalBeam()
  // and Shell::SetTapScoreStats().

  void probe_final_beam_report();
  void probe_score_term_stats_report();

}  // namespace mkfit

#endif
