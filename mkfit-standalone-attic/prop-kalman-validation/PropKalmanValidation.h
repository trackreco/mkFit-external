#ifndef RecoTracker_MkFitCore_standalone_PropKalmanValidation_h
#define RecoTracker_MkFitCore_standalone_PropKalmanValidation_h

// Analytic / synthetic validation of the propagation and the Kalman update.
//
// No MC sample, no event file: every number below is generated from an exact
// uniform-B helix, so the expected answer is known in closed form and there is
// no truth-matching ambiguity. That is the whole point -- these tests say
// whether the machinery is self-consistent, which has to be settled before any
// pull measurement against sim tracks means anything.
//
// The geometry IS needed (layer surfaces, module planes and normals, the
// material grid), so run these from inside mkFit with a geometry loaded:
//
//   cd /foo/matevz/mic-dev/current/src/standalone
//   unset DISPLAY
//   export LD_LIBRARY_PATH=.
//   echo .q | ./mkFit --geom CMS-phase2 --input-file <any .bin> --num-events 1
//        --num-thr 1 --shell
//        --shell-command 'gROOT->SetBatch(kTRUE)'
//        --shell-command 'gROOT->ProcessLine(".L ../RecoTracker/MkFitCore/standalone/test/prop-kalman-validation.C")'
//        --shell-command 'pkv_all()'
//
// (one invocation; the --shell-command options are continuations of the same
// line.)  The whole suite takes about 7 s of compute on top of the ~11 s
// geometry startup.  An input file is needed only because the shell wants one;
// none of these tests reads an event.  The GEOMETRY is needed: layer surfaces,
// module planes and normals, and the material grid all come from the loaded
// TrackerInfo (Config::TrkInfo).
//
// See test/prop-kalman-validation.C for the per-task entry points.

namespace mkfit {

  class Event;

  // Task 0: is the material contribution of the window propagation ("pea") dead?
  void pkv_task0_pea_material();

  // Task 1: propagation closure, A -> B -> A, position and covariance.
  void pkv_task1_closure();

  // Task 2: reported covariance transport vs a numerical jacobian.
  void pkv_task2_numerical_jacobian();

  // Task 3: analytic limits of kalmanOperationPlaneLocal.
  void pkv_task3_kalman_limits();

  // Task 4: synthetic hit sequence through the backward-fit kernel.
  void pkv_task4_sequence();

  // Task 5: timing A/B of the two plane-intersection solvers.
  void pkv_task5_solver_timing();

  // Task 6: pT resolution and propagation-failure rate from a forward+backward
  // fit of REAL SIM TRACKS out of a real event file.  Unlike Tasks 0-5 this one
  // reads an Event, so it is driven per event:
  //
  //   pkv_task6_reset();
  //   for each event:  s.GoToEvent(i);  pkv_task6_add_event(s.event());
  //   pkv_task6_report();
  //
  // pkv_task6_one_event() is the single-event convenience (reset+add+report).
  void pkv_task6_reset();
  void pkv_task6_add_event(const Event *ev);
  void pkv_task6_report();
  void pkv_task6_one_event(const Event *ev);
  void pkv_task6_debug(int n);
  void pkv_task6_set_max_hits_per_layer(int n);
  // Sim-track shaping bitmask: 1 = order by turn angle (path length) instead
  // of by radius, 2 = truncate at the apex (78.5 deg position-momentum
  // angle), 4 = iterative per-hit-chi2 outlier removal.  Default 7.
  void pkv_task6_set_shaping(int mask);
  void pkv_task6_debug_build(int n);
  void pkv_task6_set_outlier(double chi2_cut, int max_removals);
  void pkv_task6_set_order_eps(double eps_rad);
  // Replayable list of the bad fits: '<event> <sim_label> <category>' lines.
  void pkv_task6_write_badcases(const char *path, const char *sample_full_path);

  // ROOT output: run several solver configurations over the same tracks in one
  // process, stashing a snapshot after each, then write one .root holding the
  // per-snapshot histograms and the overlay canvases.
  //   pkv_set_solver(stable_root);   // 0/1 flag, see PropagationMPlex.h
  //   ... reset / add events / report ...
  //   pkv_task6_stash("label");
  //   ... repeat ...
  //   pkv_task6_write_root("pkv-task6-abc.root");
  void pkv_task6_stash(const char *label);
  void pkv_task6_clear_stash();
  void pkv_task6_write_root(const char *path);
  void pkv_set_solver(int stable_root);

  void pkv_run_all();

}  // namespace mkfit

#endif
