// Driver for the synthetic propagation / Kalman validation suite.
//
// Same pattern as minipropagator-unit-test.C: the work lives in the compiled
// library (MkFitCore/standalone/PropKalmanValidation.cc), this file only gives
// cling a declaration to call through, so that no Matriplex header has to be
// parsed by the interpreter.
//
// Load and run, from the standalone build directory:
//
//   cd /foo/matevz/mic-dev/current/src/standalone
//   unset DISPLAY
//   export LD_LIBRARY_PATH=.
//   echo .q | ./mkFit --geom CMS-phase2 \
//     --input-file /foo/matevz/mic-dev/TTbar_noPU_100events_..._LSTSeedsIncludingT5s.bin \
//     --num-events 1 --num-thr 1 --shell \
//     --shell-command 'gROOT->SetBatch(kTRUE)' \
//     --shell-command 'gROOT->ProcessLine(".L ../RecoTracker/MkFitCore/standalone/test/prop-kalman-validation.C")' \
//     --shell-command 'pkv_all()'
//
// Individual tasks: pkv_task0() .. pkv_task4().
//
// An input file is needed only because the shell wants one; none of these tests
// reads an event. The GEOMETRY is needed: layer surfaces, module planes and
// normals, and the material grid all come from the loaded TrackerInfo.

#include "RecoTracker/MkFitCore/standalone/PropKalmanValidation.h"

void pkv_task0() { mkfit::pkv_task0_pea_material(); }
void pkv_task1() { mkfit::pkv_task1_closure(); }
void pkv_task2() { mkfit::pkv_task2_numerical_jacobian(); }
void pkv_task3() { mkfit::pkv_task3_kalman_limits(); }
void pkv_task5() { mkfit::pkv_task5_solver_timing(); }
void pkv_task4() { mkfit::pkv_task4_sequence(); }
void pkv_all()   { mkfit::pkv_run_all(); }

// --- Task 6: needs a real Event.  Driven per event from the shell:
//   --shell-command 'gROOT->ProcessLine(".L .../prop-kalman-validation.C")'
//   --shell-command 'pkv_task6_reset()'
//   --shell-command 's.GoToEvent(1)'   --shell-command 'pkv_task6_ev(s.event())'
//   --shell-command 's.NextEvent()'    --shell-command 'pkv_task6_ev(s.event())'
//   ...
//   --shell-command 'pkv_task6_report()'
void pkv_task6_reset()                 { mkfit::pkv_task6_reset(); }
void pkv_task6_ev(mkfit::Event *ev)    { mkfit::pkv_task6_add_event(ev); }
void pkv_task6_report()                { mkfit::pkv_task6_report(); }
void pkv_task6(mkfit::Event *ev)       { mkfit::pkv_task6_one_event(ev); }
void pkv_task6_debug(int n)            { mkfit::pkv_task6_debug(n); }
void pkv_task6_badcases(const char *path, const char *sample) { mkfit::pkv_task6_write_badcases(path, sample); }
void pkv_task6_maxhpl(int n) { mkfit::pkv_task6_set_max_hits_per_layer(n); }
void pkv_task6_shaping(int m) { mkfit::pkv_task6_set_shaping(m); }
void pkv_task6_dbgbuild(int n) { mkfit::pkv_task6_debug_build(n); }
void pkv_task6_outlier(double c, int n) { mkfit::pkv_task6_set_outlier(c, n); }
void pkv_task6_ordereps(double e) { mkfit::pkv_task6_set_order_eps(e); }
void pkv_task6_stash(const char *label)       { mkfit::pkv_task6_stash(label); }
void pkv_task6_clear_stash()                  { mkfit::pkv_task6_clear_stash(); }
void pkv_task6_root(const char *path)         { mkfit::pkv_task6_write_root(path); }
void pkv_solver(int stable_root) { mkfit::pkv_set_solver(stable_root); }
