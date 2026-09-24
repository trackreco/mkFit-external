#ifndef mkfit_standalone_seedgeom_SeedGeom_h
#define mkfit_standalone_seedgeom_SeedGeom_h

// Pure-geometry pixel seeding study: can a good starting sample be extracted
// from layer-to-layer deltas alone, with no fit, no covariance and no vertex
// as an INPUT?
//
// Two independent parts, both driven per event from the shell.
//
// PART A -- truth vertices from the sample.  The .bin carries no vertex
// collection, but sim-track parameters are at the PRODUCTION POINT and every
// track of one vertex shares the identical position, so an exact group-by on
// the (x,y,z) bit pattern recovers the vertex set.  Tracks also carry
// ProdType {Signal, InTimePU, OutOfTimePU}, set in WriteMemoryFile.cc from
// sim_bunchCrossing / sim_event, so the trigger vertex is labelled.
// Off-beamline groups are secondary production points -- a free displaced
// sample, binned in radius.
//
// PART B -- geometric triplets.  Three space points are 6 measurements against
// 5 helix parameters, so there is EXACTLY ONE constraint, and it is
// longitudinal: z is linear in arc length.  That test involves neither the
// beamline nor D0, so it keeps displaced tracks.  The two windows are:
//
//   phi   w = (rb - ra)/(2 R_min)  +  D0_max (1/ra - 1/rb)
//   q     z_c predicted from a,b linearly in r; the curvature correction
//         s - r = r^3/(24 R^2) is sub-mm inside the pixel volume
//
// so the only physics inputs are two scalars, pT_min and D0_max, and both do
// nothing but widen one window.  What is measured: how many pairs each window
// admits, how much the q prediction reduces them, the purity of what survives,
// and the resolution of the z_v intercept that falls out of every surviving
// triplet.
//
// Run (one invocation; --shell-command options are continuations):
//
//   cd /foo/matevz/mic-dev/current/src/standalone
//   unset DISPLAY; export LD_LIBRARY_PATH=.
//   echo .q | ./mkFit --geom CMS-phase2 --seed-input cmssw
//        --input-file /foo/matevz/mic-dev/trackingNtuple_HLT_2026_March.bin
//        --num-events N --num-thr 1 --shell
//        --shell-command 'gROOT->SetBatch(kTRUE)'
//        --shell-command 'gROOT->ProcessLine(".L ../RecoTracker/MkFitCore/standalone/test/an-seedgeom.C")'
//        --shell-command 'sg_reset()'
//        --shell-command 's.GoToEvent(1)'  --shell-command 'sg_ev(s.event())'
//        ... repeat ...
//        --shell-command 'sg_report()'
//
// The driver script seedgeom-run.sh (beside this file) generates the loop.

namespace mkfit {
  class Event;

  void sg_reset();
  void sg_add_event(Event *ev);
  void sg_diag(Event *ev);
  void sg_rdiag(Event *ev);
  void sg_report();
  void sg_write_root(const char *prefix);

  // Knobs.  Defaults: pT_min 0.9 GeV, D0_max 0.1 cm, layers (0,1,2),
  // q tolerance 0.15 cm, part B on.
  void sg_set_ptmin(float pt_min);
  void sg_set_d0max(float d0_max_cm);
  void sg_set_layers(int la, int lb, int lc);
  void sg_set_qwin(float dz_cm);
  void sg_set_do_partb(bool on);
  void sg_set_verbose(bool on);
  void sg_set_vtx_gap(float dz_cm);
  void sg_set_bl_rcut(float dr_cm);
  void sg_set_layer4(int ld);
  void sg_set_qwin_d(float dz_cm);
  void sg_set_phiwin_d(float dphi_rad);
  void sg_set_phi_margin(float dphi_rad);
  void sg_set_nr(int n_radial_subbins);
  void sg_set_phi_lin(int mode, float margin_rad);
  void sg_set_grid(int nphi, int nz);
  void sg_set_light(bool on);
  void sg_set_vf(float win_cm, int min_votes, float match_cm);
}  // namespace mkfit

#endif
