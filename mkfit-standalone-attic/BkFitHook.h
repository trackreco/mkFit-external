#ifndef RecoTracker_MkFitCore_src_BkFitHook_h
#define RecoTracker_MkFitCore_src_BkFitHook_h

// A null-by-default hook for recording the backward fit hit by hit.
//
// bkFitFitTracksProp2Plane() computes a per-hit chi2 in `tmp_chi2` and then
// only accumulates it, so nothing downstream can see where a track's chi2 came
// from. This lets a validation build observe each step without putting any
// storage or ROOT dependency into the core library: the pointer is null unless
// something sets it, and the call costs one predictable branch per lane.
//
// The residual is given in the MODULE frame -- across the strip, along the
// strip, and along the normal -- never as a 3D distance. In this detector a 3D
// residual folds the precise direction together with one that is not measured
// at all, and reporting the sum has already produced two wrong conclusions.

namespace mkfit {

  struct BkFitHitRec {
    int label = -1;      // track label
    int layer = -1;
    int mc_hit_id = -1;  // Hit::mcHitID(), for joining to truth offline
    int step = -1;       // 0 = outermost hit, counting inward
    int fail = 0;        // propagation fail flag
    float chi2 = -1.f;
    float chi2_cum = -1.f;  // running total BEFORE this hit
    float d_xdir = 0.f;  // propagated - hit, across the strip  [cm]
    float d_ydir = 0.f;  // along the strip (unmeasured by a parallel pair)
    float d_zdir = 0.f;  // along the module normal (should be ~0)
    float pt = 0.f, eta = 0.f;
  };

  // Set to a collector to observe; leave null in production.
  extern void (*g_bkfit_hit_hook)(const BkFitHitRec &);

  // Validation knob: the factor bkFitInputTracks() applies to the input
  // covariance. It is a VARIANCE scale, so 100 means 10x in sigma. Production
  // value is 100; exposed here so the scan can be run without a rebuild.
  extern float g_bkfit_err_scale;

}  // namespace mkfit

#endif
