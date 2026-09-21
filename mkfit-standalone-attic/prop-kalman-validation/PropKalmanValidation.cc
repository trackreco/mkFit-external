#include "RecoTracker/MkFitCore/standalone/PropKalmanValidation.h"
#include "RecoTracker/MkFitCore/standalone/ConfigStandalone.h"

#include "RecoTracker/MkFitCore/interface/Config.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/MkFitCore/interface/cms_common_macros.h"

#include "RecoTracker/MkFitCore/src/Matrix.h"
#include "RecoTracker/MkFitCore/src/PropagationMPlex.h"
#include "RecoTracker/MkFitCore/src/KalmanUtilsMPlex.h"
#include "RecoTracker/MkFitCore/src/MiniPropagators.h"
#include "RecoTracker/MkFitCore/src/MkBins.h"

#include "RecoTracker/MkFitCore/interface/Hit.h"
#include "RecoTracker/MkFitCore/interface/Track.h"
#include "RecoTracker/MkFitCore/standalone/Event.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <array>
#include <map>
#include <random>
#include <string>
#include <vector>

#ifdef WITH_ROOT
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TStyle.h"
#endif

// ============================================================================
// Synthetic validation of propagation and Kalman update.  See the header for
// how to run it.  Everything here is generated from an exact uniform-B helix,
// so the expected answer is known analytically.
//
// Two conventions that everything below depends on, both taken from the code
// under test rather than invented here:
//
//  * the 6-parameter state is CCS, par = (x, y, z, 1/pT, phi_mom, theta), with
//    1/pT kept positive and the charge carried separately (see
//    kalmanCheckChargeFlip());
//  * the helix, exactly, in the turning angle alpha, with k = 1/inv_k and
//    inv_k = -0.01 * sol * B for a positive charge (MiniPropagators.h):
//      x(a) = x0 + k (px sin a - py (1 - cos a))
//      y(a) = y0 + k (py sin a + px (1 - cos a))
//      z(a) = z0 + k pz a
//      p(a) = p rotated by +a about z
//    so dr/da = k p(a) exactly and the signed 3D path length is s = k |p| a.
//    Since k < 0 for a positive charge, FORWARD motion means alpha < 0 for
//    q > 0 -- the trap that the existing minipropagator unit test documents.
//    Parametrizing by s instead of alpha removes the sign bookkeeping, so that
//    is what Helix does.
// ============================================================================

namespace mkfit {

  namespace {

    // --------------------------------------------------------------------
    // Small statistics accumulator.  Nothing fancy: these tests want a
    // median, a 90th percentile and a max per bin, and want to say how many
    // entries were non-finite rather than silently dropping them.

    // std::isfinite() is NOT usable in this build: -Ofast implies
    // -ffinite-math-only, under which the compiler may fold it to true.  That is
    // exactly why mkfit::isFinite() exists (cms_common_macros.h) and does the
    // test on the exponent bits.  Same thing for double here.  Getting this
    // wrong let NaNs into the sorted vector and produced NaN medians on the
    // first attempt -- worth knowing, because any analysis code written against
    // this build has the same trap.
    inline bool fin(double x) {
      unsigned long long u;
      std::memcpy(&u, &x, sizeof(u));
      return ((u >> 52) & 0x7ffull) != 0x7ffull;
    }
    inline bool fin(float x) { return mkfit::isFinite(x); }

    struct Stats {
      std::vector<double> v;
      int n_bad = 0;   // non-finite entries, counted and excluded
      bool sorted = false;

      void add(double x) {
        if (fin(x)) { v.push_back(x); sorted = false; }
        else ++n_bad;
      }
      size_t n() const { return v.size(); }
      void finalize() { if (!sorted) { std::sort(v.begin(), v.end()); sorted = true; } }
      double pct(double p) {
        if (v.empty()) return 0.0;
        finalize();
        double idx = 0.01 * p * (v.size() - 1);
        size_t i = (size_t) idx;
        if (i + 1 >= v.size()) return v.back();
        double f = idx - i;
        return v[i] * (1.0 - f) + v[i + 1] * f;
      }
      double med() { return pct(50.0); }
      double max() { if (v.empty()) return 0.0; finalize(); return v.back(); }
      double mean() {
        if (v.empty()) return 0.0;
        double s = 0; for (double x : v) s += x; return s / v.size();
      }
    };

    // --------------------------------------------------------------------
    // Exact uniform-B helix, parametrized by signed 3D path length s.

    struct Helix {
      double x0 = 0, y0 = 0, z0 = 0;
      double px0 = 0, py0 = 0, pz0 = 0;
      double pt = 1, pmag = 1, theta = 1.57, phi0 = 0;
      int chg = 1;
      double k = 0;      // cm/GeV, signed; k = 1/inv_k
      double bfield = Config::Bfield;

      void init(double a_pt, double eta, int charge, double a_phi0, double bf = Config::Bfield) {
        pt = a_pt; chg = charge; phi0 = a_phi0; bfield = bf;
        theta = 2.0 * std::atan(std::exp(-eta));
        px0 = pt * std::cos(phi0);
        py0 = pt * std::sin(phi0);
        pz0 = pt / std::tan(theta);
        pmag = pt / std::sin(theta);
        const double inv_k = ((charge < 0) ? 0.01 : -0.01) * (double)Const::sol * bf;
        k = 1.0 / inv_k;
        x0 = y0 = z0 = 0.0;
      }

      double alpha_of_s(double s) const { return s / (k * pmag); }

      void at_alpha(double a, double p[3], double m[3]) const {
        const double sa = std::sin(a), ca = std::cos(a);
        p[0] = x0 + k * (px0 * sa - py0 * (1.0 - ca));
        p[1] = y0 + k * (py0 * sa + px0 * (1.0 - ca));
        p[2] = z0 + k * pz0 * a;
        m[0] = px0 * ca - py0 * sa;
        m[1] = px0 * sa + py0 * ca;
        m[2] = pz0;
      }
      void at_s(double s, double p[3], double m[3]) const { at_alpha(alpha_of_s(s), p, m); }

      double r_at_s(double s) const {
        double p[3], m[3]; at_s(s, p, m); return std::hypot(p[0], p[1]);
      }
      double z_at_s(double s) const {
        double p[3], m[3]; at_s(s, p, m); return p[2];
      }

      // Fill a CCS parameter vector at path length s.
      void ccs_at_s(double s, float par[6]) const {
        double p[3], m[3]; at_s(s, p, m);
        par[0] = (float) p[0];
        par[1] = (float) p[1];
        par[2] = (float) p[2];
        par[3] = (float) (1.0 / pt);
        par[4] = (float) std::atan2(m[1], m[0]);
        par[5] = (float) theta;
      }

      // Radius of the trajectory circle and of the reachable band -- used to
      // decide whether a target radius is reachable at all.
      double rc() const { return std::abs(k) * pt; }
      double dc() const {  // distance of the circle centre from the origin
        const double cx = x0 - k * py0, cy = y0 + k * px0;
        return std::hypot(cx, cy);
      }
    };

    // --------------------------------------------------------------------
    // State-space sampling.  Deliberately a grid rather than a random draw:
    // per-pT-bin reporting is asked for, and a grid makes the bins exact.

    struct SamplePoint { double pt, eta, phi0; int chg; };

    // pT log-spaced 0.5 .. 100 GeV.
    const double kPtGrid[] = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 100.0};
    constexpr int kNPt = sizeof(kPtGrid) / sizeof(kPtGrid[0]);
    const double kEtaGrid[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
    constexpr int kNEta = sizeof(kEtaGrid) / sizeof(kEtaGrid[0]);
    // Four phis, deliberately not multiples of pi/2 -- module boundaries and
    // the atan2 branch cut are both places where a bug could hide.
    const double kPhiGrid[] = {0.3, 1.7, 3.0, -2.1};
    constexpr int kNPhi = sizeof(kPhiGrid) / sizeof(kPhiGrid[0]);

    void build_sample(std::vector<SamplePoint> &out, int nphi = kNPhi) {
      out.clear();
      for (int ip = 0; ip < kNPt; ++ip)
        for (int ie = 0; ie < kNEta; ++ie)
          for (int iq = 0; iq < 2; ++iq)
            for (int ih = 0; ih < nphi; ++ih)
              out.push_back({kPtGrid[ip], kEtaGrid[ie], kPhiGrid[ih], iq ? +1 : -1});
    }

    int pt_bin(double pt) {
      for (int i = 0; i < kNPt; ++i) if (pt <= kPtGrid[i] * 1.0001) return i;
      return kNPt - 1;
    }

    // --------------------------------------------------------------------
    // Layer crossings of the truth helix, and the module plane at each.

    struct Crossing {
      int layer = -1;
      int module_sid = -1;
      double s = 0;             // path length of the crossing with the MODULE PLANE
      double pos[3] = {0, 0, 0};
      double mom[3] = {0, 0, 0};
      bool is_barrel = true;
    };

    // Bisect f(s) = target - coordinate(s) on [s1, s2], where the two ends
    // bracket a sign change.
    double bisect_surface(const Helix &h, bool barrel, double target, double s1, double s2) {
      auto f = [&](double s) { return (barrel ? h.r_at_s(s) : h.z_at_s(s)) - target; };
      double f1 = f(s1);
      for (int i = 0; i < 80; ++i) {
        const double sm = 0.5 * (s1 + s2);
        const double fm = f(sm);
        if ((f1 < 0) == (fm < 0)) { s1 = sm; f1 = fm; } else { s2 = sm; }
      }
      return 0.5 * (s1 + s2);
    }

    // Refine a crossing onto a module's infinite plane: solve
    // n . (r(s) - modpos) = 0 by Newton, with dr/ds = phat exactly.
    // Returns false if it fails to converge (grazing incidence).
    bool refine_to_plane(const Helix &h, const SVector3 &mpos, const SVector3 &mnrm,
                         double s_in, double &s_out) {
      double s = s_in;
      for (int it = 0; it < 40; ++it) {
        double p[3], m[3]; h.at_s(s, p, m);
        const double g = mnrm[0] * (p[0] - mpos[0]) + mnrm[1] * (p[1] - mpos[1]) + mnrm[2] * (p[2] - mpos[2]);
        const double dg = (mnrm[0] * m[0] + mnrm[1] * m[1] + mnrm[2] * m[2]) / h.pmag;
        if (std::abs(dg) < 1e-4) return false;       // grazing: ill-conditioned
        const double ds = -g / dg;
        s += ds;
        if (std::abs(ds) < 1e-9) { s_out = s; return true; }
      }
      s_out = s;
      return true;
    }

    // Nearest module of a layer to a point, by 3D distance of the module centre.
    int nearest_module(const LayerInfo &li, const double p[3]) {
      int best = -1;
      double bd2 = 1e30;
      const int nm = li.n_modules();
      for (int i = 0; i < nm; ++i) {
        const ModuleInfo &mi = li.module_info(i);
        const double dx = mi.pos[0] - p[0], dy = mi.pos[1] - p[1], dz = mi.pos[2] - p[2];
        const double d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < bd2) { bd2 = d2; best = i; }
      }
      return best;
    }

    // All layer crossings of the helix, ordered by path length.
    //
    // The sequence is generated from the geometry, not from a layer plan, so it
    // is exactly the set of layers the track really traverses -- which is also
    // why the resulting hit counts (see the printout) come out at the 10-14 of a
    // real track without any tuning.
    void find_crossings(const Helix &h, const TrackerInfo &ti, std::vector<Crossing> &out,
                        double smax = 1400.0, double ds = 0.5) {
      out.clear();
      const int nsamp = (int) (smax / ds) + 1;
      std::vector<float> rs(nsamp), zs(nsamp);
      for (int i = 0; i < nsamp; ++i) {
        double p[3], m[3];
        h.at_s(i * ds, p, m);
        rs[i] = (float) std::hypot(p[0], p[1]);
        zs[i] = (float) p[2];
      }

      for (int l = 0; l < ti.n_layers(); ++l) {
        const LayerInfo &li = ti.layer(l);
        if (li.n_modules() == 0) continue;
        const bool brl = li.is_barrel();
        const double target = brl ? li.r_mean() : li.z_mean();

        // first sign change of (coord - target) along s
        int ihit = -1;
        for (int i = 0; i + 1 < nsamp; ++i) {
          const double c1 = brl ? rs[i] : zs[i];
          const double c2 = brl ? rs[i + 1] : zs[i + 1];
          if ((c1 - target) * (c2 - target) <= 0.0 && i > 0) {
            // check the other coordinate is inside the layer at this crossing
            const double sc = bisect_surface(h, brl, target, i * ds, (i + 1) * ds);
            double p[3], m[3]; h.at_s(sc, p, m);
            const double rr = std::hypot(p[0], p[1]);
            const bool inside = brl ? (p[2] > li.zmin() && p[2] < li.zmax())
                                    : (rr > li.rin() && rr < li.rout() && !li.is_in_r_hole(rr));
            if (inside) { ihit = i; break; }
          }
        }
        if (ihit < 0) continue;

        const double sc = bisect_surface(h, brl, target, ihit * ds, (ihit + 1) * ds);
        double p[3], m[3]; h.at_s(sc, p, m);
        const int sid = nearest_module(li, p);
        if (sid < 0) continue;
        const ModuleInfo &mi = li.module_info(sid);
        double sp;
        if (!refine_to_plane(h, mi.pos, mi.zdir, sc, sp)) continue;
        if (sp < 0.5) continue;                  // do not put a hit at the origin
        if (std::abs(sp - sc) > 8.0) continue;   // plane solve ran away

        Crossing c;
        c.layer = l;
        c.module_sid = sid;
        c.s = sp;
        c.is_barrel = brl;
        h.at_s(sp, c.pos, c.mom);
        out.push_back(c);
      }
      std::sort(out.begin(), out.end(), [](const Crossing &a, const Crossing &b) { return a.s < b.s; });
    }

    // --------------------------------------------------------------------
    // Per-layer-type hit resolution.
    //
    // The numbers are the nominal phase-2 sensor geometry, converted the way
    // LayerOfHits::registerHit() reads them back:  for a strip of half-length L
    // with the crossing uniform along it the covariance is (L^2/3) u u^T, so
    // sigma_along = L / sqrt(3) and registerHit's sqrt(3 * ezz) recovers
    // L * |u_z| exactly.  Feeding sigma_along = L/sqrt(3) here therefore makes
    // hit_q_half_length() come out at the geometric half-length, tilt included.
    //
    // sigma_across is pitch / sqrt(12).  The values:
    //   IT pixel        25 x 100 um             -> 7.2 um  x 28.9 um
    //   TBPS "P"        pitch 100 um, L 0.75 mm -> 28.9 um x 0.0433 cm
    //   TBPS "S"        pitch 100 um, L 12 mm   -> 28.9 um x 0.6928 cm
    //   TOB 2S          pitch  90 um, L 25 mm   -> 26.0 um x 1.4434 cm
    // Phase-2 barrel layer map (RecoTracker/CLAUDE.md): 0-3 PixB, 4-9 TBPS with
    // is_stereo = 1 on the P sensor, 10-15 TOB 2S.  The endcap strip disks are
    // radially about half PS and half 2S; there is no per-hit information here
    // to tell which, so they are taken as PS with the stereo flag choosing P/S.
    // That is an ASSUMPTION, flagged because it changes the endcap numbers --
    // but not the acceptance criteria, which are self-consistency tests against
    // whatever covariance is declared.

    struct HitRes {
      float sig_across;   // along module xdir -- the precise / phi direction
      float sig_along;    // along module ydir -- the strip direction
      float sig_normal;   // along module zdir; only the in-plane 2x2 is ever used
      const char *kind;
      int kind_id;        // 0 pix, 1 PS-P, 2 PS-S, 3 2S
    };
    constexpr int kNKind = 4;
    const char *kKindName[kNKind] = {"pix", "PS-P", "PS-S", "2S"};

    // Set non-zero to force every sensor to this sigma in both in-plane
    // directions.  Used by the "track-covariance-dominated" row of Task 4: with
    // the hit error made negligible, resErr_loc ~ psErrLoc and chi2/ndf becomes
    // a direct measurement of the PROPAGATED covariance rather than of the sum.
    float g_res_override = 0.0f;

    HitRes hit_res_for_layer(const LayerInfo &li) {
      const float p_over_sq12 = 1.0f / std::sqrt(12.0f);
      if (g_res_override > 0.0f)
        return { g_res_override, g_res_override, g_res_override, "ovr", li.is_pixel() ? 0 : 2 };
      if (li.is_pixel())
        return { 25e-4f * p_over_sq12, 100e-4f * p_over_sq12, 10e-4f, "pix", 0 };
      const int l = li.layer_id();
      const bool tob2s = (l >= 10 && l <= 15);
      if (tob2s)
        return { 90e-4f * p_over_sq12, 2.5f / std::sqrt(3.0f), 10e-4f, "2S", 3 };
      if (li.is_stereo())
        return { 100e-4f * p_over_sq12, 0.075f / std::sqrt(3.0f), 10e-4f, "PS-P", 1 };
      return { 100e-4f * p_over_sq12, 1.2f / std::sqrt(3.0f), 10e-4f, "PS-S", 2 };
    }

    // 3x3 global hit covariance from the module frame resolutions.
    // C_glob = R^T diag(sx^2, sy^2, sz^2) R, rows of R being xdir, ydir, zdir.
    void hit_cov_global(const ModuleInfo &mi, const HitRes &hr, float cov[6]) {
      const SVector3 yd = mi.calc_ydir();
      const double ax[3] = {mi.xdir[0], mi.xdir[1], mi.xdir[2]};
      const double ay[3] = {yd[0], yd[1], yd[2]};
      const double az[3] = {mi.zdir[0], mi.zdir[1], mi.zdir[2]};
      const double vx = (double) hr.sig_across * hr.sig_across;
      const double vy = (double) hr.sig_along * hr.sig_along;
      const double vz = (double) hr.sig_normal * hr.sig_normal;
      double c[3][3];
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          c[i][j] = vx * ax[i] * ax[j] + vy * ay[i] * ay[j] + vz * az[i] * az[j];
      // MPlexHS packing, lower triangle row-major: 00, 10, 11, 20, 21, 22
      cov[0] = (float) c[0][0];
      cov[1] = (float) c[1][0];
      cov[2] = (float) c[1][1];
      cov[3] = (float) c[2][0];
      cov[4] = (float) c[2][1];
      cov[5] = (float) c[2][2];
    }

    // --------------------------------------------------------------------
    // Propagation flag helpers.  These build LOCAL PropagationFlags copies and
    // never touch the global PropagationConfig, so nothing has to be restored.

    PropagationFlags pflags_with_material(const PropagationFlags &base, bool on) {
      PropagationFlags pf = base;
      pf.apply_material = on;
      return pf;
    }
    PropagationFlags pflags_with_param_bfield(const PropagationFlags &base, bool on) {
      PropagationFlags pf = base;
      pf.use_param_b_field = on;
      return pf;
    }

    const TrackerInfo &tinfo() { return Config::TrkInfo; }

    // --------------------------------------------------------------------
    // Matriplex helpers: fill every lane identically so that no lane carries
    // garbage (the propagation kernels loop to NN in several places).

    void fill_all_lanes(MPlexLV &par, const float p[6]) {
      for (int n = 0; n < NN; ++n)
        for (int i = 0; i < 6; ++i) par(n, i, 0) = p[i];
    }
    void fill_all_lanes_sym(MPlexLS &err, const float c[21]) {
      for (int n = 0; n < NN; ++n)
        for (int i = 0; i < 21; ++i) err.fArray[i * NN + n] = c[i];
    }
    void fill_all_lanes_hv(MPlexHV &v, const double p[3]) {
      for (int n = 0; n < NN; ++n)
        for (int i = 0; i < 3; ++i) v(n, i, 0) = (float) p[i];
    }
    void fill_all_lanes_hv(MPlexHV &v, const SVector3 &p) {
      for (int n = 0; n < NN; ++n)
        for (int i = 0; i < 3; ++i) v(n, i, 0) = p[i];
    }
    void fill_all_lanes_q(MPlexQI &v, int x) { for (int n = 0; n < NN; ++n) v(n, 0, 0) = x; }
    void fill_all_lanes_q(MPlexQF &v, float x) { for (int n = 0; n < NN; ++n) v(n, 0, 0) = x; }

    // MPlexLS index of (i,j), lower triangle row-major.
    inline int ls_idx(int i, int j) {
      if (i < j) std::swap(i, j);
      return i * (i + 1) / 2 + j;
    }

    // A plausible track covariance at a hit -- used where a test needs an input
    // covariance but its exact value does not matter (Tasks 0 to 2 all check
    // internal consistency, not size).  Scaled to the state so that high-pT and
    // low-pT tracks get comparably conditioned matrices.
    void nominal_cov(const float par[6], float c[21], float scale = 1.0f) {
      std::memset(c, 0, 21 * sizeof(float));
      const float s_pos = 0.02f * scale;                   // 200 um
      const float s_ipt = 0.03f * par[3] * scale;          // 3% of 1/pT
      const float s_ang = 0.002f * scale;                  // 2 mrad
      c[ls_idx(0, 0)] = s_pos * s_pos;
      c[ls_idx(1, 1)] = s_pos * s_pos;
      c[ls_idx(2, 2)] = s_pos * s_pos;
      c[ls_idx(3, 3)] = s_ipt * s_ipt;
      c[ls_idx(4, 4)] = s_ang * s_ang;
      c[ls_idx(5, 5)] = s_ang * s_ang;
      // A little correlation, so that a jacobian bug cannot hide behind a
      // diagonal input.
      c[ls_idx(0, 1)] = 0.2f * s_pos * s_pos;
      c[ls_idx(3, 4)] = 0.1f * s_ipt * s_ang;
      c[ls_idx(4, 5)] = 0.15f * s_ang * s_ang;
    }

  }  // anonymous namespace

  // ==========================================================================
  // TASK 0 -- is the material contribution of the window propagation dead?
  //
  // PropErrsArgs::do_propagation_stuff() (MkFinderV2p2Structures.cc:16) passes
  // MPlexHV dummy{0.0f} for both plPnt and plNrm because it supplies sPerp and
  // therefore needs no plane.  finding_inter_layer_pflags has PF_apply_material
  // set, so applyMaterialEffects still runs and computes
  //     invCos = p / |pt cosP nx + pt sinP ny + pz nz|
  // with n = (0,0,0), i.e. p/0.  The claim under test: the result lands only in
  // err(3,3), (3,5), (4,4), (5,5) plus outPar(3,0), and the only consumer,
  // MkBinTrackCovExtract, reads (0,0), (0,1), (1,1), (2,2) -- so the window is
  // untouched and dropping PF_apply_material here is a no-op.
  //
  // The test reproduces that exact call rather than flipping the global config,
  // so nothing needs restoring and the Kalman propagation (which shares the
  // same pflags but passes a REAL normal) is not disturbed.
  // ==========================================================================

  void pkv_task0_pea_material() {
    printf("\n"
           "============================================================================\n"
           "TASK 0 -- pea material contribution: dead or not\n"
           "============================================================================\n");

    const TrackerInfo &ti = tinfo();
    const PropagationFlags base = ti.prop_config().finding_inter_layer_pflags;
    printf("finding_inter_layer_pflags as configured: use_param_b_field=%d apply_material=%d "
           "copy_input_state_on_fail=%d  tracker_info=%p\n",
           (int) base.use_param_b_field, (int) base.apply_material,
           (int) base.copy_input_state_on_fail, (const void *) base.tracker_info);

    const PropagationFlags pf_on  = pflags_with_material(base, true);
    const PropagationFlags pf_off = pflags_with_material(base, false);

    std::vector<SamplePoint> sample;
    build_sample(sample);

    // Counters
    long n_calls = 0;
    long n_window_diff = 0;          // any of (0,0),(0,1),(1,1),(2,2) differs bitwise
    long n_dphi_diff = 0, n_dq_diff = 0;
    long n_ang_nonfinite_on = 0;     // (3,3)/(3,5)/(4,4)/(5,5) non-finite with material ON
    long n_ang_nonfinite_off = 0;
    long n_ang_diff = 0;             // angular block differs at all
    long n_ipt_diff = 0;             // outPar(3,0) differs
    long n_radl_zero = 0;            // material lookup returned an empty bin
    Stats ipt_ratio;                 // outPar(3,0)_on / _off
    Stats ang44_ratio;               // err(4,4)_on / err(4,4)_off, finite cases only

    // Track which elements ever differ, over all calls.
    bool elem_differs[21] = {false};

    for (const SamplePoint &sp : sample) {
      Helix h;
      h.init(sp.pt, sp.eta, sp.chg, sp.phi0);
      std::vector<Crossing> cx;
      find_crossings(h, ti, cx);
      if (cx.size() < 2) continue;

      for (size_t ic = 0; ic + 1 < cx.size(); ++ic) {
        // Previous hit -> next layer's module plane, exactly as the window
        // propagation is set up: tsPar/tsErr at the previous hit, propPar
        // pre-filled with the destination state, sPerp the transverse path.
        float par_in[6], par_dest[6], cov_in[21];
        h.ccs_at_s(cx[ic].s, par_in);
        h.ccs_at_s(cx[ic + 1].s, par_dest);
        nominal_cov(par_in, cov_in);

        const double alpha = h.alpha_of_s(cx[ic + 1].s) - h.alpha_of_s(cx[ic].s);
        const float sPerp_v = (float) (alpha * h.k * h.pt);   // s_T = k pT alpha

        MPlexLV tsPar;  fill_all_lanes(tsPar, par_in);
        MPlexLS tsErr;  fill_all_lanes_sym(tsErr, cov_in);
        MPlexQI tsChg;  fill_all_lanes_q(tsChg, sp.chg);
        MPlexQF sPerp;  fill_all_lanes_q(sPerp, sPerp_v);
        MPlexHV dummy{0.0f};

        MPlexLS errA, errB;
        MPlexLV parA, parB;
        fill_all_lanes(parA, par_dest);
        fill_all_lanes(parB, par_dest);
        MPlexQI ffA{0}, ffB{0};

        propagateHelixToPlaneMPlex(tsErr, tsPar, tsChg, dummy, dummy, &sPerp,
                                   errA, parA, ffA, NN, pf_on, nullptr);
        propagateHelixToPlaneMPlex(tsErr, tsPar, tsChg, dummy, dummy, &sPerp,
                                   errB, parB, ffB, NN, pf_off, nullptr);

        ++n_calls;

        // What the material lookup would have found at the destination.
        {
          const auto mat = ti.material_checked(std::abs(par_dest[2]), std::hypot(par_dest[0], par_dest[1]));
          if (mat.radl < 1e-13f) ++n_radl_zero;
        }

        // Window numbers -- bitwise.
        const int wi[4][2] = {{0, 0}, {0, 1}, {1, 1}, {2, 2}};
        bool wdiff = false;
        for (auto &ij : wi) {
          const float a = errA.constAt(0, ij[0], ij[1]);
          const float b = errB.constAt(0, ij[0], ij[1]);
          if (std::memcmp(&a, &b, sizeof(float)) != 0) wdiff = true;
        }
        if (wdiff) ++n_window_diff;

        // The two derived quantities the window is actually built from.
        // Reproduce MkBins::determine_bin_windows() for dphi_track (the
        // jacobian is evaluated at the smaller of the two radii there; here
        // both variants get the same jacobian, so any difference is the
        // covariance's) and dq_track for a barrel layer.
        {
          MkBinTrackCovExtract exA(errA), exB(errB);
          const float x = par_dest[0], y = par_dest[1];
          const float r2inv = 1.0f / (x * x + y * y);
          MPlexQF jx, jy;
          fill_all_lanes_q(jx, -y * r2inv);
          fill_all_lanes_q(jy,  x * r2inv);
          const float dphiA = 3.0f * std::sqrt(std::abs(exA.calc_err_xy(jx, jy)[0]));
          const float dphiB = 3.0f * std::sqrt(std::abs(exB.calc_err_xy(jx, jy)[0]));
          const float dqA = 3.0f * std::sqrt(std::abs(exA.m_cov_2_2[0]));
          const float dqB = 3.0f * std::sqrt(std::abs(exB.m_cov_2_2[0]));
          if (std::memcmp(&dphiA, &dphiB, sizeof(float)) != 0) ++n_dphi_diff;
          if (std::memcmp(&dqA, &dqB, sizeof(float)) != 0) ++n_dq_diff;
        }

        // The angular block: is it actually inf/nan at run time?
        const int ai[4][2] = {{3, 3}, {3, 5}, {4, 4}, {5, 5}};
        bool anyang = false;
        for (auto &ij : ai) {
          const float a = errA.constAt(0, ij[0], ij[1]);
          const float b = errB.constAt(0, ij[0], ij[1]);
          if (!fin(a)) anyang = true;
          if (std::memcmp(&a, &b, sizeof(float)) != 0) n_ang_diff += 0;  // counted below
        }
        if (anyang) ++n_ang_nonfinite_on;
        {
          bool bad = false;
          for (auto &ij : ai) if (!fin(errB.constAt(0, ij[0], ij[1]))) bad = true;
          if (bad) ++n_ang_nonfinite_off;
        }

        // Full 21-element scan, so the answer is not restricted to the four
        // elements the claim names.
        bool any_diff = false;
        for (int i = 0; i < 6; ++i)
          for (int j = 0; j <= i; ++j) {
            const float a = errA.constAt(0, i, j);
            const float b = errB.constAt(0, i, j);
            if (std::memcmp(&a, &b, sizeof(float)) != 0) { elem_differs[ls_idx(i, j)] = true; any_diff = true; }
          }
        if (any_diff) ++n_ang_diff;

        // outPar: applyMaterialEffects also rewrites ipt (energy loss).
        for (int i = 0; i < 6; ++i) {
          const float a = parA.constAt(0, i, 0), b = parB.constAt(0, i, 0);
          if (i == 3 && std::memcmp(&a, &b, sizeof(float)) != 0) ++n_ipt_diff;
        }
        if (parB.constAt(0, 3, 0) != 0.0f)
          ipt_ratio.add(parA.constAt(0, 3, 0) / parB.constAt(0, 3, 0));
        const float e44A = errA.constAt(0, 4, 4), e44B = errB.constAt(0, 4, 4);
        if (fin(e44A) && fin(e44B) && e44B != 0.0f) ang44_ratio.add(e44A / e44B);
      }
    }

    printf("\nwindow propagations exercised: %ld  (from %zu sampled tracks)\n", n_calls, sample.size());
    printf("material lookups landing in an empty/off-grid bin (radl < 1e-13): %ld (%.1f%%)\n",
           n_radl_zero, 100.0 * n_radl_zero / std::max(1L, n_calls));
    printf("\nTHE ANSWER:\n");
    printf("  MkBinTrackCovExtract inputs err(0,0),(0,1),(1,1),(2,2) differ bitwise: %ld / %ld\n",
           n_window_diff, n_calls);
    printf("  derived dphi_track differs bitwise: %ld / %ld\n", n_dphi_diff, n_calls);
    printf("  derived dq_track   differs bitwise: %ld / %ld\n", n_dq_diff, n_calls);
    printf("  any of the 21 covariance elements differs: %ld / %ld\n", n_ang_diff, n_calls);
    printf("  outPar(3,0) (1/pT, rewritten by the energy-loss term) differs: %ld / %ld\n",
           n_ipt_diff, n_calls);
    printf("\nwhich covariance elements ever differ:\n   ");
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j <= i; ++j) printf(" (%d,%d)%s", i, j, elem_differs[ls_idx(i, j)] ? "*" : " ");
      printf("\n   ");
    }
    printf("\n  ('*' = differs for at least one propagation)\n");

    printf("\nis the inf real at run time (-Ofast implies -ffinite-math-only, under which\n"
           "a guaranteed 1/0 is UB rather than a benign inf):\n");
    printf("  angular block non-finite with material ON : %ld / %ld\n", n_ang_nonfinite_on, n_calls);
    printf("  angular block non-finite with material OFF: %ld / %ld\n", n_ang_nonfinite_off, n_calls);
    if (ang44_ratio.n())
      printf("  err(4,4) ON/OFF ratio where both finite: n=%zu med=%.4g max=%.4g\n",
             ang44_ratio.n(), ang44_ratio.med(), ang44_ratio.max());
    if (ipt_ratio.n())
      printf("  outPar(3,0) ON/OFF ratio: n=%zu med=%.6g p90=%.6g max=%.6g  (1 = untouched)\n",
             ipt_ratio.n(), ipt_ratio.med(), ipt_ratio.pct(90), ipt_ratio.max());

    // ------------------------------------------------------------------
    // Direct probe of applyMaterialEffects with a zero normal.
    //
    // The scan above says nothing at all changes, which is a STRONGER claim
    // than "only the angular terms change and nobody reads them".  So call the
    // material function directly, with a hand-set radl, and watch what it does
    // to a known covariance: once with n = (0,0,0) as the window propagation
    // passes it, once with a real unit normal for comparison.
    {
      printf("\ndirect probe of applyMaterialEffects (radl and xi set by hand):\n");
      const float radl_probe = 0.05f, xi_probe = 1.0e-3f;
      // Scan the normal from exactly zero up to unit length.  If the exact-zero
      // case is a no-op while a tiny-but-finite normal blows the covariance up,
      // the difference is the compiler's, not the arithmetic's -- which is the
      // question -ffinite-math-only raises.
      const struct { float n[3]; const char *nm; } nrms[] = {
        {{0.0f,   0.0f, 0.0f}, "n=(0,0,0)      as pea passes it"},
        {{1e-30f, 0.0f, 0.0f}, "n=(1e-30,0,0)                  "},
        {{1e-10f, 0.0f, 0.0f}, "n=(1e-10,0,0)                  "},
        {{1e-4f,  0.0f, 0.0f}, "n=(1e-4,0,0)                   "},
        {{0.1f,   0.0f, 0.0f}, "n=(0.1,0,0)                    "},
        {{1.0f,   0.0f, 0.0f}, "n=(1,0,0)      a real normal   "},
      };
      for (auto &nn : nrms) {
        float par[6] = {30.0f, 5.0f, 10.0f, 0.25f, 0.2f, 1.2f};
        float cov[21]; nominal_cov(par, cov);
        MPlexLV oPar; fill_all_lanes(oPar, par);
        MPlexLS oErr; fill_all_lanes_sym(oErr, cov);
        MPlexQF rl, xi, ps;
        fill_all_lanes_q(rl, radl_probe);
        fill_all_lanes_q(xi, xi_probe);
        fill_all_lanes_q(ps, 1.0f);
        MPlexHV nrm; { double d[3] = {nn.n[0], nn.n[1], nn.n[2]}; fill_all_lanes_hv(nrm, d); }
        const float e33_0 = oErr.constAt(0, 3, 3), e44_0 = oErr.constAt(0, 4, 4);
        const float e55_0 = oErr.constAt(0, 5, 5), ipt_0 = oPar.constAt(0, 3, 0);
        applyMaterialEffects(rl, xi, ps, nrm, oErr, oPar, NN);
        (void) e33_0; (void) e55_0;
        printf("  %s: err(3,3) %.6g -> %-12.6g err(4,4) %.6g -> %-12.6g err(5,5) %.6g -> %-12.6g "
               "1/pT %.6g -> %.6g\n",
               nn.nm, e33_0, oErr.constAt(0, 3, 3), e44_0, oErr.constAt(0, 4, 4),
               e55_0, oErr.constAt(0, 5, 5), ipt_0, oPar.constAt(0, 3, 0));
      }
      // And what the two arithmetic pieces actually evaluate to, computed here
      // the same way applyMaterialEffects computes them.
      const float p_probe = 4.0f / std::sin(1.2f);
      const float denom = 0.0f;
      const float invcos = p_probe / std::abs(denom);
      printf("  invCos = p / |pt cosP nx + pt sinP ny + pz nz| with n=0 : p=%.6g denom=%.6g "
             "invCos=%.6g  isfinite=%d\n", p_probe, denom, invcos, (int) fin(invcos));
      const float radl_times = 0.05f * invcos;
      printf("  radL * invCos = %.6g  isfinite=%d ; the bailout tests (radL < 1e-13f) = %d\n",
             radl_times, (int) fin(radl_times), (int) (radl_times < 1e-13f));
    }

    // One worked example, printed in full, so the reader can see the two
    // matrices rather than only the verdict.
    {
      Helix h; h.init(4.0, 0.4, +1, 0.3);
      std::vector<Crossing> cx;
      find_crossings(h, ti, cx);
      if (cx.size() >= 3) {
        float par_in[6], par_dest[6], cov_in[21];
        h.ccs_at_s(cx[1].s, par_in);
        h.ccs_at_s(cx[2].s, par_dest);
        nominal_cov(par_in, cov_in);
        const double alpha = h.alpha_of_s(cx[2].s) - h.alpha_of_s(cx[1].s);
        MPlexLV tsPar;  fill_all_lanes(tsPar, par_in);
        MPlexLS tsErr;  fill_all_lanes_sym(tsErr, cov_in);
        MPlexQI tsChg;  fill_all_lanes_q(tsChg, +1);
        MPlexQF sPerp;  fill_all_lanes_q(sPerp, (float) (alpha * h.k * h.pt));
        MPlexHV dummy{0.0f};
        MPlexLS errA, errB; MPlexLV parA, parB; MPlexQI ffA{0}, ffB{0};
        fill_all_lanes(parA, par_dest); fill_all_lanes(parB, par_dest);
        propagateHelixToPlaneMPlex(tsErr, tsPar, tsChg, dummy, dummy, &sPerp, errA, parA, ffA, NN, pf_on, nullptr);
        propagateHelixToPlaneMPlex(tsErr, tsPar, tsChg, dummy, dummy, &sPerp, errB, parB, ffB, NN, pf_off, nullptr);
        const auto mat_we = ti.material_checked(std::abs(par_dest[2]), std::hypot(par_dest[0], par_dest[1]));
        printf("\nworked example, pT=4 eta=0.4 q=+1, layer %d -> %d:\n", cx[1].layer, cx[2].layer);
        printf("   destination (r,|z|) = (%.3f, %.3f); material bin radl=%.6g bbxi=%.6g\n",
               std::hypot(par_dest[0], par_dest[1]), std::abs(par_dest[2]), mat_we.radl, mat_we.bbxi);
        printf("   %-9s %14s %14s\n", "element", "material ON", "material OFF");
        const int show[8][2] = {{0,0},{0,1},{1,1},{2,2},{3,3},{3,5},{4,4},{5,5}};
        for (auto &ij : show)
          printf("   err(%d,%d)   %14.7g %14.7g\n", ij[0], ij[1],
                 errA.constAt(0, ij[0], ij[1]), errB.constAt(0, ij[0], ij[1]));
        printf("   par(3,0)   %14.7g %14.7g   (1/pT)\n", parA.constAt(0, 3, 0), parB.constAt(0, 3, 0));
      }
    }
  }


  // ==========================================================================
  // TASK 4 -- synthetic hit sequence through the backward-fit kernel.
  //
  // What is reproduced, and what is not.  MkFinder::bkFitFitTracksProp2Plane()
  // (MkFinder.cc:2356) is a HoTNode walk wrapped around exactly two calls:
  //
  //     propagateHelixToPlaneMPlex(m_Err[iC], m_Par[iC], m_Chg, plPnt, plNrm,
  //                                nullptr, m_Err[iP], m_Par[iP], m_FailFlag,
  //                                N_proc, backward_fit_pflags, nullptr);
  //     kalmanOperationPlaneLocal(KFO_Calculate_Chi2 | KFO_Update_Params |
  //                               KFO_Local_Cov, m_Err[iP], m_Par[iP], m_Chg,
  //                               m_msErr, m_msPar, plNrm, plDir, plPnt,
  //                               m_Err[iC], m_Par[iC], tmp_chi2, N_proc);
  //     kalmanCheckChargeFlip(m_Par[iC], m_Chg, N_proc);
  //     m_Chi2.add(tmp_chi2);
  //
  // Those four lines, in that order, with the same flags, are what runs below.
  // Everything else in that function is hit navigation and copy-out, which a
  // synthetic sequence supplies directly.  Driving the real function instead
  // would require a fabricated EventOfHits with registered hits and TrackCands
  // carrying HoTNode chains; that buys no additional coverage of the fit itself
  // and a lot of scaffolding that could itself be wrong.
  //
  // NOTE on the two "recovery branches": in the current source they are NOT
  // dprintf-only, they are commented out entirely (MkFinder.cc, the block after
  // "// Fixup for failed propagation." inside /* ... */, plus the
  // "PROP-FAIL-ENABLE We do not check for pfailed propagation here" comment).
  // So there is no live branch to make visible and nothing to count in situ.
  // The two conditions are evaluated and counted here instead.
  // ==========================================================================

  namespace {

    enum MatMode {
      MM_None = 0,        // material off entirely
      MM_AllBefore,       // PF_apply_material on the propagation -- what runs today
      MM_AllBeforeManual, // same, applied by hand after a material-off propagation
      MM_HalfHalf         // half before the update, half after
    };

    const char *matmode_name(MatMode m) {
      switch (m) {
        case MM_None: return "material OFF";
        case MM_AllBefore: return "all before update (current code)";
        case MM_AllBeforeManual: return "all before update (by hand)";
        case MM_HalfHalf: return "half before / half after";
      }
      return "?";
    }

    // One hit of a synthetic sequence.  The module plane is a REAL phase-2
    // module (nearest to the helix crossing), so normals carry the actual TBPS
    // tilt; the hit sits at the exact analytic helix-plane intersection.
    struct SeqHit {
      int layer = -1, sid = -1;
      double s = 0;
      double pos_true[3] = {0, 0, 0};
      SVector3 plpnt, plnrm, pldir, plydir;
      HitRes hr = {0, 0, 0, "?", 0};
      float cov[6] = {0, 0, 0, 0, 0, 0};   // MPlexHS packing, global frame
      // |cos| of the angle between the truth momentum and the module normal.
      // 1 = normal incidence.  Task 2 finds the reported CCS position
      // covariance blind along the MOMENTUM while the state it describes is
      // constrained to the module PLANE; those two subspaces coincide only at
      // normal incidence, so this is the variable to bin on if that mismatch is
      // what drives the excess chi2.
      float cos_inc = 1.0f;
    };

    // |cos(incidence)| bins.
    const double kIncEdges[] = {0.0, 0.5, 0.7, 0.85, 0.93, 0.97, 0.99, 1.0001};
    constexpr int kNInc = sizeof(kIncEdges) / sizeof(kIncEdges[0]) - 1;
    int inc_bin(double c) {
      for (int i = 0; i < kNInc; ++i) if (c < kIncEdges[i + 1]) return i;
      return kNInc - 1;
    }

    // Per-hit, per-lane hit positions: hp[ihit][3*NN], laid out (i, lane).
    using HitPosArr = std::vector<std::vector<float>>;

    struct FitOut {
      MPlexQF chi2_total{0.0f};
      std::vector<MPlexQF> chi2_hit;     // one per hit, in fit order
      int n_prop_fail = 0;               // FailFlag set, any lane
      std::vector<int> fail_hits;        // fit-order indices where it happened
      int n_chi2_over_200 = 0;           // the commented-out recovery condition
      int n_chi2_negative = 0;
      int n_chi2_nonfinite = 0;
      MPlexLV par_end;
      MPlexLS err_end;
      MPlexQI chg_end{0};
    };

    // Apply a fraction of the destination bin's material as one thin scatterer.
    // Halving radl (and xi) is the physically right way to split a scatterer in
    // two: thetaMSC2 goes as radL * (1 + 0.038 log radL)^2, so half the
    // thickness twice is not exactly the whole thickness once -- which is the
    // point, that is what a two-slab model says.
    void apply_material_fraction(const TrackerInfo &ti, float frac,
                                 const MPlexHV &plNrm, const MPlexLV &refPar,
                                 MPlexLS &err, MPlexLV &par, int n_proc) {
      if (frac <= 0.0f) return;
      MPlexQF rl, xi, ps;
      for (int n = 0; n < NN; ++n) {
        const float z = par.constAt(n, 2, 0);
        const float r = std::hypot(par.constAt(n, 0, 0), par.constAt(n, 1, 0));
        const auto mat = ti.material_checked(std::abs(z), r);
        rl(n, 0, 0) = frac * mat.radl;
        xi(n, 0, 0) = frac * mat.bbxi;
        // propSign: +1 if the step went along the momentum. The propagation
        // itself uses sign(pathL); reconstruct it from the displacement, which
        // is the same thing and does not need pathL plumbed out.
        const float dx = par.constAt(n, 0, 0) - refPar.constAt(n, 0, 0);
        const float dy = par.constAt(n, 1, 0) - refPar.constAt(n, 1, 0);
        const float dz = par.constAt(n, 2, 0) - refPar.constAt(n, 2, 0);
        const float ph = refPar.constAt(n, 4, 0), th = refPar.constAt(n, 5, 0);
        const float d = dx * std::cos(ph) * std::sin(th) + dy * std::sin(ph) * std::sin(th) + dz * std::cos(th);
        ps(n, 0, 0) = (d > 0.0f) ? 1.0f : -1.0f;
      }
      applyMaterialEffects(rl, xi, ps, plNrm, err, par, n_proc);
    }

    // The fit.  hits are in FIT order, i.e. outermost first.
    void run_bkfit(const TrackerInfo &ti,
                   const std::vector<SeqHit> &hits, const HitPosArr &hp,
                   const MPlexLV &par0, const MPlexLS &err0, const MPlexQI &chg0,
                   MatMode mm, const PropagationFlags &pf_base, FitOut &out) {
      const int nh = (int) hits.size();
      out.chi2_hit.assign(nh, MPlexQF{0.0f});
      out.chi2_total = MPlexQF{0.0f};
      out.n_prop_fail = out.n_chi2_over_200 = out.n_chi2_negative = out.n_chi2_nonfinite = 0;
      out.fail_hits.clear();

      const PropagationFlags pf = pflags_with_material(pf_base, mm == MM_AllBefore);

      MPlexLS errC = err0, errP;
      MPlexLV parC = par0, parP;
      MPlexQI chg = chg0;

      for (int ih = 0; ih < nh; ++ih) {
        const SeqHit &h = hits[ih];
        MPlexHV plPnt, plNrm, plDir;
        fill_all_lanes_hv(plPnt, h.plpnt);
        fill_all_lanes_hv(plNrm, h.plnrm);
        fill_all_lanes_hv(plDir, h.pldir);

        MPlexHV msPar;
        MPlexHS msErr;
        for (int n = 0; n < NN; ++n) {
          for (int i = 0; i < 3; ++i) msPar(n, i, 0) = hp[ih][i * NN + n];
          for (int i = 0; i < 6; ++i) msErr.fArray[i * NN + n] = h.cov[i];
        }

        MPlexQI ff{0};
        propagateHelixToPlaneMPlex(errC, parC, chg, plPnt, plNrm, nullptr,
                                   errP, parP, ff, NN, pf, nullptr);
        for (int n = 0; n < NN; ++n) if (ff(n, 0, 0)) { ++out.n_prop_fail; out.fail_hits.push_back(ih); break; }

        if (mm == MM_AllBeforeManual) apply_material_fraction(ti, 1.0f, plNrm, parC, errP, parP, NN);
        else if (mm == MM_HalfHalf)   apply_material_fraction(ti, 0.5f, plNrm, parC, errP, parP, NN);

        MPlexQF chi2{0.0f};
        kalmanOperationPlaneLocal(KFO_Calculate_Chi2 | KFO_Update_Params | KFO_Local_Cov,
                                  errP, parP, chg, msErr, msPar, plNrm, plDir, plPnt,
                                  errC, parC, chi2, NN);
        kalmanCheckChargeFlip(parC, chg, NN);

        if (mm == MM_HalfHalf) apply_material_fraction(ti, 0.5f, plNrm, parP, errC, parC, NN);

        out.chi2_hit[ih] = chi2;
        for (int n = 0; n < NN; ++n) {
          const float c = chi2(n, 0, 0);
          if (!fin(c)) ++out.n_chi2_nonfinite;
          else {
            if (c < 0.0f) ++out.n_chi2_negative;
            if (c > 200.0f) ++out.n_chi2_over_200;
          }
          out.chi2_total(n, 0, 0) += c;
        }
      }
      out.par_end = parC;
      out.err_end = errC;
      out.chg_end = chg;
    }

    // Build the hit sequence for one helix. Returns hits in OUTWARD order.
    bool build_sequence(const Helix &h, const TrackerInfo &ti, std::vector<SeqHit> &out) {
      std::vector<Crossing> cx;
      find_crossings(h, ti, cx);
      out.clear();
      for (const Crossing &c : cx) {
        const LayerInfo &li = ti.layer(c.layer);
        const ModuleInfo &mi = li.module_info(c.module_sid);
        SeqHit sh;
        sh.layer = c.layer;
        sh.sid = c.module_sid;
        sh.s = c.s;
        for (int i = 0; i < 3; ++i) sh.pos_true[i] = c.pos[i];
        sh.plpnt = mi.pos;
        sh.plnrm = mi.zdir;
        sh.pldir = mi.xdir;
        sh.plydir = mi.calc_ydir();
        sh.hr = hit_res_for_layer(li);
        hit_cov_global(mi, sh.hr, sh.cov);
        {
          const double pm = std::sqrt(c.mom[0] * c.mom[0] + c.mom[1] * c.mom[1] + c.mom[2] * c.mom[2]);
          sh.cos_inc = (float) std::abs((c.mom[0] * mi.zdir[0] + c.mom[1] * mi.zdir[1] +
                                         c.mom[2] * mi.zdir[2]) / pm);
        }
        out.push_back(sh);
      }
      return out.size() >= 4;
    }

    // Fill hit positions: exact, or smeared in the module plane by the very
    // covariance the hit declares (so criterion (b) is a statement about the
    // fit, not about the smearing model).
    void fill_hit_positions(const std::vector<SeqHit> &hits, bool smear,
                            std::mt19937_64 &rng, HitPosArr &hp) {
      std::normal_distribution<float> g(0.0f, 1.0f);
      hp.assign(hits.size(), std::vector<float>(3 * NN, 0.0f));
      for (size_t ih = 0; ih < hits.size(); ++ih) {
        const SeqHit &h = hits[ih];
        for (int n = 0; n < NN; ++n) {
          double p[3] = {h.pos_true[0], h.pos_true[1], h.pos_true[2]};
          if (smear) {
            const float a = h.hr.sig_across * g(rng);
            const float b = h.hr.sig_along * g(rng);
            for (int i = 0; i < 3; ++i) p[i] += a * h.pldir[i] + b * h.plydir[i];
          }
          for (int i = 0; i < 3; ++i) hp[ih][i * NN + n] = (float) p[i];
        }
      }
    }

  }  // anonymous namespace

  namespace {
    // pT x eta grid of accumulators, so every number below can be read against
    // where in the state space it came from.
    struct Grid {
      Stats c[kNPt][kNEta];
      void add(double pt, double eta, double x) {
        int ip = pt_bin(pt);
        int ie = 0;
        for (int i = 0; i < kNEta; ++i) if (std::abs(eta - kEtaGrid[i]) < 1e-6) ie = i;
        c[ip][ie].add(x);
      }
      void print(const char *what, double (Stats::*f)()) {
        printf("  %-16s", what);
        for (int ie = 0; ie < kNEta; ++ie) printf("  eta=%.1f", kEtaGrid[ie]);
        printf("\n");
        for (int ip = 0; ip < kNPt; ++ip) {
          printf("  pT=%-6.4g       ", kPtGrid[ip]);
          for (int ie = 0; ie < kNEta; ++ie) {
            if (c[ip][ie].n()) printf(" %8.3g", (c[ip][ie].*f)());
            else                printf(" %8s", "-");
          }
          printf("\n");
        }
      }
    };

    // |dalpha| bins for the step between consecutive hits -- the quantity the
    // plane solve's convergence depends on (Config::nSStepsInProp2Plane = 2, so
    // one refinement iteration beyond the initial quadratic guess).
    const double kDalEdges[] = {0.0, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8, 1e9};
    constexpr int kNDal = sizeof(kDalEdges) / sizeof(kDalEdges[0]) - 1;
    int dal_bin(double d) {
      for (int i = 0; i < kNDal; ++i) if (d < kDalEdges[i + 1]) return i;
      return kNDal - 1;
    }
  }  // anonymous namespace

  void pkv_task4_sequence() {
    printf("\n"
           "============================================================================\n"
           "TASK 4 -- synthetic hit sequence through the backward-fit kernel\n"
           "============================================================================\n");

    const TrackerInfo &ti = tinfo();
    const PropagationFlags bf_base = ti.prop_config().backward_fit_pflags;
    printf("backward_fit_pflags as configured: use_param_b_field=%d apply_material=%d\n",
           (int) bf_base.use_param_b_field, (int) bf_base.apply_material);
    printf("Config::usePropToPlane=%d usePtMultScat=%d nSStepsInProp2Plane=%d\n",
           (int) Config::usePropToPlane, (int) Config::usePtMultScat, Config::nSStepsInProp2Plane);
    printf("truth trajectory: exact uniform-B helix, B = %.4f T; hits at the exact\n"
           "analytic intersection with a REAL phase-2 module plane (nearest module to the\n"
           "layer crossing), so the module tilt is carried by plNrm/plDir as in the fit.\n",
           Config::Bfield);
    {
      // Demonstrate the -ffinite-math-only trap, since every number here depends
      // on getting the finiteness test right.
      volatile float zero = 0.0f, one = 1.0f;
      const float inf_rt = one / zero;
      printf("build check: a runtime 1/0 gives %g ; std::isfinite says %d, "
             "mkfit::isFinite says %d\n", inf_rt, (int) std::isfinite(inf_rt), (int) mkfit::isFinite(inf_rt));
    }

    // The truth is a uniform-B helix, so the reference configuration must also
    // use a uniform field: use_param_b_field is OFF for (a) and (b), and its
    // cost is measured separately at the end.
    const PropagationFlags pf_uniform = pflags_with_param_bfield(bf_base, false);
    const PropagationFlags pf_paramB  = pflags_with_param_bfield(bf_base, true);

    std::vector<SamplePoint> sample;
    build_sample(sample);

    printf("\nhit resolutions used (sigma_across = pitch/sqrt(12) along module xdir,\n"
           "sigma_along = L/sqrt(3) along ydir, L the sensor half-length):\n");
    {
      const int show[] = {0, 4, 5, 10, 16, 28, 29};
      for (int l : show) {
        if (l >= ti.n_layers()) continue;
        const LayerInfo &li = ti.layer(l);
        const HitRes hr = hit_res_for_layer(li);
        printf("  layer %2d %-5s pix=%d stereo=%d brl=%d : sigma_across=%9.2e cm  sigma_along=%9.2e cm"
               "  -> q_half_len would read %.4f cm\n",
               l, hr.kind, (int) li.is_pixel(), (int) li.is_stereo(), (int) li.is_barrel(),
               hr.sig_across, hr.sig_along,
               (li.is_pixel() ? 3.0f : std::sqrt(3.0f)) * hr.sig_along);
      }
    }

    auto make_prior = [](const float par[6], float cov[21]) {
      // Loose but not absurd: sigma_pos 1 cm, sigma_ang 50 mrad, sigma_ipt 50%.
      std::memset(cov, 0, 21 * sizeof(float));
      const float sp = 1.0f, sa = 0.05f, si = 0.5f * par[3];
      cov[ls_idx(0, 0)] = sp * sp;
      cov[ls_idx(1, 1)] = sp * sp;
      cov[ls_idx(2, 2)] = sp * sp;
      cov[ls_idx(3, 3)] = si * si;
      cov[ls_idx(4, 4)] = sa * sa;
      cov[ls_idx(5, 5)] = sa * sa;
    };

    std::mt19937_64 rng(20260910u);
    std::normal_distribution<float> gauss(0.0f, 1.0f);

    // ====================================================================
    // (a) hits EXACTLY on the helix, material OFF, uniform B -> chi2 ~ 0
    // ====================================================================
    {
      Grid g_tot, g_max;
      Stats by_seq[24], by_seq_hipt[24];
      Stats by_dal[kNDal];
      long n_fits = 0, n_hits_tot = 0, n_over200 = 0, n_neg = 0, n_nonfin = 0, n_propfail = 0;
      Stats tot_all;
      struct Worst { double chi2, pt, eta, phi; int chg, ihit, layer; double dal; };
      std::vector<Worst> worst;

      for (const SamplePoint &sp : sample) {
        Helix h; h.init(sp.pt, sp.eta, sp.chg, sp.phi0);
        std::vector<SeqHit> seq;
        if (!build_sequence(h, ti, seq)) continue;
        std::vector<SeqHit> fs(seq.rbegin(), seq.rend());
        const int nh = (int) fs.size();

        float par0[6], cov0[21];
        h.ccs_at_s(fs.front().s, par0);
        make_prior(par0, cov0);
        MPlexLV P;  fill_all_lanes(P, par0);
        MPlexLS E;  fill_all_lanes_sym(E, cov0);
        MPlexQI Q;  fill_all_lanes_q(Q, sp.chg);

        HitPosArr hp;
        fill_hit_positions(fs, false, rng, hp);
        FitOut fo;
        run_bkfit(ti, fs, hp, P, E, Q, MM_None, pf_uniform, fo);

        ++n_fits; n_hits_tot += nh;
        n_over200 += fo.n_chi2_over_200; n_neg += fo.n_chi2_negative;
        n_nonfin += fo.n_chi2_nonfinite; n_propfail += fo.n_prop_fail;
        const double tot = fo.chi2_total(0, 0, 0);
        tot_all.add(tot);
        g_tot.add(sp.pt, sp.eta, tot);
        double mx = 0;
        for (int ih = 0; ih < nh; ++ih) {
          const double c = fo.chi2_hit[ih](0, 0, 0);
          if (fin(c) && c > mx) mx = c;
          if (ih < 24) {
            by_seq[ih].add(c);
            if (sp.pt >= 4.0) by_seq_hipt[ih].add(c);
          }
          // |dalpha| of the step that led into this hit
          double dal = 0;
          if (ih > 0) dal = std::abs(h.alpha_of_s(fs[ih].s) - h.alpha_of_s(fs[ih - 1].s));
          by_dal[dal_bin(dal)].add(c);
          if (fin(c) && c > 1.0)
            worst.push_back({c, sp.pt, sp.eta, sp.phi0, sp.chg, ih, fs[ih].layer, dal});
        }
        g_max.add(sp.pt, sp.eta, mx);
      }

      printf("\n--- (a) hits EXACTLY on the helix, material OFF, uniform B -------------------\n");
      printf("ACCEPTANCE: every per-hit chi2 increment must be ~0 (float precision only).\n");
      printf("fits=%ld  hits=%ld  <hits/fit>=%.1f  (phase-2 has 16 barrel layers, so a\n"
             "barrel track legitimately crosses 16 -- these are geometric crossings, not\n"
             "a reconstruction efficiency)\n",
             n_fits, n_hits_tot, (double) n_hits_tot / std::max(1L, n_fits));
      printf("total chi2 over the sequence: med=%.4g p90=%.4g p99=%.4g max=%.4g  (n_bad=%d)\n",
             tot_all.med(), tot_all.pct(90), tot_all.pct(99), tot_all.max(), tot_all.n_bad);
      printf("counters over all hits: chi2>200 %ld, chi2<0 %ld, non-finite %ld, prop-fail %ld\n",
             n_over200, n_neg, n_nonfin, n_propfail);

      printf("\ntotal chi2, MEDIAN, by (pT, eta):\n");
      g_tot.print("", &Stats::med);
      printf("\nworst single-hit chi2 in the sequence, by (pT, eta):\n");
      g_max.print("", &Stats::max);

      printf("\nper-hit chi2 along the fit sequence (0 = outermost hit):\n");
      printf("  all pT   median:");
      for (int i = 0; i < 17; ++i) if (by_seq[i].n()) printf(" %d:%.2g", i, by_seq[i].med());
      printf("\n  all pT   max   :");
      for (int i = 0; i < 17; ++i) if (by_seq[i].n()) printf(" %d:%.2g", i, by_seq[i].max());
      printf("\n  pT>=4    median:");
      for (int i = 0; i < 17; ++i) if (by_seq_hipt[i].n()) printf(" %d:%.2g", i, by_seq_hipt[i].med());
      printf("\n  pT>=4    max   :");
      for (int i = 0; i < 17; ++i) if (by_seq_hipt[i].n()) printf(" %d:%.2g", i, by_seq_hipt[i].max());
      printf("\n");

      printf("\nper-hit chi2 vs |dalpha| of the step into that hit:\n");
      printf("  %-16s %8s %10s %10s %10s\n", "|dalpha| bin", "n", "median", "p99", "max");
      for (int i = 0; i < kNDal; ++i) {
        if (!by_dal[i].n()) continue;
        char lab[64];
        if (i == kNDal - 1) snprintf(lab, sizeof(lab), ">%.2f", kDalEdges[i]);
        else snprintf(lab, sizeof(lab), "%.2f-%.2f", kDalEdges[i], kDalEdges[i + 1]);
        printf("  %-16s %8zu %10.3g %10.3g %10.3g\n", lab, by_dal[i].n(),
               by_dal[i].med(), by_dal[i].pct(99), by_dal[i].max());
      }

      std::sort(worst.begin(), worst.end(), [](const Worst &a, const Worst &b) { return a.chi2 > b.chi2; });
      printf("\nworst offenders (exact hits should give chi2 = 0), top 12 of %zu with chi2 > 1:\n",
             worst.size());
      printf("  %10s %6s %5s %6s %3s %5s %6s %9s\n", "chi2", "pT", "eta", "phi", "q", "ihit", "layer", "|dalpha|");
      for (size_t i = 0; i < worst.size() && i < 12; ++i)
        printf("  %10.3g %6.4g %5.1f %6.2f %3d %5d %6d %9.4f\n", worst[i].chi2, worst[i].pt,
               worst[i].eta, worst[i].phi, worst[i].chg, worst[i].ihit, worst[i].layer, worst[i].dal);
    }

    // ====================================================================
    // (b) hits smeared by their own covariance; start state DRAWN from its
    //     own prior, so the expectation is exactly chi2 with ndf = 2*n_hits.
    // ====================================================================
    struct BRes { double med_ndf, med_ndf_hipt, frac_over200; };

    auto run_pass_b = [&](const char *label, MatMode mm, const PropagationFlags &pf, int nrep) -> BRes {
      Stats ndf_all, ndf_hipt;
      Grid g_ndf;
      Stats by_seq[24], by_seq_hipt[24];
      Stats by_nh[24];
      Stats by_kind[kNKind], by_dal[kNDal], by_inc[kNInc];
      Stats tot_chi2;
      long n_tot_gt_200 = 0, n_tot_gt_500 = 0, n_tot_gt_2000 = 0;
      long n_hits_tot = 0, n_over200 = 0, n_neg = 0, n_nonfin = 0, n_propfail = 0, n_fits = 0;
      for (const SamplePoint &sp : sample) {
        Helix h; h.init(sp.pt, sp.eta, sp.chg, sp.phi0);
        std::vector<SeqHit> seq;
        if (!build_sequence(h, ti, seq)) continue;
        std::vector<SeqHit> fs(seq.rbegin(), seq.rend());
        const int nh = (int) fs.size();

        float par_true[6], cov0[21];
        h.ccs_at_s(fs.front().s, par_true);
        make_prior(par_true, cov0);

        std::vector<double> dal(nh, 0.0);
        for (int ih = 1; ih < nh; ++ih)
          dal[ih] = std::abs(h.alpha_of_s(fs[ih].s) - h.alpha_of_s(fs[ih - 1].s));

        for (int rep = 0; rep < nrep; ++rep) {
          MPlexLV P;
          MPlexLS E;  fill_all_lanes_sym(E, cov0);
          MPlexQI Q;  fill_all_lanes_q(Q, sp.chg);
          for (int n = 0; n < NN; ++n) {
            for (int i = 0; i < 6; ++i) {
              const float sig = std::sqrt(cov0[ls_idx(i, i)]);
              P(n, i, 0) = par_true[i] + sig * gauss(rng);
            }
            if (P(n, 3, 0) <= 1e-4f) P(n, 3, 0) = 1e-4f;   // 1/pT is kept positive
          }
          HitPosArr hp;
          fill_hit_positions(fs, true, rng, hp);
          FitOut fo;
          run_bkfit(ti, fs, hp, P, E, Q, mm, pf, fo);

          n_fits += NN; n_hits_tot += (long) NN * nh;
          n_over200 += fo.n_chi2_over_200; n_neg += fo.n_chi2_negative;
          n_nonfin += fo.n_chi2_nonfinite; n_propfail += fo.n_prop_fail;
          for (int n = 0; n < NN; ++n) {
            const double r = fo.chi2_total(n, 0, 0) / (2.0 * nh);
            ndf_all.add(r);
            g_ndf.add(sp.pt, sp.eta, r);
            if (sp.pt >= 4.0) ndf_hipt.add(r);
            if (nh < 24) by_nh[nh].add(r);
            const double tc = fo.chi2_total(n, 0, 0);
            tot_chi2.add(tc);
            if (fin(tc)) {
              if (tc > 200.0) ++n_tot_gt_200;
              if (tc > 500.0) ++n_tot_gt_500;
              if (tc > 2000.0) ++n_tot_gt_2000;
            }
            for (int ih = 0; ih < nh && ih < 24; ++ih) {
              const double c = fo.chi2_hit[ih](n, 0, 0);
              by_seq[ih].add(c);
              if (sp.pt >= 4.0) by_seq_hipt[ih].add(c);
              by_kind[fs[ih].hr.kind_id].add(c);
              by_dal[dal_bin(dal[ih])].add(c);
              by_inc[inc_bin(fs[ih].cos_inc)].add(c);
            }
          }
        }
      }
      printf("\n--- (b) %s ---\n", label);
      printf("mode: %s ; use_param_b_field=%d ; fits=%ld hits=%ld\n",
             matmode_name(mm), (int) pf.use_param_b_field, n_fits, n_hits_tot);
      printf("chi2/ndf (ndf = 2*n_hits), expectation 1:\n");
      printf("  all pT : med=%.4g  p90=%.4g  p99=%.4g  max=%.4g  n_nonfinite=%d\n",
             ndf_all.med(), ndf_all.pct(90), ndf_all.pct(99), ndf_all.max(), ndf_all.n_bad);
      printf("  pT>=4  : med=%.4g  p90=%.4g  p99=%.4g  max=%.4g\n",
             ndf_hipt.med(), ndf_hipt.pct(90), ndf_hipt.pct(99), ndf_hipt.max());
      printf("counters: chi2>200 %ld (%.3f%% of hits), chi2<0 %ld, non-finite %ld, prop-fail %ld\n",
             n_over200, 100.0 * n_over200 / std::max(1L, n_hits_tot), n_neg, n_nonfin, n_propfail);
      printf("chi2/ndf MEDIAN by (pT, eta):\n");
      g_ndf.print("", &Stats::med);
      printf("per-hit chi2 along the sequence (2-dof chi2 has median 1.386, p90 4.605):\n");
      printf("  all pT median:");
      for (int i = 0; i < 17; ++i) if (by_seq[i].n()) printf(" %d:%.3g", i, by_seq[i].med());
      printf("\n  all pT p90   :");
      for (int i = 0; i < 17; ++i) if (by_seq[i].n()) printf(" %d:%.3g", i, by_seq[i].pct(90));
      printf("\n  pT>=4 median :");
      for (int i = 0; i < 17; ++i) if (by_seq_hipt[i].n()) printf(" %d:%.3g", i, by_seq_hipt[i].med());
      printf("\n  pT>=4 p90    :");
      for (int i = 0; i < 17; ++i) if (by_seq_hipt[i].n()) printf(" %d:%.3g", i, by_seq_hipt[i].pct(90));
      printf("\nchi2/ndf median vs n_hits (accumulation with hit count?):\n  ");
      for (int i = 0; i < 24; ++i) if (by_nh[i].n() > 200) printf(" n=%d:%.3g", i, by_nh[i].med());
      printf("\nTOTAL track chi2 (the number quoted as O(2000) over 10-12 hits):\n");
      printf("  med=%.4g p90=%.4g p99=%.4g max=%.4g ;  >200: %.3f%%  >500: %.3f%%  >2000: %.3f%%\n",
             tot_chi2.med(), tot_chi2.pct(90), tot_chi2.pct(99), tot_chi2.max(),
             100.0 * n_tot_gt_200 / std::max(1L, n_fits), 100.0 * n_tot_gt_500 / std::max(1L, n_fits),
             100.0 * n_tot_gt_2000 / std::max(1L, n_fits));
      printf("per-hit chi2 by SENSOR TYPE (separates 'position in the sequence' from\n"
             "'hit precision' -- the outer hits of a barrel track are the coarse 2S ones):\n");
      for (int k = 0; k < kNKind; ++k) {
        if (!by_kind[k].n()) continue;
        printf("  %-5s n=%-8zu med=%-8.3g p90=%-9.3g p99=%-10.3g max=%.3g\n", kKindName[k],
               by_kind[k].n(), by_kind[k].med(), by_kind[k].pct(90), by_kind[k].pct(99), by_kind[k].max());
      }
      printf("per-hit chi2 vs |dalpha| of the step into that hit:\n");
      for (int i = 0; i < kNDal; ++i) {
        if (!by_dal[i].n()) continue;
        char lab[64];
        if (i == kNDal - 1) snprintf(lab, sizeof(lab), ">%.2f", kDalEdges[i]);
        else snprintf(lab, sizeof(lab), "%.2f-%.2f", kDalEdges[i], kDalEdges[i + 1]);
        printf("  %-12s n=%-8zu med=%-8.3g p90=%-9.3g p99=%-10.3g max=%.3g\n", lab,
               by_dal[i].n(), by_dal[i].med(), by_dal[i].pct(90), by_dal[i].pct(99), by_dal[i].max());
      }
      printf("per-hit chi2 vs |cos(incidence)| -- 1.0 = the track hits the module head on,\n"
             "which is where the reported covariance's null direction (the momentum) and the\n"
             "state's actual constraint (the module plane) coincide:\n");
      for (int i = 0; i < kNInc; ++i) {
        if (!by_inc[i].n()) continue;
        char lab[64];
        snprintf(lab, sizeof(lab), "%.2f-%.2f", kIncEdges[i], std::min(1.0, kIncEdges[i + 1]));
        printf("  %-12s n=%-8zu med=%-8.3g p90=%-9.3g p99=%-10.3g max=%.3g\n", lab,
               by_inc[i].n(), by_inc[i].med(), by_inc[i].pct(90), by_inc[i].pct(99), by_inc[i].max());
      }
      return {ndf_all.med(), ndf_hipt.med(), 100.0 * n_over200 / std::max(1L, n_hits_tot)};
    };

    const int nrep = 8;
    printf("\n============================================================================\n");
    printf("ACCEPTANCE CRITERION (b): chi2/ndf must be ~1 for the material-OFF row.\n");
    printf("============================================================================\n");
    const BRes b_none = run_pass_b("smeared hits, material OFF, uniform B   [THE CRITERION]",
                                   MM_None, pf_uniform, nrep);

    printf("\n============================================================================\n");
    printf("THE MATERIAL ORDERING A/B, and the field-model cost, same sample\n");
    printf("============================================================================\n");
    const BRes b_before = run_pass_b("ALL material before the update (what the code does today)",
                                     MM_AllBefore, pf_uniform, nrep);
    const BRes b_manual = run_pass_b("all material before, applied BY HAND (self-check of the A/B)",
                                     MM_AllBeforeManual, pf_uniform, nrep);
    const BRes b_half   = run_pass_b("HALF before / HALF after the update",
                                     MM_HalfHalf, pf_uniform, nrep);
    const BRes b_paramB = run_pass_b("material OFF, PARAMETRIZED B field (truth is uniform B)",
                                     MM_None, pf_paramB, nrep);

    // Config::usePtMultScat is FALSE in standalone (Config.cc:8) and TRUE in
    // CMSSW (MkFitGeometryESProducer.cc:563); CMS-phase2.cc sets
    // Config::usePropToPlane but not this one, and mkFit has --use-ptms to set
    // it.  With it false, multiple scattering is added ONLY to err(4,4) and
    // err(5,5) -- so the standalone build puts no scattering into 1/pT at all,
    // which is precisely the element a momentum-resolution or chi2 study cares
    // about.  Repeat the ordering A/B under the CMSSW setting.  Restored below.
    const bool ptms_saved = Config::usePtMultScat;
    Config::usePtMultScat = true;
    printf("\n============================================================================\n");
    printf("SAME A/B with Config::usePtMultScat = true (the CMSSW setting; standalone\n"
           "defaults to false, so multiple scattering normally misses err(3,3)/(3,5))\n");
    printf("============================================================================\n");
    const BRes b_before_ms = run_pass_b("ALL material before update, usePtMultScat=true",
                                        MM_AllBefore, pf_uniform, nrep);
    const BRes b_half_ms   = run_pass_b("HALF before / HALF after, usePtMultScat=true",
                                        MM_HalfHalf, pf_uniform, nrep);
    Config::usePtMultScat = ptms_saved;
    printf("\nConfig::usePtMultScat restored to %d\n", (int) Config::usePtMultScat);

    // Track-covariance-dominated row.  Task 2 finds the analytic error
    // propagation under-reporting the in-plane ALONG-STRIP position variance by
    // ~16%, pT-independent, with the discrepancy a rank-1 term along the
    // direction of motion -- i.e. the free-propagation jacobian without the
    // ds/d(inPar) term the plane constraint implies.  Whether that reaches the
    // Kalman is a separate question: kalmanOperationPlaneLocal transforms CCS ->
    // curvilinear -> local with jacCurv2Loc's cosz terms, which is where CMSSW
    // puts the surface-crossing correction, so it may be recovered there.
    //
    // This row settles it observationally.  With every sensor at 1 um in both
    // in-plane directions the hit error is negligible against the track error,
    // so resErr_loc ~ psErrLoc and chi2/ndf measures the PROPAGATED covariance
    // directly.  1.0 => the correction is recovered; ~1.16 => it is not.
    printf("\n============================================================================\n");
    printf("TRACK-COVARIANCE-DOMINATED row: all sensors forced to 1 um x 1 um, so chi2\n"
           "measures the propagated covariance and not the hit error\n");
    printf("============================================================================\n");
    g_res_override = 1e-4f;
    const BRes b_fine_none = run_pass_b("1 um hits, material OFF, uniform B", MM_None, pf_uniform, nrep);
    const BRes b_fine_bef  = run_pass_b("1 um hits, all material before update", MM_AllBefore, pf_uniform, nrep);
    g_res_override = 0.0f;

    printf("\n--- summary -------------------------------------------------------------\n");
    printf("  %-46s %10s %10s %12s\n", "configuration", "med c2/ndf", "pT>=4", "% hits >200");
    printf("  %-46s %10.4g %10.4g %12.3f\n", "material OFF  [criterion (b), want 1]",
           b_none.med_ndf, b_none.med_ndf_hipt, b_none.frac_over200);
    printf("  %-46s %10.4g %10.4g %12.3f\n", "all material before update (current code)",
           b_before.med_ndf, b_before.med_ndf_hipt, b_before.frac_over200);
    printf("  %-46s %10.4g %10.4g %12.3f\n", "all material before, by hand (self-check)",
           b_manual.med_ndf, b_manual.med_ndf_hipt, b_manual.frac_over200);
    printf("  %-46s %10.4g %10.4g %12.3f\n", "half before / half after",
           b_half.med_ndf, b_half.med_ndf_hipt, b_half.frac_over200);
    printf("  %-46s %10.4g %10.4g %12.3f\n", "material OFF, parametrized B field",
           b_paramB.med_ndf, b_paramB.med_ndf_hipt, b_paramB.frac_over200);
    printf("  %-46s %10.4g %10.4g %12.3f\n", "all before, usePtMultScat=true (CMSSW)",
           b_before_ms.med_ndf, b_before_ms.med_ndf_hipt, b_before_ms.frac_over200);
    printf("  %-46s %10.4g %10.4g %12.3f\n", "half/half, usePtMultScat=true (CMSSW)",
           b_half_ms.med_ndf, b_half_ms.med_ndf_hipt, b_half_ms.frac_over200);
    printf("  %-46s %10.4g %10.4g %12.3f\n", "1um hits, material OFF (cov-dominated)",
           b_fine_none.med_ndf, b_fine_none.med_ndf_hipt, b_fine_none.frac_over200);
    printf("  %-46s %10.4g %10.4g %12.3f\n", "1um hits, all material before",
           b_fine_bef.med_ndf, b_fine_bef.med_ndf_hipt, b_fine_bef.frac_over200);

    printf("\nHOOK -- the pull test against MC truth attaches here.\n"
           "It is deliberately NOT implemented: it needs a clean single-muon sample with\n"
           "tuned Geant4 physics, and sim tracks from the existing ttbar samples are\n"
           "actively wrong for it -- filtering to 'tracks consistent with a helix'\n"
           "selects on the very null hypothesis the pull test measures and biases pulls\n"
           "low, in a study whose whole purpose is to find out whether the covariance is\n"
           "too tight.  What attaches here: replace build_sequence()'s synthetic hits\n"
           "with the sim hits of one muon and fill_hit_positions() with the real hit\n"
           "positions, keeping run_bkfit() and the per-hit chi2 bookkeeping unchanged;\n"
           "then histogram (fitted - sim) / sqrt(C_ii) at each layer.\n");
  }

  // ==========================================================================
  // Shared helpers for Tasks 1 to 3.
  // ==========================================================================

  namespace {

    // Angle difference, unwrapped -- phi is squashed into (-pi, pi] by
    // squashPhiMPlex() at the end of every propagation, so a raw subtraction
    // can pick up 2 pi and poison a finite difference or a closure residual.
    inline double dphi_unwrap(double a, double b) { return std::remainder(a - b, 2.0 * M_PI); }

    // Normalized covariance difference: |C'_ij - C_ij| / sqrt(C_ii C_jj).
    // Dividing by the geometric mean of the diagonals makes every element
    // dimensionless and comparable, and makes a pure SCALE error on any one
    // element visible at its own size -- which is what the rout/rin jacobian bug
    // recorded in CLAUDE.md was.  An absolute or Frobenius norm hides that,
    // because the six CCS parameters differ by ten orders of magnitude in units.
    double cov_rel_diff(const MPlexLS &a, const MPlexLS &b, int lane, int i, int j) {
      const double d = (double) a.constAt(lane, i, j) - (double) b.constAt(lane, i, j);
      const double s = std::sqrt(std::abs((double) b.constAt(lane, i, i) * (double) b.constAt(lane, j, j)));
      return (s > 0.0) ? std::abs(d) / s : 0.0;
    }

    // Eigenvalues of the CORRELATION matrix, by cyclic Jacobi.  Working on the
    // correlation matrix rather than the covariance makes the six eigenvalues
    // comparable despite the six parameters spanning ten orders of magnitude in
    // units, and puts them on a scale where "how negative is negative" means
    // something (they sum to 6).
    //
    // A CCS covariance out of kalmanOperationPlaneLocal is built as
    // jacLoc2CCS (6x5) * C_loc(5x5) * jacLoc2CCS^T, so it is rank 5 BY
    // CONSTRUCTION and ONE ~zero eigenvalue is expected and correct.  A
    // significantly NEGATIVE eigenvalue is not.
    // Returns false if any diagonal element is non-positive or non-finite.
    bool corr_eigenvalues(const MPlexLS &m, int lane, double ev[6]) {
      double a[6][6];
      for (int i = 0; i < 6; ++i) {
        const double d = m.constAt(lane, i, i);
        if (!fin(d) || d <= 0.0) return false;
      }
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) {
          const double v = m.constAt(lane, i, j);
          if (!fin(v)) return false;
          a[i][j] = v / std::sqrt((double) m.constAt(lane, i, i) * (double) m.constAt(lane, j, j));
        }
      for (int sweep = 0; sweep < 60; ++sweep) {
        double off = 0;
        for (int i = 0; i < 6; ++i) for (int j = i + 1; j < 6; ++j) off += a[i][j] * a[i][j];
        if (off < 1e-24) break;
        for (int i = 0; i < 6; ++i)
          for (int j = i + 1; j < 6; ++j) {
            if (std::abs(a[i][j]) < 1e-30) continue;
            const double th = 0.5 * std::atan2(2.0 * a[i][j], a[i][i] - a[j][j]);
            const double c = std::cos(th), sn = std::sin(th);
            for (int k = 0; k < 6; ++k) {
              const double aik = a[i][k], ajk = a[j][k];
              a[i][k] = c * aik + sn * ajk;
              a[j][k] = -sn * aik + c * ajk;
            }
            for (int k = 0; k < 6; ++k) {
              const double aki = a[k][i], akj = a[k][j];
              a[k][i] = c * aki + sn * akj;
              a[k][j] = -sn * aki + c * akj;
            }
          }
      }
      for (int i = 0; i < 6; ++i) ev[i] = a[i][i];
      std::sort(ev, ev + 6);
      return true;
    }

    // Rotate a vector about a unit axis by angle b (Rodrigues).
    SVector3 rotate_about(const SVector3 &v, const SVector3 &axis, double b) {
      const double c = std::cos(b), s = std::sin(b);
      const double kx = axis[0], ky = axis[1], kz = axis[2];
      const double dot = kx * v[0] + ky * v[1] + kz * v[2];
      const double cx = ky * v[2] - kz * v[1];
      const double cy = kz * v[0] - kx * v[2];
      const double cz = kx * v[1] - ky * v[0];
      SVector3 out;
      out[0] = (float) (v[0] * c + cx * s + kx * dot * (1 - c));
      out[1] = (float) (v[1] * c + cy * s + ky * dot * (1 - c));
      out[2] = (float) (v[2] * c + cz * s + kz * dot * (1 - c));
      return out;
    }

    // A (state, two planes) test case, drawn from the truth helix so that both
    // planes are real modules the track actually crosses.
    struct PairCase {
      SamplePoint sp;
      float par_a[6], cov_a[21];
      SVector3 pnt_a, nrm_a, dir_a;
      SVector3 pnt_b, nrm_b, dir_b;
      double r_a, r_b, z_a, z_b;
      double dalpha;
      int layer_a, layer_b;
      // The crossing BEFORE A, used to condition the input state onto the
      // manifold each propagator preserves.  See the note in pkv_task1_closure().
      float par_o[6], cov_o[21];
    };

    // Build one PairCase per sampled track: A is the crossing nearest the middle
    // of the sequence, B the next one out.  Middle rather than first so that the
    // step is a representative layer-to-layer one and not the seed region.
    void build_pair_cases(const TrackerInfo &ti, std::vector<PairCase> &out, int stride = 1) {
      std::vector<SamplePoint> sample;
      build_sample(sample);
      out.clear();
      for (const SamplePoint &sp : sample) {
        Helix h; h.init(sp.pt, sp.eta, sp.chg, sp.phi0);
        std::vector<Crossing> cx;
        find_crossings(h, ti, cx);
        if ((int) cx.size() < 2 + stride) continue;
        const size_t ia = cx.size() / 2;
        const size_t ib = ia + stride;
        if (ib >= cx.size() || ia < 1) continue;
        PairCase pc;
        pc.sp = sp;
        h.ccs_at_s(cx[ia].s, pc.par_a);
        nominal_cov(pc.par_a, pc.cov_a);
        const ModuleInfo &ma = ti.layer(cx[ia].layer).module_info(cx[ia].module_sid);
        const ModuleInfo &mb = ti.layer(cx[ib].layer).module_info(cx[ib].module_sid);
        pc.pnt_a = ma.pos; pc.nrm_a = ma.zdir; pc.dir_a = ma.xdir;
        pc.pnt_b = mb.pos; pc.nrm_b = mb.zdir; pc.dir_b = mb.xdir;
        pc.r_a = std::hypot(cx[ia].pos[0], cx[ia].pos[1]);
        pc.r_b = std::hypot(cx[ib].pos[0], cx[ib].pos[1]);
        pc.z_a = cx[ia].pos[2];
        pc.z_b = cx[ib].pos[2];
        pc.dalpha = std::abs(h.alpha_of_s(cx[ib].s) - h.alpha_of_s(cx[ia].s));
        pc.layer_a = cx[ia].layer;
        pc.layer_b = cx[ib].layer;
        h.ccs_at_s(cx[ia - 1].s, pc.par_o);
        nominal_cov(pc.par_o, pc.cov_o);
        out.push_back(pc);
      }
    }

    // Per-pT-bin reporting block.
    struct PtBins {
      Stats s[kNPt];
      void add(double pt, double x) { s[pt_bin(pt)].add(x); }
      void print(const char *name) {
        printf("  %-26s", name);
        for (int i = 0; i < kNPt; ++i) printf(" %9.4g", kPtGrid[i]);
        printf("   <- pT [GeV]\n");
        printf("  %-26s", "   median");
        for (int i = 0; i < kNPt; ++i) { if (s[i].n()) printf(" %9.3g", s[i].med()); else printf(" %9s", "-"); }
        printf("\n  %-26s", "   max");
        for (int i = 0; i < kNPt; ++i) { if (s[i].n()) printf(" %9.3g", s[i].max()); else printf(" %9s", "-"); }
        printf("\n");
      }
    };

  }  // anonymous namespace

  // ==========================================================================
  // TASK 1 -- propagation closure, A -> B -> A.
  //
  // Material OFF and uniform field: the propagation is a diffeomorphism, so the
  // round trip must return the input to float precision in the parameters AND
  // the covariance. The covariance is the interesting half -- it is where a
  // jacobian SCALE error lives -- so its residual is reported per element,
  // normalized by sqrt(C_ii C_jj).
  //
  // use_param_b_field is OFF here by necessity: the field is evaluated at each
  // propagation's START point, so A->B uses B(A) and B->A uses B(B), and the
  // round trip cannot close by construction. Its residual then measures the
  // field gradient over the step, not the propagator. Both are reported.
  //
  // A propagate-to-SURFACE constrains its output to that surface, so its 6x6
  // jacobian has rank 5 and the output covariance is degenerate along the
  // normal (Task 3(d) expects the matching zero Cholesky pivot). A round trip
  // from a full-rank 6x6 therefore cannot close however correct the propagator
  // is. Conditioned around it: the input state is itself produced by one
  // propagation onto A's surface, which is what the real code always feeds in.
  // ==========================================================================

  void pkv_task1_closure() {
    printf("\n"
           "============================================================================\n"
           "TASK 1 -- propagation closure, A -> B -> A\n"
           "============================================================================\n");

    const TrackerInfo &ti = tinfo();
    const PropagationFlags base = ti.prop_config().finding_inter_layer_pflags;
    const PropagationFlags pf = pflags_with_param_bfield(pflags_with_material(base, false), false);
    const PropagationFlags pf_pb = pflags_with_param_bfield(pflags_with_material(base, false), true);

    std::vector<PairCase> cases;
    build_pair_cases(ti, cases);
    printf("cases: %zu (one per sampled track; A = middle crossing, B = the next one out)\n",
           cases.size());
    printf("material OFF throughout; uniform B unless stated.\n");

    // ---------------- does the propagation actually land ON the plane? -------
    //
    // Not a closure test, but it belongs here and nothing else measures it.
    // helixAtPlane_impl() solves for the path length with an initial quadratic
    // guess (getS) plus Config::nSStepsInProp2Plane - 1 = ONE refinement
    // iteration, and there is no convergence test.  So the propagated position
    // does not sit exactly on the target plane, and the size of that miss is a
    // hard floor on everything downstream: kalmanOperationPlaneLocal discards
    // the normal component of the propagated position (lxu(2) = 0), so the
    // update silently PROJECTS the state onto the plane -- which is why the
    // hit-covariance-to-infinity limit in Task 3 does not converge to a no-op.
    {
      PtBins d_plane;
      Stats by_dal[kNDal], all;
      for (const PairCase &pc : cases) {
        MPlexQI q;   fill_all_lanes_q(q, pc.sp.chg);
        MPlexHV pntB, nrmB;
        fill_all_lanes_hv(pntB, pc.pnt_b); fill_all_lanes_hv(nrmB, pc.nrm_b);
        MPlexLV p0;  fill_all_lanes(p0, pc.par_a);
        MPlexLS e0;  fill_all_lanes_sym(e0, pc.cov_a);
        MPlexLS e1;  MPlexLV p1;  MPlexQI f1{0};
        propagateHelixToPlaneMPlex(e0, p0, q, pntB, nrmB, nullptr, e1, p1, f1, NN, pf, nullptr);
        const double d = std::abs((p1.constAt(0, 0, 0) - pc.pnt_b[0]) * pc.nrm_b[0] +
                                  (p1.constAt(0, 1, 0) - pc.pnt_b[1]) * pc.nrm_b[1] +
                                  (p1.constAt(0, 2, 0) - pc.pnt_b[2]) * pc.nrm_b[2]);
        d_plane.add(pc.sp.pt, d);
        by_dal[dal_bin(pc.dalpha)].add(d);
        all.add(d);
      }
      printf("\n--- distance from the PROPAGATED position to the target module plane, cm ---\n");
      printf("  all cases: med=%.4g p90=%.4g p99=%.4g max=%.4g\n",
             all.med(), all.pct(90), all.pct(99), all.max());
      d_plane.print("|n.(x_prop - p_mod)| [cm]");
      printf("  by |dalpha| of the step:\n");
      printf("  %-12s %8s %12s %12s %12s\n", "|dalpha|", "n", "median", "p90", "max");
      for (int i = 0; i < kNDal; ++i) {
        if (!by_dal[i].n()) continue;
        char lab[64];
        if (i == kNDal - 1) snprintf(lab, sizeof(lab), ">%.2f", kDalEdges[i]);
        else snprintf(lab, sizeof(lab), "%.2f-%.2f", kDalEdges[i], kDalEdges[i + 1]);
        printf("  %-12s %8zu %12.3g %12.3g %12.3g\n", lab, by_dal[i].n(),
               by_dal[i].med(), by_dal[i].pct(90), by_dal[i].max());
      }
    }

    // ---------------- propagateHelixToPlaneMPlex ----------------
    auto plane_closure = [&](const PropagationFlags &flags, const char *label) {
      PtBins d_pos, d_ipt, d_phi, d_tht, d_cov_max, d_cov_diag;
      Stats by_dal_pos[kNDal], by_dal_cov[kNDal];
      long n_fail = 0;
      for (const PairCase &pc : cases) {
        MPlexQI q;   fill_all_lanes_q(q, pc.sp.chg);
        MPlexHV pntB, nrmB, pntA, nrmA;
        fill_all_lanes_hv(pntB, pc.pnt_b); fill_all_lanes_hv(nrmB, pc.nrm_b);
        fill_all_lanes_hv(pntA, pc.pnt_a); fill_all_lanes_hv(nrmA, pc.nrm_a);

        // Conditioning: O -> plane A, so the input covariance is already
        // degenerate along A's normal, as it always is in the real code.
        MPlexLV po;  fill_all_lanes(po, pc.par_o);
        MPlexLS eo;  fill_all_lanes_sym(eo, pc.cov_o);
        MPlexLS e0;  MPlexLV p0;  MPlexQI f0{0};
        propagateHelixToPlaneMPlex(eo, po, q, pntA, nrmA, nullptr, e0, p0, f0, NN, flags, nullptr);
        if (f0(0, 0, 0)) { ++n_fail; continue; }

        MPlexLS e1, e2;  MPlexLV p1, p2;  MPlexQI f1{0}, f2{0};
        propagateHelixToPlaneMPlex(e0, p0, q, pntB, nrmB, nullptr, e1, p1, f1, NN, flags, nullptr);
        propagateHelixToPlaneMPlex(e1, p1, q, pntA, nrmA, nullptr, e2, p2, f2, NN, flags, nullptr);
        if (f1(0, 0, 0) || f2(0, 0, 0)) ++n_fail;

        const double dp = std::sqrt(std::pow(p2.constAt(0, 0, 0) - p0.constAt(0, 0, 0), 2) +
                                    std::pow(p2.constAt(0, 1, 0) - p0.constAt(0, 1, 0), 2) +
                                    std::pow(p2.constAt(0, 2, 0) - p0.constAt(0, 2, 0), 2));
        d_pos.add(pc.sp.pt, dp);
        d_ipt.add(pc.sp.pt, std::abs(p2.constAt(0, 3, 0) - p0.constAt(0, 3, 0)) / p0.constAt(0, 3, 0));
        d_phi.add(pc.sp.pt, std::abs(dphi_unwrap(p2.constAt(0, 4, 0), p0.constAt(0, 4, 0))));
        d_tht.add(pc.sp.pt, std::abs(p2.constAt(0, 5, 0) - p0.constAt(0, 5, 0)));

        double mx = 0, mdiag = 0;
        for (int i = 0; i < 6; ++i)
          for (int j = 0; j <= i; ++j) {
            const double r = cov_rel_diff(e2, e0, 0, i, j);
            if (fin(r)) {
              mx = std::max(mx, r);
              if (i == j) mdiag = std::max(mdiag, r);
            }
          }
        d_cov_max.add(pc.sp.pt, mx);
        d_cov_diag.add(pc.sp.pt, mdiag);
        by_dal_pos[dal_bin(pc.dalpha)].add(dp);
        by_dal_cov[dal_bin(pc.dalpha)].add(mx);
      }
      printf("\n%s\n", label);
      printf("  prop-fail flags set: %ld / %zu\n", n_fail, cases.size());
      d_pos.print("|d pos| [cm]");
      d_ipt.print("|d(1/pT)|/(1/pT)");
      d_phi.print("|d phi| [rad]");
      d_tht.print("|d theta| [rad]");
      d_cov_max.print("cov: max |dC|/sqrt(CiiCjj)");
      d_cov_diag.print("cov: max on the diagonal");
      printf("  by |dalpha| of the A->B step:\n");
      printf("  %-12s %8s %12s %12s %12s %12s\n", "|dalpha|", "n", "med |dpos|", "max |dpos|", "med cov", "max cov");
      for (int i = 0; i < kNDal; ++i) {
        if (!by_dal_pos[i].n()) continue;
        char lab[64];
        if (i == kNDal - 1) snprintf(lab, sizeof(lab), ">%.2f", kDalEdges[i]);
        else snprintf(lab, sizeof(lab), "%.2f-%.2f", kDalEdges[i], kDalEdges[i + 1]);
        printf("  %-12s %8zu %12.3g %12.3g %12.3g %12.3g\n", lab, by_dal_pos[i].n(),
               by_dal_pos[i].med(), by_dal_pos[i].max(), by_dal_cov[i].med(), by_dal_cov[i].max());
      }
    };

    plane_closure(pf, "--- propagateHelixToPlaneMPlex, uniform B ---");
    plane_closure(pf_pb, "--- propagateHelixToPlaneMPlex, PARAMETRIZED B (cannot close: the\n"
                         "    field is evaluated at the start point of each leg) ---");

    // Per-element covariance closure, so a single bad element cannot hide in a max.
    {
      printf("\n--- propagateHelixToPlaneMPlex: per-element covariance closure, uniform B ---\n");
      printf("    median |dC_ij| / sqrt(C_ii C_jj) over all %zu cases\n", cases.size());
      Stats el[21];
      for (const PairCase &pc : cases) {
        MPlexQI q;   fill_all_lanes_q(q, pc.sp.chg);
        MPlexHV pntB, nrmB, pntA, nrmA;
        fill_all_lanes_hv(pntB, pc.pnt_b); fill_all_lanes_hv(nrmB, pc.nrm_b);
        fill_all_lanes_hv(pntA, pc.pnt_a); fill_all_lanes_hv(nrmA, pc.nrm_a);
        MPlexLV po;  fill_all_lanes(po, pc.par_o);
        MPlexLS eoo; fill_all_lanes_sym(eoo, pc.cov_o);
        MPlexLS e0;  MPlexLV p0;  MPlexQI f0{0};
        propagateHelixToPlaneMPlex(eoo, po, q, pntA, nrmA, nullptr, e0, p0, f0, NN, pf, nullptr);
        MPlexLS e1, e2;  MPlexLV p1, p2;  MPlexQI f1{0}, f2{0};
        propagateHelixToPlaneMPlex(e0, p0, q, pntB, nrmB, nullptr, e1, p1, f1, NN, pf, nullptr);
        propagateHelixToPlaneMPlex(e1, p1, q, pntA, nrmA, nullptr, e2, p2, f2, NN, pf, nullptr);
        for (int i = 0; i < 6; ++i)
          for (int j = 0; j <= i; ++j) el[ls_idx(i, j)].add(cov_rel_diff(e2, e0, 0, i, j));
      }
      const char *pn[6] = {"x", "y", "z", "1/pT", "phi", "theta"};
      for (int i = 0; i < 6; ++i) {
        printf("    ");
        for (int j = 0; j <= i; ++j) printf(" %9.2e", el[ls_idx(i, j)].med());
        printf("   %s\n", pn[i]);
      }
    }

    // ---------------- propagateHelixToRMPlex / ZMPlex ----------------
    {
      PtBins r_pos, r_cov, z_pos, z_cov;
      long nr = 0, nz = 0, nr_fail = 0, nz_fail = 0;
      for (const PairCase &pc : cases) {
        MPlexLV po;  fill_all_lanes(po, pc.par_o);
        MPlexLS eoo; fill_all_lanes_sym(eoo, pc.cov_o);
        MPlexQI q;   fill_all_lanes_q(q, pc.sp.chg);

        if (std::abs(pc.r_b - pc.r_a) > 0.5) {
          MPlexQF rb, ra;
          fill_all_lanes_q(rb, (float) pc.r_b);
          fill_all_lanes_q(ra, (float) pc.r_a);
          MPlexLS e0;  MPlexLV p0;  MPlexQI f0{0};
          propagateHelixToRMPlex(eoo, po, q, ra, e0, p0, f0, NN, pf, nullptr);   // conditioning
          MPlexLS e1, e2;  MPlexLV p1, p2;  MPlexQI f1{0}, f2{0};
          propagateHelixToRMPlex(e0, p0, q, rb, e1, p1, f1, NN, pf, nullptr);
          propagateHelixToRMPlex(e1, p1, q, ra, e2, p2, f2, NN, pf, nullptr);
          if (f0(0, 0, 0) || f1(0, 0, 0) || f2(0, 0, 0)) ++nr_fail;
          r_pos.add(pc.sp.pt, std::sqrt(std::pow(p2.constAt(0, 0, 0) - p0.constAt(0, 0, 0), 2) +
                                        std::pow(p2.constAt(0, 1, 0) - p0.constAt(0, 1, 0), 2) +
                                        std::pow(p2.constAt(0, 2, 0) - p0.constAt(0, 2, 0), 2)));
          double mx = 0;
          for (int i = 0; i < 6; ++i)
            for (int j = 0; j <= i; ++j) { const double r = cov_rel_diff(e2, e0, 0, i, j); if (fin(r)) mx = std::max(mx, r); }
          r_cov.add(pc.sp.pt, mx);
          ++nr;
        }
        if (std::abs(pc.z_b - pc.z_a) > 0.5) {
          MPlexQF zb, za;
          fill_all_lanes_q(zb, (float) pc.z_b);
          fill_all_lanes_q(za, (float) pc.z_a);
          MPlexLS e0;  MPlexLV p0;  MPlexQI f0{0};
          propagateHelixToZMPlex(eoo, po, q, za, e0, p0, f0, NN, pf, nullptr);   // conditioning
          MPlexLS e1, e2;  MPlexLV p1, p2;  MPlexQI f1{0}, f2{0};
          propagateHelixToZMPlex(e0, p0, q, zb, e1, p1, f1, NN, pf, nullptr);
          propagateHelixToZMPlex(e1, p1, q, za, e2, p2, f2, NN, pf, nullptr);
          if (f0(0, 0, 0) || f1(0, 0, 0) || f2(0, 0, 0)) ++nz_fail;
          z_pos.add(pc.sp.pt, std::sqrt(std::pow(p2.constAt(0, 0, 0) - p0.constAt(0, 0, 0), 2) +
                                        std::pow(p2.constAt(0, 1, 0) - p0.constAt(0, 1, 0), 2) +
                                        std::pow(p2.constAt(0, 2, 0) - p0.constAt(0, 2, 0), 2)));
          double mx = 0;
          for (int i = 0; i < 6; ++i)
            for (int j = 0; j <= i; ++j) { const double r = cov_rel_diff(e2, e0, 0, i, j); if (fin(r)) mx = std::max(mx, r); }
          z_cov.add(pc.sp.pt, mx);
          ++nz;
        }
      }
      printf("\n--- propagateHelixToRMPlex, uniform B, r_A -> r_B -> r_A (%ld cases, %ld fail-flagged) ---\n",
             nr, nr_fail);
      r_pos.print("|d pos| [cm]");
      r_cov.print("cov: max |dC|/sqrt(CiiCjj)");
      printf("\n--- propagateHelixToZMPlex, uniform B, z_A -> z_B -> z_A (%ld cases, %ld fail-flagged) ---\n",
             nz, nz_fail);
      z_pos.print("|d pos| [cm]");
      z_cov.print("cov: max |dC|/sqrt(CiiCjj)");
    }

    // ---------------- mini-propagators (no covariance) ----------------
    {
      namespace mp = mini_propagators;
      PtBins mr_pos, mr_mom, mz_pos, mz_mom;
      long nr = 0, nz = 0, nr_fail = 0, nz_fail = 0;
      for (const PairCase &pc : cases) {
        MPlexLV p0;  fill_all_lanes(p0, pc.par_a);
        MPlexQI q;   fill_all_lanes_q(q, pc.sp.chg);
        const mp::InitialStatePlex isp(p0, q);

        if (std::abs(pc.r_b - pc.r_a) > 0.5) {
          MPlexQF rb, ra;
          fill_all_lanes_q(rb, (float) pc.r_b);
          fill_all_lanes_q(ra, (float) pc.r_a);
          mp::StatePlex s1, s2;
          const int ff1 = isp.propagate_to_r(mp::PA_Exact, rb, s1, true, NN);
          const mp::InitialStatePlex isp1(s1, isp);
          const int ff2 = isp1.propagate_to_r(mp::PA_Exact, ra, s2, true, NN);
          if (ff1 || ff2) ++nr_fail;
          mr_pos.add(pc.sp.pt, std::sqrt(std::pow(s2.x[0] - isp.x[0], 2) + std::pow(s2.y[0] - isp.y[0], 2) +
                                         std::pow(s2.z[0] - isp.z[0], 2)));
          mr_mom.add(pc.sp.pt, std::sqrt(std::pow(s2.px[0] - isp.px[0], 2) + std::pow(s2.py[0] - isp.py[0], 2) +
                                         std::pow(s2.pz[0] - isp.pz[0], 2)) / pc.sp.pt);
          ++nr;
        }
        if (std::abs(pc.z_b - pc.z_a) > 0.5) {
          MPlexQF zb, za;
          fill_all_lanes_q(zb, (float) pc.z_b);
          fill_all_lanes_q(za, (float) pc.z_a);
          mp::StatePlex s1, s2;
          const int ff1 = isp.propagate_to_z(mp::PA_Exact, zb, s1, true, NN);
          const mp::InitialStatePlex isp1(s1, isp);
          const int ff2 = isp1.propagate_to_z(mp::PA_Exact, za, s2, true, NN);
          if (ff1 || ff2) ++nz_fail;
          mz_pos.add(pc.sp.pt, std::sqrt(std::pow(s2.x[0] - isp.x[0], 2) + std::pow(s2.y[0] - isp.y[0], 2) +
                                         std::pow(s2.z[0] - isp.z[0], 2)));
          mz_mom.add(pc.sp.pt, std::sqrt(std::pow(s2.px[0] - isp.px[0], 2) + std::pow(s2.py[0] - isp.py[0], 2) +
                                         std::pow(s2.pz[0] - isp.pz[0], 2)) / pc.sp.pt);
          ++nz;
        }
      }
      printf("\n--- mini_propagators::InitialStatePlex::propagate_to_r, round trip (%ld cases, %ld flagged) ---\n",
             nr, nr_fail);
      mr_pos.print("|d pos| [cm]");
      mr_mom.print("|d p| / pT");
      printf("\n--- mini_propagators::InitialStatePlex::propagate_to_z, round trip (%ld cases, %ld flagged) ---\n",
             nz, nz_fail);
      mz_pos.print("|d pos| [cm]");
      mz_mom.print("|d p| / pT");
    }
  }

  // ==========================================================================
  // TASK 2 -- reported covariance transport vs a NUMERICAL jacobian.
  //
  // Perturb each of the six CCS input parameters by +/- delta, re-propagate to
  // the SAME module plane (sPerp = nullptr, so the plane is re-solved every
  // time and the implicit dependence of the path length on the input state is
  // included), and finite-difference.  Then compare J Sigma J^T against the
  // covariance the propagator reports, element by element.
  //
  // delta is chosen by a convergence scan rather than hardcoded: the propagator
  // is single precision, so the usable window is bounded below by cancellation
  // (~eps |f| / delta) and above by truncation (~delta^2 f'''/6).  For each
  // parameter the scan halves delta and keeps the value where successive
  // central-difference estimates agree best; that delta and the residual
  // disagreement are both reported, so the reader can see how much of any
  // mismatch below is numerical.
  // ==========================================================================

  void pkv_task2_numerical_jacobian() {
    printf("\n"
           "============================================================================\n"
           "TASK 2 -- covariance transport vs numerical jacobian\n"
           "============================================================================\n");

    const TrackerInfo &ti = tinfo();
    const PropagationFlags base = ti.prop_config().finding_inter_layer_pflags;
    const PropagationFlags pf = pflags_with_param_bfield(pflags_with_material(base, false), false);

    std::vector<PairCase> cases;
    build_pair_cases(ti, cases);
    printf("cases: %zu ; material OFF, uniform B, plane re-solved for every perturbation\n",
           cases.size());

    // Natural scale per parameter, from which the scanned deltas are built.
    auto scale_of = [](const float par[6], int i) -> double {
      switch (i) {
        case 0: case 1: case 2: return 0.05;                    // cm
        case 3: return 0.02 * par[3];                           // 2% of 1/pT
        default: return 5e-4;                                   // rad
      }
    };

    // One propagation, returning the 6 output parameters (phi unwrapped
    // relative to a reference so the finite difference never sees a 2 pi jump).
    auto prop_once = [&](const PairCase &pc, const float par[6], double out[6]) {
      MPlexLV p;  fill_all_lanes(p, par);
      MPlexLS e;  fill_all_lanes_sym(e, pc.cov_a);
      MPlexQI q;  fill_all_lanes_q(q, pc.sp.chg);
      MPlexHV pnt, nrm;
      fill_all_lanes_hv(pnt, pc.pnt_b); fill_all_lanes_hv(nrm, pc.nrm_b);
      MPlexLS eo;  MPlexLV po;  MPlexQI ff{0};
      propagateHelixToPlaneMPlex(e, p, q, pnt, nrm, nullptr, eo, po, ff, NN, pf, nullptr);
      for (int i = 0; i < 6; ++i) out[i] = po.constAt(0, i, 0);
    };

    const int kNScan = 9;
    Stats scan_pick[6];        // chosen delta / scale, per parameter
    Stats scan_resid[6];       // relative disagreement between the two best estimates
    Stats ratio_diag[6];       // (J S J^T)_ii / reported_ii
    Stats ratio_offd;          // off-diagonal, normalized
    Stats worst_diag;
    PtBins pt_worst_diag;
    long n_cases_used = 0;

    // Also keep a per-element median ratio table.
    Stats el_ratio[21];
    // In-plane (module frame) position covariance -- the rank-safe comparison.
    Stats inplane_ratio[2], inplane_offd;
    PtBins pt_inplane_x, pt_inplane_y;
    // Which direction is each position covariance blind to?
    Stats q_rep_mom, q_rep_nrm, q_num_mom, q_num_nrm, cos_mom_nrm;

    for (const PairCase &pc : cases) {
      // Reference propagation
      double f0[6];
      prop_once(pc, pc.par_a, f0);

      MPlexLV p;  fill_all_lanes(p, pc.par_a);
      MPlexLS e;  fill_all_lanes_sym(e, pc.cov_a);
      MPlexQI q;  fill_all_lanes_q(q, pc.sp.chg);
      MPlexHV pnt, nrm;
      fill_all_lanes_hv(pnt, pc.pnt_b); fill_all_lanes_hv(nrm, pc.nrm_b);
      MPlexLS eo;  MPlexLV po;  MPlexQI ff{0};
      propagateHelixToPlaneMPlex(e, p, q, pnt, nrm, nullptr, eo, po, ff, NN, pf, nullptr);
      if (ff(0, 0, 0)) continue;

      double J[6][6];
      bool bad = false;
      for (int i = 0; i < 6 && !bad; ++i) {
        const double s = scale_of(pc.par_a, i);
        double col[kNScan][6];
        for (int k = 0; k < kNScan; ++k) {
          const double d = s * std::pow(0.5, k);
          float pp[6], pm[6];
          for (int m = 0; m < 6; ++m) { pp[m] = pc.par_a[m]; pm[m] = pc.par_a[m]; }
          pp[i] = (float) (pc.par_a[i] + d);
          pm[i] = (float) (pc.par_a[i] - d);
          double fp[6], fm[6];
          prop_once(pc, pp, fp);
          prop_once(pc, pm, fm);
          for (int m = 0; m < 6; ++m) {
            const double df = (m == 4) ? dphi_unwrap(fp[m], fm[m]) : (fp[m] - fm[m]);
            col[k][m] = df / (2.0 * d);
          }
        }
        // Pick the delta whose column agrees best with the next-finer one,
        // measured relative to the column's own magnitude.
        int best = 0;
        double best_r = 1e30;
        for (int k = 0; k + 1 < kNScan; ++k) {
          double num = 0, den = 0;
          for (int m = 0; m < 6; ++m) {
            const double sc = std::max(std::abs(col[k][m]), std::abs(col[k + 1][m]));
            if (sc <= 0) continue;
            num += std::pow((col[k][m] - col[k + 1][m]) / sc, 2);
            den += 1.0;
          }
          const double r = (den > 0) ? std::sqrt(num / den) : 1e30;
          if (fin(r) && r < best_r) { best_r = r; best = k; }
        }
        if (!fin(best_r) || best_r > 0.5) bad = true;
        scan_pick[i].add(std::pow(0.5, best));
        scan_resid[i].add(best_r);
        for (int m = 0; m < 6; ++m) J[m][i] = col[best][m];
      }
      if (bad) continue;
      ++n_cases_used;

      // C = J Sigma J^T
      double S[6][6], C[6][6];
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) S[i][j] = e.constAt(0, i, j);
      double T[6][6];
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) { double v = 0; for (int k = 0; k < 6; ++k) v += J[i][k] * S[k][j]; T[i][j] = v; }
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) { double v = 0; for (int k = 0; k < 6; ++k) v += T[i][k] * J[j][k]; C[i][j] = v; }

      // The position block of BOTH matrices is degenerate along the plane
      // normal (rank 5, see Task 1), so a raw diagonal ratio in x/y/z is
      // ill-conditioned: it is a ratio of two small differences of larger
      // numbers, and for an endcap disk the z variance is exactly zero.  The
      // well-conditioned comparison is the 2x2 position covariance projected
      // onto the module's own in-plane axes (xdir, ydir) -- which is also
      // exactly the block the Kalman update consumes (psErrLoc(3,3), (3,4),
      // (4,4) after RotateResidualsOnPlane).  Both are reported.
      {
        const SVector3 yd = ROOT::Math::Cross(pc.nrm_b, pc.dir_b);
        const double ax[3] = {pc.dir_b[0], pc.dir_b[1], pc.dir_b[2]};
        const double ay[3] = {yd[0], yd[1], yd[2]};
        double gn[2][2] = {{0, 0}, {0, 0}}, gr[2][2] = {{0, 0}, {0, 0}};
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) {
            const double cn = C[i][j], cr = eo.constAt(0, i, j);
            gn[0][0] += ax[i] * ax[j] * cn;  gr[0][0] += ax[i] * ax[j] * cr;
            gn[0][1] += ax[i] * ay[j] * cn;  gr[0][1] += ax[i] * ay[j] * cr;
            gn[1][1] += ay[i] * ay[j] * cn;  gr[1][1] += ay[i] * ay[j] * cr;
          }
        if (gr[0][0] > 0) inplane_ratio[0].add(gn[0][0] / gr[0][0]);
        if (gr[1][1] > 0) inplane_ratio[1].add(gn[1][1] / gr[1][1]);
        const double nz = std::sqrt(std::abs(gr[0][0] * gr[1][1]));
        if (nz > 0) inplane_offd.add((gn[0][1] - gr[0][1]) / nz);
        if (gr[0][0] > 0) pt_inplane_x.add(pc.sp.pt, gn[0][0] / gr[0][0]);
        if (gr[1][1] > 0) pt_inplane_y.add(pc.sp.pt, gn[1][1] / gr[1][1]);
      }

      // Structure: WHICH DIRECTION is each position covariance blind to?
      //
      // The reported CCS covariance is built as jacCurv2CCS * C_curv * ..., and
      // jacCurv2CCS feeds the three position rows from curvilinear coordinates
      // 3 and 4 only, whose basis vectors are u = (-sinP, cosP, 0) and
      // v = (-cosT cosP, -cosT sinP, sinT) -- the two directions PERPENDICULAR
      // TO THE MOMENTUM.  So the reported position block is degenerate along the
      // MOMENTUM by construction (the standard curvilinear convention).
      //
      // The true covariance of the state that propagate-to-plane produces is
      // degenerate along the PLANE NORMAL instead, because the output position
      // lies on the plane.  Those are different null directions unless the plane
      // happens to be perpendicular to the momentum -- so the two matrices are
      // not even describing the same 2-dimensional subspace, and a
      // direction-by-direction ratio between them is only meaningful where the
      // subspaces overlap.  Measured as normalized quadratic forms below.
      {
        double Mr[3][3], Mn[3][3];
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) { Mr[i][j] = eo.constAt(0, i, j); Mn[i][j] = C[i][j]; }
        const double ph = po.constAt(0, 4, 0), th = po.constAt(0, 5, 0);
        const double m[3] = {std::cos(ph) * std::sin(th), std::sin(ph) * std::sin(th), std::cos(th)};
        const double nr[3] = {pc.nrm_b[0], pc.nrm_b[1], pc.nrm_b[2]};
        auto qform = [](const double M[3][3], const double u[3]) {
          double v = 0;
          for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) v += u[i] * M[i][j] * u[j];
          return v;
        };
        const double tr_r = Mr[0][0] + Mr[1][1] + Mr[2][2];
        const double tr_n = Mn[0][0] + Mn[1][1] + Mn[2][2];
        if (tr_r > 0) {
          q_rep_mom.add(qform(Mr, m) / tr_r);
          q_rep_nrm.add(qform(Mr, nr) / tr_r);
        }
        if (tr_n > 0) {
          q_num_mom.add(qform(Mn, m) / tr_n);
          q_num_nrm.add(qform(Mn, nr) / tr_n);
        }
        double mn = 0;
        for (int i = 0; i < 3; ++i) mn += m[i] * nr[i];
        cos_mom_nrm.add(std::abs(mn));
      }

      double wd = 0;
      for (int i = 0; i < 6; ++i) {
        const double rep = eo.constAt(0, i, i);
        if (rep > 0 && C[i][i] > 0) {
          const double r = C[i][i] / rep;
          ratio_diag[i].add(r);
          wd = std::max(wd, std::abs(r - 1.0));
        }
      }
      worst_diag.add(wd);
      pt_worst_diag.add(pc.sp.pt, wd);
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j <= i; ++j) {
          const double rep = eo.constAt(0, i, j);
          const double nrmz = std::sqrt(std::abs((double) eo.constAt(0, i, i) * (double) eo.constAt(0, j, j)));
          if (nrmz > 0) el_ratio[ls_idx(i, j)].add((C[i][j] - rep) / nrmz);
          if (i != j && nrmz > 0) ratio_offd.add(std::abs((C[i][j] - rep) / nrmz));
        }
    }

    printf("cases used: %ld (a case is dropped if the delta scan does not converge)\n", n_cases_used);
    printf("\n--- delta convergence scan (delta = scale * 2^-k; scale = 0.05 cm, 2%% of 1/pT, 5e-4 rad) ---\n");
    const char *pn[6] = {"x", "y", "z", "1/pT", "phi", "theta"};
    printf("  %-8s %14s %14s %14s\n", "par", "median 2^-k", "median resid", "p90 resid");
    for (int i = 0; i < 6; ++i)
      printf("  %-8s %14.4g %14.3g %14.3g\n", pn[i], scan_pick[i].med(), scan_resid[i].med(), scan_resid[i].pct(90));
    printf("  (resid = relative disagreement between the two best central differences;\n"
           "   it bounds how much of the mismatch below is numerical, not the propagator's)\n");

    printf("\n--- IN-PLANE POSITION COVARIANCE (module frame): (J Sigma J^T) / reported ---\n");
    printf("  This is the rank-safe comparison and the block the Kalman actually uses.\n");
    printf("  %-14s %10s %10s %10s %10s %10s\n", "axis", "n", "median", "p10", "p90", "max");
    printf("  %-14s %10zu %10.5f %10.5f %10.5f %10.4g\n", "xdir (precise)", inplane_ratio[0].n(),
           inplane_ratio[0].med(), inplane_ratio[0].pct(10), inplane_ratio[0].pct(90), inplane_ratio[0].max());
    printf("  %-14s %10zu %10.5f %10.5f %10.5f %10.4g\n", "ydir (strip)", inplane_ratio[1].n(),
           inplane_ratio[1].med(), inplane_ratio[1].pct(10), inplane_ratio[1].pct(90), inplane_ratio[1].max());
    printf("  %-14s %10zu %10.3g %10.3g %10.3g %10.4g   (normalized difference)\n", "xy correlation",
           inplane_offd.n(), inplane_offd.med(), inplane_offd.pct(10), inplane_offd.pct(90), inplane_offd.max());
    printf("  in-plane xdir variance ratio by pT:\n");
    pt_inplane_x.print("  (JSJ^T)/reported, xdir");
    printf("  in-plane ydir variance ratio by pT:\n");
    pt_inplane_y.print("  (JSJ^T)/reported, ydir");

    printf("\n--- DIAGONAL: (J Sigma J^T)_ii / reported_ii, per element ---\n");
    printf("  NOTE x/y/z here are ILL-CONDITIONED (see above); read the momentum rows\n"
           "  (1/pT, phi, theta) and the in-plane table, not these three.\n");
    printf("  %-8s %10s %10s %10s %10s %10s\n", "element", "n", "median", "p10", "p90", "max|r-1|");
    for (int i = 0; i < 6; ++i) {
      if (!ratio_diag[i].n()) continue;
      const double mx = std::max(std::abs(ratio_diag[i].max() - 1.0), std::abs(ratio_diag[i].pct(0) - 1.0));
      printf("  %-8s %10zu %10.5f %10.5f %10.5f %10.3g\n", pn[i], ratio_diag[i].n(),
             ratio_diag[i].med(), ratio_diag[i].pct(10), ratio_diag[i].pct(90), mx);
    }

    printf("\n--- ALL 21 ELEMENTS: median (J Sigma J^T - reported) / sqrt(rep_ii rep_jj) ---\n");
    for (int i = 0; i < 6; ++i) {
      printf("    ");
      for (int j = 0; j <= i; ++j) printf(" %9.2e", el_ratio[ls_idx(i, j)].med());
      printf("   %s\n", pn[i]);
    }
    printf("  off-diagonal |difference| normalized: med=%.3g p90=%.3g max=%.3g\n",
           ratio_offd.med(), ratio_offd.pct(90), ratio_offd.max());

    printf("\n--- WHICH DIRECTION IS EACH POSITION COVARIANCE BLIND TO ---\n");
    printf("  normalized quadratic form u^T M u / trace(M); ~0 means M has no variance along u\n");
    printf("  %-34s %12s %12s %12s\n", "", "median", "p90", "max");
    printf("  %-34s %12.3e %12.3e %12.3e\n", "reported, along MOMENTUM",
           q_rep_mom.med(), q_rep_mom.pct(90), q_rep_mom.max());
    printf("  %-34s %12.3e %12.3e %12.3e\n", "reported, along PLANE NORMAL",
           q_rep_nrm.med(), q_rep_nrm.pct(90), q_rep_nrm.max());
    printf("  %-34s %12.3e %12.3e %12.3e\n", "numerical, along MOMENTUM",
           q_num_mom.med(), q_num_mom.pct(90), q_num_mom.max());
    printf("  %-34s %12.3e %12.3e %12.3e\n", "numerical, along PLANE NORMAL",
           q_num_nrm.med(), q_num_nrm.pct(90), q_num_nrm.max());
    printf("  |cos(momentum, plane normal)| : med=%.4f p10=%.4f  (1 = plane perpendicular to p)\n",
           cos_mom_nrm.med(), cos_mom_nrm.pct(10));
    printf("  If the first row is ~0 and the fourth is ~0 while the second and third are\n"
           "  not, the two matrices are degenerate along DIFFERENT directions -- the\n"
           "  reported one in the curvilinear convention (blind along the momentum), the\n"
           "  numerical one describing a state constrained to the module plane.  The\n"
           "  surface-crossing correction then has to be supplied downstream, which is\n"
           "  what jacCurv2Loc's cosz terms do inside kalmanOperationPlaneLocal.\n");

    printf("\nworst diagonal deviation |ratio - 1| per pT bin:\n");
    pt_worst_diag.print("max|r-1| over diag");
  }

  // ==========================================================================
  // TASK 3 -- analytic limits of kalmanOperationPlaneLocal.
  //
  //  (a) hit covariance -> infinity : the update must be a no-op and chi2 -> 0
  //  (b) hit covariance -> zero     : the updated state must lie ON the hit
  //  (c) rotating plDir about plNrm must leave chi2 and the updated GLOBAL
  //      state unchanged -- the local frame is a choice, not physics
  //  (d) the filtered covariance must be positive semi-definite, and its
  //      diagonal must not exceed the propagated one
  // ==========================================================================

  void pkv_task3_kalman_limits() {
    printf("\n"
           "============================================================================\n"
           "TASK 3 -- Kalman update limits\n"
           "============================================================================\n");

    const TrackerInfo &ti = tinfo();
    const PropagationFlags base = ti.prop_config().finding_inter_layer_pflags;
    const PropagationFlags pf = pflags_with_param_bfield(pflags_with_material(base, false), false);

    std::vector<PairCase> cases;
    build_pair_cases(ti, cases);

    // Build, for each case, a propagated state at plane B and a hit on that
    // plane offset from it by roughly one sigma of the sensor, so that the
    // residual is realistic rather than zero.
    struct Prep {
      MPlexLS P; MPlexLV pp; MPlexQI q;
      MPlexHS msErr; MPlexHV msPar, plNrm, plDir, plPnt;
      SVector3 nrm, dir, ydir, pnt;
      HitRes hr;
      int layer;
      double pt;
      bool ok = false;
    };
    std::vector<Prep> preps;
    preps.reserve(cases.size());
    for (const PairCase &pc : cases) {
      Prep pr;
      pr.layer = pc.layer_b;
      pr.pt = pc.sp.pt;
      MPlexLV p0;  fill_all_lanes(p0, pc.par_a);
      MPlexLS e0;  fill_all_lanes_sym(e0, pc.cov_a);
      fill_all_lanes_q(pr.q, pc.sp.chg);
      fill_all_lanes_hv(pr.plPnt, pc.pnt_b);
      fill_all_lanes_hv(pr.plNrm, pc.nrm_b);
      fill_all_lanes_hv(pr.plDir, pc.dir_b);
      MPlexQI ff{0};
      propagateHelixToPlaneMPlex(e0, p0, pr.q, pr.plPnt, pr.plNrm, nullptr, pr.P, pr.pp, ff, NN, pf, nullptr);
      if (ff(0, 0, 0)) continue;
      pr.nrm = pc.nrm_b; pr.dir = pc.dir_b; pr.pnt = pc.pnt_b;
      pr.ydir = ROOT::Math::Cross(pc.nrm_b, pc.dir_b);
      const LayerInfo &li = ti.layer(pc.layer_b);
      pr.hr = hit_res_for_layer(li);
      ModuleInfo mi(pc.pnt_b, pc.nrm_b, pc.dir_b, 0, 0);
      float cov[6];
      hit_cov_global(mi, pr.hr, cov);
      for (int n = 0; n < NN; ++n) for (int i = 0; i < 6; ++i) pr.msErr.fArray[i * NN + n] = cov[i];
      // Hit = propagated position pushed onto the plane and offset by 1 sigma
      // in each in-plane direction.
      for (int n = 0; n < NN; ++n) {
        const double px = pr.pp.constAt(n, 0, 0), py = pr.pp.constAt(n, 1, 0), pz = pr.pp.constAt(n, 2, 0);
        const double dn = (px - pr.pnt[0]) * pr.nrm[0] + (py - pr.pnt[1]) * pr.nrm[1] + (pz - pr.pnt[2]) * pr.nrm[2];
        const double hx = px - dn * pr.nrm[0] + pr.hr.sig_across * pr.dir[0] + pr.hr.sig_along * pr.ydir[0];
        const double hy = py - dn * pr.nrm[1] + pr.hr.sig_across * pr.dir[1] + pr.hr.sig_along * pr.ydir[1];
        const double hz = pz - dn * pr.nrm[2] + pr.hr.sig_across * pr.dir[2] + pr.hr.sig_along * pr.ydir[2];
        pr.msPar(n, 0, 0) = (float) hx;
        pr.msPar(n, 1, 0) = (float) hy;
        pr.msPar(n, 2, 0) = (float) hz;
      }
      pr.ok = true;
      preps.push_back(pr);
    }
    printf("prepared cases: %zu of %zu\n", preps.size(), cases.size());

    const int kfop = KFO_Calculate_Chi2 | KFO_Update_Params | KFO_Local_Cov;

    // ---------------- (a) hit covariance -> infinity ----------------
    {
      printf("\n--- (a) hit covariance -> infinity: update must be a no-op, chi2 -> 0 ---\n");
      printf("  %-10s %11s %11s %11s %11s %11s %11s\n", "scale", "med chi2", "max chi2",
             "med |dpos|", "p99 |dpos|", "max |dpos|", "n>1um");
      for (double sc : {1.0, 1e4, 1e8, 1e12}) {
        Stats chi2s, dpar, dpos;
        long n_big = 0;
        for (const Prep &pr : preps) {
          MPlexHS me = pr.msErr;
          for (int k = 0; k < me.kTotSize; ++k) me.fArray[k] *= (float) sc;
          MPlexLS C;  MPlexLV up;  MPlexQF chi2{0.0f};
          MPlexQI q = pr.q;
          kalmanOperationPlaneLocal(kfop, pr.P, pr.pp, q, me, pr.msPar, pr.plNrm, pr.plDir, pr.plPnt,
                                    C, up, chi2, NN);
          chi2s.add(chi2(0, 0, 0));
          // Compare the updated state to the propagated one in units of the
          // propagated covariance -- the only scale-free comparison available.
          double mx = 0;
          for (int i = 0; i < 6; ++i) {
            const double sig = std::sqrt(std::abs((double) pr.P.constAt(0, i, i)));
            if (sig <= 0) continue;
            const double d = (i == 4) ? dphi_unwrap(up.constAt(0, i, 0), pr.pp.constAt(0, i, 0))
                                      : (up.constAt(0, i, 0) - pr.pp.constAt(0, i, 0));
            mx = std::max(mx, std::abs(d) / sig);
          }
          dpar.add(mx);
          // Absolute position move, which unlike a sigma ratio cannot be
          // inflated by a near-degenerate direction of the propagated matrix.
          const double dp = std::sqrt(std::pow(up.constAt(0, 0, 0) - pr.pp.constAt(0, 0, 0), 2) +
                                      std::pow(up.constAt(0, 1, 0) - pr.pp.constAt(0, 1, 0), 2) +
                                      std::pow(up.constAt(0, 2, 0) - pr.pp.constAt(0, 2, 0), 2));
          dpos.add(dp);
          if (dp > 1e-4) ++n_big;
        }
        printf("  %-10.0g %11.4g %11.4g %11.4g %11.4g %11.4g %11ld\n", sc, chi2s.med(), chi2s.max(),
               dpos.med(), dpos.pct(99), dpos.max(), n_big);
      }
      printf("  (|dpos| in cm: how far the update moves the state. n>1um counts cases\n"
             "   moving more than 1 um. A sigma-normalized version is not used here: the\n"
             "   propagated position covariance is degenerate along the momentum, so\n"
             "   dividing by sqrt(P_ii) can inflate a harmless move without bound.)\n");
    }

    // ---------------- (b) hit covariance -> zero ----------------
    {
      printf("\n--- (b) hit covariance -> zero: updated state must lie ON the hit ---\n");
      printf("  in-plane residual |x_upd - x_hit| projected on (xdir, ydir), in cm\n");
      printf("  %-10s %14s %14s %14s\n", "scale", "med |res|", "p99 |res|", "max |res|");
      for (double sc : {1.0, 1e-4, 1e-8, 1e-12}) {
        Stats res;
        for (const Prep &pr : preps) {
          MPlexHS me = pr.msErr;
          for (int k = 0; k < me.kTotSize; ++k) me.fArray[k] *= (float) sc;
          MPlexLS C;  MPlexLV up;  MPlexQF chi2{0.0f};
          MPlexQI q = pr.q;
          kalmanOperationPlaneLocal(kfop, pr.P, pr.pp, q, me, pr.msPar, pr.plNrm, pr.plDir, pr.plPnt,
                                    C, up, chi2, NN);
          const double dx = up.constAt(0, 0, 0) - pr.msPar.constAt(0, 0, 0);
          const double dy = up.constAt(0, 1, 0) - pr.msPar.constAt(0, 1, 0);
          const double dz = up.constAt(0, 2, 0) - pr.msPar.constAt(0, 2, 0);
          const double a = dx * pr.dir[0] + dy * pr.dir[1] + dz * pr.dir[2];
          const double b = dx * pr.ydir[0] + dy * pr.ydir[1] + dz * pr.ydir[2];
          res.add(std::hypot(a, b));
        }
        printf("  %-10.0g %14.4g %14.4g %14.4g\n", sc, res.med(), res.pct(99), res.max());
      }
    }

    // ---------------- (c) local-frame rotation invariance ----------------
    {
      printf("\n--- (c) rotate plDir about plNrm: chi2 and the updated GLOBAL state must not move ---\n");
      printf("  %-10s %14s %14s %14s %14s\n", "beta [rad]", "med |dchi2|/chi2", "max |dchi2|/chi2",
             "med |dpos| cm", "max |dpos| cm");
      // Reference
      std::vector<double> chi2_ref(preps.size(), 0.0);
      std::vector<std::array<double, 3>> pos_ref(preps.size());
      for (size_t k = 0; k < preps.size(); ++k) {
        const Prep &pr = preps[k];
        MPlexLS C;  MPlexLV up;  MPlexQF chi2{0.0f};  MPlexQI q = pr.q;
        kalmanOperationPlaneLocal(kfop, pr.P, pr.pp, q, pr.msErr, pr.msPar, pr.plNrm, pr.plDir, pr.plPnt,
                                  C, up, chi2, NN);
        chi2_ref[k] = chi2(0, 0, 0);
        pos_ref[k] = {up.constAt(0, 0, 0), up.constAt(0, 1, 0), up.constAt(0, 2, 0)};
      }
      for (double beta : {0.1, 0.7854, 1.5708, 2.5}) {
        Stats dchi2, dpos;
        for (size_t k = 0; k < preps.size(); ++k) {
          const Prep &pr = preps[k];
          const SVector3 nd = rotate_about(pr.dir, pr.nrm, beta);
          MPlexHV plDir2;  fill_all_lanes_hv(plDir2, nd);
          MPlexLS C;  MPlexLV up;  MPlexQF chi2{0.0f};  MPlexQI q = pr.q;
          kalmanOperationPlaneLocal(kfop, pr.P, pr.pp, q, pr.msErr, pr.msPar, pr.plNrm, plDir2, pr.plPnt,
                                    C, up, chi2, NN);
          if (chi2_ref[k] > 0) dchi2.add(std::abs(chi2(0, 0, 0) - chi2_ref[k]) / chi2_ref[k]);
          dpos.add(std::sqrt(std::pow(up.constAt(0, 0, 0) - pos_ref[k][0], 2) +
                             std::pow(up.constAt(0, 1, 0) - pos_ref[k][1], 2) +
                             std::pow(up.constAt(0, 2, 0) - pos_ref[k][2], 2)));
        }
        printf("  %-10.4f %14.4g %14.4g %14.4g %14.4g\n", beta, dchi2.med(), dchi2.max(),
               dpos.med(), dpos.max());
      }
      printf("  NOTE: msErr is rotated into the local frame inside the function too, so a\n"
             "  rotation of plDir is a pure reparametrization and must cancel exactly.\n");
    }

    // ---------------- (d) positive definiteness, and C <= P on the diagonal ----------------
    {
      printf("\n--- (d) filtered covariance: positive semi-definite, and diag(C) <= diag(P) ---\n");
      long n = 0, n_neg_eig = 0, n_nonfinite = 0, n_chi2_neg = 0, n_bad_diag = 0;
      long n_diag_grew[6] = {0, 0, 0, 0, 0, 0};
      Stats min_eig, second_eig, diag_ratio[6], diag_grow_size[6];
      for (const Prep &pr : preps) {
        MPlexLS C;  MPlexLV up;  MPlexQF chi2{0.0f};  MPlexQI q = pr.q;
        kalmanOperationPlaneLocal(kfop, pr.P, pr.pp, q, pr.msErr, pr.msPar, pr.plNrm, pr.plDir, pr.plPnt,
                                  C, up, chi2, NN);
        ++n;
        if (chi2(0, 0, 0) < 0) ++n_chi2_neg;
        bool nf = false;
        for (int i = 0; i < 6; ++i)
          for (int j = 0; j <= i; ++j) if (!fin(C.constAt(0, i, j))) nf = true;
        if (nf) { ++n_nonfinite; continue; }
        double ev[6];
        if (!corr_eigenvalues(C, 0, ev)) ++n_bad_diag;
        else {
          min_eig.add(ev[0]);
          second_eig.add(ev[1]);
          if (ev[0] < -1e-6) ++n_neg_eig;
        }
        for (int i = 0; i < 6; ++i) {
          const double p = pr.P.constAt(0, i, i), c = C.constAt(0, i, i);
          if (p > 0) {
            diag_ratio[i].add(c / p);
            if (c > p * (1.0 + 1e-5)) { ++n_diag_grew[i]; diag_grow_size[i].add(c / p - 1.0); }
          }
        }
      }
      printf("  cases: %ld ; non-finite covariance: %ld ; non-positive diagonal: %ld ; negative chi2: %ld\n",
             n, n_nonfinite, n_bad_diag, n_chi2_neg);
      printf("  correlation-matrix eigenvalues (they sum to 6; ONE ~0 is expected because\n"
             "  the output is jacLoc2CCS(6x5) * C_loc(5x5) * ..., i.e. rank 5):\n");
      printf("    smallest : med=%.4g  p1=%.4g  min=%.4g\n", min_eig.med(), min_eig.pct(1), min_eig.pct(0));
      printf("    2nd smallest: med=%.4g  p1=%.4g  min=%.4g\n",
             second_eig.med(), second_eig.pct(1), second_eig.pct(0));
      printf("    cases with the smallest eigenvalue < -1e-6 (a genuine defect): %ld / %zu\n",
             n_neg_eig, min_eig.n());
      const char *pn[6] = {"x", "y", "z", "1/pT", "phi", "theta"};
      printf("  %-8s %10s %10s %10s %14s %14s\n", "element", "med C/P", "p90", "max", "n grew", "med growth");
      for (int i = 0; i < 6; ++i)
        printf("  %-8s %10.5f %10.5f %10.4g %14ld %14.3g\n", pn[i], diag_ratio[i].med(),
               diag_ratio[i].pct(90), diag_ratio[i].max(), n_diag_grew[i],
               diag_grow_size[i].n() ? diag_grow_size[i].med() : 0.0);
      printf("  NOTE on 'grew': C is the covariance at the UPDATED parameter values and P\n"
             "  at the propagated ones, and the CCS <-> curvilinear <-> local chain is a\n"
             "  state-dependent reparametrization, so a small diagonal increase is not\n"
             "  strictly a violation of 'an update cannot add information'. The size is\n"
             "  what matters, and it is reported.\n");
    }

    printf("\nSIDE NOTE from reading the code: KFO_Local_Cov has no effect on\n"
           "kalmanOperationPlaneLocal at all -- the only two places it is tested in that\n"
           "function (KalmanUtilsMPlex.cc:1772, :1888) are inside #ifdef DEBUG print\n"
           "blocks. Every caller passes it; none of them gets anything for it.\n");
  }

  // ==========================================================================
  // Entry points
  // ==========================================================================


  // ==========================================================================
  // TASK 5 -- timing A/B of the two getS root forms.
  //
  // Same inputs, same call, same binary: only g_getS_stable_root differs.
  // Inputs are real (state, module plane) pairs from the synthetic crossings, so
  // the mix of incidence angles and step sizes is representative.
  // ==========================================================================
  void pkv_task5_solver_timing() {
    const TrackerInfo &ti = Config::TrkInfo;
    printf("\n============================================================================\n");
    printf("TASK 5 -- getS root-form timing A/B (stable vs original)\n");
    printf("============================================================================\n");

    std::vector<SamplePoint> sample;
    build_sample(sample);

    // Collect representative (state, plane) pairs.
    std::vector<MPlexLV> vpar;
    std::vector<MPlexQI> vchg;
    std::vector<MPlexHV> vpnt, vnrm;
    std::vector<MPlexLS> verr;
    for (const SamplePoint &sp : sample) {
      Helix h;
      h.init(sp.pt, sp.eta, sp.chg, sp.phi0);
      std::vector<SeqHit> seq;
      if (!build_sequence(h, ti, seq)) continue;
      for (size_t ih = 1; ih < seq.size(); ++ih) {
        float par[6], cov[21];
        h.ccs_at_s(seq[ih - 1].s, par);
        std::memset(cov, 0, sizeof(cov));
        for (int i = 0; i < 6; ++i) cov[ls_idx(i, i)] = (i < 3) ? 1e-2f : 1e-4f;
        MPlexLV P; MPlexLS E; MPlexQI Q; MPlexHV PN, NR;
        fill_all_lanes(P, par); fill_all_lanes_sym(E, cov); fill_all_lanes_q(Q, sp.chg);
        fill_all_lanes_hv(PN, seq[ih].plpnt); fill_all_lanes_hv(NR, seq[ih].plnrm);
        vpar.push_back(P); verr.push_back(E); vchg.push_back(Q); vpnt.push_back(PN); vnrm.push_back(NR);
      }
    }
    const int nb = (int) vpar.size();
    const PropagationFlags base = ti.prop_config().finding_inter_layer_pflags;
    const PropagationFlags pf = pflags_with_param_bfield(pflags_with_material(base, false), false);
    printf("batches: %d  (x %d lanes = %d propagations per pass)\n", nb, NN, nb * NN);

    auto bench = [&](bool hermite, int nrep) {
      g_getS_stable_root = hermite;
      // warm up
      for (int b = 0; b < nb; ++b) {
        MPlexLS eo; MPlexLV po; MPlexQI ff{0};
        propagateHelixToPlaneMPlex(verr[b], vpar[b], vchg[b], vpnt[b], vnrm[b], nullptr, eo, po, ff, NN, pf, nullptr);
      }
      auto t0 = std::chrono::steady_clock::now();
      double sink = 0.0;
      for (int r = 0; r < nrep; ++r)
        for (int b = 0; b < nb; ++b) {
          MPlexLS eo; MPlexLV po; MPlexQI ff{0};
          propagateHelixToPlaneMPlex(verr[b], vpar[b], vchg[b], vpnt[b], vnrm[b], nullptr, eo, po, ff, NN, pf, nullptr);
          sink += po(0, 0, 0);
        }
      auto t1 = std::chrono::steady_clock::now();
      const double sec = std::chrono::duration<double>(t1 - t0).count();
      return std::make_pair(sec, sink);
    };

    const int nrep = 200;
    const double nprop = double(nrep) * nb * NN;

    // Two full passes; report the second, so neither solver pays the cold cache.
    double t_getS = 0.0, t_herm = 0.0;
    for (int pass = 0; pass < 2; ++pass) {
      auto [sg, kg] = bench(false, nrep);
      auto [sh, kh] = bench(true, nrep);
      t_getS = sg; t_herm = sh;
      if (pass == 0) continue;
      printf("%-34s %12s %14s %10s\n", "solver", "total [s]", "ns / prop", "relative");
      printf("%-34s %12.4f %14.2f %10.3f\n", "getS quadratic iteration", t_getS, 1e9 * t_getS / nprop, 1.0);
      printf("%-34s %12.4f %14.2f %10.3f\n", "Hermite + Newton + bracket", t_herm, 1e9 * t_herm / nprop,
             t_herm / t_getS);
      printf("(checksums %.6g / %.6g -- differ only by the solve, as intended)\n", kg, kh);
    }
    g_getS_stable_root = true;
  }


  // ==========================================================================
  // TASK 6 -- pT resolution and propagation-failure rate from a forward +
  //           backward fit of REAL SIM TRACKS.
  //
  // Structurally Task 4: same run_bkfit() kernel, SeqHit / HitPosArr / FitOut.
  // Only the hit source changes -- build_sequence()'s synthetic helix-plane
  // intersections give way to a sim track's rec hits. This is the attachment
  // the HOOK note at the end of Task 4 describes.
  //
  // The fit mirrors MkBuilder::fit_tracks() -> MkFitter::fwdFitFitTracks()
  // followed by MkFitter::bkReFitFitTracks():
  //   * forward pass: hits in OUTWARD order (innermost first), starting from
  //     the sim state at the production vertex with a deliberately loose prior,
  //     so the first propagation lands on the innermost module plane;
  //   * backward pass: the same hits in INWARD order, seeded from the forward
  //     pass's end state / covariance / charge.
  // Reported: d(1/pT)/(1/pT) = (ipt_fit - ipt_sim)/ipt_sim at the END of the
  // backward pass, against sim truth at the vertex. Robust statistics only --
  // the distribution has heavy tails, so an RMS is meaningless; width is half
  // the 16-84 percentile spread.
  // ==========================================================================

  namespace {

    const double kT6PtEdges[]  = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 1e9};
    constexpr int kT6NPt = sizeof(kT6PtEdges) / sizeof(kT6PtEdges[0]) - 1;
    const double kT6EtaEdges[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 1e9};
    constexpr int kT6NEta = sizeof(kT6EtaEdges) / sizeof(kT6EtaEdges[0]) - 1;

    int t6_pt_bin(double pt) {
      for (int i = 0; i < kT6NPt; ++i) if (pt < kT6PtEdges[i + 1]) return i;
      return kT6NPt - 1;
    }
    int t6_eta_bin(double e) {
      e = std::abs(e);
      for (int i = 0; i < kT6NEta; ++i) if (e < kT6EtaEdges[i + 1]) return i;
      return kT6NEta - 1;
    }

    // Minimum number of valid rec hits on the sim track.  8 is the maintainer's
    // suggested value and is also roughly what a phase-2 track needs before a
    // 5-parameter fit has any redundancy at all.
    constexpr int kT6MinHits = 8;
    constexpr double kT6MinPt = 0.5;

    // --- SIM-TRACK SHAPING ---------------------------------------------
    //
    // The raw hit list of a seed-labelled sim track is a zoo: loopers come
    // back inward, back-scattered tracks re-enter layers they already
    // crossed, and a few hits are simply nowhere near the trajectory.  Any
    // ordering by radius is then WRONG, and the fit diverges (chi2 1e14).
    // The previous stop-gap -- reject any track with more than two rec hits
    // in one mkFit layer -- threw away ~50% of the sample.
    //
    // The shaping below replaces it with three stages, all keyed on the
    // IDEAL HELIX of the sim track's own production state (uniform B, no
    // material).  Everything it computes is a function of a hit POSITION
    // alone, so no stage depends on the list order we are trying to
    // establish in the first place:
    //
    //   centre  C  = (x0 - k py0, y0 + k px0),  R_c = |k| pT
    //   v(a) = P(a) - C = R(+a) v(0)   -- the turn angle 'a' of a hit is the
    //                                     CCW angle from v(0) to its own v.
    //   m(a) = (1/k) R90 v(a)          -- the transverse MOMENTUM direction
    //                                     at a hit follows from its position.
    //
    // Stage 1  order by turn angle 'a', i.e. by path length, not by radius.
    // Stage 2  truncate at the apex, using the production criterion
    //          (MkBuilder.cc:693-703): stop once the transverse angle
    //          between position and momentum reaches pi/2 - 0.2 rad, i.e.
    //          cos(momPhi - posPhi) < sin(0.2).  Dot-product form, no wrap.
    // Stage 3  drop severe outliers by an iterative per-hit chi2 from a
    //          forward-only first pass.
    //
    // Shaping is bit-selectable so the before/after comparison is one flag.
    enum { kT6ShOrder = 1, kT6ShApex = 2, kT6ShOutlier = 4 };

    const double kT6SinApex = std::sin(0.2);   // 0.19867, = cos(78.54 deg)
    // Defaults chosen by a scan (see CLAUDE.md, Task 6 shaping): 25 is 4.5
    // sigma per coordinate at 2 d.o.f. and sits on a plateau -- 25 and 100
    // give the same resolution -- while 10 removals rather than 3 is what
    // actually matters, because a track can carry several stray hits.
    double kT6OutlierChi2 = 25.0;
    int    kT6MaxOutliers = 10;

    // Ideal helix of a sim track, from its production state.
    struct T6Helix {
      double cx = 0, cy = 0, z0 = 0, kk = 1, rc = 1, kpz = 0, v0x = 1, v0y = 0;
      bool ok = false;

      void init(const Track &st) {
        const auto &pp = st.state().parameters;
        const double x0 = pp[0], y0 = pp[1];
        z0 = pp[2];
        const double pt  = (pp[3] != 0.0f) ? std::abs(1.0 / (double) pp[3]) : 0.0;
        const double phi = pp[4], th = pp[5];
        const double px0 = pt * std::cos(phi), py0 = pt * std::sin(phi);
        const double pz0 = pt / std::tan(th);
        const double inv_k = ((st.charge() < 0) ? 0.01 : -0.01) *
                             (double) Const::sol * (double) Config::Bfield;
        kk  = 1.0 / inv_k;
        cx  = x0 - kk * py0;
        cy  = y0 + kk * px0;
        v0x = x0 - cx;  v0y = y0 - cy;
        rc  = std::abs(kk) * pt;
        kpz = std::abs(kk) * pz0;   // z = z0 + kpz * alpha, alpha as returned below
        ok  = fin((float) rc) && rc > 1e-3 && fin((float) cx) && fin((float) cy) &&
              fin((float) kpz) && (std::hypot(v0x, v0y) > 1e-6);
      }

      // Turn angle of a hit, measured from the production point and folded
      // with sign(k) so that it INCREASES along the direction of travel for
      // both charges (v(a) = R(+a) v(0) with a = s/(k |p|), and k is signed).
      // Range [0, 2pi): a track turning by more than a full circle wraps, and
      // such a hit then lands early in the sorted list -- where the apex test
      // truncates it, which is the conservative direction.
      double alpha(const double p[3]) const {
        const double vx = p[0] - cx, vy = p[1] - cy;
        double a = std::atan2(v0x * vy - v0y * vx, v0x * vx + v0y * vy) * ((kk > 0) ? 1.0 : -1.0);
        if (a < -0.1) a += 2.0 * M_PI;   // hits are at a >= 0; only noise is negative
        return a;
      }
      // cos of the transverse angle between the momentum and position vectors.
      // 1 = heading straight out, 0 = apex, < 0 = coming back in.
      double cos_apex(const double p[3]) const {
        const double vx = p[0] - cx, vy = p[1] - cy;
        const double vn = std::hypot(vx, vy), r = std::hypot(p[0], p[1]);
        if (vn < 1e-9 || r < 1e-9) return 1.0;
        const double sg = (kk > 0) ? 1.0 : -1.0;
        return (sg * (-vy) * p[0] + sg * vx * p[1]) / (vn * r);
      }
      double dcirc(const double p[3]) const { return std::hypot(p[0] - cx, p[1] - cy) - rc; }
      double dz(const double p[3], double a) const { return p[2] - (z0 + kpz * a); }  // kpz already |k| pz0
    };

    // Per-hit shaping metadata, parallel to the SeqHit vector.
    struct T6HitInfo { double alpha = 0, cosap = 1, dcirc = 0, dz = 0, r3 = 0; };

    struct T6Acc {
      Stats res[kT6NPt][kT6NEta];          // d(1/pT)/(1/pT), fail-flagged EXCLUDED
      long  n_ok[kT6NPt][kT6NEta];         // fits entering res
      long  n_fail[kT6NPt][kT6NEta];       // fits with a propagation failure
      Stats res_all;
      Stats nhits, chi2ndf_fwd, chi2ndf_bkw;
      Stats hit_plane_dist, hit_mod_dist;
      Stats res_loose;                     // prior-sensitivity cross-check

      long n_events = 0;
      long n_sim = 0, n_seeds = 0, n_seed_lab = 0, n_uniq_lab = 0;
      long n_after_hits = 0, n_after_loop = 0, n_after_pt = 0, n_fitted = 0;
      long n_fail_tot = 0, n_nonfinite = 0, n_gross = 0;
      long n_fail_fwd = 0, n_fail_bkw = 0, n_fail_firststep = 0;
      long n_rej_multihit = 0;
      long maxperlayer[16] = {};
      long n_badstate = 0, n_outlier = 0;
      std::vector<std::array<int, 3>> badcases;   // (event, label, category)
      long faillayer[64] = {};

      // --- RAW per-fit / per-hit values, kept so the SHAPE of the
      // distributions can be written out as ROOT histograms.  The text report
      // is all medians and robust widths, which is exactly what hides a tail.
      // Parallel vectors, one entry per CLEAN fit (v_hitchi2 is per hit, over
      // every fit that ran, clean or not).
      std::vector<float> v_res, v_pt, v_eta, v_c2f, v_c2b;
      std::vector<int> v_evt, v_lbl;      // identity, so two snapshots can be PAIRED
      std::vector<float> v_hitchi2;

      // --- shaping bookkeeping -----------------------------------------
      long n_after_apex = 0, n_after_outl = 0;
      long n_tr_truncated = 0, n_hits_apex = 0;      // apex stage
      long n_tr_outl = 0, n_hits_outl = 0;           // outlier stage
      long n_lost_apex = 0, n_lost_outl = 0;         // dropped below kT6MinHits
      long n_nohelix = 0;                            // sim state unusable
      long n_order_swaps = 0, n_order_tracks = 0;    // alpha vs radius ordering
      long maxperlayer_shaped[16] = {};
      Stats st_alpha_span, st_dcirc, st_dz, st_cosap, st_outchi2;
      // The subset that the OLD blunt veto (<=2 raw hits per layer) also kept,
      // so shaped and unshaped numbers can be compared on identical tracks.
      Stats res_legacy, c2_legacy;
      long n_leg_fit = 0, n_leg_bad = 0;
      // Multi-hit-layer anatomy: are 3+ hits in one layer a turn-around, or a
      // genuine multi-module crossing, or duplicate/stray hits?
      Stats st_mhl_aspan, st_mhl_dcirc;
      long n_mhl = 0;
      long n_apex_by_ptbin[kT6NPt] = {};
      long n_tot_by_ptbin[kT6NPt] = {};

      void reset() { *this = T6Acc(); }
    };

    T6Acc g_t6;
    int g_t6_debug = 0;
    int g_t6_maxhits_per_layer = 2;
    // Default: apex + outlier.  The turn-angle ORDERING bit is off, and that
    // is a measured choice, not an oversight: once the apex truncation has
    // removed the turn-arounds the track is monotonically outgoing, so 3D
    // radius IS the path order and is numerically the more robust of the two
    // (pure turn-angle ordering scrambles near-degenerate multi-module
    // crossings in the forward disks -- measured chi2/ndf p90 29.3 against
    // 4.6).  Bit 1 is kept so the comparison can be re-made.
    int g_t6_shaping = kT6ShApex | kT6ShOutlier;
    int g_t6_dbgbuild = 0;
    // Turn-angle bucket for the ordering: hits whose alpha differs by less than
    // this are treated as simultaneous and ordered by 3D radius instead.  Pure
    // alpha ordering (eps = 0) scrambles the near-degenerate multi-module
    // crossings that dominate the forward disks.
    double g_t6_order_eps = 0.05;
    long g_t6_badmod = 0, g_t6_zeronrm = 0;

    // Build the SeqHit sequence of one sim track out of its rec hits, in
    // OUTWARD order (sorted by 3D radius -- monotone for a track leaving the
    // interaction region, and it does not assume the file's hit ordering).
    bool t6_build_hits(const Event &ev, const TrackerInfo &ti, const Track &st,
                       const T6Helix &hx, std::vector<SeqHit> &out,
                       std::vector<T6HitInfo> &info) {
      struct Ent { double r3; SeqHit sh; T6HitInfo hi; };
      std::vector<Ent> v;
      const int nl = (int) ev.layerHits_.size();
      for (int i = 0; i < st.nTotalHits(); ++i) {
        const HitOnTrack hot = st.getHitOnTrack(i);
        if (hot.index < 0 || hot.layer < 0 || hot.layer >= nl) continue;
        if (hot.index >= (int) ev.layerHits_[hot.layer].size()) continue;
        if (hot.layer >= ti.n_layers()) continue;
        const Hit &h = ev.layerHits_[hot.layer][hot.index];
        const LayerInfo &li = ti.layer(hot.layer);
        const unsigned int mid = h.detIDinLayer();
        if ((int) mid >= li.n_modules()) { ++g_t6_badmod; continue; }
        Ent e;
        e.sh.layer = hot.layer;
        e.sh.sid = (int) mid;
        const ModuleInfo &mi = li.module_info(mid);
        e.sh.plpnt = mi.pos;
        e.sh.plnrm = mi.zdir;
        e.sh.pldir = mi.xdir;
        e.sh.plydir = mi.calc_ydir();
        bool ok = true;
        for (int k = 0; k < 3; ++k) {
          const float p = h.posArray()[k];
          if (!fin(p)) ok = false;
          e.sh.pos_true[k] = p;
        }
        for (int k = 0; k < 6; ++k) {
          const float c = h.errArray()[k];
          if (!fin(c)) ok = false;
          e.sh.cov[k] = c;
        }
        if (!ok) continue;
        {
          const double nn2 = mi.zdir[0]*mi.zdir[0] + mi.zdir[1]*mi.zdir[1] + mi.zdir[2]*mi.zdir[2];
          if (nn2 < 0.25) { ++g_t6_zeronrm; continue; }
        }
        {
          double d = 0, dt = 0;
          for (int k = 0; k < 3; ++k) d += mi.zdir[k] * (e.sh.pos_true[k] - mi.pos[k]);
          for (int k = 0; k < 3; ++k) { const double t = e.sh.pos_true[k] - mi.pos[k]; dt += t * t; }
          g_t6.hit_plane_dist.add(std::abs(d));
          g_t6.hit_mod_dist.add(std::sqrt(dt));
        }
        e.r3 = std::sqrt(e.sh.pos_true[0] * e.sh.pos_true[0] +
                         e.sh.pos_true[1] * e.sh.pos_true[1] +
                         e.sh.pos_true[2] * e.sh.pos_true[2]);
        e.hi.r3 = e.r3;
        if (hx.ok) {
          e.hi.alpha = hx.alpha(e.sh.pos_true);
          e.hi.cosap = hx.cos_apex(e.sh.pos_true);
          e.hi.dcirc = hx.dcirc(e.sh.pos_true);
          e.hi.dz    = hx.dz(e.sh.pos_true, e.hi.alpha);
        }
        v.push_back(e);
      }
      // Stage 1: order by turn angle (== path length) rather than by radius.
      // For a track that turns around, radius is not monotone along the
      // trajectory and the resulting sequence is simply not the order the
      // hits were crossed in.  Turn angle is, and it is order-free.
      const bool by_alpha = hx.ok && (g_t6_shaping & kT6ShOrder);
      if (by_alpha) {
        const double eps = g_t6_order_eps;
        if (eps > 0.0)
          std::sort(v.begin(), v.end(), [eps](const Ent &a, const Ent &b) {
            const long ba = (long) std::floor(a.hi.alpha / eps), bb = (long) std::floor(b.hi.alpha / eps);
            return (ba != bb) ? (ba < bb) : (a.r3 < b.r3);
          });
        else
          std::sort(v.begin(), v.end(), [](const Ent &a, const Ent &b) { return a.hi.alpha < b.hi.alpha; });
      }
      else
        std::sort(v.begin(), v.end(), [](const Ent &a, const Ent &b) { return a.r3 < b.r3; });
      // How often do the two orderings disagree?  (Adjacent-transposition
      // count between the alpha order and the radius order.)
      if (hx.ok && v.size() > 1) {
        std::vector<int> byr(v.size());
        for (size_t i = 0; i < v.size(); ++i) byr[i] = (int) i;
        std::sort(byr.begin(), byr.end(), [&](int a, int b) { return v[a].r3 < v[b].r3; });
        long sw = 0;
        for (size_t i = 0; i < byr.size(); ++i) for (size_t j = i + 1; j < byr.size(); ++j)
          if (byr[i] > byr[j]) ++sw;
        g_t6.n_order_swaps += sw;
        ++g_t6.n_order_tracks;
      }
      out.clear();  info.clear();
      out.reserve(v.size());  info.reserve(v.size());
      for (const Ent &e : v) { out.push_back(e.sh); info.push_back(e.hi); }
      if (g_t6_dbgbuild > 0 && v.size() >= 8) {
        --g_t6_dbgbuild;
        printf("T6BUILD pt=%.3g eta=%.3g chg=%d nh=%d  C=(%.4g,%.4g) Rc=%.4g k=%.4g z0=%.4g kpz=%.4g\n",
               st.pT(), st.momEta(), st.charge(), (int) v.size(), hx.cx, hx.cy, hx.rc, hx.kk, hx.z0, hx.kpz);
        printf("   vtx=(%.4g,%.4g,%.4g) phi=%.4g th=%.4g\n", st.state().parameters[0],
               st.state().parameters[1], st.state().parameters[2],
               st.state().parameters[4], st.state().parameters[5]);
        for (size_t i = 0; i < v.size(); ++i)
          printf("   L%-3d r3=%8.3f  a=%8.4f cos=%7.4f dcirc=%9.3g dz=%10.4g  z=%8.3f\n",
                 v[i].sh.layer, v[i].r3, v[i].hi.alpha, v[i].hi.cosap, v[i].hi.dcirc, v[i].hi.dz,
                 v[i].sh.pos_true[2]);
      }
      return true;
    }

    void t6_prior(const float par[6], float cov[21], float ipt_rel) {
      std::memset(cov, 0, 21 * sizeof(float));
      const float sp = 1.0f, sa = 0.05f, si = ipt_rel * par[3];
      cov[ls_idx(0, 0)] = sp * sp;
      cov[ls_idx(1, 1)] = sp * sp;
      cov[ls_idx(2, 2)] = sp * sp;
      cov[ls_idx(3, 3)] = si * si;
      cov[ls_idx(4, 4)] = sa * sa;
      cov[ls_idx(5, 5)] = sa * sa;
    }

    // Outcome of one forward+backward fit.
    struct T6Out {
      double ipt = 0;
      bool propfail = false;     // a flagged propagation OTHER than the structural first backward step
      bool first_step = false;   // only the structural first backward step was flagged
      bool nonfinite = false;    // non-finite end state
      bool badstate = false;     // finite but non-physical end state
      double chi2ndf_fwd = 0, chi2ndf_bkw = 0;
      std::vector<float> hitchi2;   // per-hit chi2, forward pass then backward
    };

    // Returns the fitted 1/pT; sets fail if any propagation flagged or the
    // result is not finite.  nh_used is the hit count actually fitted.
    double t6_fwd_then_bkw(const TrackerInfo &ti, const Track &st,
                           const std::vector<SeqHit> &outward, const HitPosArr &hp_out,
                           const PropagationFlags &pf, MatMode mm, float ipt_rel,
                           T6Out &res) {
      const int nh = (int) outward.size();
      float par0[6], cov0[21];
      for (int i = 0; i < 6; ++i) par0[i] = st.state().parameters[i];
      t6_prior(par0, cov0, ipt_rel);

      MPlexLV P;  fill_all_lanes(P, par0);
      MPlexLS E;  fill_all_lanes_sym(E, cov0);
      MPlexQI Q;  fill_all_lanes_q(Q, st.charge());

      FitOut fo_f;
      run_bkfit(ti, outward, hp_out, P, E, Q, mm, pf, fo_f);

      // Backward: same hits, reversed, seeded from the forward end state.
      //
      // The OUTERMOST hit is dropped from the backward list.  The forward pass
      // ends with its state ON that module's plane, so propagating to it again
      // is a zero-length step: the Hermite span collapses, sl = 0, and the sign
      // bracket in helixAtPlane_impl reports a propagation failure -- correctly,
      // there is nothing to solve.  Measured before this was done: it fired on
      // essentially every track and made the FailFlag rate uninformative.
      std::vector<SeqHit> inward(outward.rbegin(), outward.rend() - 1);
      HitPosArr hp_in(hp_out.rbegin(), hp_out.rend() - 1);
      FitOut fo_b;
      run_bkfit(ti, inward, hp_in, fo_f.par_end, fo_f.err_end, fo_f.chg_end, mm, pf, fo_b);

      res.hitchi2.clear();
      res.hitchi2.reserve(fo_f.chi2_hit.size() + fo_b.chi2_hit.size());
      for (const auto &c : fo_f.chi2_hit) res.hitchi2.push_back(c(0, 0, 0));
      for (size_t q = 1; q < fo_b.chi2_hit.size(); ++q)   // skip the structural first backward step
        res.hitchi2.push_back(fo_b.chi2_hit[q](0, 0, 0));

      const double c2f = fo_f.chi2_total(0, 0, 0);
      const double c2b = fo_b.chi2_total(0, 0, 0);
      res.chi2ndf_fwd = c2f / (2.0 * nh);
      res.chi2ndf_bkw = c2b / (2.0 * std::max(1, (int) inward.size()));

      // The FIRST backward propagation is structurally degenerate whenever the
      // two outermost hits are a stereo pair / overlap: the forward pass ends
      // on one module and the target plane is ~2 mm away, so the Hermite span
      // collapses and the sign bracket reports a failure.  It is real (the
      // production bkReFit has it too) but it is not a statement about the
      // track, so it is counted separately.
      int n_real = fo_f.n_prop_fail;
      for (int q : fo_b.fail_hits) { if (q == 0) res.first_step = true; else ++n_real; }
      res.propfail = (n_real > 0);
      g_t6.n_fail_fwd += (fo_f.n_prop_fail > 0);
      g_t6.n_fail_bkw += (n_real > fo_f.n_prop_fail);
      g_t6.n_fail_firststep += res.first_step;
      for (int q : fo_f.fail_hits) g_t6.faillayer[std::min(63, outward[q].layer)]++;
      for (int q : fo_b.fail_hits) if (q != 0) g_t6.faillayer[std::min(63, inward[q].layer)]++;

      if (g_t6_debug > 0) {
        --g_t6_debug;
        printf("T6DBG nh=%d fwd_fail(%d):", nh, fo_f.n_prop_fail);
        for (int q : fo_f.fail_hits) printf(" %d/L%d", q, outward[q].layer);
        printf("  bkw_fail(%d):", fo_b.n_prop_fail);
        for (int q : fo_b.fail_hits) printf(" %d/L%d", q, inward[q].layer);
        printf("  c2f=%.4g c2b=%.4g ipt=%.6g\n", res.chi2ndf_fwd, res.chi2ndf_bkw, fo_b.par_end(0,3,0));
        printf("   fwd chi2/hit:");
        for (int q = 0; q < nh; ++q) printf(" %.3g", fo_f.chi2_hit[q](0,0,0));
        printf("\n");
      }
      const double ipt = fo_b.par_end(0, 3, 0);
      res.ipt = ipt;
      for (int i = 0; i < 6; ++i) if (!fin(fo_b.par_end(0, i, 0))) res.nonfinite = true;
      if (!fin(ipt)) res.nonfinite = true;
      if (!res.nonfinite) {
        const double rr = std::hypot(fo_b.par_end(0, 0, 0), fo_b.par_end(0, 1, 0));
        const double zz = std::abs(fo_b.par_end(0, 2, 0));
        const double pt = (ipt != 0.0) ? std::abs(1.0 / ipt) : 1e30;
        // Non-physical end state: outside the tracker volume, or a pT the
        // finding loop would have thrown away (minPtCut-equivalent), or one so
        // stiff it is meaningless.
        if (rr > 150.0 || zz > 350.0 || pt < 0.1 || pt > 1.0e4) res.badstate = true;
      }
      return ipt;
    }

    // ------------------------------------------------------------------
    // Stage 2: truncate at the apex.
    //
    // Production criterion, MkBuilder::find_tracks_unroll_candidates()
    // (MkBuilder.cc:693-703): stop a candidate once the transverse angle
    // between its position and momentum vectors reaches pi/2 - 0.2 rad
    // (78.54 deg).  That file spells it as |posPhi - momPhi| against a pair
    // of bounds whose upper one has to be DERIVED as TwoPI - kMaxAngPosMom;
    // here the dot-product form is used instead, which is the same test with
    // no branch cut anywhere: cos(momPhi - posPhi) < sin(0.2).
    //
    // The momentum direction is the sim track's IDEAL-HELIX momentum at the
    // hit's own position (not a chord between consecutive hits, which would
    // depend on the ordering this stage exists to make trustworthy, and not
    // a per-hit sim momentum, which the ntuple does not carry -- MCHitInfo
    // has only the track/layer ids).
    //
    // Unlike production this is NOT gated on pT < 1.2 GeV or r > 25 cm: with
    // sim truth the test simply cannot fire on a stiff track, so the gates
    // would only hide real turn-arounds at small radius.
    int t6_truncate_at_apex(std::vector<SeqHit> &seq, std::vector<T6HitInfo> &info) {
      const int n0 = (int) seq.size();
      int keep = n0;
      for (int i = 0; i < n0; ++i)
        if (info[i].cosap < kT6SinApex) { keep = i; break; }
      if (keep < n0) {
        seq.resize(keep);
        info.resize(keep);
      }
      return n0 - keep;
    }

    // ------------------------------------------------------------------
    // Stage 3: severe-outlier removal by iterative per-hit chi2.
    //
    // A forward-only first pass over the surviving hits gives a per-hit chi2
    // (2 d.o.f.) from the very machinery the fit uses, material model
    // included -- which a geometric residual to the ideal helix would not,
    // and which matters because energy loss makes a low-pT track curl inside
    // its own vertex helix by more than any fixed geometric tolerance.
    // Drop the single worst hit if it exceeds kT6OutlierChi2, refit, repeat.
    // One at a time, because one wrecked hit poisons every chi2 after it.
    int t6_remove_outliers(const TrackerInfo &ti, const Track &st,
                           std::vector<SeqHit> &seq, std::vector<T6HitInfo> &info,
                           const PropagationFlags &pf, MatMode mm,
                           std::mt19937_64 &rng) {
      int nrem = 0;
      for (int it = 0; it < kT6MaxOutliers; ++it) {
        const int nh = (int) seq.size();
        if (nh <= kT6MinHits) break;
        HitPosArr hp;
        fill_hit_positions(seq, false, rng, hp);
        float par0[6], cov0[21];
        for (int i = 0; i < 6; ++i) par0[i] = st.state().parameters[i];
        t6_prior(par0, cov0, 0.5f);
        MPlexLV P;  fill_all_lanes(P, par0);
        MPlexLS E;  fill_all_lanes_sym(E, cov0);
        MPlexQI Q;  fill_all_lanes_q(Q, st.charge());
        FitOut fo;
        run_bkfit(ti, seq, hp, P, E, Q, mm, pf, fo);
        int iw = -1;
        double cw = kT6OutlierChi2;
        for (int ih = 0; ih < nh; ++ih) {
          const double c = fo.chi2_hit[ih](0, 0, 0);
          if (!fin(c)) { iw = ih; cw = 1e30; break; }   // non-finite is the worst kind
          if (c > cw) { cw = c; iw = ih; }
        }
        if (iw < 0) break;
        g_t6.st_outchi2.add(cw);
        seq.erase(seq.begin() + iw);
        info.erase(info.begin() + iw);
        ++nrem;
      }
      return nrem;
    }

    void t6_print_grid(const char *title, const char *what,
                       double (*get)(int, int), bool as_pct) {
      printf("\n  %s\n", title);
      printf("    %-14s", what);
      for (int ie = 0; ie < kT6NEta; ++ie) {
        if (kT6EtaEdges[ie + 1] > 100) printf("  |eta|>%.1f", kT6EtaEdges[ie]);
        else printf("  %.1f-%.1f ", kT6EtaEdges[ie], kT6EtaEdges[ie + 1]);
      }
      printf("\n");
      for (int ip = 0; ip < kT6NPt; ++ip) {
        char lab[32];
        if (kT6PtEdges[ip + 1] > 1e8) snprintf(lab, sizeof(lab), "pT>%.0f", kT6PtEdges[ip]);
        else snprintf(lab, sizeof(lab), "pT %.3g-%.3g", kT6PtEdges[ip], kT6PtEdges[ip + 1]);
        printf("    %-14s", lab);
        for (int ie = 0; ie < kT6NEta; ++ie) {
          const double v = get(ip, ie);
          if (v == -1e30) printf(" %9s", "-");
          else if (as_pct) printf(" %8.2f%%", 100.0 * v);
          else printf(" %9.3g", v);
        }
        printf("\n");
      }
    }

    double t6_get_n(int ip, int ie)    { return (double) g_t6.n_ok[ip][ie]; }
    double t6_get_med(int ip, int ie)  { return g_t6.n_ok[ip][ie] ? g_t6.res[ip][ie].med() : -1e30; }
    double t6_get_wid(int ip, int ie)  {
      if (g_t6.n_ok[ip][ie] < 4) return -1e30;
      return 0.5 * (g_t6.res[ip][ie].pct(84.0) - g_t6.res[ip][ie].pct(16.0));
    }
    double t6_get_fail(int ip, int ie) {
      const long tot = g_t6.n_ok[ip][ie] + g_t6.n_fail[ip][ie];
      return tot ? (double) g_t6.n_fail[ip][ie] / tot : -1e30;
    }
    double t6_get_tot(int ip, int ie)  { return (double) (g_t6.n_ok[ip][ie] + g_t6.n_fail[ip][ie]); }

  }  // anonymous namespace

  void pkv_task6_reset() { g_t6.reset(); }
  void pkv_task6_debug(int n) { g_t6_debug = n; }
  void pkv_task6_set_max_hits_per_layer(int n) { g_t6_maxhits_per_layer = n; }
  void pkv_task6_set_shaping(int m) { g_t6_shaping = m; }
  void pkv_task6_debug_build(int n) { g_t6_dbgbuild = n; }
  void pkv_task6_set_outlier(double chi2, int nmax) { kT6OutlierChi2 = chi2; kT6MaxOutliers = nmax; }
  void pkv_task6_set_order_eps(double e) { g_t6_order_eps = e; }

  void pkv_task6_add_event(const Event *evp) {
    if (evp == nullptr) { printf("pkv_task6_add_event: null event\n"); return; }
    const Event &ev = *evp;
    const TrackerInfo &ti = tinfo();

    // MkBuilder::fit_tracks() builds its own flags rather than reading
    // PropagationConfig: PF_use_param_b_field | PF_apply_material.
    PropagationFlags pf = ti.prop_config().backward_fit_pflags;
    pf.use_param_b_field = true;
    const MatMode mm = MM_AllBefore;   // what the code does today

    ++g_t6.n_events;
    g_t6.n_sim += (long) ev.simTracks_.size();
    g_t6.n_seeds += (long) ev.seedTracks_.size();

    // ---- selection: labels from the seed tracks -------------------------
    std::vector<char> want(ev.simTracks_.size(), 0);
    for (const Track &sd : ev.seedTracks_) {
      const int l = sd.label();
      if (l < 0) continue;
      ++g_t6.n_seed_lab;
      if (l < (int) ev.simTracks_.size()) want[l] = 1;
    }
    long n_uniq = 0;
    for (char c : want) n_uniq += c;
    g_t6.n_uniq_lab += n_uniq;

    std::mt19937_64 rng(0xB0BAu);

    for (size_t il = 0; il < want.size(); ++il) {
      if (!want[il]) continue;
      const Track &st = ev.simTracks_[il];

      T6Helix hx;
      hx.init(st);
      if (!hx.ok) ++g_t6.n_nohelix;

      std::vector<SeqHit> seq;
      std::vector<T6HitInfo> info;
      t6_build_hits(ev, ti, st, hx, seq, info);
      if ((int) seq.size() < kT6MinHits) continue;
      ++g_t6.n_after_hits;

      // Raw-list diagnostics, BEFORE any shaping: the max-hits-per-layer
      // distribution that the old blunt veto acted on, and the geometric
      // residuals of every hit to the sim track's ideal helix.
      bool legacy_ok = true;
      {
        std::map<int, int> per_layer;
        int mx = 0, mxl = -1;
        for (const SeqHit &q : seq) { const int c = ++per_layer[q.layer]; if (c > mx) { mx = c; mxl = q.layer; } }
        ++g_t6.maxperlayer[std::min(15, mx)];
        legacy_ok = (mx <= 2);
        if (mx >= 3) {
          double amin = 1e30, amax = -1e30, dcmax = 0;
          for (size_t q = 0; q < seq.size(); ++q) if (seq[q].layer == mxl) {
            amin = std::min(amin, info[q].alpha); amax = std::max(amax, info[q].alpha);
            dcmax = std::max(dcmax, std::abs(info[q].dcirc));
          }
          ++g_t6.n_mhl;
          g_t6.st_mhl_aspan.add(amax - amin);
          g_t6.st_mhl_dcirc.add(dcmax);
        }
      }
      if (hx.ok) {
        for (const T6HitInfo &q : info) {
          g_t6.st_dcirc.add(std::abs(q.dcirc));
          g_t6.st_dz.add(std::abs(q.dz));
          g_t6.st_cosap.add(q.cosap);
        }
        g_t6.st_alpha_span.add(info.back().alpha - info.front().alpha);
      }

      const double pt_sim = st.pT();
      const double eta_sim = st.momEta();
      const int ipb = t6_pt_bin(pt_sim);
      ++g_t6.n_tot_by_ptbin[ipb];

      // --- Stage 2: apex truncation ---
      if (hx.ok && (g_t6_shaping & kT6ShApex)) {
        const int ncut = t6_truncate_at_apex(seq, info);
        if (ncut > 0) {
          ++g_t6.n_tr_truncated;
          g_t6.n_hits_apex += ncut;
          ++g_t6.n_apex_by_ptbin[ipb];
        }
      }
      if ((int) seq.size() < kT6MinHits) { ++g_t6.n_lost_apex; continue; }
      ++g_t6.n_after_apex;

      // --- Stage 3: severe-outlier removal ---
      if (g_t6_shaping & kT6ShOutlier) {
        const int nrem = t6_remove_outliers(ti, st, seq, info, pf, mm, rng);
        if (nrem > 0) { ++g_t6.n_tr_outl; g_t6.n_hits_outl += nrem; }
      }
      if ((int) seq.size() < kT6MinHits) { ++g_t6.n_lost_outl; continue; }
      ++g_t6.n_after_outl;

      // <= 2 hits per mkFit layer -- kept, but now as a CONSEQUENCE of the
      // shaping rather than a track-level veto: a track with 3+ in a layer
      // because it came back is truncated at the apex above and survives
      // here, and only what is still multiply-hit afterwards is dropped.
      {
        std::map<int, int> per_layer;
        int mx = 0;
        for (const SeqHit &q : seq) mx = std::max(mx, ++per_layer[q.layer]);
        ++g_t6.maxperlayer_shaped[std::min(15, mx)];
        if (mx > g_t6_maxhits_per_layer) { ++g_t6.n_rej_multihit; continue; }
      }
      ++g_t6.n_after_loop;

      if (!(pt_sim >= kT6MinPt) || !fin((float) pt_sim) || !fin((float) eta_sim)) continue;
      ++g_t6.n_after_pt;

      HitPosArr hp;
      fill_hit_positions(seq, false, rng, hp);

      T6Out r;
      const double ipt = t6_fwd_then_bkw(ti, st, seq, hp, pf, mm, 0.5f, r);
      ++g_t6.n_fitted;
      for (float c : r.hitchi2) g_t6.v_hitchi2.push_back(c);

      const int ip = t6_pt_bin(pt_sim), ie = t6_eta_bin(eta_sim);
      const double ipt_sim = 1.0 / pt_sim;
      const int lbl = (int) il;

      if (legacy_ok) ++g_t6.n_leg_fit;
      const bool bad = r.propfail || r.nonfinite || r.badstate;
      if (bad) {
        if (legacy_ok) ++g_t6.n_leg_bad;
        ++g_t6.n_fail[ip][ie];
        ++g_t6.n_fail_tot;
        int cat = 0;                      // propfail
        if (r.badstate) { cat = 2; ++g_t6.n_badstate; }
        if (r.nonfinite) { cat = 1; ++g_t6.n_nonfinite; }
        g_t6.badcases.push_back({ev.evtID(), lbl, cat});
        continue;
      }

      const double d = (ipt - ipt_sim) / ipt_sim;
      if (std::abs(d) > 0.5) { ++g_t6.n_outlier; g_t6.badcases.push_back({ev.evtID(), lbl, 3}); }
      g_t6.res[ip][ie].add(d);
      if (legacy_ok) { g_t6.res_legacy.add(d); g_t6.c2_legacy.add(r.chi2ndf_bkw); }
      ++g_t6.n_ok[ip][ie];
      g_t6.res_all.add(d);
      g_t6.v_res.push_back((float) d);
      g_t6.v_pt.push_back((float) pt_sim);
      g_t6.v_eta.push_back((float) eta_sim);
      g_t6.v_c2f.push_back((float) r.chi2ndf_fwd);
      g_t6.v_c2b.push_back((float) r.chi2ndf_bkw);
      g_t6.v_evt.push_back(ev.evtID());
      g_t6.v_lbl.push_back(lbl);
      g_t6.nhits.add((double) seq.size());
      g_t6.chi2ndf_fwd.add(r.chi2ndf_fwd);
      g_t6.chi2ndf_bkw.add(r.chi2ndf_bkw);
      if (std::abs(d) > 1.0) ++g_t6.n_gross;

      // Prior sensitivity: same fit with a 10x looser 1/pT prior.  If the
      // numbers below move, the truth-seeded prior is leaking into the answer.
      {
        T6Out r2;
        const double ipt2 = t6_fwd_then_bkw(ti, st, seq, hp, pf, mm, 5.0f, r2);
        if (!(r2.propfail || r2.nonfinite || r2.badstate))
          g_t6.res_loose.add((ipt2 - ipt_sim) / ipt_sim);
      }
    }
  }

  void pkv_task6_write_badcases(const char *path, const char *sample) {
    static const char *kCat[4] = {"propfail", "nonfinite", "badstate", "outlier"};
    auto &v = g_t6.badcases;
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    FILE *f = std::fopen(path, "w");
    if (!f) { printf("pkv_task6: CANNOT open %s for writing\n", path); return; }
    std::fprintf(f, "# %s\n", sample ? sample : "(sample path not given)");
    std::fprintf(f, "# events processed: %ld ; columns: event  sim_label  category\n", g_t6.n_events);
    std::fprintf(f, "# categories: propfail nonfinite badstate outlier\n");
    long cnt[4] = {0, 0, 0, 0};
    for (const auto &e : v) { std::fprintf(f, "%d %d %s\n", e[0], e[1], kCat[e[2]]); ++cnt[e[2]]; }
    std::fclose(f);
    printf("\nBAD-CASE LIST written to %s : %ld lines\n", path, (long) v.size());
    for (int i = 0; i < 4; ++i) printf("    %-10s %6ld\n", kCat[i], cnt[i]);
  }

  void pkv_task6_report() {
    printf("\n"
           "============================================================================\n"
           "TASK 6 -- pT resolution from a forward+backward fit of REAL SIM TRACKS\n"
           "============================================================================\n");
    printf("fit: MkFitter::fwdFitFitTracks then bkReFitFitTracks, reproduced with the\n"
           "     suite's run_bkfit() kernel (propagateHelixToPlaneMPlex +\n"
           "     kalmanOperationPlaneLocal + kalmanCheckChargeFlip).\n");
    printf("material mode: %s ; use_param_b_field=1 (as MkBuilder::fit_tracks)\n",
           matmode_name(MM_AllBefore));
    printf("Config::usePropToPlane=%d usePtMultScat=%d nSStepsInProp2Plane=%d Bfield=%.4f\n",
           (int) Config::usePropToPlane, (int) Config::usePtMultScat,
           Config::nSStepsInProp2Plane, Config::Bfield);
    printf("quantity: d(1/pT)/(1/pT) = (ipt_fit - ipt_sim)/ipt_sim, sim truth at vertex.\n");
    printf("seed prior at the vertex: sigma_pos 1 cm, sigma_ang 50 mrad, sigma(1/pT) 50%%.\n");

    printf("\n  hit-to-its-own-module-plane distance |n.(hit-plPnt)| [cm]: med %.4g p90 %.4g max %.4g (n=%ld)\n",
           g_t6.hit_plane_dist.med(), g_t6.hit_plane_dist.pct(90.0), g_t6.hit_plane_dist.max(),
           (long) g_t6.hit_plane_dist.n());
    printf("  hit-to-module-centre distance [cm]: med %.4g p90 %.4g max %.4g\n",
           g_t6.hit_mod_dist.med(), g_t6.hit_mod_dist.pct(90.0), g_t6.hit_mod_dist.max());
    printf("\n  hits dropped: module id out of range %ld, degenerate module normal %ld\n",
           g_t6_badmod, g_t6_zeronrm);
    printf("\nSELECTION FUNNEL (%ld events)\n", g_t6.n_events);
    printf("  sim tracks in file                         %10ld\n", g_t6.n_sim);
    printf("  seed tracks                                %10ld\n", g_t6.n_seeds);
    printf("  seeds with a label >= 0                    %10ld\n", g_t6.n_seed_lab);
    printf("  -> distinct sim labels from seeds          %10ld\n", g_t6.n_uniq_lab);
    printf("  -> after n_rec_hits >= %d                   %10ld\n", kT6MinHits, g_t6.n_after_hits);
    printf("  -> after APEX truncation                   %10ld  (lost %ld to <%d hits)\n",
           g_t6.n_after_apex, g_t6.n_lost_apex, kT6MinHits);
    printf("  -> after OUTLIER removal                   %10ld  (lost %ld to <%d hits)\n",
           g_t6.n_after_outl, g_t6.n_lost_outl, kT6MinHits);
    printf("  -> after <=%d hits per layer                %10ld  (rejected %ld)\n",
           g_t6_maxhits_per_layer, g_t6.n_after_loop, g_t6.n_rej_multihit);
    printf("  -> after pT >= %.1f GeV                     %10ld\n", kT6MinPt, g_t6.n_after_pt);

    printf("\nSHAPING (mode 0x%x: order=%d apex=%d outlier=%d)\n", g_t6_shaping,
           (g_t6_shaping & kT6ShOrder) ? 1 : 0, (g_t6_shaping & kT6ShApex) ? 1 : 0,
           (g_t6_shaping & kT6ShOutlier) ? 1 : 0);
    printf("  sim tracks with an unusable production state %8ld\n", g_t6.n_nohelix);
    printf("  tracks TRUNCATED at the apex               %10ld  (%ld hits removed)\n",
           g_t6.n_tr_truncated, g_t6.n_hits_apex);
    printf("  tracks with an OUTLIER removed             %10ld  (%ld hits removed)\n",
           g_t6.n_tr_outl, g_t6.n_hits_outl);
    if (g_t6.st_outchi2.n())
      printf("     per-hit chi2 of the removed hits: med %.4g p90 %.4g max %.4g (cut %.0f)\n",
             g_t6.st_outchi2.med(), g_t6.st_outchi2.pct(90.0), g_t6.st_outchi2.max(), kT6OutlierChi2);
  printf("  outlier cut: per-hit chi2 > %.0f, at most %d removals per track\n", kT6OutlierChi2, kT6MaxOutliers);
    printf("  apex truncations by pT bin (fired / tracks):");
    for (int i = 0; i < kT6NPt; ++i)
      printf(" [%.3g-%.3g]%ld/%ld", kT6PtEdges[i], std::min(kT6PtEdges[i+1], 1.0e3),
             g_t6.n_apex_by_ptbin[i], g_t6.n_tot_by_ptbin[i]);
    printf("\n");
    printf("  max hits in any one layer, BEFORE shaping:");
    for (int i = 1; i < 16; ++i) if (g_t6.maxperlayer[i]) printf(" %d:%ld", i, g_t6.maxperlayer[i]);
    printf("\n  max hits in any one layer, AFTER  shaping:");
    for (int i = 1; i < 16; ++i) if (g_t6.maxperlayer_shaped[i]) printf(" %d:%ld", i, g_t6.maxperlayer_shaped[i]);
    printf("\n");
    if (g_t6.n_order_tracks)
      printf("  ordering: alpha (path length) vs 3D radius -- %.3g inversions per track (%ld tracks)\n",
             (double) g_t6.n_order_swaps / g_t6.n_order_tracks, g_t6.n_order_tracks);
    printf("  ordering bucket eps = %.4g rad (0 = pure turn angle)\n", g_t6_order_eps);
    if (g_t6.st_cosap.n())
      printf("  cos(momPhi-posPhi) over ALL raw hits: p1 %.4g p10 %.4g med %.4g   (apex cut < %.5f)\n",
             g_t6.st_cosap.pct(1.0), g_t6.st_cosap.pct(10.0), g_t6.st_cosap.med(), kT6SinApex);
    if (g_t6.st_dcirc.n())
      printf("  |hit - ideal-helix circle| [cm]: med %.4g p90 %.4g p99 %.4g max %.4g\n",
             g_t6.st_dcirc.med(), g_t6.st_dcirc.pct(90.0), g_t6.st_dcirc.pct(99.0), g_t6.st_dcirc.max());
    if (g_t6.st_dz.n())
      printf("  |dz to ideal helix|       [cm]: med %.4g p90 %.4g p99 %.4g max %.4g\n",
             g_t6.st_dz.med(), g_t6.st_dz.pct(90.0), g_t6.st_dz.pct(99.0), g_t6.st_dz.max());
    if (g_t6.res_legacy.n())
      printf("  LEGACY SUBSET (raw max hits/layer <= 2, i.e. what the old cut kept):\n"
             "     fits %ld, bad %ld (%.2f%%), clean %ld, median %.4g, width %.4g, chi2/ndf bkw med %.4g\n",
             g_t6.n_leg_fit, g_t6.n_leg_bad,
             g_t6.n_leg_fit ? 100.0 * g_t6.n_leg_bad / g_t6.n_leg_fit : 0.0,
             (long) g_t6.res_legacy.n(), g_t6.res_legacy.med(),
             0.5 * (g_t6.res_legacy.pct(84.0) - g_t6.res_legacy.pct(16.0)), g_t6.c2_legacy.med());
    if (g_t6.n_mhl)
      printf("  tracks with 3+ raw hits in one layer: %ld ; within that layer, turn-angle span"
             " med %.4g p90 %.4g max %.4g rad ; max |dcirc| med %.4g p90 %.4g cm\n",
             g_t6.n_mhl, g_t6.st_mhl_aspan.med(), g_t6.st_mhl_aspan.pct(90.0), g_t6.st_mhl_aspan.max(),
             g_t6.st_mhl_dcirc.med(), g_t6.st_mhl_dcirc.pct(90.0));
    if (g_t6.st_alpha_span.n())
      printf("  turn-angle span of the raw hit list [rad]: med %.4g p90 %.4g max %.4g\n",
             g_t6.st_alpha_span.med(), g_t6.st_alpha_span.pct(90.0), g_t6.st_alpha_span.max());
    printf("  -> fits run                                %10ld\n", g_t6.n_fitted);
    printf("     of which BAD (any category)             %10ld  (%.2f%%)\n",
           g_t6.n_fail_tot, g_t6.n_fitted ? 100.0 * g_t6.n_fail_tot / g_t6.n_fitted : 0.0);
    printf("        nonfinite end state                  %10ld\n", g_t6.n_nonfinite);
    printf("        badstate  end state                  %10ld\n", g_t6.n_badstate);
    printf("        propfail  (the rest)                 %10ld\n",
           g_t6.n_fail_tot - g_t6.n_nonfinite - g_t6.n_badstate);
    printf("     fwd pass flagged / bkw pass flagged     %10ld %9ld\n", g_t6.n_fail_fwd, g_t6.n_fail_bkw);
    printf("     structural first-backward-step flags    %10ld   (NOT counted as bad)\n",
           g_t6.n_fail_firststep);
    printf("     outliers |d(1/pT)/(1/pT)|>0.5 (kept)    %10ld\n", g_t6.n_outlier);
    printf("     clean fits entering the resolution      %10ld\n", (long) g_t6.res_all.n());
    printf("     |d(1/pT)/(1/pT)| > 1 among those        %10ld\n", g_t6.n_gross);

    if (g_t6.res_all.n()) {
      printf("\nOVERALL (clean fits only)\n");
      printf("  median d(1/pT)/(1/pT)   %10.4g\n", g_t6.res_all.med());
      printf("  robust width (p84-p16)/2 %9.4g\n",
             0.5 * (g_t6.res_all.pct(84.0) - g_t6.res_all.pct(16.0)));
      printf("  p16 / p84               %10.4g %10.4g\n", g_t6.res_all.pct(16.0), g_t6.res_all.pct(84.0));
      printf("  p1  / p99               %10.4g %10.4g\n", g_t6.res_all.pct(1.0), g_t6.res_all.pct(99.0));
      printf("  median n_hits fitted    %10.4g\n", g_t6.nhits.med());
      printf("  chi2/ndf fwd  med/p90   %10.4g %10.4g\n", g_t6.chi2ndf_fwd.med(), g_t6.chi2ndf_fwd.pct(90.0));
      printf("  chi2/ndf bkw  med/p90   %10.4g %10.4g\n", g_t6.chi2ndf_bkw.med(), g_t6.chi2ndf_bkw.pct(90.0));
      if (g_t6.res_loose.n())
        printf("  prior x10 looser: median %9.4g  width %9.4g  (n=%ld) -- compare above\n",
               g_t6.res_loose.med(),
               0.5 * (g_t6.res_loose.pct(84.0) - g_t6.res_loose.pct(16.0)),
               (long) g_t6.res_loose.n());
    }

    {
      printf("\n  flagged propagations by LAYER (fwd+bkw, both prior settings):\n   ");
      int shown = 0;
      for (int l = 0; l < 64; ++l) if (g_t6.faillayer[l]) {
        printf(" L%d:%ld", l, g_t6.faillayer[l]);
        if (++shown % 10 == 0) printf("\n   ");
      }
      printf("\n");
    }
    t6_print_grid("n entries in the resolution (clean fits)", "n", t6_get_n, false);
    t6_print_grid("MEDIAN d(1/pT)/(1/pT)  (bias)", "median", t6_get_med, false);
    t6_print_grid("ROBUST WIDTH (p84-p16)/2 of d(1/pT)/(1/pT)", "width", t6_get_wid, false);
    t6_print_grid("n fits attempted (clean + flagged)", "n", t6_get_tot, false);
    t6_print_grid("FRACTION with a propagation FailFlag", "fail frac", t6_get_fail, true);
    printf("\n");
  }

  // ==========================================================================
  // ROOT OUTPUT AND THE THREE-WAY SOLVER COMPARISON
  //
  // Why plain ROOT histograms and not AnRun/CanvasGroup: CanvasGroup::Add()
  // takes an ROOT::RDF::RResultPtr<TH1> and nothing else, and AnRun is built
  // around an RDataFrame over the per-Event trace collections.  Task 6 is a
  // scalar C++ loop over sim tracks that produces plain doubles and never
  // creates a trace record, so wrapping it in an RDF would mean inventing a
  // dataframe purely to satisfy the container.  What IS kept is the output
  // CONVENTION: <prefix>.root of canvases + histograms next to a <prefix>.txt
  // of the text report, which is exactly what `show-anrun <prefix>.root`
  // consumes.
  //
  // Snapshots: pkv_task6_stash("label") copies the current accumulator aside so
  // several solver configurations can be run over the SAME tracks in ONE
  // process and then overlaid.
  // ==========================================================================

  namespace {
    struct T6Snap { std::string label; T6Acc acc; };
    std::vector<T6Snap> g_t6_snaps;

#ifdef WITH_ROOT
    const int kT6Col[6] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7};

    TH1D *t6_h1(const std::string &pfx, const char *name, const char *title,
                int nb, double lo, double hi) {
      TH1D *h = new TH1D((pfx + "_" + name).c_str(), title, nb, lo, hi);
      h->SetDirectory(nullptr);
      return h;
    }

    // One overlay canvas out of the same-named histogram from every snapshot.
    TCanvas *t6_overlay(const char *cname, const char *ctitle,
                        const std::vector<std::vector<TH1D *>> &per_snap, int which,
                        const std::vector<std::string> &labels, bool logy) {
      TCanvas *c = new TCanvas(cname, ctitle, 900, 650);
      if (logy) c->SetLogy();
      double ymax = 0;
      for (const auto &v : per_snap) ymax = std::max(ymax, v[which]->GetMaximum());
      TLegend *lg = new TLegend(0.60, 0.70, 0.98, 0.92);
      lg->SetFillStyle(0);
      for (size_t i = 0; i < per_snap.size(); ++i) {
        TH1D *h = per_snap[i][which];
        h->SetLineColor(kT6Col[i % 6]);
        h->SetLineWidth(2);
        h->SetStats(0);
        h->SetMaximum(1.35 * std::max(ymax, 1.0));
        if (logy) h->SetMinimum(0.5);
        h->Draw(i == 0 ? "HIST" : "HIST SAME");
        lg->AddEntry(h, Form("%s  (n=%.0f)", labels[i].c_str(), h->GetEntries()), "l");
      }
      lg->Draw();
      return c;
    }
#endif
  }  // anonymous namespace

  void pkv_task6_stash(const char *label) {
    g_t6_snaps.push_back({label ? label : "unnamed", g_t6});
    printf("\npkv_task6_stash: snapshot '%s' kept (%ld clean fits, %ld bad, %ld per-hit chi2)\n",
           g_t6_snaps.back().label.c_str(), (long) g_t6.v_res.size(),
           g_t6.n_fail_tot, (long) g_t6.v_hitchi2.size());
  }

  void pkv_task6_clear_stash() { g_t6_snaps.clear(); }

  void pkv_task6_write_root(const char *path) {
#ifndef WITH_ROOT
    printf("pkv_task6_write_root: built without WITH_ROOT, nothing written.\n");
    (void) path;
#else
    if (g_t6_snaps.empty()) {
      printf("pkv_task6_write_root: no snapshots stashed; stashing the current one as 'current'.\n");
      pkv_task6_stash("current");
    }
    TFile f(path, "RECREATE");
    if (f.IsZombie()) { printf("pkv_task6_write_root: CANNOT open %s\n", path); return; }

    std::vector<std::string> labels;
    std::vector<std::vector<TH1D *>> H;          // [snapshot][which]
    std::vector<std::vector<TH1D *>> Hcell;      // [snapshot][ip*kT6NEta+ie]

    enum { kDipt = 0, kDiptWide, kC2fLog, kC2bLog, kHitChi2Log, kNH };

    for (const T6Snap &sn : g_t6_snaps) {
      labels.push_back(sn.label);
      const std::string &p = sn.label;
      std::vector<TH1D *> v(kNH, nullptr);
      v[kDipt]     = t6_h1(p, "dipt", "d(1/pT)/(1/pT);(ipt_{fit}-ipt_{sim})/ipt_{sim};fits", 200, -0.1, 0.1);
      v[kDiptWide] = t6_h1(p, "dipt_wide", "d(1/pT)/(1/pT), full range;(ipt_{fit}-ipt_{sim})/ipt_{sim};fits", 200, -1.0, 1.0);
      v[kC2fLog]   = t6_h1(p, "log10_chi2ndf_fwd", "forward pass;log_{10}(#chi^{2}/ndf);fits", 120, -3.0, 6.0);
      v[kC2bLog]   = t6_h1(p, "log10_chi2ndf_bkw", "backward pass;log_{10}(#chi^{2}/ndf);fits", 120, -3.0, 6.0);
      v[kHitChi2Log] = t6_h1(p, "log10_chi2_hit", "per-hit #chi^{2};log_{10}(#chi^{2}_{hit});hits", 160, -6.0, 10.0);
      for (float d : sn.acc.v_res) { v[kDipt]->Fill(d); v[kDiptWide]->Fill(d); }
      for (float c : sn.acc.v_c2f) v[kC2fLog]->Fill(std::log10(std::max(1e-30f, c)));
      for (float c : sn.acc.v_c2b) v[kC2bLog]->Fill(std::log10(std::max(1e-30f, c)));
      for (float c : sn.acc.v_hitchi2) if (fin(c) && c > 0) v[kHitChi2Log]->Fill(std::log10(c));
      H.push_back(v);

      // per (pT, eta) cell
      std::vector<TH1D *> vc(kT6NPt * kT6NEta, nullptr);
      for (int ip = 0; ip < kT6NPt; ++ip)
        for (int ie = 0; ie < kT6NEta; ++ie) {
          char nm[64], ti[160];
          snprintf(nm, sizeof(nm), "dipt_pt%d_eta%d", ip, ie);
          snprintf(ti, sizeof(ti), "pT %.3g-%.3g, |eta| %.1f-%.1f;d(1/pT)/(1/pT);fits",
                   kT6PtEdges[ip], std::min(kT6PtEdges[ip + 1], 1.0e3),
                   kT6EtaEdges[ie], std::min(kT6EtaEdges[ie + 1], 1.0e2));
          vc[ip * kT6NEta + ie] = t6_h1(p, nm, ti, 60, -0.15, 0.15);
        }
      for (size_t i = 0; i < sn.acc.v_res.size(); ++i) {
        const int ip = t6_pt_bin(sn.acc.v_pt[i]), ie = t6_eta_bin(sn.acc.v_eta[i]);
        vc[ip * kT6NEta + ie]->Fill(sn.acc.v_res[i]);
      }
      Hcell.push_back(vc);
    }

    // --- per-snapshot directories ---------------------------------------
    for (size_t i = 0; i < g_t6_snaps.size(); ++i) {
      TDirectory *d = f.mkdir(g_t6_snaps[i].label.c_str());
      d->cd();
      for (TH1D *h : H[i]) h->Write();
      for (TH1D *h : Hcell[i]) h->Write();

      // propagation-failure counts and fractions per cell
      const T6Acc &a = g_t6_snaps[i].acc;
      TH2D *hf = new TH2D(("fail_frac_" + g_t6_snaps[i].label).c_str(),
                          "fraction of fits with a propagation FailFlag;pT bin;|eta| bin",
                          kT6NPt, -0.5, kT6NPt - 0.5, kT6NEta, -0.5, kT6NEta - 0.5);
      TH2D *hn = new TH2D(("fail_n_" + g_t6_snaps[i].label).c_str(),
                          "n fits with a propagation FailFlag;pT bin;|eta| bin",
                          kT6NPt, -0.5, kT6NPt - 0.5, kT6NEta, -0.5, kT6NEta - 0.5);
      TH2D *ht = new TH2D(("fit_n_" + g_t6_snaps[i].label).c_str(),
                          "n fits attempted;pT bin;|eta| bin",
                          kT6NPt, -0.5, kT6NPt - 0.5, kT6NEta, -0.5, kT6NEta - 0.5);
      for (int ip = 0; ip < kT6NPt; ++ip) {
        char lab[32];
        if (kT6PtEdges[ip + 1] > 1e8) snprintf(lab, sizeof(lab), "pT>%.0f", kT6PtEdges[ip]);
        else snprintf(lab, sizeof(lab), "%.3g-%.3g", kT6PtEdges[ip], kT6PtEdges[ip + 1]);
        hf->GetXaxis()->SetBinLabel(ip + 1, lab);
        hn->GetXaxis()->SetBinLabel(ip + 1, lab);
        ht->GetXaxis()->SetBinLabel(ip + 1, lab);
      }
      for (int ie = 0; ie < kT6NEta; ++ie) {
        char lab[32];
        if (kT6EtaEdges[ie + 1] > 100) snprintf(lab, sizeof(lab), ">%.1f", kT6EtaEdges[ie]);
        else snprintf(lab, sizeof(lab), "%.1f-%.1f", kT6EtaEdges[ie], kT6EtaEdges[ie + 1]);
        hf->GetYaxis()->SetBinLabel(ie + 1, lab);
        hn->GetYaxis()->SetBinLabel(ie + 1, lab);
        ht->GetYaxis()->SetBinLabel(ie + 1, lab);
      }
      for (int ip = 0; ip < kT6NPt; ++ip)
        for (int ie = 0; ie < kT6NEta; ++ie) {
          const long tot = a.n_ok[ip][ie] + a.n_fail[ip][ie];
          hn->SetBinContent(ip + 1, ie + 1, (double) a.n_fail[ip][ie]);
          ht->SetBinContent(ip + 1, ie + 1, (double) tot);
          if (tot) hf->SetBinContent(ip + 1, ie + 1, (double) a.n_fail[ip][ie] / tot);
        }
      hf->SetDirectory(nullptr); hn->SetDirectory(nullptr); ht->SetDirectory(nullptr);
      hf->Write(); hn->Write(); ht->Write();

      // flagged-propagation count by layer
      TH1D *hl = t6_h1(g_t6_snaps[i].label, "faillayer",
                       "flagged propagations by mkFit layer;layer;count", 64, -0.5, 63.5);
      for (int l = 0; l < 64; ++l) hl->SetBinContent(l + 1, (double) a.faillayer[l]);
      hl->Write();
      f.cd();
    }

    // --- overlay canvases -------------------------------------------------
    TCanvas *c1 = t6_overlay("c_dipt", "d(1/pT)/(1/pT)", H, kDipt, labels, false);
    TCanvas *c2 = t6_overlay("c_dipt_wide", "d(1/pT)/(1/pT), wide", H, kDiptWide, labels, true);
    TCanvas *c3 = t6_overlay("c_chi2ndf_bkw", "backward-pass chi2/ndf", H, kC2bLog, labels, true);
    TCanvas *c4 = t6_overlay("c_chi2ndf_fwd", "forward-pass chi2/ndf", H, kC2fLog, labels, true);
    TCanvas *c5 = t6_overlay("c_chi2_hit", "per-hit chi2", H, kHitChi2Log, labels, true);
    c1->Write(); c2->Write(); c3->Write(); c4->Write(); c5->Write();

    TCanvas *cc = new TCanvas("c_dipt_cells", "d(1/pT)/(1/pT) per (pT, eta) cell", 1600, 1100);
    cc->Divide(kT6NEta, kT6NPt);
    for (int ip = 0; ip < kT6NPt; ++ip)
      for (int ie = 0; ie < kT6NEta; ++ie) {
        cc->cd(ip * kT6NEta + ie + 1);
        double ymax = 0;
        for (size_t i = 0; i < Hcell.size(); ++i) ymax = std::max(ymax, Hcell[i][ip * kT6NEta + ie]->GetMaximum());
        for (size_t i = 0; i < Hcell.size(); ++i) {
          TH1D *h = Hcell[i][ip * kT6NEta + ie];
          h->SetLineColor(kT6Col[i % 6]);
          h->SetStats(0);
          h->SetMaximum(1.25 * std::max(ymax, 1.0));
          h->Draw(i == 0 ? "HIST" : "HIST SAME");
        }
      }
    cc->Write();

    TCanvas *cf = new TCanvas("c_fail_frac", "propagation-failure fraction per cell",
                              450 * (int) g_t6_snaps.size(), 450);
    cf->Divide((int) g_t6_snaps.size(), 1);
    for (size_t i = 0; i < g_t6_snaps.size(); ++i) {
      cf->cd((int) i + 1);
      TH2D *hf = (TH2D *) f.Get((g_t6_snaps[i].label + "/fail_frac_" + g_t6_snaps[i].label).c_str());
      if (hf) { hf->SetStats(0); hf->SetTitle(g_t6_snaps[i].label.c_str()); hf->Draw("COLZ TEXT"); }
    }
    cf->Write();

    f.Close();

    // --- PAIRED comparison: same (event, sim label), snapshot i vs snapshot 0.
    // This is the question "is the fix visible at all" asked per track rather
    // than on an aggregate, so it is not limited by the ~2%/sqrt(2n) noise on a
    // robust width.  Only tracks that produced a clean fit in BOTH snapshots
    // can be paired.
    if (g_t6_snaps.size() > 1) {
      printf("\nPAIRED per-track comparison of d(1/pT)/(1/pT), vs snapshot '%s'\n",
             g_t6_snaps[0].label.c_str());
      const T6Acc &a0 = g_t6_snaps[0].acc;
      std::map<std::pair<int, int>, std::pair<float, float>> m0;   // res, chi2ndf_bkw
      for (size_t i = 0; i < a0.v_res.size(); ++i)
        m0[{a0.v_evt[i], a0.v_lbl[i]}] = {a0.v_res[i], a0.v_c2b[i]};
      for (size_t k = 1; k < g_t6_snaps.size(); ++k) {
        const T6Acc &a = g_t6_snaps[k].acc;
        Stats dif, rel;
        long npair = 0, nsame = 0, c2_better = 0, c2_worse = 0;
        Stats c2ratio;
        for (size_t i = 0; i < a.v_res.size(); ++i) {
          auto it = m0.find({a.v_evt[i], a.v_lbl[i]});
          if (it == m0.end()) continue;
          ++npair;
          const double d = a.v_res[i] - it->second.first;
          if (d == 0.0) ++nsame;
          dif.add(std::abs(d));
          rel.add(std::abs(a.v_res[i]) - std::abs(it->second.first));  // did it get closer to 0?
          const double c_new = a.v_c2b[i], c_old = it->second.second;
          if (c_new < c_old) ++c2_better; else if (c_new > c_old) ++c2_worse;
          if (c_old > 0 && c_new > 0) c2ratio.add(std::log10(c_new / c_old));
        }
        printf("  %-18s paired %5ld ; bit-identical %5ld (%.1f%%) ; |delta| med %.3g p90 %.3g max %.3g ;"
               " median(|d_new|-|d_old|) %+.3g\n",
               g_t6_snaps[k].label.c_str(), npair, nsame,
               npair ? 100.0 * nsame / npair : 0.0,
               dif.med(), dif.pct(90.0), dif.max(), rel.med());
        printf("  %-18s   chi2/ndf bkw PAIRED: better %ld / worse %ld (sign test) ;"
               " log10(new/old) med %+.3g p10 %+.3g p90 %+.3g\n",
               "", c2_better, c2_worse, c2ratio.med(), c2ratio.pct(10.0), c2ratio.pct(90.0));
        // Aggregates restricted to the PAIRED subset, so a change cannot be an
        // artefact of the two snapshots having kept different tracks.
        Stats ro, rn, co, cn;
        for (size_t i = 0; i < a.v_res.size(); ++i) {
          auto it = m0.find({a.v_evt[i], a.v_lbl[i]});
          if (it == m0.end()) continue;
          ro.add(it->second.first);  rn.add(a.v_res[i]);
          co.add(it->second.second); cn.add(a.v_c2b[i]);
        }
        printf("  %-18s   on the PAIRED subset ONLY: width old %.4g -> new %.4g ;"
               " chi2/ndf med %.4g -> %.4g ; p90 %.4g -> %.4g ; p99 %.4g -> %.4g\n",
               "", 0.5 * (ro.pct(84.0) - ro.pct(16.0)), 0.5 * (rn.pct(84.0) - rn.pct(16.0)),
               co.med(), cn.med(), co.pct(90.0), cn.pct(90.0), co.pct(99.0), cn.pct(99.0));
      }
    }

    printf("\npkv_task6_write_root: wrote %s with %d snapshot(s):",
           path, (int) g_t6_snaps.size());
    for (const auto &l : labels) printf(" %s", l.c_str());
    printf("\n  view with:  show-anrun %s   (put the text report in the matching .txt)\n", path);
#endif
  }

  // Runtime getS root-form selection, so the A/B needs no rebuild.
  void pkv_set_solver(int stable_root) {
    g_getS_stable_root = (stable_root != 0);
    printf("\npkv_set_solver: g_getS_stable_root=%d -> getS quadratic, %s roots\n",
           (int) g_getS_stable_root, g_getS_stable_root ? "STABLE" : "ORIGINAL");
  }


  void pkv_task6_one_event(const Event *ev) {
    pkv_task6_reset();
    pkv_task6_add_event(ev);
    pkv_task6_report();
  }

  void pkv_run_all() {
    pkv_task0_pea_material();
    pkv_task4_sequence();
    pkv_task1_closure();
    pkv_task2_numerical_jacobian();
    pkv_task3_kalman_limits();
  }

}  // namespace mkfit
