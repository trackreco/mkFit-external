#include "RecoTracker/MkFitCore/standalone/RdfTrace/ValProp.h"
#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
#include "RecoTracker/MkFitCore/standalone/ConfigStandalone.h"

#include "RecoTracker/MkFitCore/interface/Config.h"
#include "RecoTracker/MkFitCore/interface/cms_common_macros.h"
#include "RecoTracker/MkFitCore/interface/PropagationConfig.h"

#include "RecoTracker/MkFitCore/src/Matrix.h"
#include "RecoTracker/MkFitCore/src/PropagationMPlex.h"
#include "RecoTracker/MkFitCore/src/KalmanUtilsMPlex.h"
#include "RecoTracker/MkFitCore/src/MkFinder.h"
#include "RecoTracker/MkFitCore/src/MkFinderV2p2.h"
#include "RecoTracker/MkFitCore/src/V2p2Score.h"
#include "RecoTracker/MkFitCore/src/MkBins.h"

#include "RecoTracker/MkFitCore/interface/Hit.h"
#include "RecoTracker/MkFitCore/interface/Track.h"
#include "RecoTracker/MkFitCore/standalone/Event.h"
#include "RecoTracker/MkFitCore/standalone/TrackExtra.h"
#include "RecoTracker/MkFitCore/standalone/ConfigStandalone.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/MkFitCore/interface/IterationConfig.h"

#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <array>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixDSym.h"

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

namespace mkfit {
  namespace V2p2 = Config::V2p2;

  // Lite trace: keep only pixel-barrel hits and skip the exact double-precision
  // helix-plane solve. That solve plus the sheer volume (~142k records/event) is
  // what makes val_search_event ~45 s/event; with this it is seconds, which is
  // what makes a covariance SCAN affordable.
  static int g_val_search_lite = 0;   // bit0 skip exact solve, bit1 pixel barrel only

  namespace {

    // ----------------------------------------------------------------------
    // Exact uniform-B helix from a production vertex, parametrized by signed
    // 3D path length s. This is the TRUTH: it defines where the planes are,
    // and nothing in it comes from the code under test.

    struct ValHelix {
      double x0 = 0, y0 = 0, z0 = 0;
      double px0 = 0, py0 = 0, pz0 = 0;
      double pt = 1, pmag = 1, theta = 1.57;
      double k = 0;  // cm/GeV, signed
      int chg = 1;

      void init(double a_pt, double eta, int charge, double phi0, double d0, double zv) {
        pt = a_pt;
        chg = charge;
        theta = 2.0 * std::atan(std::exp(-eta));
        px0 = pt * std::cos(phi0);
        py0 = pt * std::sin(phi0);
        pz0 = pt / std::tan(theta);
        pmag = pt / std::sin(theta);
        const double inv_k = ((charge < 0) ? 0.01 : -0.01) * (double)Const::sol * Config::Bfield;
        k = 1.0 / inv_k;
        // Displace the vertex PERPENDICULAR to the momentum in the transverse
        // plane, so d0 is a signed impact parameter and not just a shift along
        // the track (which would be no displacement at all).
        x0 = -d0 * std::sin(phi0);
        y0 = d0 * std::cos(phi0);
        z0 = zv;
      }

      double alpha_of_s(double s) const { return s / (k * pmag); }
      double s_of_alpha(double a) const { return a * k * pmag; }

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
        double p[3], m[3];
        at_s(s, p, m);
        return std::hypot(p[0], p[1]);
      }

      // CCS parameter vector (x, y, z, 1/pT, phi, theta) at path length s.
      void ccs_at_s(double s, float par[6]) const {
        double p[3], m[3];
        at_s(s, p, m);
        par[0] = (float)p[0];
        par[1] = (float)p[1];
        par[2] = (float)p[2];
        par[3] = (float)(1.0 / pt);
        par[4] = (float)std::atan2(m[1], m[0]);
        par[5] = (float)theta;
      }

      // First crossing of radius R at s > 0, or false if the track never gets
      // there. Scan-then-bisect, capped at half a turn: beyond |alpha| = pi a
      // track from the origin is coming back in, and a closure test has no
      // business out there.
      bool first_crossing_r(double R, double &s_out) const {
        const double s_max = std::fabs(M_PI * k * pmag);
        const int N = 4000;
        double s_prev = 0.0, f_prev = r_at_s(0.0) - R;
        for (int i = 1; i <= N; ++i) {
          const double s = s_max * i / N;
          const double f = r_at_s(s) - R;
          if ((f_prev < 0.0) != (f < 0.0)) {
            double lo = s_prev, hi = s, flo = f_prev;
            for (int it = 0; it < 80; ++it) {
              const double sm = 0.5 * (lo + hi);
              const double fm = r_at_s(sm) - R;
              if ((flo < 0.0) == (fm < 0.0)) { lo = sm; flo = fm; } else { hi = sm; }
            }
            s_out = 0.5 * (lo + hi);
            return true;
          }
          s_prev = s;
          f_prev = f;
        }
        return false;
      }
    };

    // ----------------------------------------------------------------------
    // Plane construction.
    //
    // A plane is a point and a normal. The propagator solves n.(x(s) - p) = 0,
    // so the derivative of that equation along the track is n.dx/ds = |p^ . n|
    // -- the INCIDENCE. At incidence 1 the plane is normal to the track and the
    // solve is as well conditioned as it can be; as incidence -> 0 the track
    // runs parallel to the plane, touches rather than crosses, and the root
    // becomes double. That last case cannot occur in a detector (a module a
    // track grazes produces no hit), so it is excluded rather than measured.
    //
    // The tilt is applied in the plane spanned by the momentum and the radial
    // direction, which is where real module tilt lives.
    void make_plane(const double p[3], const double m[3], double cos_inc,
                    float pnt[3], float nrm[3]) {
      const double pm = std::sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
      const double ph[3] = {m[0]/pm, m[1]/pm, m[2]/pm};
      double n[3] = {ph[0], ph[1], ph[2]};
      if (cos_inc < 0.999999) {
        // Rotation axis: perpendicular to both the momentum and the radial
        // direction, so the tilt happens in the (p^, r^) plane.
        const double rr = std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        double rh[3] = {p[0]/rr, p[1]/rr, p[2]/rr};
        double u[3] = {ph[1]*rh[2] - ph[2]*rh[1],
                       ph[2]*rh[0] - ph[0]*rh[2],
                       ph[0]*rh[1] - ph[1]*rh[0]};
        double un = std::sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
        if (un < 1e-6) {  // momentum parallel to the radius: any perpendicular will do
          const double a[3] = {0.0, 0.0, 1.0};
          u[0] = ph[1]*a[2] - ph[2]*a[1];
          u[1] = ph[2]*a[0] - ph[0]*a[2];
          u[2] = ph[0]*a[1] - ph[1]*a[0];
          un = std::sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
        }
        for (int i = 0; i < 3; ++i) u[i] /= un;
        const double th = std::acos(cos_inc), c = std::cos(th), sn = std::sin(th);
        const double dot = u[0]*ph[0] + u[1]*ph[1] + u[2]*ph[2];
        const double cx[3] = {u[1]*ph[2] - u[2]*ph[1],
                              u[2]*ph[0] - u[0]*ph[2],
                              u[0]*ph[1] - u[1]*ph[0]};
        for (int i = 0; i < 3; ++i) n[i] = ph[i]*c + cx[i]*sn + u[i]*dot*(1.0 - c);
      }
      const double nn = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
      for (int i = 0; i < 3; ++i) { pnt[i] = (float)p[i]; nrm[i] = (float)(n[i]/nn); }
    }

    // ----------------------------------------------------------------------
    // Matriplex helpers. Every lane is filled identically and lane 0 is read
    // back: wasteful, and it keeps the call identical to production.

    inline int ls_idx(int i, int j) {
      if (i < j) { const int t = i; i = j; j = t; }
      return i * (i + 1) / 2 + j;
    }

    void fill_lanes(MPlexLV &par, const float p[6]) {
      for (int n = 0; n < NN; ++n)
        for (int i = 0; i < 6; ++i) par(n, i, 0) = p[i];
    }
    void fill_lanes_sym(MPlexLS &err, const float c[21]) {
      for (int n = 0; n < NN; ++n)
        for (int i = 0; i < 21; ++i) err.fArray[i * NN + n] = c[i];
    }
    void fill_lanes_hv(MPlexHV &v, const float p[3]) {
      for (int n = 0; n < NN; ++n)
        for (int i = 0; i < 3; ++i) v(n, i, 0) = p[i];
    }
    void fill_lanes_q(MPlexQI &v, int x) {
      for (int n = 0; n < NN; ++n) v(n, 0, 0) = x;
    }
    void read_lane0(const MPlexLV &par, const MPlexLS &err, float p[6], float c[21]) {
      for (int i = 0; i < 6; ++i) p[i] = par.constAt(0, i, 0);
      for (int i = 0; i < 21; ++i) c[i] = err.fArray[i * NN];
    }

    // A plausible input covariance. Its exact value does not matter -- closure
    // is an internal-consistency test -- but it must be non-degenerate and
    // have off-diagonal terms, so that a jacobian bug cannot hide behind a
    // diagonal input.
    void nominal_cov(const float par[6], float c[21]) {
      std::memset(c, 0, 21 * sizeof(float));
      const float s_pos = 0.02f;             // 200 um
      const float s_ipt = 0.03f * par[3];    // 3% of 1/pT
      const float s_ang = 0.002f;            // 2 mrad
      c[ls_idx(0, 0)] = s_pos * s_pos;
      c[ls_idx(1, 1)] = s_pos * s_pos;
      c[ls_idx(2, 2)] = s_pos * s_pos;
      c[ls_idx(3, 3)] = s_ipt * s_ipt;
      c[ls_idx(4, 4)] = s_ang * s_ang;
      c[ls_idx(5, 5)] = s_ang * s_ang;
      c[ls_idx(0, 1)] = 0.2f * s_pos * s_pos;
      c[ls_idx(3, 4)] = 0.1f * s_ipt * s_ang;
      c[ls_idx(4, 5)] = 0.15f * s_ang * s_ang;
    }

    // ----------------------------------------------------------------------
    // The configurations. The two controls, scanned directly:
    //   dalpha  -- how far the outbound step turns  [rad]
    //   cos_inc -- how squarely it meets plane B
    // Nothing here is a radius. The radius reached is an output.
    const double kDalpha[]  = {0.02, 0.05, 0.10, 0.20, 0.40, 0.80};
    const double kCosInc[]  = {1.00, 0.70, 0.50, 0.30, 0.15};
    constexpr int kNDal = sizeof(kDalpha) / sizeof(kDalpha[0]);
    constexpr int kNInc = sizeof(kCosInc) / sizeof(kCosInc[0]);
    constexpr int kNCfg = kNDal * kNInc;

    // Tracker envelope. A configuration landing outside it is not measured:
    // it is not a propagation the detector ever asks for, and letting it in is
    // how plane "r = 30 cm" ended up at z = 8 m.
    constexpr double kMaxR = 120.0, kMaxZ = 300.0;

  }  // anonymous namespace

  //============================================================================

  void val_gen_closure(const char *out_file, int n_pt, int n_eta, int n_phi,
                       unsigned seed, float d0, float z0) {
    // Material OFF, uniform field: the propagation is then a diffeomorphism and
    // the exact closure residual is ZERO. Anything measured is float
    // arithmetic, or a defect.
    const PropagationFlags pf;  // all flags false

    const double kL0 = std::log10(0.3), kL1 = std::log10(200.0);
    const double dL = (kL1 - kL0) / n_pt;
    const double kE0 = -4.0, kE1 = 4.0;
    const double dE = (kE1 - kE0) / n_eta;
    const double dP = 2.0 * M_PI / n_phi;

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);

    TFile f(out_file, "RECREATE");
    TTree *t_cfg = new TTree("cfgs", "step configurations");
    ValCfgInfo ci;
    t_cfg->Branch("g", &ci);
    for (int ia = 0; ia < kNDal; ++ia)
      for (int ic = 0; ic < kNInc; ++ic) {
        ci.cfg = ia * kNInc + ic;
        ci.dalpha = (float)kDalpha[ia];
        ci.cos_inc = (float)kCosInc[ic];
        t_cfg->Fill();
      }

    TTree *t = new TTree("closure", "one-step-there-and-back propagation closure");
    ValClosure c;
    t->Branch("c", &c);

    int id = 0, track_id = 0;
    long n_try = 0, n_outside = 0;

    for (int ip = 0; ip < n_pt; ++ip)
      for (int ie = 0; ie < n_eta; ++ie)
        for (int iph = 0; iph < n_phi; ++iph)
          for (int iq = 0; iq < 2; ++iq, ++track_id) {
            // Jittered stratification: one draw inside each cell, so every cell
            // is populated but no two tracks share a pT, an eta or a phi.
            const double a_pt = std::pow(10.0, kL0 + (ip + U(rng)) * dL);
            const double eta = kE0 + (ie + U(rng)) * dE;
            const double phi0 = (iph + U(rng)) * dP;
            const int chg = iq ? +1 : -1;

            ValHelix h;
            h.init(a_pt, eta, chg, phi0, d0, z0);

            // PLANE A: at the production point, normal along the momentum there.
            // Incidence 1, so the return leg is as well conditioned as a plane
            // can be and adds no error of its own -- everything measured is a
            // property of the OUTBOUND step.
            double pa[3], ma[3];
            h.at_s(0.0, pa, ma);
            float par_o[6], cov_o[21];
            h.ccs_at_s(0.0, par_o);
            nominal_cov(par_o, cov_o);

            float pntA[3], nrmA[3];
            make_plane(pa, ma, 1.0, pntA, nrmA);

            for (int ia = 0; ia < kNDal; ++ia) {
              // B is placed at a chosen TURN ANGLE, and the radius it reaches
              // is an output. alpha runs with the sign of k, so fold it.
              const double alpha = (h.k < 0 ? -1.0 : 1.0) * kDalpha[ia];
              const double s_b = h.s_of_alpha(alpha);
              double pb[3], mb[3];
              h.at_s(s_b, pb, mb);
              const double r_b = std::hypot(pb[0], pb[1]);

              for (int ic = 0; ic < kNInc; ++ic) {
                ++n_try;
                // Outside the tracker envelope is not a propagation the
                // detector ever asks for.
                if (r_b > kMaxR || std::fabs(pb[2]) > kMaxZ) { ++n_outside; continue; }

                c = ValClosure{};
                c.id = id++;
                c.track_id = track_id;
                c.cfg = ia * kNInc + ic;
                c.pt = (float)a_pt; c.eta = (float)eta; c.phi0 = (float)phi0;
                c.d0 = d0; c.z0 = z0; c.charge = chg;
                c.dalpha_target = (float)kDalpha[ia];
                c.cosinc_target = (float)kCosInc[ic];
                c.r_a = (float)std::hypot(pa[0], pa[1]);  c.z_a = (float)pa[2];
                c.r_b = (float)r_b;                        c.z_b = (float)pb[2];
                c.s_ab = (float)s_b;
                c.dalpha_ab = (float)kDalpha[ia];

                for (int q = 0; q < 3; ++q) { c.pnt_a[q] = pntA[q]; c.nrm_a[q] = nrmA[q]; }
                make_plane(pb, mb, kCosInc[ic], c.pnt_b, c.nrm_b);
                // Keep the exact crossing in double -- see ValStructs.h.
                c.ref_pnt_b[0] = pb[0]; c.ref_pnt_b[1] = pb[1]; c.ref_pnt_b[2] = pb[2];
                c.ref_s_b = s_b;
                {
                  const double na = std::sqrt(ma[0]*ma[0] + ma[1]*ma[1] + ma[2]*ma[2]);
                  const double nb = std::sqrt(mb[0]*mb[0] + mb[1]*mb[1] + mb[2]*mb[2]);
                  c.cos_inc_a = (float)std::fabs((ma[0]*c.nrm_a[0] + ma[1]*c.nrm_a[1] + ma[2]*c.nrm_a[2]) / na);
                  c.cos_inc_b = (float)std::fabs((mb[0]*c.nrm_b[0] + mb[1]*c.nrm_b[1] + mb[2]*c.nrm_b[2]) / nb);
                }

                MPlexQI q_;  fill_lanes_q(q_, chg);
                MPlexHV mpA, mnA, mpB, mnB;
                fill_lanes_hv(mpA, c.pnt_a); fill_lanes_hv(mnA, c.nrm_a);
                fill_lanes_hv(mpB, c.pnt_b); fill_lanes_hv(mnB, c.nrm_b);

                MPlexLV p0;  fill_lanes(p0, par_o);
                MPlexLS e0;  fill_lanes_sym(e0, cov_o);

                // Leg O -> A. The state is already on plane A, so this is a
                // zero-length step; its job is to give the covariance the
                // rank-5 degeneracy the real code always carries. Needed for
                // the covariance closure, harmless for the position one.
                MPlexLS eA, eB, eC;  MPlexLV pA, pB_, pC;  MPlexQI fA{0}, fB{0}, fC{0};
                propagateHelixToPlaneMPlex(e0, p0, q_, mpA, mnA, eA, pA, fA, NN, pf, nullptr);
                // Leg A -> B: the step under test.
                propagateHelixToPlaneMPlex(eA, pA, q_, mpB, mnB, eB, pB_, fB, NN, pf, nullptr);
                // Leg B -> A': the return, onto the same well-conditioned plane A.
                propagateHelixToPlaneMPlex(eB, pB_, q_, mpA, mnA, eC, pC, fC, NN, pf, nullptr);

                c.fail_oa = fA.constAt(0, 0, 0);
                c.fail_ab = fB.constAt(0, 0, 0);
                c.fail_ba = fC.constAt(0, 0, 0);
                read_lane0(pA, eA, c.par_a, c.err_a);
                read_lane0(pB_, eB, c.par_b, c.err_b);
                read_lane0(pC, eC, c.par_c, c.err_c);

                t->Fill();
              }
            }
          }

    t_cfg->Write();
    t->Write();
    f.Close();

    printf("val_gen_closure: %d tracks x %d cfgs = %ld trials, %ld outside the tracker\n"
           "                 envelope (r < %g, |z| < %g cm), %d written to '%s'\n",
           track_id, kNCfg, n_try, n_outside, kMaxR, kMaxZ, id, out_file);
    printf("  one step there and back: plane A at the production point with its normal ALONG\n"
           "  the momentum (incidence 1, so the return adds no error of its own); plane B at a\n"
           "  chosen TURN ANGLE with a chosen incidence. Radius is an output, not a control.\n");
    printf("  dalpha scan [rad] :");
    for (int i = 0; i < kNDal; ++i) printf(" %g", kDalpha[i]);
    printf("\n  incidence scan    :");
    for (int i = 0; i < kNInc; ++i) printf(" %g", kCosInc[i]);
    printf("\n  sample: log10 pT flat over [0.3, 200] GeV in %d cells, eta flat over [-4, +4] in\n"
           "          %d cells, phi in %d cells, both charges, one jittered draw per cell; seed %u\n",
           n_pt, n_eta, n_phi, seed);
    printf("  vertex: d0 = %g cm, z0 = %g cm;  material OFF, uniform B = %g T\n",
           d0, z0, Config::Bfield);
  }

  //============================================================================
  // What turn angles does a REAL track present between consecutive hits?
  //
  // A synthetic scan chooses the turn angle; a detector does not. Without this
  // distribution, a response measured against turn angle is being weighted
  // uniformly over a variable the experiment is not uniform in -- which is how
  // a scan can spend its statistics on a region that never occurs.
  //============================================================================

  namespace {
    std::vector<ValStep> g_val_steps;
    int g_val_evt = 0;
  }

  void val_dalpha_reset() { g_val_steps.clear(); g_val_evt = 0; }

  void val_dalpha_add_event(const Event *ev, int which, float pt_min) {
    if (ev == nullptr) { printf("val_dalpha_add_event: null event\n"); return; }
    const TrackVec &tv = (which == 1) ? ev->simTracks_
                       : (which == 2) ? ev->candidateTracks_
                                      : ev->cmsswTracks_;
    // R_c [cm] = pT / (0.01 * sol * B)
    const double inv_k = 0.01 * (double)Const::sol * Config::Bfield;
    const int evt = g_val_evt++;

    for (int it = 0; it < (int)tv.size(); ++it) {
      const Track &t = tv[it];
      if (t.pT() < pt_min) continue;
      const double Rc = t.pT() / inv_k;
      int prev_lay = -1, prev_idx = -1, step = 0;
      for (int ih = 0; ih < t.nTotalHits(); ++ih) {
        const int idx = t.getHitIdx(ih), lay = t.getHitLyr(ih);
        if (idx < 0 || lay < 0 || lay >= (int)ev->layerHits_.size()) continue;
        if (idx >= (int)ev->layerHits_[lay].size()) continue;
        if (prev_idx >= 0) {
          const Hit &a = ev->layerHits_[prev_lay][prev_idx];
          const Hit &b = ev->layerHits_[lay][idx];
          const double dx = b.x() - a.x(), dy = b.y() - a.y(), dz = b.z() - a.z();
          const double dperp = std::sqrt(dx * dx + dy * dy);
          ValStep vs;
          vs.event = evt; vs.track = it; vs.step = step++;
          vs.lay_a = prev_lay; vs.lay_b = lay;
          vs.pt = t.pT(); vs.eta = t.momEta();
          vs.r_a = (float)std::hypot(a.x(), a.y()); vs.z_a = a.z();
          vs.r_b = (float)std::hypot(b.x(), b.y()); vs.z_b = b.z();
          vs.d_perp = (float)dperp;
          vs.d_3d = (float)std::sqrt(dx * dx + dy * dy + dz * dz);
          // chord -> turn angle. Beyond 2*R_c the two hits cannot be on one
          // circle of this radius; clamp and let the analysis see it.
          vs.dalpha = (float)(2.0 * std::asin(std::min(1.0, dperp / (2.0 * Rc))));
          g_val_steps.push_back(vs);
        }
        prev_lay = lay; prev_idx = idx;
      }
    }
  }

  void val_dalpha_write(const char *out_file) {
    TFile f(out_file, "RECREATE");
    TTree *t = new TTree("steps", "turn angle between consecutive hits of real tracks");
    ValStep vs;
    t->Branch("s", &vs);
    for (const ValStep &v : g_val_steps) { vs = v; t->Fill(); }
    t->Write();
    f.Close();
    printf("val_dalpha_write: %zu steps from %d events -> '%s'\n",
           g_val_steps.size(), g_val_evt, out_file);
  }

  //============================================================================
  // Geometry compatibility of a sample -- see the header.
  //============================================================================

  namespace {
    std::vector<double> g_gc_d;
    long g_gc_bad_sid = 0, g_gc_hits = 0, g_gc_evt = 0;
  }

  void val_geom_check_reset() {
    g_gc_d.clear(); g_gc_bad_sid = 0; g_gc_hits = 0; g_gc_evt = 0;
  }

  void val_geom_check_event(const Event *ev) {
    if (ev == nullptr) { printf("val_geom_check_event: null event\n"); return; }
    const TrackerInfo &ti = Config::TrkInfo;
    ++g_gc_evt;
    for (int lay = 0; lay < (int)ev->layerHits_.size(); ++lay) {
      if (lay >= ti.n_layers()) break;
      const LayerInfo &li = ti.layer(lay);
      for (const Hit &h : ev->layerHits_[lay]) {
        ++g_gc_hits;
        const unsigned int sid = h.detIDinLayer();
        if ((int)sid >= li.n_modules()) { ++g_gc_bad_sid; continue; }
        const ModuleInfo &mi = li.module_info(sid);
        const double d = std::fabs((h.x() - mi.pos[0]) * mi.zdir[0] +
                                   (h.y() - mi.pos[1]) * mi.zdir[1] +
                                   (h.z() - mi.pos[2]) * mi.zdir[2]);
        g_gc_d.push_back(d);
      }
    }
  }

  void val_geom_check_report(const char *sample_name) {
    printf("\n================ geometry compatibility ================\n");
    printf("sample : %s\n", sample_name);
    printf("events : %ld,  hits : %ld\n", g_gc_evt, g_gc_hits);
    printf("module short_id out of range : %ld  (%.4f %%)\n",
           g_gc_bad_sid, g_gc_hits ? 100.0 * g_gc_bad_sid / g_gc_hits : 0.0);
    if (g_gc_d.empty()) { printf("no usable hits\n"); return; }
    std::sort(g_gc_d.begin(), g_gc_d.end());
    auto q = [&](double p) { return g_gc_d[(size_t)(0.01 * p * (g_gc_d.size() - 1))]; };
    printf("|n.(hit - module plane point)|, cm -- a hit is ON its module by construction\n");
    printf("  median %10.4g   p90 %10.4g   p99 %10.4g   max %10.4g\n",
           q(50), q(90), q(99), g_gc_d.back());
    long over = 0; for (double x : g_gc_d) if (x > 0.01) ++over;
    printf("  fraction over 100 um : %.4f %%\n", 100.0 * over / g_gc_d.size());
    printf("\nreference: 2.9e-9 cm median = COMPATIBLE;  0.26 cm median = NOT\n");
    printf("=======================================================\n");
  }

  void val_sample_info(const Event *ev) {
    if (ev == nullptr) { printf("val_sample_info: null event\n"); return; }
    auto describe = [](const char *nm, const TrackVec &tv) {
      if (tv.empty()) { printf("  %-18s empty\n", nm); return; }
      std::map<int, int> algo;
      double ptmin = 1e9, ptmax = -1e9;
      long nh = 0, nlab = 0;
      for (const Track &t : tv) {
        algo[t.algoint()]++;
        ptmin = std::min(ptmin, (double)t.pT());
        ptmax = std::max(ptmax, (double)t.pT());
        nh += t.nFoundHits();
        if (t.label() >= 0) ++nlab;
      }
      printf("  %-18s n=%-7zu  pT %.3g .. %.4g  <nhits> %.1f  labelled %.0f%%  algos:",
             nm, tv.size(), ptmin, ptmax, (double)nh / tv.size(), 100.0 * nlab / tv.size());
      for (auto &a : algo) printf(" %d:%d", a.first, a.second);
      printf("\n");
    };
    printf("\n---- sample contents, one event ----\n");
    printf("  %-18s %zu layers\n", "layerHits_", ev->layerHits_.size());
    long nh = 0; for (auto &v : ev->layerHits_) nh += v.size();
    printf("  %-18s %ld hits\n", "total hits", nh);
    printf("  %-18s %zu\n", "simHitsInfo_", ev->simHitsInfo_.size());
    printf("  %-18s %zu   <-- if this equals the number of SIM HITS there is a\n",
           "simTrackStates_", ev->simTrackStates_.size());
    printf("  %-18s     true state at every hit, and the bias measurement needs no fit\n", "");
    {  // how many hits actually carry usable MC info, and do labels line up?
      long with_mc = 0, foreign = 0, checked = 0;
      for (int lay = 0; lay < (int)ev->layerHits_.size(); ++lay)
        for (const Hit &h : ev->layerHits_[lay]) {
          const int mid = h.mcHitID();
          if (mid >= 0 && mid < (int)ev->simHitsInfo_.size()) ++with_mc;
        }
      for (const Track &t : ev->simTracks_) {
        for (int ih = 0; ih < t.nTotalHits(); ++ih) {
          const int idx = t.getHitIdx(ih), lay = t.getHitLyr(ih);
          if (idx < 0 || lay < 0 || lay >= (int)ev->layerHits_.size()) continue;
          if (idx >= (int)ev->layerHits_[lay].size()) continue;
          const Hit &h = ev->layerHits_[lay][idx];
          const int mid = h.mcHitID();
          if (mid < 0 || mid >= (int)ev->simHitsInfo_.size()) continue;
          ++checked;
          if (ev->simHitsInfo_[mid].mcTrackID() != t.label()) ++foreign;
        }
      }
      printf("  %-18s %ld of %ld hits carry MC info\n", "mc coverage", with_mc, nh);
      printf("  %-18s %ld of %ld hits ON a sim track belong to a DIFFERENT particle"
             " (%.2f%%)\n", "foreign hits", foreign, checked,
             checked ? 100.0 * foreign / checked : 0.0);
    }
    describe("simTracks_", ev->simTracks_);
    describe("seedTracks_", ev->seedTracks_);
    describe("cmsswTracks_", ev->cmsswTracks_);
  }

  // ==========================================================================
  // S11 step 1 -- ENSEMBLE COVARIANCE TRANSPORT, ONE STEP.
  //
  // Draw N samples from C_A, push each through the production propagator, and
  // compare their empirical covariance S against the ONE covariance C_B that
  // the same propagator produces from C_A. No truth, no material, no MC: the
  // only thing under test is the jacobian.
  //
  // Two things this sees that nothing measured so far could:
  //  - a wrong CORRELATION. Six marginal ratios are blind to it; the
  //    eigenvalues of the whitened empirical covariance are not.
  //  - a wrong RANK. Propagating to a surface is a rank-5 operation in 6
  //    dimensions, and a missing surface-crossing term is exactly a rank/null
  //    space defect.
  //
  // The control is `cov_scale`: correct LINEAR transport is scale-invariant, so
  // a scale dependence separates "the jacobian is wrong" from "the step is too
  // nonlinear for any jacobian".
  void val_cov_transport(const char *out_file, int n_samp, unsigned seed) {
    const double kPt[]   = {1.0, 10.0};
    // Extended to 2.5 on 2026-09-17: the measured q under-estimate lives at
    // |eta| > 1.6 and this test stopped at 1.5, so it had never probed the
    // region where the covariance is claimed to be wrong.
    const double kEta[]  = {0.0, 1.5, 2.0, 2.5};
    const double kDal[]  = {0.02, 0.10, 0.40, 0.80, 1.50};
    const double kInc[]  = {1.00, 0.50};
    const float  kScale[] = {1.0f, 100.0f};

    PropagationFlags pf;   // no material, uniform B -- transport only
    std::vector<ValCovXport> out;
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> gauss(0.0, 1.0);

    int cfg = 0;
    for (double pt : kPt) for (double eta : kEta)
    for (double dal : kDal) for (double inc : kInc) for (float scale : kScale) {
      ValCovXport v;
      v.cfg = cfg++;  v.pt = (float)pt;  v.eta = (float)eta;
      v.dalpha = (float)dal;  v.cos_inc = (float)inc;  v.cov_scale = scale;

      ValHelix h;  h.init(pt, eta, -1, 0.3, 0.0, 0.0);

      // State A: the helix's first crossing of r = 30 cm. Plane B: dalpha of
      // turn further along, tilted to the requested incidence.
      double sA = 0.0;
      if (!h.first_crossing_r(30.0, sA)) continue;
      const double sB = h.s_of_alpha(h.alpha_of_s(sA) + dal);
      double pB3[3], mB3[3];  h.at_s(sB, pB3, mB3);
      if (std::hypot(pB3[0], pB3[1]) > kMaxR || std::fabs(pB3[2]) > kMaxZ) continue;
      float plp[3], pln[3];   make_plane(pB3, mB3, inc, plp, pln);

      float parA[6], covA[21];
      h.ccs_at_s(sA, parA);
      nominal_cov(parA, covA);
      for (int i = 0; i < 21; ++i) covA[i] *= scale;

      // --- the reference: one propagation of the mean state with C_A.
      MPlexLV p0, pB;  MPlexLS e0, eB;  MPlexQI q_, fB;  MPlexHV mp, mn;
      fill_lanes(p0, parA);  fill_lanes_sym(e0, covA);
      fill_lanes_q(q_, h.chg);  fill_lanes_hv(mp, plp);  fill_lanes_hv(mn, pln);
      propagateHelixToPlaneMPlex(e0, p0, q_, mp, mn, eB, pB, fB, NN, pf, nullptr);
      if (fB.constAt(0, 0, 0) != 0) { out.push_back(v); continue; }
      float parB[6], covB[21];
      read_lane0(pB, eB, parB, covB);
      {
        const double st = std::sin((double)parB[5]);
        const double ph[3] = {st*std::cos((double)parB[4]), st*std::sin((double)parB[4]),
                              std::cos((double)parB[5])};
        v.cos_inc_meas = (float)std::fabs(ph[0]*pln[0] + ph[1]*pln[1] + ph[2]*pln[2]);
        double cs = 0; for (int i = 0; i < 21; ++i) cs += std::fabs((double)covB[i]);
        v.covb_sum = (float)cs;
      }
      v.d_plane_ref = (float)std::fabs((double)(parB[0]-plp[0])*pln[0] +
                                       (double)(parB[1]-plp[1])*pln[1] +
                                       (double)(parB[2]-plp[2])*pln[2]);

      // --- Cholesky of C_A, in double. C_A is constructed PD, so a failure
      //     here is a bug in nominal_cov, not a result.
      double L[6][6] = {};
      {
        double A[6][6];
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) A[i][j] = covA[ls_idx(i, j)];
        bool ok = true;
        for (int i = 0; i < 6 && ok; ++i) {
          for (int j = 0; j <= i; ++j) {
            double sum = A[i][j];
            for (int k = 0; k < j; ++k) sum -= L[i][k] * L[j][k];
            if (i == j) { if (sum <= 0) { ok = false; break; } L[i][i] = std::sqrt(sum); }
            else L[i][j] = sum / L[j][j];
          }
        }
        if (!ok) { printf("val_cov_transport: cfg %d C_A not PD\n", v.cfg); out.push_back(v); continue; }
      }

      // --- the ensemble. NN samples per propagator call; the covariance
      //     argument is transported too and simply not read.
      std::vector<std::array<double,6>> acc, accv;
      std::vector<double> dpl;
      acc.reserve(n_samp);  accv.reserve(n_samp);  dpl.reserve(n_samp);
      double pref[3];
      { const double st = std::sin((double)parB[5]);
        pref[0] = st*std::cos((double)parB[4]); pref[1] = st*std::sin((double)parB[4]);
        pref[2] = std::cos((double)parB[5]); }
      int n_fail = 0;
      for (int base = 0; base < n_samp; base += NN) {
        const int nlane = std::min(NN, n_samp - base);
        MPlexLV ps;  MPlexLS es;  MPlexQI fs;  MPlexLV po;  MPlexLS eo;
        fill_lanes_sym(es, covA);
        float smp[NN][6];
        for (int n = 0; n < NN; ++n) {
          double z[6];
          for (int i = 0; i < 6; ++i) z[i] = gauss(rng);
          for (int i = 0; i < 6; ++i) {
            double d = 0.0;
            for (int k = 0; k <= i; ++k) d += L[i][k] * z[k];
            smp[n][i] = (float)(parA[i] + d);
          }
          for (int i = 0; i < 6; ++i) ps(n, i, 0) = smp[n][i];
        }
        propagateHelixToPlaneMPlex(es, ps, q_, mp, mn, eo, po, fs, NN, pf, nullptr);
        for (int n = 0; n < nlane; ++n) {
          if (fs.constAt(n, 0, 0) != 0) { ++n_fail; continue; }
          std::array<double,6> a;
          for (int i = 0; i < 6; ++i) a[i] = po.constAt(n, i, 0);
          dpl.push_back(std::fabs((a[0]-plp[0])*pln[0] + (a[1]-plp[1])*pln[1] + (a[2]-plp[2])*pln[2]));
          // phi is squashed into (-pi, pi] by the propagator; unwrap it against
          // the reference before averaging, or a track near the branch cut
          // manufactures a variance of order pi^2.
          double dphi = a[4] - parB[4];
          while (dphi >  M_PI) dphi -= 2.0 * M_PI;
          while (dphi <= -M_PI) dphi += 2.0 * M_PI;
          a[4] = parB[4] + dphi;
          acc.push_back(a);
          // Curvilinear projection: slide the sample along ITS OWN momentum
          // until it reaches the plane through the reference point normal to
          // the REFERENCE momentum. Momentum is unchanged by the slide, so only
          // the three position components move.
          {
            const double st = std::sin(a[5]);
            const double ph[3] = {st*std::cos(a[4]), st*std::sin(a[4]), std::cos(a[5])};
            const double den = ph[0]*pref[0] + ph[1]*pref[1] + ph[2]*pref[2];
            std::array<double,6> b = a;
            if (std::fabs(den) > 1e-6) {
              const double tt = (((double)parB[0]-a[0])*pref[0] +
                                 ((double)parB[1]-a[1])*pref[1] +
                                 ((double)parB[2]-a[2])*pref[2]) / den;
              for (int i = 0; i < 3; ++i) b[i] = a[i] + tt*ph[i];
            }
            accv.push_back(b);
          }
        }
      }
      v.n_samp = (int)acc.size();
      v.n_fail = n_fail;
      if (!dpl.empty()) {
        std::sort(dpl.begin(), dpl.end());
        v.d_plane_med = (float)dpl[dpl.size()/2];
        v.d_plane_p90 = (float)dpl[(size_t)(0.90*(dpl.size()-1))];
        v.d_plane_max = (float)dpl.back();
      }
      if (v.n_samp < 100) { out.push_back(v); continue; }

      // --- empirical mean and covariance
      double mean[6] = {};
      for (auto &a : acc) for (int i = 0; i < 6; ++i) mean[i] += a[i];
      for (int i = 0; i < 6; ++i) mean[i] /= v.n_samp;
      double S[6][6] = {};
      for (auto &a : acc)
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j)
          S[i][j] += (a[i] - mean[i]) * (a[j] - mean[j]);
      for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) S[i][j] /= (v.n_samp - 1);

      for (int i = 0; i < 6; ++i) {
        const double cb = covB[ls_idx(i, i)];
        v.ratio[i]     = cb > 0 ? (float)std::sqrt(S[i][i] / cb) : -1.f;
        v.mean_pull[i] = cb > 0 ? (float)((mean[i] - parB[i]) / std::sqrt(cb)) : 0.f;
      }

      // --- whiten inside C_B's column space.
      //
      // NON-DIMENSIONALISE FIRST. C_B mixes cm^2 (position), rad^2 (two angles)
      // and GeV^-2 (1/pT), so its raw eigenvalues are not comparable to one
      // another and "lambda > eps * lambda_max" is not a rank test at all -- it
      // is a statement about the units the state happens to be written in.
      // Scaling each coordinate by its own sigma gives the CORRELATION matrix:
      // dimensionless, unit diagonal, eigenvalues summing to 6. The rank gap is
      // then a real gap. (Getting this wrong makes the whitening divide by
      // float32 dust and reports lam_max of 1e9 -- a measurement artefact, not
      // a propagator defect.)
      double sc[6];
      for (int i = 0; i < 6; ++i) {
        const double cb = covB[ls_idx(i, i)];
        sc[i] = cb > 0 ? 1.0 / std::sqrt(cb) : 0.0;
      }
      TMatrixDSym R(6), Sn(6);
      for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
        R(i, j)  = covB[ls_idx(i, j)] * sc[i] * sc[j];
        Sn(i, j) = S[i][j]            * sc[i] * sc[j];
      }
      TMatrixDSymEigen eig(R);
      TVectorD d = eig.GetEigenValues();     // descending, dimensionless
      TMatrixD U = eig.GetEigenVectors();
      for (int i = 0; i < 6; ++i) v.cb_spec[i] = (float)d[i];
      // C_B is held in float32, so a structurally-zero eigenvalue comes back at
      // ~1e-7 of the unit scale. 1e-5 sits two decades above that dust and many
      // decades below any real direction.
      int rank = 0;
      while (rank < 6 && d[rank] > 1e-5) ++rank;
      v.rank_cb = rank;
      if (rank == 0) { out.push_back(v); continue; }

      TMatrixD W(rank, 6);
      for (int r = 0; r < rank; ++r) {
        const double f = 1.0 / std::sqrt(d[r]);
        for (int i = 0; i < 6; ++i) W(r, i) = f * U(i, r);
      }
      TMatrixD WS(W, TMatrixD::kMult, Sn);
      TMatrixD Wt(TMatrixD::kTransposed, W);
      TMatrixD Z(WS, TMatrixD::kMult, Wt);
      TMatrixDSym Zs(rank);
      for (int i = 0; i < rank; ++i) for (int j = 0; j < rank; ++j) Zs(i, j) = 0.5 * (Z(i, j) + Z(j, i));
      TMatrixDSymEigen ze(Zs);
      TVectorD zl = ze.GetEigenValues();
      v.lam_max = (float)zl[0];  v.lam_min = (float)zl[rank - 1];
      for (int i = 0; i < rank && i < 6; ++i) v.lam[i] = (float)zl[i];

      // --- the same comparison on the curvilinear surface. Off normal
      //     incidence this is the only fair one; at cos_inc = 1 it must agree
      //     with the block above, which is its own cross-check.
      {
        double mv[6] = {};
        for (auto &a : accv) for (int i = 0; i < 6; ++i) mv[i] += a[i];
        for (int i = 0; i < 6; ++i) mv[i] /= accv.size();
        TMatrixDSym Sv(6);
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) Sv(i, j) = 0.0;
        for (auto &a : accv)
          for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j)
            Sv(i, j) += (a[i]-mv[i])*(a[j]-mv[j]);
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j)
          Sv(i, j) *= sc[i]*sc[j] / (accv.size()-1);
        TMatrixD Zv(TMatrixD(W, TMatrixD::kMult, Sv), TMatrixD::kMult, Wt);
        TMatrixDSym Zvs(rank);
        for (int i = 0; i < rank; ++i) for (int j = 0; j < rank; ++j)
          Zvs(i, j) = 0.5*(Zv(i, j) + Zv(j, i));
        TVectorD vl = TMatrixDSymEigen(Zvs).GetEigenValues();
        v.rank_cv = rank;
        v.lam_max_cv = (float)vl[0];  v.lam_min_cv = (float)vl[rank-1];
      }

      // --- the null direction: the samples' spread along the direction C_B
      //     says has no width, in units of the smallest direction it does keep.
      if (rank < 6) {
        double nv = 0.0;
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) nv += U(i, rank) * Sn(i, j) * U(j, rank);
        v.null_sig = (float)(std::sqrt(std::max(nv, 0.0)) / std::sqrt(d[rank - 1]));
      }
      // --- FULL CHAIN. Offer each sample's landing point to the production
      //     chi2 as a hit with a negligible error, against the reference state
      //     and its transported covariance. This runs jacCurv2Loc and the local
      //     projection, so unlike the bare comparison above it is well posed at
      //     any incidence.
      {
        float pdir[3];
        {   // any in-plane direction will do: chi2 is invariant under rotation
            // about the normal (verified independently as the suite's finding 11)
          double a[3] = {0.0, 0.0, 1.0};
          double d[3] = {a[1]*pln[2]-a[2]*pln[1], a[2]*pln[0]-a[0]*pln[2], a[0]*pln[1]-a[1]*pln[0]};
          double dn = std::sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
          if (dn < 1e-6) { d[0]=1; d[1]=0; d[2]=0; dn=1; }
          for (int i = 0; i < 3; ++i) pdir[i] = (float)(d[i]/dn);
        }
        MPlexLS pe;  MPlexLV pp;  MPlexHS me;  MPlexHV mv, nn, dd, pt3;  MPlexQF c2;
        fill_lanes(pp, parB);  fill_lanes_sym(pe, covB);
        fill_lanes_hv(nn, pln);  fill_lanes_hv(dd, pdir);  fill_lanes_hv(pt3, plp);
        const float tiny = 1e-10f;   // (0.1 um)^2 -- present only to keep the
                                     // 2-D matrix invertible, 4 orders below
                                     // any track covariance here
        for (int n = 0; n < NN; ++n)
          for (int i = 0; i < 6; ++i) me.fArray[i*NN + n] = (i==0||i==2||i==5) ? tiny : 0.f;
        std::vector<double> c2v;  c2v.reserve(acc.size());
        double c2sum = 0.0;
        for (size_t base = 0; base < acc.size(); base += NN) {
          const size_t nlane = std::min((size_t)NN, acc.size() - base);
          for (int n = 0; n < NN; ++n) {
            const auto &a = acc[std::min(base + n, acc.size() - 1)];
            for (int i = 0; i < 3; ++i) mv(n, i, 0) = (float)a[i];
          }
          kalmanComputeChi2Plane(pe, pp, q_, me, mv, nn, dd, pt3, c2, NN);
          for (size_t n = 0; n < nlane; ++n) {
            const float c = c2.constAt(n, 0, 0);
            // NOT std::isfinite: -Ofast implies -ffinite-math-only, under
            // which it returns true for an inf. isFinite is bit-based.
            if (!isFinite(c)) continue;
            c2v.push_back(c);  c2sum += c;
          }
        }
        if (c2v.size() > 100) {
          std::sort(c2v.begin(), c2v.end());
          v.chi2_med  = (float)c2v[c2v.size()/2];
          v.chi2_p90  = (float)c2v[(size_t)(0.90*(c2v.size()-1))];
          v.chi2_p99  = (float)c2v[(size_t)(0.99*(c2v.size()-1))];
          v.chi2_mean = (float)(c2sum / c2v.size());
          v.chi2_med_ratio = v.chi2_med / 1.38629436f;
        }
      }

      out.push_back(v);
    }

    TFile f(out_file, "RECREATE");
    TTree *t = new TTree("covxport", "ensemble covariance transport, one step");
    ValCovXport v;  t->Branch("v", &v);
    for (const auto &x : out) { v = x; t->Fill(); }
    t->Write();  f.Close();
    printf("val_cov_transport: %zu configs, %d samples each -> '%s'\n",
           out.size(), n_samp, out_file);
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
    std::vector<ValSearchHit> g_sh;
    std::vector<ValSearchMiss> g_sm;
    std::vector<ValCovStep>    g_cs;

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

  // Truth association for val_eff and anything else on setMCTrackIDInfo. false
  // is quality-val's rule, 2*mccount >= nCandHits over the non-seed hits. true
  // is CMSSW MTV's quickTrackAssociatorByHits: shared / valid reco hits > 0.75,
  // strictly, over all hits, seed hits included. Only the association: the
  // --mtv-like-val option also changes the sim selection, this does not.
  //
  // MTV has no minimum track length, so MTV mode also sets nMinFoundHits = 0
  // (as --mtv-like-val does): otherwise modifyRefTrackID() turns a correctly
  // matched track with fewer than 10 non-seed hits into code -2, which
  // val_eff then counted as a fake -- and a finder that builds longer tracks
  // looked cleaner for it. The previous value is restored on switching off.
  void val_assoc_mtv(bool on) {
    static int saved_min_found_hits = Config::nMinFoundHits;
    if (on && !Config::mtvLikeValidation)
      saved_min_found_hits = Config::nMinFoundHits;
    Config::mtvLikeValidation = on;
    Config::nMinFoundHits = on ? 0 : saved_min_found_hits;
    printf("val_assoc_mtv: %s, nMinFoundHits %d\n",
           on ? "MTV, 4*mccount > 3*nCandHits over all hits" : "quality-val, 2*mccount >= nCandHits over non-seed hits",
           Config::nMinFoundHits);
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

  void v2p2_final_beam_diag_reset();
  void v2p2_final_beam_diag_report();
  void val_final_beam(bool on) {
    V2p2::Diag::final_beam_purity = on;
    v2p2_final_beam_diag_reset();
    printf("val_final_beam: %d (end-of-search beam purity, reset)\n", (int) on);
  }
  void val_final_beam_report() { v2p2_final_beam_diag_report(); }

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

  // Mean ln(rho) and mean ln(det V) per hit taken, so the constants for
  // val_score_terms() are measured. Accumulation is not thread safe; a trace build serialises the
  // in-event loops, which is where this is meant to run.
  void val_score_term_stats(bool on) {
    if (on) {
      V2p2::ScoreStats::n_hits = 0;
      V2p2::ScoreStats::sum_log_rho = V2p2::ScoreStats::sum_log_detv = 0.0;
      V2p2::ScoreStats::accum = true;
      printf("val_score_term_stats: accumulating.\n");
      return;
    }
    V2p2::ScoreStats::accum = false;
    const long n = V2p2::ScoreStats::n_hits;
    printf("val_score_term_stats: %ld hits scored; mean ln(rho) = %.4f, "
           "mean ln(det V) = %.4f\n", n,
           n ? V2p2::ScoreStats::sum_log_rho / n : 0.0,
           n ? V2p2::ScoreStats::sum_log_detv / n : 0.0);
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

  void val_kinfo_reset() { g_ki = KIAcc(); }

  void val_kinfo_event(const Event *ev) {
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
    printf("val_kinfo_event: needs MKFIT_TRACE_KALMAN_DEBUG\n");
#endif
  }

  void val_kinfo_report() {
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

  // bit 0 = skip the exact double-precision helix-plane solve (the expensive
  // part); bit 1 = keep only pixel-barrel hits. 3 = both, which is what a
  // pixel-barrel covariance scan wants; 1 = whole detector but still fast,
  // which is what a phi-everywhere measurement wants.
  void val_search_lite(int mode) {
    g_val_search_lite = mode;
    printf("val_search_lite: mode %d (skip_solve=%d pixb_only=%d)\n",
           mode, mode & 1, (mode >> 1) & 1);
  }


  // ==========================================================================
  // Track-level true-hit efficiency of the inward search: of the sim track's
  // hits in the pixel barrel and the forward disks, how many did the FINAL
  // candidate actually collect? This is the quantity "the backward search gets
  // true hits down to the first pixel layers with 75% efficiency" refers to,
  // and unlike any per-layer-search measure it is not survivorship-biased: a
  // candidate that dies early simply has few matched hits.
  namespace {
    struct TEff { long n_trk = 0;
                  long sim_pix = 0, got_pix = 0, sim_dsk = 0, got_dsk = 0;
                  long deepest_ok = 0;
                  // "any true pixel-barrel hit" is only a fair question for a
                  // track that HAS one to find. High-|eta| tracks largely miss
                  // the pixel barrel (~1 hit each against ~3 for barrel-only),
                  // so counting them in the denominator understates the rate.
                  long n_trk_with_pix = 0;
                  // Distinct LAYERS with a sim hit, against the hit count. The
                  // search took at most one hit per layer before the in-layer
                  // search, so layers/hits was a hard ceiling on the efficiency above --
                  // and phase-2 disks carry module overlaps, so it is well below 1.
                  long sim_pix_lay = 0, sim_dsk_lay = 0; };
    TEff g_te[2];   // [0] barrel-only sim track, [1] touches the disks
    // Split by how many pixel-barrel hits the SIM track offers. A plain T5 is a
    // track whose pixel match FAILED upstream (pT5 = LST T5 + Patatrack pixels),
    // so the sharp question is: when the pixel hits are genuinely there, does
    // the inward search get them?
    TEff g_tp[3];   // spix 0 / 1-2 / >=3
    long g_te_ntrk = 0, g_te_nolbl = 0;
  }

  void val_track_eff_reset() { g_te[0] = g_te[1] = TEff();
                               g_tp[0] = g_tp[1] = g_tp[2] = TEff();
                               g_te_ntrk = g_te_nolbl = 0; }

  void val_track_eff_event(const Event *ev) {
    if (ev == nullptr) return;
    for (const Track &c : ev->candidateTracks_) {
      ++g_te_ntrk;
      auto si = ev->simInfoForTrack(c);
      if (!si.is_set()) { ++g_te_nolbl; continue; }
      const int L = si.label;
      const Track &st = ev->simTracks_[L];
      // the sim track's own hit content, by region
      long spix = 0, sdsk = 0;
      std::set<int> lpix, ldsk;
      for (int i = 0; i < st.nTotalHits(); ++i) {
        const HitOnTrack hot = st.getHitOnTrack(i);
        if (hot.index < 0) continue;
        if (hot.layer < 4) { ++spix; lpix.insert(hot.layer); }
        else if (hot.layer >= 16) { ++sdsk; ldsk.insert(hot.layer); }
      }
      if (spix == 0 && sdsk == 0) continue;
      const int cls = (sdsk > 0) ? 1 : 0;
      // what the candidate actually holds there, TRUTH-checked per hit
      long gpix = 0, gdsk = 0; int deepest = -1;
      for (int i = 0; i < c.nTotalHits(); ++i) {
        const HitOnTrack hot = c.getHitOnTrack(i);
        if (hot.index < 0 || hot.layer < 0) continue;
        const Hit &h = ev->layerHits_[hot.layer][hot.index];
        const int hl = ev->simHitsInfo_[h.mcHitID()].mcTrackID();
        if (hl != L) continue;
        if (hot.layer < 4)  { ++gpix; if (hot.layer > deepest) deepest = hot.layer; }
        else if (hot.layer >= 16) ++gdsk;
      }
      TEff &tp = g_tp[spix == 0 ? 0 : (spix < 3 ? 1 : 2)];
      ++tp.n_trk;  tp.sim_pix += spix;  tp.got_pix += gpix;
      tp.sim_dsk += sdsk;  tp.got_dsk += gdsk;
      if (spix > 0) { ++tp.n_trk_with_pix; if (gpix > 0) ++tp.deepest_ok; }
      tp.sim_pix_lay += (long) lpix.size();

      TEff &t = g_te[cls];
      ++t.n_trk;  t.sim_pix += spix;  t.got_pix += gpix;
      if (spix > 0) ++t.n_trk_with_pix;
      t.sim_pix_lay += (long) lpix.size();  t.sim_dsk_lay += (long) ldsk.size();
      t.sim_dsk += sdsk;  t.got_dsk += gdsk;
      if (gpix > 0) ++t.deepest_ok;
    }
  }

  void val_track_eff_report(const char *tag) {
    printf("\n=== TRACK-LEVEL TRUE-HIT EFFICIENCY [%s] ===\n", tag);
    printf("%ld candidate tracks, %ld with no sim label.\n", g_te_ntrk, g_te_nolbl);
    printf("  %-14s %7s | %9s %9s %7s | %9s %9s %7s | %8s %9s\n",
           "sim track","tracks","sim pixB","got","eff","sim disks","got","eff",
           "w/ pixB","any pixB*");
    printf("  (* fraction of tracks that HAVE a pixel-barrel hit and got at least one)\n");
    printf("  CEILING = distinct layers / hits. It was a hard ceiling while the search\n");
    printf("  added at most ONE hit per layer (the best-hit path); with the in-layer\n");
    printf("  combinatorial search it is a soft one -- a path may take\n");
    printf("  several hits in a layer, so read these as a FRACTION OF CEILING, not against\n");
    printf("  100%%. Overlaps are what put the ceiling well below 100%% in the first place.\n");
    const char *nm[2] = {"barrel-only", "touches disks"};
    for (int c = 0; c < 2; ++c) {
      const TEff &t = g_te[c];
      if (t.n_trk == 0) continue;
      printf("  %-14s %7ld | %9ld %9ld %6.1f%% | %9ld %9ld %6.1f%% | %8ld %8.1f%%\n",
             nm[c], t.n_trk, t.sim_pix, t.got_pix,
             t.sim_pix ? 100.0*t.got_pix/t.sim_pix : 0.0,
             t.sim_dsk, t.got_dsk, t.sim_dsk ? 100.0*t.got_dsk/t.sim_dsk : 0.0,
             t.n_trk_with_pix,
             t.n_trk_with_pix ? 100.0*t.deepest_ok/t.n_trk_with_pix : 0.0);
      printf("  %-14s %7s | CEILING pixB %5.1f%% (of %ld) | CEILING disks %5.1f%% (of %ld)\n",
             "", "", t.sim_pix ? 100.0*t.sim_pix_lay/t.sim_pix : 0.0, t.sim_pix_lay,
             t.sim_dsk ? 100.0*t.sim_dsk_lay/t.sim_dsk : 0.0, t.sim_dsk_lay);
    }
    printf("\n  SPLIT BY HOW MANY PIXEL-BARREL HITS THE SIM TRACK OFFERS:\n");
    printf("  %-14s %7s | %8s %8s %7s %9s | %9s\n",
           "sim pixB hits","tracks","sim","got","eff","of ceiling","any pixB");
    const char *pn[3] = {"0 (none to find)", "1 - 2", ">= 3"};
    for (int c = 0; c < 3; ++c) {
      const TEff &t = g_tp[c];
      if (t.n_trk == 0) continue;
      const double eff = t.sim_pix ? 100.0*t.got_pix/t.sim_pix : 0.0;
      const double cei = t.sim_pix ? 100.0*t.sim_pix_lay/t.sim_pix : 0.0;
      printf("  %-14s %7ld | %8ld %8ld %6.1f%% %8.0f%% | %8.1f%%\n",
             pn[c], t.n_trk, t.sim_pix, t.got_pix, eff, cei > 0 ? 100.0*eff/cei : 0.0,
             t.n_trk_with_pix ? 100.0*t.deepest_ok/t.n_trk_with_pix : 0.0);
    }
  }

  void val_search_reset() { g_sh.clear(); g_sm.clear(); g_cs.clear(); g_sh_evt = -1; }

  void val_search_event(const Event *ev, int event_idx) {
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
      if ((g_val_search_lite & 2) && (hm.layer < 0 || hm.layer >= 4))
        continue;                                   // pixel barrel only
      ValSearchHit v;
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
      if (!(g_val_search_lite & 1) &&
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
      ValCovStep c;
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
      ValSearchMiss m;
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

    printf("val_search_event: ev %d -- %d hit matches, %d searches, %d kalman updates\n",
           event_idx, n_add, n_ls, n_ku);
  }

  void val_search_write(const char *out_file) {
    TFile f(out_file, "RECREATE");
    TTree *t = new TTree("search", "backward-search hit matches");
    ValSearchHit v;  t->Branch("h", &v);
    for (const auto &x : g_sh) { v = x; t->Fill(); }
    TTree *t3 = new TTree("cov", "one row per Kalman update: covariance before/after");
    ValCovStep c;  t3->Branch("c", &c);
    for (const auto &x : g_cs) { c = x; t3->Fill(); }

    TTree *t2 = new TTree("miss", "one row per layer-search: where was the sim hit");
    ValSearchMiss m;  t2->Branch("m", &m);
    for (const auto &x : g_sm) { m = x; t2->Fill(); }
    t->Write();  t2->Write();  t3->Write();  f.Close();
    printf("val_search_write: %zu hit matches, %zu searches, %zu kalman updates -> '%s'\n",
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

  void val_bkfit_trace_reset() { g_bt.clear(); }

  void val_bkfit_trace_event(const Event *ev) {
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

  void val_bkfit_trace_report(const char *prefix) {
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



  // Toggle PF_apply_material on the backward-fit propagation only. Returns the
  // previous state. Diagnostic: material inflates the propagated covariance, so
  // turning it off must RAISE chi2; how much says how big the material term is.
  void val_bkfit_material(bool on) {
    auto &pc = const_cast<PropagationConfig&>(Config::TrkInfo.prop_config());
    pc.backward_fit_pflags.apply_material = on;
    printf("val_bkfit_material: backward_fit_pflags.apply_material = %d\n", (int) on);
  }

}  // namespace mkfit
