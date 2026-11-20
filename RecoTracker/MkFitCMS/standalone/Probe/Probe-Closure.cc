// ===========================================================================
// Probe-Closure.cc -- CLOSURE probes: synthetic tracks with an exact answer.
//   probe_gen_closure     one step there and back through the plane propagator
//   probe_cov_transport   ensemble test of the covariance transport
//   probe_dalpha_*        the turn angles between consecutive hits of real
//                         tracks, which weight the two scans above; it reads
//                         events
//
// No detector geometry is used. Closure is a property of the PROPAGATOR, not
// of the detector, so the planes here are constructed analytically at chosen
// radii. That also removes a whole class of contamination -- a synthetic plane
// cannot be geometry-incompatible with a sample.
//
// GENERATION needs the mkFit build (the propagator). ANALYSIS does not: the
// output is a plain TTree of dictionary-ed Probe* structs, which RDataFrame
// reads in a bare root.exe session after
//     gSystem->Load("libMkFitRootDataFormats.so")
// ===========================================================================
#include "RecoTracker/MkFitCMS/standalone/Probe/Probe.h"
#include "RecoTracker/MkFitCore/standalone/DataFormats/ProbeStructs.h"

#include "RecoTracker/MkFitCore/interface/Config.h"
#include "RecoTracker/MkFitCore/interface/cms_common_macros.h"
#include "RecoTracker/MkFitCore/interface/PropagationConfig.h"
#include "RecoTracker/MkFitCore/interface/Track.h"
#include "RecoTracker/MkFitCore/src/Matrix.h"
#include "RecoTracker/MkFitCore/src/PropagationMPlex.h"
#include "RecoTracker/MkFitCore/src/KalmanUtilsMPlex.h"
#include "RecoTracker/MkFitCore/standalone/Event.h"

#include "TFile.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

namespace mkfit {

  namespace {

    // ----------------------------------------------------------------------
    // Exact uniform-B helix from a production vertex, parametrized by signed
    // 3D path length s. This is the TRUTH: it defines where the planes are,
    // and nothing in it comes from the code under test.

    struct ProbeHelix {
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

  void probe_gen_closure(const char *out_file, int n_pt, int n_eta, int n_phi,
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
    ProbeCfgInfo ci;
    t_cfg->Branch("g", &ci);
    for (int ia = 0; ia < kNDal; ++ia)
      for (int ic = 0; ic < kNInc; ++ic) {
        ci.cfg = ia * kNInc + ic;
        ci.dalpha = (float)kDalpha[ia];
        ci.cos_inc = (float)kCosInc[ic];
        t_cfg->Fill();
      }

    TTree *t = new TTree("closure", "one-step-there-and-back propagation closure");
    ProbeClosure c;
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

            ProbeHelix h;
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

                c = ProbeClosure{};
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
                // Keep the exact crossing in double -- see ProbeStructs.h.
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

    printf("probe_gen_closure: %d tracks x %d cfgs = %ld trials, %ld outside the tracker\n"
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
    std::vector<ProbeStep> g_probe_steps;
    int g_probe_evt = 0;
  }

  void probe_dalpha_reset() { g_probe_steps.clear(); g_probe_evt = 0; }

  void probe_dalpha_add_event(const Event *ev, int which, float pt_min) {
    if (ev == nullptr) { printf("probe_dalpha_add_event: null event\n"); return; }
    const TrackVec &tv = (which == 1) ? ev->simTracks_
                       : (which == 2) ? ev->candidateTracks_
                                      : ev->cmsswTracks_;
    // R_c [cm] = pT / (0.01 * sol * B)
    const double inv_k = 0.01 * (double)Const::sol * Config::Bfield;
    const int evt = g_probe_evt++;

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
          ProbeStep vs;
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
          g_probe_steps.push_back(vs);
        }
        prev_lay = lay; prev_idx = idx;
      }
    }
  }

  void probe_dalpha_write(const char *out_file) {
    TFile f(out_file, "RECREATE");
    TTree *t = new TTree("steps", "turn angle between consecutive hits of real tracks");
    ProbeStep vs;
    t->Branch("s", &vs);
    for (const ProbeStep &v : g_probe_steps) { vs = v; t->Fill(); }
    t->Write();
    f.Close();
    printf("probe_dalpha_write: %zu steps from %d events -> '%s'\n",
           g_probe_steps.size(), g_probe_evt, out_file);
  }

  // ==========================================================================
  // ENSEMBLE COVARIANCE TRANSPORT, ONE STEP.
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
  void probe_cov_transport(const char *out_file, int n_samp, unsigned seed) {
    const double kPt[]   = {1.0, 10.0};
    // Extended to 2.5 on 2026-09-17: the measured q under-estimate lives at
    // |eta| > 1.6 and this test stopped at 1.5, so it had never probed the
    // region where the covariance is claimed to be wrong.
    const double kEta[]  = {0.0, 1.5, 2.0, 2.5};
    const double kDal[]  = {0.02, 0.10, 0.40, 0.80, 1.50};
    const double kInc[]  = {1.00, 0.50};
    const float  kScale[] = {1.0f, 100.0f};

    PropagationFlags pf;   // no material, uniform B -- transport only
    std::vector<ProbeCovXport> out;
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> gauss(0.0, 1.0);

    int cfg = 0;
    for (double pt : kPt) for (double eta : kEta)
    for (double dal : kDal) for (double inc : kInc) for (float scale : kScale) {
      ProbeCovXport v;
      v.cfg = cfg++;  v.pt = (float)pt;  v.eta = (float)eta;
      v.dalpha = (float)dal;  v.cos_inc = (float)inc;  v.cov_scale = scale;

      ProbeHelix h;  h.init(pt, eta, -1, 0.3, 0.0, 0.0);

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
        if (!ok) { printf("probe_cov_transport: cfg %d C_A not PD\n", v.cfg); out.push_back(v); continue; }
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
    ProbeCovXport v;  t->Branch("v", &v);
    for (const auto &x : out) { v = x; t->Fill(); }
    t->Write();  f.Close();
    printf("probe_cov_transport: %zu configs, %d samples each -> '%s'\n",
           out.size(), n_samp, out_file);
  }

}  // namespace mkfit
