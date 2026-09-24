#ifndef mkfit_seeding_SeedMargins_h
#define mkfit_seeding_SeedMargins_h

// Every cut of one quad, evaluated in one arithmetic, with its signed margin.
//
// This is the difference tool's core.  When two arithmetics give different
// quad lists, each differing quad is evaluated in both: the cuts whose pass
// flag differs are the ones that flipped, and their margin in the reference
// arithmetic says how close to the threshold the quad was.  A flip within a
// stated epsilon of the threshold is rounding; a flip far from it is a bug.
//
// The pair cuts (a-b, b-c) are the b-major finder's expressions; the triplet
// and quad cuts call helix_triplet<A>() and quad_cuts<A>(), the functions
// finish_triplets<A>() itself runs, so they cannot drift.

#include "SeedFinder.h"
#include "SeedMath.h"

namespace mkfit::seeding {

  enum MarginCut {
    MC_ab_r,
    MC_ab_phi,
    MC_bc_r,
    MC_bc_phi,
    MC_bc_lin,
    MC_c_z,
    MC_c3,
    MC_reach,
    MC_dfetch_phi,  // informational: the fetch rounds outward to whole bins
    MC_dfetch_z,    // informational
    MC_cd_r,
    MC_d_cross,
    MC_d_phi,
    MC_d_z,
    MC_N
  };

  struct MarginCutInfo {
    const char *name;
    bool rad;       // unit: rad, else cm
    bool decisive;  // part of the accept decision
  };

  inline const MarginCutInfo &margin_cut_info(int i) {
    static const MarginCutInfo info[MC_N] = {
        {"ab_r", false, true},   {"ab_phi", true, true},   {"bc_r", false, true},     {"bc_phi", true, true},
        {"bc_lin", true, true},  {"c_z", false, true},     {"c3", false, true},       {"reach", false, true},
        {"dfetch_phi", true, false}, {"dfetch_z", false, false}, {"cd_r", false, true}, {"d_cross", false, true},
        {"d_phi", true, true},   {"d_z", false, true}};
    return info[i];
  }

  struct QuadEval {
    CutVal c[MC_N];
    double R = 0, cot = 0;  // for the printout: pT = 0.0114 R [GeV, cm], and cot(theta)
    bool pass() const {
      for (int i = 0; i < MC_N; ++i)
        if (margin_cut_info(i).decisive && c[i].valid && !c[i].pass)
          return false;
      // a decisive cut left invalid means an earlier one failed
      for (int i = 0; i < MC_N; ++i)
        if (margin_cut_info(i).decisive && !c[i].valid && i != MC_bc_phi && i != MC_bc_lin)
          return false;
      return true;
    }
  };

  // ka..kd are BIN-ORDER indices in the four layers.
  template <class A, typename LA, typename LB, typename LC, typename LD>
  void eval_quad(const SeedParams &P,
                 const LA &ga,
                 const LB &gb,
                 const LC &gc,
                 const LD &gd,
                 unsigned int ka,
                 unsigned int kb,
                 unsigned int kc,
                 unsigned int kd,
                 QuadEval &e) {
    using detail::wrap_pi;
    using T = typename A::real;
    e = QuadEval();
    constexpr float kBfield = 3.8f;
    const double Rmin = P.pt_min / (0.003f * kBfield);
    const float inv2R = (float)(1.0 / (2 * Rmin)), d0m = P.d0_max, marg = P.phi_margin;
    auto wphi_w = [=](float r_in, float inv_in, float r_out, float inv_out) {
      return (r_out - r_in) * inv2R + d0m * (inv_in - inv_out) + marg;
    };
    const float pa = ga.phi_[ka], za = ga.z_[ka], rra = ga.r_[ka], inva = ga.invr_[ka];
    const float pb = gb.phi_[kb], zb = gb.z_[kb], rrb = gb.r_[kb], invb = gb.invr_[kb];
    const float pc = gc.phi_[kc], zc = gc.z_[kc], rrc = gc.r_[kc], invc = gc.invr_[kc];

    //---- a-b
    e.c[MC_ab_r] = {double(rrb - (rra + 0.1f)), !(rrb <= rra + 0.1f), true};
    {
      const float w = wphi_w(rra, inva, rrb, invb), x = std::abs(wrap_pi(pb - pa));
      e.c[MC_ab_phi] = {double(w) - double(x), !(x > w), true};
    }
    //---- b-c
    e.c[MC_bc_r] = {double(rrc - (rrb + 0.1f)), !(rrc <= rrb + 0.1f), true};
    if (P.phi_lin != 1) {
      const float w = wphi_w(rrb, invb, rrc, invc), x = std::abs(wrap_pi(pc - pb));
      e.c[MC_bc_phi] = {double(w) - double(x), !(x > w), true};
    }
    if (P.phi_lin) {
      const float slope = wrap_pi(pb - pa) / (rrb - rra);
      const float kab = 1.0f / (rrb - rra);
      const float dinv_ab = invb - inva;
      const float wl = d0m * std::abs(invc - invb - (rrc - rrb) * kab * dinv_ab) + P.phi_lin_marg;
      const float x = std::abs(wrap_pi(pc - pb - slope * (rrc - rrb)));
      e.c[MC_bc_lin] = {double(wl) - double(x), !(x > wl), true};
    }
    {
      const double cot = (zb - za) / (rrb - rra);
      const double zpred = za + cot * (rrc - rra);
      const double dz = zc - zpred;
      e.c[MC_c_z] = {double(P.qwin) - std::abs(dz), !(std::abs(dz) > P.qwin), true};
      e.cot = cot;
    }
    //---- triplet helix and the 4th-layer window
    TripletHelix<A> h;
    helix_triplet<A>(P.qwin_d, P.phiwin_d, ga.x_[ka], ga.y_[ka], za, gb.x_[kb], gb.y_[kb], gc.x_[kc], gc.y_[kc], zc,
                     T(gd.rlo_), T(gd.rhi_), h);
    e.R = h.R;
    e.c[MC_c3] = h.c3;
    e.c[MC_reach] = h.reach;
    if (h.ok) {
      const float x = std::abs(wrap_pi((float)(gd.phi_[kd] - (float)h.phid)));
      e.c[MC_dfetch_phi] = {double(h.phidh) - double(x), !(x > h.phidh), true};
      const double zd = gd.z_[kd];
      const double m = std::min(zd - double(h.zlo), double(h.zhi) - zd);
      e.c[MC_dfetch_z] = {m, m >= 0, true};
    }
    //---- quad
    QuadCuts qc;
    quad_cuts<A>(P.qwin_d, P.phiwin_d, h, zc, rrc, gd.r_[kd], gd.phi_[kd], gd.z_[kd], &qc);
    e.c[MC_cd_r] = qc.cd_r;
    e.c[MC_d_cross] = qc.d_cross;
    e.c[MC_d_phi] = qc.d_phi;
    e.c[MC_d_z] = qc.d_z;
  }

}  // namespace mkfit::seeding

#endif
