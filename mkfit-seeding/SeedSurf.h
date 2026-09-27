#ifndef mkfit_seeding_SeedSurf_h
#define mkfit_seeding_SeedSurf_h

// Quadruplet seeds on any four pixel surfaces, barrel cylinders or discs:
// the scalar REFERENCE finder for the patterns beyond the pixel barrel
// (B1 B2 B3 F1, B1 B2 F1 F2, B1 F1 F2 F3, F_k .. F_k+3).  Double precision,
// written for being right and readable, not fast; the kernel finder for the
// barrel (SeedFinderTile.h) is where speed comes from, later.
//
// Every prediction is made in the TARGET layer's own coordinates: a barrel
// layer is parametrised by r (qbar) and measures z (q), a disc by z (qbar) and
// measures r (q).  A hit is predicted at its own qbar, so the layer's
// thickness never widens a window.
//
//   b, doublet   |dphi| < (r_b - r_a) / 2 R_min + D0 (1/r_a - 1/r_b) + marg_b,
//                the geometric bound of the barrel finder; and the r-z line
//                through a and b reaches r = 0 within zv of the beam spot z
//   c, triplet   the a-b line in (u, q) and (u, phi), with u the target's qbar.
//                For D0 = 0 the azimuth of a point is phi0 + s / 2R and z is
//                linear in s, so phi is EXACTLY linear in z: on a disc the
//                phi extrapolation has no curvature term.  r is linear in s
//                up to s^3 / 24 R^2.  Barrel targets use u = r as the barrel
//                finder does.  Tolerance: D0 times the departure of 1/r from
//                linear in u (what D0 adds to phi), plus phi_c; and q_c.
//   d, quad      the circle through a, b, c, in its curvature form (signed k,
//                point and tangent at c), with z linear in arc length.  Disc:
//                move along the arc to the hit's z, compare r and phi.
//                Barrel: solve |P(s)| = r_hit by Newton, compare z and phi.
//                Tolerances phi_d, q_d.
//
// eval() computes all of it for four given hits, so the truth tools measure
// the residuals of true quads on the same arithmetic the finder cuts on.

#include "SeedLayer.h"
#include "SeedFinder.h"  // Quad, SeedCounters

#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"

#include <algorithm>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <array>

namespace mkfit::seeding {

  struct SurfParams {
    float pt_min = 0.9f;   // GeV
    float d0_max = 0.1f;   // cm
    float zv = 25.0f;      // cm, |z0 - z_beamspot| of the a-b line
    float bs_z = 0.0f;     // cm
    float marg_b = 0.002f; // rad
    float phi_c = 0.004f;  // rad
    float q_c = 0.10f;     // cm
    float phi_d = 0.004f;  // rad
    float q_d = 0.10f;     // cm
    int no_cut = 0;        // 1: accept every candidate the fetch returns (for window studies only)
    // d windows growing as 1/pT_est: w = a + b / max(pT_est, pt_min); on when b_phi_d or b_q_d > 0.
    // phi_d and q_d are then the a terms.  The fetch uses pT_est = pt_min, the widest.
    float b_phi_d = 0;  // rad GeV
    float b_q_d = 0;    // cm GeV
    // s_ref > 0: the b terms scale with the candidate's own path length s from hit c
    // (3D, cm) as s / s_ref -- the lever arm of the scattering between c and d.
    float s_ref = 0;
    double lever(double s) const { return s_ref > 0 && s > 0 ? s / s_ref : 1.0; }
    // d windows per |eta| slice of the triplet's helix; a slice replaces phi_d, b_phi_d, q_d, b_q_d
    struct EtaWin {
      float lo, hi, aphi, bphi, aq, bq;
    };
    static constexpr int kMaxEtaWin = 8;
    EtaWin eta_win[kMaxEtaWin];
    int n_eta_win = 0;
    void add_eta_win(const EtaWin &w) {
      if (n_eta_win < kMaxEtaWin)
        eta_win[n_eta_win++] = w;
    }
    // the d windows for a triplet at |eta| ae: this, or a copy with the slice's a and b
    SurfParams for_eta(double ae) const {
      SurfParams Q = *this;
      for (int i = 0; i < n_eta_win; ++i)
        if (ae >= eta_win[i].lo && ae < eta_win[i].hi) {
          Q.phi_d = eta_win[i].aphi, Q.b_phi_d = eta_win[i].bphi, Q.q_d = eta_win[i].aq, Q.b_q_d = eta_win[i].bq;
          break;
        }
      return Q;
    }
    double wphi_d(double pte, double s = -1) const { return phi_d + lever(s) * b_phi_d / std::max((double)pt_min, pte); }
    double wq_d(double pte, double s = -1) const { return q_d + lever(s) * b_q_d / std::max((double)pt_min, pte); }
  };

  // A layer as the surface finder sees it.
  struct SurfLayer {
    using AxPhi = axis_pow2_u1<float, unsigned short, 16, 8>;
    using AxQ = axis<float, unsigned short, 16, 8>;
    using L = SeedLayer<AxPhi, AxQ>;
    int id = -1;
    bool disc = false;
    double qbar_lo = 0, qbar_hi = 0;  // r range (barrel) or z range (disc)
    double q_lo = 0, q_hi = 0;        // z range (barrel) or r range (disc)
    L sl;

    static unsigned int nq(double lo, double hi, double bin) {
      return std::max(1u, (unsigned int)std::ceil((hi - lo) / bin));
    }
    SurfLayer(int id_, const LayerInfo &li, double qbin)
        : id(id_),
          disc(!li.is_barrel()),
          qbar_lo(li.is_barrel() ? li.rin() : li.zmin()),
          qbar_hi(li.is_barrel() ? li.rout() : li.zmax()),
          q_lo(li.is_barrel() ? li.zmin() : li.rin()),
          q_hi(li.is_barrel() ? li.zmax() : li.rout()),
          sl((float)q_lo, (float)q_hi, nq(li.is_barrel() ? li.zmin() : li.rin(), li.is_barrel() ? li.zmax() : li.rout(), qbin),
             !li.is_barrel()) {}

    // r and phi of each hit as P3 computes them from the layer's float x, y, in double
    std::vector<double> pr_, pphi_;
    void fill(const HitVec &hits) {
      sl.fill(hits);
      pr_.resize(sl.n());
      pphi_.resize(sl.n());
      for (unsigned int k = 0; k < sl.n(); ++k) {
        pr_[k] = std::hypot((double)sl.x_[k], (double)sl.y_[k]);
        pphi_[k] = std::atan2((double)sl.y_[k], (double)sl.x_[k]);
      }
    }
    // the hit's own qbar and q
    double qbar(unsigned int k) const { return disc ? sl.z_[k] : sl.r_[k]; }
    double q(unsigned int k) const { return disc ? sl.r_[k] : sl.z_[k]; }
  };

  namespace surf {
    constexpr double kTwoPi = 2 * 3.14159265358979323846;
    inline double wrap(double d) {
      while (d > kTwoPi / 2)
        d -= kTwoPi;
      while (d < -kTwoPi / 2)
        d += kTwoPi;
      return d;
    }

    // A point with its transverse radius and azimuth computed once, in double, at
    // construction: r() and phi() are called many times per hit.
    struct P3 {
      double x = 0, y = 0, z = 0, rr = 0, pp = 0;
      P3() = default;
      P3(double x_, double y_, double z_) : x(x_), y(y_), z(z_), rr(std::hypot(x_, y_)), pp(std::atan2(y_, x_)) {}
      // with r and phi already computed by the same expressions (SurfLayer's cache)
      static P3 cached(double x_, double y_, double z_, double r_, double phi_) {
        P3 p;
        p.x = x_, p.y = y_, p.z = z_, p.rr = r_, p.pp = phi_;
        return p;
      }
      double r() const { return rr; }
      double phi() const { return pp; }
    };

    // The circle through a, b, c in the transverse plane: signed curvature k
    // (> 0 turning left, i.e. counter-clockwise), unit tangent at c, and the
    // arc lengths a->b, b->c; z is taken linear in arc length.
    struct Helix {
      double k = 0, tx = 0, ty = 0, s_ab = 0, s_bc = 0, cot = 0;
      P3 c{};
      bool ok = false;

      static double arc(double chord, double k) {
        const double h = 0.5 * std::abs(k) * chord;
        if (h < 1e-9)
          return chord;
        return h >= 1 ? kTwoPi / 2 / std::abs(k) : 2 * std::asin(h) / std::abs(k);
      }

      Helix(const P3 &a, const P3 &b, const P3 &c_) : c(c_) {
        const double ux = b.x - a.x, uy = b.y - a.y, vx = c.x - b.x, vy = c.y - b.y;
        const double lu = std::hypot(ux, uy), lv = std::hypot(vx, vy), lw = std::hypot(c.x - a.x, c.y - a.y);
        if (lu < 1e-6 || lv < 1e-6)
          return;
        k = 2 * (ux * vy - uy * vx) / (lu * lv * lw);
        // the tangent at c is the chord b->c turned by half its central angle
        const double half = std::asin(std::clamp(0.5 * k * lv, -1.0, 1.0));
        const double cx = vx / lv, cy = vy / lv, cs = std::cos(half), sn = std::sin(half);
        tx = cx * cs - cy * sn;
        ty = cx * sn + cy * cs;
        s_ab = arc(lu, k);
        s_bc = arc(lv, k);
        cot = (c.z - a.z) / (s_ab + s_bc);
        ok = true;
      }

      // the point at arc length ds past c
      void at(double ds, double &x, double &y) const {
        const double al = k * ds;
        double A, B;  // sin(al) / k, (1 - cos(al)) / k
        if (std::abs(al) < 1e-4) {
          A = ds * (1 - al * al / 6);
          B = ds * al / 2 * (1 - al * al / 12);
        } else {
          A = std::sin(al) / k;
          B = (1 - std::cos(al)) / k;
        }
        x = c.x + A * tx - B * ty;
        y = c.y + A * ty + B * tx;
      }
      void dir(double ds, double &dx, double &dy) const {
        const double al = k * ds, cs = std::cos(al), sn = std::sin(al);
        dx = tx * cs - ty * sn;
        dy = tx * sn + ty * cs;
      }

      // ds > 0 with |P(ds)| = r, by Newton from the straight-line solution
      bool cross_r(double r, double &ds) const {
        const double ct = c.x * tx + c.y * ty, c2 = c.x * c.x + c.y * c.y;
        const double disc = ct * ct - c2 + r * r;
        if (disc < 0)
          return false;
        ds = -ct + std::sqrt(disc);
        for (int it = 0; it < 6; ++it) {
          double x, y, dx, dy;
          at(ds, x, y);
          dir(ds, dx, dy);
          const double f = x * x + y * y - r * r, fp = 2 * (x * dx + y * dy);
          if (std::abs(fp) < 1e-12)
            return false;
          const double step = f / fp;
          ds -= step;
          if (std::abs(step) < 1e-9)
            break;
        }
        double x, y;
        at(ds, x, y);
        return ds > 0 && std::abs(std::hypot(x, y) - r) < 1e-5;
      }

      // Prediction on a target surface at qbar u: disc (u = z) gives r and
      // phi, barrel (u = r) gives z and phi; s3, if given, the 3D path length from c.
      bool predict(bool disc, double u, double &q, double &phi, double *s3 = nullptr) const {
        double ds;
        if (disc) {
          if (std::abs(cot) < 1e-9)
            return false;
          ds = (u - c.z) / cot;
          if (ds <= 0)
            return false;
        } else if (!cross_r(u, ds))
          return false;
        double x, y;
        at(ds, x, y);
        phi = std::atan2(y, x);
        q = disc ? std::hypot(x, y) : c.z + cot * ds;
        if (s3)
          *s3 = ds * std::sqrt(1 + cot * cot);
        return true;
      }

      double pt() const { return std::abs(k) > 1e-12 ? 0.0114 / std::abs(k) : 1e9; }
    };

    // The a-b line in the target's (u, q) and (u, phi), and the D0 term.
    struct Line {
      double ua, ub, qa, qb, pa, dp, ia, ib;
      Line(bool disc, const P3 &a, const P3 &b) {
        ua = disc ? a.z : a.r();
        ub = disc ? b.z : b.r();
        qa = disc ? a.r() : a.z;
        qb = disc ? b.r() : b.z;
        pa = a.phi();
        dp = wrap(b.phi() - pa);
        ia = 1 / a.r();
        ib = 1 / b.r();
      }
      // fraction along a->b; > 1 past b
      double t(double u) const { return (u - ua) / (ub - ua); }
      double q(double u) const { return qa + (qb - qa) * t(u); }
      double phi(double u) const { return pa + dp * t(u); }
      // |1/r - linear in u|, what D0 adds to phi beyond the line
      double d0term(double u, double r) const { return std::abs(1 / r - (ia + (ib - ia) * t(u))); }
    };
  }  // namespace surf

  // Everything the cuts look at, for four given hits.  NaN where a stage could
  // not be evaluated.
  struct SurfEval {
    double z0 = NAN, dphi_b = NAN, w_b = NAN;       // b: residual and its geometric window (no margin)
    double dphi_c = NAN, wd0_c = NAN, dq_c = NAN;   // c: residuals and the D0 term (no phi_c)
    double dphi_d = NAN, dq_d = NAN, pte = NAN;     // d
    double s_cd = NAN;                              // 3D path length from c to d on the helix
  };

  inline double surf_wb(const SurfParams &P, double ra, double rb) {
    const double Rmin = P.pt_min / (0.003 * 3.8);
    return (rb - ra) / (2 * Rmin) + P.d0_max * (1 / ra - 1 / rb);
  }

  inline surf::P3 surf_p3(const SurfLayer &L, unsigned int k) {
    return surf::P3::cached(L.sl.x_[k], L.sl.y_[k], L.sl.z_[k], L.pr_[k], L.pphi_[k]);
  }

  inline void surf_eval(const SurfParams &P, const SurfLayer *Ls[4], const surf::P3 h[4], SurfEval &e) {
    using namespace surf;
    const double ra = h[0].r(), rb = h[1].r();
    e.z0 = h[0].z - ra * (h[1].z - h[0].z) / (rb - ra);
    e.dphi_b = wrap(h[1].phi() - h[0].phi());
    e.w_b = surf_wb(P, ra, rb);
    const bool dc = Ls[2]->disc;
    const Line ln(dc, h[0], h[1]);
    const double uc = dc ? h[2].z : h[2].r();
    if (ln.t(uc) > 1) {
      e.dphi_c = wrap(h[2].phi() - ln.phi(uc));
      e.wd0_c = P.d0_max * ln.d0term(uc, h[2].r());
      e.dq_c = (dc ? h[2].r() : h[2].z) - ln.q(uc);
    }
    const Helix hx(h[0], h[1], h[2]);
    if (!hx.ok)
      return;
    e.pte = hx.pt();
    const bool dd = Ls[3]->disc;
    double q, phi;
    if (hx.predict(dd, dd ? h[3].z : h[3].r(), q, phi, &e.s_cd)) {
      e.dphi_d = wrap(h[3].phi() - phi);
      e.dq_d = (dd ? h[3].r() : h[3].z) - q;
    }
  }

  namespace surf {
    inline auto phi_bins(const SurfLayer &S, double c, double w) {
      w = std::min(w, 0.9 * kTwoPi / 2);
      return S.sl.phi_range((float)(c - w), (float)(c + w));
    }
    // q bins over [lo, hi], widened by a float ULP margin
    inline auto q_bins(const SurfLayer &S, double lo, double hi) {
      return S.sl.q_range((float)(lo - 1e-4), (float)(hi + 1e-4));
    }
  }  // namespace surf

  // Stage b: every hit of layer B that makes a doublet with ha -- the geometric
  // phi bound and the a-b line reaching r = 0 within zv of the beam spot.
  // Returns false if B cannot be reached from ha at all.  f(kb, hb).
  // Stage b's fetch: the phi and q bin ranges of layer Bl for the a-hit ha; false if empty.
  struct SurfFetch {
    SurfLayer::AxPhi::I_pair p;
    SurfLayer::AxQ::I_pair q;
  };
  inline bool surf_b_fetch(const SurfParams &P, const surf::P3 &ha, const SurfLayer &Bl, SurfFetch &fe) {
    using namespace surf;
    const double zlo = P.bs_z - P.zv, zhi = P.bs_z + P.zv;
    const double ra = ha.r(), pa = ha.phi();
    // phi from the geometric bound at the layer's largest r; q from the beam-line window
    const double rbmax = Bl.disc ? Bl.q_hi : Bl.qbar_hi;
    const double w_ab = surf_wb(P, ra, rbmax) + P.marg_b;
    double bq_lo, bq_hi;
    if (!Bl.disc) {
      // z_b = z_a + (z_a - z0) (r_b - r_a) / r_a, bilinear: take the corners
      double zz[4];
      int n = 0;
      for (double rb : {Bl.qbar_lo, Bl.qbar_hi})
        for (double z0 : {zlo, zhi})
          zz[n++] = ha.z + (ha.z - z0) * (rb - ra) / ra;
      bq_lo = *std::min_element(zz, zz + 4);
      bq_hi = *std::max_element(zz, zz + 4);
    } else {
      // r_b = r_a + r_a (z_b - z_a) / (z_a - z0); in s z, s = side of the disc
      const double s = (Bl.qbar_lo + Bl.qbar_hi) > 0 ? 1 : -1;
      const double za = s * ha.z, zb0 = std::min(s * Bl.qbar_lo, s * Bl.qbar_hi), zb1 = std::max(s * Bl.qbar_lo, s * Bl.qbar_hi);
      const double z0min = std::min(s * zlo, s * zhi), z0max = std::max(s * zlo, s * zhi);
      if (za - z0min <= 0)
        return false;
      bq_lo = ra + ra * (zb0 - za) / (za - z0min);
      bq_hi = za - z0max > 0 ? ra + ra * (zb1 - za) / (za - z0max) : Bl.q_hi;
    }
    if (bq_hi < Bl.q_lo || bq_lo > Bl.q_hi)
      return false;
    fe.p = phi_bins(Bl, pa, w_ab);
    fe.q = q_bins(Bl, bq_lo, bq_hi);
    return true;
  }

  template <typename F>
  inline void surf_stage_b(const SurfParams &P, const surf::P3 &ha, const SurfLayer &Bl, SeedCounters &cnt, F &&f) {
    using namespace surf;
    SurfFetch fe;
    if (!surf_b_fetch(P, ha, Bl, fe))
      return;
    const double zlo = P.bs_z - P.zv, zhi = P.bs_z + P.zv;
    const double ra = ha.r(), pa = ha.phi();
    Bl.sl.for_each_in(fe.p, fe.q, [&](unsigned int kb) {
      const P3 hb = surf_p3(Bl, kb);
      const double rb = hb.r();
      if (rb <= ra + 0.1)
        return;
      if (!P.no_cut && std::abs(wrap(hb.phi() - pa)) > surf_wb(P, ra, rb) + P.marg_b)
        return;
      const double z0 = ha.z - ra * (hb.z - ha.z) / (rb - ra);
      if (!P.no_cut && (z0 < zlo || z0 > zhi))
        return;
      cnt.doublets++;
      f(kb, hb);
    });
  }

  // Stage b in float on the layer's arrays (K2a): the same fetch as surf_stage_b, and
  // per candidate the same three cuts with no division -- the phi bound from the
  // cached 1/r, z0 in [zlo, zhi] multiplied through by dr = r_b - r_a > 0. The
  // survivors are handed on as surf_stage_b does (double P3 from the cache).
  template <typename F>
  inline void surf_stage_b_fast(const SurfParams &P, const surf::P3 &ha, const SurfLayer &Bl, SeedCounters &cnt, F &&f) {
    using namespace surf;
    SurfFetch fe;
    if (!surf_b_fetch(P, ha, Bl, fe))
      return;
    constexpr float kPi = 3.14159265358979f, k2Pi = 6.28318530717959f;
    const float ra = ha.r(), pa = ha.phi(), za = ha.z, inva = 1.0f / ra;
    const float inv2R = 0.003f * 3.8f / (2.0f * P.pt_min), d0 = P.d0_max, marg = P.marg_b;
    const float zlo = P.bs_z - P.zv, zhi = P.bs_z + P.zv;
    const float *phi = Bl.sl.phi_.data(), *r = Bl.sl.r_.data(), *z = Bl.sl.z_.data(), *ir = Bl.sl.invr_.data();
    Bl.sl.for_each_run(fe.p, fe.q, [&](unsigned int b, unsigned int e) {
      for (unsigned int i0 = b; i0 < e; i0 += 64) {
        const unsigned int n = std::min(64u, e - i0);
        unsigned char m[64];
        for (unsigned int j = 0; j < n; ++j) {
          const unsigned int i = i0 + j;
          const float dr = r[i] - ra;
          float dp = phi[i] - pa;
          dp = dp > kPi ? dp - k2Pi : dp;
          dp = dp < -kPi ? dp + k2Pi : dp;
          const float w = dr * inv2R + d0 * (inva - ir[i]) + marg;
          const float num = za * dr - ra * (z[i] - za);
          m[j] = (dr > 0.1f) & (std::abs(dp) <= w) & (num >= zlo * dr) & (num <= zhi * dr);
        }
        for (unsigned int j = 0; j < n; ++j)
          if (m[j]) {
            cnt.doublets++;
            f(i0 + j, surf_p3(Bl, i0 + j));
          }
      }
    });
  }

  // Stage c: every hit of layer C on the line through h1, h2 in the target's
  // (qbar, q) and (qbar, phi); phi_c, q_c from P.  f(kc, hc).
  template <typename F>
  inline void surf_stage_c(const SurfParams &P, const surf::P3 &h1, const surf::P3 &h2, const SurfLayer &C, SeedCounters &cnt, F &&f) {
    using namespace surf;
    const Line ln(C.disc, h1, h2);
    const double u0 = C.qbar_lo, u1 = C.qbar_hi;
    if (std::max(ln.t(u0), ln.t(u1)) <= 1)
      return;
    const double cq0 = ln.q(u0), cq1 = ln.q(u1), cp0 = ln.phi(u0), cp1 = ln.phi(u1);
    // the D0 term is largest at the far edge; r there from the line (barrel: the edge itself)
    const double ufar = std::abs(ln.t(u0)) > std::abs(ln.t(u1)) ? u0 : u1;
    const double rfar = std::max(0.5, C.disc ? ln.q(ufar) : ufar);
    const double wc_phi = P.d0_max * ln.d0term(ufar, rfar) + P.phi_c;
    const auto qc = q_bins(C, std::min(cq0, cq1) - P.q_c, std::max(cq0, cq1) + P.q_c);
    const double cmid = cp0 + 0.5 * wrap(cp1 - cp0), chalf = 0.5 * std::abs(wrap(cp1 - cp0));
    C.sl.for_each_in(phi_bins(C, cmid, chalf + wc_phi), qc, [&](unsigned int kc) {
      cnt.c_touched++;
      const P3 hc = surf_p3(C, kc);
      const double uc = C.qbar(kc);
      if (ln.t(uc) <= 1)
        return;
      if (!P.no_cut) {
        if (std::abs((C.disc ? hc.r() : hc.z) - ln.q(uc)) > P.q_c)
          return;
        if (std::abs(wrap(hc.phi() - ln.phi(uc))) > P.d0_max * ln.d0term(uc, hc.r()) + P.phi_c)
          return;
      }
      cnt.triplets++;
      f(kc, hc);
    });
  }

  // Stage d: every hit of layer D on the helix hx; windows P.wphi_d, P.wq_d.  f(kd, hd).
  template <typename F>
  inline void surf_stage_d(const SurfParams &P0, const surf::Helix &hx, const SurfLayer &D, SeedCounters &cnt, F &&f) {
    using namespace surf;
    if (!hx.ok)
      return;
    SurfParams Pe;
    const SurfParams *pp = &P0;
    if (P0.n_eta_win > 0) {
      Pe = P0.for_eta(std::asinh(std::abs(hx.cot)));
      pp = &Pe;
    }
    const SurfParams &P = *pp;
    double dq0, dp0, dq1, dp1, s0 = -1, s1 = -1;
    const bool ok0 = hx.predict(D.disc, D.qbar_lo, dq0, dp0, &s0);
    const bool ok1 = hx.predict(D.disc, D.qbar_hi, dq1, dp1, &s1);
    if (!ok0 && !ok1)
      return;
    if (!ok0)
      dq0 = dq1, dp0 = dp1;
    if (!ok1)
      dq1 = dq0, dp1 = dp0;
    // fetch with the window at the longer end of the layer's slab; cut per hit on its own path length
    const double pte = hx.pt(), smax = std::max(s0, s1), wpd = P.wphi_d(pte, smax), wqd = P.wq_d(pte, smax);
    const auto qd = q_bins(D, std::min(dq0, dq1) - wqd, std::max(dq0, dq1) + wqd);
    const double dmid = dp0 + 0.5 * wrap(dp1 - dp0), dhalf = 0.5 * std::abs(wrap(dp1 - dp0));
    D.sl.for_each_in(phi_bins(D, dmid, dhalf + wpd), qd, [&](unsigned int kd) {
      cnt.d_touched++;
      const P3 hd = surf_p3(D, kd);
      double q, phi, s = -1;
      if (!hx.predict(D.disc, D.qbar(kd), q, phi, &s))
        return;
      const double wq = P.s_ref > 0 ? P.wq_d(pte, s) : wqd, wp = P.s_ref > 0 ? P.wphi_d(pte, s) : wpd;
      if (!P.no_cut && (std::abs((D.disc ? hd.r() : hd.z) - q) > wq || std::abs(wrap(hd.phi() - phi)) > wp))
        return;
      cnt.quads++;
      f(kd, hd);
    });
  }

  // The finder for one pattern.  L[0..3]: the pattern's layers, filled.
  inline void find_quads_surf(const SurfParams &P, const SurfLayer *L[4], std::vector<Quad> &out, SeedCounters &cnt) {
    const SurfLayer &A = *L[0], &Bl = *L[1], &C = *L[2], &D = *L[3];
    for (unsigned int ka = 0; ka < A.sl.n(); ++ka) {
      const surf::P3 ha = surf_p3(A, ka);
      surf_stage_b(P, ha, Bl, cnt, [&](unsigned int kb, const surf::P3 &hb) {
        surf_stage_c(P, ha, hb, C, cnt, [&](unsigned int kc, const surf::P3 &hc) {
          const surf::Helix hx(ha, hb, hc);
          surf_stage_d(P, hx, D, cnt, [&](unsigned int kd, const surf::P3 &) {
            out.push_back({A.sl.orig_[ka], Bl.sl.orig_[kb], C.sl.orig_[kc], D.sl.orig_[kd]});
          });
        });
      });
    }
  }

  // Phase-space ownership: which pattern a quad belongs to, from its own r-z
  // line.  The line through the a- and d-hits (z0, cot) is intersected with
  // every pixel layer's envelope; a layer is DEFINITE if the crossing is inside
  // it by more than delta (cm: in z for a barrel layer, in r for a disc),
  // MAYBE within delta of an edge.  A pattern owns the quad if its four layers
  // can be the first four crossed: every pattern layer is definite or maybe,
  // and no definite layer crossed before the pattern's last one is missing
  // from it.  delta = 0 gives one owner per quad (up to exact edges); a larger
  // delta lets neighbouring patterns overlap by that much; delta < 0 turns the
  // test off.
  struct SurfOwnership {
    struct Env {
      int id;
      bool disc;
      double pos, pos_lo, pos_hi;  // qbar: barrel r (mid, rin, rout); disc z (mid, zmin, zmax)
      double lo, hi;               // q range: barrel z; disc r
    };
    std::vector<Env> env;
    double delta = -1;
    int max_skip = 0;  // definite layers crossed before the pattern's last one and not in it

    // the pixel layers, plus the non-pixel layers patterns use (e.g. OT1-P)
    template <typename C>
    void setup(const TrackerInfo &ti, const C &extra) {
      env.clear();
      for (int l = 0; l < ti.n_layers(); ++l) {
        const bool pix = l <= 3 || (l >= 16 && l <= 27) || (l >= 38 && l <= 49);
        if (!pix && !extra.count(l))
          continue;
        const LayerInfo &li = ti[l];
        if (li.is_barrel())
          env.push_back({l, false, 0.5 * (li.rin() + li.rout()), li.rin(), li.rout(), li.zmin(), li.zmax()});
        else
          env.push_back({l, true, 0.5 * (li.zmin() + li.zmax()), li.zmin(), li.zmax(), li.rin(), li.rout()});
      }
    }

    // 0 no, 1 maybe, 2 definite; s: the path parameter (r) at the layer's mid qbar.
    // The layer is a slab in qbar (a barrel layer's two shells span ~1 cm in r):
    // the line's q over the whole slab is compared with the q range, so a track
    // that meets the inner shell inside the layer and the outer one outside it
    // is a MAYBE, not a no.
    int cross(const Env &e, double z0, double cot, double &s) const {
      double x1, x2;
      if (!e.disc) {
        s = e.pos;
        x1 = z0 + cot * e.pos_lo;
        x2 = z0 + cot * e.pos_hi;
      } else {
        if (std::abs(cot) < 1e-9)
          return 0;
        s = (e.pos - z0) / cot;
        if (s <= 0)
          return 0;
        x1 = (e.pos_lo - z0) / cot;
        x2 = (e.pos_hi - z0) / cot;
      }
      const double xl = std::min(x1, x2), xh = std::max(x1, x2);
      if (xl > e.lo + delta && xh < e.hi - delta)
        return 2;
      if (xh > e.lo - delta && xl < e.hi + delta)
        return 1;
      return 0;
    }

    // one line per layer the line meets: id, state (1 maybe, 2 definite), s, coordinate
    void print(double z0, double cot) const {
      for (const Env &e : env) {
        double s;
        const int st = cross(e, z0, cot, s);
        if (st)
          printf(" %d:%d@%.1f", e.id, st, s);
      }
      printf("\n");
    }

    // max_skip_ot: the allowance for a pattern with a non-pixel layer (< 0: max_skip)
    int max_skip_ot = -1;
    static bool is_pix(int l) { return l <= 3 || (l >= 16 && l <= 27) || (l >= 38 && l <= 49); }

    // Definite crossings of layers NOT in the pattern: before the pattern's first
    // layer (lead) and between its first and last (inner).  The chain with a late
    // start takes a quad with lead <= start_holes and inner == 0.
    void skips(const int lay[4], double z0, double cot, int &lead, int &inner) const {
      double s_first = 1e30, s_last = -1;
      std::vector<double> other;
      for (const Env &e : env) {
        double s;
        const int st = cross(e, z0, cot, s);
        const bool inpat = e.id == lay[0] || e.id == lay[1] || e.id == lay[2] || e.id == lay[3];
        if (inpat)
          s_first = std::min(s_first, s), s_last = std::max(s_last, s);
        else if (st == 2)
          other.push_back(s);
      }
      lead = inner = 0;
      for (double s : other) {
        lead += s < s_first;
        inner += s > s_first && s < s_last;
      }
    }

    bool owns(const int lay[4], double z0, double cot) const {
      if (delta < 0)
        return true;
      const bool ot = !is_pix(lay[0]) || !is_pix(lay[1]) || !is_pix(lay[2]) || !is_pix(lay[3]);
      const int allow = ot && max_skip_ot >= 0 ? max_skip_ot : max_skip;
      std::vector<std::pair<double, int>> v;  // (s, state) per layer, pattern membership in the sign
      double s_last = -1;
      int n_in = 0;
      for (const Env &e : env) {
        double s;
        const int st = cross(e, z0, cot, s);
        const bool inpat = e.id == lay[0] || e.id == lay[1] || e.id == lay[2] || e.id == lay[3];
        if (inpat) {
          if (st == 0)
            return false;
          ++n_in;
          s_last = std::max(s_last, s);
        } else if (st == 2)
          v.push_back({s, 0});
      }
      if (n_in != 4)
        return false;
      int skipped = 0;
      for (const auto &p : v)
        skipped += p.first < s_last;
      return skipped <= allow;
    }
  };

  // The feed-forward chain: one pass over the layers of one z side in crossing
  // order (barrel B1..B4, the discs by |z|, then OT1-P), instead of a list of
  // four-layer patterns.
  //
  //   start   doublets on layer pairs that a line from the beam region can cross
  //           first and second, or with up to max_holes crossed layers skipped
  //           (derived at setup by scanning lines, not listed by hand)
  //   extend  a candidate's NEXT layer is the next one its own line (first to
  //           last hit) crosses; it is queued there.  Layers are processed in
  //           order, a candidate only moves to later layers, so the pass is
  //           strictly feed-forward.
  //   miss    no compatible hit in the target: the candidate moves on to the
  //           following crossed layer, charged a hole if the target was a
  //           DEFINITE crossing (a MAYBE costs nothing), dropped past max_holes.
  //
  // Windows come from the per-pattern tables (win_c by the three layers, win_d
  // by the four), defaults from P where a combination has none.
  struct SurfChain {
    struct DW {
      float aphi, bphi, aq, bq, sref = 0;
      std::vector<SurfParams::EtaWin> eta;
    };
    std::map<std::array<int, 3>, std::pair<float, float>> win_c;
    std::map<std::array<int, 4>, DW> win_d;
    SurfParams P;
    int max_holes = 1;
    int hole_always = 0;  // 1: also forward a candidate that found hits (a hole in competition)
    int start_holes = -1; // holes allowed in the start doublet (< 0: max_holes); 0 = start on the first two crossed
    int lead_only = 0;    // 1: the start's holes may only LEAD (a later start), a and b consecutive crossings
    int known_only = 1;   // search a target only for a layer combination with a window table; else move on
    int fast = 0;         // 1: the float kernels (K2) where they exist; 0: the double reference
    int max_holes_ot = 0; // holes a candidate may carry into an outer-tracker layer: 0 = OT1-P only
                          // where the pixels have a geometric gap, not after a missed hit
    int side = 1;
    const SurfOwnership *own = nullptr;
    std::vector<int> order;                  // mkFit layer ids, crossing order
    std::vector<int> env_of;                 // per chain position: index into own->env
    std::vector<std::pair<int, int>> starts; // chain positions (a, b)
    long n_forwarded = 0, n_dropped = 0;

    static int mirror(int l) { return (l >= 16 && l <= 27) ? l + 22 : l; }

    // order: B1..B4, the side's discs, OT1-P if present
    void setup(const SurfOwnership &o, int side_, const std::set<int> &have) {
      own = &o;
      side = side_;
      order.clear();
      for (int l : {0, 1, 2, 3})
        order.push_back(l);
      for (int l = 16; l <= 27; ++l)
        order.push_back(side > 0 ? l : l + 22);
      if (have.count(4))
        order.push_back(4);
      env_of.assign(order.size(), -1);
      for (int p = 0; p < (int)order.size(); ++p)
        for (int e = 0; e < (int)o.env.size(); ++e)
          if (o.env[e].id == order[p])
            env_of[p] = e;
      // start pairs: scan lines z0 in the beam region, |eta| 0-4.2 on this side
      std::set<std::pair<int, int>> sp;
      for (double z0 = P.bs_z - P.zv; z0 <= P.bs_z + P.zv + 1e-9; z0 += 0.5)
        for (double eta = 0; eta < 4.2; eta += 0.002) {
          const double cot = side * std::sinh(eta);
          std::vector<int> seq;  // pixel chain positions crossed (maybe or definite)
          for (int p = 0; p < (int)order.size(); ++p) {
            if (!SurfOwnership::is_pix(order[p]) || env_of[p] < 0)
              continue;
            double sdum;
            if (o.cross(o.env[env_of[p]], z0, cot, sdum))
              seq.push_back(p);
          }
          const int hs = start_holes < 0 ? max_holes : start_holes;
          for (int i = 0; i <= hs && i < (int)seq.size(); ++i)
            for (int j = i + 1; j <= (lead_only ? i + 1 : i + 1 + hs - i) && j < (int)seq.size(); ++j)
              sp.insert({seq[i], seq[j]});
        }
      starts.assign(sp.begin(), sp.end());
    }

    // crossing state of chain position p for the line through h1, h2
    // the r-z line through two hits, as the crossing test takes it
    struct LineRZ {
      double z0 = 0, cot = 0;
      LineRZ() = default;
      LineRZ(const surf::P3 &h1, const surf::P3 &h2) {
        const double r1 = h1.r(), r2 = h2.r();
        cot = (h2.z - h1.z) / (r2 - r1);
        z0 = h1.z - cot * r1;
      }
    };
    int state(int p, const LineRZ &l) const {
      if (env_of[p] < 0)
        return 0;
      double sdum;
      return own->cross(own->env[env_of[p]], l.z0, l.cot, sdum);
    }

    // A candidate holds its hits as (chain position, index in the layer); the points
    // come back from the layer's cache (surf_p3), identical to storing them.
    struct Cand {
      int nh = 0, holes = 0;
      int pos[4];
      unsigned int k[4];
      LineRZ ln;  // through its first and last hit
    };

    // Window parameters per chain-position combination, resolved from win_c / win_d
    // once (at the first run, after the tables are filled): an index into par_,
    // -1 where the combination has no table; par_[0] is P.
    std::vector<SurfParams> par_;
    std::vector<int> idx_c_, idx_d_;  // (i, j, p) and (i, j, k, p), n positions each
    bool par_built_ = false;
    void build_params() {
      const int n = order.size();
      par_.assign(1, P);
      idx_c_.assign(n * n * n, -1);
      idx_d_.assign(n * n * n * n, -1);
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          for (int p = 0; p < n; ++p) {
            auto it = win_c.find({order[i], order[j], order[p]});
            if (it != win_c.end()) {
              SurfParams Q = P;
              Q.phi_c = it->second.first, Q.q_c = it->second.second;
              par_.push_back(Q);
              idx_c_[(i * n + j) * n + p] = par_.size() - 1;
            }
            for (int k = 0; k < n; ++k) {
              auto jt = win_d.find({order[i], order[j], order[k], order[p]});
              if (jt == win_d.end())
                continue;
              SurfParams Q = P;
              Q.phi_d = jt->second.aphi, Q.b_phi_d = jt->second.bphi, Q.q_d = jt->second.aq, Q.b_q_d = jt->second.bq,
              Q.s_ref = jt->second.sref;
              for (const auto &w : jt->second.eta)
                Q.add_eta_win(w);
              par_.push_back(Q);
              idx_d_[((i * n + j) * n + k) * n + p] = par_.size() - 1;
            }
          }
      par_built_ = true;
    }
    int idx_c(const Cand &c, int p) const {
      const int n = order.size();
      return idx_c_[(c.pos[0] * n + c.pos[1]) * n + p];
    }
    int idx_d(const Cand &c, int p) const {
      const int n = order.size();
      return idx_d_[((c.pos[0] * n + c.pos[1]) * n + c.pos[2]) * n + p];
    }
    // the per-target-layer queues, kept across events for their capacity
    std::vector<std::vector<Cand>> Qc_;
    std::vector<std::vector<int>> Qst_;

    // the next chain position after p its line crosses, and its state; -1 if none
    int next(const Cand &c, int p, int &st) const {
      for (int q = p + 1; q < (int)order.size(); ++q) {
        // OT1-P only as a triplet's 4th hit
        if (!SurfOwnership::is_pix(order[q]) && (c.nh != 3 || c.holes > max_holes_ot))
          continue;
        st = state(q, c.ln);
        if (st)
          return q;
      }
      return -1;
    }

    // L: mkFit layer id -> filled layer.  out: (layer ids, original hit indices)
    void run(const std::map<int, const SurfLayer *> &L, std::vector<std::pair<std::array<int, 4>, Quad>> &out,
             SeedCounters &cnt) {
      const int n = order.size();
      if (!par_built_)
        build_params();
      Qc_.resize(n);
      Qst_.resize(n);
      auto &Q = Qc_;
      auto &Qst = Qst_;  // the target's crossing state, per queued candidate
      std::vector<const SurfLayer *> lay(n);
      for (int p = 0; p < n; ++p) {
        auto it = L.find(order[p]);
        lay[p] = it == L.end() ? nullptr : it->second;
      }
      auto layer = [&](int p) -> const SurfLayer * { return lay[p]; };
      auto hit = [&](const Cand &c, int i) { return surf_p3(*lay[c.pos[i]], c.k[i]); };
      auto push_next = [&](Cand &c, int p) {
        int st = 0;
        const int q = next(c, p, st);
        if (q < 0 || !layer(q)) {
          ++n_dropped;
          return;
        }
        Q[q].push_back(c);
        Qst[q].push_back(st);
      };
      // start doublets
      for (const auto &se : starts) {
        const SurfLayer *A = layer(se.first), *B = layer(se.second);
        if (!A || !B)
          continue;
        for (unsigned int ka = 0; ka < A->sl.n(); ++ka) {
          const surf::P3 ha = surf_p3(*A, ka);
          auto on_b = [&](unsigned int kb, const surf::P3 &hb) {
            // the side: a doublet belongs to the side its line goes to
            if ((hb.z - ha.z) * side < 0 || (hb.z == ha.z && side < 0))
              return;
            // holes: definite crossings before b that the doublet does not use
            const LineRZ lab(ha, hb);
            int holes = 0, between = 0;
            for (int q = 0; q < se.second && holes <= 8; ++q)
              if (q != se.first && SurfOwnership::is_pix(order[q]) && state(q, lab) == 2) {
                ++holes;
                between += q > se.first;
              }
            if (lead_only && between)
              return;
            if (holes > (start_holes < 0 ? max_holes : start_holes) || !state(se.first, lab) || !state(se.second, lab))
              return;
            Cand c;
            c.nh = 2;
            c.holes = holes;
            c.pos[0] = se.first, c.pos[1] = se.second;
            c.k[0] = ka, c.k[1] = kb;
            c.ln = lab;
            push_next(c, se.second);
          };
          if (fast)
            surf_stage_b_fast(P, ha, *B, cnt, on_b);
          else
            surf_stage_b(P, ha, *B, cnt, on_b);
        }
      }
      // the forward pass
      for (int p = 0; p < n; ++p) {
        const SurfLayer *T = layer(p);
        if (!T)
          continue;
        // Q[p] is not appended to while it is walked: next() only queues at positions after p
        for (size_t i = 0; i < Q[p].size(); ++i) {
          const Cand &c = Q[p][i];
          const int st = Qst[p][i];
          bool found = false;
          const int ix = c.nh == 2 ? idx_c(c, p) : idx_d(c, p);
          const bool known = !known_only || ix >= 0;
          if (!known) {
            // no windows for this combination: not searched, the candidate moves on as after a miss
          } else if (c.nh == 2) {
            const SurfParams &Pc = par_[ix >= 0 ? ix : 0];
            const surf::P3 h0 = hit(c, 0), h1 = hit(c, 1);
            surf_stage_c(Pc, h0, h1, *T, cnt, [&](unsigned int kc, const surf::P3 &hc) {
              found = true;
              Cand t = c;
              t.nh = 3, t.pos[2] = p, t.k[2] = kc;
              t.ln = LineRZ(h0, hc);
              push_next(t, p);
            });
          } else if (c.nh == 3) {
            const SurfParams &Pd = par_[ix >= 0 ? ix : 0];
            const surf::Helix hx(hit(c, 0), hit(c, 1), hit(c, 2));
            surf_stage_d(Pd, hx, *T, cnt, [&](unsigned int kd, const surf::P3 &) {
              found = true;
              std::array<int, 4> ids{order[c.pos[0]], order[c.pos[1]], order[c.pos[2]], order[p]};
              const SurfLayer *La = layer(c.pos[0]), *Lb = layer(c.pos[1]), *Lc = layer(c.pos[2]);
              out.push_back({ids, {La->sl.orig_[c.k[0]], Lb->sl.orig_[c.k[1]], Lc->sl.orig_[c.k[2]], T->sl.orig_[kd]}});
            });
          }
          if (!found || hole_always) {
            if (c.holes + (st == 2) <= max_holes) {
              Cand f = c;
              f.holes += st == 2;
              ++n_forwarded;
              push_next(f, p);
            } else
              ++n_dropped;
          }
        }
        Q[p].clear();
        Qst[p].clear();
      }
    }
  };

}  // namespace mkfit::seeding

#endif
