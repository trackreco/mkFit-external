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

#include "RecoTracker/MkFitCore/interface/SeedStructures.h"  // SeedLayerOfHits, SeedQuad, SeedCounters
#include "RecoTracker/MkFitCore/interface/SeedChain.h"       // SeedingParams, SeedLayerEnvelopes, SeedChain

#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <array>

namespace mkfit::seeding {

  // the parameters, the crossing envelopes and the chain's configuration are in MkFitCore since 2026-10-01
  using SurfParams = SeedingParams;

  constexpr float kPi = 3.14159265358979323846f;

  // A layer as the surface finder sees it: the seeder's layer of hits, in MkFitCore since 2026-10-01.
  using SurfLayer = SeedLayerOfHits;
  using Quad = SeedQuad;

  namespace surf {
    using seedchain::kTwoPi;
    inline double wrap(double d) { return seedchain::wrap_d(d); }
    using seedchain::P3;

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

      // The first crossing ds > 0 with |P(ds)| = r in closed form: the helix circle
      // (centre C = c + n / k, n the left normal of the tangent) against |P| = r,
      // by the radical line, then the arc from c to the point. d^2 - rho^2 is taken
      // as |c|^2 + 2 (c.n) / k, which does not cancel for a stiff track. No iteration
      // and no sin/cos; for |k| below 1e-7 the straight-line solution.
      bool cross_r_cf(double r, double &ds, double &px, double &py) const {
        const double nx = -ty, ny = tx;
        if (std::abs(k) < 1e-7) {
          const double ct = c.x * tx + c.y * ty, c2 = c.x * c.x + c.y * c.y, disc = ct * ct - c2 + r * r;
          if (disc < 0)
            return false;
          ds = -ct + std::sqrt(disc);
          px = c.x + ds * tx, py = c.y + ds * ty;
          return ds > 0;
        }
        const double ik = 1 / k, Cx = c.x + nx * ik, Cy = c.y + ny * ik;
        const double d2 = Cx * Cx + Cy * Cy, d = std::sqrt(d2);
        if (d < 1e-12)
          return false;
        const double dmr = c.x * c.x + c.y * c.y + 2 * (c.x * nx + c.y * ny) * ik;  // d^2 - rho^2
        const double a = (r * r + dmr) / (2 * d), h2 = r * r - a * a;
        if (h2 < 0)
          return false;
        const double h = std::sqrt(h2), ux = Cx / d, uy = Cy / d;
        // the two points, and the arc from c to each in the direction of motion
        const double qx[2] = {a * ux - h * uy, a * ux + h * uy}, qy[2] = {a * uy + h * ux, a * uy - h * ux};
        const double rcx = c.x - Cx, rcy = c.y - Cy;
        double best = -1;
        for (int j = 0; j < 2; ++j) {
          const double rpx = qx[j] - Cx, rpy = qy[j] - Cy;
          double th = std::atan2(rcx * rpy - rcy * rpx, rcx * rpx + rcy * rpy);  // CCW angle c -> P about C
          double s = th * ik;
          if (s <= 0)
            s += kTwoPi * std::abs(ik);
          if (best < 0 || s < best) {
            best = s;
            px = qx[j], py = qy[j];
          }
        }
        ds = best;
        return ds > 0;
      }
      // predict() with the closed-form barrel crossing; the disc branch is predict()'s.
      bool predict_cf(bool disc, double u, double &q, double &phi, double *s3 = nullptr) const {
        if (disc)
          return predict(true, u, q, phi, s3);
        double ds, x, y;
        if (!cross_r_cf(u, ds, x, y))
          return false;
        phi = std::atan2(y, x);
        q = c.z + cot * ds;
        if (s3)
          *s3 = ds * std::sqrt(1 + cot * cot);
        return true;
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

  inline double surf_wb(const SurfParams &P, double ra, double rb) { return seedchain::b_window(P, ra, rb); }

  inline surf::P3 surf_p3(const SurfLayer &L, unsigned int k) { return seedchain::p3(L, k); }

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
    inline auto phi_bins(const SurfLayer &S, double c, double w) { return seedchain::phi_bins_d(S, c, w); }
    inline auto q_bins(const SurfLayer &S, double lo, double hi) { return seedchain::q_bins_d(S, lo, hi); }
  }  // namespace surf

  // Stage b: every hit of layer B that makes a doublet with ha -- the geometric
  // phi bound and the a-b line reaching r = 0 within zv of the beam spot.
  // Returns false if B cannot be reached from ha at all.  f(kb, hb).
  // Stage b's fetch: the phi and q bin ranges of layer Bl for the a-hit ha; false if empty.
  using SurfFetch = seedchain::BFetch;
  inline bool surf_b_fetch(const SurfParams &P, const surf::P3 &ha, const SurfLayer &Bl, SurfFetch &fe) {
    return seedchain::b_fetch(P, ha, Bl, fe);
  }

  template <typename F>
  inline void surf_stage_b(const SurfParams &P, const surf::P3 &ha, const SurfLayer &Bl, SeedCounters &cnt, F &&f) {
    using namespace surf;
    SurfFetch fe;
    if (!surf_b_fetch(P, ha, Bl, fe))
      return;
    const double zlo = P.bs_z - P.zv, zhi = P.bs_z + P.zv;
    const double ra = ha.r(), pa = ha.phi();
    Bl.for_each_in(fe.p, fe.q, [&](unsigned int kb) {
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
    const float *phi = Bl.phi_.data(), *r = Bl.r_.data(), *z = Bl.z_.data(), *ir = Bl.invr_.data();
    Bl.for_each_run(fe.p, fe.q, [&](unsigned int b, unsigned int e) {
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
    C.for_each_in(phi_bins(C, cmid, chalf + wc_phi), qc, [&](unsigned int kc) {
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
    D.for_each_in(phi_bins(D, dmid, dhalf + wpd), qd, [&](unsigned int kd) {
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

  // --chain-fast-check: the closed-form nodes against Newton, and the quadratic
  // against the exact prediction at every fetched candidate (relative to its window).
  struct SurfFastCheck {
    bool on = false;
    long nodes = 0, cands = 0, node_fail_mismatch = 0;
    double max_node_dq = 0, max_node_dphi = 0, max_rel_q = 0, max_rel_phi = 0;
    // the batched float finder (SeedSurfBatch.h): its stage d prediction at every fetched hit
    long b_cands = 0, b_fail_mismatch = 0;
    double b_max_dq = 0, b_max_dphi = 0, b_max_rel_q = 0, b_max_rel_phi = 0;
  };
  inline SurfFastCheck g_surf_fast_check;

  // Stage d as a pre-filter plus confirm (K2c). The helix is predicted at the
  // layer's two qbar edges and its middle with the closed-form crossing
  // (predict_cf, within 1e-8 cm of Newton), which also give the fetch. q, phi and
  // the path length are quadratics in qbar through the three points; every fetched
  // candidate is tested against them in float with the window widened by
  // kSlackD. The worst quadratic error measured (--chain-fast-check, 5 events,
  // 4.3 M candidates) is 2.0 % of the window in q and 0.42 % in phi, so the widened
  // test keeps every hit the exact cut accepts. The survivors are confirmed with
  // surf_stage_d's own prediction and double cut. Falls back to surf_stage_d if a
  // node prediction fails.
  constexpr float kSlackD = 1.05f;
  template <typename F>
  inline void surf_stage_d_fast(const SurfParams &P0, const surf::Helix &hx, const SurfLayer &D, SeedCounters &cnt, F &&f) {
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
    const double u0 = D.qbar_lo, u1 = D.qbar_hi, um = 0.5 * (u0 + u1), hh = 0.5 * (u1 - u0);
    double q0, p0, q1, p1, qm, pm, s0 = -1, s1 = -1, sm = -1;
    const bool ok0 = hx.predict_cf(D.disc, u0, q0, p0, &s0);
    const bool ok1 = hx.predict_cf(D.disc, u1, q1, p1, &s1);
    const bool okm = hx.predict_cf(D.disc, um, qm, pm, &sm);
    if (!(ok0 && ok1 && okm) || hh <= 0) {
      surf_stage_d(P0, hx, D, cnt, f);
      return;
    }
    SurfFastCheck &CK = g_surf_fast_check;
    if (CK.on) {
      const double uu[3] = {u0, u1, um}, qq[3] = {q0, q1, qm}, ph[3] = {p0, p1, pm};
      for (int j = 0; j < 3; ++j) {
        double qn, pn;
        if (!hx.predict(D.disc, uu[j], qn, pn)) {
          ++CK.node_fail_mismatch;
          continue;
        }
        ++CK.nodes;
        CK.max_node_dq = std::max(CK.max_node_dq, std::abs(qn - qq[j]));
        CK.max_node_dphi = std::max(CK.max_node_dphi, std::abs(wrap(pn - ph[j])));
      }
    }
    // the same fetch as surf_stage_d
    const double pte = hx.pt(), smax = std::max(s0, s1), wpd = P.wphi_d(pte, smax), wqd = P.wq_d(pte, smax);
    const auto qd = q_bins(D, std::min(q0, q1) - wqd, std::max(q0, q1) + wqd);
    const double dmid = p0 + 0.5 * wrap(p1 - p0), dhalf = 0.5 * std::abs(wrap(p1 - p0));
    // quadratics in x = u - um: v(x) = vm + b1 x + b2 x^2
    const double d0p = wrap(p0 - pm), d1p = wrap(p1 - pm);
    const float bq1 = (q1 - q0) / (2 * hh), bq2 = (q1 + q0 - 2 * qm) / (2 * hh * hh);
    const float bp1 = (d1p - d0p) / (2 * hh), bp2 = (d1p + d0p) / (2 * hh * hh);
    const float bs1 = (s1 - s0) / (2 * hh), bs2 = (s1 + s0 - 2 * sm) / (2 * hh * hh);
    const float fum = um, fqm = qm, fpm = pm, fsm = sm;
    // the windows: a + lever(s) b / pT, in float
    const float ipt = 1.0f / (float)std::max((double)P.pt_min, pte);
    const float aq = P.q_d, bq = P.b_q_d * ipt, ap = P.phi_d, bp = P.b_phi_d * ipt, isr = P.s_ref > 0 ? 1.0f / P.s_ref : 0.0f;
    constexpr float kPi = 3.14159265358979f, k2Pi = 6.28318530717959f;
    const float *hphi = D.phi_.data(), *hq = D.disc ? D.r_.data() : D.z_.data(),
                *hu = D.disc ? D.z_.data() : D.r_.data();
    D.for_each_run(phi_bins(D, dmid, dhalf + wpd), qd, [&](unsigned int b, unsigned int e) {
      for (unsigned int i0 = b; i0 < e; i0 += 64) {
        const unsigned int n = std::min(64u, e - i0);
        unsigned char m[64];
        for (unsigned int j = 0; j < n; ++j) {
          const unsigned int i = i0 + j;
          const float x = hu[i] - fum;
          const float qp = fqm + (bq1 + bq2 * x) * x;
          float dp = hphi[i] - fpm;
          dp = dp > kPi ? dp - k2Pi : dp;
          dp = dp < -kPi ? dp + k2Pi : dp;
          dp -= (bp1 + bp2 * x) * x;
          const float lev = isr > 0 ? (fsm + (bs1 + bs2 * x) * x) * isr : 1.0f;
          m[j] = (std::abs(hq[i] - qp) <= kSlackD * (aq + lev * bq)) & (std::abs(dp) <= kSlackD * (ap + lev * bp));
          if (CK.on) {
            double qe, pe;
            if (hx.predict(D.disc, D.qbar(i), qe, pe)) {
              ++CK.cands;
              const double pq = fqm + (bq1 + bq2 * x) * x, pphi = pm + (bp1 + bp2 * x) * x;
              CK.max_rel_q = std::max(CK.max_rel_q, std::abs(pq - qe) / (aq + lev * bq));
              CK.max_rel_phi = std::max(CK.max_rel_phi, std::abs(wrap(pphi - pe)) / (ap + lev * bp));
            }
          }
        }
        for (unsigned int j = 0; j < n; ++j) {
          cnt.d_touched++;
          if (!m[j])
            continue;
          // confirm: surf_stage_d's prediction and cut, in double
          const unsigned int kd = i0 + j;
          const P3 hd = surf_p3(D, kd);
          double q, phi, s = -1;
          if (!hx.predict(D.disc, D.qbar(kd), q, phi, &s))
            continue;
          const double wq = P.s_ref > 0 ? P.wq_d(pte, s) : wqd, wp = P.s_ref > 0 ? P.wphi_d(pte, s) : wpd;
          if (!P.no_cut && (std::abs((D.disc ? hd.r() : hd.z) - q) > wq || std::abs(wrap(hd.phi() - phi)) > wp))
            continue;
          cnt.quads++;
          f(kd, hd);
        }
      }
    });
  }

  // The finder for one pattern.  L[0..3]: the pattern's layers, filled.
  inline void find_quads_surf(const SurfParams &P, const SurfLayer *L[4], std::vector<Quad> &out, SeedCounters &cnt) {
    const SurfLayer &A = *L[0], &Bl = *L[1], &C = *L[2], &D = *L[3];
    for (unsigned int ka = 0; ka < A.n(); ++ka) {
      const surf::P3 ha = surf_p3(A, ka);
      surf_stage_b(P, ha, Bl, cnt, [&](unsigned int kb, const surf::P3 &hb) {
        surf_stage_c(P, ha, hb, C, cnt, [&](unsigned int kc, const surf::P3 &hc) {
          const surf::Helix hx(ha, hb, hc);
          surf_stage_d(P, hx, D, cnt, [&](unsigned int kd, const surf::P3 &) {
            out.push_back({A.orig_[ka], Bl.orig_[kb], C.orig_[kc], D.orig_[kd]});
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
  struct SurfOwnership : SeedLayerEnvelopes {
    // the envelopes (Env, env, delta), setup(), cross() and is_pix() are SeedLayerEnvelopes'
    int max_skip = 0;  // definite layers crossed before the pattern's last one and not in it

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
  // The double-precision reference finder of the chain: SeedChain's configuration (MkFitCore) and the run.
  struct SurfChain : SeedChain {
    int fast = 0;  // 1: the float kernels (K2) where they exist; 0: the double reference
    long n_forwarded = 0, n_dropped = 0;
    double t_start = 0;                       // seconds in the start doublets
    unsigned long long cyc_c = 0, cyc_d = 0;  // rdtsc cycles in stage c (with its pushes) and stage d
    long n_c_cand = 0, n_d_cand = 0;          // queued candidates searched with stage c, stage d

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
      using sclk = std::chrono::steady_clock;
      auto tsec = [](sclk::time_point a, sclk::time_point b) { return std::chrono::duration<double>(b - a).count(); };
      const auto t0 = sclk::now();
      // start doublets
      for (const auto &se : starts) {
        const SurfLayer *A = layer(se.first), *B = layer(se.second);
        if (!A || !B)
          continue;
        for (unsigned int ka = 0; ka < A->n(); ++ka) {
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
      t_start += tsec(t0, sclk::now());
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
            const unsigned long long tc0 = phases ? __builtin_ia32_rdtsc() : 0;
            ++n_c_cand;
            const SurfParams &Pc = par_[ix >= 0 ? ix : 0];
            const surf::P3 h0 = hit(c, 0), h1 = hit(c, 1);
            surf_stage_c(Pc, h0, h1, *T, cnt, [&](unsigned int kc, const surf::P3 &hc) {
              found = true;
              Cand t = c;
              t.nh = 3, t.pos[2] = p, t.k[2] = kc;
              t.ln = LineRZ(h0, hc);
              push_next(t, p);
            });
            if (phases)
              cyc_c += __builtin_ia32_rdtsc() - tc0;
          } else if (c.nh == 3) {
            const unsigned long long td0 = phases ? __builtin_ia32_rdtsc() : 0;
            ++n_d_cand;
            const SurfParams &Pd = par_[ix >= 0 ? ix : 0];
            const surf::Helix hx(hit(c, 0), hit(c, 1), hit(c, 2));
            auto on_d = [&](unsigned int kd, const surf::P3 &) {
              found = true;
              std::array<int, 4> ids{order[c.pos[0]], order[c.pos[1]], order[c.pos[2]], order[p]};
              const SurfLayer *La = layer(c.pos[0]), *Lb = layer(c.pos[1]), *Lc = layer(c.pos[2]);
              out.push_back({ids, {La->orig_[c.k[0]], Lb->orig_[c.k[1]], Lc->orig_[c.k[2]], T->orig_[kd]}});
            };
            if (fast)
              surf_stage_d_fast(Pd, hx, *T, cnt, on_d);
            else
              surf_stage_d(Pd, hx, *T, cnt, on_d);
            if (phases)
              cyc_d += __builtin_ia32_rdtsc() - td0;
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
