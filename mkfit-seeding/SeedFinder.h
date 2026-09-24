#ifndef mkfit_seeding_SeedFinder_h
#define mkfit_seeding_SeedFinder_h

// Pure-geometry quadruplet finder on the mkFit binnor.
//
// STEP 1 of milestone 1: a scalar port of the prototype's search
// (mkfit-standalone-seedgeom/SeedGeom.cc, part_b) in which the per-hit cuts are
// copied expression for expression -- same float/double types, same order of
// operations -- and only the FETCH is new.  The acceptance test is that the
// quad list equals the prototype's in cover mode, as a set, for identical
// parameters.  Optimisation comes after that, one step at a time, each step
// re-checked against the same list.

#include "SeedLayer.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace mkfit::seeding {

  struct SeedParams {
    float pt_min = 0.9f;       // GeV
    float d0_max = 0.1f;       // cm
    float qwin = 0.035f;       // cm, 3rd-layer z tolerance
    float qwin_d = 0.025f;     // cm, 4th-layer z tolerance
    float phiwin_d = 0.002f;   // rad, 4th-layer phi tolerance
    float phi_margin = 0.002f; // rad
    int phi_lin = 0;           // 0 off, 1 linear only, 2 linear AND generic band
    float phi_lin_marg = 0.003f;
  };

  struct SeedCounters {
    long doublets = 0, c_touched = 0, triplets = 0, d_touched = 0, quads = 0;
    void add(const SeedCounters &o) {
      doublets += o.doublets;
      c_touched += o.c_touched;
      triplets += o.triplets;
      d_touched += o.d_touched;
      quads += o.quads;
    }
  };

  // (ia, ib, ic, id) in ORIGINAL hit indices within each layer's HitVec
  using Quad = std::array<unsigned int, 4>;

  namespace detail {
    constexpr float kTwoPi = 2.0f * kPi;
    constexpr float kCoverEps = 1e-4f;  // cm, as the prototype's cover mode

    inline float wrap_pi(float d) {
      d -= kTwoPi * (d > kPi);
      d += kTwoPi * (d < -kPi);
      return d;
    }

    // verbatim from the prototype, including the float parameters: callers pass
    // doubles and they are narrowed here, exactly as there
    inline bool circle3(float x1, float y1, float x2, float y2, float x3, float y3, double &cx, double &cy, double &R) {
      double A = x2 - x1, B = y2 - y1, C = x3 - x1, D = y3 - y1;
      double E = A * (x1 + x2) + B * (y1 + y2);
      double F = C * (x1 + x3) + D * (y1 + y3);
      double G = 2.0 * (A * (y3 - y2) - B * (x3 - x2));
      if (std::abs(G) < 1e-12)
        return false;
      cx = (D * E - B * F) / G;
      cy = (A * F - C * E) / G;
      R = std::hypot(x1 - cx, y1 - cy);
      return R > 1e-6 && std::isfinite(R);
    }

    inline double arc(double chord, double R) {
      double h = chord / (2.0 * R);
      if (h > 1.0)
        h = 1.0;
      return 2.0 * R * std::asin(h);
    }

    inline bool circle_cross_r(
        double cx, double cy, double R, double rt, double xref, double yref, double &px, double &py) {
      const double dC = std::hypot(cx, cy);
      if (dC < 1e-9)
        return false;
      const double aa = (rt * rt + dC * dC - R * R) / (2 * dC);
      const double b2 = rt * rt - aa * aa;
      if (b2 < 0)
        return false;
      const double bb = std::sqrt(b2);
      const double ux = cx / dC, uy = cy / dC;
      const double vx = -uy, vy = ux;
      const double x1 = aa * ux + bb * vx, y1 = aa * uy + bb * vy;
      const double x2 = aa * ux - bb * vx, y2 = aa * uy - bb * vy;
      if (std::hypot(x1 - xref, y1 - yref) <= std::hypot(x2 - xref, y2 - yref)) {
        px = x1;
        py = y1;
      } else {
        px = x2;
        py = y2;
      }
      return true;
    }
  }  // namespace detail

  // LA: layer a, iterated whole.  LB: doublet target, queried in phi only.
  // LC, LD: queried in (phi, q).
  template <typename LA, typename LB, typename LC, typename LD>
  void find_quads(const SeedParams &P,
                  const LA &ga,
                  const LB &gb,
                  const LC &gc,
                  const LD &gd,
                  std::vector<Quad> &out,
                  SeedCounters &cnt) {
    using namespace detail;
    constexpr float kBfield = 3.8f;
    const double Rmin = P.pt_min / (0.003f * kBfield);
    const float inv2R = (float)(1.0 / (2 * Rmin)), d0m = P.d0_max, marg = P.phi_margin;
    const double rbhi = gb.rhi_, rclo = gc.rlo_, rchi = gc.rhi_, rdlo = gd.rlo_, rdhi = gd.rhi_;
    const float inv_rbhi = 1.0f / (float)rbhi, inv_rchi = 1.0f / (float)rchi;
    auto wphi_w = [=](float r_in, float inv_in, float r_out, float inv_out) {
      return (r_out - r_in) * inv2R + d0m * (inv_in - inv_out) + marg;
    };
    // a window of half-width w below pi never produces the ambiguous begin == end
    auto phi_bins = [](const auto &L, float c, float w) {
      w = std::min(w, 0.9f * kPi);
      return L.phi_range(c - w, c + w);
    };

    for (unsigned int ka = 0; ka < ga.n(); ++ka) {
      const float pa = ga.phi_[ka], za = ga.z_[ka], rra = ga.r_[ka];
      const float inva = ga.invr_[ka];
      const double xa = rra * std::cos(pa), ya = rra * std::sin(pa);
      const float w_ab = wphi_w(rra, inva, (float)rbhi, inv_rbhi);

      gb.for_each_in(phi_bins(gb, pa, w_ab), gb.q_all(), [&](unsigned int kb) {
        const float rrb = gb.r_[kb];
        if (rrb <= rra + 0.1f)
          return;
        const float pbph = gb.phi_[kb], zb = gb.z_[kb];
        const float invb = gb.invr_[kb];
        if (std::abs(wrap_pi(pbph - pa)) > wphi_w(rra, inva, rrb, invb))
          return;
        cnt.doublets++;

        const double cot = (zb - za) / (rrb - rra);
        const float slope = wrap_pi(pbph - pa) / (rrb - rra);
        const float kab = 1.0f / (rrb - rra);
        const float dinv_ab = invb - inva;
        const float mlin = P.phi_lin_marg;
        auto wlin = [d0m, invb, rrb, kab, dinv_ab, mlin](float r_c, float inv_c) {
          return d0m * std::abs(inv_c - invb - (r_c - rrb) * kab * dinv_ab) + mlin;
        };
        float w_bc, phic_ctr;
        if (P.phi_lin) {
          const float rmid = 0.5f * (float)(rclo + rchi);
          phic_ctr = pbph + slope * (rmid - rrb);
          w_bc = std::max(wlin((float)rclo, 1.0f / (float)rclo), wlin((float)rchi, inv_rchi)) +
                 std::abs(slope) * 0.5f * (float)(rchi - rclo);
        } else {
          phic_ctr = pbph;
          w_bc = wphi_w(rrb, invb, (float)rchi, inv_rchi);
        }
        // z range over the layer's radial extent, linear in r, plus the tolerance
        const double z0 = za + cot * (rclo - rra), z1 = za + cot * (rchi - rra);
        const double zh = P.qwin + kCoverEps;
        const auto qr = gc.q_range((float)(std::min(z0, z1) - zh), (float)(std::max(z0, z1) + zh));

        gc.for_each_in(phi_bins(gc, phic_ctr, w_bc), qr, [&](unsigned int kc) {
          cnt.c_touched++;
          const float rrc = gc.r_[kc], pcph = gc.phi_[kc];
          if (rrc <= rrb + 0.1f)
            return;
          if (P.phi_lin) {
            if (std::abs(wrap_pi(pcph - pbph - slope * (rrc - rrb))) > wlin(rrc, gc.invr_[kc]))
              return;
          }
          if (P.phi_lin != 1 && std::abs(wrap_pi(pcph - pbph)) > wphi_w(rrb, invb, rrc, gc.invr_[kc]))
            return;
          const double zpred = za + cot * (rrc - rra);
          const double dz = gc.z_[kc] - zpred;
          if (std::abs(dz) > P.qwin)
            return;
          cnt.triplets++;

          const double xb = rrb * std::cos(pbph), yb = rrb * std::sin(pbph);
          const double xc = rrc * std::cos(pcph), yc = rrc * std::sin(pcph);
          double cx = 0, cy = 0, R;
          if (!circle3(xa, ya, xb, yb, xc, yc, cx, cy, R))
            R = 1e6;
          const double s_ab = arc(std::hypot(xb - xa, yb - ya), R);
          const double s_bc = arc(std::hypot(xc - xb, yc - yb), R);
          const double s_ac = s_ab + s_bc;
          const double cot_s = (s_ac > 1e-6) ? (gc.z_[kc] - za) / s_ac : 0.0;

          // 4th layer: phi and z ranges from the predictions at the layer's two
          // radial edges, so every hit the per-hit cut can accept is fetched
          double p0x, p0y, p1x, p1y;
          const bool ok0 = circle_cross_r(cx, cy, R, rdlo, xc, yc, p0x, p0y);
          const bool ok1 = circle_cross_r(cx, cy, R, rdhi, xc, yc, p1x, p1y);
          if (!ok0 && !ok1)
            return;
          const double f0 = ok0 ? std::atan2(p0y, p0x) : std::atan2(p1y, p1x);
          const double f1 = ok1 ? std::atan2(p1y, p1x) : f0;
          const double dfh = 0.5 * wrap_pi((float)(f1 - f0));
          const double phid = f0 + dfh;
          const double phid_half = P.phiwin_d + std::abs(dfh) + kCoverEps;
          const double zz0 = ok0 ? gc.z_[kc] + cot_s * arc(std::hypot(p0x - xc, p0y - yc), R)
                                 : gc.z_[kc] + cot_s * arc(std::hypot(p1x - xc, p1y - yc), R);
          const double zz1 = ok1 ? gc.z_[kc] + cot_s * arc(std::hypot(p1x - xc, p1y - yc), R) : zz0;
          const double zdh = P.qwin_d + kCoverEps;
          const auto qd = gd.q_range((float)(std::min(zz0, zz1) - zdh), (float)(std::max(zz0, zz1) + zdh));

          gd.for_each_in(phi_bins(gd, (float)phid, (float)phid_half), qd, [&](unsigned int kd) {
            cnt.d_touched++;
            const float rrd = gd.r_[kd];
            if (rrd <= rrc + 0.1f)
              return;
            double px2, py2;
            if (!circle_cross_r(cx, cy, R, rrd, xc, yc, px2, py2))
              return;
            const double dphi_d = wrap_pi((float)(gd.phi_[kd] - std::atan2(py2, px2)));
            const double zd2 = gc.z_[kc] + cot_s * arc(std::hypot(px2 - xc, py2 - yc), R);
            const double dz_d = gd.z_[kd] - zd2;
            if (std::abs(dphi_d) > P.phiwin_d || std::abs(dz_d) > P.qwin_d)
              return;
            cnt.quads++;
            out.push_back({ga.orig_[ka], gb.orig_[kb], gc.orig_[kc], gd.orig_[kd]});
          });
        });
      });
    }
  }

}  // namespace mkfit::seeding

#endif
