#ifndef mkfit_seeding_SeedSurfBatch_h
#define mkfit_seeding_SeedSurfBatch_h

// The feed-forward chain in float, batched: SurfChain's logic (SeedSurf.h) with
// its configuration -- the chain order, the start pairs, the window tables and
// the crossing envelopes are all taken from a set-up SurfChain, which stays the
// double-precision reference.  Accepted by `seedsurf --margins` against a
// reference list, not by list identity.
//
// What differs from SurfChain::run:
//   - a candidate is 28 bytes: three hit indices, three chain positions, holes,
//     the target's crossing state and its r-z line (z0, cot) in float;
//   - each stage takes the candidates of one target layer in blocks, and the
//     r-z crossing tests (holes before b, the next crossed layer) run over the
//     block as lanes, one layer at a time;
//   - stages b, c and d cut in float on the layer's struct-of-arrays, with no
//     double confirm.
//
// Stage d predicts the helix at each hit's own qbar directly, from one point:
// the circle through a, b, c as curvature k, point c and the tangent at c.
//   barrel, |P| = r   the radical line of |P| = r and the helix circle,
//                     multiplied through by k: P.(k c + n) = k (r^2 + |c|^2) / 2 + c.n,
//                     with n the left normal of the tangent.  Every term is O(1)
//                     for any k, so nothing cancels for a stiff track, and the
//                     helix centre (1/k away) never appears.  Of the two points,
//                     the one ahead of c with the shorter chord L; the arc is
//                     2 asin(|k| L / 2) / |k|, and z = z_c + cot * arc.
//   disc, z = u       the arc length is (u - z_c) / cot; the point is c plus the
//                     chord ds sinc(h) along the tangent turned by h = k ds / 2.
// The phi cut is |dphi| < w taken as dot > 0 and cross^2 < sin^2(w) |P|^2 |h|^2,
// on the hit's own x, y, so no atan2 per hit.

#include "SeedSurf.h"

#include <chrono>
#include <cmath>
#include <vector>

namespace mkfit::seeding {

  namespace surfb {
    constexpr float k2Pi = 6.28318530717959f, kInv2Pi = 0.159154943091895f;  // kPi from SeedLayer.h
    inline float wrap(float d) { return d - k2Pi * std::floor(d * kInv2Pi + 0.5f); }
    // asin(x) / x, for 0 <= x <= 1; the series below x = 0.2 (next term < 1e-10)
    inline float asin_ox(float x) {
      const float x2 = x * x;
      if (x2 < 0.04f)
        return 1 + x2 * (1.f / 6 + x2 * (3.f / 40 + x2 * (5.f / 112 + x2 * (35.f / 1152 + x2 * (63.f / 2816)))));
      return x >= 1 ? kPi / 2 : std::asin(x) / x;
    }
    // sin(h) / h, cos(h), sin(h)
    inline void sinc_cs(float h, float &sc, float &c, float &s) {
      const float h2 = h * h;
      if (h2 < 0.04f) {
        sc = 1 - h2 * (1.f / 6 - h2 * (1.f / 120 - h2 * (1.f / 5040)));
        c = 1 - h2 * (0.5f - h2 * (1.f / 24 - h2 * (1.f / 720 - h2 * (1.f / 40320))));
        s = h * sc;
      } else {
        s = std::sin(h), c = std::cos(h);
        sc = s / h;
      }
    }

    // The helix at hit c, from the circle through a, b, c (as surf::Helix, in float).
    struct HelixF {
      float k = 0, tx = 0, ty = 0, cx = 0, cy = 0, cz = 0, cot = 0;
      // the barrel solve: unit m = (k c + n) / |k c + n|, 1 / |k c + n|, |c|^2, c.n
      float ux = 0, uy = 0, im = 0, c2 = 0, cn = 0;
      bool ok = false;

      void make(float ax, float ay, float az, float bx, float by, float bz, float cx_, float cy_, float cz_) {
        (void)bz;
        cx = cx_, cy = cy_, cz = cz_;
        const float vx0 = bx - ax, vy0 = by - ay, vx = cx - bx, vy = cy - by;
        const float lu = std::sqrt(vx0 * vx0 + vy0 * vy0), lv = std::sqrt(vx * vx + vy * vy);
        const float wx = cx - ax, wy = cy - ay, lw = std::sqrt(wx * wx + wy * wy);
        ok = lu >= 1e-6f && lv >= 1e-6f;
        if (!ok)
          return;
        k = 2 * (vx0 * vy - vy0 * vx) / (lu * lv * lw);
        // the tangent at c: the chord b->c turned by half its central angle, sin(half) = k lv / 2
        const float sn = std::clamp(0.5f * k * lv, -1.0f, 1.0f), cs = std::sqrt(1 - sn * sn);
        const float ex = vx / lv, ey = vy / lv;
        tx = ex * cs - ey * sn;
        ty = ex * sn + ey * cs;
        const float ak = std::abs(k);
        const float s_ab = lu * asin_ox(std::min(1.0f, 0.5f * ak * lu)), s_bc = lv * asin_ox(std::min(1.0f, 0.5f * ak * lv));
        cot = (cz - az) / (s_ab + s_bc);
        const float mx = k * cx - ty, my = k * cy + tx;
        im = 1.0f / std::sqrt(mx * mx + my * my);
        ux = mx * im, uy = my * im;
        c2 = cx * cx + cy * cy;
        cn = -cx * ty + cy * tx;
      }
      float pt() const { return std::abs(k) > 1e-12f ? 0.0114f / std::abs(k) : 1e9f; }

      // barrel: the first crossing ahead of c with |P| = r; q = z, s3 = 3D path length
      bool at_r(float r, float &px, float &py, float &q, float &s3) const {
        const float A = (0.5f * k * (r * r + c2) + cn) * im, g2 = r * r - A * A;
        if (g2 < 0)
          return false;
        const float g = std::sqrt(g2);
        const float p1x = A * ux - g * uy, p1y = A * uy + g * ux, p2x = A * ux + g * uy, p2y = A * uy - g * ux;
        const float d1x = p1x - cx, d1y = p1y - cy, d2x = p2x - cx, d2y = p2y - cy;
        const bool a1 = d1x * tx + d1y * ty > 0, a2 = d2x * tx + d2y * ty > 0;
        const float L1 = d1x * d1x + d1y * d1y, L2 = d2x * d2x + d2y * d2y;
        float L2s;
        if (a1 && (!a2 || L1 <= L2))
          px = p1x, py = p1y, L2s = L1;
        else if (a2)
          px = p2x, py = p2y, L2s = L2;
        else
          return false;
        const float L = std::sqrt(L2s), s = L * asin_ox(std::min(1.0f, 0.5f * std::abs(k) * L));
        q = cz + cot * s;
        s3 = s * std::sqrt(1 + cot * cot);
        return s > 0;
      }
      // disc: the point at z = u; q = r
      bool at_z(float u, float &px, float &py, float &q, float &s3) const {
        if (std::abs(cot) < 1e-9f)
          return false;
        const float ds = (u - cz) / cot;
        if (ds <= 0)
          return false;
        float sc, ch, sh;
        sinc_cs(0.5f * k * ds, sc, ch, sh);
        const float f = ds * sc;
        px = cx + f * (ch * tx - sh * ty);
        py = cy + f * (sh * tx + ch * ty);
        q = std::sqrt(px * px + py * py);
        s3 = ds * std::sqrt(1 + cot * cot);
        return true;
      }
      bool at(bool disc, float u, float &px, float &py, float &q, float &s3) const {
        return disc ? at_z(u, px, py, q, s3) : at_r(u, px, py, q, s3);
      }

      // --chain-batch-d 1: Hermite3D's one-point mode, the Taylor cubic of the
      // helix about c in the transverse arc s,
      //   P(s) = c + t (s - k^2 s^3 / 6) + n k s^2 / 2,
      // truncation ~ R_c (k s)^4 / 24. Barrel: |P(s)| = r by Newton from the
      // straight line; disc: s from z, which is exact (z is linear in s).
      void cubic1(float s, float &x, float &y, float &dx, float &dy) const {
        const float nx = -ty, ny = tx, ks = k * s;
        const float a = s * (1 - ks * ks * (1.f / 6)), b = 0.5f * ks * s;
        x = cx + a * tx + b * nx, y = cy + a * ty + b * ny;
        const float da = 1 - 0.5f * ks * ks;
        dx = da * tx + ks * nx, dy = da * ty + ks * ny;
      }
      bool at1(bool disc, float u, float &px, float &py, float &q, float &s3) const {
        float s, dx, dy;
        if (disc) {
          if (std::abs(cot) < 1e-9f)
            return false;
          s = (u - cz) / cot;
          if (s <= 0)
            return false;
          cubic1(s, px, py, dx, dy);
          q = std::sqrt(px * px + py * py);
        } else {
          const float ct = cx * tx + cy * ty, d = ct * ct - c2 + u * u;
          if (d < 0)
            return false;
          s = -ct + std::sqrt(d);
          for (int it = 0; it < 3; ++it) {
            cubic1(s, px, py, dx, dy);
            s -= (px * px + py * py - u * u) / (2 * (px * dx + py * dy));
          }
          cubic1(s, px, py, dx, dy);
          q = cz + cot * s;
        }
        s3 = s * std::sqrt(1 + cot * cot);
        return s > 0;
      }
    };

    // --chain-batch-d 2: Hermite3D's two-point mode across the target slab. The
    // helix at the slab's two qbar edges (at(), exact), the cubic in t through
    // both points with the tangents scaled by the transverse arc L between them,
    // truncation ~ R_c (k L)^4 / 384. z and the arc are linear in t. Barrel:
    // |H(t)| = r by Newton from t linear in r; disc: t from z, exact.
    struct Hermite2F {
      float ax[4], ay[4], z0 = 0, dz = 0, s0 = 0, L = 0, u0 = 0, du = 0, sq = 1;
      bool ok = false;
      void make(const HelixF &H, bool disc, float u0_, float u1_) {
        float x0, y0, q0, a0, x1, y1, q1, a1;
        ok = H.ok && H.at(disc, u0_, x0, y0, q0, a0) && H.at(disc, u1_, x1, y1, q1, a1);
        if (!ok)
          return;
        sq = std::sqrt(1 + H.cot * H.cot);
        s0 = a0 / sq, L = (a1 - a0) / sq;
        ok = L > 1e-4f;
        if (!ok)
          return;
        // tangents: t turned by k s
        float sc, c0, sn0, c1, sn1;
        sinc_cs(H.k * s0, sc, c0, sn0);
        sinc_cs(H.k * (s0 + L), sc, c1, sn1);
        const float t0x = H.tx * c0 - H.ty * sn0, t0y = H.tx * sn0 + H.ty * c0;
        const float t1x = H.tx * c1 - H.ty * sn1, t1y = H.tx * sn1 + H.ty * c1;
        auto coef = [&](float p0, float d0, float p1, float d1, float *a) {
          a[0] = p0, a[1] = d0, a[2] = 3 * (p1 - p0) - 2 * d0 - d1, a[3] = 2 * (p0 - p1) + d0 + d1;
        };
        coef(x0, L * t0x, x1, L * t1x, ax);
        coef(y0, L * t0y, y1, L * t1y, ay);
        z0 = H.cz + H.cot * s0, dz = H.cot * L;
        u0 = u0_, du = u1_ - u0_;
      }
      void eval(float t, float &x, float &y, float &dx, float &dy) const {
        x = ((ax[3] * t + ax[2]) * t + ax[1]) * t + ax[0];
        y = ((ay[3] * t + ay[2]) * t + ay[1]) * t + ay[0];
        dx = (3 * ax[3] * t + 2 * ax[2]) * t + ax[1];
        dy = (3 * ay[3] * t + 2 * ay[2]) * t + ay[1];
      }
      bool at(bool disc, float u, float &px, float &py, float &q, float &s3) const {
        float t, dx, dy;
        if (disc) {
          t = (u - z0) / dz;
          eval(t, px, py, dx, dy);
          q = std::sqrt(px * px + py * py);
        } else {
          t = (u - u0) / du;
          for (int it = 0; it < 3; ++it) {
            eval(t, px, py, dx, dy);
            t -= (px * px + py * py - u * u) / (2 * (px * dx + py * dy));
          }
          eval(t, px, py, dx, dy);
          q = z0 + dz * t;
        }
        s3 = (s0 + L * t) * sq;
        return s3 > 0;
      }
    };
  }  // namespace surfb

  struct SurfChainBatch {
    struct Cand {
      unsigned int k[3];
      unsigned char pos[3], holes, st;
      float z0, cot;  // the r-z line through its first and last hit
      float sc;       // a triplet: the stage c part of the residual score (fk_score)
    };
    // the crossing envelope of one chain position, thresholds with the margin folded in
    struct EnvF {
      bool ok = false, disc = false, pix = false;
      float pos = 0, plo = 0, phi = 0, lo_in = 0, hi_in = 0, lo_out = 0, hi_out = 0;
    };
    // one entry of SurfChain::par_, in float
    struct ParF {
      float phi_c, q_c, aphi, bphi, aq, bq, sref;
      int e0, ne;  // |eta| slices in etaw_
    };

    SurfChain *C = nullptr;
    int n = 0;
    float d0_max = 0, pt_min = 0;
    std::vector<EnvF> env_;
    std::vector<ParF> par_;
    std::vector<SurfParams::EtaWin> etaw_;
    std::vector<std::vector<std::pair<int, int>>> hole_pos_;  // per start pair: (position, between a and b)
    std::vector<std::vector<Cand>> Q2_, Q3_;                  // doublets, triplets, per target position

    long n_forwarded = 0, n_dropped = 0, n_c_cand = 0, n_d_cand = 0;
    double t_start = 0;
    unsigned long long cyc_c = 0, cyc_d = 0;
    bool phases = false;
    int d_mode = 0;  // stage d prediction: 0 direct from c, 1 one-point cubic, 2 two-point Hermite across the slab

    // Fake rejection (README "Fakes"), each off by default.
    // fk_score > 0: a quad needs its residual score, (dq_c/w)^2 + (dphi_c/w)^2 + (dphi_d/w)^2 + (dq_d/w)^2
    //   with each residual over its own window, below fk_score. Stage c drops a triplet whose c part
    //   alone reaches it, and stage d adds the d part.
    float fk_score = 0;
    // fk_shape: the cluster length along z (Hit::spanCols()) of every hit on a barrel pixel layer (0-3)
    //   within the band of true hits for the line's |cot theta|: at stage b on the a-b line, at c and d
    //   on the line the candidate already has (a-b, a-c).
    struct ShapeTab {
      float inv_bw = 0;
      std::vector<int> lo, hi;  // per |cot| bin
      bool ok() const { return !lo.empty(); }
      int bin(float acot) const { return std::min((int)(acot * inv_bw), (int)lo.size() - 1); }
      int pass(int span, float acot) const {
        const int b = bin(acot);
        return (span >= lo[b]) & (span <= hi[b]);
      }
    };
    ShapeTab shape_[4];
    bool fk_shape = false;
    // fk_ot2 > 0: a quad whose d is on OT1-P (layer 4) needs a hit on OT2-P (layer 6) for the helix through
    //   b, c, d: the best hit by the score of next_hit() within fk_ot2 x (a + b / pT) in phi and z. Quads
    //   whose helix does not reach OT2-P, or reaches it outside |z| < zmax - 2 cm, pass.
    float fk_ot2 = 0;
    // q97 of true quads, a + b / pT in rad and cm, events 0-39 (windows-D121 sample)
    float ot2_aphi = -8.1e-4f, ot2_bphi = 5.92e-3f, ot2_aq = 0.3084f, ot2_bq = 0.0909f;
    long n_fk_shape = 0, n_fk_ot2 = 0, n_ot2_tested = 0;

    static constexpr int kBlk = 256;

    void setup(SurfChain &c) {
      C = &c;
      if (!c.par_built_)
        c.build_params();
      n = c.order.size();
      d0_max = c.P.d0_max, pt_min = c.P.pt_min;
      phases = c.phases;
      const SurfOwnership &o = *c.own;
      const double dl = o.delta;
      env_.assign(n, EnvF{});
      for (int p = 0; p < n; ++p) {
        EnvF &e = env_[p];
        e.pix = SurfOwnership::is_pix(c.order[p]);
        if (c.env_of[p] < 0)
          continue;
        const SurfOwnership::Env &E = o.env[c.env_of[p]];
        e.ok = true, e.disc = E.disc, e.pos = E.pos, e.plo = E.pos_lo, e.phi = E.pos_hi;
        e.lo_in = E.lo + dl, e.hi_in = E.hi - dl, e.lo_out = E.lo - dl, e.hi_out = E.hi + dl;
      }
      par_.clear();
      etaw_.clear();
      for (const SurfParams &Q : c.par_) {
        par_.push_back({Q.phi_c, Q.q_c, Q.phi_d, Q.b_phi_d, Q.q_d, Q.b_q_d, Q.s_ref, (int)etaw_.size(), Q.n_eta_win});
        for (int i = 0; i < Q.n_eta_win; ++i)
          etaw_.push_back(Q.eta_win[i]);
      }
      hole_pos_.assign(c.starts.size(), {});
      for (size_t si = 0; si < c.starts.size(); ++si)
        for (int q = 0; q < c.starts[si].second; ++q)
          if (q != c.starts[si].first && env_[q].pix && env_[q].ok)
            hole_pos_[si].push_back({q, q > c.starts[si].first});
      Q2_.assign(n, {});
      Q3_.assign(n, {});
    }

    // The work arrays of a block: one r-z line per lane, as (z0, cot) for a barrel
    // crossing and as (1/cot, -z0/cot) for a disc crossing (0, 0 for |cot| < 1e-9),
    // so the crossing tests run over the lanes with no division and no branch.
    // Everything per lane is float or int32, which the vectorizer takes with -mavx.
    std::vector<float> w_z0, w_cot, w_ic, w_zr;
    std::vector<int> w_allow, w_st, w_holes, w_bt, w_s, w_sa, w_sb, w_nq;
    std::vector<unsigned int> w_i0, w_i1;  // per lane: hit or candidate indices (meaning per stage)
    void work(size_t m) {
      if (w_z0.size() < m) {
        for (auto *v : {&w_z0, &w_cot, &w_ic, &w_zr})
          v->resize(m);
        for (auto *v : {&w_allow, &w_st, &w_holes, &w_bt, &w_s, &w_sa, &w_sb, &w_nq})
          v->resize(m);
        w_i0.resize(m), w_i1.resize(m);
      }
    }
    void prep_lines(int m) {
      float *__restrict ic = w_ic.data(), *__restrict zr = w_zr.data();
      const float *__restrict z0 = w_z0.data(), *__restrict ct = w_cot.data();
      for (int i = 0; i < m; ++i) {
        const float c = ct[i], cs = std::abs(c) >= 1e-9f ? c : 1.0f, v = std::abs(c) >= 1e-9f ? 1.0f / cs : 0.0f;
        ic[i] = v, zr[i] = -z0[i] * v;
      }
    }
    // 0 no, 1 maybe, 2 definite; the definite band lies inside the maybe band (margin >= 0)
    static int cls(const EnvF &e, float x1, float x2) {
      const float xl = std::min(x1, x2), xh = std::max(x1, x2);
      return ((xh > e.lo_out) & (xl < e.hi_out)) + ((xl > e.lo_in) & (xh < e.hi_in));
    }
    // crossing state of position e for the m work lines
    void states(const EnvF &e, int m, int *__restrict st) const {
      const float *__restrict z0 = w_z0.data(), *__restrict ct = w_cot.data(), *__restrict ic = w_ic.data(),
                  *__restrict zr = w_zr.data();
      const float plo = e.plo, phi = e.phi, pos = e.pos;
      if (!e.disc) {
        for (int i = 0; i < m; ++i)
          st[i] = cls(e, z0[i] + ct[i] * plo, z0[i] + ct[i] * phi);
      } else {
        for (int i = 0; i < m; ++i) {
          const int valid = (ic[i] != 0) & (pos * ic[i] + zr[i] > 0);
          st[i] = valid * cls(e, plo * ic[i] + zr[i], phi * ic[i] + zr[i]);
        }
      }
    }

    // For the m work lines leaving position p: the next crossed position w_nq (-1:
    // none) and its state w_st; w_allow: a non-pixel position may be taken
    // (SurfChain::next()). One position at a time over the lanes, until every lane has one.
    void route(int p, int m) {
      prep_lines(m);
      int *__restrict nq = w_nq.data(), *__restrict nst = w_st.data(), *__restrict s = w_s.data();
      const int *__restrict allow = w_allow.data();
      for (int i = 0; i < m; ++i)
        nq[i] = -1, nst[i] = 0;
      for (int q = p + 1; q < n && m > 0; ++q) {
        const EnvF &e = env_[q];
        if (!e.ok)
          continue;
        states(e, m, s);
        const int pix = e.pix;
        int open = 0;
        for (int i = 0; i < m; ++i) {
          const int take = (nq[i] < 0) & (s[i] > 0) & (pix | allow[i]);
          nq[i] = take ? q : nq[i];
          nst[i] = take ? s[i] : nst[i];
        }
        for (int i = 0; i < m; ++i)
          open += nq[i] < 0;
        if (!open)
          return;
      }
    }

    // The hit on barrel layer LP closest to the helix H (through b, c, d), by the score
    // (dphi / sphi)^2 + (dz / sz)^2 with sphi = 0.0005 + 3.2e-3 / pT and sz = 0.075 + 0.0316 / pT (pT >= 0.5),
    // among the hits fetched within fphi and fz of the prediction at the layer's two edges. Returns the hit (in the
    // layer's sorted order), -1 if none, -2 if H does not reach the layer or reaches it outside its z range; zmid: the helix's z at the middle
    // of the two edges; dphi, dz, score: the best hit's.
    static int next_hit(const SurfLayer &LP, const surfb::HelixF &H, float pte, float fphi, float fz, float &zmid,
                        float &dphi, float &dz, float &score) {
      using namespace surfb;
      float x0, y0, q0, s0, x1, y1, q1, s1;
      const bool ok0 = H.ok && H.at_r(LP.qbar_lo, x0, y0, q0, s0), ok1 = H.ok && H.at_r(LP.qbar_hi, x1, y1, q1, s1);
      if (!ok0 && !ok1)
        return -2;
      if (!ok0)
        x0 = x1, y0 = y1, q0 = q1;
      if (!ok1)
        x1 = x0, y1 = y0, q1 = q0;
      zmid = 0.5f * (q0 + q1);
      if (std::min(q0, q1) > LP.q_hi || std::max(q0, q1) < LP.q_lo)
        return -2;  // outside the layer at both edges
      const float pt = std::max(0.5f, pte), sphi = 0.0005f + 3.2e-3f / pt, sz = 0.075f + 0.0316f / pt;
      const float p0 = std::atan2(y0, x0), p1 = std::atan2(y1, x1), dpp = wrap(p1 - p0);
      float best = 1e30f;
      int kb = -1;
      LP.sl.for_each_in(surf::phi_bins(LP, p0 + 0.5f * dpp, 0.5f * std::abs(dpp) + fphi),
                        surf::q_bins(LP, std::min(q0, q1) - fz, std::max(q0, q1) + fz), [&](unsigned int k) {
                          float px, py, qp, s3;
                          if (!H.at_r(LP.sl.r_[k], px, py, qp, s3))
                            return;
                          const float dp = wrap(LP.sl.phi_[k] - std::atan2(py, px)), dq = LP.sl.z_[k] - qp;
                          const float sc = (dp / sphi) * (dp / sphi) + (dq / sz) * (dq / sz);
                          if (sc < best)
                            best = sc, kb = k, dphi = dp, dz = dq;
                        });
      score = best;
      return kb;
    }

    void run(const std::map<int, const SurfLayer *> &L, std::vector<std::pair<std::array<int, 4>, Quad>> &out,
             SeedCounters &cnt) {
      using namespace surfb;
      const SurfChain &Ch = *C;
      std::vector<const SurfLayer *> lay(n);
      for (int p = 0; p < n; ++p) {
        auto it = L.find(Ch.order[p]);
        lay[p] = it == L.end() ? nullptr : it->second;
      }
      const int max_holes = Ch.max_holes, max_holes_ot = Ch.max_holes_ot;
      // the shape table of position p, or null
      auto shape_of = [&](int p) -> const ShapeTab * {
        const int l = Ch.order[p];
        return fk_shape && l >= 0 && l < 4 && shape_[l].ok() ? &shape_[l] : nullptr;
      };
      const float fks = fk_score > 0 ? fk_score : 1e30f;
      const bool fk_on = fk_score > 0 || fk_shape;
      const SurfLayer *lay_ot2 = nullptr;
      if (fk_ot2 > 0)
        if (auto it = L.find(6); it != L.end())
          lay_ot2 = it->second;
      const int nn = n;
      auto idx_c = [&](const Cand &c, int p) { return Ch.idx_c_[(c.pos[0] * nn + c.pos[1]) * nn + p]; };
      auto idx_d = [&](const Cand &c, int p) { return Ch.idx_d_[((c.pos[0] * nn + c.pos[1]) * nn + c.pos[2]) * nn + p]; };
      // push the m lanes of the work arrays routed from p: w_nq, w_st from route()
      auto push = [&](std::vector<std::vector<Cand>> &Q, int m, auto &&make) {
        for (int i = 0; i < m; ++i) {
          const int q = w_nq[i];
          if (q < 0 || !lay[q]) {
            ++n_dropped;
            continue;
          }
          Cand c = make(i);
          c.st = w_st[i];
          Q[q].push_back(c);
        }
      };

      using sclk = std::chrono::steady_clock;
      const auto t0 = sclk::now();

      // ---- start doublets
      // the stage b survivors of the current block: hits, and their r-z line computed in the mask loop
      static constexpr int kCap = kBlk + 64;
      alignas(32) unsigned int bka[kCap], bkb[kCap];
      alignas(32) float bz0[kCap], bct[kCap];
      const SurfParams &P = Ch.P;
      const int hs = Ch.start_holes < 0 ? max_holes : Ch.start_holes;
      long n_doublets = 0;
      for (size_t si = 0; si < Ch.starts.size(); ++si) {
        const int pa_ = Ch.starts[si].first, pb_ = Ch.starts[si].second;
        const SurfLayer *A = lay[pa_], *B = lay[pb_];
        if (!A || !B)
          continue;
        const auto &hp = hole_pos_[si];
        const ShapeTab *shA = shape_of(pa_), *shB = shape_of(pb_);
        int nb = 0;
        auto flush = [&]() {
          const int m = nb;
          if (!m)
            return;
          work(m);
          std::copy(bz0, bz0 + m, w_z0.data());
          std::copy(bct, bct + m, w_cot.data());
          prep_lines(m);
          int *__restrict holes = w_holes.data(), *__restrict bt = w_bt.data(), *__restrict s = w_s.data();
          for (int i = 0; i < m; ++i)
            holes[i] = 0, bt[i] = 0;
          // holes: definite crossings before b that the doublet does not use
          for (const auto &h : hp) {
            states(env_[h.first], m, s);
            const int btw = h.second;
            for (int i = 0; i < m; ++i) {
              const int d = s[i] == 2;
              holes[i] += d;
              bt[i] |= d & btw;
            }
          }
          states(env_[pa_], m, w_sa.data());
          states(env_[pb_], m, w_sb.data());
          const int lead_only = Ch.lead_only;
          int *__restrict shp = w_bt.data();  // bt is consumed here, reuse it
          for (int i = 0; i < m; ++i)
            shp[i] = !(lead_only & bt[i]);
          if (shA || shB) {
            for (int i = 0; i < m; ++i) {
              const float ac = std::abs(w_cot[i]);
              const int ok = (!shA || shA->pass(A->span_[bka[i]], ac)) & (!shB || shB->pass(B->span_[bkb[i]], ac));
              n_fk_shape += shp[i] & !ok;
              shp[i] &= ok;
            }
          }
          int g = 0;
          for (int i = 0; i < m; ++i) {
            const int good = shp[i] & (holes[i] <= hs) & (w_sa[i] > 0) & (w_sb[i] > 0);
            // compact the good lanes in place (g <= i)
            w_z0[g] = w_z0[i], w_cot[g] = w_cot[i], holes[g] = holes[i], w_allow[g] = 0;
            w_i0[g] = bka[i], w_i1[g] = bkb[i];
            g += good;
          }
          route(pb_, g);
          push(Q2_, g, [&](int i) {
            Cand c;
            c.k[0] = w_i0[i], c.k[1] = w_i1[i], c.k[2] = 0;
            c.pos[0] = pa_, c.pos[1] = pb_, c.pos[2] = 0;
            c.holes = w_holes[i];
            c.z0 = w_z0[i], c.cot = w_cot[i], c.sc = 0;
            return c;
          });
          nb = 0;
        };
        const float inv2R = 0.003f * 3.8f / (2.0f * P.pt_min), d0 = P.d0_max, marg = P.marg_b;
        const float zlo = P.bs_z - P.zv, zhi = P.bs_z + P.zv, side = Ch.side;
        const float *bphi = B->sl.phi_.data(), *br = B->sl.r_.data(), *bz = B->sl.z_.data(), *bir = B->sl.invr_.data();
        for (unsigned int ka = 0; ka < A->sl.n(); ++ka) {
          const surf::P3 ha = surf_p3(*A, ka);
          SurfFetch fe;
          if (!surf_b_fetch(P, ha, *B, fe))
            continue;
          // the float cuts of surf_stage_b_fast, the side, and the survivors' line
          const float ra = ha.r(), pa = ha.phi(), za = ha.z, inva = 1.0f / ra;
          B->sl.for_each_run(fe.p, fe.q, [&](unsigned int b, unsigned int e) {
            for (unsigned int i0 = b; i0 < e; i0 += 64) {
              const unsigned int nk = std::min(64u, e - i0);
              if (nb + (int)nk > kCap)
                flush();
              alignas(32) int msk[64];
              alignas(32) float lz0[64], lct[64];
              for (unsigned int j = 0; j < nk; ++j) {
                const unsigned int i = i0 + j;
                const float dr = br[i] - ra, dz = bz[i] - za;
                float dp = bphi[i] - pa;
                dp = dp > kPi ? dp - k2Pi : dp;
                dp = dp < -kPi ? dp + k2Pi : dp;
                const float w = dr * inv2R + d0 * (inva - bir[i]) + marg;
                const float num = za * dr - ra * dz;
                const int cut = (dr > 0.1f) & (std::abs(dp) <= w) & (num >= zlo * dr) & (num <= zhi * dr);
                // a doublet belongs to the side its line goes to
                const int sd = !((dz * side < 0) | ((dz == 0) & (side < 0)));
                msk[j] = cut + ((cut & sd) << 1);
                const float cot = dz / (dr > 0.1f ? dr : 1.0f);
                lct[j] = cot, lz0[j] = za - cot * ra;
              }
              int nd = 0, g = nb;
              for (unsigned int j = 0; j < nk; ++j) {
                nd += msk[j] & 1;
                bka[g] = ka, bkb[g] = i0 + j, bz0[g] = lz0[j], bct[g] = lct[j];
                g += msk[j] >> 1;
              }
              nb = g;
              n_doublets += nd;
            }
          });
        }
        flush();
      }
      cnt.doublets += n_doublets;
      t_start += std::chrono::duration<double>(sclk::now() - t0).count();

      // ---- the forward pass
      std::vector<unsigned int> t_j, t_k;  // stage c survivors: (candidate in block, hit)
      std::vector<float> t_s;              // ... and their stage c score
      std::vector<unsigned int> f_j;       // candidates forwarded as after a miss
      std::vector<HelixF> hx(kBlk);
      SurfFastCheck &CK = g_surf_fast_check;
      long n_ct = 0, n_tr = 0, n_dt = 0, n_qd = 0;
      for (int p = 0; p < n; ++p) {
        const SurfLayer *T = lay[p];
        if (!T)
          continue;
        const bool disc = T->disc;
        const float *hphi = T->sl.phi_.data(), *hz = T->sl.z_.data(), *hr = T->sl.r_.data(), *hir = T->sl.invr_.data();
        const float *hx_ = T->sl.x_.data(), *hy_ = T->sl.y_.data();
        const float *hu = disc ? hz : hr, *hq = disc ? hr : hz;
        const float u0 = T->qbar_lo, u1 = T->qbar_hi;
        const ShapeTab *shT = shape_of(p);
        const int *hsp = T->span_.data();

        // -- stage c: the doublets queued at p
        auto &Qd = Q2_[p];
        for (size_t b0 = 0; b0 < Qd.size(); b0 += kBlk) {
          const unsigned long long tc0 = phases ? __builtin_ia32_rdtsc() : 0;
          const size_t b1 = std::min(Qd.size(), b0 + kBlk);
          t_j.clear(), t_k.clear(), t_s.clear(), f_j.clear();
          for (size_t j = b0; j < b1; ++j) {
            const Cand &c = Qd[j];
            const int ix = idx_c(c, p);
            bool found = false;
            if (!Ch.known_only || ix >= 0) {
              ++n_c_cand;
              const ParF &w = par_[ix >= 0 ? ix : 0];
              const SurfLayer &La = *lay[c.pos[0]], &Lb = *lay[c.pos[1]];
              const unsigned ka = c.k[0], kb = c.k[1];
              const float ua = disc ? La.sl.z_[ka] : La.sl.r_[ka], ub = disc ? Lb.sl.z_[kb] : Lb.sl.r_[kb];
              const float qa = disc ? La.sl.r_[ka] : La.sl.z_[ka], qb = disc ? Lb.sl.r_[kb] : Lb.sl.z_[kb];
              const float pa = La.sl.phi_[ka], dp = wrap(Lb.sl.phi_[kb] - pa);
              const float ia = La.sl.invr_[ka], ib = Lb.sl.invr_[kb], idu = 1.0f / (ub - ua);
              const float t0_ = (u0 - ua) * idu, t1_ = (u1 - ua) * idu;
              if (std::max(t0_, t1_) > 1) {
                const float cq0 = qa + (qb - qa) * t0_, cq1 = qa + (qb - qa) * t1_;
                const float cp0 = pa + dp * t0_, cp1 = pa + dp * t1_;
                const bool far0 = std::abs(t0_) > std::abs(t1_);
                const float ufar = far0 ? u0 : u1, tfar = far0 ? t0_ : t1_;
                const float rfar = std::max(0.5f, disc ? qa + (qb - qa) * tfar : ufar);
                const float wcphi = d0_max * std::abs(1 / rfar - (ia + (ib - ia) * tfar)) + w.phi_c;
                const auto qc = surf::q_bins(*T, std::min(cq0, cq1) - w.q_c, std::max(cq0, cq1) + w.q_c);
                const float dcp = wrap(cp1 - cp0), cmid = cp0 + 0.5f * dcp, chalf = 0.5f * std::abs(dcp);
                const float qcw = w.q_c, pcw = w.phi_c, dqab = qb - qa, dib = ib - ia, iqcw = 1.0f / qcw;
                // the shape band of hit c, on the a-b line (barrel pixels only)
                int slo = 0, shi = 1 << 30;
                if (shT) {
                  const int sb = shT->bin(std::abs(dqab * idu));
                  slo = shT->lo[sb], shi = shT->hi[sb];
                }
                T->sl.for_each_run(surf::phi_bins(*T, cmid, chalf + wcphi + 1e-6f), qc, [&](unsigned int b, unsigned int e) {
                  for (unsigned int i0 = b; i0 < e; i0 += 64) {
                    const unsigned int nk = std::min(64u, e - i0);
                    unsigned char msk[64];
                    for (unsigned int jj = 0; jj < nk; ++jj) {
                      const unsigned int i = i0 + jj;
                      const float t = (hu[i] - ua) * idu;
                      const float dq = hq[i] - (qa + dqab * t);
                      const float dph = wrap(wrap(hphi[i] - pa) - dp * t);
                      const float d0t = std::abs(hir[i] - (ia + dib * t));
                      msk[jj] = (t > 1) & (std::abs(dq) <= qcw) & (std::abs(dph) <= d0_max * d0t + pcw);
                    }
                    n_ct += nk;
                    for (unsigned int jj = 0; jj < nk; ++jj)
                      if (msk[jj]) {
                        // the fake cuts on the survivors: the score, recomputed, and the shape
                        const unsigned int i = i0 + jj;
                        float sc = 0;
                        if (fk_on) {
                          const float t = (hu[i] - ua) * idu;
                          const float rq = (hq[i] - (qa + dqab * t)) * iqcw;
                          const float rp =
                              wrap(wrap(hphi[i] - pa) - dp * t) / (d0_max * std::abs(hir[i] - (ia + dib * t)) + pcw);
                          sc = rq * rq + rp * rp;
                          if (sc >= fks || hsp[i] < slo || hsp[i] > shi)
                            continue;
                        }
                        t_j.push_back(j), t_k.push_back(i), t_s.push_back(sc);
                        found = true;
                        ++n_tr;
                      }
                  }
                });
              }
            }
            if (!found || Ch.hole_always) {
              if (c.holes + (c.st == 2) <= max_holes) {
                ++n_forwarded;
                f_j.push_back(j);
              } else
                ++n_dropped;
            }
          }
          // the triplets: their a-c line, routed from p
          {
            const int m = t_j.size();
            work(m);
            for (int i = 0; i < m; ++i) {
              const Cand &c = Qd[t_j[i]];
              const SurfLayer &La = *lay[c.pos[0]];
              const float za = La.sl.z_[c.k[0]], ra = La.sl.r_[c.k[0]], zc = hz[t_k[i]], rc = hr[t_k[i]];
              const float cot = (zc - za) / (rc - ra);
              w_cot[i] = cot, w_z0[i] = za - cot * ra;
              w_allow[i] = c.holes <= max_holes_ot;
            }
            route(p, m);
            push(Q3_, m, [&](int i) {
              Cand c = Qd[t_j[i]];
              c.k[2] = t_k[i], c.pos[2] = p;
              c.z0 = w_z0[i], c.cot = w_cot[i], c.sc = t_s[i];
              return c;
            });
          }
          // the misses, on their own line
          {
            const int m = f_j.size();
            work(m);
            for (int i = 0; i < m; ++i) {
              const Cand &c = Qd[f_j[i]];
              w_z0[i] = c.z0, w_cot[i] = c.cot, w_allow[i] = 0, w_holes[i] = c.holes + (c.st == 2);
            }
            route(p, m);
            push(Q2_, m, [&](int i) {
              Cand c = Qd[f_j[i]];
              c.holes = w_holes[i];
              return c;
            });
          }
          if (phases)
            cyc_c += __builtin_ia32_rdtsc() - tc0;
        }
        Qd.clear();

        // -- stage d: the triplets queued at p
        auto &Qt = Q3_[p];
        for (size_t b0 = 0; b0 < Qt.size(); b0 += kBlk) {
          const unsigned long long td0 = phases ? __builtin_ia32_rdtsc() : 0;
          const size_t b1 = std::min(Qt.size(), b0 + kBlk);
          // the helices of the block
          for (size_t j = b0; j < b1; ++j) {
            const Cand &c = Qt[j];
            const SurfLayer &La = *lay[c.pos[0]], &Lb = *lay[c.pos[1]], &Lc = *lay[c.pos[2]];
            const unsigned ka = c.k[0], kb = c.k[1], kc = c.k[2];
            hx[j - b0].make(La.sl.x_[ka], La.sl.y_[ka], La.sl.z_[ka], Lb.sl.x_[kb], Lb.sl.y_[kb], Lb.sl.z_[kb],
                            Lc.sl.x_[kc], Lc.sl.y_[kc], Lc.sl.z_[kc]);
          }
          f_j.clear();
          for (size_t j = b0; j < b1; ++j) {
            const Cand &c = Qt[j];
            const int ix = idx_d(c, p);
            bool found = false;
            const HelixF &H = hx[j - b0];
            if ((!Ch.known_only || ix >= 0) && H.ok) {
              ++n_d_cand;
              const ParF &w = par_[ix >= 0 ? ix : 0];
              float aphi = w.aphi, bphi = w.bphi, aq = w.aq, bq = w.bq;
              if (w.ne > 0) {
                const float ae = std::asinh(std::abs(H.cot));
                for (int i = 0; i < w.ne; ++i) {
                  const SurfParams::EtaWin &E = etaw_[w.e0 + i];
                  if (ae >= E.lo && ae < E.hi) {
                    aphi = E.aphi, bphi = E.bphi, aq = E.aq, bq = E.bq;
                    break;
                  }
                }
              }
              Hermite2F H2;
              if (d_mode == 2)
                H2.make(H, disc, u0, u1);
              // the prediction at qbar u in the chosen mode; the two-point mode falls back to the direct
              // one where its span is degenerate (an edge not reached)
              auto pred = [&](float u, float &px, float &py, float &q, float &s3) {
                return d_mode == 1 ? H.at1(disc, u, px, py, q, s3)
                       : d_mode == 2 && H2.ok ? H2.at(disc, u, px, py, q, s3)
                                             : H.at(disc, u, px, py, q, s3);
              };
              float x0, y0, q0, s0 = -1, x1, y1, q1, s1 = -1;
              const bool ok0 = pred(u0, x0, y0, q0, s0), ok1 = pred(u1, x1, y1, q1, s1);
              if (ok0 || ok1) {
                if (!ok0)
                  x0 = x1, y0 = y1, q0 = q1;
                if (!ok1)
                  x1 = x0, y1 = y0, q1 = q0;
                const float p0 = std::atan2(y0, x0), p1 = std::atan2(y1, x1);
                const float pte = H.pt(), ipt = 1.0f / std::max(pt_min, pte);
                const float isr = w.sref > 0 ? 1.0f / w.sref : 0.0f;
                const float smax = std::max(s0, s1), lmax = isr > 0 && smax > 0 ? smax * isr : 1.0f;
                const float wpd = aphi + lmax * bphi * ipt, wqd = aq + lmax * bq * ipt;
                const auto qd = surf::q_bins(*T, std::min(q0, q1) - wqd, std::max(q0, q1) + wqd);
                const float dpp = wrap(p1 - p0), dmid = p0 + 0.5f * dpp, dhalf = 0.5f * std::abs(dpp);
                const float bqi = bq * ipt, bpi = bphi * ipt;
                const float sw = std::sin(wpd), sw2 = sw * sw;
                // the rest of the score, and the shape band of hit d on the a-c line
                const float fkd = fks - c.sc;
                int slo = 0, shi = 1 << 30;
                if (shT) {
                  const int sb = shT->bin(std::abs(c.cot));
                  slo = shT->lo[sb], shi = shT->hi[sb];
                }
                T->sl.for_each_run(surf::phi_bins(*T, dmid, dhalf + wpd + 1e-6f), qd, [&](unsigned int b, unsigned int e) {
                  for (unsigned int i0 = b; i0 < e; i0 += 64) {
                    const unsigned int nk = std::min(64u, e - i0);
                    unsigned char msk[64];
                    float dqv[64], wqv[64], c2v[64], snv[64];
                    for (unsigned int jj = 0; jj < nk; ++jj) {
                      const unsigned int i = i0 + jj;
                      float px, py, qp, s3 = -1;
                      bool ok = pred(hu[i], px, py, qp, s3);
                      float wq = wqd, sp2 = sw2;
                      if (isr > 0) {
                        const float lv = s3 > 0 ? s3 * isr : 1.0f, wp = aphi + lv * bpi, sp = std::sin(wp);
                        wq = aq + lv * bqi, sp2 = sp * sp;
                      }
                      const float hxx = hx_[i], hyy = hy_[i];
                      const float cr = px * hyy - py * hxx, dt = px * hxx + py * hyy;
                      const float n2 = (px * px + py * py) * (hxx * hxx + hyy * hyy);
                      const float dq = hq[i] - qp, c2 = cr * cr, sn = sp2 * n2;
                      msk[jj] = ok & (std::abs(dq) <= wq) & (dt > 0) & (c2 <= sn);
                      // for the score of the survivors: (dq / wq)^2 + sin^2(dphi) / sin^2(wp)
                      dqv[jj] = dq, wqv[jj] = wq, c2v[jj] = c2, snv[jj] = sn;
                      if (CK.on)
                        check_d(Qt[j], lay, disc, hu[i], ok, px, py, qp, wq, isr > 0 ? std::asin(std::sqrt(sp2)) : wpd, CK);
                    }
                    n_dt += nk;
                    for (unsigned int jj = 0; jj < nk; ++jj)
                      if (msk[jj]) {
                        if (fk_on) {
                          const float rq = dqv[jj] / wqv[jj];
                          if (rq * rq + c2v[jj] / snv[jj] >= fkd || hsp[i0 + jj] < slo || hsp[i0 + jj] > shi)
                            continue;
                        }
                        if (lay_ot2 && Ch.order[p] == 4) {
                          // OT2-P for the helix through b, c, d
                          const SurfLayer &Lb = *lay[c.pos[1]], &Lc = *lay[c.pos[2]];
                          const unsigned kb = c.k[1], kc = c.k[2], kd = i0 + jj;
                          HelixF H3;
                          H3.make(Lb.sl.x_[kb], Lb.sl.y_[kb], Lb.sl.z_[kb], Lc.sl.x_[kc], Lc.sl.y_[kc], Lc.sl.z_[kc],
                                  hx_[kd], hy_[kd], hz[kd]);
                          float zm = 0, dpb = 0, dzb = 0, scb = 0;
                          const float ip2 = 1.0f / std::max(0.9f, pte);
                          const float wp2 = fk_ot2 * (ot2_aphi + ot2_bphi * ip2), wz2 = fk_ot2 * (ot2_aq + ot2_bq * ip2);
                          const int kn = next_hit(*lay_ot2, H3, pte, wp2, wz2, zm, dpb, dzb, scb);
                          if (kn != -2 && std::abs(zm) < lay_ot2->q_hi - 2) {
                            ++n_ot2_tested;
                            const bool pass = kn >= 0 && std::abs(dpb) < wp2 && std::abs(dzb) < wz2;
                            if (!pass) {
                              ++n_fk_ot2;
                              continue;
                            }
                          }
                        }
                        found = true;
                        ++n_qd;
                        const std::array<int, 4> ids{Ch.order[c.pos[0]], Ch.order[c.pos[1]], Ch.order[c.pos[2]], Ch.order[p]};
                        out.push_back({ids,
                                       {lay[c.pos[0]]->sl.orig_[c.k[0]], lay[c.pos[1]]->sl.orig_[c.k[1]],
                                        lay[c.pos[2]]->sl.orig_[c.k[2]], T->sl.orig_[i0 + jj]}});
                      }
                  }
                });
              }
            }
            if (!found || Ch.hole_always) {
              if (c.holes + (c.st == 2) <= max_holes) {
                ++n_forwarded;
                f_j.push_back(j);
              } else
                ++n_dropped;
            }
          }
          {
            const int m = f_j.size();
            work(m);
            for (int i = 0; i < m; ++i) {
              const Cand &c = Qt[f_j[i]];
              const int h = c.holes + (c.st == 2);
              w_z0[i] = c.z0, w_cot[i] = c.cot, w_holes[i] = h, w_allow[i] = h <= max_holes_ot;
            }
            route(p, m);
            push(Q3_, m, [&](int i) {
              Cand c = Qt[f_j[i]];
              c.holes = w_holes[i];
              return c;
            });
          }
          if (phases)
            cyc_d += __builtin_ia32_rdtsc() - td0;
        }
        Qt.clear();
      }
      cnt.c_touched += n_ct, cnt.triplets += n_tr, cnt.d_touched += n_dt, cnt.quads += n_qd;
    }

    // --chain-fast-check: the float single-point prediction against surf::Helix::predict in double
    static void check_d(const Cand &c, const std::vector<const SurfLayer *> &lay, bool disc, float u, bool ok, float px,
                        float py, float qp, float wq, float wp, SurfFastCheck &CK) {
      const surf::Helix H(surf_p3(*lay[c.pos[0]], c.k[0]), surf_p3(*lay[c.pos[1]], c.k[1]), surf_p3(*lay[c.pos[2]], c.k[2]));
      double qe, pe;
      const bool oke = H.ok && H.predict(disc, u, qe, pe);
      if (oke != ok) {
        ++CK.b_fail_mismatch;
        return;
      }
      if (!ok)
        return;
      ++CK.b_cands;
      const double dq = std::abs(qp - qe), dp = std::abs(surf::wrap(std::atan2((double)py, (double)px) - pe));
      CK.b_max_dq = std::max(CK.b_max_dq, dq), CK.b_max_dphi = std::max(CK.b_max_dphi, dp);
      CK.b_max_rel_q = std::max(CK.b_max_rel_q, dq / wq), CK.b_max_rel_phi = std::max(CK.b_max_rel_phi, dp / wp);
    }
  };

}  // namespace mkfit::seeding

#endif
