#ifndef mkfit_seeding_SeedFinderStaged_h
#define mkfit_seeding_SeedFinderStaged_h

// The same search as find_quads(), restructured as a pipeline of flat loops.
//
// find_quads() walks a deep call tree per doublet: form the layer-c window,
// walk its bin runs, and test each hit behind a chain of data-dependent
// branches.  Measured on one event: ~600 instructions and 5.5 branch
// mispredictions per doublet, with the data all in L2.  Here each stage is one
// loop over a flat array that the previous stage filled, the tests append
// branch-free (`out[n] = x; n += pass;`), and the arithmetic of every stage
// sits in a loop without calls, so it can be vectorised ACROSS doublets or
// candidates rather than within one.
//
//   1  doublets     a-hit x layer-b phi run          -> (ka, kb)
//   2  windows      per doublet, branch-free         -> layer-c phi / z window
//   3  candidates   per doublet, the bin runs        -> (doublet, kc)
//   4  triplets     per candidate, branch-free       -> (doublet, kc)
//   5  helix        per triplet: circle, 4th window  -> (cx, cy, R, cot_s, ...)
//   6  candidates   per triplet, the bin runs        -> (triplet, kd)
//   7  quads        per candidate                    -> Quad
//
// The a-hits go through in BLOCKS so the per-stage arrays stay in L2.
//
// The per-candidate arithmetic is find_quads()'s, expression for expression,
// so the quad list must equal it exactly.

#include "SeedFinder.h"

namespace mkfit::seeding {

  struct StagedWork {
    // stage 1: doublets
    std::vector<unsigned int> d_ka, d_kb;
    // stage 2: per-doublet quantities and the layer-c window
    std::vector<double> d_cot;
    std::vector<float> d_slope, d_kab, d_dinv, d_phic, d_wbc, d_zlo, d_zhi;
    // stage 3/4: (doublet, c-hit) candidates, then triplets
    std::vector<unsigned int> c_d, c_kc, t_d, t_kc;
    // stage 5: per triplet
    std::vector<double> t_cx, t_cy, t_R, t_cots, t_xc, t_yc, t_phid, t_phidh, t_zlo, t_zhi;
    std::vector<unsigned char> t_ok;
    // stage 6: (triplet, d-hit) candidates
    std::vector<unsigned int> e_t, e_kd;

    template <typename V>
    static void fit(V &v, size_t n) {
      if (v.size() < n)
        v.resize(std::max(n, 2 * v.size()));
    }
  };

  // Stages 5-7, shared by every variant: the helix per triplet and the 4th
  // layer.  The triplets are (W.t_d, W.t_kc)[0, nt), with W.d_ka / W.d_kb
  // holding the doublets they index.
  template <typename LA, typename LB, typename LC, typename LD, typename Tick>
  void finish_triplets(const SeedParams &P,
                       const LA &ga,
                       const LB &gb,
                       const LC &gc,
                       const LD &gd,
                       unsigned int nt,
                       StagedWork &W,
                       std::vector<Quad> &out,
                       SeedCounters &cnt,
                       Tick &&tick) {
    using namespace detail;
    const double rdlo = gd.rlo_, rdhi = gd.rhi_;
    auto phi_bins = [](const auto &L, float c, float w) {
      w = std::min(w, 0.9f * kPi);
      return L.phi_range(c - w, c + w);
    };
    //---- 5: helix per triplet, and the 4th-layer window
    for (auto *v : {&W.t_cx, &W.t_cy, &W.t_R, &W.t_cots, &W.t_xc, &W.t_yc, &W.t_phid, &W.t_phidh, &W.t_zlo, &W.t_zhi})
      StagedWork::fit(*v, nt);
    StagedWork::fit(W.t_ok, nt);
    for (unsigned int t = 0; t < nt; ++t) {
      const unsigned int d = W.t_d[t], kc = W.t_kc[t];
      const unsigned int ka = W.d_ka[d], kb = W.d_kb[d];
      const float za = ga.z_[ka];
      // precomputed at fill with the same float expression, r * cos(phi)
      const double xa = ga.x_[ka], ya = ga.y_[ka];
      const double xb = gb.x_[kb], yb = gb.y_[kb];
      const double xc = gc.x_[kc], yc = gc.y_[kc];
      double cx = 0, cy = 0, R;
      if (!circle3(xa, ya, xb, yb, xc, yc, cx, cy, R))
        R = 1e6;
      const double s_ab = arc(std::hypot(xb - xa, yb - ya), R);
      const double s_bc = arc(std::hypot(xc - xb, yc - yb), R);
      const double s_ac = s_ab + s_bc;
      const double cot_s = (s_ac > 1e-6) ? (gc.z_[kc] - za) / s_ac : 0.0;
      double p0x, p0y, p1x, p1y;
      const bool ok0 = circle_cross_r(cx, cy, R, rdlo, xc, yc, p0x, p0y);
      const bool ok1 = circle_cross_r(cx, cy, R, rdhi, xc, yc, p1x, p1y);
      W.t_ok[t] = ok0 || ok1;
      if (!(ok0 || ok1))
        continue;
      const double f0 = ok0 ? std::atan2(p0y, p0x) : std::atan2(p1y, p1x);
      const double f1 = ok1 ? std::atan2(p1y, p1x) : f0;
      const double dfh = 0.5 * wrap_pi((float)(f1 - f0));
      const double zz0 = ok0 ? gc.z_[kc] + cot_s * arc(std::hypot(p0x - xc, p0y - yc), R)
                             : gc.z_[kc] + cot_s * arc(std::hypot(p1x - xc, p1y - yc), R);
      const double zz1 = ok1 ? gc.z_[kc] + cot_s * arc(std::hypot(p1x - xc, p1y - yc), R) : zz0;
      const double zdh = P.qwin_d + kCoverEps;
      W.t_cx[t] = cx;
      W.t_cy[t] = cy;
      W.t_R[t] = R;
      W.t_cots[t] = cot_s;
      W.t_xc[t] = xc;
      W.t_yc[t] = yc;
      W.t_phid[t] = f0 + dfh;
      W.t_phidh[t] = P.phiwin_d + std::abs(dfh) + kCoverEps;
      W.t_zlo[t] = std::min(zz0, zz1) - zdh;
      W.t_zhi[t] = std::max(zz0, zz1) + zdh;
    }

    tick(4);
    //---- 6: (triplet, d-hit) candidates
    unsigned int ne = 0;
    for (unsigned int t = 0; t < nt; ++t) {
      if (!W.t_ok[t])
        continue;
      gd.for_each_run(phi_bins(gd, (float)W.t_phid[t], (float)W.t_phidh[t]),
                      gd.q_range((float)W.t_zlo[t], (float)W.t_zhi[t]),
                      [&](unsigned int e0, unsigned int e1) {
                        StagedWork::fit(W.e_t, ne + (e1 - e0));
                        StagedWork::fit(W.e_kd, ne + (e1 - e0));
                        for (unsigned int kd = e0; kd < e1; ++kd, ++ne) {
                          W.e_t[ne] = t;
                          W.e_kd[ne] = kd;
                        }
                      });
    }
    cnt.d_touched += ne;

    tick(5);
    //---- 7: quads
    for (unsigned int i = 0; i < ne; ++i) {
      const unsigned int t = W.e_t[i], kd = W.e_kd[i];
      const unsigned int kc = W.t_kc[t];
      const float rrc = gc.r_[kc], rrd = gd.r_[kd];
      if (rrd <= rrc + 0.1f)
        continue;
      const double cx = W.t_cx[t], cy = W.t_cy[t], R = W.t_R[t], xc = W.t_xc[t], yc = W.t_yc[t];
      double px2, py2;
      if (!circle_cross_r(cx, cy, R, rrd, xc, yc, px2, py2))
        continue;
      const double dphi_d = wrap_pi((float)(gd.phi_[kd] - std::atan2(py2, px2)));
      const double zd2 = gc.z_[kc] + W.t_cots[t] * arc(std::hypot(px2 - xc, py2 - yc), R);
      const double dz_d = gd.z_[kd] - zd2;
      if (std::abs(dphi_d) > P.phiwin_d || std::abs(dz_d) > P.qwin_d)
        continue;
      cnt.quads++;
      const unsigned int d = W.t_d[t];
      out.push_back({ga.orig_[W.d_ka[d]], gb.orig_[W.d_kb[d]], gc.orig_[kc], gd.orig_[kd]});
    }
    tick(6);
  }

  template <typename LA, typename LB, typename LC, typename LD>
  void find_quads_staged(const SeedParams &P,
                         const LA &ga,
                         const LB &gb,
                         const LC &gc,
                         const LD &gd,
                         std::vector<Quad> &out,
                         SeedCounters &cnt,
                         StagedWork &W,
                         unsigned int block = 64,
                         bool fuse = false) {
    using namespace detail;
    constexpr float kBfield = 3.8f;
    const double Rmin = P.pt_min / (0.003f * kBfield);
    const float inv2R = (float)(1.0 / (2 * Rmin)), d0m = P.d0_max, marg = P.phi_margin;
    const double rbhi = gb.rhi_, rclo = gc.rlo_, rchi = gc.rhi_, rdlo = gd.rlo_, rdhi = gd.rhi_;
    const float inv_rbhi = 1.0f / (float)rbhi, inv_rchi = 1.0f / (float)rchi;
    const float mlin = P.phi_lin_marg;
    const float rmid = 0.5f * (float)(rclo + rchi);
    const float inv_rclo_f = 1.0f / (float)rclo;
    auto wphi_w = [=](float r_in, float inv_in, float r_out, float inv_out) {
      return (r_out - r_in) * inv2R + d0m * (inv_in - inv_out) + marg;
    };
    auto phi_bins = [](const auto &L, float c, float w) {
      w = std::min(w, 0.9f * kPi);
      return L.phi_range(c - w, c + w);
    };

    using clk = std::chrono::steady_clock;
    for (unsigned int ka0 = 0; ka0 < ga.n(); ka0 += block) {
      const unsigned int ka1 = std::min(ga.n(), ka0 + block);
      auto tprev = clk::now();
      auto tick = [&](int st) {
        const auto now = clk::now();
        cnt.t_stage[st] += std::chrono::duration<double>(now - tprev).count();
        tprev = now;
      };

      //---- 1: doublets
      unsigned int nd = 0;
      for (unsigned int ka = ka0; ka < ka1; ++ka) {
        const float pa = ga.phi_[ka], rra = ga.r_[ka], inva = ga.invr_[ka];
        const float w_ab = wphi_w(rra, inva, (float)rbhi, inv_rbhi);
        gb.for_each_run(phi_bins(gb, pa, w_ab), gb.q_all(), [&](unsigned int b0, unsigned int b1) {
          StagedWork::fit(W.d_ka, nd + (b1 - b0));
          StagedWork::fit(W.d_kb, nd + (b1 - b0));
          unsigned int *__restrict oka = W.d_ka.data();
          unsigned int *__restrict okb = W.d_kb.data();
          for (unsigned int kb = b0; kb < b1; ++kb) {
            const float rrb = gb.r_[kb];
            const bool pass = !(rrb <= rra + 0.1f) &&
                              !(std::abs(wrap_pi(gb.phi_[kb] - pa)) > wphi_w(rra, inva, rrb, gb.invr_[kb]));
            oka[nd] = ka;
            okb[nd] = kb;
            nd += pass;
          }
        });
      }
      cnt.doublets += nd;

      tick(0);
      //---- 2: per-doublet quantities and the layer-c window
      for (auto *v : {&W.d_slope, &W.d_kab, &W.d_dinv, &W.d_phic, &W.d_wbc, &W.d_zlo, &W.d_zhi})
        StagedWork::fit(*v, nd);
      StagedWork::fit(W.d_cot, nd);
      for (unsigned int d = 0; d < nd; ++d) {
        const unsigned int ka = W.d_ka[d], kb = W.d_kb[d];
        const float pa = ga.phi_[ka], za = ga.z_[ka], rra = ga.r_[ka], inva = ga.invr_[ka];
        const float pbph = gb.phi_[kb], zb = gb.z_[kb], rrb = gb.r_[kb], invb = gb.invr_[kb];
        const double cot = (zb - za) / (rrb - rra);
        const float slope = wrap_pi(pbph - pa) / (rrb - rra);
        const float kab = 1.0f / (rrb - rra);
        const float dinv_ab = invb - inva;
        auto wlin = [&](float r_c, float inv_c) {
          return d0m * std::abs(inv_c - invb - (r_c - rrb) * kab * dinv_ab) + mlin;
        };
        float w_bc, phic_ctr;
        if (P.phi_lin) {
          phic_ctr = pbph + slope * (rmid - rrb);
          w_bc = std::max(wlin((float)rclo, inv_rclo_f), wlin((float)rchi, inv_rchi)) +
                 std::abs(slope) * 0.5f * (float)(rchi - rclo);
        } else {
          phic_ctr = pbph;
          w_bc = wphi_w(rrb, invb, (float)rchi, inv_rchi);
        }
        const double z0 = za + cot * (rclo - rra), z1 = za + cot * (rchi - rra);
        const double zh = P.qwin + kCoverEps;
        W.d_cot[d] = cot;
        W.d_slope[d] = slope;
        W.d_kab[d] = kab;
        W.d_dinv[d] = dinv_ab;
        W.d_phic[d] = phic_ctr;
        W.d_wbc[d] = w_bc;
        W.d_zlo[d] = (float)(std::min(z0, z1) - zh);
        W.d_zhi[d] = (float)(std::max(z0, z1) + zh);
      }

      tick(1);
      //---- 3: (doublet, c-hit) candidates
      unsigned int nc = 0, nt = 0;
      if (fuse) {
        // 3+4 FUSED: the doublet's quantities stay in registers and the
        // branch-free test runs straight over each contiguous c-run, so no
        // candidate is materialised and nothing is gathered through indices
        for (unsigned int d = 0; d < nd; ++d) {
          const unsigned int ka = W.d_ka[d], kb = W.d_kb[d];
          const float za = ga.z_[ka], rra = ga.r_[ka];
          const float pbph = gb.phi_[kb], rrb = gb.r_[kb], invb = gb.invr_[kb];
          const double cot = W.d_cot[d];
          const float slope = W.d_slope[d], kab = W.d_kab[d], dinv_ab = W.d_dinv[d];
          gc.for_each_run(phi_bins(gc, W.d_phic[d], W.d_wbc[d]), gc.q_range(W.d_zlo[d], W.d_zhi[d]),
                          [&](unsigned int c0, unsigned int c1) {
                            nc += c1 - c0;
                            StagedWork::fit(W.t_d, nt + (c1 - c0));
                            StagedWork::fit(W.t_kc, nt + (c1 - c0));
                            unsigned int *__restrict otd = W.t_d.data();
                            unsigned int *__restrict otc = W.t_kc.data();
                            const float *__restrict cr = gc.r_.data();
                            const float *__restrict cp = gc.phi_.data();
                            const float *__restrict ci = gc.invr_.data();
                            const float *__restrict cz = gc.z_.data();
                            for (unsigned int kc = c0; kc < c1; ++kc) {
                              const float rrc = cr[kc], pcph = cp[kc], invc = ci[kc];
                              bool pass = !(rrc <= rrb + 0.1f);
                              if (P.phi_lin) {
                                const float wl = d0m * std::abs(invc - invb - (rrc - rrb) * kab * dinv_ab) + mlin;
                                pass &= !(std::abs(wrap_pi(pcph - pbph - slope * (rrc - rrb))) > wl);
                              }
                              if (P.phi_lin != 1)
                                pass &= !(std::abs(wrap_pi(pcph - pbph)) > wphi_w(rrb, invb, rrc, invc));
                              const double zpred = za + cot * (rrc - rra);
                              const double dz = cz[kc] - zpred;
                              pass &= !(std::abs(dz) > P.qwin);
                              otd[nt] = d;
                              otc[nt] = kc;
                              nt += pass;
                            }
                          });
        }
      } else
      for (unsigned int d = 0; d < nd; ++d) {
        gc.for_each_run(phi_bins(gc, W.d_phic[d], W.d_wbc[d]),
                        gc.q_range(W.d_zlo[d], W.d_zhi[d]),
                        [&](unsigned int c0, unsigned int c1) {
                          StagedWork::fit(W.c_d, nc + (c1 - c0));
                          StagedWork::fit(W.c_kc, nc + (c1 - c0));
                          for (unsigned int kc = c0; kc < c1; ++kc, ++nc) {
                            W.c_d[nc] = d;
                            W.c_kc[nc] = kc;
                          }
                        });
      }
      cnt.c_touched += nc;

      tick(2);
      //---- 4: triplets
      if (!fuse) {
      StagedWork::fit(W.t_d, nc);
      StagedWork::fit(W.t_kc, nc);
        unsigned int *__restrict otd = W.t_d.data();
        unsigned int *__restrict otc = W.t_kc.data();
        for (unsigned int i = 0; i < nc; ++i) {
          const unsigned int d = W.c_d[i], kc = W.c_kc[i];
          const unsigned int ka = W.d_ka[d], kb = W.d_kb[d];
          const float za = ga.z_[ka], rra = ga.r_[ka];
          const float pbph = gb.phi_[kb], rrb = gb.r_[kb], invb = gb.invr_[kb];
          const float rrc = gc.r_[kc], pcph = gc.phi_[kc], invc = gc.invr_[kc];
          bool pass = !(rrc <= rrb + 0.1f);
          if (P.phi_lin) {
            const float slope = W.d_slope[d], kab = W.d_kab[d], dinv_ab = W.d_dinv[d];
            const float wl = d0m * std::abs(invc - invb - (rrc - rrb) * kab * dinv_ab) + mlin;
            pass &= !(std::abs(wrap_pi(pcph - pbph - slope * (rrc - rrb))) > wl);
          }
          if (P.phi_lin != 1)
            pass &= !(std::abs(wrap_pi(pcph - pbph)) > wphi_w(rrb, invb, rrc, invc));
          const double zpred = za + W.d_cot[d] * (rrc - rra);
          const double dz = gc.z_[kc] - zpred;
          pass &= !(std::abs(dz) > P.qwin);
          otd[nt] = d;
          otc[nt] = kc;
          nt += pass;
        }
      }
      cnt.triplets += nt;

      tick(3);
      finish_triplets(P, ga, gb, gc, gd, nt, W, out, cnt, tick);
    }
  }

}  // namespace mkfit::seeding

#endif
