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
#include "SeedMath.h"

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
    // stage 5, curvature form: the triplet's hits as struct-of-arrays
    std::vector<float> g_xa, g_ya, g_za, g_xb, g_yb, g_xc, g_yc, g_zc;
    // stage 7, curvature form: each (triplet, d-hit) candidate as struct-of-arrays
    std::vector<float> q_k, q_ux, q_uy, q_xc, q_yc, q_cs, q_zc, q_rc, q_rd, q_pd, q_zd;
    std::vector<unsigned char> q_ok;
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
  // holding the doublets they index.  The arithmetic is A's (SeedMath.h), and
  // it is the same code the margin evaluator runs.
  template <class A, typename LA, typename LB, typename LC, typename LD, typename Tick>
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
    using T = typename A::real;
    const T rdlo = gd.rlo_, rdhi = gd.rhi_;
    auto phi_bins = [](const auto &L, float c, float w) {
      w = std::min(w, 0.9f * kPi);
      return L.phi_range(c - w, c + w);
    };
    //---- 5: helix per triplet, and the 4th-layer window
    for (auto *v : {&W.t_cx, &W.t_cy, &W.t_R, &W.t_cots, &W.t_xc, &W.t_yc, &W.t_phid, &W.t_phidh, &W.t_zlo, &W.t_zhi})
      StagedWork::fit(*v, nt);
    StagedWork::fit(W.t_ok, nt);
    if constexpr (A::curvature_form) {
      // gather the three hits into struct-of-arrays, then the helix across
      // triplets in one simd loop
      for (auto *v : {&W.g_xa, &W.g_ya, &W.g_za, &W.g_xb, &W.g_yb, &W.g_xc, &W.g_yc, &W.g_zc})
        StagedWork::fit(*v, nt);
      for (unsigned int t = 0; t < nt; ++t) {
        const unsigned int d = W.t_d[t], kc = W.t_kc[t];
        const unsigned int ka = W.d_ka[d], kb = W.d_kb[d];
        W.g_xa[t] = ga.x_[ka];
        W.g_ya[t] = ga.y_[ka];
        W.g_za[t] = ga.z_[ka];
        W.g_xb[t] = gb.x_[kb];
        W.g_yb[t] = gb.y_[kb];
        W.g_xc[t] = gc.x_[kc];
        W.g_yc[t] = gc.y_[kc];
        W.g_zc[t] = gc.z_[kc];
      }
      const float *__restrict xa = W.g_xa.data(), *__restrict ya = W.g_ya.data(), *__restrict za = W.g_za.data();
      const float *__restrict xb = W.g_xb.data(), *__restrict yb = W.g_yb.data();
      const float *__restrict xc = W.g_xc.data(), *__restrict yc = W.g_yc.data(), *__restrict zc = W.g_zc.data();
      double *__restrict ocx = W.t_cx.data(), *__restrict ocy = W.t_cy.data(), *__restrict oR = W.t_R.data();
      double *__restrict ocs = W.t_cots.data(), *__restrict oxc = W.t_xc.data(), *__restrict oyc = W.t_yc.data();
      double *__restrict opd = W.t_phid.data(), *__restrict oph = W.t_phidh.data();
      double *__restrict ozl = W.t_zlo.data(), *__restrict ozh = W.t_zhi.data();
      unsigned char *__restrict ook = W.t_ok.data();
      const T rlo = rdlo, rhi = rdhi;
      const float qd = P.qwin_d, pd = P.phiwin_d;
#pragma omp simd
      for (unsigned int t = 0; t < nt; ++t) {
        HelixK<A> o;
        helix_k<A>(qd, pd, xa[t], ya[t], za[t], xb[t], yb[t], xc[t], yc[t], zc[t], rlo, rhi, o);
        ook[t] = o.ok;
        ocx[t] = o.ux;  // the curvature form stores (k, ux, uy) where the centre form has (R, cx, cy)
        ocy[t] = o.uy;
        oR[t] = o.k;
        ocs[t] = o.cots;
        oxc[t] = xc[t];
        oyc[t] = yc[t];
        opd[t] = o.phid;
        oph[t] = o.phidh;
        ozl[t] = o.zlo;
        ozh[t] = o.zhi;
      }
    } else
    for (unsigned int t = 0; t < nt; ++t) {
      const unsigned int d = W.t_d[t], kc = W.t_kc[t];
      const unsigned int ka = W.d_ka[d], kb = W.d_kb[d];
      TripletHelix<A> h;
      // x, y precomputed at fill with the same float expression, r * cos(phi)
      helix_triplet<A>(P.qwin_d, P.phiwin_d, ga.x_[ka], ga.y_[ka], ga.z_[ka], gb.x_[kb], gb.y_[kb], gc.x_[kc],
                       gc.y_[kc], gc.z_[kc], rdlo, rdhi, h);
      W.t_ok[t] = h.ok;
      if (!h.ok)
        continue;
      // the curvature form stores (k, ux, uy) where the centre form has (R, cx, cy)
      W.t_cx[t] = A::curvature_form ? h.ux : h.cx;
      W.t_cy[t] = A::curvature_form ? h.uy : h.cy;
      W.t_R[t] = A::curvature_form ? h.k : h.R;
      W.t_cots[t] = h.cots;
      W.t_xc[t] = h.xc;
      W.t_yc[t] = h.yc;
      W.t_phid[t] = h.phid;
      W.t_phidh[t] = h.phidh;
      W.t_zlo[t] = h.zlo;
      W.t_zhi[t] = h.zhi;
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
    if constexpr (A::curvature_form) {
      // gather each candidate's triplet and d-hit into struct-of-arrays, the
      // test in one simd loop, then compact
      for (auto *v : {&W.q_k, &W.q_ux, &W.q_uy, &W.q_xc, &W.q_yc, &W.q_cs, &W.q_zc, &W.q_rc, &W.q_rd, &W.q_pd, &W.q_zd})
        StagedWork::fit(*v, ne);
      StagedWork::fit(W.q_ok, ne);
      for (unsigned int i = 0; i < ne; ++i) {
        const unsigned int t = W.e_t[i], kd = W.e_kd[i], kc = W.t_kc[t];
        W.q_k[i] = W.t_R[t];
        W.q_ux[i] = W.t_cx[t];
        W.q_uy[i] = W.t_cy[t];
        W.q_xc[i] = W.t_xc[t];
        W.q_yc[i] = W.t_yc[t];
        W.q_cs[i] = W.t_cots[t];
        W.q_zc[i] = gc.z_[kc];
        W.q_rc[i] = gc.r_[kc];
        W.q_rd[i] = gd.r_[kd];
        W.q_pd[i] = gd.phi_[kd];
        W.q_zd[i] = gd.z_[kd];
      }
      const float *__restrict qk = W.q_k.data(), *__restrict qux = W.q_ux.data(), *__restrict quy = W.q_uy.data();
      const float *__restrict qxc = W.q_xc.data(), *__restrict qyc = W.q_yc.data(), *__restrict qcs = W.q_cs.data();
      const float *__restrict qzc = W.q_zc.data(), *__restrict qrc = W.q_rc.data(), *__restrict qrd = W.q_rd.data();
      const float *__restrict qpd = W.q_pd.data(), *__restrict qzd = W.q_zd.data();
      unsigned char *__restrict qok = W.q_ok.data();
      const float qd = P.qwin_d, pd = P.phiwin_d;
#pragma omp simd
      for (unsigned int i = 0; i < ne; ++i) {
        QuadK<A> o;
        quad_k<A>(qd, pd, qk[i], qux[i], quy[i], qxc[i], qyc[i], qcs[i], qzc[i], qrc[i], qrd[i], qpd[i], qzd[i], o);
        qok[i] = o.pass;
      }
      for (unsigned int i = 0; i < ne; ++i) {
        if (!qok[i])
          continue;
        const unsigned int t = W.e_t[i], kd = W.e_kd[i];
        cnt.quads++;
        const unsigned int d = W.t_d[t];
        out.push_back({ga.orig_[W.d_ka[d]], gb.orig_[W.d_kb[d]], gc.orig_[W.t_kc[t]], gd.orig_[kd]});
      }
    } else
    for (unsigned int i = 0; i < ne; ++i) {
      const unsigned int t = W.e_t[i], kd = W.e_kd[i];
      const unsigned int kc = W.t_kc[t];
      TripletHelix<A> h;
      if constexpr (A::curvature_form) {
        h.ux = (T)W.t_cx[t];
        h.uy = (T)W.t_cy[t];
        h.k = (T)W.t_R[t];
      } else {
        h.cx = (T)W.t_cx[t];
        h.cy = (T)W.t_cy[t];
        h.R = (T)W.t_R[t];
      }
      h.cots = (T)W.t_cots[t];
      h.xc = (T)W.t_xc[t];
      h.yc = (T)W.t_yc[t];
      if (!quad_cuts<A>(P.qwin_d, P.phiwin_d, h, gc.z_[kc], gc.r_[kc], gd.r_[kd], gd.phi_[kd], gd.z_[kd], nullptr))
        continue;
      cnt.quads++;
      const unsigned int d = W.t_d[t];
      out.push_back({ga.orig_[W.d_ka[d]], gb.orig_[W.d_kb[d]], gc.orig_[kc], gd.orig_[kd]});
    }
    tick(6);
  }

  template <class A = ArithRef, typename LA, typename LB, typename LC, typename LD>
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
      finish_triplets<A>(P, ga, gb, gc, gd, nt, W, out, cnt, tick);
    }
  }

}  // namespace mkfit::seeding

#endif
