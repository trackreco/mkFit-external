#ifndef mkfit_seeding_SeedFinderBMajor_h
#define mkfit_seeding_SeedFinderBMajor_h

// The same search with the b-hit as the OUTER loop.
//
// The layer-c phi test of the generic window depends on the b-hit alone: it is
// centred on phi_b, with a half-width set by r_b and the c-hit's own r.  So the
// c-hits that pass it, and the r ordering test, are the same for every doublet
// through that b-hit.  They are collected once per b-hit into a list sorted in
// z, and each doublet then needs only a binary search for its z window plus the
// exact z test on the few hits inside it.  In find_quads_staged() the same work
// was redone per doublet, which costs ~70 ns in bin-run bookkeeping for ~4 hits.
//
// The doublets are found from the b side.  The doublet half-width
// (r_b - r_a)/2R + D0 (1/r_a - 1/r_b) + margin falls with r_a, so the widest a
// window of a b-hit is at layer a's SMALLEST radius, and fetching at that width
// covers the per-pair cut.
//
// phi_lin mode 2 (linear AND generic) works unchanged: the generic test is
// applied when the list is built, the linear one per candidate.  Mode 1 drops
// the generic test and is not supported here.
//
// The cuts are find_quads()'s, expression for expression, so the quad list
// must equal it exactly.

#include "SeedFinderStaged.h"

#include <algorithm>

namespace mkfit::seeding {

  struct BMajorWork {
    std::vector<unsigned int> b_dbeg;  // per b-hit of the block: first doublet
    std::vector<float> l_z;            // per b-hit: c-hits passing the phi and r tests, z-sorted
    std::vector<unsigned int> l_kc;
    std::vector<std::pair<float, unsigned int>> l_tmp;
    std::vector<unsigned int> bk_start;  // z buckets of the list, CSR
    std::vector<unsigned int> l_tmp2;    // fill cursors
  };

  namespace detail {
    // Branch-free lower / upper bound: the comparison feeds an index update
    // (a conditional move), not a jump, so it cannot mispredict.  First index
    // with v[i] >= x, resp. v[i] > x, in [0, n).
    inline unsigned int lower_bound_bf(const float *v, unsigned int n, float x) {
      const float *base = v;
      while (n > 1) {
        const unsigned int half = n / 2;
        base = (base[half - 1] < x) ? base + half : base;
        n -= half;
      }
      return (base - v) + (n == 1 && *base < x);
    }
    inline unsigned int upper_bound_bf(const float *v, unsigned int n, float x) {
      const float *base = v;
      while (n > 1) {
        const unsigned int half = n / 2;
        base = (base[half - 1] <= x) ? base + half : base;
        n -= half;
      }
      return (base - v) + (n == 1 && *base <= x);
    }
  }  // namespace detail

  template <class A = ArithRef, typename LA, typename LB, typename LC, typename LD>
  void find_quads_bmajor(const SeedParams &P,
                         const LA &ga,
                         const LB &gb,
                         const LC &gc,
                         const LD &gd,
                         std::vector<Quad> &out,
                         SeedCounters &cnt,
                         StagedWork &W,
                         BMajorWork &BW,
                         unsigned int block = 32) {
    using namespace detail;
    if (P.phi_lin == 1) {
      fprintf(stderr, "find_quads_bmajor: phi_lin mode 1 not supported\n");
      return;
    }
    constexpr float kBfield = 3.8f;
    const double Rmin = P.pt_min / (0.003f * kBfield);
    const float inv2R = (float)(1.0 / (2 * Rmin)), d0m = P.d0_max, marg = P.phi_margin;
    const double ralo = ga.rlo_, rclo = gc.rlo_, rchi = gc.rhi_;
    const float inv_ralo = 1.0f / (float)ralo, inv_rchi = 1.0f / (float)rchi;
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
    constexpr float kFetchEps = 1e-6f;  // rad, float rounding of the window width
    const float bz0 = gc.q_min(), binv = 1.0f / P.lbin;
    const unsigned int nbk = std::max(1u, (unsigned int)std::ceil((gc.q_max() - bz0) * binv));
    BW.bk_start.resize(nbk + 1);
    auto bucket = [=](float z) -> unsigned int {
      const int b = (int)((z - bz0) * binv);
      return b < 0 ? 0u : (b >= (int)nbk ? nbk - 1 : (unsigned int)b);
    };

    using clk = std::chrono::steady_clock;
    for (unsigned int kb0 = 0; kb0 < gb.n(); kb0 += block) {
      const unsigned int kb1 = std::min(gb.n(), kb0 + block);
      auto tprev = clk::now();
      auto tick = [&](int st) {
        const auto now = clk::now();
        cnt.t_stage[st] += std::chrono::duration<double>(now - tprev).count();
        tprev = now;
      };

      //---- 1: doublets, from the b side
      unsigned int nd = 0;
      BW.b_dbeg.resize(kb1 - kb0 + 1);
      for (unsigned int kb = kb0; kb < kb1; ++kb) {
        BW.b_dbeg[kb - kb0] = nd;
        const float pbph = gb.phi_[kb], rrb = gb.r_[kb], invb = gb.invr_[kb];
        const float w = wphi_w((float)ralo, inv_ralo, rrb, invb) + kFetchEps;
        ga.for_each_run(phi_bins(ga, pbph, w), ga.q_all(), [&](unsigned int a0, unsigned int a1) {
          StagedWork::fit(W.d_ka, nd + (a1 - a0));
          StagedWork::fit(W.d_kb, nd + (a1 - a0));
          unsigned int *__restrict oka = W.d_ka.data();
          unsigned int *__restrict okb = W.d_kb.data();
          for (unsigned int ka = a0; ka < a1; ++ka) {
            const float rra = ga.r_[ka];
            const bool pass = !(rrb <= rra + 0.1f) &&
                              !(std::abs(wrap_pi(pbph - ga.phi_[ka])) > wphi_w(rra, ga.invr_[ka], rrb, invb));
            oka[nd] = ka;
            okb[nd] = kb;
            nd += pass;
          }
        });
      }
      BW.b_dbeg[kb1 - kb0] = nd;
      cnt.doublets += nd;

      tick(0);
      //---- 2: per-doublet quantities and the layer-c z window
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
        const double z0 = za + cot * (rclo - rra), z1 = za + cot * (rchi - rra);
        const double zh = P.qwin + kCoverEps;
        W.d_cot[d] = cot;
        W.d_slope[d] = slope;
        W.d_kab[d] = kab;
        W.d_dinv[d] = dinv_ab;
        W.d_zlo[d] = (float)(std::min(z0, z1) - zh);
        W.d_zhi[d] = (float)(std::max(z0, z1) + zh);
      }

      tick(1);
      //---- 3+4: per b-hit, the phi- and r-passing c-hits sorted in z; per
      //          doublet, a binary search and the exact z test
      unsigned int nt = 0;
      for (unsigned int kb = kb0; kb < kb1; ++kb) {
        const unsigned int db = BW.b_dbeg[kb - kb0], de = BW.b_dbeg[kb - kb0 + 1];
        if (db == de)
          continue;
        const float pbph = gb.phi_[kb], rrb = gb.r_[kb], invb = gb.invr_[kb];
        const float w_bc = wphi_w(rrb, invb, (float)rchi, inv_rchi);
        BW.l_tmp.clear();
        gc.for_each_run(phi_bins(gc, pbph, w_bc), gc.q_all(), [&](unsigned int c0, unsigned int c1) {
          for (unsigned int kc = c0; kc < c1; ++kc) {
            const float rrc = gc.r_[kc];
            if (rrc <= rrb + 0.1f)
              continue;
            if (std::abs(wrap_pi(gc.phi_[kc] - pbph)) > wphi_w(rrb, invb, rrc, gc.invr_[kc]))
              continue;
            BW.l_tmp.emplace_back(gc.z_[kc], kc);
          }
        });
        // counting sort into fixed z buckets: a doublet's window then maps to a
        // bucket range in O(1).  The range is a slight superset of the window,
        // which is harmless since every candidate gets the exact z test.
        const unsigned int nl = BW.l_tmp.size();
        std::fill(BW.bk_start.begin(), BW.bk_start.end(), 0u);
        for (const auto &e : BW.l_tmp)
          ++BW.bk_start[bucket(e.first) + 1];
        for (unsigned int k = 0; k < nbk; ++k)
          BW.bk_start[k + 1] += BW.bk_start[k];
        BW.l_z.resize(nl);
        BW.l_kc.resize(nl);
        {
          // fill through a cursor copy of the starts
          BW.l_tmp2.assign(BW.bk_start.begin(), BW.bk_start.end() - 1);
          for (const auto &e : BW.l_tmp) {
            const unsigned int pos = BW.l_tmp2[bucket(e.first)]++;
            BW.l_z[pos] = e.first;
            BW.l_kc[pos] = e.second;
          }
        }
        const unsigned int *bks = BW.bk_start.data();
        tick(2);

        for (unsigned int d = db; d < de; ++d) {
          const unsigned int i0 = bks[bucket(W.d_zlo[d])];
          const unsigned int i1 = bks[bucket(W.d_zhi[d]) + 1];
          cnt.c_touched += i1 - i0;
          if (i0 == i1)
            continue;
          const unsigned int ka = W.d_ka[d];
          const float za = ga.z_[ka], rra = ga.r_[ka];
          const double cot = W.d_cot[d];
          const float slope = W.d_slope[d], kab = W.d_kab[d], dinv_ab = W.d_dinv[d];
          StagedWork::fit(W.t_d, nt + (i1 - i0));
          StagedWork::fit(W.t_kc, nt + (i1 - i0));
          for (unsigned int i = i0; i < i1; ++i) {
            const unsigned int kc = BW.l_kc[i];
            const float rrc = gc.r_[kc];
            bool pass = true;
            if (P.phi_lin) {
              const float pcph = gc.phi_[kc], invc = gc.invr_[kc];
              const float wl = d0m * std::abs(invc - invb - (rrc - rrb) * kab * dinv_ab) + mlin;
              pass &= !(std::abs(wrap_pi(pcph - pbph - slope * (rrc - rrb))) > wl);
            }
            const double zpred = za + cot * (rrc - rra);
            const double dz = gc.z_[kc] - zpred;
            pass &= !(std::abs(dz) > P.qwin);
            W.t_d[nt] = d;
            W.t_kc[nt] = kc;
            nt += pass;
          }
        }
        tick(3);
      }
      cnt.triplets += nt;
      (void)rmid;
      (void)inv_rclo_f;

      tprev = clk::now();
      finish_triplets<A>(P, ga, gb, gc, gd, nt, W, out, cnt, tick);
    }
  }

}  // namespace mkfit::seeding

#endif
