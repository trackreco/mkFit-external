#ifndef mkfit_seeding_SeedFinderTile_h
#define mkfit_seeding_SeedFinderTile_h

// The per-pair stages as three kernels per b-hit, each a loop over contiguous
// hits against one broadcast b-hit (step B, README "Step B0"):
//
//   K1  the a-window (one contiguous run, two at the phi wrap): the doublet
//       tests, and the r-z slope s_a = (z_b - z_a)/(r_b - r_a) of each;
//   K3  the c-window, from a layer c binned in PHI ONLY so it is one run too:
//       the b-c tests, the slope s_c = (z_c - z_b)/(r_c - r_b) and the
//       tolerance t_c = q_win/(r_c - r_b);
//   K4  |s_a - s_c| <= t_c over doublets x list.  This IS the layer-c z test,
//       divided by dr_a dr_c > 0.  Run as a conservative pre-filter (t_c
//       carries a slack for the float rounding); every hit is re-tested with
//       the finder's exact c_z expression.  Two forms:
//         brute    the dense tile, 8 list entries per AVX instruction.  Measured
//                  ~9 cycles per 8 pairs on Zen+, and 67 entries per doublet
//                  for 0.058 hits: the volume is the cost.
//         buckets  (default) the list counting-sorted into slope buckets of
//                  width 0.5 over +-16, clamped at the edges.  With t_c <=
//                  t_max ~0.01 a doublet tests the ONE contiguous range
//                  [bucket(s_a - t_max), bucket(s_a + t_max)], ~3 entries, in
//                  one masked 8-wide step.  The bucket map is monotone, so the
//                  range holds every entry within t_max: still conservative.
//
// K1 and K3 compute the b-major finder's expressions exactly, as vectorisable
// mask loops followed by a scalar compaction, so doublets and list are the
// same.  K4 plus its confirm give the same triplets.  Stages 5-7 are
// finish_triplets<A>(), with the phi-only layer c as its layer c: they only
// read per-hit arrays through the index, and orig_ maps back to the same hits.
//
// phi_lin 0 (the generic b-c band) and 2 (that band AND the layer-c phi
// predicted linearly in r from the doublet's own phi slope).  Mode 2 applies
// the linear cut in the confirm step, next to the exact c_z test, in the
// b-major finder's float expression, so the two give the same triplets.

#include "SeedFinderBMajor.h"

#if defined(__AVX__)
#include <immintrin.h>
#endif

namespace mkfit::seeding {

  struct TileWork {
    std::vector<unsigned int> b_dbeg;
    std::vector<float> d_sa;                    // per doublet of the block: slope
    std::vector<unsigned short> d_blo, d_bhi;   // per doublet: slope-bucket range [blo, bhi)
    std::vector<float> t_s, t_t;                // K1/K3 scratch per fetched hit: slope, tolerance
    std::vector<unsigned short> t_b0, t_b1;     // K1/K3 scratch: bucket indices
    std::vector<unsigned char> t_m;             // K1/K3 scratch: pass mask
    std::vector<float> l_sc, l_tc;              // per b-hit list, compacted
    std::vector<unsigned int> l_kc;
    std::vector<unsigned short> l_bk;
    std::vector<float> s_sc, s_tc;              // the list sorted into slope buckets, padded
    std::vector<unsigned int> s_kc, bk_start, bk_cur;
    std::vector<unsigned int> h_d, h_j;         // K4 phase 1: (doublet, first list slot) of a nonzero mask
    std::vector<unsigned char> h_m;             //             and the mask
    bool brute = false;
    long k4_cand = 0;  // pre-filter hits, for the counters
  };

  namespace detail {
    // slope slack of the K4 pre-filter [dimensionless]: float rounding of s_a,
    // s_c and t_c is ~1e-6 at |s| ~ 20; the exact test decides after it
    constexpr float kSlopeSlack = 1e-4f;

    // One masked 8-wide slope step over list slots [j, j + 8): bit k is set if
    // |sa - s_c| <= t_c for slot j + k; lanes at or past `end` are cleared.
    inline unsigned int k4_step(const float *__restrict ssc,
                                const float *__restrict stc,
                                float sa,
                                unsigned int j,
                                unsigned int end) {
      const unsigned int left = std::min(end - j, 8u);
#if defined(__AVX__)
      const __m256 vabs = _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
      const __m256 diff = _mm256_and_ps(_mm256_sub_ps(_mm256_set1_ps(sa), _mm256_loadu_ps(ssc + j)), vabs);
      const unsigned int m = _mm256_movemask_ps(_mm256_cmp_ps(diff, _mm256_loadu_ps(stc + j), _CMP_LE_OQ));
#else
      unsigned int m = 0;
      for (unsigned int k = 0; k < 8; ++k)
        m |= (unsigned int)(std::abs(sa - ssc[j + k]) <= stc[j + k]) << k;
#endif
      return m & (0xffu >> (8 - left));
    }

    // The same step without the lane mask, for the first step of a doublet in
    // bucket mode.  The mask is not needed there: the slots past the range
    // end hold higher slope buckets, whose slopes exceed sa + t_max >= sa + t_c
    // (the bucket map is monotone), and past the list end come 8 padding
    // entries that never pass.  So those lanes are false anyway.  The slope
    // is broadcast straight from memory.
    inline unsigned int k4_step_unmasked(const float *__restrict ssc,
                                         const float *__restrict stc,
                                         const float *sa,
                                         unsigned int j) {
#if defined(__AVX__)
      const __m256 vabs = _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
      const __m256 diff = _mm256_and_ps(_mm256_sub_ps(_mm256_broadcast_ss(sa), _mm256_loadu_ps(ssc + j)), vabs);
      return _mm256_movemask_ps(_mm256_cmp_ps(diff, _mm256_loadu_ps(stc + j), _CMP_LE_OQ));
#else
      unsigned int m = 0;
      for (unsigned int k = 0; k < 8; ++k)
        m |= (unsigned int)(std::abs(*sa - ssc[j + k]) <= stc[j + k]) << k;
      return m;
#endif
    }

    // The steps after the first, for a range longer than 8 slots.  Rare in
    // bucket mode, so it lives out of line and its values do not compete for
    // registers with the hot loop.
    __attribute__((noinline, cold)) inline unsigned int k4_more_steps(const float *__restrict ssc,
                                                                      const float *__restrict stc,
                                                                      float sa,
                                                                      unsigned int d,
                                                                      unsigned int j,
                                                                      unsigned int i1,
                                                                      unsigned int *__restrict hd,
                                                                      unsigned int *__restrict hj,
                                                                      unsigned char *__restrict hm,
                                                                      unsigned int nh) {
      for (; j < i1; j += 8) {
        const unsigned int m = k4_step(ssc, stc, sa, j, i1);
        hd[nh] = d;
        hj[nh] = j;
        hm[nh] = (unsigned char)m;
        nh += m != 0;
      }
      return nh;
    }

    // K4 phase 1 over the doublets [db, de) of one b-hit: per doublet, one
    // masked 8-wide slope test over its list range, [bs[blo], bs[bhi]) with
    // slope buckets or the whole list [0, nl) for the brute tile.  Records
    // (doublet, first slot, mask) where the mask is nonzero and returns how
    // many; hd/hj/hm must hold (de - db) * ((nl + 7) / 8 + 1) entries.
    //
    // A function of its own, with every array a __restrict argument and every
    // counter a local.  Inlined into find_quads_tile, the byte store to hm may
    // alias TileWork's members, so GCC reloaded every array pointer from
    // TileWork per doublet and kept the loop counter on the stack.
    template <bool Brute>
    __attribute__((noinline)) unsigned int k4_phase1(unsigned int db,
                                                     unsigned int de,
                                                     unsigned int nl,
                                                     const unsigned int *__restrict bs,
                                                     const unsigned short *__restrict blo,
                                                     const unsigned short *__restrict bhi,
                                                     const float *__restrict dsa,
                                                     const float *__restrict ssc,
                                                     const float *__restrict stc,
                                                     unsigned int *__restrict hd,
                                                     unsigned int *__restrict hj,
                                                     unsigned char *__restrict hm) {
      unsigned int nh = 0;
      for (unsigned int d = db; d < de; ++d) {
        const unsigned int i0 = Brute ? 0 : bs[blo[d]];
        const unsigned int i1 = Brute ? nl : bs[bhi[d]];
        const float sa = dsa[d];
        const unsigned int m = Brute ? k4_step(ssc, stc, sa, i0, i1) : k4_step_unmasked(ssc, stc, dsa + d, i0);
        // Stored only where the mask is nonzero, 5.9 % of doublets.  Storing
        // always, at slot nh, made every store address wait for the previous
        // doublet's mask, at the end of the longest chain in the loop.
        if (m) {
          hd[nh] = d;
          hj[nh] = i0;
          hm[nh] = (unsigned char)m;
          ++nh;
        }
        if (__builtin_expect(i0 + 8 < i1, 0))
          nh = k4_more_steps(ssc, stc, sa, d, i0 + 8, i1, hd, hj, hm, nh);
      }
      return nh;
    }
  }  // namespace detail

  template <class A = ArithFastK, typename LA, typename LB, typename LC, typename LD>
  void find_quads_tile(const SeedParams &P,
                       const LA &ga,
                       const LB &gb,
                       const LC &gc,  // phi-only binning
                       const LD &gd,
                       std::vector<Quad> &out,
                       SeedCounters &cnt,
                       StagedWork &W,
                       TileWork &TW,
                       unsigned int block = 32) {
    using namespace detail;
    if (P.phi_lin == 1) {
      fprintf(stderr, "find_quads_tile: phi_lin mode 1 not supported (the c-list is the generic band)\n");
      return;
    }
    constexpr float kBfield = 3.8f;
    const double Rmin = P.pt_min / (0.003f * kBfield);
    const float inv2R = (float)(1.0 / (2 * Rmin)), d0m = P.d0_max, marg = P.phi_margin;
    const double ralo = ga.rlo_, rchi = gc.rhi_;
    const float inv_ralo = 1.0f / (float)ralo, inv_rchi = 1.0f / (float)rchi;
    const float qwin = P.qwin;
    const float mlin = P.phi_lin_marg;
    auto wphi_w = [=](float r_in, float inv_in, float r_out, float inv_out) {
      return (r_out - r_in) * inv2R + d0m * (inv_in - inv_out) + marg;
    };
    auto phi_bins = [](const auto &L, float c, float w) {
      w = std::min(w, 0.9f * kPi);
      return L.phi_range(c - w, c + w);
    };
    constexpr float kFetchEps = 1e-6f;
    // slope buckets over +-kSMax, clamped; t_c <= t_max from the smallest r_c - r_b
    constexpr float kSMax = 16.0f, kBw = 0.5f;
    constexpr unsigned int kNb = (unsigned int)(2 * kSMax / kBw);
    const float t_max = qwin / std::max(0.1f, (float)(gc.rlo_ - gb.rhi_)) + kSlopeSlack;
    TW.bk_start.resize(kNb + 1);
    TW.bk_cur.resize(kNb);

    // [begin, end) runs of a phi window over all q: one, or two at the wrap
    auto runs = [](const auto &L, auto p, unsigned int r[4]) -> int {
      int k = 0;
      L.for_each_run(p, L.q_all(), [&](unsigned int b, unsigned int e) {
        r[2 * k] = b;
        r[2 * k + 1] = e;
        ++k;
      });
      return k;
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

      //---- K1: doublets, their slopes and slope-bucket ranges
      unsigned int nd = 0;
      TW.b_dbeg.resize(kb1 - kb0 + 1);
      for (unsigned int kb = kb0; kb < kb1; ++kb) {
        TW.b_dbeg[kb - kb0] = nd;
        const float pbph = gb.phi_[kb], zb = gb.z_[kb], rrb = gb.r_[kb], invb = gb.invr_[kb];
        const float w = wphi_w((float)ralo, inv_ralo, rrb, invb) + kFetchEps;
        unsigned int rr[4];
        const int nr = runs(ga, phi_bins(ga, pbph, w), rr);
        for (int ir = 0; ir < nr; ++ir) {
          const unsigned int a0 = rr[2 * ir], n = rr[2 * ir + 1] - a0;
          StagedWork::fit(TW.t_s, n);
          StagedWork::fit(TW.t_m, n);
          StagedWork::fit(TW.t_b0, n);
          StagedWork::fit(TW.t_b1, n);
          const float *__restrict ar = ga.r_.data() + a0;
          const float *__restrict ap = ga.phi_.data() + a0;
          const float *__restrict ai = ga.invr_.data() + a0;
          const float *__restrict az = ga.z_.data() + a0;
          float *__restrict ts = TW.t_s.data();
          unsigned char *__restrict tm = TW.t_m.data();
          unsigned short *__restrict tb0 = TW.t_b0.data();
          unsigned short *__restrict tb1 = TW.t_b1.data();
#pragma omp simd
          for (unsigned int i = 0; i < n; ++i) {
            const float rra = ar[i];
            const bool pass =
                !(rrb <= rra + 0.1f) & !(std::abs(wrap_pi(pbph - ap[i])) > wphi_w(rra, ai[i], rrb, invb));
            tm[i] = pass;
            const float sa = (zb - az[i]) / (rrb - rra);
            ts[i] = sa;
            const float y0 = std::min(std::max((sa - t_max + kSMax) * (1.0f / kBw), 0.0f), float(kNb - 1));
            const float y1 = std::min(std::max((sa + t_max + kSMax) * (1.0f / kBw), 0.0f), float(kNb - 1));
            tb0[i] = (unsigned short)(int)y0;
            tb1[i] = (unsigned short)((int)y1 + 1);
          }
          StagedWork::fit(W.d_ka, nd + n);
          StagedWork::fit(W.d_kb, nd + n);
          StagedWork::fit(TW.d_sa, nd + n);
          StagedWork::fit(TW.d_blo, nd + n);
          StagedWork::fit(TW.d_bhi, nd + n);
          unsigned int *__restrict oka = W.d_ka.data();
          unsigned int *__restrict okb = W.d_kb.data();
          float *__restrict osa = TW.d_sa.data();
          unsigned short *__restrict olo = TW.d_blo.data();
          unsigned short *__restrict ohi = TW.d_bhi.data();
          unsigned int m = nd;  // a local count: a captured one is kept in memory
          for (unsigned int i = 0; i < n; ++i) {
            oka[m] = a0 + i;
            osa[m] = ts[i];
            olo[m] = tb0[i];
            ohi[m] = tb1[i];
            m += tm[i];
          }
          std::fill(okb + nd, okb + m, kb);
          nd = m;
        }
      }
      TW.b_dbeg[kb1 - kb0] = nd;
      cnt.doublets += nd;
      tick(0);

      //---- K3 + K4 per b-hit
      unsigned int nt = 0;
      for (unsigned int kb = kb0; kb < kb1; ++kb) {
        const unsigned int db = TW.b_dbeg[kb - kb0], de = TW.b_dbeg[kb - kb0 + 1];
        if (db == de)
          continue;
        const float pbph = gb.phi_[kb], zb = gb.z_[kb], rrb = gb.r_[kb], invb = gb.invr_[kb];
        const float w_bc = wphi_w(rrb, invb, (float)rchi, inv_rchi);
        unsigned int nl = 0;
        unsigned int rr[4];
        const int nr = runs(gc, phi_bins(gc, pbph, w_bc), rr);
        for (int ir = 0; ir < nr; ++ir) {
          const unsigned int c0 = rr[2 * ir], n = rr[2 * ir + 1] - c0;
          StagedWork::fit(TW.t_s, n);
          StagedWork::fit(TW.t_t, n);
          StagedWork::fit(TW.t_m, n);
          StagedWork::fit(TW.t_b0, n);
          const float *__restrict cr = gc.r_.data() + c0;
          const float *__restrict cp = gc.phi_.data() + c0;
          const float *__restrict ci = gc.invr_.data() + c0;
          const float *__restrict cz = gc.z_.data() + c0;
          float *__restrict ts = TW.t_s.data();
          float *__restrict tt = TW.t_t.data();
          unsigned char *__restrict tm = TW.t_m.data();
          unsigned short *__restrict tb = TW.t_b0.data();
#pragma omp simd
          for (unsigned int i = 0; i < n; ++i) {
            const float rrc = cr[i];
            const bool pass =
                !(rrc <= rrb + 0.1f) & !(std::abs(wrap_pi(cp[i] - pbph)) > wphi_w(rrb, invb, rrc, ci[i]));
            tm[i] = pass;
            const float idr = 1.0f / (rrc - rrb);
            const float sc = (cz[i] - zb) * idr;
            ts[i] = sc;
            tt[i] = qwin * idr + kSlopeSlack;
            tb[i] = (unsigned short)(int)std::min(std::max((sc + kSMax) * (1.0f / kBw), 0.0f), float(kNb - 1));
          }
          StagedWork::fit(TW.l_sc, nl + n);
          StagedWork::fit(TW.l_tc, nl + n);
          StagedWork::fit(TW.l_kc, nl + n);
          StagedWork::fit(TW.l_bk, nl + n);
          float *__restrict osc = TW.l_sc.data();
          float *__restrict otc = TW.l_tc.data();
          unsigned int *__restrict okc = TW.l_kc.data();
          unsigned short *__restrict obk = TW.l_bk.data();
          unsigned int m = nl;
          for (unsigned int i = 0; i < n; ++i) {
            osc[m] = ts[i];
            otc[m] = tt[i];
            okc[m] = c0 + i;
            obk[m] = tb[i];
            m += tm[i];
          }
          nl = m;
        }
        cnt.c_touched += nl;

        // the list, sorted into slope buckets (or as is, for the brute tile),
        // padded with 8 entries that can never pass
        StagedWork::fit(TW.s_sc, nl + 8);
        StagedWork::fit(TW.s_tc, nl + 8);
        StagedWork::fit(TW.s_kc, nl + 8);
        float *__restrict ssc = TW.s_sc.data();
        float *__restrict stc = TW.s_tc.data();
        unsigned int *__restrict skc = TW.s_kc.data();
        unsigned int *__restrict bs = TW.bk_start.data();
        if (!TW.brute) {
          const unsigned short *__restrict lbk = TW.l_bk.data();
          unsigned int *__restrict cur = TW.bk_cur.data();
          // histogram, and the prefix sum over the occupied buckets only
          for (unsigned int k = 0; k <= kNb; ++k)
            bs[k] = 0;
          unsigned int bmin = kNb, bmax = 0;
          for (unsigned int j = 0; j < nl; ++j) {
            const unsigned int b = lbk[j];
            ++bs[b + 1];
            bmin = std::min(bmin, b);
            bmax = std::max(bmax, b);
          }
          if (nl == 0)
            bmin = bmax = 0;
          for (unsigned int k = bmin; k <= bmax; ++k) {
            bs[k + 1] += bs[k];
            cur[k] = bs[k];
          }
          for (unsigned int k = bmax + 2; k <= kNb; ++k)
            bs[k] = nl;
          for (unsigned int j = 0; j < nl; ++j) {
            const unsigned int pos = cur[lbk[j]]++;
            ssc[pos] = TW.l_sc[j];
            stc[pos] = TW.l_tc[j];
            skc[pos] = TW.l_kc[j];
          }
        } else {
          std::copy(TW.l_sc.begin(), TW.l_sc.begin() + nl, ssc);
          std::copy(TW.l_tc.begin(), TW.l_tc.begin() + nl, stc);
          std::copy(TW.l_kc.begin(), TW.l_kc.begin() + nl, skc);
        }
        for (unsigned int j = nl; j < nl + 8; ++j) {
          ssc[j] = 1e30f;
          stc[j] = -1.0f;
        }
        tick(2);

        //---- K4: the slope pre-filter, then the exact c_z test on its hits
        auto confirm = [&](unsigned int d, unsigned int j) {
          ++TW.k4_cand;
          const unsigned int ka = W.d_ka[d], kc = skc[j];
          const float za = ga.z_[ka], rra = ga.r_[ka];
          const double cot = (zb - za) / (rrb - rra);
          const double zpred = za + cot * (gc.r_[kc] - rra);
          const double dz = gc.z_[kc] - zpred;
          if (std::abs(dz) > P.qwin)
            return;
          if (P.phi_lin) {
            // the layer-c phi, linear in r from the doublet's phi slope
            const float rrc = gc.r_[kc], invc = gc.invr_[kc], inva = ga.invr_[ka];
            const float kab = 1.0f / (rrb - rra);
            const float slope = wrap_pi(pbph - ga.phi_[ka]) / (rrb - rra);
            const float wl = d0m * std::abs(invc - invb - (rrc - rrb) * kab * (invb - inva)) + mlin;
            if (std::abs(wrap_pi(gc.phi_[kc] - pbph - slope * (rrc - rrb))) > wl)
              return;
          }
          StagedWork::fit(W.t_d, nt + 1);
          StagedWork::fit(W.t_kc, nt + 1);
          W.t_d[nt] = d;
          W.t_kc[nt] = kc;
          ++nt;
        };
        // phase 1: per doublet, record (d, j, mask) only where the mask is nonzero
        const unsigned int hcap = (de - db) * ((nl + 7) / 8 + 1);
        StagedWork::fit(TW.h_d, hcap);
        StagedWork::fit(TW.h_j, hcap);
        StagedWork::fit(TW.h_m, hcap);
        const unsigned int nh =
            (TW.brute ? k4_phase1<true> : k4_phase1<false>)(db, de, nl, bs, TW.d_blo.data(), TW.d_bhi.data(),
                                                             TW.d_sa.data(), ssc, stc, TW.h_d.data(),
                                                             TW.h_j.data(), TW.h_m.data());
        // phase 2: the exact z test on each candidate
        for (unsigned int h = 0; h < nh; ++h) {
          unsigned int m = TW.h_m[h];
          while (m) {
            confirm(TW.h_d[h], TW.h_j[h] + __builtin_ctz(m));
            m &= m - 1;
          }
        }
        tick(3);
      }
      cnt.triplets += nt;

      tprev = clk::now();
      finish_triplets<A>(P, ga, gb, gc, gd, nt, W, out, cnt, tick);
    }
  }

}  // namespace mkfit::seeding

#endif
