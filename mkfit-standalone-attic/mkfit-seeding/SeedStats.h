#ifndef mkfit_seeding_SeedStats_h
#define mkfit_seeding_SeedStats_h

// Step B0: the numbers the vectorised per-pair stages have to be designed
// around.  Walks the b-major finder's windows and cuts (same expressions) and
// records, instead of finding quads:
//
//   trip counts   what each loop would iterate over: a-hits per b window, the
//                 number of contiguous runs, doublets per b-hit, c-hits per b
//                 window and per b list, bucket candidates per doublet, and
//                 the dense tile (doublets x list) a b-hit would present;
//   value ranges  every delta the fixed-point tests would take, RELATIVE TO THE
//                 b-HIT: dz, dr, dphi of the a side (per doublet) and of the c
//                 side (per list entry, and per bucket candidate), and the z
//                 determinant with its threshold.
//
// Not timed, not on the fast path.  Its doublet and triplet counts must equal
// the finder's; seedfind prints both.

#include "SeedFinderBMajor.h"

#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

namespace mkfit::seeding {

  struct Dist {
    std::string name, unit;
    std::vector<double> v;
    void add(double x) { v.push_back(x); }
    void print(FILE *f = stdout) {
      if (v.empty()) {
        fprintf(f, "     %-22s %6s  (empty)\n", name.c_str(), unit.c_str());
        return;
      }
      std::sort(v.begin(), v.end());
      double s = 0;
      for (double x : v)
        s += x;
      auto q = [&](double p) { return v[std::min(v.size() - 1, (size_t)(p * v.size()))]; };
      fprintf(f, "     %-22s %6s %10zu %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n", name.c_str(), unit.c_str(),
              v.size(), s / v.size(), v.front(), q(0.5), q(0.9), q(0.99), q(0.999), v.back());
    }
  };

  struct SeedStats {
    long doublets = 0, triplets = 0;
    // trip counts
    Dist a_fetch{"a fetched / b-hit", "hits"}, a_runs{"a runs / b-hit", "runs"}, d_per_b{"doublets / b-hit", "n"};
    Dist c_fetch{"c fetched / b-hit", "hits"}, c_list{"c list / b-hit", "hits"}, c_runs{"c runs / b-hit", "runs"};
    Dist cand_per_d{"bucket cands / doublet", "hits"}, trip_per_d{"triplets / doublet", "n"};
    Dist tile{"tile nd x nl / b-hit", "pairs"}, bucket_sum{"bucket cands / b-hit", "pairs"};
    // value ranges, relative to the b-hit
    Dist dz_a{"|zb - za| doublet", "cm"}, dr_a{"rb - ra doublet", "cm"}, dphi_a{"|dphi_ab| doublet", "rad"};
    Dist dz_c_list{"|zc - zb| list", "cm"}, dz_c_cand{"|zc - zb| cand", "cm"}, dr_c{"rc - rb list", "cm"};
    Dist dphi_c{"|dphi_bc| list", "rad"}, dphi_c_f{"|dphi_bc| fetched", "rad"};
    Dist det_z{"|det_z| cand", "cm2"}, thr_z{"qwin*dr_a doublet", "cm2"};
    Dist z_abs_a{"|z| layer a", "cm"}, z_abs_c{"|z| layer c", "cm"};
    // the z test as a band test on the r-z slope relative to b:
    //   |s_c - s_a| < qwin / dr_c,  s = dz / dr   (exactly the determinant / (dr_a dr_c))
    Dist s_a{"|s_a| doublet", "1"}, s_c{"|s_c| list", "1"}, tol_c{"qwin/dr_c list", "1"};
    Dist cand_slope{"slope cands / doublet", "hits"}, cand_slope_b{"slope cands / b-hit", "pairs"};
    long trip_slope = 0;

    template <typename F>
    void each(F &&f) {
      for (Dist *d : {&a_fetch, &a_runs, &d_per_b, &c_fetch, &c_runs, &c_list, &cand_per_d, &trip_per_d, &tile,
                      &bucket_sum, &dz_a, &dr_a, &dphi_a, &dz_c_list, &dz_c_cand, &dr_c, &dphi_c, &dphi_c_f,
                      &det_z, &thr_z, &z_abs_a, &z_abs_c, &s_a, &s_c, &tol_c, &cand_slope, &cand_slope_b})
        f(*d);
    }
    void report() {
      printf("[stats] doublets %ld, triplets %ld (must equal the finder's); slope-form triplets %ld\n", doublets,
             triplets, trip_slope);
      printf("     %-22s %6s %10s %10s %10s %10s %10s %10s %10s %10s\n", "quantity", "unit", "n", "mean", "min",
             "p50", "p90", "p99", "p99.9", "max");
      each([](Dist &d) { d.print(); });
    }
  };

  template <typename LA, typename LB, typename LC, typename LD>
  void seed_stats(const SeedParams &P, const LA &ga, const LB &gb, const LC &gc, const LD &gd, SeedStats &S) {
    using namespace detail;
    constexpr float kBfield = 3.8f;
    const double Rmin = P.pt_min / (0.003f * kBfield);
    const float inv2R = (float)(1.0 / (2 * Rmin)), d0m = P.d0_max, marg = P.phi_margin;
    const double ralo = ga.rlo_, rclo = gc.rlo_, rchi = gc.rhi_;
    const float inv_ralo = 1.0f / (float)ralo, inv_rchi = 1.0f / (float)rchi;
    auto wphi_w = [=](float r_in, float inv_in, float r_out, float inv_out) {
      return (r_out - r_in) * inv2R + d0m * (inv_in - inv_out) + marg;
    };
    auto phi_bins = [](const auto &L, float c, float w) {
      w = std::min(w, 0.9f * kPi);
      return L.phi_range(c - w, c + w);
    };
    constexpr float kFetchEps = 1e-6f;
    const float bz0 = gc.q_min(), binv = 1.0f / P.lbin;
    const unsigned int nbk = std::max(1u, (unsigned int)std::ceil((gc.q_max() - bz0) * binv));
    auto bucket = [=](float z) -> unsigned int {
      const int b = (int)((z - bz0) * binv);
      return b < 0 ? 0u : (b >= (int)nbk ? nbk - 1 : (unsigned int)b);
    };
    for (unsigned int k = 0; k < ga.n(); ++k)
      S.z_abs_a.add(std::abs(ga.z_[k]));
    for (unsigned int k = 0; k < gc.n(); ++k)
      S.z_abs_c.add(std::abs(gc.z_[k]));

    std::vector<unsigned int> da;  // doublets of this b-hit: a indices
    std::vector<unsigned int> lk;  // list of this b-hit: c indices
    for (unsigned int kb = 0; kb < gb.n(); ++kb) {
      const float pbph = gb.phi_[kb], zb = gb.z_[kb], rrb = gb.r_[kb], invb = gb.invr_[kb];
      //---- stage 1
      da.clear();
      unsigned int nfa = 0, nra = 0;
      const float w = wphi_w((float)ralo, inv_ralo, rrb, invb) + kFetchEps;
      ga.for_each_run(phi_bins(ga, pbph, w), ga.q_all(), [&](unsigned int a0, unsigned int a1) {
        nfa += a1 - a0;
        nra += a1 > a0;
        for (unsigned int ka = a0; ka < a1; ++ka) {
          const float rra = ga.r_[ka];
          if (!(rrb <= rra + 0.1f) &&
              !(std::abs(wrap_pi(pbph - ga.phi_[ka])) > wphi_w(rra, ga.invr_[ka], rrb, invb)))
            da.push_back(ka);
        }
      });
      S.a_fetch.add(nfa);
      S.a_runs.add(nra);
      S.d_per_b.add(da.size());
      S.doublets += da.size();
      if (da.empty())
        continue;
      for (unsigned int ka : da) {
        S.dz_a.add(std::abs(zb - ga.z_[ka]));
        S.dr_a.add(rrb - ga.r_[ka]);
        S.dphi_a.add(std::abs(wrap_pi(pbph - ga.phi_[ka])));
        S.thr_z.add(P.qwin * (rrb - ga.r_[ka]));
      }
      //---- stage 3
      lk.clear();
      unsigned int nfc = 0, nrc = 0;
      const float w_bc = wphi_w(rrb, invb, (float)rchi, inv_rchi);
      gc.for_each_run(phi_bins(gc, pbph, w_bc), gc.q_all(), [&](unsigned int c0, unsigned int c1) {
        nfc += c1 - c0;
        nrc += c1 > c0;
        for (unsigned int kc = c0; kc < c1; ++kc) {
          S.dphi_c_f.add(std::abs(wrap_pi(gc.phi_[kc] - pbph)));
          const float rrc = gc.r_[kc];
          if (rrc <= rrb + 0.1f)
            continue;
          if (std::abs(wrap_pi(gc.phi_[kc] - pbph)) > wphi_w(rrb, invb, rrc, gc.invr_[kc]))
            continue;
          lk.push_back(kc);
        }
      });
      S.c_fetch.add(nfc);
      S.c_runs.add(nrc);
      S.c_list.add(lk.size());
      S.tile.add(double(da.size()) * lk.size());
      for (unsigned int kc : lk) {
        S.dz_c_list.add(std::abs(gc.z_[kc] - zb));
        S.dr_c.add(gc.r_[kc] - rrb);
        S.dphi_c.add(std::abs(wrap_pi(gc.phi_[kc] - pbph)));
      }
      //---- stage 4, bucketed as the finder does
      long bsum = 0;
      for (unsigned int ka : da) {
        const float za = ga.z_[ka], rra = ga.r_[ka];
        const double cot = (zb - za) / (rrb - rra);
        const double z0 = za + cot * (rclo - rra), z1 = za + cot * (rchi - rra);
        const double zh = P.qwin + kCoverEps;
        const unsigned int blo = bucket((float)(std::min(z0, z1) - zh)), bhi = bucket((float)(std::max(z0, z1) + zh));
        unsigned int nc = 0, nt = 0;
        for (unsigned int kc : lk) {
          const unsigned int bk = bucket(gc.z_[kc]);
          if (bk < blo || bk > bhi)
            continue;
          ++nc;
          const float rrc = gc.r_[kc], zc = gc.z_[kc];
          S.dz_c_cand.add(std::abs(zc - zb));
          // the determinant form of the z test, relative to b
          const double det = double(zc - zb) * (rrb - rra) - double(zb - za) * (rrc - rrb);
          S.det_z.add(std::abs(det));
          const double dz = zc - (za + cot * (rrc - rra));
          nt += !(std::abs(dz) > P.qwin);
        }
        S.cand_per_d.add(nc);
        S.trip_per_d.add(nt);
        S.triplets += nt;
        bsum += nc;
      }
      S.bucket_sum.add(bsum);
      //---- stage 4 as a slope band join: candidates within the list's LARGEST tolerance
      double tmax = 0;
      for (unsigned int kc : lk) {
        const double sc = double(gc.z_[kc] - zb) / (gc.r_[kc] - rrb), tc = P.qwin / double(gc.r_[kc] - rrb);
        S.s_c.add(std::abs(sc));
        S.tol_c.add(tc);
        tmax = std::max(tmax, tc);
      }
      long ssum = 0;
      for (unsigned int ka : da) {
        const double sa = double(zb - ga.z_[ka]) / (rrb - ga.r_[ka]);
        S.s_a.add(std::abs(sa));
        unsigned int nc = 0;
        for (unsigned int kc : lk) {
          const double sc = double(gc.z_[kc] - zb) / (gc.r_[kc] - rrb);
          if (std::abs(sc - sa) < tmax)
            ++nc;
          if (std::abs(sc - sa) < P.qwin / double(gc.r_[kc] - rrb))
            ++S.trip_slope;
        }
        S.cand_slope.add(nc);
        ssum += nc;
      }
      S.cand_slope_b.add(ssum);
    }
  }

}  // namespace mkfit::seeding

#endif
