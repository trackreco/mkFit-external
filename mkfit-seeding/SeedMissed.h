#ifndef mkfit_seeding_SeedMissed_h
#define mkfit_seeding_SeedMissed_h

// Why a findable track has no quad: `seedfind --why-missed`.
//
// For every findable track (SeedTruth's definition) that no quad of the run
// matches with all four hits, evaluate every cut on the track's OWN hits with
// eval_quad<A>(), the evaluator the difference tool uses, which calls the
// finder's own functions.  A track can have more than one hit in a layer
// (module overlaps); every combination is evaluated and the one that gets
// furthest down the cut sequence is kept.  The first decisive cut that fails,
// in the finder's order, is the reason.  If none fails, the finder should have
// found it and did not: counted as "no cut fails".

#include "SeedMargins.h"
#include "SeedTruth.h"

#include <array>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace mkfit::seeding {

  struct SeedMissed {
    // pipeline order of the decisive cuts; bc_lin is off in the default window
    static constexpr int kOrder[] = {MC_ab_r, MC_ab_phi, MC_bc_r, MC_bc_phi, MC_bc_lin, MC_c_z, MC_c3,
                                     MC_reach, MC_cd_r, MC_d_cross, MC_d_phi, MC_d_z};
    static constexpr int kNo = sizeof(kOrder) / sizeof(int);
    static constexpr int kNpt = 4;
    static constexpr double kPtEdge[kNpt + 1] = {0, 1.2, 2.0, 5.0, 1e9};
    static constexpr int kNeta = 3;
    static constexpr double kEtaEdge[kNeta + 1] = {0, 0.6, 1.0, 1.6};

    double n_missed = 0, n_findable = 0, n_events = 0;
    // [reason][pt bin], [reason][eta bin]; reason kNo = no cut fails
    double by_pt[kNo + 1][kNpt] = {}, by_eta[kNo + 1][kNeta] = {};
    std::vector<double> margin[kNo + 1];  // failing margin, cm or rad (negative = outside)
    double fetch_only = 0;                // all decisive cuts pass but the d fetch window does not
    double findable_pt[kNpt] = {};
    // fine pT bins for efficiency versus pT: findable and missed per bin
    static constexpr double kFine[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.35, 1.5, 1.75,
                                       2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0, 50.0};
    static constexpr int kNfine = sizeof(kFine) / sizeof(double) - 1;
    double fine_findable[kNfine] = {}, fine_missed[kNfine] = {};
    static int fine_bin(double pt) {
      if (pt < kFine[0] || pt >= kFine[kNfine])
        return -1;
      int i = 0;
      while (pt >= kFine[i + 1])
        ++i;
      return i;
    }
    // for the fixed-window cuts (c_z, d_phi, d_z): |residual| / window of the
    // first failing cut, with the track's pT bin; > 1 by construction
    std::vector<std::pair<int, double>> ratio[kNo + 1];

    template <class A, typename L>
    void event(const Event &ev, int la, int lb, int lc, int ld, const SeedParams &P, const L &ga, const L &gb,
               const L &gc, const L &gd, const std::vector<Quad> &quads) {
      n_events += 1;
      const MCHitInfoVec &mc = ev.simHitsInfo_;
      const HitVec *H[4] = {&ev.layerHits_[la], &ev.layerHits_[lb], &ev.layerHits_[lc], &ev.layerHits_[ld]};
      auto label = [&mc](const Hit &h) -> int {
        const int mid = h.mcHitID();
        return (mid >= 0 && mid < (int)mc.size()) ? mc[mid].mcTrackID() : -1;
      };
      // label -> original hit indices per layer
      std::unordered_map<int, std::array<std::vector<unsigned int>, 4>> hits;
      for (int i = 0; i < 4; ++i)
        for (unsigned int k = 0; k < H[i]->size(); ++k) {
          const int l = label((*H[i])[k]);
          if (l >= 0)
            hits[l][i].push_back(k);
        }
      // labels a run's quads find with all four hits
      std::unordered_map<int, bool> found;
      for (const auto &q : quads) {
        const int l0 = label((*H[0])[q[0]]);
        if (l0 >= 0 && l0 == label((*H[1])[q[1]]) && l0 == label((*H[2])[q[2]]) && l0 == label((*H[3])[q[3]]))
          found[l0] = true;
      }
      // original -> bin-order index
      std::vector<unsigned int> io[4];
      const L *G[4] = {&ga, &gb, &gc, &gd};
      for (int i = 0; i < 4; ++i) {
        io[i].resize(G[i]->n());
        for (unsigned int k = 0; k < G[i]->n(); ++k)
          io[i][G[i]->orig_[k]] = k;
      }
      for (const auto &kv : hits) {
        const int l = kv.first;
        const auto &hl = kv.second;
        if (hl[0].empty() || hl[1].empty() || hl[2].empty() || hl[3].empty() || l >= (int)ev.simTracks_.size())
          continue;
        const Track &t = ev.simTracks_[l];
        if (t.pT() < P.pt_min || std::hypot(t.x() - ev.beamSpot_.x, t.y() - ev.beamSpot_.y) > P.d0_max ||
            std::abs(t.momEta()) > SeedTruth::kEtaMax)
          continue;
        n_findable += 1;
        int ip = 0;
        while (ip < kNpt - 1 && t.pT() >= kPtEdge[ip + 1])
          ++ip;
        findable_pt[ip] += 1;
        const int fb = fine_bin(t.pT());
        if (fb >= 0)
          fine_findable[fb] += 1;
        if (found.count(l))
          continue;
        if (fb >= 0)
          fine_missed[fb] += 1;
        n_missed += 1;
        // best combination: the one whose first failing cut is latest in the order
        int best_stage = -1;
        double best_m = 0;
        bool best_fetch = false;
        for (unsigned int a : hl[0])
          for (unsigned int b : hl[1])
            for (unsigned int c : hl[2])
              for (unsigned int d : hl[3]) {
                QuadEval e;
                eval_quad<A>(P, ga, gb, gc, gd, io[0][a], io[1][b], io[2][c], io[3][d], e);
                int stage = kNo;
                double m = 0;
                for (int s = 0; s < kNo; ++s) {
                  const CutVal &cv = e.c[kOrder[s]];
                  if (cv.valid && !cv.pass) {
                    stage = s;
                    m = cv.m;
                    break;
                  }
                  // a cut left invalid after the earlier ones passed means the
                  // circle could not be formed or crossed: count it there
                  if (!cv.valid && kOrder[s] != MC_bc_lin && kOrder[s] != MC_bc_phi && kOrder[s] != MC_c3) {
                    stage = s;
                    m = 0;
                    break;
                  }
                }
                const bool fetch_fail = stage == kNo && ((e.c[MC_dfetch_phi].valid && !e.c[MC_dfetch_phi].pass) ||
                                                         (e.c[MC_dfetch_z].valid && !e.c[MC_dfetch_z].pass));
                if (stage > best_stage) {
                  best_stage = stage;
                  best_m = m;
                  best_fetch = fetch_fail;
                }
              }
        int ie = 0;
        const double ae = std::abs(t.momEta());
        while (ie < kNeta - 1 && ae >= kEtaEdge[ie + 1])
          ++ie;
        by_pt[best_stage][ip] += 1;
        by_eta[best_stage][ie] += 1;
        margin[best_stage].push_back(best_m);
        if (best_stage < kNo) {
          const int c = kOrder[best_stage];
          const double w = c == MC_c_z ? P.qwin : c == MC_d_phi ? P.phiwin_d : c == MC_d_z ? P.qwin_d : 0;
          if (w > 0)
            ratio[best_stage].push_back({ip, 1.0 - best_m / w});
        }
        if (best_fetch)
          fetch_only += 1;
      }
    }

    void report() const {
      auto name = [](int s) -> std::string { return s == kNo ? "no cut fails" : margin_cut_info(kOrder[s]).name; };
      printf("   why missed: %.1f of %.1f findable tracks /ev have no 4-of-4 quad (%.4f)\n", n_missed / n_events,
             n_findable / n_events, n_missed / std::max(1.0, n_findable));
      printf("   %-14s %8s %7s |  pT <1.2  1.2-2   2-5    >5 | |eta| <0.6 0.6-1 1-1.6 | median fail margin\n", "first failing",
             "/ev", "share");
      for (int s = 0; s <= kNo; ++s) {
        double tot = 0;
        for (int i = 0; i < kNpt; ++i)
          tot += by_pt[s][i];
        if (tot == 0)
          continue;
        std::vector<double> m = margin[s];
        std::sort(m.begin(), m.end());
        const bool rad = s < kNo && margin_cut_info(kOrder[s]).rad;
        printf("   %-14s %8.2f %7.3f | %6.2f %6.2f %6.2f %6.2f | %6.2f %6.2f %6.2f | %+.3g %s\n", name(s).c_str(),
               tot / n_events, tot / std::max(1.0, n_missed), by_pt[s][0] / n_events, by_pt[s][1] / n_events,
               by_pt[s][2] / n_events, by_pt[s][3] / n_events, by_eta[s][0] / n_events, by_eta[s][1] / n_events,
               by_eta[s][2] / n_events, m.empty() ? 0.0 : m[m.size() / 2], s == kNo ? "" : rad ? "rad" : "cm");
      }
      printf("   (of the 'no cut fails': %.2f /ev fail only the d-hit fetch window)\n", fetch_only / n_events);
      printf("   missed fraction by pT:");
      for (int i = 0; i < kNpt; ++i) {
        double m = 0;
        for (int s = 0; s <= kNo; ++s)
          m += by_pt[s][i];
        printf("  [%g,%g) %.2f/%.2f = %.4f", kPtEdge[i], kPtEdge[i + 1] > 1e8 ? 99 : kPtEdge[i + 1], m / n_events,
               findable_pt[i] / n_events, m / std::max(1.0, findable_pt[i]));
      }
      printf("\n   efficiency vs pT, for plotting: effpt lo hi findable missed (totals over all events)\n");
      for (int i = 0; i < kNfine; ++i)
        if (fine_findable[i] > 0)
          printf("   effpt %g %g %.0f %.0f\n", kFine[i], kFine[i + 1], fine_findable[i], fine_missed[i]);
      printf("   |residual| / window of the failing cut (fixed windows only), by pT bin:\n");
      printf("   %-8s %-10s %6s %6s %6s %6s %7s %7s\n", "cut", "pT", "n/ev", "p50", "p90", "max", "f(>2)", "f(>5)");
      for (int s = 0; s < kNo; ++s) {
        if (ratio[s].empty())
          continue;
        for (int i = -1; i < kNpt; ++i) {
          std::vector<double> v;
          for (auto &pr : ratio[s])
            if (i < 0 || pr.first == i)
              v.push_back(pr.second);
          if (v.empty())
            continue;
          std::sort(v.begin(), v.end());
          double f2 = 0, f5 = 0;
          for (double x : v) {
            f2 += x > 2;
            f5 += x > 5;
          }
          char pt[32];
          if (i < 0)
            snprintf(pt, sizeof pt, "all");
          else
            snprintf(pt, sizeof pt, "%g-%g", kPtEdge[i], kPtEdge[i + 1] > 1e8 ? 99 : kPtEdge[i + 1]);
          printf("   %-8s %-10s %6.2f %6.2f %6.2f %6.1f %7.3f %7.3f\n", i < 0 ? name(s).c_str() : "", pt,
                 v.size() / n_events, v[v.size() / 2], v[(size_t)(0.9 * (v.size() - 1))], v.back(), f2 / v.size(),
                 f5 / v.size());
        }
      }
    }
  };

}  // namespace mkfit::seeding

#endif
