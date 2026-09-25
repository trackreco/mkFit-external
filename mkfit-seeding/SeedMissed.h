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
        if (found.count(l))
          continue;
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
        int ip = 0, ie = 0;
        while (ip < kNpt - 1 && t.pT() >= kPtEdge[ip + 1])
          ++ip;
        const double ae = std::abs(t.momEta());
        while (ie < kNeta - 1 && ae >= kEtaEdge[ie + 1])
          ++ie;
        by_pt[best_stage][ip] += 1;
        by_eta[best_stage][ie] += 1;
        margin[best_stage].push_back(best_m);
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
    }
  };

}  // namespace mkfit::seeding

#endif
