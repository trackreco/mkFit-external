#ifndef mkfit_seeding_SeedTruth_h
#define mkfit_seeding_SeedTruth_h

// Seeding efficiency, fake rate and duplicate rate of a quad list, from the
// hits' truth links, in bins of |eta|.  `seedfind --truth OUT.txt`.
//
// The definitions are the study's (../mkfit-standalone-seedgeom/SeedGeom.cc),
// so the numbers are comparable with its README:
//
//   truth label of a hit  Hit::mcHitID() -> simHitsInfo_ -> mcTrackID(), -1 if
//                         the hit has no link (rec -> sim direction)
//   findable              a sim track with at least one hit carrying its label
//                         in EACH of the four layers, pT > pT_min, and its
//                         production point within D0_max of the beam spot in xy
//   found                 at least one quad whose four hits all carry its label
//   quad outcome          true (four equal valid labels), fake (four valid
//                         labels that disagree), undecidable (>= 1 label is -1)
//   duplicates            per found track, true quads beyond the first
//
// and, alongside, a LOOSE matching, "3 of 4": a quad matches the track that
// at least three of its hits carry; found3 needs one such quad.  A quad is
// fake3 if no label could reach three hits even with every unlabelled hit
// given to the most frequent label (max count + unlabelled < 3), and
// undecidable3 if it could.  Efficiency is quoted over the strict findable
// set AND over the wider set of tracks with a labelled hit in at least three
// of the four layers, which a 3-of-4 match can find.
//
// The strict undecidable class above is the study's and is lenient: a quad
// whose valid labels already disagree cannot be true whatever its unlabelled
// hits are.  fake4c / undecidable4c apply the same rule as the loose matching
// to 4 of 4: fake if max count + unlabelled < 4, undecidable if the valid
// labels agree and some hit is unlabelled.
//
// Efficiency and duplicates are binned in the sim track's |eta| (momentum
// direction at production); quad outcomes in the quad's own |eta|, from the r-z
// line through its a- and d-hits, so no truth is needed to place it.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

namespace mkfit::seeding {

  struct SeedTruth {
    static constexpr int kNb = 16;
    static constexpr double kEtaMax = 1.6;
    // per |eta| bin
    double den[kNb] = {}, found[kNb] = {}, dup_extra[kNb] = {}, found_dup[kNb] = {};
    double q_all[kNb] = {}, q_true[kNb] = {}, q_fake[kNb] = {}, q_undec[kNb] = {};
    // true quads of tracks outside the findable denominator, by reason
    double q_true_lowpt[kNb] = {}, q_true_displ[kNb] = {};
    double n_events = 0;
    // tracks passing the pT and D0 cuts with hits in exactly 3 of the 4 layers,
    // by the layer that has none: what a skip-one-layer variant could recover
    double miss_one[4] = {}, miss_more = 0;
    // loose, 3-of-4 matching
    double found3[kNb] = {}, dup3_extra[kNb] = {}, q_true3[kNb] = {}, q_fake3[kNb] = {}, q_undec3[kNb] = {};
    double den3w[kNb] = {}, found3w[kNb] = {};
    double q_fake4c[kNb] = {}, q_undec4c[kNb] = {};  // strict, with the consistent rule  // tracks with >= 3 of 4 layers, and found by a 3-of-4 quad

    static int bin(double aeta) {
      const int b = (int)(aeta / kEtaMax * kNb);
      return b < 0 ? 0 : (b >= kNb ? -1 : b);  // beyond the range is dropped
    }

    // quads: original hit indices within the four layers, as the finders return them
    template <typename QuadVec>
    void event(const Event &ev, int la, int lb, int lc, int ld, float pt_min, float d0_max, const QuadVec &quads) {
      n_events += 1;
      const MCHitInfoVec &mc = ev.simHitsInfo_;
      const HitVec *L[4] = {&ev.layerHits_[la], &ev.layerHits_[lb], &ev.layerHits_[lc], &ev.layerHits_[ld]};
      auto label = [&mc](const Hit &h) -> int {
        const int mid = h.mcHitID();
        return (mid >= 0 && mid < (int)mc.size()) ? mc[mid].mcTrackID() : -1;
      };

      // findable: a hit with the label in all four layers
      std::unordered_map<int, unsigned> present;
      for (int i = 0; i < 4; ++i)
        for (const Hit &h : *L[i]) {
          const int l = label(h);
          if (l >= 0)
            present[l] |= 1u << i;
        }
      std::unordered_map<int, int> n_true;   // findable label -> true quads
      std::unordered_map<int, int> n_true3;  // findable label -> 3-of-4 quads
      std::unordered_map<int, int> n_wide3;  // label with >= 3 of 4 layers -> 3-of-4 quads
      for (const auto &kv : present) {
        if (kv.first >= (int)ev.simTracks_.size())
          continue;
        const Track &t = ev.simTracks_[kv.first];
        if (t.pT() < pt_min || std::hypot(t.x() - ev.beamSpot_.x, t.y() - ev.beamSpot_.y) > d0_max ||
            std::abs(t.momEta()) > kEtaMax)
          continue;
        const int nl = __builtin_popcount(kv.second);
        if (nl >= 3)
          n_wide3[kv.first] = 0;
        if (kv.second == 15u) {
          n_true[kv.first] = 0;
          n_true3[kv.first] = 0;
          continue;
        }
        if (nl == 3)
          miss_one[__builtin_ctz(~kv.second & 15u)] += 1;
        else if (nl == 2)
          miss_more += 1;
      }

      for (const auto &q : quads) {
        const Hit &ha = (*L[0])[q[0]], &hd = (*L[3])[q[3]];
        const double cot = (hd.z() - ha.z()) / (hd.r() - ha.r());
        const int b = bin(std::abs(std::asinh(cot)));
        const int l0 = label(ha), l1 = label((*L[1])[q[1]]), l2 = label((*L[2])[q[2]]), l3 = label(hd);
        const bool undec = l0 < 0 || l1 < 0 || l2 < 0 || l3 < 0;
        const bool tru = !undec && l0 == l1 && l0 == l2 && l0 == l3;
        if (b >= 0) {
          q_all[b] += 1;
          if (undec)
            q_undec[b] += 1;
          else if (tru)
            q_true[b] += 1;
          else
            q_fake[b] += 1;
        }
        {  // loose: the label carried by at least three hits, if any
          const int ls[4] = {l0, l1, l2, l3};
          int l3of4 = -1, nun = 0, maxc = 0;
          for (int u = 0; u < 4; ++u) {
            nun += ls[u] < 0;
            if (ls[u] < 0)
              continue;
            int c = 0;
            for (int w = 0; w < 4; ++w)
              c += ls[w] == ls[u];
            maxc = std::max(maxc, c);
            if (c >= 3)
              l3of4 = ls[u];
          }
          if (b >= 0 && maxc < 4) {
            if (maxc + nun < 4)
              q_fake4c[b] += 1;
            else
              q_undec4c[b] += 1;
          }
          if (b >= 0) {
            if (l3of4 >= 0)
              q_true3[b] += 1;
            else if (maxc + nun < 3)
              q_fake3[b] += 1;
            else
              q_undec3[b] += 1;
          }
          if (l3of4 >= 0) {
            if (auto it = n_true3.find(l3of4); it != n_true3.end())
              ++it->second;
            if (auto it = n_wide3.find(l3of4); it != n_wide3.end())
              ++it->second;
          }
        }
        if (tru) {
          if (auto it = n_true.find(l0); it != n_true.end())
            ++it->second;
          else if (b >= 0 && l0 < (int)ev.simTracks_.size()) {
            if (ev.simTracks_[l0].pT() < pt_min)
              q_true_lowpt[b] += 1;
            else
              q_true_displ[b] += 1;
          }
        }
      }

      for (const auto &kv : n_wide3) {
        const int b = bin(std::abs(ev.simTracks_[kv.first].momEta()));
        if (b < 0)
          continue;
        den3w[b] += 1;
        if (kv.second > 0)
          found3w[b] += 1;
      }
      for (const auto &kv : n_true3) {
        const int b = bin(std::abs(ev.simTracks_[kv.first].momEta()));
        if (b >= 0 && kv.second > 0) {
          found3[b] += 1;
          dup3_extra[b] += kv.second - 1;
        }
      }
      for (const auto &kv : n_true) {
        const int b = bin(std::abs(ev.simTracks_[kv.first].momEta()));
        if (b < 0)
          continue;
        den[b] += 1;
        if (kv.second > 0) {
          found[b] += 1;
          dup_extra[b] += kv.second - 1;
          if (kv.second > 1)
            found_dup[b] += 1;
        }
      }
    }

    bool write(const std::string &fn, float pt_min, float d0_max) const {
      FILE *f = fopen(fn.c_str(), "w");
      if (!f)
        return false;
      fprintf(f, "# seedfind --truth; pt_min %.3f GeV, d0_max %.3f cm, %.0f events\n", pt_min, d0_max, n_events);
      fprintf(f, "# lo hi  findable found extra_true_quads found_with_dup  quads true fake undecidable"
                 "  true_below_ptmin true_displaced  found3 extra3 true3 fake3 undecidable3  findable3w found3w  fake4c undecidable4c\n");
      for (int b = 0; b < kNb; ++b)
        fprintf(f, "%.2f %.2f  %.0f %.0f %.0f %.0f  %.0f %.0f %.0f %.0f  %.0f %.0f  %.0f %.0f %.0f %.0f %.0f  %.0f %.0f  %.0f %.0f\n",
                b * kEtaMax / kNb, (b + 1) * kEtaMax / kNb, den[b], found[b], dup_extra[b], found_dup[b], q_all[b],
                q_true[b], q_fake[b], q_undec[b], q_true_lowpt[b], q_true_displ[b], found3[b], dup3_extra[b],
                q_true3[b], q_fake3[b], q_undec3[b], den3w[b], found3w[b], q_fake4c[b], q_undec4c[b]);
      fclose(f);
      return true;
    }

    void report() const {
      double D = 0, F = 0, X = 0, A = 0, T = 0, K = 0, U = 0, TL = 0, TD = 0;
      double F3 = 0, X3 = 0, T3 = 0, K3 = 0, U3 = 0, D3W = 0, F3W = 0, K4C = 0, U4C = 0;
      for (int b = 0; b < kNb; ++b) {
        F3 += found3[b], X3 += dup3_extra[b], T3 += q_true3[b], K3 += q_fake3[b], U3 += q_undec3[b];
        D3W += den3w[b], F3W += found3w[b];
        K4C += q_fake4c[b], U4C += q_undec4c[b];
        D += den[b], F += found[b], X += dup_extra[b], TL += q_true_lowpt[b], TD += q_true_displ[b];
        A += q_all[b], T += q_true[b], K += q_fake[b], U += q_undec[b];
      }
      const double ne = n_events > 0 ? n_events : 1;
      printf("   truth: findable %.1f /ev, efficiency %.4f; extra true quads per found track %.3f\n", D / ne,
             F / std::max(1.0, D), X / std::max(1.0, F));
      printf("   truth: quads %.1f /ev: true %.4f, fake %.4f, undecidable %.4f; fake among decidable %.4f\n", A / ne,
             T / std::max(1.0, A), K / std::max(1.0, A), U / std::max(1.0, A), K / std::max(1.0, T + K));
      printf("   truth 3of4: efficiency %.4f; extra matched quads per found track %.3f; quads matched %.4f, fake %.4f, "
             "undecidable %.4f; fake among decidable %.4f\n",
             F3 / std::max(1.0, D), X3 / std::max(1.0, F3), T3 / std::max(1.0, A), K3 / std::max(1.0, A),
             U3 / std::max(1.0, A), K3 / std::max(1.0, T3 + K3));
      printf("   truth 4of4 consistent: quads fake %.4f, undecidable %.4f; fake among decidable %.4f\n",
             K4C / std::max(1.0, A), U4C / std::max(1.0, A), K4C / std::max(1.0, T + K4C));
      printf("   truth 3of4: over tracks with >= 3 of 4 layers (%.1f /ev): efficiency %.4f\n", D3W / ne,
             F3W / std::max(1.0, D3W));
      printf("   truth: tracks passing pT/D0 with 3 of 4 layers, missing a/b/c/d: %.1f %.1f %.1f %.1f /ev; 2 of 4: %.1f /ev\n",
             miss_one[0] / ne, miss_one[1] / ne, miss_one[2] / ne, miss_one[3] / ne, miss_more / ne);
      printf("   truth: true quads of non-findable tracks %.1f /ev: below pT_min %.1f, produced beyond D0_max %.1f\n",
             (TL + TD) / ne, TL / ne, TD / ne);
    }
  };

}  // namespace mkfit::seeding

#endif
